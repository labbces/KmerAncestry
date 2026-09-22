"""
Scan a genome FASTA with sliding windows and compare each window to references.

High-level idea
---------------
This is the final low-disk pipeline:

1. load four precomputed branchwater/sourmash reference ``.sig`` or ``.zip`` files;
2. choose one shared ``scaled`` value, like a branchwater batch comparison;
3. stream the genome FASTA record by record in bounded sequence blocks;
4. hash each block once through the sourmash native API;
5. keep only the sampled hashes that belong to the current sliding window; and
6. write one CSV row per window-vs-reference comparison.

Why the script does not write per-window sketches
-------------------------------------------------
Branchwater is excellent for fast sketching/searching from files, but running
``singlesketch`` for every window would create avoidable I/O and process
overhead. Here the window sketch is represented only by the current in-memory
hash catalog, and the only file written by the pipeline is the final CSV.

Coordinate convention
---------------------
Coordinates in the output are 1-based and inclusive. The step ``d`` advances
the start of the next window. For example, ``w=4`` and ``d=2`` gives windows
``1-4``, ``3-6``, ``5-8``, and so on.

Expected reference files
------------------------
By default, the script searches the reference directory for these species:

* ``officinarum``
* ``robustum``
* ``spontaneum``
* ``barberi``

For each species and k-mer size, it accepts a few common filename variants,
including ``officinarum_31.zip``, ``officinarum31.zip``, and
``officinarum_k31.zip``. Each reference file must contain exactly one DNA
signature for the requested ``k``.

When multiple scaled copies are stored in the same directory, pass
``--scaled``. For example, ``--scaled 500`` selects only files ending in
``_scaled500.sig`` or ``_scaled500.zip`` and also verifies the signature's
internal scaled metadata.

Output metrics
--------------
The final CSV includes both directions of containment, percentage forms of
containment, hash counts, intersection/union sizes, abundance summaries, and
Jaccard, Sorensen-Dice, Bray-Curtis, and cosine similarities plus their
complementary distances. ``average_reference_abundance`` is the mean abundance
of reference hashes shared with a window, analogous to the ``avg_abund`` value
reported by abundance-aware sourmash searches.

Checkpoint and job logs
-----------------------
Use ``--log "description"`` to place a free-text label at the top of both
scheduler ``.o`` and ``.e`` files. Every subsequent script-generated log line
starts with the exact overall completion percentage.

If a job stops, rerun with the same ``--output`` path and pass its prior log to
``--resume-log``. The existing CSV is validated and only missing
window-reference pairs are appended; complete rows are never overwritten.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import sys
from collections import Counter, deque
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Iterable, Iterator, TextIO

# KaMcEs: camada de escrita Parquet com retomada por partes.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "lib"))
import parquetio


DEFAULT_REFERENCE_NAMES = ("officinarum", "robustum", "spontaneum", "barberi")
MINHASH_MAX_HASH = 0xFFFFFFFFFFFFFFFF
DEFAULT_CHECKPOINT_EVERY = 100
DEFAULT_HASH_CHUNK_SIZE = 1_000_000
DEFAULT_OUTPUT_BUFFER_SIZE = 1024 * 1024
DEFAULT_PART_WINDOWS = 5000
CHECKPOINT_COMPLETED_PATTERN = re.compile(r"\bcompleted_windows=(?P<count>\d+)\b")


@dataclass(frozen=True)
class FastaRecordInfo:
    """Identifier and length of one FASTA record, collected in a light first pass."""

    record_id: str
    length: int


@dataclass(frozen=True)
class ResumeCursor:
    """First window that still needs processing in a contiguous completed prefix."""

    record_id: str
    start_1based: int
    completed_windows_before: int


@dataclass
class OutputState:
    """Rows already present in an output CSV, used to resume safely after interruption."""

    written_references_by_window: dict[str, set[str]]
    completed_window_ids: set[str]
    row_count: int


@dataclass
class ProgressLogger:
    """Mirror progress-prefixed lines to the scheduler's stdout and stderr logs."""

    # Cache da decisao de _saidas(). Sem anotacao de tipo DE PROPOSITO: um
    # atributo anotado viraria campo do dataclass.
    _saidas_cache = None

    total_windows: int | None = None
    completed_windows: int = 0
    streams: tuple[TextIO, TextIO] = (sys.stdout, sys.stderr)

    def percentage(self) -> float:
        """Return overall completion percentage, including work restored from a checkpoint."""
        if not self.total_windows:
            return 0.0
        return min(100.0, 100.0 * self.completed_windows / self.total_windows)

    def _saidas(self) -> tuple[TextIO, ...]:
        """Os streams em que realmente escrever.

        O sliding.job usa `#$ -j y`, e com isso o SGE aponta o .o e o .e para o
        MESMO arquivo. Escrever nos dois duplicava TODAS as linhas do log - o
        que apareceu como um log inteiro em dobro, ilegivel e do dobro do
        tamanho. Num terminal e a mesma coisa: stdout e stderr sao o mesmo tty.
        Se os dois descritores apontam para o mesmo arquivo (dev, inode),
        escreve uma vez; se sao arquivos diferentes (sem -j y), escreve nos dois,
        que e a intencao original de a linha aparecer onde voce olhar.
        """
        if self._saidas_cache is not None:
            return self._saidas_cache
        fluxos = tuple(self.streams)
        escolha = fluxos
        if len(fluxos) >= 2:
            try:
                a = os.fstat(fluxos[0].fileno())
                b = os.fstat(fluxos[1].fileno())
                if (a.st_dev, a.st_ino) == (b.st_dev, b.st_ino):
                    escolha = (fluxos[0],)
            except (OSError, ValueError, AttributeError):
                pass          # sem fileno (StringIO nos testes): escreve nos dois
        self._saidas_cache = escolha
        return escolha

    def write_banner(self, message: str | None) -> None:
        """Write the user label followed by exactly two blank lines in both job logs."""
        if message is None:
            return

        safe_message = message.replace("\r", " ").replace("\n", " ")
        for stream in self._saidas():
            print(safe_message, file=stream, flush=True)
            print(file=stream, flush=True)
            print(file=stream, flush=True)

    def write(self, stage: str, message: str, level: str = "INFO") -> None:
        """Write one percentage-prefixed line to both scheduler capture streams."""
        line = f"[{self.percentage():6.2f}%] {level} [{stage}] {message}"
        for stream in self._saidas():
            print(line, file=stream, flush=True)


_LOGGER = ProgressLogger()


@dataclass(frozen=True)
class ReferenceSketch:
    """
    One loaded reference sketch, normalized for fast comparison.

    ``hash_to_weight`` is intentionally kept as a dictionary because the custom
    metrics need direct access to the hash vectors. For non-abundance sketches,
    every hash receives weight 1.0.
    """

    name: str
    path: Path
    ksize: int
    original_scaled: int
    scaled: int
    seed: int | None
    moltype: str | None
    hash_to_weight: dict[int, float]
    uses_abundance: bool

    @property
    def hash_count(self) -> int:
        return len(self.hash_to_weight)

    @cached_property
    def total_weight(self) -> float:
        # O HashCatalog ja soma na construcao; somar 50 milhoes de floats em
        # Python de novo custaria mais que carregar o arquivo.
        pronto = getattr(self.hash_to_weight, "total_weight", None)
        if pronto is not None:
            return pronto
        return sum(self.hash_to_weight.values())

    @cached_property
    def squared_weight_sum(self) -> float:
        pronto = getattr(self.hash_to_weight, "squared_weight_sum", None)
        if pronto is not None:
            return pronto
        return sum(weight * weight for weight in self.hash_to_weight.values())


class HashCatalog:
    """Catalogo de hashes de uma REFERENCIA, em array numpy ordenado.

    POR QUE ISSO EXISTE
    -------------------
    A referencia era um dict[int, float]. Num dict do CPython cada entrada
    custa a chave (int grande: 32 bytes), o ponteiro do valor e o slot da
    tabela: da ordem de 80-100 bytes por hash. Um catalogo de especie de cana
    em scaled=200 tem dezenas de milhoes de hashes - so ele ja passa de 3 GB,
    e sao QUATRO especies carregadas ao mesmo tempo, antes da primeira janela
    ser processada. E o "MemoryError" em 0.00%.

    Aqui os mesmos hashes ficam num uint64 ordenado: 8 bytes cada, mais 8 do
    peso quando ha abundancia. Os mesmos 50 milhoes de hashes passam de ~4 GB
    para 400 MB.

    A busca deixa de ser hash de tabela e passa a ser busca binaria, o que
    seria mais lento se fosse feita um hash de cada vez. Nao e: o compute_
    metrics faz UMA busca vetorizada com todos os hashes da janela de uma vez.
    """

    __slots__ = ("hashes", "weights", "_total", "_sq")

    def __init__(self, hashes, weights=None):
        import numpy as np
        h = np.asarray(hashes, dtype=np.uint64)
        if weights is None:
            ordem = np.argsort(h, kind="stable")
            self.hashes = h[ordem]
            self.weights = None
            self._total = float(self.hashes.size)
            self._sq = float(self.hashes.size)
        else:
            w = np.asarray(weights, dtype=np.float64)
            ordem = np.argsort(h, kind="stable")
            self.hashes = h[ordem]
            self.weights = w[ordem]
            self._total = float(self.weights.sum())
            self._sq = float((self.weights * self.weights).sum())

    def __len__(self):
        return int(self.hashes.size)

    def __bool__(self):
        return self.hashes.size > 0

    def __contains__(self, hash_value):
        import numpy as np
        i = np.searchsorted(self.hashes, np.uint64(hash_value))
        return bool(i < self.hashes.size and self.hashes[i] == np.uint64(hash_value))

    def get(self, hash_value, default=None):
        import numpy as np
        alvo = np.uint64(hash_value)
        i = np.searchsorted(self.hashes, alvo)
        if i >= self.hashes.size or self.hashes[i] != alvo:
            return default
        return 1.0 if self.weights is None else float(self.weights[i])

    def values(self):
        import numpy as np
        return (np.ones(self.hashes.size, dtype=np.float64)
                if self.weights is None else self.weights)

    def items(self):
        """Compatibilidade. Materializa - so use em diagnostico."""
        pesos = self.values()
        for h, w in zip(self.hashes.tolist(), pesos.tolist()):
            yield int(h), float(w)

    @property
    def total_weight(self):
        return self._total

    @property
    def squared_weight_sum(self):
        return self._sq

    def lookup(self, alvos):
        """Busca VETORIZADA. Devolve (achou, pesos) para o array de hashes."""
        import numpy as np
        alvos = np.asarray(alvos, dtype=np.uint64)
        if self.hashes.size == 0 or alvos.size == 0:
            return (np.zeros(alvos.size, dtype=bool),
                    np.zeros(alvos.size, dtype=np.float64))
        i = np.searchsorted(self.hashes, alvos)
        i_seguro = np.minimum(i, self.hashes.size - 1)
        achou = self.hashes[i_seguro] == alvos
        if self.weights is None:
            pesos = achou.astype(np.float64)
        else:
            pesos = np.where(achou, self.weights[i_seguro], 0.0)
        return achou, pesos


@dataclass(frozen=True)
class WindowSketch:
    """Hash catalog for one genome window."""

    record_id: str
    start_1based: int
    end_1based: int
    hash_to_weight: dict[int, float]
    uses_abundance: bool

    @cached_property
    def arrays(self):
        """(hashes, pesos) da janela como arrays numpy, calculado UMA vez.

        A mesma janela e comparada com as 4 referencias; sem este cache os
        arrays seriam remontados 4 vezes por janela, e a montagem e o custo
        dominante numa janela de poucas dezenas de hashes."""
        import numpy as np
        n = len(self.hash_to_weight)
        return (np.fromiter(self.hash_to_weight.keys(), dtype=np.uint64, count=n),
                np.fromiter(self.hash_to_weight.values(), dtype=np.float64, count=n))

    @property
    def window_id(self) -> str:
        return f"{self.record_id}_{self.start_1based}_{self.end_1based}"

    @property
    def hash_count(self) -> int:
        return len(self.hash_to_weight)

    @cached_property
    def total_weight(self) -> float:
        return sum(self.hash_to_weight.values())

    @cached_property
    def squared_weight_sum(self) -> float:
        return sum(weight * weight for weight in self.hash_to_weight.values())


class PipelineError(RuntimeError):
    """Raised when an expected pipeline condition is not met."""


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the final sliding-window workflow."""
    parser = argparse.ArgumentParser(
        description=(
            "Stream a genome FASTA, sketch each sliding window in memory with "
            "sourmash-compatible hashes, and compare against four reference zips."
        )
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Input genome FASTA path. Plain-text and .gz files are supported.",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        required=True,
        help="Directory containing officinarum/robustum/spontaneum/barberi sig sketches.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Final CSV output path.",
    )
    parser.add_argument(
        "-w",
        "--window-size",
        type=int,
        required=True,
        help="Window size in bases.",
    )
    parser.add_argument(
        "-d",
        "--step",
        type=int,
        required=True,
        help="Step between consecutive window starts in bases.",
    )
    parser.add_argument(
        "-k",
        "--ksize",
        type=int,
        required=True,
        help="k-mer size used for sourmash hashing and reference selection.",
    )
    parser.add_argument(
        "--references",
        default=",".join(DEFAULT_REFERENCE_NAMES),
        help=(
            "Comma-separated reference labels. The default is "
            "officinarum,robustum,spontaneum,barberi."
        ),
    )
    parser.add_argument(
        "--scaled",
        type=int,
        default=None,
        help=(
            "Exact scaled value to use. When provided, the script selects only "
            "reference files named with _scaled<value> and verifies that the "
            "signature metadata has that same scaled value."
        ),
    )
    parser.add_argument(
        "--strict-dna",
        action="store_true",
        help=(
            "Fail when a k-mer contains characters outside A/C/G/T. By default, "
            "those k-mers are skipped, matching sourmash add_sequence(force=True)."
        ),
    )
    parser.add_argument(
        "--log-every",
        type=int,
        default=1000,
        help="Print progress to both job logs every N scanned windows. Use 0 to disable.",
    )
    parser.add_argument(
        "--log",
        default=None,
        help=(
            "Free-text label written at the top of both scheduler .o and .e logs, "
            "followed by two blank lines. Quote the string when it contains spaces."
        ),
    )
    parser.add_argument(
        "--resume-log",
        type=Path,
        default=None,
        help=(
            "Optional .o or .e log from an interrupted run. It enables safe append "
            "mode for an existing --output CSV and is inspected for the last checkpoint."
        ),
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=DEFAULT_CHECKPOINT_EVERY,
        help=(
            "Flush and synchronize the output CSV, then write a checkpoint log line, "
            "after this many newly completed windows."
        ),
    )
    parser.add_argument(
        "--part-windows",
        type=int,
        default=DEFAULT_PART_WINDOWS,
        help=(
            "Parquet: janelas por arquivo de parte antes de rotacionar "
            "(padrao: %d). Uma parte so vira legivel quando fechada, entao "
            "este numero e o teto de retrabalho se o job morrer. Valores "
            "muito baixos geram milhares de arquivinhos."
        ) % DEFAULT_PART_WINDOWS,
    )
    parser.add_argument(
        "--compression",
        default="zstd",
        choices=["zstd", "snappy", "gzip", "none"],
        help="Compressao do Parquet (padrao: zstd).",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=9,
        help="Nivel de compressao, quando o codec aceita (padrao: 9).",
    )
    parser.add_argument(
        "--no-sort",
        action="store_true",
        help=(
            "Nao ordenar na consolidacao. Por padrao a saida sai ordenada por "
            "(record_id, start_1based, ordem das referencias), o que garante "
            "que as 4 especies de uma janela fiquem consecutivas mesmo apos "
            "varias retomadas."
        ),
    )
    parser.add_argument(
        "--hash-chunk-size",
        type=int,
        default=DEFAULT_HASH_CHUNK_SIZE,
        help=(
            "Number of FASTA bases hashed per sourmash API call. Larger blocks are "
            "usually faster but use more memory (default: 1000000)."
        ),
    )
    return parser.parse_args()


def configure_logging(message: str | None, streams: tuple[TextIO, TextIO] | None = None) -> None:
    """Reset the process logger and emit the optional user-provided log banner."""
    global _LOGGER
    _LOGGER = ProgressLogger(streams=streams or (sys.stdout, sys.stderr))
    _LOGGER.write_banner(message)


def set_progress(total_windows: int | None = None, completed_windows: int | None = None) -> None:
    """Update the shared progress state used by every subsequent log line."""
    if total_windows is not None:
        _LOGGER.total_windows = total_windows
    if completed_windows is not None:
        _LOGGER.completed_windows = completed_windows


def log(stage: str, message: str, level: str = "INFO") -> None:
    """Write a cluster-friendly status line to both .o and .e scheduler logs."""
    _LOGGER.write(stage, message, level)


def import_sourmash():
    """Import sourmash with a clear cluster/module installation message."""
    try:
        import sourmash
    except ImportError as exc:
        raise PipelineError(
            "sourmash is required. Install it with 'pip install sourmash' or "
            "load the sourmash/branchwater module in the cluster environment."
        ) from exc

    return sourmash


def validate_arguments(args: argparse.Namespace) -> list[str]:
    """Validate paths and numeric arguments before any expensive work starts."""
    if not args.fasta.exists():
        raise PipelineError(f"Input FASTA does not exist: {args.fasta}")
    if not args.fasta.is_file():
        raise PipelineError(f"Input FASTA is not a file: {args.fasta}")
    if not args.reference_dir.exists():
        raise PipelineError(f"Reference directory does not exist: {args.reference_dir}")
    if not args.reference_dir.is_dir():
        raise PipelineError(f"Reference path is not a directory: {args.reference_dir}")
    if args.window_size <= 0:
        raise PipelineError("--window-size/-w must be > 0")
    if args.step <= 0:
        raise PipelineError("--step/-d must be > 0")
    if args.ksize <= 0:
        raise PipelineError("--ksize/-k must be > 0")
    if args.ksize > args.window_size:
        raise PipelineError(
            f"k-mer size ({args.ksize}) cannot be larger than window size "
            f"({args.window_size})."
        )
    if args.scaled is not None and args.scaled <= 0:
        raise PipelineError("--scaled must be > 0 when provided")
    if args.log_every < 0:
        raise PipelineError("--log-every must be >= 0")
    if args.checkpoint_every <= 0:
        raise PipelineError("--checkpoint-every must be > 0")
    if args.hash_chunk_size <= 0:
        raise PipelineError("--hash-chunk-size must be > 0")
    if args.resume_log is not None:
        if not args.resume_log.exists():
            raise PipelineError(f"Resume log does not exist: {args.resume_log}")
        if not args.resume_log.is_file():
            raise PipelineError(f"Resume log is not a file: {args.resume_log}")

    reference_names = [item.strip() for item in args.references.split(",") if item.strip()]
    if not reference_names:
        raise PipelineError("At least one reference name must be provided.")
    if len(set(reference_names)) != len(reference_names):
        raise PipelineError("Reference names in --references must be unique.")
    return reference_names


def open_text_maybe_gzip(path: Path):
    """Open a plain-text or gzipped FASTA file transparently."""
    with path.open("rb") as raw_handle:
        magic = raw_handle.read(2)

    opener = gzip.open if magic == b"\x1f\x8b" else open
    return opener(path, "rt", encoding="utf-8")


def parse_record_id(header_line: str) -> str:
    """Extract the FASTA record ID from a header line without the leading ``>``."""
    record_id = header_line.strip().split()[0] if header_line.strip() else ""
    if not record_id:
        raise PipelineError("Found a FASTA header without a record identifier.")
    return record_id


def iter_fasta_bases(path: Path) -> Iterator[tuple[str, str, int | None, str]]:
    """
    Stream FASTA bases without storing full records or windows on disk.

    Yields
    ------
    ("record_start", record_id, None, "")
        Marks the beginning of a new FASTA record.
    ("base", record_id, position_1based, base)
        One base from the active FASTA record.

    Notes
    -----
    Biopython's ``SeqIO.parse`` is a good iterator over records, but each
    ``SeqRecord`` still holds its sequence. This small parser keeps only the
    active k-mer/window buffers, which is cheaper for chromosome-scale FASTA.
    """
    record_id: str | None = None
    position = 0
    saw_header = False

    with open_text_maybe_gzip(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                record_id = parse_record_id(line[1:])
                position = 0
                saw_header = True
                yield ("record_start", record_id, None, "")
                continue

            if record_id is None:
                raise PipelineError(
                    f"Found sequence data before the first FASTA header at line {line_number}."
                )

            for base in line.upper():
                position += 1
                yield ("base", record_id, position, base)

    if not saw_header:
        raise PipelineError(f"No FASTA records were found in {path}.")


def iter_fasta_sequence_chunks(
    path: Path,
    chunk_size: int,
) -> Iterator[tuple[str, str, int, str]]:
    """
    Stream FASTA records in bounded sequence blocks.

    Hashing one full block per sourmash API call avoids a Python-to-native call
    and a new k-mer string for every base. The yielded sequence coordinates are
    1-based, and records are never joined across FASTA headers.
    """

    current_record_id: str | None = None
    emitted_length = 0
    buffered_parts: list[str] = []
    buffered_length = 0
    saw_header = False

    def take_chunks(final: bool = False) -> Iterator[tuple[str, str, int, str]]:
        nonlocal emitted_length, buffered_parts, buffered_length
        if current_record_id is None or not buffered_parts:
            return

        sequence = "".join(buffered_parts)
        offset = 0
        available = len(sequence)
        while available - offset >= chunk_size:
            block = sequence[offset : offset + chunk_size]
            yield ("sequence", current_record_id, emitted_length + 1, block)
            emitted_length += len(block)
            offset += len(block)

        remainder = sequence[offset:]
        if final and remainder:
            yield ("sequence", current_record_id, emitted_length + 1, remainder)
            emitted_length += len(remainder)
            remainder = ""

        buffered_parts = [remainder] if remainder else []
        buffered_length = len(remainder)

    with open_text_maybe_gzip(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                yield from take_chunks(final=True)
                current_record_id = parse_record_id(line[1:])
                emitted_length = 0
                buffered_parts = []
                buffered_length = 0
                saw_header = True
                yield ("record_start", current_record_id, 0, "")
                continue

            if current_record_id is None:
                raise PipelineError(
                    f"Found sequence data before the first FASTA header at line {line_number}."
                )

            buffered_parts.append(line.upper())
            buffered_length += len(line)
            if buffered_length >= chunk_size:
                yield from take_chunks()

    yield from take_chunks(final=True)
    if not saw_header:
        raise PipelineError(f"No FASTA records were found in {path}.")


def read_fasta_record_info(path: Path) -> list[FastaRecordInfo]:
    """
    Return FASTA record IDs and lengths without materializing their sequences.

    This preflight pass makes the progress percentage exact. It is intentionally
    line based, so it only keeps one record length in memory at a time.
    """

    records: list[FastaRecordInfo] = []
    current_record_id: str | None = None
    current_length = 0

    with open_text_maybe_gzip(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_record_id is not None:
                    records.append(FastaRecordInfo(current_record_id, current_length))
                current_record_id = parse_record_id(line[1:])
                current_length = 0
                continue

            if current_record_id is None:
                raise PipelineError(
                    f"Found sequence data before the first FASTA header at line {line_number}."
                )
            current_length += len(line)

    if current_record_id is not None:
        records.append(FastaRecordInfo(current_record_id, current_length))

    if not records:
        raise PipelineError(f"No FASTA records were found in {path}.")

    record_id_counts = Counter(record.record_id for record in records)
    duplicate_ids = {
        record_id for record_id, count in record_id_counts.items() if count > 1
    }
    if duplicate_ids:
        raise PipelineError(
            "FASTA record IDs must be unique for checkpointing. Duplicates: "
            + ", ".join(sorted(duplicate_ids))
        )

    return records


def window_count_for_record(sequence_length: int, window_size: int, step: int) -> int:
    """Return the number of complete sliding windows available in one FASTA record."""
    if sequence_length < window_size:
        return 0
    return 1 + ((sequence_length - window_size) // step)


def count_complete_windows(
    records: Iterable[FastaRecordInfo],
    window_size: int,
    step: int,
) -> int:
    """Return the exact number of windows that the streaming scan will emit."""
    return sum(
        window_count_for_record(record.length, window_size, step)
        for record in records
    )


def find_resume_cursor(
    records: Iterable[FastaRecordInfo],
    window_size: int,
    step: int,
    completed_window_ids: set[str],
) -> ResumeCursor | None:
    """
    Find the first missing window when existing results form a contiguous prefix.

    A contiguous prefix lets the scanner skip hashing all sequence before that
    point. If a CSV contains completed windows after a gap, the caller receives
    ``None`` and falls back to the slower, but still safe, duplicate-free scan.
    """

    if not completed_window_ids:
        return None

    seen_completed = 0
    for record in records:
        count = window_count_for_record(record.length, window_size, step)
        for index in range(count):
            start_1based = 1 + (index * step)
            end_1based = start_1based + window_size - 1
            window_id = f"{record.record_id}_{start_1based}_{end_1based}"
            if window_id in completed_window_ids:
                seen_completed += 1
                continue

            if seen_completed != len(completed_window_ids):
                return None
            return ResumeCursor(
                record_id=record.record_id,
                start_1based=start_1based,
                completed_windows_before=seen_completed,
            )

    return None


def get_max_hash_for_scaled(scaled: int) -> int:
    """Convert sourmash ``scaled`` into the corresponding hash threshold."""
    if scaled <= 0:
        raise PipelineError("scaled must be > 0 for this pipeline")
    if scaled == 1:
        return MINHASH_MAX_HASH
    return min(int(round(MINHASH_MAX_HASH / scaled, 0)), MINHASH_MAX_HASH)


def make_sourmash_minhash(sourmash, ksize: int, scaled: int, seed: int | None):
    """Create a DNA MinHash using the same seed as the references."""
    kwargs = {"n": 0, "ksize": ksize, "scaled": scaled}
    if seed is not None:
        kwargs["seed"] = seed

    try:
        return sourmash.MinHash(**kwargs)
    except TypeError as exc:
        if seed is None:
            raise
        raise PipelineError(
            "The installed sourmash version did not accept the reference seed "
            f"parameter ({seed}) when creating query sketches."
        ) from exc


def build_reference_candidates(
    reference_dir: Path,
    reference_name: str,
    ksize: int,
    scaled: int | None,
) -> list[Path]:
    """
    Return accepted filename variants for one reference, k-mer size, and scaled value.

    When ``scaled`` is set, only files with the matching ``_scaled<value>``
    suffix are candidates. This is what allows scaled 200, 500, and 1000
    reference collections to live safely in the same directory.
    """

    stems = [
        f"{reference_name}_{ksize}",
        f"{reference_name}{ksize}",
        f"{reference_name}_k{ksize}",
        f"{reference_name}k{ksize}",
    ]
    # Aceita _scaled200 (canonico) e _s200 (nome que o singlesketch usava
    # antes), para nao invalidar catalogos ja gerados.
    suffixes = [f"_scaled{scaled}", f"_s{scaled}"] if scaled is not None else [""]
    return [
        reference_dir / f"{stem}{suffix}{extension}"
        for stem in stems
        for suffix in suffixes
        for extension in (".sig", ".zip")
    ]


def find_available_scaled_values(reference_dir: Path, reference_name: str, ksize: int) -> list[int]:
    """List scaled suffixes available for one reference so ambiguity errors are actionable."""

    name_pattern = re.escape(reference_name)
    ksize_pattern = re.escape(str(ksize))
    pattern = re.compile(
        rf"^{name_pattern}(?:_{ksize_pattern}|{ksize_pattern}|_k{ksize_pattern}|k{ksize_pattern})"
        rf"_s(?:caled)?(?P<scaled>\d+)\.(?:sig|zip)$",
        re.IGNORECASE,
    )
    available = []
    for path in reference_dir.iterdir():
        if not path.is_file():
            continue
        match = pattern.fullmatch(path.name)
        if match:
            available.append(int(match.group("scaled")))
    return sorted(set(available))


def find_reference_file(
    reference_dir: Path,
    reference_name: str,
    ksize: int,
    scaled: int | None,
) -> Path:
    """Find exactly one reference file, using ``--scaled`` to resolve matching copies."""
    candidates = build_reference_candidates(reference_dir, reference_name, ksize, scaled)
    existing = [path for path in candidates if path.exists() and path.is_file()]

    if len(existing) == 1:
        return existing[0]

    candidate_text = ", ".join(path.name for path in candidates)
    if not existing:
        available_scaled = find_available_scaled_values(reference_dir, reference_name, ksize)
        scaled_hint = (
            f" Available scaled values are: {', '.join(map(str, available_scaled))}. "
            "Pass --scaled with one of those values."
            if scaled is None and available_scaled
            else ""
        )
        raise PipelineError(
            f"Could not find reference '{reference_name}' for k={ksize}. "
            f"Tried: {candidate_text}.{scaled_hint}"
        )

    existing_text = ", ".join(path.name for path in existing)
    raise PipelineError(
        f"More than one filename variant exists for reference '{reference_name}' "
        f"and k={ksize}, scaled={scaled}: {existing_text}. Keep only one to avoid ambiguity."
    )


def load_one_reference_signature(sourmash, path: Path, ksize: int):
    """
    Load exactly one DNA signature for ``ksize`` from a sourmash/branchwater file.

    A reference zip with multiple matching signatures would produce more than
    one comparison per species, so the script fails early instead of silently
    merging or choosing the wrong sketch.
    """
    try:
        signatures = list(
            sourmash.load_file_as_signatures(
                str(path),
                ksize=ksize,
                select_moltype="DNA",
            )
        )
    except Exception as exc:
        raise PipelineError(f"Could not load sourmash signatures from {path}: {exc}") from exc

    if not signatures:
        raise PipelineError(f"No DNA signature with k={ksize} found in {path}.")
    if len(signatures) > 1:
        raise PipelineError(
            f"More than one DNA signature with k={ksize} found in {path}. "
            "This pipeline expects one sketch per reference file."
        )
    return signatures[0]


def get_minhash_scaled(minhash, path: Path) -> int:
    """Return and validate the ``scaled`` value from a sourmash MinHash object."""
    scaled = int(getattr(minhash, "scaled", 0) or 0)
    if scaled <= 0:
        raise PipelineError(
            f"Reference {path} is not a scaled/FracMinHash signature. "
            "Branchwater-compatible comparisons require scaled sketches."
        )
    return scaled


def get_minhash_seed(minhash) -> int | None:
    """Return the sourmash seed when exposed by the installed version."""
    seed = getattr(minhash, "seed", None)
    return int(seed) if seed not in (None, "") else None


def get_minhash_moltype(minhash) -> str | None:
    """Return a normalized molecule type when exposed by sourmash."""
    moltype = getattr(minhash, "moltype", None)
    return str(moltype).upper() if moltype not in (None, "") else None


def minhash_weight_map(minhash) -> tuple[dict[int, float], bool]:
    """Convert a sourmash MinHash ``hashes`` object into ``hash -> weight``."""
    raw_hashes = getattr(minhash, "hashes", None)
    if raw_hashes is None:
        raise PipelineError("Loaded MinHash object does not expose a hashes attribute.")

    track_abundance = getattr(minhash, "track_abundance", False)
    if callable(track_abundance):
        track_abundance = track_abundance()

    # Constroi o HashCatalog DIRETO, sem passar por um dict intermediario: o
    # dict e justamente o que nao cabe na memoria num catalogo de especie.
    import numpy as np

    if hasattr(raw_hashes, "items") and track_abundance:
        pares = list(raw_hashes.items())
        h = np.fromiter((int(a) for a, _ in pares), dtype=np.uint64,
                        count=len(pares))
        w = np.fromiter((float(b) for _, b in pares), dtype=np.float64,
                        count=len(pares))
        del pares
        uses_abundance = bool((w != 1.0).any())
        if not uses_abundance:
            return HashCatalog(h), False
        return HashCatalog(h, w), True

    chaves = list(raw_hashes)
    h = np.fromiter((int(v) for v in chaves), dtype=np.uint64, count=len(chaves))
    del chaves
    return HashCatalog(h), False


def downsample_minhash(minhash, scaled: int):
    """Downsample a MinHash to the common scaled value with clear errors."""
    current_scaled = int(getattr(minhash, "scaled", 0) or 0)
    if current_scaled == scaled:
        return minhash
    if current_scaled > scaled:
        raise PipelineError(
            f"Cannot upsample reference from scaled={current_scaled} to {scaled}."
        )

    try:
        return minhash.downsample(scaled=scaled)
    except Exception as exc:
        raise PipelineError(
            f"Could not downsample reference from scaled={current_scaled} "
            f"to scaled={scaled}: {exc}"
        ) from exc


def load_reference_metadata(
    sourmash,
    reference_dir: Path,
    reference_names: Iterable[str],
    ksize: int,
    requested_scaled: int | None,
) -> list[tuple[str, Path, object, int]]:
    """Load raw reference signatures so the common scaled value can be chosen."""
    loaded = []
    for reference_name in reference_names:
        path = find_reference_file(
            reference_dir,
            reference_name,
            ksize,
            requested_scaled,
        )
        signature = load_one_reference_signature(sourmash, path, ksize)
        minhash = signature.minhash
        scaled = get_minhash_scaled(minhash, path)
        if requested_scaled is not None and scaled != requested_scaled:
            raise PipelineError(
                f"Reference {path.name} was selected for --scaled {requested_scaled}, "
                f"but its signature metadata says scaled={scaled}."
            )
        loaded.append((reference_name, path, minhash, scaled))
        log("references", f"loaded {reference_name} from {path.name} with scaled={scaled}")
    return loaded


def normalize_references(
    loaded_references: list[tuple[str, Path, object, int]],
    common_scaled: int,
    ksize: int,
) -> list[ReferenceSketch]:
    """Downsample references to one shared scaled value and convert to maps."""
    references: list[ReferenceSketch] = []
    seeds = set()
    moltypes = set()

    for reference_name, path, minhash, original_scaled in loaded_references:
        if int(getattr(minhash, "ksize", 0) or 0) != ksize:
            raise PipelineError(
                f"Reference {path} has k={getattr(minhash, 'ksize', None)}, "
                f"expected k={ksize}."
            )

        normalized = downsample_minhash(minhash, common_scaled)
        seed = get_minhash_seed(normalized)
        moltype = get_minhash_moltype(normalized)
        hash_to_weight, uses_abundance = minhash_weight_map(normalized)

        if seed is not None:
            seeds.add(seed)
        if moltype is not None:
            moltypes.add(moltype)

        references.append(
            ReferenceSketch(
                name=reference_name,
                path=path,
                ksize=ksize,
                original_scaled=original_scaled,
                scaled=common_scaled,
                seed=seed,
                moltype=moltype,
                hash_to_weight=hash_to_weight,
                uses_abundance=uses_abundance,
            )
        )

    if len(seeds) > 1:
        raise PipelineError(f"Reference sketches use incompatible seeds: {sorted(seeds)}")
    if len(moltypes) > 1:
        raise PipelineError(
            f"Reference sketches use incompatible molecule types: {sorted(moltypes)}"
        )

    return references


def decrement_counter(counter: Counter[int], hash_value: int) -> None:
    """Remove one occurrence from a Counter and delete empty keys."""
    counter[hash_value] -= 1
    if counter[hash_value] <= 0:
        del counter[hash_value]


def sampled_hashes_for_kmer(
    hasher,
    kmer: str,
    max_hash: int,
    strict_dna: bool,
) -> tuple[list[int], bool]:
    """
    Hash one k-mer through sourmash and apply the scaled sampling threshold.

    When ``strict_dna`` is false, invalid DNA k-mers are skipped. This mirrors
    the practical behavior of ``add_sequence(..., force=True)`` and is safer for
    real assemblies containing ``N`` runs.
    """
    try:
        hash_values = hasher.seq_to_hashes(kmer)
    except ValueError as exc:
        if strict_dna:
            raise PipelineError(f"Invalid DNA k-mer '{kmer}': {exc}") from exc
        return [], True

    return [
        int(hash_value)
        for hash_value in hash_values
        if int(hash_value) <= max_hash
    ], False


def make_window_hash_to_weight(
    active_hash_counts: Counter[int],
    use_abundance: bool,
) -> dict[int, float]:
    """Snapshot the active window hashes for metric calculation."""
    if use_abundance:
        return {
            int(hash_value): float(count)
            for hash_value, count in active_hash_counts.items()
            if count > 0
        }

    return {
        int(hash_value): 1.0
        for hash_value, count in active_hash_counts.items()
        if count > 0
    }


def iter_window_sketches(
    sourmash,
    fasta_path: Path,
    window_size: int,
    step: int,
    ksize: int,
    common_scaled: int,
    seed: int | None,
    strict_dna: bool,
    use_abundance: bool,
    resume_cursor: ResumeCursor | None = None,
    hash_chunk_size: int = DEFAULT_HASH_CHUNK_SIZE,
) -> Iterator[WindowSketch]:
    """
    Yield in-memory sketches for every complete sliding window in the FASTA.

    The function never stores a full FASTA record. It keeps:
    * one bounded sequence block and the k-1 bases that precede it;
    * sampled hashes for k-mers still inside the current window; and
    * counters needed to remove hashes when the window moves.
    """
    hasher = make_sourmash_minhash(
        sourmash=sourmash,
        ksize=ksize,
        scaled=common_scaled,
        seed=seed,
    )
    max_hash = get_max_hash_for_scaled(common_scaled)
    active_kmers: deque[tuple[int, int]] = deque()
    active_hash_counts: Counter[int] = Counter()
    sequence_tail = ""
    current_record_id: str | None = None
    next_window_start = 1
    skipped_invalid_kmers = 0
    resume_reached = resume_cursor is None
    process_current_record = resume_cursor is None

    for event, record_id, position, sequence in iter_fasta_sequence_chunks(
        fasta_path,
        hash_chunk_size,
    ):
        if event == "record_start":
            if current_record_id is not None and skipped_invalid_kmers:
                log(
                    "scan",
                    f"{current_record_id}: skipped {skipped_invalid_kmers} invalid k-mer(s)",
                )
            current_record_id = record_id
            skipped_invalid_kmers = 0
            active_kmers.clear()
            active_hash_counts.clear()
            sequence_tail = ""

            if not resume_reached:
                if record_id != resume_cursor.record_id:
                    process_current_record = False
                    continue
                resume_reached = True
                process_current_record = True
                next_window_start = resume_cursor.start_1based
                log(
                    "resume",
                    f"resumed hashing at {record_id}:{next_window_start} after "
                    f"{resume_cursor.completed_windows_before} completed window(s)",
                )
            else:
                process_current_record = True
                next_window_start = 1
                log("scan", f"started FASTA record {current_record_id}")
            continue

        if current_record_id is None:
            raise PipelineError("Internal FASTA parser state became inconsistent.")
        if not process_current_record:
            continue

        chunk_start = int(position)
        chunk_end = chunk_start + len(sequence) - 1
        if chunk_end < next_window_start:
            continue
        if chunk_start < next_window_start:
            sequence = sequence[next_window_start - chunk_start :]
            chunk_start = next_window_start
            sequence_tail = ""

        hash_input = sequence_tail + sequence
        hash_input_start = chunk_start - len(sequence_tail)
        if len(hash_input) < ksize:
            sequence_tail = hash_input
            continue

        if strict_dna:
            invalid_match = re.search(r"[^ACGT]", hash_input)
            if invalid_match:
                invalid_position = hash_input_start + invalid_match.start()
                raise PipelineError(
                    f"Invalid DNA base {hash_input[invalid_match.start()]!r} at "
                    f"{current_record_id}:{invalid_position}."
                )

        try:
            hash_values = hasher.seq_to_hashes(
                hash_input,
                force=not strict_dna,
                bad_kmers_as_zeroes=not strict_dna,
            )
        except ValueError as exc:
            raise PipelineError(
                f"Could not hash {current_record_id}:{hash_input_start}-"
                f"{hash_input_start + len(hash_input) - 1}: {exc}"
            ) from exc

        expected_hash_count = len(hash_input) - ksize + 1
        if len(hash_values) != expected_hash_count:
            raise PipelineError(
                "sourmash returned an unexpected number of DNA hashes: "
                f"{len(hash_values)}, expected {expected_hash_count}."
            )

        for offset, raw_hash_value in enumerate(hash_values):
            kmer_start = hash_input_start + offset
            kmer_end = kmer_start + ksize - 1
            hash_value = int(raw_hash_value)

            if not hash_value:
                skipped_invalid_kmers += 1
            elif hash_value <= max_hash:
                active_kmers.append((kmer_start, hash_value))
                active_hash_counts[hash_value] += 1

            while kmer_end >= next_window_start + window_size - 1:
                while active_kmers and active_kmers[0][0] < next_window_start:
                    _old_start, old_hash = active_kmers.popleft()
                    decrement_counter(active_hash_counts, old_hash)

                end_1based = next_window_start + window_size - 1
                yield WindowSketch(
                    record_id=current_record_id,
                    start_1based=next_window_start,
                    end_1based=end_1based,
                    hash_to_weight=make_window_hash_to_weight(
                        active_hash_counts,
                        use_abundance,
                    ),
                    uses_abundance=use_abundance,
                )
                next_window_start += step

        sequence_tail = hash_input[-(ksize - 1) :] if ksize > 1 else ""

    if current_record_id is not None and skipped_invalid_kmers:
        log("scan", f"{current_record_id}: skipped {skipped_invalid_kmers} invalid k-mer(s)")


def sum_min_weights(left, right) -> float:
    """Sum ``min(weight_left, weight_right)`` across shared hashes."""
    smaller = left.hash_to_weight
    larger = right.hash_to_weight
    if len(smaller) > len(larger):
        smaller, larger = larger, smaller

    total = 0.0
    for hash_value, weight in smaller.items():
        if hash_value in larger:
            total += min(weight, larger[hash_value])
    return total


def bray_curtis_similarity(left, right) -> float:
    """Compute Bray-Curtis similarity from the active hash weights."""
    denominator = left.total_weight + right.total_weight
    if denominator == 0.0:
        return 0.0
    return (2.0 * sum_min_weights(left, right)) / denominator


def cosine_similarity(left, right) -> float:
    """Compute cosine similarity between hash-weight vectors."""
    if not left.hash_to_weight or not right.hash_to_weight:
        return 0.0

    smaller = left.hash_to_weight
    larger = right.hash_to_weight
    if len(smaller) > len(larger):
        smaller, larger = larger, smaller

    dot_product = 0.0
    for hash_value, weight in smaller.items():
        dot_product += weight * larger.get(hash_value, 0.0)

    magnitude_left = math.sqrt(
        sum(weight * weight for weight in left.hash_to_weight.values())
    )
    magnitude_right = math.sqrt(
        sum(weight * weight for weight in right.hash_to_weight.values())
    )
    if magnitude_left == 0.0 or magnitude_right == 0.0:
        return 0.0
    return dot_product / (magnitude_left * magnitude_right)


def shared_weight_sum(source, other) -> float:
    """Return the total source abundance carried by hashes shared with ``other``."""
    return sum(
        weight
        for hash_value, weight in source.hash_to_weight.items()
        if hash_value in other.hash_to_weight
    )


def _interseccao(window, reference):
    """Soma os quatro acumuladores da interseccao janela x referencia.

    Devolve (n, peso_janela, peso_referencia, min, produto).

    Dois caminhos, com o MESMO resultado:
      - vetorizado, quando a referencia e um HashCatalog. Os hashes da janela
        (dezenas a milhares) sao procurados de uma vez so, com busca binaria
        em array. E o caminho normal do pipeline.
      - dicionario, para qualquer outro caso. E a implementacao original,
        mantida como referencia de corretude e usada nos testes para provar
        que o caminho vetorizado da exatamente o mesmo numero.
    """
    janela = window.hash_to_weight
    refer = reference.hash_to_weight

    lookup = getattr(refer, "lookup", None)
    if lookup is not None:
        import numpy as np
        n_j = len(janela)
        if n_j == 0 or len(refer) == 0:
            return 0, 0.0, 0.0, 0.0, 0.0
        prontos = getattr(window, "arrays", None)
        if prontos is not None:
            hj, wj = prontos
        else:
            hj = np.fromiter(janela.keys(), dtype=np.uint64, count=n_j)
            wj = np.fromiter(janela.values(), dtype=np.float64, count=n_j)
        achou, wr = lookup(hj)
        if not achou.any():
            return 0, 0.0, 0.0, 0.0, 0.0
        wj = wj[achou]
        wr = wr[achou]
        return (int(achou.sum()), float(wj.sum()), float(wr.sum()),
                float(np.minimum(wj, wr).sum()), float((wj * wr).sum()))

    menor, maior = janela, refer
    menor_e_janela = True
    if len(menor) > len(maior):
        menor, maior = maior, menor
        menor_e_janela = False

    n = 0
    pj = pr = pmin = prod = 0.0
    for hash_value, peso_menor in menor.items():
        peso_maior = maior.get(hash_value)
        if peso_maior is None:
            continue
        n += 1
        if menor_e_janela:
            w_j, w_r = peso_menor, peso_maior
        else:
            w_j, w_r = peso_maior, peso_menor
        pj += w_j
        pr += w_r
        pmin += min(w_j, w_r)
        prod += w_j * w_r
    return n, pj, pr, pmin, prod


def compute_metrics(window: WindowSketch, reference: ReferenceSketch) -> dict[str, float | int]:
    """Compute overlap, abundance, similarity, and distance metrics for one pair."""
    (intersection_count, shared_window_weight, shared_reference_weight,
     shared_min_weight, dot_product) = _interseccao(window, reference)

    union_count = window.hash_count + reference.hash_count - intersection_count

    window_containment = (
        intersection_count / window.hash_count if window.hash_count else 0.0
    )
    reference_containment = (
        intersection_count / reference.hash_count if reference.hash_count else 0.0
    )

    jaccard = intersection_count / union_count if union_count else 0.0

    total_weight = window.total_weight + reference.total_weight
    bray_curtis = (
        (2.0 * shared_min_weight) / total_weight if total_weight else 0.0
    )
    squared_magnitude_product = window.squared_weight_sum * reference.squared_weight_sum
    cosine = (
        dot_product / math.sqrt(squared_magnitude_product)
        if squared_magnitude_product
        else 0.0
    )

    return {
        "intersection_hash_count": intersection_count,
        "union_hash_count": union_count,
        "window_containment_in_reference": window_containment,
        "reference_containment_in_window": reference_containment,
        "max_containment": max(window_containment, reference_containment),
        "shared_window_weight": shared_window_weight,
        "shared_reference_weight": shared_reference_weight,
        "average_window_abundance": (
            shared_window_weight / intersection_count if intersection_count else 0.0
        ),
        "average_reference_abundance": (
            shared_reference_weight / intersection_count if intersection_count else 0.0
        ),
        "jaccard_similarity": jaccard,
        "jaccard_distance": 1.0 - jaccard,
        "bray_curtis_similarity": bray_curtis,
        "bray_curtis_distance": 1.0 - bray_curtis,
        "cosine_similarity": cosine,
        "cosine_distance": 1.0 - cosine,
    }


def format_float(value: float) -> float:
    """Devolve a metrica como float64.

    A versao CSV formatava com f"{value:.10f}", o que em Parquet criaria 19
    colunas de TEXTO - e ai o log2 do R e a entropia teriam que converter tudo
    de volta. Em Parquet o tipo e float64 nativo, com precisao total.
    Na saida CSV o formato de 10 casas continua sendo aplicado, agora na hora
    de escrever (parquetio.RowWriter.float_format).
    """
    return float(value)


def write_csv_header(writer) -> None:
    """Cabecalho: no CSV e escrito pelo RowWriter na abertura; no Parquet o
    schema faz esse papel. Mantido para nao mudar as chamadas existentes."""
    writer.writeheader()


def build_output_row(
    window: WindowSketch,
    reference: ReferenceSketch,
    window_size: int,
    step: int,
    ksize: int,
) -> dict[str, str | int]:
    """Build one CSV row with the full comparison and abundance summary."""
    metrics = compute_metrics(window, reference)
    return {
        "window_id": window.window_id,
        "record_id": window.record_id,
        "window_size": window_size,
        "step": step,
        "ksize": ksize,
        "start_1based": window.start_1based,
        "end_1based": window.end_1based,
        "reference_name": reference.name,
        "reference_file": str(reference.path),
        "reference_scaled_original": reference.original_scaled,
        "common_scaled": reference.scaled,
        "window_hash_count": window.hash_count,
        "reference_hash_count": reference.hash_count,
        "window_total_abundance": format_float(window.total_weight),
        "reference_total_abundance": format_float(reference.total_weight),
        "intersection_hash_count": metrics["intersection_hash_count"],
        "union_hash_count": metrics["union_hash_count"],
        "window_containment_in_reference": format_float(
            metrics["window_containment_in_reference"]
        ),
        "reference_containment_in_window": format_float(
            metrics["reference_containment_in_window"]
        ),
        "window_containment_percent": format_float(
            100.0 * metrics["window_containment_in_reference"]
        ),
        "reference_containment_percent": format_float(
            100.0 * metrics["reference_containment_in_window"]
        ),
        "max_containment": format_float(metrics["max_containment"]),
        "shared_window_weight": format_float(metrics["shared_window_weight"]),
        "shared_reference_weight": format_float(metrics["shared_reference_weight"]),
        "average_window_abundance": format_float(metrics["average_window_abundance"]),
        "average_reference_abundance": format_float(
            metrics["average_reference_abundance"]
        ),
        "uses_abundance": window.uses_abundance or reference.uses_abundance,
        "jaccard_similarity": format_float(metrics["jaccard_similarity"]),
        "jaccard_distance": format_float(metrics["jaccard_distance"]),
        "bray_curtis_similarity": format_float(metrics["bray_curtis_similarity"]),
        "bray_curtis_distance": format_float(metrics["bray_curtis_distance"]),
        "cosine_similarity": format_float(metrics["cosine_similarity"]),
        "cosine_distance": format_float(metrics["cosine_distance"]),
    }


OUTPUT_FIELDNAMES = [
    "window_id",
    "record_id",
    "window_size",
    "step",
    "ksize",
    "start_1based",
    "end_1based",
    "reference_name",
    "reference_file",
    "reference_scaled_original",
    "common_scaled",
    "window_hash_count",
    "reference_hash_count",
    "window_total_abundance",
    "reference_total_abundance",
    "intersection_hash_count",
    "union_hash_count",
    "window_containment_in_reference",
    "reference_containment_in_window",
    "window_containment_percent",
    "reference_containment_percent",
    "max_containment",
    "shared_window_weight",
    "shared_reference_weight",
    "average_window_abundance",
    "average_reference_abundance",
    "uses_abundance",
    "jaccard_similarity",
    "jaccard_distance",
    "bray_curtis_similarity",
    "bray_curtis_distance",
    "cosine_similarity",
    "cosine_distance",
]


def read_last_checkpoint_from_log(path: Path) -> int | None:
    """Read the final ``completed_windows`` value written by a prior job run."""
    last_count: int | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = CHECKPOINT_COMPLETED_PATTERN.search(line)
            if match:
                last_count = int(match.group("count"))
    return last_count


def inspect_existing_output(
    output_path: Path,
    fieldnames: list[str],
    reference_names: set[str],
    args: argparse.Namespace,
    common_scaled: int,
) -> OutputState:
    """
    Validate an existing CSV and index rows already written by an interrupted job.

    Every ``(window_id, reference_name)`` pair must occur at most once. A window
    is marked complete only after all expected references are present, allowing a
    resume to fill an incomplete final window without rewriting any CSV rows.
    """

    if parquetio.is_parquet(output_path):
        # A comparacao de configuracao e feita como INTEIRO. A versao CSV
        # comparava str(...) porque tudo em CSV e texto; em Parquet as colunas
        # voltam tipadas e a comparacao por string reprovaria todas as linhas.
        expected_config = {
            "window_size": args.window_size,
            "step": args.step,
            "ksize": args.ksize,
            "common_scaled": common_scaled,
        }
        try:
            # Devolve so as janelas INCOMPLETAS no dicionario; as completas vem
            # num set separado. Guardar um set de referencias por janela
            # completa custava mais de 1 GB num cromossomo grande - era isso
            # que estourava a memoria antes da primeira janela ser processada.
            partial, completed, row_count = parquetio.inspect_existing(
                output_path, fieldnames, reference_names, expected_config, log=log)
        except RuntimeError as exc:
            raise PipelineError(f"Cannot resume {output_path}: {exc}") from exc
        return OutputState(
            written_references_by_window=partial,
            completed_window_ids=completed,
            row_count=row_count,
        )

    if not output_path.is_file():
        raise PipelineError(f"Existing output path is not a file: {output_path}")

    written_references_by_window: dict[str, set[str]] = {}
    seen_pairs: set[tuple[str, str]] = set()
    row_count = 0
    expected_config = {
        "window_size": str(args.window_size),
        "step": str(args.step),
        "ksize": str(args.ksize),
        "common_scaled": str(common_scaled),
    }

    with output_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != fieldnames:
            raise PipelineError(
                f"Cannot resume {output_path}: CSV columns do not match this script version. "
                "Choose a new --output path for a fresh run."
            )

        for line_number, row in enumerate(reader, start=2):
            row_count += 1
            window_id = row.get("window_id", "")
            reference_name = row.get("reference_name", "")
            if not window_id or not reference_name:
                raise PipelineError(
                    f"Cannot resume {output_path}: missing window_id or reference_name "
                    f"on CSV line {line_number}."
                )
            if reference_name not in reference_names:
                raise PipelineError(
                    f"Cannot resume {output_path}: reference '{reference_name}' on CSV "
                    "does not belong to the current --references selection."
                )
            for column, expected_value in expected_config.items():
                if row.get(column) != expected_value:
                    raise PipelineError(
                        f"Cannot resume {output_path}: CSV line {line_number} has "
                        f"{column}={row.get(column)!r}, expected {expected_value!r}."
                    )

            pair = (window_id, reference_name)
            if pair in seen_pairs:
                raise PipelineError(
                    f"Cannot resume {output_path}: duplicate comparison {pair} on "
                    f"or before CSV line {line_number}."
                )
            seen_pairs.add(pair)
            written_references_by_window.setdefault(window_id, set()).add(reference_name)

    completed_window_ids = {
        window_id
        for window_id, written_references in written_references_by_window.items()
        if written_references == reference_names
    }
    return OutputState(
        written_references_by_window=written_references_by_window,
        completed_window_ids=completed_window_ids,
        row_count=row_count,
    )


def open_output_for_run(
    args: argparse.Namespace,
    fieldnames: list[str],
    reference_names: set[str],
    common_scaled: int,
):
    """Abre a saida (Parquet em partes, ou CSV) e indexa o que ja existe.

    Devolve (writer, writer, state): o mesmo objeto faz o papel do handle e do
    DictWriter, entao run_pipeline segue chamando handle.close() e
    writer.writerow() sem alteracao.
    """
    args.output.parent.mkdir(parents=True, exist_ok=True)
    parts_dir = args.output.with_name(args.output.name + parquetio.PARTS_SUFFIX)
    has_previous = args.output.exists() or parts_dir.is_dir()

    if has_previous and args.resume_log is None:
        raise PipelineError(
            f"Output already exists: {args.output}. Pass --resume-log <old .o/.e log> "
            "to append only missing results, or choose a different --output path."
        )

    state = OutputState({}, set(), 0)
    if has_previous:
        state = inspect_existing_output(
            output_path=args.output,
            fieldnames=fieldnames,
            reference_names=reference_names,
            args=args,
            common_scaled=common_scaled,
        )

    writer = parquetio.RowWriter(
        args.output,
        fieldnames,
        part_rows=max(1, args.part_windows) * max(1, len(reference_names)),
        compression=args.compression,
        compression_level=args.compression_level,
    )
    return writer, writer, state


def synchronize_output(handle) -> None:
    """Ponto de checkpoint duravel.

    No CSV: flush + fsync, como antes. No Parquet: esvazia o buffer e fecha a
    parte corrente quando ela ja esta cheia - so depois de fechada ela tem
    rodape e pode ser lida por uma retomada.
    """
    try:
        handle.checkpoint()
    except OSError as exc:
        raise PipelineError(f"Could not synchronize checkpoint output: {exc}") from exc


def _rss_mb() -> float:
    """Memoria residente do processo, em MB. 0 se a plataforma nao disser."""
    try:
        import resource
        kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return kb / 1024.0            # Linux devolve KB
    except Exception:
        return 0.0


def _mem_log(references) -> None:
    total_hashes = 0
    for r in references:
        n = r.hash_count
        total_hashes += n
        com_peso = getattr(r.hash_to_weight, "weights", "ausente") is not None
        bytes_por = 16 if com_peso else 8
        log("memoria",
            f"{r.name}: {n:,} hashes (~{n * bytes_por / 1e6:.0f} MB)")
    log("memoria",
        f"{len(references)} referencia(s), {total_hashes:,} hashes no total; "
        f"RSS do processo: {_rss_mb():.0f} MB")


def run_pipeline(args: argparse.Namespace) -> int:
    """Run the complete workflow after argument parsing."""
    reference_names = validate_arguments(args)
    sourmash = import_sourmash()

    log(
        "start",
        (
            f"fasta={args.fasta} reference_dir={args.reference_dir} "
            f"w={args.window_size} d={args.step} k={args.ksize}"
        ),
    )

    loaded_references = load_reference_metadata(
        sourmash=sourmash,
        reference_dir=args.reference_dir,
        reference_names=reference_names,
        ksize=args.ksize,
        requested_scaled=args.scaled,
    )
    reference_scaled_values = [scaled for _name, _path, _minhash, scaled in loaded_references]
    common_scaled = args.scaled if args.scaled is not None else max(reference_scaled_values)
    log("references", f"using common scaled={common_scaled}")

    references = normalize_references(
        loaded_references=loaded_references,
        common_scaled=common_scaled,
        ksize=args.ksize,
    )
    reference_seeds = {reference.seed for reference in references if reference.seed is not None}
    common_seed = next(iter(reference_seeds), None) if reference_seeds else None
    use_abundance = any(reference.uses_abundance for reference in references)
    if common_seed is not None:
        log("references", f"using reference seed={common_seed}")
    if use_abundance:
        log("references", "at least one reference stores abundance; window counts will be used")

    fasta_records = read_fasta_record_info(args.fasta)
    # Contabilidade de memoria, dita em voz alta. Quando um job morre com
    # MemoryError, o log passa a dizer QUANTO cada peca custava - em vez de
    # deixar todo mundo adivinhando entre catalogo grande, janelas demais e
    # vizinho no mesmo no.
    _mem_log(references)

    total_window_count = count_complete_windows(
        fasta_records,
        args.window_size,
        args.step,
    )
    if total_window_count == 0:
        raise PipelineError(
            f"No complete windows of size {args.window_size} were found in {args.fasta}."
        )
    set_progress(total_windows=total_window_count, completed_windows=0)
    log("preflight", f"found {total_window_count} complete sliding window(s)")

    expected_reference_names = {reference.name for reference in references}
    prior_log_checkpoint = (
        read_last_checkpoint_from_log(args.resume_log)
        if args.resume_log is not None
        else None
    )
    handle, writer, output_state = open_output_for_run(
        args=args,
        fieldnames=OUTPUT_FIELDNAMES,
        reference_names=expected_reference_names,
        common_scaled=common_scaled,
    )

    if len(output_state.completed_window_ids) > total_window_count:
        handle.close()
        raise PipelineError(
            f"Output has {len(output_state.completed_window_ids)} completed windows, "
            f"but this FASTA configuration contains only {total_window_count}."
        )

    set_progress(completed_windows=len(output_state.completed_window_ids))
    if args.resume_log is not None:
        prior_checkpoint_text = (
            str(prior_log_checkpoint) if prior_log_checkpoint is not None else "not found"
        )
        log(
            "resume",
            f"loaded {output_state.row_count} existing row(s), including "
            f"{len(output_state.completed_window_ids)} completed window(s); "
            f"last checkpoint in {args.resume_log.name}: {prior_checkpoint_text}",
        )

    resume_cursor = find_resume_cursor(
        fasta_records,
        args.window_size,
        args.step,
        output_state.completed_window_ids,
    )
    if (
        output_state.completed_window_ids
        and len(output_state.completed_window_ids) < total_window_count
        and resume_cursor is None
    ):
        log(
            "resume",
            "existing completed windows are not a contiguous prefix; running a full "
            "safe scan and writing only missing window-reference pairs",
        )

    window_count_this_run = 0
    row_count = output_state.row_count
    zero_hash_windows = 0
    windows_since_checkpoint = 0

    try:
        if len(output_state.completed_window_ids) == total_window_count:
            # TODAS as janelas ja foram calculadas - e isso NAO quer dizer que o
            # arquivo final exista. Quando o trabalho esta nas PARTES, sair aqui
            # deixava a saida SEM CONSOLIDAR: o .parquet nunca aparecia, o
            # checar_sliding classificava como RESTO (partes sem final), o
            # refazer re-submetia, e a tarefa gastava dezenas de GB de RAM e
            # minutos de leitura para dizer "nada a acrescentar" outra vez, em
            # looping, sem nunca gravar nada. Consolidar aqui e o que fecha o
            # ciclo: a consolidacao junta as partes (mais um final preexistente,
            # se houver) e remove a pasta de partes.
            summary = handle.consolidate(
                reference_order=[reference.name for reference in references],
                sort=not args.no_sort,
            )
            if summary.get("parts"):
                log(
                    "consolidate",
                    f"{summary['parts']} parte(s) -> {args.output} "
                    f"({summary['rows']} linhas"
                    + (", ordenadas)" if summary.get("sorted") else ")"),
                )
                log("done", "nada a acrescentar: todas as janelas ja estavam "
                            f"calculadas; saida consolidada em {args.output}")
            else:
                log("done", "a saida ja contem todas as comparacoes esperadas; "
                            "nada a acrescentar e nada a consolidar")
            return 0

        for window in iter_window_sketches(
            sourmash=sourmash,
            fasta_path=args.fasta,
            window_size=args.window_size,
            step=args.step,
            ksize=args.ksize,
            common_scaled=common_scaled,
            seed=common_seed,
            strict_dna=args.strict_dna,
            use_abundance=use_abundance,
            resume_cursor=resume_cursor,
            hash_chunk_size=args.hash_chunk_size,
        ):
            window_count_this_run += 1
            # Janela ja completa nao entra no dicionario (ver inspect_existing):
            # a consulta tem de vir ANTES do setdefault, senao ele criaria um
            # set vazio e a janela seria recalculada.
            if window.window_id in output_state.completed_window_ids:
                continue
            written_references = output_state.written_references_by_window.setdefault(
                window.window_id,
                set(),
            )
            if written_references == expected_reference_names:
                continue

            if window.hash_count == 0:
                zero_hash_windows += 1
                if zero_hash_windows <= 10:
                    log(
                        "window",
                        f"{window.window_id} has no sampled hashes; writing zero metrics",
                    )
                elif zero_hash_windows == 11:
                    log(
                        "window",
                        "more zero-hash windows found; suppressing repeated warnings",
                    )

            for reference in references:
                if reference.name in written_references:
                    continue
                writer.writerow(
                    build_output_row(
                        window=window,
                        reference=reference,
                        window_size=args.window_size,
                        step=args.step,
                        ksize=args.ksize,
                    )
                )
                written_references.add(reference.name)
                row_count += 1

            if written_references == expected_reference_names:
                output_state.completed_window_ids.add(window.window_id)
                # Tira a janela do dicionario assim que ela fecha. Sem isto o
                # dicionario cresce ate ter uma entrada com um set por janela do
                # cromossomo - mais de 1 GB num cromossomo de 1,1 milhao de
                # janelas, acumulado ao longo da propria execucao. Agora so as
                # janelas em andamento ficam nele.
                output_state.written_references_by_window.pop(window.window_id, None)
                set_progress(completed_windows=len(output_state.completed_window_ids))
                windows_since_checkpoint += 1

            if windows_since_checkpoint >= args.checkpoint_every:
                synchronize_output(handle)
                log(
                    "checkpoint",
                    f"completed_windows={len(output_state.completed_window_ids)} "
                    f"total_windows={total_window_count} csv_rows={row_count}",
                )
                windows_since_checkpoint = 0

            if args.log_every and window_count_this_run % args.log_every == 0:
                log(
                    "progress",
                    f"scanned {window_count_this_run} window(s) in this run and "
                    f"wrote {row_count} total CSV row(s)",
                )

        expected_windows_this_run = (
            total_window_count - resume_cursor.completed_windows_before
            if resume_cursor is not None
            else total_window_count
        )
        if window_count_this_run != expected_windows_this_run:
            raise PipelineError(
                f"Expected to scan {expected_windows_this_run} window(s), but scanned "
                f"{window_count_this_run}. The FASTA changed during execution or the "
                "checkpoint no longer matches it."
            )
        if len(output_state.completed_window_ids) != total_window_count:
            raise PipelineError(
                f"Checkpoint ended with {len(output_state.completed_window_ids)} completed "
                f"windows, expected {total_window_count}."
            )

        synchronize_output(handle)
        log(
            "checkpoint",
            f"completed_windows={len(output_state.completed_window_ids)} "
            f"total_windows={total_window_count} csv_rows={row_count}",
        )
        log("done", f"completed {total_window_count} window(s)")

        # Consolidacao: junta as partes num unico arquivo com o nome pedido.
        # So acontece no caminho de sucesso - um job que morre deixa as partes
        # intactas para a proxima retomada.
        summary = handle.consolidate(
            reference_order=[reference.name for reference in references],
            sort=not args.no_sort,
        )
        if summary.get("parts"):
            log(
                "consolidate",
                f"{summary['parts']} parte(s) -> {args.output} "
                f"({summary['rows']} linhas"
                + (", ordenadas)" if summary.get("sorted") else ")"),
            )
        log("done", f"wrote {row_count} comparison rows to {args.output}")
        if zero_hash_windows:
            log("done", f"{zero_hash_windows} newly processed window(s) had no sampled hashes")
        return 0
    finally:
        handle.close()


def main() -> int:
    """CLI entrypoint with cluster-friendly error reporting."""
    args = parse_args()
    configure_logging(args.log)
    try:
        return run_pipeline(args)
    except PipelineError as exc:
        log("pipeline", str(exc), level="ERROR")
        return 1
    except KeyboardInterrupt:
        log("pipeline", "interrupted by user or scheduler signal", level="ERROR")
        return 130
    except Exception as exc:
        log("unexpected", f"{type(exc).__name__}: {exc}", level="ERROR")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
