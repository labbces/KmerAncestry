"""
Generate validation FASTA subsets from a genome using a sliding-window scan.

Overview
--------
This script was written to help validate a future TCC pipeline that will scan
genomic sequences with windows of size ``w`` and step ``d``. For now, the goal
is not to exhaustively export every possible window, but to confirm that the
logic works correctly on a real FASTA file.

Default validation behavior
---------------------------
By default, the script:

1. reads a genome FASTA with Biopython;
2. scans the sequence(s) with a sliding window step of 10 bp;
3. collects 3 windows for each target size:
   * 10^2 bp
   * 10^4 bp
   * 10^6 bp
4. writes one output FASTA per window size; and
5. writes a TSV manifest with the exact coordinates of every exported subset.

Coordinate convention
---------------------
Internally, Python slices use the usual ``start_0based`` and
``end_0based_exclusive`` coordinates. For the output metadata, the script also
stores:

* ``start_1based``
* ``end_1based_inclusive``

This avoids ambiguity when checking the exact length of each subset. For
example, a 10,000 bp window that starts at the 3rd base has:

* ``start_1based = 3``
* ``end_1based_inclusive = 10002``
* ``end_0based_exclusive = 10002``

The FASTA header keeps the concise part in the record ID itself, for example:

``>SP803280_w10000_3_10002 ...``

Example
-------
Run with the default validation sizes and default step:

``python genome_sliding_windows.py --fasta /path/to/genome.fasta``

If you want to target only one record inside a multi-record FASTA:

``python genome_sliding_windows.py --fasta genome.fasta --record-id SP803280``
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path


DEFAULT_WINDOW_SIZES_TEXT = "10^2,10^4,10^6"
WINDOW_POWER_PATTERN = re.compile(r"10(?:\^|\*\*)(?P<exponent>\d+)$")


@dataclass(frozen=True)
class FastaSequence:
    """
    In-memory representation of one FASTA record.

    Attributes
    ----------
    record_id
        Short record identifier from the FASTA header.
    description
        Full FASTA description line, useful for provenance in the manifest.
    sequence
        Nucleotide sequence exactly as Biopython returned it, converted to a
        plain Python string for easier slicing.
    """

    record_id: str
    description: str
    sequence: str


@dataclass(frozen=True)
class WindowSubset:
    """
    One exported sliding-window subset.

    Attributes
    ----------
    source_record_id
        FASTA record from which this window was extracted.
    source_description
        Full description line of the source FASTA record.
    window_size
        Number of base pairs in the subset.
    subset_index
        Running index within each window size group, starting at 1.
    start_0based
        Inclusive start coordinate in Python's 0-based coordinate system.
    end_0based_exclusive
        Exclusive end coordinate in Python's slice convention.
    sequence
        Extracted subsequence.
    """

    source_record_id: str
    source_description: str
    window_size: int
    subset_index: int
    start_0based: int
    end_0based_exclusive: int
    sequence: str

    @property
    def start_1based(self) -> int:
        """Return the inclusive 1-based start coordinate."""
        return self.start_0based + 1

    @property
    def end_1based_inclusive(self) -> int:
        """Return the inclusive 1-based end coordinate."""
        return self.end_0based_exclusive

    @property
    def length(self) -> int:
        """Return the number of bases in the exported subset."""
        return len(self.sequence)

    @property
    def header_id(self) -> str:
        """
        Build a compact FASTA identifier with size and genomic coordinates.

        Example
        -------
        ``SP803280_w10000_3_10002``
        """

        return (
            f"{self.source_record_id}_w{self.window_size}_"
            f"{self.start_1based}_{self.end_1based_inclusive}"
        )

    @property
    def header_description(self) -> str:
        """Return a verbose description to keep all coordinate systems visible."""
        return (
            f"subset_index={self.subset_index} "
            f"window_size={self.window_size} "
            f"start_0based={self.start_0based} "
            f"end_0based_exclusive={self.end_0based_exclusive} "
            f"start_1based={self.start_1based} "
            f"end_1based_inclusive={self.end_1based_inclusive} "
            f"length={self.length}"
        )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the sliding-window validation workflow."""
    parser = argparse.ArgumentParser(
        description=(
            "Read a genome FASTA with Biopython and export validation subsets "
            "generated with a sliding window."
        )
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Input genome FASTA path. Plain-text and .gz files are supported.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("sliding_window_subsets"),
        help="Directory where the subset FASTA files and manifest will be written.",
    )
    parser.add_argument(
        "--window-sizes",
        default=DEFAULT_WINDOW_SIZES_TEXT,
        help=(
            "Comma-separated window sizes. Values can be written as integers "
            "(100,10000,1000000) or powers of ten (10^2,10^4,10^6)."
        ),
    )
    parser.add_argument(
        "--step",
        type=int,
        default=10,
        help="Sliding-window step in base pairs.",
    )
    parser.add_argument(
        "--windows-per-size",
        type=int,
        default=3,
        help="How many subsets should be exported for each window size.",
    )
    parser.add_argument(
        "--record-id",
        default=None,
        help=(
            "Optional FASTA record ID to extract from a multi-record file. "
            "If omitted, records are scanned in file order."
        ),
    )
    parser.add_argument(
        "--manifest-name",
        default="window_manifest.tsv",
        help="Filename of the TSV manifest written inside --output-dir.",
    )
    return parser.parse_args()


def parse_window_size_token(token: str) -> int:
    """
    Convert one size token into an integer number of bases.

    Accepted forms
    --------------
    * plain integers, e.g. ``10000``
    * power notation, e.g. ``10^4`` or ``10**4``
    """

    normalized = token.strip().replace("_", "").replace(" ", "")
    if not normalized:
        raise ValueError("Window size tokens cannot be empty.")

    if normalized.isdigit():
        value = int(normalized)
    else:
        match = WINDOW_POWER_PATTERN.fullmatch(normalized)
        if not match:
            raise ValueError(
                f"Unsupported window size token: '{token}'. "
                "Use integers such as 10000 or power notation such as 10^4."
            )
        value = 10 ** int(match.group("exponent"))

    if value <= 0:
        raise ValueError("Window sizes must be positive integers.")
    return value


def parse_window_sizes(raw_value: str) -> list[int]:
    """Parse a comma-separated window-size list while preserving the user order."""
    sizes: list[int] = []
    seen: set[int] = set()

    for token in raw_value.split(","):
        size = parse_window_size_token(token)
        if size not in seen:
            sizes.append(size)
            seen.add(size)

    if not sizes:
        raise ValueError("At least one window size must be provided.")

    return sizes


def open_text_maybe_gzip(path: Path):
    """Open a plain-text or gzipped FASTA file transparently."""
    with path.open("rb") as raw_handle:
        magic = raw_handle.read(2)

    opener = gzip.open if magic == b"\x1f\x8b" else open
    return opener(path, "rt", encoding="utf-8")


def import_biopython_seqio():
    """Import ``Bio.SeqIO`` with a clearer installation message if it is missing."""
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise ImportError(
            "Biopython is required to read the FASTA file. "
            "Install it with 'pip install biopython' or load the corresponding "
            "module in the cluster environment."
        ) from exc

    return SeqIO


def read_fasta_records(path: Path, selected_record_id: str | None = None) -> list[FastaSequence]:
    """
    Read FASTA records with Biopython and convert them into plain dataclasses.

    Notes
    -----
    If ``selected_record_id`` is provided, only the matching FASTA record is
    returned. Otherwise, every record is loaded in file order.
    """

    SeqIO = import_biopython_seqio()
    records: list[FastaSequence] = []

    with open_text_maybe_gzip(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if selected_record_id is not None and record.id != selected_record_id:
                continue

            records.append(
                FastaSequence(
                    record_id=record.id,
                    description=record.description,
                    sequence=str(record.seq),
                )
            )

    if selected_record_id is not None and not records:
        raise ValueError(
            f"Record '{selected_record_id}' was not found in FASTA file {path}."
        )

    if not records:
        raise ValueError(f"No FASTA records were found in {path}.")

    return records


def iter_window_starts(sequence_length: int, window_size: int, step: int):
    """
    Yield all valid 0-based window starts for one sequence.

    A start is valid when the complete window fits inside the sequence without
    truncation.
    """

    if window_size <= 0:
        raise ValueError("window_size must be > 0")
    if step <= 0:
        raise ValueError("step must be > 0")

    max_start = sequence_length - window_size
    if max_start < 0:
        return

    for start_0based in range(0, max_start + 1, step):
        yield start_0based


def collect_validation_windows(
    records: list[FastaSequence],
    window_size: int,
    step: int,
    windows_per_size: int,
) -> list[WindowSubset]:
    """
    Collect the first ``windows_per_size`` windows available for one size.

    The records are scanned in the order they appear in the FASTA file. This is
    useful for validation because the logic stays deterministic and easy to
    inspect.
    """

    if windows_per_size <= 0:
        raise ValueError("windows_per_size must be > 0")

    collected: list[WindowSubset] = []

    for record in records:
        for start_0based in iter_window_starts(len(record.sequence), window_size, step):
            end_0based_exclusive = start_0based + window_size
            collected.append(
                WindowSubset(
                    source_record_id=record.record_id,
                    source_description=record.description,
                    window_size=window_size,
                    subset_index=len(collected) + 1,
                    start_0based=start_0based,
                    end_0based_exclusive=end_0based_exclusive,
                    sequence=record.sequence[start_0based:end_0based_exclusive],
                )
            )

            if len(collected) == windows_per_size:
                return collected

    raise ValueError(
        f"Could not collect {windows_per_size} window(s) of size {window_size} "
        f"with step {step}. Only {len(collected)} valid window(s) were found "
        f"across {len(records)} FASTA record(s)."
    )


def wrap_sequence(sequence: str, width: int = 80):
    """Yield FASTA-friendly sequence lines with a fixed width."""
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def write_windows_to_fasta(windows: list[WindowSubset], output_path: Path) -> None:
    """Write a list of subsets to a multi-record FASTA file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        for window in windows:
            handle.write(f">{window.header_id} {window.header_description}\n")
            for line in wrap_sequence(window.sequence):
                handle.write(f"{line}\n")


def build_manifest_rows(windows: list[WindowSubset], output_path: Path) -> list[dict[str, str | int]]:
    """Convert exported windows into tabular rows for the TSV manifest."""
    rows: list[dict[str, str | int]] = []

    for window in windows:
        rows.append(
            {
                "output_fasta": str(output_path),
                "subset_header": window.header_id,
                "source_record_id": window.source_record_id,
                "source_description": window.source_description,
                "window_size": window.window_size,
                "subset_index": window.subset_index,
                "start_0based": window.start_0based,
                "end_0based_exclusive": window.end_0based_exclusive,
                "start_1based": window.start_1based,
                "end_1based_inclusive": window.end_1based_inclusive,
                "length": window.length,
            }
        )

    return rows


def write_manifest(rows: list[dict[str, str | int]], output_path: Path) -> None:
    """Write the final manifest as a tab-separated table."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "output_fasta",
        "subset_header",
        "source_record_id",
        "source_description",
        "window_size",
        "subset_index",
        "start_0based",
        "end_0based_exclusive",
        "start_1based",
        "end_1based_inclusive",
        "length",
    ]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def strip_fasta_suffix(path: Path) -> str:
    """
    Remove common FASTA and compression suffixes to build cleaner output names.

    Examples
    --------
    * ``genome.fasta`` -> ``genome``
    * ``genome.fa.gz`` -> ``genome``
    """

    name = path.name
    if name.endswith(".gz"):
        name = Path(name).stem
    return Path(name).stem


def make_output_fasta_path(output_dir: Path, fasta_path: Path, window_size: int) -> Path:
    """Build the FASTA filename for one window size group."""
    prefix = strip_fasta_suffix(fasta_path)
    return output_dir / f"{prefix}_w{window_size}_subsets.fasta"


def validate_arguments(step: int, windows_per_size: int) -> None:
    """Validate the numeric parameters that control the sliding-window scan."""
    if step <= 0:
        raise ValueError("--step must be > 0")
    if windows_per_size <= 0:
        raise ValueError("--windows-per-size must be > 0")


def main() -> int:
    """
    Orchestrate the full subset-generation workflow.

    Processing sequence
    -------------------
    1. Parse command-line arguments.
    2. Read the FASTA file with Biopython.
    3. Collect validation windows for each requested size.
    4. Write one FASTA file per size.
    5. Write a TSV manifest with the exact coordinates.
    """

    args = parse_args()

    if not args.fasta.exists():
        print(f"Input FASTA does not exist: {args.fasta}", file=sys.stderr)
        return 1

    try:
        validate_arguments(args.step, args.windows_per_size)
        window_sizes = parse_window_sizes(args.window_sizes)
        records = read_fasta_records(args.fasta, args.record_id)
    except (ImportError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    manifest_rows: list[dict[str, str | int]] = []

    try:
        for window_size in window_sizes:
            windows = collect_validation_windows(
                records=records,
                window_size=window_size,
                step=args.step,
                windows_per_size=args.windows_per_size,
            )
            output_fasta_path = make_output_fasta_path(
                output_dir=args.output_dir,
                fasta_path=args.fasta,
                window_size=window_size,
            )
            write_windows_to_fasta(windows, output_fasta_path)
            manifest_rows.extend(build_manifest_rows(windows, output_fasta_path))
            print(
                f"Wrote {len(windows)} subset(s) of size {window_size} to "
                f"{output_fasta_path}"
            )
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    manifest_path = args.output_dir / args.manifest_name
    write_manifest(manifest_rows, manifest_path)
    print(f"Wrote manifest with {len(manifest_rows)} row(s) to {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
