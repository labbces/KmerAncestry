#!/usr/bin/env python3
"""
Scan a genome FASTA with sliding windows and compare each window to references.

High-level idea
---------------
This is the final low-disk pipeline:

1. load four precomputed branchwater/sourmash reference zip files;
2. choose one shared ``scaled`` value, like a branchwater batch comparison;
3. stream the genome FASTA record by record and base by base;
4. hash each k-mer once with the sourmash Python API;
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
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


DEFAULT_REFERENCE_NAMES = ("officinarum", "robustum", "spontaneum", "barberi")
MINHASH_MAX_HASH = 0xFFFFFFFFFFFFFFFF


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
    scaled: int
    seed: int | None
    moltype: str | None
    hash_to_weight: dict[int, float]
    uses_abundance: bool

    @property
    def hash_count(self) -> int:
        return len(self.hash_to_weight)

    @property
    def total_weight(self) -> float:
        return sum(self.hash_to_weight.values())


@dataclass(frozen=True)
class WindowSketch:
    """Hash catalog for one genome window."""

    record_id: str
    start_1based: int
    end_1based: int
    hash_to_weight: dict[int, float]

    @property
    def window_id(self) -> str:
        return f"{self.record_id}_{self.start_1based}_{self.end_1based}"

    @property
    def hash_count(self) -> int:
        return len(self.hash_to_weight)

    @property
    def total_weight(self) -> float:
        return sum(self.hash_to_weight.values())


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
        help="Directory containing officinarum/robustum/spontaneum/barberi zip sketches.",
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
            "Optional minimum scaled value for query windows. The actual common "
            "scaled is max(--scaled, reference scaled values)."
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
        help="Print progress to stderr every N completed windows. Use 0 to disable.",
    )
    return parser.parse_args()


def log(stage: str, message: str) -> None:
    """Write a cluster-friendly status line to stderr."""
    print(f"INFO [{stage}] {message}", file=sys.stderr, flush=True)


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

    reference_names = [item.strip() for item in args.references.split(",") if item.strip()]
    if not reference_names:
        raise PipelineError("At least one reference name must be provided.")
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


def build_reference_candidates(reference_dir: Path, reference_name: str, ksize: int) -> list[Path]:
    """Return accepted filename variants for one reference and k-mer size."""
    filenames = [
        f"{reference_name}_{ksize}.zip",
        f"{reference_name}{ksize}.zip",
        f"{reference_name}_k{ksize}.zip",
        f"{reference_name}k{ksize}.zip",
    ]
    return [reference_dir / filename for filename in filenames]


def find_reference_file(reference_dir: Path, reference_name: str, ksize: int) -> Path:
    """Find one reference zip file, accepting common naming variants."""
    candidates = build_reference_candidates(reference_dir, reference_name, ksize)
    existing = [path for path in candidates if path.exists() and path.is_file()]

    if len(existing) == 1:
        return existing[0]

    candidate_text = ", ".join(path.name for path in candidates)
    if not existing:
        raise PipelineError(
            f"Could not find reference '{reference_name}' for k={ksize}. "
            f"Tried: {candidate_text}"
        )

    existing_text = ", ".join(path.name for path in existing)
    raise PipelineError(
        f"More than one filename variant exists for reference '{reference_name}' "
        f"and k={ksize}: {existing_text}. Keep only one to avoid ambiguity."
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

    if hasattr(raw_hashes, "items"):
        if not track_abundance:
            return {int(hash_value): 1.0 for hash_value in raw_hashes}, False

        hash_to_weight = {
            int(hash_value): float(weight)
            for hash_value, weight in raw_hashes.items()
        }
        uses_abundance = any(weight != 1.0 for weight in hash_to_weight.values())
        return hash_to_weight, uses_abundance

    return {int(hash_value): 1.0 for hash_value in raw_hashes}, False


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
) -> list[tuple[str, Path, object, int]]:
    """Load raw reference signatures so the common scaled value can be chosen."""
    loaded = []
    for reference_name in reference_names:
        path = find_reference_file(reference_dir, reference_name, ksize)
        signature = load_one_reference_signature(sourmash, path, ksize)
        minhash = signature.minhash
        scaled = get_minhash_scaled(minhash, path)
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

    for reference_name, path, minhash, _scaled in loaded_references:
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
) -> Iterator[WindowSketch]:
    """
    Yield in-memory sketches for every complete sliding window in the FASTA.

    The function never stores a full FASTA record. It keeps:
    * the current k-mer buffer;
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
    kmer_chars: deque[str] = deque(maxlen=ksize)
    active_kmers: deque[tuple[int, int]] = deque()
    active_hash_counts: Counter[int] = Counter()
    current_record_id: str | None = None
    next_window_start = 1
    skipped_invalid_kmers = 0

    for event, record_id, position, base in iter_fasta_bases(fasta_path):
        if event == "record_start":
            if current_record_id is not None and skipped_invalid_kmers:
                log(
                    "scan",
                    f"{current_record_id}: skipped {skipped_invalid_kmers} invalid k-mer(s)",
                )
            current_record_id = record_id
            next_window_start = 1
            skipped_invalid_kmers = 0
            kmer_chars.clear()
            active_kmers.clear()
            active_hash_counts.clear()
            log("scan", f"started FASTA record {current_record_id}")
            continue

        if current_record_id is None or position is None:
            raise PipelineError("Internal FASTA parser state became inconsistent.")

        kmer_chars.append(base)
        if len(kmer_chars) == ksize:
            kmer = "".join(kmer_chars)
            kmer_start = int(position) - ksize + 1
            sampled, was_invalid = sampled_hashes_for_kmer(
                hasher=hasher,
                kmer=kmer,
                max_hash=max_hash,
                strict_dna=strict_dna,
            )
            if was_invalid:
                skipped_invalid_kmers += 1
            for hash_value in sampled:
                active_kmers.append((kmer_start, hash_value))
                active_hash_counts[hash_value] += 1

        while int(position) >= next_window_start + window_size - 1:
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
            )
            next_window_start += step

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


def compute_metrics(window: WindowSketch, reference: ReferenceSketch) -> dict[str, float]:
    """Compute the requested similarity metrics in the requested order."""
    window_hashes = set(window.hash_to_weight)
    reference_hashes = set(reference.hash_to_weight)
    intersection_count = len(window_hashes & reference_hashes)
    union_count = len(window_hashes | reference_hashes)
    dice_denominator = window.hash_count + reference.hash_count

    jaccard = intersection_count / union_count if union_count else 0.0
    sorensen_dice = (
        (2.0 * intersection_count) / dice_denominator
        if dice_denominator
        else 0.0
    )

    return {
        "jaccard_similarity": jaccard,
        "sorensen_dice_similarity": sorensen_dice,
        "bray_curtis_similarity": bray_curtis_similarity(window, reference),
        "cosine_similarity": cosine_similarity(window, reference),
    }


def format_float(value: float) -> str:
    """Format floating-point output with stable precision for CSV export."""
    return f"{value:.10f}"


def write_csv_header(writer: csv.DictWriter) -> None:
    """Write the CSV header once, before the streaming loop starts."""
    writer.writeheader()


def build_output_row(
    window: WindowSketch,
    reference: ReferenceSketch,
    window_size: int,
    step: int,
) -> dict[str, str | int]:
    """Build one final CSV row for a window/reference comparison."""
    metrics = compute_metrics(window, reference)
    return {
        "window_id": window.window_id,
        "window_size": window_size,
        "step": step,
        "start_1based": window.start_1based,
        "end_1based": window.end_1based,
        "comparison_file": reference.name,
        "jaccard_similarity": format_float(metrics["jaccard_similarity"]),
        "sorensen_dice_similarity": format_float(metrics["sorensen_dice_similarity"]),
        "bray_curtis_similarity": format_float(metrics["bray_curtis_similarity"]),
        "cosine_similarity": format_float(metrics["cosine_similarity"]),
    }


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
    )
    reference_scaled_values = [scaled for _name, _path, _minhash, scaled in loaded_references]
    common_scaled = max(reference_scaled_values + ([args.scaled] if args.scaled else []))
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

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "window_id",
        "window_size",
        "step",
        "start_1based",
        "end_1based",
        "comparison_file",
        "jaccard_similarity",
        "sorensen_dice_similarity",
        "bray_curtis_similarity",
        "cosine_similarity",
    ]

    window_count = 0
    row_count = 0
    zero_hash_windows = 0

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        write_csv_header(writer)

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
        ):
            window_count += 1
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
                writer.writerow(
                    build_output_row(
                        window=window,
                        reference=reference,
                        window_size=args.window_size,
                        step=args.step,
                    )
                )
                row_count += 1

            if args.log_every and window_count % args.log_every == 0:
                log(
                    "progress",
                    f"processed {window_count} windows and wrote {row_count} CSV rows",
                )

    if window_count == 0:
        raise PipelineError(
            f"No complete windows of size {args.window_size} were found in {args.fasta}."
        )

    log("done", f"processed {window_count} windows")
    log("done", f"wrote {row_count} comparison rows to {args.output}")
    if zero_hash_windows:
        log("done", f"{zero_hash_windows} window(s) had no sampled hashes")
    return 0


def main() -> int:
    """CLI entrypoint with cluster-friendly error reporting."""
    args = parse_args()
    try:
        return run_pipeline(args)
    except PipelineError as exc:
        print(f"ERROR [pipeline] {exc}", file=sys.stderr, flush=True)
        return 1
    except KeyboardInterrupt:
        print("ERROR [pipeline] interrupted by user or scheduler signal", file=sys.stderr)
        return 130
    except Exception as exc:
        print(f"ERROR [unexpected] {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
