#!/usr/bin/env python3
"""
Compare precomputed FracKMC/sourmash-compatible sketches and export a CSV.

High-level idea
---------------
This script compares one set of sample sketches against one genome sketch for
each k-mer size. It was designed for cases where:

1. the sample sketches already exist on disk and should not be recomputed;
2. the sample sketches may be stored either directly as ``.sig`` files or
   wrapped inside branchwater-style ``.zip`` containers;
3. the genome sketch is available as a normal ``.sig`` file;
4. multiple k-mer sizes exist in the same directory; and
5. the output must be a single CSV containing one row per sample-vs-genome
   comparison.

Why this script is a little more careful than a raw set comparison
------------------------------------------------------------------
FracMinHash/sourmash sketches can only be compared safely after they are
harmonized to a common ``scaled`` value. Branchwater does this internally
before computing overlap-based quantities such as containment and Jaccard.

To mimic that behavior more closely, this script:

1. reads the metadata stored inside each signature;
2. extracts the original ``scaled`` value;
3. chooses a common scaled value for the comparison;
4. downsamples the hashes so both sides live in the same hash space; and only
   then
5. computes the requested metrics.

Output
------
For each valid pair (sample sketch, genome sketch) at the same k-mer size and
software label, the
script writes a CSV row containing:

* original and harmonized sketch sizes;
* containment values;
* Jaccard similarity/distance;
* cosine similarity/distance;
* Sorensen-Dice similarity/distance; and
* Bray-Curtis similarity/distance.

If you want a more tutorial-style explanation of each function, metric, and
processing step, see the companion README file:
``compare_precomputed_sketches_README.md``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import re
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path


K_PATTERN = re.compile(r"_k(?P<k>\d+)(?:\D|$)")
SOFTWARE_PATTERN = re.compile(
    r"_k\d+_(?P<software>sour|sourmash|branch|branchwater)$",
    re.IGNORECASE,
)
MINHASH_MAX_HASH = 0xFFFFFFFFFFFFFFFF


@dataclass(frozen=True)
class SignatureData:
    """
    In-memory representation of one selected signature.

    A sourmash/FracKMC file can contain metadata plus one or more signatures.
    After the script selects the exact signature for a given k-mer size, it is
    normalized into this dataclass so the rest of the code can work with a
    predictable Python object instead of raw JSON.

    Attributes
    ----------
    path
        Original file path where the signature came from.
    label
        Human-readable name stored in the signature metadata, or a fallback
        derived from the filename.
    ksize
        K-mer size of the selected signature.
    scaled
        FracMinHash scaled value used by the sketch.
    max_hash
        Maximum hash value implied by ``scaled``. Sourmash uses ``scaled`` and
        ``max_hash`` as two related ways of expressing the same sampling rule.
    seed
        Hash seed used by the sketch, when available.
    moltype
        Molecule type recorded in the signature, e.g. DNA.
    hash_to_weight
        Mapping from hash value to abundance/weight. For non-abundance sketches
        every retained hash receives weight 1.0.
    uses_abundance
        Whether the signature stores abundance information.
    """
    path: Path
    label: str
    ksize: int
    scaled: int
    max_hash: int
    seed: int | None
    moltype: str | None
    hash_to_weight: dict[int, float]
    uses_abundance: bool

    @property
    def hashes(self) -> set[int]:
        return set(self.hash_to_weight)

    @property
    def hash_count(self) -> int:
        return len(self.hash_to_weight)

    @property
    def total_weight(self) -> float:
        return sum(self.hash_to_weight.values())

    def downsample(self, new_scaled: int) -> "SignatureData":
        """
        Return a new signature downsampled to ``new_scaled``.

        Sourmash compares scaled sketches after moving both sides to the same
        effective sampling density. For scaled signatures, "downsampling" means
        keeping only hashes whose numeric value is below the threshold implied
        by the new scaled value.

        Important rule:
        * a larger scaled value keeps fewer hashes;
        * therefore we can downsample from 1000 to 2000, but not "upsample"
          from 2000 back to 1000.
        """
        if new_scaled < self.scaled:
            raise ValueError(
                f"Cannot upsample signature from scaled={self.scaled} to {new_scaled}"
            )

        if new_scaled == self.scaled:
            return self

        new_max_hash = get_max_hash_for_scaled(new_scaled)
        filtered_hash_to_weight = {
            hash_value: weight
            for hash_value, weight in self.hash_to_weight.items()
            if hash_value <= new_max_hash
        }

        return SignatureData(
            path=self.path,
            label=self.label,
            ksize=self.ksize,
            scaled=new_scaled,
            max_hash=new_max_hash,
            seed=self.seed,
            moltype=self.moltype,
            hash_to_weight=filtered_hash_to_weight,
            uses_abundance=self.uses_abundance,
        )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the sketch comparison workflow."""
    parser = argparse.ArgumentParser(
        description=(
            "Compare precomputed scaled .sig files from samples against a genome "
            "sketch with the same k-mer size and software label."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory that contains the input sketch files.",
    )
    parser.add_argument(
        "--sample-pattern",
        default="SRR*_k*.sig,SRR*_k*.zip",
        help="Glob pattern(s) used to find sample sketches. Separate patterns with commas.",
    )
    parser.add_argument(
        "--genome-pattern",
        default="SP803280_subset_k*.sig",
        help=(
            "Glob pattern used to find genome sketches. Genome filenames must "
            "also include the software suffix, e.g. '_sour' or '_branch'."
        ),
    )
    parser.add_argument(
        "--common-scaled-mode",
        choices=("per-k-max", "pairwise-max"),
        default="per-k-max",
        help=(
            "How to harmonize scaled values: 'per-k-max' mimics a single "
            "branchwater-like run per k, while 'pairwise-max' harmonizes each "
            "sample/genome pair independently."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("sample_vs_genome_similarity.csv"),
        help="Output CSV path.",
    )
    return parser.parse_args()


def extract_ksize_from_name(path: Path) -> int:
    """Extract ``k`` from filenames such as ``SRR123_k21.sig``."""
    match = K_PATTERN.search(path.stem)
    if not match:
        raise ValueError(f"Could not extract k-mer size from file name: {path.name}")
    return int(match.group("k"))


def strip_k_suffix(path: Path) -> str:
    """Remove the ``_k...`` suffix so filenames become cleaner sample labels."""
    return re.sub(r"_k\d+.*$", "", path.stem)


def strip_known_sketch_extensions(path: Path) -> str:
    """
    Remove known sketch/container extensions from a filename.

    This helper exists because ``Path.stem`` only removes the last suffix. For
    files such as ``example.sig.gz`` we want to strip both suffixes before
    parsing metadata encoded in the basename.
    """
    name = path.name
    for suffix in (".sig.gz", ".sig", ".zip"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def extract_software_from_name(path: Path) -> str:
    """
    Extract the software label encoded in a sample filename.

    Expected patterns include:

    * ``SRR..._k21_sour.sig``
    * ``SRR..._k21_branch.zip``

    For convenience, longer aliases are normalized as follows:

    * ``sourmash`` -> ``sour``
    * ``branchwater`` -> ``branch``
    """
    basename = strip_known_sketch_extensions(path)
    match = SOFTWARE_PATTERN.search(basename)
    if not match:
        raise ValueError(
            "Could not extract software from sample filename. Expected a suffix "
            f"such as '_sour' or '_branch' in {path.name}"
        )

    software = match.group("software").lower()
    if software == "sourmash":
        return "sour"
    if software == "branchwater":
        return "branch"
    return software


def get_max_hash_for_scaled(scaled: int) -> int:
    """
    Convert sourmash ``scaled`` into the corresponding maximum accepted hash.

    Sourmash retains hashes according to a threshold in the 64-bit hash space.
    Internally, the two most useful representations are:

    * ``scaled``: human-friendly sampling density;
    * ``max_hash``: exact numeric cutoff in hash space.

    This function mirrors the documented sourmash conversion so our manual
    downsampling behaves like sourmash/branchwater.
    """
    if scaled < 0:
        raise ValueError("scaled must be >= 0")
    if scaled == 0:
        return 0
    if scaled == 1:
        return MINHASH_MAX_HASH
    return min(int(round(MINHASH_MAX_HASH / scaled, 0)), MINHASH_MAX_HASH)


def get_scaled_for_max_hash(max_hash: int) -> int:
    """Invert ``get_max_hash_for_scaled`` when only ``max_hash`` is present."""
    if max_hash < 0:
        raise ValueError("max_hash must be >= 0")
    if max_hash == 0:
        return 0
    return min(int(round(MINHASH_MAX_HASH / max_hash, 0)), MINHASH_MAX_HASH)


def open_maybe_gzip(path: Path):
    """Open a plain-text or gzipped signature file transparently."""
    with path.open("rb") as raw_handle:
        magic = raw_handle.read(2)

    opener = gzip.open if magic == b"\x1f\x8b" else open
    return opener(path, "rt", encoding="utf-8")


def load_json(path: Path):
    """Load JSON from either a direct signature file or a branchwater zip."""
    if path.suffix.lower() == ".zip":
        return load_json_from_zip(path)

    with open_maybe_gzip(path) as handle:
        return json.load(handle)


def load_json_from_zip(path: Path):
    """
    Read the single ``.sig``/``.sig.gz`` payload from a branchwater-style zip.

    Branchwater often stores signatures in a zip containing:

    * one compressed signature under ``signatures/...sig.gz``; and
    * a manifest CSV.

    This helper locates that one signature payload and returns its decoded JSON
    content.
    """
    with zipfile.ZipFile(path) as archive:
        signature_members = [
            member_info
            for member_info in archive.infolist()
            if not member_info.is_dir()
            and (
                member_info.filename.endswith(".sig")
                or member_info.filename.endswith(".sig.gz")
            )
        ]

        if not signature_members:
            raise ValueError(f"No .sig or .sig.gz file found inside {path}")

        if len(signature_members) > 1:
            member_names = ", ".join(member.filename for member in signature_members)
            raise ValueError(
                f"More than one .sig-like file found inside {path}: {member_names}"
            )

        member_name = signature_members[0].filename
        raw_bytes = archive.read(member_name)
        if member_name.endswith(".gz"):
            raw_bytes = gzip.decompress(raw_bytes)

        with io.TextIOWrapper(io.BytesIO(raw_bytes), encoding="utf-8") as handle:
            return json.load(handle)


def iter_signature_entries(payload, path: Path):
    """
    Iterate over all signature-like records in the decoded JSON payload.

    Sourmash JSON can appear in slightly different shapes depending on how it
    was written. This helper normalizes those structures so the rest of the
    script can simply iterate over ``(label, signature_dict)`` pairs.
    """
    if isinstance(payload, dict):
        entries = [payload]
    elif isinstance(payload, list):
        entries = payload
    else:
        raise ValueError(f"Unsupported JSON structure in {path}")

    for entry in entries:
        entry_label = entry.get("name") or entry.get("filename") or strip_k_suffix(path)
        signatures = entry.get("signatures")

        if signatures is None and ("mins" in entry or "hashes" in entry):
            signatures = [entry]

        if signatures is None:
            continue

        for signature in signatures:
            yield entry_label, signature


def build_hash_weight_map(signature: dict) -> tuple[dict[int, float], bool]:
    """
    Convert raw sourmash fields into a single ``hash -> weight`` mapping.

    Supported cases:
    * plain ``mins`` list with no abundance tracking;
    * ``mins`` plus ``abundances`` list; and
    * dictionary-like structures found in some serialized variants.
    """
    mins = signature.get("mins")
    abundances = signature.get("abundances")
    hashes = signature.get("hashes")

    if isinstance(mins, dict):
        return {int(hash_value): float(weight) for hash_value, weight in mins.items()}, True

    if isinstance(hashes, dict):
        hash_to_weight = {}
        uses_abundance = False
        for hash_value, weight in hashes.items():
            numeric_weight = float(weight)
            hash_to_weight[int(hash_value)] = numeric_weight
            if numeric_weight != 1.0:
                uses_abundance = True
        return hash_to_weight, uses_abundance

    if mins is None:
        raise ValueError("Signature does not contain 'mins' or 'hashes'.")

    if abundances is None:
        return {int(hash_value): 1.0 for hash_value in mins}, False

    if len(mins) != len(abundances):
        raise ValueError("Signature has different numbers of mins and abundances.")

    hash_to_weight = {}
    for hash_value, abundance in zip(mins, abundances):
        hash_to_weight[int(hash_value)] = float(abundance)
    return hash_to_weight, True


def build_scaled_metadata(signature: dict) -> tuple[int, int]:
    """
    Extract and validate the scaled metadata from one signature.

    A scaled sourmash signature may store:

    * ``scaled`` directly;
    * ``max_hash`` directly; or
    * both.

    This function reconciles those fields and ensures they are internally
    consistent.
    """
    scaled_raw = signature.get("scaled")
    max_hash_raw = signature.get("max_hash")

    scaled = int(scaled_raw) if scaled_raw not in (None, "") else 0
    max_hash = int(max_hash_raw) if max_hash_raw not in (None, "") else 0

    if scaled and max_hash:
        expected_max_hash = get_max_hash_for_scaled(scaled)
        if max_hash != expected_max_hash:
            raise ValueError(
                "Signature contains inconsistent 'scaled' and 'max_hash' values."
            )
        return scaled, max_hash

    if scaled:
        return scaled, get_max_hash_for_scaled(scaled)

    if max_hash:
        return get_scaled_for_max_hash(max_hash), max_hash

    raise ValueError(
        "Only scaled/FracMinHash signatures are supported. Missing 'scaled' or 'max_hash'."
    )


def load_signature(path: Path, expected_ksize: int) -> SignatureData:
    """
    Load exactly one signature for the requested ``k`` from a file.

    The script intentionally refuses ambiguous cases where multiple signatures
    with the same k-mer size live in the same file. Failing early is safer than
    silently picking the wrong record.
    """
    payload = load_json(path)

    matches = []
    for label, signature in iter_signature_entries(payload, path):
        signature_ksize = signature.get("ksize")
        if signature_ksize is None:
            continue
        if int(signature_ksize) != expected_ksize:
            continue
        matches.append((label, signature))

    if not matches:
        raise ValueError(f"No signature with k={expected_ksize} found in file {path}")

    if len(matches) > 1:
        raise ValueError(
            f"More than one signature with k={expected_ksize} found in file {path}"
        )

    label, signature = matches[0]
    hash_to_weight, uses_abundance = build_hash_weight_map(signature)
    scaled, max_hash = build_scaled_metadata(signature)

    moltype = signature.get("molecule") or signature.get("moltype")
    seed_raw = signature.get("seed")
    seed = int(seed_raw) if seed_raw not in (None, "") else None

    return SignatureData(
        path=path,
        label=label,
        ksize=expected_ksize,
        scaled=scaled,
        max_hash=max_hash,
        seed=seed,
        moltype=moltype,
        hash_to_weight=hash_to_weight,
        uses_abundance=uses_abundance,
    )


def ensure_compatible(sig1: SignatureData, sig2: SignatureData) -> None:
    """
    Check metadata that should match before two signatures are compared.

    Even if two sketches share the same filename pattern, they should not be
    compared if fundamental sketch settings differ, such as k-mer size, seed,
    or molecule type.
    """
    if sig1.ksize != sig2.ksize:
        raise ValueError(
            f"Incompatible ksize: {sig1.path.name} has k={sig1.ksize}, "
            f"but {sig2.path.name} has k={sig2.ksize}"
        )

    if sig1.seed is not None and sig2.seed is not None and sig1.seed != sig2.seed:
        raise ValueError(
            f"Incompatible seed: {sig1.path.name} has seed={sig1.seed}, "
            f"but {sig2.path.name} has seed={sig2.seed}"
        )

    if sig1.moltype and sig2.moltype and sig1.moltype != sig2.moltype:
        raise ValueError(
            f"Incompatible molecule type: {sig1.path.name} has {sig1.moltype}, "
            f"but {sig2.path.name} has {sig2.moltype}"
        )


def count_intersection_hashes(sig1: SignatureData, sig2: SignatureData) -> int:
    """Count shared hash identities between two already-harmonized sketches."""
    return len(sig1.hashes & sig2.hashes)


def count_union_hashes(sig1: SignatureData, sig2: SignatureData) -> int:
    """Count distinct hash identities seen in at least one of the two sketches."""
    return len(sig1.hashes | sig2.hashes)


def sum_min_weights(sig1: SignatureData, sig2: SignatureData) -> float:
    """
    Sum ``min(weight1, weight2)`` across shared hashes.

    This is the weighted-overlap term used in Bray-Curtis-style calculations.
    """
    smaller = sig1.hash_to_weight
    larger = sig2.hash_to_weight
    if len(smaller) > len(larger):
        smaller, larger = larger, smaller

    total = 0.0
    for hash_value, weight in smaller.items():
        if hash_value in larger:
            total += min(weight, larger[hash_value])
    return total


def cosine_similarity(sig1: SignatureData, sig2: SignatureData) -> float:
    """
    Compute cosine similarity using the weight vector of each sketch.

    If abundance tracking is absent, every hash has weight 1 and the metric
    reduces to cosine similarity on binary presence/absence vectors.
    """
    if not sig1.hash_to_weight or not sig2.hash_to_weight:
        return 0.0

    smaller = sig1.hash_to_weight
    larger = sig2.hash_to_weight
    if len(smaller) > len(larger):
        smaller, larger = larger, smaller

    dot_product = 0.0
    for hash_value, weight in smaller.items():
        dot_product += weight * larger.get(hash_value, 0.0)

    magnitude1 = math.sqrt(sum(weight * weight for weight in sig1.hash_to_weight.values()))
    magnitude2 = math.sqrt(sum(weight * weight for weight in sig2.hash_to_weight.values()))

    if magnitude1 == 0.0 or magnitude2 == 0.0:
        return 0.0

    return dot_product / (magnitude1 * magnitude2)


def jaccard_similarity(sig1: SignatureData, sig2: SignatureData) -> float:
    """Compute unweighted Jaccard similarity from shared vs total hashes."""
    union_size = count_union_hashes(sig1, sig2)
    if union_size == 0:
        return 0.0
    return count_intersection_hashes(sig1, sig2) / union_size


def sorensen_dice_similarity(sig1: SignatureData, sig2: SignatureData) -> float:
    """Compute the binary Sorensen-Dice coefficient."""
    denominator = sig1.hash_count + sig2.hash_count
    if denominator == 0:
        return 0.0
    return (2.0 * count_intersection_hashes(sig1, sig2)) / denominator


def bray_curtis_similarity(sig1: SignatureData, sig2: SignatureData) -> float:
    """
    Compute a Bray-Curtis-style similarity using hash abundances.

    Classical Bray-Curtis is often presented as a dissimilarity:

        BC_distance = 1 - (2 * sum(min) / (sum_x + sum_y))

    Here we expose both forms:
    * this function returns the similarity term;
    * the CSV also reports the complementary distance.
    """
    denominator = sig1.total_weight + sig2.total_weight
    if denominator == 0.0:
        return 0.0
    return (2.0 * sum_min_weights(sig1, sig2)) / denominator


def compute_pair_metrics(
    sample_sig: SignatureData,
    genome_sig: SignatureData,
) -> dict[str, float | int]:
    """
    Compute all final metrics for one harmonized sample/genome pair.

    At this point both signatures must already:
    * refer to the same k-mer size;
    * be metadata-compatible; and
    * share the same effective ``scaled`` value.
    """
    intersection_hash_count = count_intersection_hashes(sample_sig, genome_sig)
    union_hash_count = count_union_hashes(sample_sig, genome_sig)

    sample_containment_in_genome = (
        intersection_hash_count / sample_sig.hash_count if sample_sig.hash_count else 0.0
    )
    genome_containment_in_sample = (
        intersection_hash_count / genome_sig.hash_count if genome_sig.hash_count else 0.0
    )
    max_containment = max(sample_containment_in_genome, genome_containment_in_sample)

    jaccard = (
        intersection_hash_count / union_hash_count if union_hash_count else 0.0
    )
    cosine = cosine_similarity(sample_sig, genome_sig)
    sorensen_dice = sorensen_dice_similarity(sample_sig, genome_sig)
    bray_curtis_sim = bray_curtis_similarity(sample_sig, genome_sig)

    return {
        "intersection_hash_count": intersection_hash_count,
        "union_hash_count": union_hash_count,
        "sample_containment_in_genome": sample_containment_in_genome,
        "genome_containment_in_sample": genome_containment_in_sample,
        "max_containment": max_containment,
        "jaccard_similarity": jaccard,
        "jaccard_distance": 1.0 - jaccard,
        "cosine_similarity": cosine,
        "cosine_distance": 1.0 - cosine,
        "sorensen_dice_similarity": sorensen_dice,
        "sorensen_dice_distance": 1.0 - sorensen_dice,
        "bray_curtis_similarity": bray_curtis_sim,
        "bray_curtis_distance": 1.0 - bray_curtis_sim,
    }


def collect_files(input_dir: Path, pattern: str) -> list[Path]:
    """Resolve one or more comma-separated glob patterns inside ``input_dir``."""
    patterns = [item.strip() for item in pattern.split(",") if item.strip()]
    matched_paths = set()
    for single_pattern in patterns:
        matched_paths.update(
            path for path in input_dir.glob(single_pattern) if path.is_file()
        )
    return sorted(matched_paths)


def build_genome_index(genome_files: list[Path]) -> dict[tuple[int, str], Path]:
    """
    Map each ``(k, software)`` pair to its unique genome sketch file.

    This is what guarantees that a sample generated by one software is only
    compared to a genome sketch generated by that same software.
    """
    genome_by_k = {}
    for genome_file in genome_files:
        ksize = extract_ksize_from_name(genome_file)
        software = extract_software_from_name(genome_file)
        key = (ksize, software)
        if key in genome_by_k:
            raise ValueError(
                f"More than one genome sketch found for k={ksize}, software={software}: "
                f"{genome_by_k[key].name} and {genome_file.name}"
            )
        genome_by_k[key] = genome_file
    return genome_by_k


def group_sample_files_by_key(sample_files: list[Path]) -> dict[tuple[int, str], list[Path]]:
    """
    Group sample files by ``(k, software)``.

    This ensures that the comparison stage keeps software families separated.
    """
    grouped = {}
    for sample_file in sample_files:
        key = (
            extract_ksize_from_name(sample_file),
            extract_software_from_name(sample_file),
        )
        grouped.setdefault(key, []).append(sample_file)
    return grouped


def choose_common_scaled(
    sample_signatures: list[SignatureData],
    genome_signature: SignatureData,
    mode: str,
) -> int | None:
    """
    Choose the scaled value used for harmonization.

    Modes
    -----
    per-k-max
        Use the maximum scaled value seen across the genome plus all samples for
        that specific k-mer size. This is the closest match to a single batch
        search workflow such as branchwater.
    pairwise-max
        Use the maximum scaled value separately for each individual
        sample/genome pair.
    """
    if mode != "per-k-max":
        return None
    return max([genome_signature.scaled] + [sig.scaled for sig in sample_signatures])


def format_float(value: float) -> str:
    """Format floating-point output with a stable precision for CSV export."""
    return f"{value:.10f}"


def write_output(rows: list[dict], output_path: Path) -> None:
    """Write the final comparison table to CSV with an explicit column order."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "ksize",
        "software",
        "sample_id",
        "sample_file",
        "genome_id",
        "genome_file",
        "sample_scaled_original",
        "genome_scaled_original",
        "common_scaled",
        "sample_hash_count_original",
        "genome_hash_count_original",
        "sample_hash_count",
        "genome_hash_count",
        "intersection_hash_count",
        "union_hash_count",
        "sample_containment_in_genome",
        "genome_containment_in_sample",
        "max_containment",
        "uses_abundance",
        "cosine_similarity",
        "cosine_distance",
        "jaccard_similarity",
        "jaccard_distance",
        "sorensen_dice_similarity",
        "sorensen_dice_distance",
        "bray_curtis_similarity",
        "bray_curtis_distance",
    ]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    """
    Orchestrate the full workflow.

    Processing sequence
    -------------------
    1. Read command-line arguments.
    2. Locate sample and genome sketch files.
    3. Group sample sketches by k-mer size and software label.
    4. Load the genome and all samples for each ``(k, software)`` group.
    5. Validate metadata compatibility.
    6. Harmonize scaled values.
    7. Compute metrics.
    8. Write one CSV row per comparison.
    """
    args = parse_args()
    input_dir = args.input_dir

    if not input_dir.exists():
        print(f"Input directory does not exist: {input_dir}", file=sys.stderr)
        return 1

    sample_files = collect_files(input_dir, args.sample_pattern)
    genome_files = collect_files(input_dir, args.genome_pattern)

    if not sample_files:
        print(
            f"No sample sketches found with pattern '{args.sample_pattern}' in {input_dir}",
            file=sys.stderr,
        )
        return 1

    if not genome_files:
        print(
            f"No genome sketches found with pattern '{args.genome_pattern}' in {input_dir}",
            file=sys.stderr,
        )
        return 1

    genome_by_k = build_genome_index(genome_files)
    sample_files_by_key = group_sample_files_by_key(sample_files)

    rows = []
    skipped = []

    for ksize, software in sorted(sample_files_by_key):
        genome_file = genome_by_k.get((ksize, software))
        if genome_file is None:
            skipped.extend(path.name for path in sample_files_by_key[(ksize, software)])
            continue

        genome_signature = load_signature(genome_file, ksize)
        sample_signatures = [
            load_signature(sample_file, ksize)
            for sample_file in sample_files_by_key[(ksize, software)]
        ]

        for sample_signature in sample_signatures:
            ensure_compatible(sample_signature, genome_signature)

        group_common_scaled = choose_common_scaled(
            sample_signatures,
            genome_signature,
            args.common_scaled_mode,
        )
        group_genome_signature = (
            genome_signature.downsample(group_common_scaled)
            if group_common_scaled is not None
            else None
        )

        for sample_signature in sample_signatures:
            common_scaled = (
                max(sample_signature.scaled, genome_signature.scaled)
                if args.common_scaled_mode == "pairwise-max"
                else group_common_scaled
            )
            if common_scaled is None:
                raise RuntimeError("common_scaled unexpectedly missing")

            compared_sample_signature = sample_signature.downsample(common_scaled)
            compared_genome_signature = (
                group_genome_signature
                if group_genome_signature is not None
                else genome_signature.downsample(common_scaled)
            )

            metrics = compute_pair_metrics(
                compared_sample_signature,
                compared_genome_signature,
            )

            rows.append(
                {
                    "ksize": ksize,
                    "software": software,
                    "sample_id": sample_signature.label or strip_k_suffix(sample_signature.path),
                    "sample_file": str(sample_signature.path),
                    "genome_id": genome_signature.label or strip_k_suffix(genome_signature.path),
                    "genome_file": str(genome_signature.path),
                    "sample_scaled_original": sample_signature.scaled,
                    "genome_scaled_original": genome_signature.scaled,
                    "common_scaled": common_scaled,
                    "sample_hash_count_original": sample_signature.hash_count,
                    "genome_hash_count_original": genome_signature.hash_count,
                    "sample_hash_count": compared_sample_signature.hash_count,
                    "genome_hash_count": compared_genome_signature.hash_count,
                    "intersection_hash_count": metrics["intersection_hash_count"],
                    "union_hash_count": metrics["union_hash_count"],
                    "sample_containment_in_genome": format_float(
                        metrics["sample_containment_in_genome"]
                    ),
                    "genome_containment_in_sample": format_float(
                        metrics["genome_containment_in_sample"]
                    ),
                    "max_containment": format_float(metrics["max_containment"]),
                    "uses_abundance": sample_signature.uses_abundance
                    or genome_signature.uses_abundance,
                    "cosine_similarity": format_float(metrics["cosine_similarity"]),
                    "cosine_distance": format_float(metrics["cosine_distance"]),
                    "jaccard_similarity": format_float(metrics["jaccard_similarity"]),
                    "jaccard_distance": format_float(metrics["jaccard_distance"]),
                    "sorensen_dice_similarity": format_float(
                        metrics["sorensen_dice_similarity"]
                    ),
                    "sorensen_dice_distance": format_float(
                        metrics["sorensen_dice_distance"]
                    ),
                    "bray_curtis_similarity": format_float(
                        metrics["bray_curtis_similarity"]
                    ),
                    "bray_curtis_distance": format_float(
                        metrics["bray_curtis_distance"]
                    ),
                }
            )

    rows.sort(key=lambda row: (int(row["ksize"]), row["sample_id"]))
    write_output(rows, args.output)

    print(f"Wrote {len(rows)} comparisons to {args.output}")
    if skipped:
        print(
            "Skipped sample sketches without a matching genome k: "
            + ", ".join(skipped),
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
