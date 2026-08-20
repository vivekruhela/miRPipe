#!/usr/bin/env python3
"""Perform deterministic local sequence, reverse-complement, locus, and seed annotation."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import shutil
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def normalize_sequence(sequence: str) -> str:
    return "".join(sequence.upper().replace("U", "T").split())


def reverse_complement(sequence: str) -> str:
    return normalize_sequence(sequence).translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else path.open(encoding="utf-8", errors="replace")


def read_fasta(path: Path | None) -> tuple[dict[str, str], dict[str, list[str]]]:
    names: dict[str, str] = {}
    by_sequence: dict[str, list[str]] = defaultdict(list)
    if path is None:
        return names, by_sequence
    current = ""
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequence = normalize_sequence("".join(chunks))
                    names[current] = sequence
                    by_sequence[sequence].append(current)
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current:
        sequence = normalize_sequence("".join(chunks))
        names[current] = sequence
        by_sequence[sequence].append(current)
    return names, by_sequence


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
        else:
            continue
        result[key.strip()] = unquote(value.strip().strip('"'))
    return result


def read_gff(path: Path | None) -> dict[str, list[tuple[str, int, int, str]]]:
    loci: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
    if path is None:
        return loci
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = parse_attributes(fields[8])
            candidates = [attrs.get(key, "") for key in ("Name", "ID", "Alias", "Derives_from")]
            for candidate in candidates:
                for name in filter(None, candidate.split(",")):
                    loci[name].append((fields[0], int(fields[3]), int(fields[4]), fields[6]))
    return loci


def coordinate_status(row: dict[str, str], names: list[str], loci: dict[str, list[tuple[str, int, int, str]]], tolerance: int) -> str:
    available = []
    for name in names:
        available.extend(loci.get(name, []))
    if not available:
        return "not_available"
    try:
        start, end = int(row.get("start", "")), int(row.get("end", ""))
    except ValueError:
        return "not_available"
    chrom = row.get("chrom", "")
    strand = row.get("strand", ".")
    for ref_chrom, ref_start, ref_end, ref_strand in available:
        strand_ok = strand in ("", ".") or ref_strand in ("", ".") or strand == ref_strand
        if chrom == ref_chrom and strand_ok and abs(start - ref_start) <= tolerance and abs(end - ref_end) <= tolerance:
            return "yes"
    return "no"


def numeric(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "0"))
    except ValueError:
        return 0.0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--known-mature-fasta", required=True, type=Path)
    parser.add_argument("--known-gff", required=True, type=Path)
    parser.add_argument("--external-mature-fasta", type=Path)
    parser.add_argument("--external-gff", type=Path)
    parser.add_argument("--locus-tolerance", type=int, default=2)
    parser.add_argument("--min-samples", type=int, default=2)
    parser.add_argument("--min-total-count", type=int, default=10)
    parser.add_argument("--min-caller-score", type=float, default=5.0)
    parser.add_argument("--output-counts", required=True, type=Path)
    parser.add_argument("--output-annotation", required=True, type=Path)
    parser.add_argument("--output-fasta", required=True, type=Path)
    args = parser.parse_args()

    known_names, known_by_sequence = read_fasta(args.known_mature_fasta)
    external_names, external_by_sequence = read_fasta(args.external_mature_fasta)
    known_loci = read_gff(args.known_gff)
    external_loci = read_gff(args.external_gff)

    known_seed: dict[str, list[str]] = defaultdict(list)
    for name, sequence in known_names.items():
        if len(sequence) >= 8:
            known_seed[sequence[1:8]].append(name)

    with args.candidates.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        input_fields = reader.fieldnames or []

    extra_fields = [
        "canonical_seed_2_8",
        "seed_cluster",
        "reference_match",
        "reference_match_type",
        "coordinate_concordance",
        "known_seed_family",
        "external_match",
        "candidate_class",
        "evidence_score",
        "confidence_tier",
        "review_flags",
    ]

    fasta_records: list[tuple[str, str]] = []
    for row in rows:
        sequence = normalize_sequence(row.get("mature_sequence", ""))
        precursor = normalize_sequence(row.get("precursor_sequence", ""))
        seed = sequence[1:8] if len(sequence) >= 8 else ""
        exact_names = sorted(known_by_sequence.get(sequence, []))
        reverse_names = sorted(known_by_sequence.get(reverse_complement(sequence), [])) if sequence else []
        external_exact = sorted(external_by_sequence.get(sequence, []))
        external_reverse = sorted(external_by_sequence.get(reverse_complement(sequence), [])) if sequence else []

        match_names = exact_names or reverse_names
        match_type = "exact" if exact_names else ("reverse_complement" if reverse_names else "none")
        coord = coordinate_status(row, match_names, known_loci, args.locus_tolerance) if match_names else "not_applicable"
        ext_names = external_exact or external_reverse
        ext_coord = coordinate_status(row, ext_names, external_loci, args.locus_tolerance) if ext_names else "not_applicable"

        flags: list[str] = []
        if match_names and coord in ("yes", "not_available"):
            candidate_class = "known_miRNA_reference" if match_type == "exact" else "known_miRNA_reverse_complement"
            confidence = "reference"
            evidence_score = 5
        elif match_names:
            candidate_class = "known_sequence_paralogue_candidate"
            confidence = "medium"
            evidence_score = 3
            flags.append("known_sequence_coordinate_mismatch")
        elif ext_names and ext_coord in ("yes", "not_available"):
            candidate_class = "external_database_match"
            confidence = "medium"
            evidence_score = 3
        else:
            candidate_class = "novel_miRNA_candidate"
            evidence_score = 0
            evidence_score += numeric(row, "max_caller_score") >= args.min_caller_score
            evidence_score += numeric(row, "total_star_reads") > 0
            evidence_score += numeric(row, "n_samples_detected") >= args.min_samples
            evidence_score += numeric(row, "total_count") >= args.min_total_count
            confidence = "high" if evidence_score == 4 else ("medium" if evidence_score >= 2 else "low")
            if numeric(row, "total_star_reads") == 0:
                flags.append("no_star_read_evidence")
            if numeric(row, "n_samples_detected") < args.min_samples:
                flags.append("limited_sample_recurrence")
            if not precursor:
                flags.append("missing_precursor_sequence")

        row.update(
            {
                "canonical_seed_2_8": seed,
                "seed_cluster": f"seed_{hashlib.sha1(seed.encode()).hexdigest()[:10]}" if seed else "",
                "reference_match": ";".join(match_names),
                "reference_match_type": match_type,
                "coordinate_concordance": coord,
                "known_seed_family": ";".join(sorted(known_seed.get(seed, []))),
                "external_match": ";".join(ext_names),
                "candidate_class": candidate_class,
                "evidence_score": int(evidence_score),
                "confidence_tier": confidence,
                "review_flags": ";".join(flags),
            }
        )
        if candidate_class in ("novel_miRNA_candidate", "known_sequence_paralogue_candidate") and precursor:
            fasta_records.append((row["feature_id"], precursor))

    shutil.copyfile(args.counts, args.output_counts)
    with args.output_annotation.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=input_fields + extra_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with args.output_fasta.open("w", encoding="utf-8") as handle:
        for feature_id, precursor in fasta_records:
            handle.write(f">{feature_id}\n{precursor}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
