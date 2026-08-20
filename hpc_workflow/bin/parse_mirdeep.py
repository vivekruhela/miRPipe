#!/usr/bin/env python3
"""Convert miRDeep* result tables into a stable, per-sample feature table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path


OUTPUT_FIELDS = [
    "feature_id",
    "name",
    "biotype",
    "chrom",
    "start",
    "end",
    "strand",
    "mature_sequence",
    "precursor_sequence",
    "caller_score",
    "star_read_count",
    "count",
    "sample",
]


def normalized(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (value or "").lower())


def pick(row: dict[str, str], exact: tuple[str, ...], contains: tuple[str, ...] = ()) -> str:
    keyed = {normalized(key): value.strip() for key, value in row.items() if key is not None}
    for key in exact:
        if normalized(key) in keyed:
            return keyed[normalized(key)]
    for key, value in keyed.items():
        if all(token in key for token in contains):
            return value
    return ""


def as_number(value: str, default: float = 0.0) -> float:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def coordinates(value: str) -> tuple[str, str]:
    numbers = re.findall(r"\d+", value or "")
    return (numbers[0], numbers[1]) if len(numbers) >= 2 else ("", "")


def feature_key(name: str, chrom: str, start: str, end: str, strand: str, sequence: str, known: bool) -> str:
    locus = f"{chrom or 'NA'}:{start or 'NA'}-{end or 'NA'}:{strand or '.'}"
    if known:
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
        return f"known:{safe_name}|{locus}"
    digest = hashlib.sha1(sequence.encode("ascii", errors="ignore")).hexdigest()[:12]
    return f"novel:{locus}:{digest}"


def parse_file(path: Path, sample: str, source_known: bool) -> list[dict[str, object]]:
    parsed: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header in miRDeep* result: {path}")
        first_column = reader.fieldnames[0]
        for line_number, row in enumerate(reader, start=2):
            name = (row.get(first_column) or "").strip()
            if not name:
                continue
            chrom = pick(row, ("chr", "chromosome"))
            locus = pick(row, ("mature_loci", "mature loci"))
            start, end = coordinates(locus)
            strand = pick(row, ("strand",)) or "."
            mature = pick(row, ("mature miR", "mature_miR", "mature sequence"))
            precursor = pick(row, ("sequence", "precursor sequence"))
            count_text = pick(
                row,
                ("expression(number of mature reads)", "expression", "mature reads"),
                contains=("expression", "mature", "reads"),
            )
            if not count_text:
                raise ValueError(f"Cannot identify mature-read count column in {path}:{line_number}")
            score_text = pick(row, ("score", "miRDeep score"), contains=("score",))
            star_text = pick(row, ("star reads", "star read count"), contains=("star", "read"))
            known = source_known or name.lower().startswith(("hsa-mir", "hsa-let"))
            mature = mature.upper().replace("U", "T")
            parsed.append(
                {
                    "feature_id": feature_key(name, chrom, start, end, strand, mature, known),
                    "name": name,
                    "biotype": "known_miRNA" if known else "novel_miRNA_candidate",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "mature_sequence": mature,
                    "precursor_sequence": precursor.upper().replace("U", "T"),
                    "caller_score": score_text,
                    "star_read_count": int(round(as_number(star_text))),
                    "count": int(round(as_number(count_text))),
                    "sample": sample,
                }
            )
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    files = sorted(args.results_dir.rglob("*.result"))
    files = [path for path in files if path.is_file()]
    if not files:
        print("ERROR: miRDeep* produced no .result files", file=sys.stderr)
        return 2

    merged: dict[str, dict[str, object]] = {}
    try:
        for path in files:
            source_known = ".known_miR.result" in path.name
            for row in parse_file(path, args.sample, source_known):
                key = str(row["feature_id"])
                if key not in merged:
                    merged[key] = row
                else:
                    # miRDeep* can repeat a known entry in both its general and
                    # known-miRNA tables. Retain one quantified call rather than
                    # double counting the same sample/locus.
                    merged[key]["count"] = max(int(merged[key]["count"]), int(row["count"]))
                    merged[key]["star_read_count"] = max(
                        int(merged[key]["star_read_count"]), int(row["star_read_count"])
                    )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(merged):
            writer.writerow(merged[key])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
