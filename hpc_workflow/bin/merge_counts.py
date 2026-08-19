#!/usr/bin/env python3
"""Merge per-sample miRNA or piRNA tables without relying on column order."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


BASE_ANNOTATION_FIELDS = [
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
]


def numeric(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True, type=Path)
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    args = parser.parse_args()

    with args.samplesheet.open(newline="", encoding="utf-8-sig") as handle:
        sample_rows = list(csv.DictReader(handle))
    samples = [row["sample"] for row in sample_rows]
    sample_set = set(samples)

    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    annotations: dict[str, dict[str, str]] = {}
    max_score: dict[str, float] = defaultdict(lambda: float("-inf"))
    total_star: dict[str, int] = defaultdict(int)

    for input_path in sorted(args.inputs):
        with input_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or not {"feature_id", "sample", "count"}.issubset(reader.fieldnames):
                raise ValueError(f"Invalid count table: {input_path}")
            for row in reader:
                feature_id = row["feature_id"]
                sample = row["sample"]
                if sample not in sample_set:
                    raise ValueError(f"Sample {sample!r} in {input_path} is absent from the sample sheet.")
                counts[feature_id][sample] += int(round(numeric(row.get("count", "0"))))
                max_score[feature_id] = max(max_score[feature_id], numeric(row.get("caller_score", "")))
                total_star[feature_id] += int(round(numeric(row.get("star_read_count", "0"))))
                if feature_id not in annotations:
                    annotations[feature_id] = {field: row.get(field, "") for field in BASE_ANNOTATION_FIELDS}
                else:
                    for field in BASE_ANNOTATION_FIELDS:
                        incoming = row.get(field, "")
                        if not annotations[feature_id].get(field) and incoming:
                            annotations[feature_id][field] = incoming

    if not counts:
        raise ValueError("No feature counts were produced.")

    def sort_key(feature_id: str) -> tuple[str, int, str]:
        annotation = annotations[feature_id]
        try:
            start = int(annotation.get("start", ""))
        except ValueError:
            start = 0
        return annotation.get("chrom", ""), start, feature_id

    features = sorted(counts, key=sort_key)
    with args.counts.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature_id", *samples])
        for feature_id in features:
            writer.writerow([feature_id, *(counts[feature_id].get(sample, 0) for sample in samples)])

    annotation_fields = BASE_ANNOTATION_FIELDS + ["n_samples_detected", "total_count", "max_caller_score", "total_star_reads"]
    with args.annotation.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=annotation_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for feature_id in features:
            row = dict(annotations[feature_id])
            row["n_samples_detected"] = sum(counts[feature_id].get(sample, 0) > 0 for sample in samples)
            row["total_count"] = sum(counts[feature_id].values())
            row["max_caller_score"] = "" if max_score[feature_id] == float("-inf") else f"{max_score[feature_id]:.6g}"
            row["total_star_reads"] = total_star[feature_id]
            writer.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
