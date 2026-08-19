#!/usr/bin/env python3
"""Validate and normalize a miRPipe-HPC single-end sample sheet."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path


REQUIRED_COLUMNS = ("sample", "fastq", "condition")
SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def fail(message: str) -> None:
    raise ValueError(message)


def validate(input_path: Path) -> tuple[list[str], list[dict[str, str]], dict[str, object]]:
    input_path = input_path.resolve()
    if not input_path.is_file():
        fail(f"Sample sheet does not exist: {input_path}")

    with input_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            fail("Sample sheet is empty or has no header.")
        fieldnames = [name.strip() for name in reader.fieldnames]
        missing = [name for name in REQUIRED_COLUMNS if name not in fieldnames]
        if missing:
            fail(f"Missing required sample-sheet columns: {', '.join(missing)}")

        rows: list[dict[str, str]] = []
        for line_number, raw in enumerate(reader, start=2):
            row = {(key or "").strip(): (value or "").strip() for key, value in raw.items()}
            if not any(row.values()):
                continue
            sample = row.get("sample", "")
            condition = row.get("condition", "")
            fastq_value = row.get("fastq", "")
            if not sample:
                fail(f"Line {line_number}: sample is empty.")
            if not SAMPLE_RE.fullmatch(sample):
                fail(
                    f"Line {line_number}: invalid sample ID {sample!r}; use letters, digits, '.', '_', or '-'."
                )
            if not condition:
                fail(f"Line {line_number}: condition is empty for sample {sample!r}.")
            if not fastq_value:
                fail(f"Line {line_number}: fastq is empty for sample {sample!r}.")

            fastq = Path(fastq_value).expanduser()
            if not fastq.is_absolute():
                fastq = input_path.parent / fastq
            fastq = fastq.resolve()
            if not fastq.is_file():
                fail(f"Line {line_number}: FASTQ does not exist: {fastq}")
            if not str(fastq).lower().endswith(FASTQ_SUFFIXES):
                fail(f"Line {line_number}: unsupported FASTQ suffix: {fastq.name}")
            row["fastq"] = str(fastq)
            rows.append(row)

    if not rows:
        fail("Sample sheet contains no data rows.")

    samples = [row["sample"] for row in rows]
    duplicate_samples = sorted(name for name, n in Counter(samples).items() if n > 1)
    if duplicate_samples:
        fail(f"Duplicate sample IDs: {', '.join(duplicate_samples)}")

    fastqs = [row["fastq"] for row in rows]
    duplicate_fastqs = sorted(name for name, n in Counter(fastqs).items() if n > 1)
    if duplicate_fastqs:
        fail("Each FASTQ must occur once; duplicates: " + ", ".join(duplicate_fastqs))

    conditions = Counter(row["condition"] for row in rows)
    batches = Counter(row.get("batch", "") or "not_provided" for row in rows)
    summary = {
        "samples": len(rows),
        "conditions": dict(sorted(conditions.items())),
        "batches": dict(sorted(batches.items())),
        "columns": fieldnames,
        "source": str(input_path),
    }
    return fieldnames, rows, summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    try:
        fieldnames, rows, summary = validate(args.input)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
