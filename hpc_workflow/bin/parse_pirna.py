#!/usr/bin/env python3
"""Convert BEDTools GFF coverage output into miRPipe's generic count schema."""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


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


def attributes(text: str) -> dict[str, str]:
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    merged: dict[str, dict[str, object]] = {}
    with args.input.open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                raise ValueError(f"Expected GFF3 plus count at {args.input}:{line_number}")
            chrom, _, feature_type, start, end, _, strand, _, attr_text = fields[:9]
            count = int(round(float(fields[-1])))
            attrs = attributes(attr_text)
            name = next(
                (attrs[key] for key in ("Name", "ID", "Alias", "piRNA", "gene_id") if attrs.get(key)),
                "",
            )
            if not name:
                digest = hashlib.sha1(line.encode("utf-8")).hexdigest()[:12]
                name = f"piRNA_{chrom}_{start}_{end}_{digest}"
            feature_id = f"piRNA:{name}|{chrom}:{start}-{end}:{strand}"
            if feature_id not in merged:
                merged[feature_id] = {
                    "feature_id": feature_id,
                    "name": name,
                    "biotype": feature_type or "piRNA",
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "mature_sequence": "",
                    "precursor_sequence": "",
                    "caller_score": "",
                    "star_read_count": 0,
                    "count": 0,
                    "sample": args.sample,
                }
            merged[feature_id]["count"] = int(merged[feature_id]["count"]) + count

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(merged):
            writer.writerow(merged[key])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
