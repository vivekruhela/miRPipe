#!/usr/bin/env python3
"""Attach the best Infernal/Rfam hit to each candidate annotation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_hits(path: Path) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(maxsplit=17)
            if len(fields) < 17:
                continue
            target, accession, query = fields[0], fields[1], fields[2]
            score, evalue = fields[14], fields[15]
            current = best.get(query)
            if current is None or float(evalue) < float(current["rfam_evalue"]):
                best[query] = {
                    "rfam_name": target,
                    "rfam_accession": accession,
                    "rfam_score": score,
                    "rfam_evalue": evalue,
                }
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--tblout", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    hits = read_hits(args.tblout)
    with args.annotation.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []

    extra = ["rfam_name", "rfam_accession", "rfam_score", "rfam_evalue", "rfam_interpretation"]
    for row in rows:
        hit = hits.get(row["feature_id"], {})
        row.update({field: hit.get(field, "") for field in extra[:-1]})
        name = hit.get("rfam_name", "").lower()
        if not name:
            interpretation = "no_significant_hit"
        elif "mir" in name or "let-7" in name:
            interpretation = "miRNA_family_homology"
        else:
            interpretation = "other_ncRNA_hit_review_or_reject"
            flags = filter(None, [row.get("review_flags", ""), "rfam_other_ncRNA_hit"])
            row["review_flags"] = ";".join(flags)
            if row.get("candidate_class") == "novel_miRNA_candidate":
                row["confidence_tier"] = "rejected_other_ncRNA"
        row["rfam_interpretation"] = interpretation

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields + extra, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
