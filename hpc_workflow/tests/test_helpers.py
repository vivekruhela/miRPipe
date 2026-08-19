from __future__ import annotations

import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"


def run_script(name: str, *args: object) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["python3", str(BIN / name), *(str(arg) for arg in args)],
        check=True,
        text=True,
        capture_output=True,
    )


class HelperTests(unittest.TestCase):
    def test_samplesheet_validation_resolves_relative_fastq(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "sample.fastq.gz").write_bytes(b"")
            (root / "samples.csv").write_text(
                "sample,fastq,condition,batch\nS1,sample.fastq.gz,case,b1\n", encoding="utf-8"
            )
            output = root / "validated.csv"
            summary = root / "summary.json"
            run_script(
                "check_samplesheet.py",
                "--input", root / "samples.csv",
                "--output", output,
                "--summary", summary,
            )
            with output.open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(row["sample"], "S1")
            self.assertEqual(row["fastq"], str((root / "sample.fastq.gz").resolve()))

    def test_parse_and_merge_mirdeep_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            header = (
                "pre-miR name\tchr\tmature_loci\tstrand\tmature miR\t"
                "expression(number of mature reads)\tscore\tstar reads\tfoo\tbar\tsequence\n"
            )
            row = "hsa-miR-test\tchr1\t100-121\t+\tAUGCAUGCA\t12\t8.5\t2\tx\ty\tAUGCAUGCAAAAA\n"
            # miRDeep* may repeat a known call in both result tables; it must
            # not be counted twice.
            (root / "S1.result").write_text(header + row, encoding="utf-8")
            (root / "S1.known_miR.result").write_text(header + row, encoding="utf-8")
            parsed = root / "S1.tsv"
            run_script(
                "parse_mirdeep.py",
                "--sample", "S1",
                "--results-dir", root,
                "--output", parsed,
            )
            (root / "samples.csv").write_text(
                "sample,fastq,condition\nS1,/tmp/S1.fastq.gz,case\n", encoding="utf-8"
            )
            counts = root / "counts.tsv"
            annotation = root / "annotation.tsv"
            run_script(
                "merge_counts.py",
                "--inputs", parsed,
                "--samplesheet", root / "samples.csv",
                "--counts", counts,
                "--annotation", annotation,
            )
            self.assertIn("\t12\n", counts.read_text(encoding="utf-8"))
            self.assertIn("known_miRNA", annotation.read_text(encoding="utf-8"))

    def test_local_reverse_complement_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "counts.tsv").write_text("feature_id\tS1\nf1\t10\n", encoding="utf-8")
            fields = [
                "feature_id", "name", "biotype", "chrom", "start", "end", "strand",
                "mature_sequence", "precursor_sequence", "caller_score", "star_read_count",
                "n_samples_detected", "total_count", "max_caller_score", "total_star_reads",
            ]
            with (root / "candidates.tsv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerow(
                    {
                        "feature_id": "f1", "name": "novel", "biotype": "novel_miRNA_candidate",
                        "chrom": "chr1", "start": "100", "end": "108", "strand": "-",
                        "mature_sequence": "TGCAT", "precursor_sequence": "TGCATAAAAA",
                        "caller_score": "8", "star_read_count": "1", "n_samples_detected": "2",
                        "total_count": "10", "max_caller_score": "8", "total_star_reads": "1",
                    }
                )
            (root / "known.fa").write_text(">hsa-miR-test\nATGCA\n", encoding="utf-8")
            (root / "known.gff3").write_text(
                "chr1\tmiRBase\tmiRNA\t100\t108\t.\t-\t.\tID=x;Name=hsa-miR-test\n", encoding="utf-8"
            )
            output_annotation = root / "annotated.tsv"
            run_script(
                "annotate_mirna.py",
                "--counts", root / "counts.tsv",
                "--candidates", root / "candidates.tsv",
                "--known-mature-fasta", root / "known.fa",
                "--known-gff", root / "known.gff3",
                "--output-counts", root / "counts.out.tsv",
                "--output-annotation", output_annotation,
                "--output-fasta", root / "novel.fa",
            )
            self.assertIn("known_miRNA_reverse_complement", output_annotation.read_text(encoding="utf-8"))

    def test_parse_pirna_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            coverage = root / "coverage.gff3"
            coverage.write_text(
                "chr1\tpiRNAdb\tpiRNA\t10\t35\t.\t+\t.\tID=piR-test;Name=piR-test\t7\n",
                encoding="utf-8",
            )
            output = root / "pi.tsv"
            run_script("parse_pirna.py", "--sample", "S1", "--input", coverage, "--output", output)
            self.assertIn("\t7\tS1\n", output.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
