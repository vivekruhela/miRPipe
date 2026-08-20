# Reproducibility checklist

Archive this checklist with each analysis. A checked box should correspond to a retained file, checksum, or Methods statement.

## Inputs and experimental design

- [ ] Raw FASTQ files are immutable and have recorded SHA-256 checksums.
- [ ] Sample IDs are unique and map to biological samples, not ambiguous lane names.
- [ ] The adapter comes from the library protocol or adapter-detection QC.
- [ ] Biological conditions, batches, pairs, and exclusions are documented before testing.
- [ ] The sample sheet used by the run is archived exactly as executed.
- [ ] Controlled data and credentials are excluded from Git and public issue reports.

## Software and workflow

- [ ] The miRPipe Git commit is recorded with `git rev-parse HEAD`.
- [ ] The Docker image digest or Apptainer SIF SHA-256 is recorded; `latest` alone is not used as an archival identifier.
- [ ] Docker, Nextflow, Java, Apptainer/Singularity, and scheduler versions are recorded as applicable.
- [ ] The exact launch command and working directory are recorded.
- [ ] Local workflow/configuration changes are committed or saved as a patch.
- [ ] The executed notebook or `params.yml` is archived.

## References

- [ ] Genome FASTA/build, chromosome naming scheme, and checksum are recorded.
- [ ] Bowtie index was built from the same genome and its version/checksums are recorded.
- [ ] miRNA FASTA/GFF database name, release, download source/date, and checksums are recorded.
- [ ] piRNA FASTA/GFF database name, release, download source/date, and checksums are recorded.
- [ ] Optional Rfam/external database releases and index checksums are recorded.
- [ ] No hg19/hg38 or incompatible chromosome-name mixture exists.

## Parameters and quality control

- [ ] Adapter, trimming quality, read-length intervals, alignment settings, strandedness, and random seeds are recorded.
- [ ] Enabled branches and all non-default parameters are reported.
- [ ] MultiQC/FastQC and trimming logs were reviewed before downstream interpretation.
- [ ] Every expected sample appears once in the merged count matrices.
- [ ] Count filtering, design formula, contrast direction, alpha, and multiple-testing method are reported.
- [ ] Exclusions or reruns are justified and traceable.

## Outputs and interpretation

- [ ] Raw counts are retained separately from normalized values and significant subsets.
- [ ] HPC `pipeline_info/software_versions.txt`, trace, report, timeline, and DAG are archived.
- [ ] Docker notebook choices, sample ordering, image digest, and core output files are archived.
- [ ] Novel miRNA IDs are described as provisional unless accepted by an authoritative database.
- [ ] piRNAdb overlap is described as annotation, not automatically as active piRNA biogenesis.
- [ ] Comparisons between Docker and HPC hold references and biological thresholds constant and disclose implementation differences.
- [ ] A small truth-set or positive/negative control was included when accuracy claims are made.

For the rationale behind stronger biological validation, see the [miRPipe2 roadmap](../hpc_workflow/docs/MIRPIPE2_ROADMAP.md).
