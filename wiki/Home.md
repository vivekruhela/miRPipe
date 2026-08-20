# miRPipe user guide

This guide covers both supported ways to run miRPipe:

| Workflow | Best for | Interface | Execution model |
|---|---|---|---|
| [Docker notebook](Docker-Workflow.md) | First runs, small experiments, and comparison with the published notebook | Jupyter notebook | Interactive, one host |
| [HPC workflow](HPC-Workflow.md) | Cohorts, shared clusters, unattended runs, and reproducible reruns | Nextflow configuration | Parallel, resumable SLURM jobs |

Both workflows analyze single-end small-RNA sequencing data. The Docker notebook is the legacy implementation used for the original miRPipe analyses. The HPC workflow is a non-interactive Nextflow DSL2 implementation with input validation, per-sample parallelism, containers, execution reports, and `-resume`.

> **Do not compare or combine results across incompatible references.** The genome assembly, miRNA catalogue, piRNA annotation, chromosome names, and aligner indexes must describe the same reference build.

## Start here

### Docker notebook in five commands

```bash
git clone https://github.com/vivekruhela/miRPipe.git
cd miRPipe
mkdir -p data
# Copy FASTQ/FASTQ.GZ files and sample_list.csv into data/.
docker pull docker.io/vivekruhela/mirpipe:latest
docker run --rm --name mirpipe \
  -p 127.0.0.1:8880:8888 \
  -e PASSWORD='choose-a-password' \
  -e USE_HTTP=1 \
  -v "$(pwd)/data:/miRPipe/data" \
  docker.io/vivekruhela/mirpipe:latest
```

Open <http://localhost:8880/mirpipe>, enter the password, and run `mirpipe_pipeline.ipynb` from top to bottom. See the [complete Docker procedure](Docker-Workflow.md) before using research data.

### HPC/SLURM in five commands

```bash
git clone https://github.com/vivekruhela/miRPipe.git
cd miRPipe/hpc_workflow
cp assets/samplesheet.example.csv samplesheet.csv
cp assets/params.example.yml params.yml
# Edit samplesheet.csv and params.yml, then submit from this directory.
sbatch submit_mirpipe.slurm
```

The submission wrapper requires Nextflow plus Apptainer/Singularity on `PATH`. See the [complete HPC procedure](HPC-Workflow.md) for reference configuration, site-specific SLURM options, monitoring, and recovery.

## Example and supporting pages

- [Worked four-sample example](Worked-Example.md): the same experiment expressed in both input formats.
- [Troubleshooting](Troubleshooting.md): common Docker, Jupyter, Nextflow, SLURM, container, and data errors.
- [Reproducibility checklist](Reproducibility-Checklist.md): what to pin and archive for a defensible analysis.
- [miRPipe2 roadmap](../hpc_workflow/docs/MIRPIPE2_ROADMAP.md): proposed accuracy, speed, automation, and biological-evidence improvements.

## Which workflow should I report?

Record the exact workflow in the Methods section. “miRPipe” alone is ambiguous because the implementations differ in orchestration, length handling, annotation, and provenance.

- For the notebook, report the repository commit, Docker image digest, notebook choices, adapter, reference versions, and sample-sheet order.
- For HPC, report the repository commit, Nextflow version, container digest or SIF checksum, `params.yml`, validated sample sheet, reference checksums, and files in `results/pipeline_info/`.

## Citation

If miRPipe contributes to a publication, cite:

Ruhela V, Gupta A, Krishnamachari S, Ahuja G, Kaur G, Gupta R. *miRPipe: A Unified Computational Framework for a Robust, Reliable, and Reproducible Identification of Novel miRNAs from the RNA Sequencing Data.* Frontiers in Bioinformatics (2022). <https://doi.org/10.3389/fbinf.2022.842051>
