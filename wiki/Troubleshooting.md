# Troubleshooting

Start with the first failing stage. Preserve logs and partial outputs until the cause is understood.

## Docker and Jupyter

### Docker cannot connect to the daemon

```bash
docker version
systemctl status docker
```

Start Docker using your operating system's supported method. If the socket is permission denied, use the institution's approved Docker group or rootless setup; do not weaken socket permissions.

### Port 8880 is already in use

Use another host port while keeping the container port at 8888:

```bash
docker run --rm -p 127.0.0.1:8890:8888 ...
```

Then open <http://localhost:8890/mirpipe>.

### Jupyter opens but data is missing

Confirm the host path and container mount:

```bash
pwd
find data -maxdepth 1 -type f -print
docker inspect mirpipe --format '{{json .Mounts}}'
```

The prepared host directory must be mounted at `/miRPipe/data`. On SELinux hosts, use the distribution's documented bind-mount labeling option if required.

### Results are owned by root

Do not recursively change permissions on shared storage without checking policy. A cluster or lab administrator may recommend Docker user mapping or an ACL for the specific data directory.

### Notebook rerun mixes old and new results

The legacy notebook is not cache-aware. Move the previous result directory to a dated archive, confirm which intermediate files remain, and rerun into a clean mounted directory. Never merge outputs produced with different adapters or reference choices.

### Differential expression has mislabelled samples

Check all three together:

- `sample_list.csv` is sorted alphanumerically by `sample`;
- sample IDs and filenames are exact and unique;
- count-matrix columns occur in the same order expected by the notebook.

This ordering constraint applies to the legacy Docker notebook, not the name-matched HPC workflow.

## Nextflow, SLURM, and Apptainer

### `nextflow: command not found` or Java is incompatible

```bash
module avail nextflow java
module load java nextflow
java -version
nextflow -version
```

miRPipe-HPC requires Nextflow 24.04.0 or newer and Java 11 or newer.

### Driver exits because Apptainer/Singularity is missing

Load the site's container module before `sbatch`, and ensure modules are available inside the submitted driver environment. The wrapper accepts either `apptainer` or `singularity` on `PATH`.

### FASTQ path is reported missing

The `fastq` value is absolute or relative to `samplesheet.csv`, not necessarily the shell's current directory. Confirm every path from a compute node and ensure the storage is mounted there.

### SLURM rejects child jobs

Typical causes are an invalid account, partition, constraint, memory, time, or submission quota. Read the scheduler's rejection message, obtain correct values from the cluster administrators, and place them in `cluster.config` rather than hard-coding a site into the shared workflow.

### Container cannot see a reference

Use shared absolute host paths and verify Apptainer auto-mount policy. Paths such as `/miRPipe/...` refer to files bundled inside the legacy-compatible image; host paths must be bind-mounted and visible on compute nodes.

### `MD.jar`, Bowtie index, or annotation files do not match

Check `genome`, `mirdeep_genome_dir`, `known_mature_fasta`, `known_mirna_gff`, `bowtie_index`, and `pirna_gff` as one reference bundle. Mixing hg19 and hg38 can finish computationally while producing invalid biology.

### A task exits 137, 143, or 247

These often indicate memory, termination, or scheduler resource pressure; the SLURM profile retries selected failures up to two times. Inspect `.command.err`, the SLURM reason, and peak RSS in `execution_trace.tsv`. If memory is genuinely insufficient, override the affected process label in `cluster.config` and rerun with `-resume`.

### Resume recomputes more tasks than expected

Nextflow invalidates cached tasks when commands, inputs, parameters, container identity, or relevant workflow code change. Use the same launch directory and work directory, retain `.nextflow/`, and keep `-resume`. Recalculation after a meaningful change is correct behavior.

### Differential-expression stage is empty or unstable

Confirm there are biological replicates in each contrast level, the design columns exist without missing values, the contrast levels are spelled exactly, and features pass `min_count`/`min_samples`. A four-sample example is instructional; larger, well-designed cohorts are needed for robust inference.

## Scientific checks

### Too few reads remain after trimming

Verify the actual library adapter, inspect raw FastQC overrepresented sequences, review Cutadapt logs, and check length/quality thresholds. Do not choose an adapter solely because it appears in an example.

### Unexpectedly many piRNA calls

Length and database overlap are not sufficient evidence of piRNA biogenesis. Review multi-mapping, strandedness, repeat context, cluster structure, 1U/10A biases, ping-pong overlap, phasing, and PIWI expression before making biological claims.

### A novel miRNA name looks official

Legacy names are putative, and HPC `novel:<locus>:<hash>` IDs are deliberately provisional. Validate structure, expression, reproducibility, contamination/Rfam screens, locus context, and current database status before applying official nomenclature.

If a failure persists, collect the repository commit, Nextflow/container versions, redacted parameters, validated sample sheet, failed process name, `.command.sh`, `.command.err`, and relevant scheduler log before opening a GitHub issue. Do not upload controlled FASTQs or credentials.
