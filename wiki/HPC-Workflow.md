# Automated HPC/SLURM workflow

The `hpc_workflow/` directory provides a non-interactive Nextflow DSL2 implementation for single-end small-RNA cohorts. It validates metadata, runs samples and stages in parallel, uses Docker or Apptainer/Singularity, retries selected resource failures, and supports cache-based `-resume`.

## 1. Check cluster prerequisites

You need:

- Nextflow 24.04.0 or newer and Java 11 or newer;
- SLURM for the `slurm` profile;
- Apptainer or Singularity on compute nodes;
- shared storage visible from the login/driver node and all compute nodes;
- compatible miRDeep\*, Bowtie, genome, miRNA, and piRNA references.

Load your site's modules, then verify the commands without launching an analysis:

```bash
module load java nextflow apptainer  # Names vary by cluster.
java -version
nextflow -version
apptainer --version || singularity --version
sbatch --version
```

Ask the cluster administrators which account, partition/queue, storage area, and wall-time limits to use. Do not run the full workflow on a login node.

## 2. Clone and enter the workflow

```bash
git clone https://github.com/vivekruhela/miRPipe.git
cd miRPipe/hpc_workflow
git rev-parse HEAD | tee mirpipe-commit.txt
```

Run all following commands from `hpc_workflow/` unless a command says otherwise.

## 3. Create the HPC sample sheet

```bash
cp assets/samplesheet.example.csv samplesheet.csv
```

Required columns are `sample`, `fastq`, and `condition`. Extra covariates such as `batch`, `sex`, `pair`, or `RIN` are retained for the DESeq2 design.

```csv
sample,fastq,condition,batch
CLL_T01,/shared/project/fastq/CLL_T01.fastq.gz,treated,batch1
CLL_T02,/shared/project/fastq/CLL_T02.fastq.gz,treated,batch1
CLL_U01,/shared/project/fastq/CLL_U01.fastq.gz,untreated,batch1
CLL_U02,/shared/project/fastq/CLL_U02.fastq.gz,untreated,batch1
```

Rules:

- one row per biological sample;
- unique sample IDs using letters, numbers, `.`, `_`, or `-`;
- one single-end FASTQ/FASTQ.GZ per row;
- absolute paths are safest on a cluster; relative paths resolve from the sample sheet;
- merge technical lanes before this alpha workflow;
- use biological replicates for differential expression.

The HPC workflow matches counts to metadata by sample name, so rows do not need legacy alphanumeric ordering. Do not use the Docker column name `file` here.

## 4. Configure parameters and references

```bash
cp assets/params.example.yml params.yml
```

Edit at least `samplesheet`, `outdir`, `adapter`, `genome`, the reference paths, `design`, and `contrast`:

```yaml
samplesheet: /shared/project/mirpipe/samplesheet.csv
outdir: /shared/project/mirpipe/results

run_mirna: true
run_pirna: true
run_deseq2: true
run_rfam: false

adapter: TGGAATTCTCGGGTGCCAAGG
genome: hg38

mirdeep_jar: /miRPipe/Tools/MDS_command_line_v38/MDS_command_line/MD.jar
mirdeep_genome_dir: /miRPipe/Tools/MDS_command_line_v38/MDS_command_line/genome
known_mature_fasta: /miRPipe/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase/mature.fa
known_mirna_gff: /miRPipe/Tools/MDS_command_line_v38/MDS_command_line/genome/hg38/miRBase/knownMiR.gff3
bowtie_index: /miRPipe/refs/hg38/hg38
pirna_gff: ../refs/pirnadb/pirnadb.hg38.gff3.gz

design: "~ batch + condition"
contrast: condition,treated,untreated
min_count: 10
min_samples: 2
alpha: 0.05

container: docker.io/vivekruhela/mirpipe:latest
```

Paths beginning `/miRPipe` above are inside the legacy-compatible container. Host reference paths must be visible on compute nodes and automatically bind-mounted by Apptainer. Verify that every genome-dependent file uses the selected assembly and chromosome naming convention.

For an archival run, replace `latest` with a fixed image digest or a locally managed SIF and record its SHA-256 checksum.

## 5. Add site-specific SLURM settings

The bundled `conf/slurm.config` submits one task per process with bounded submission rate and retry behavior. Keep site policy in a separate file:

```bash
cp conf/columbia.example.config cluster.config
```

Example:

```groovy
process {
    queue = 'cpu'
    clusterOptions = '--account=my_allocation'
}
```

Replace the examples with values supplied by your cluster. You can also override resources for one label without editing the workflow:

```groovy
process {
    withLabel: 'large' {
        cpus = 12
        memory = '32 GB'
        time = '36h'
    }
}
```

## 6. Choose persistent work and cache locations

Nextflow needs its launch metadata and work directory to resume. Container caching prevents repeated image downloads.

```bash
export NXF_WORK=/shared/project/mirpipe/work
export NXF_APPTAINER_CACHEDIR=/shared/project/mirpipe/apptainer-cache
mkdir -p "$NXF_WORK" "$NXF_APPTAINER_CACHEDIR"
```

Do not place these on short-lived node-local scratch if you intend to use `-resume`.

## 7. Submit the workflow

### Bundled driver job

After editing `params.yml`, submit:

```bash
sbatch submit_mirpipe.slurm
```

The driver verifies Nextflow and Apptainer/Singularity, then runs the `slurm` profile with `-resume`. To use `cluster.config`, either add `-c cluster.config` to the final `nextflow run` command in a site-local copy of the wrapper, or submit a driver script containing the direct command below.

### Direct command from an allocated driver session

```bash
nextflow run main.nf \
  -profile slurm \
  -c cluster.config \
  -params-file params.yml \
  -work-dir "$NXF_WORK" \
  -resume
```

For a workstation smoke run, use `-profile docker`; for local Apptainer, use `-profile apptainer`.

## 8. Monitor and resume

```bash
squeue -u "$USER"
tail -f miRPipe_driver.<job-id>.out
nextflow log
```

If the driver or a task fails:

1. Read the driver log and the failed task's `.command.err` in its work directory.
2. Correct the parameter, path, reference, container, or resource problem.
3. Rerun the same command with the same launch directory and work directory.
4. Keep `-resume`; successful tasks will be reused when their inputs and command signatures still match.

Do not delete `work/` or `.nextflow/` until the run is finalized.

## 9. Review outputs

```text
results/
├── preprocessing/                 # trimmed and length-split reads
├── qc/
│   ├── fastqc/
│   └── multiqc/multiqc_report.html
├── mirna/
│   ├── mirdeep_star/
│   ├── counts/
│   ├── annotation/
│   └── differential_expression/
├── pirna/
│   ├── alignment/
│   ├── counts/
│   └── differential_expression/
└── pipeline_info/
    ├── software_versions.txt
    ├── execution_trace.tsv
    ├── execution_report.html
    ├── execution_timeline.html
    └── workflow_dag.html
```

Start with MultiQC, then confirm sample completeness in merged counts and inspect model/contrast outputs before biological interpretation. A piRNAdb overlap is a **piRNA annotation**, not proof of active piRNA biogenesis. Novel miRNA evidence tiers prioritize candidates but do not create official miRBase names.

Next: run the [worked example](Worked-Example.md), check [HPC troubleshooting](Troubleshooting.md#nextflow-slurm-and-apptainer), and archive the [reproducibility checklist](Reproducibility-Checklist.md).
