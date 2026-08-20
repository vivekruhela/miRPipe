# Worked four-sample example

This example shows one two-group experiment represented for both workflows. It is a template: obtain appropriately consented data, replace the placeholder FASTQs, and verify the real library adapter and reference build.

The example uses two treated and two untreated biological samples:

```text
miRPipe/
├── example-data/
│   ├── CLL_T01.fastq.gz
│   ├── CLL_T02.fastq.gz
│   ├── CLL_U01.fastq.gz
│   └── CLL_U02.fastq.gz
└── hpc_workflow/
```

The adapter below, `TGGAATTCTCGGGTGCCAAGG`, is the adapter documented for the historical CLL demonstration. Confirm it from your own library protocol or adapter-detection QC before analysis.

## A. Run the example with Docker

### A1. Add the legacy sample sheet

Create `example-data/sample_list.csv`:

```csv
sample,file,condition
CLL_T01,CLL_T01.fastq.gz,treated
CLL_T02,CLL_T02.fastq.gz,treated
CLL_U01,CLL_U01.fastq.gz,untreated
CLL_U02,CLL_U02.fastq.gz,untreated
```

Rows are alphanumerically sorted because the notebook's differential-expression stage depends on legacy sample/count ordering.

### A2. Launch the notebook

From the repository root:

```bash
docker pull docker.io/vivekruhela/mirpipe:latest
docker run --rm --name mirpipe \
  -p 127.0.0.1:8880:8888 \
  -e PASSWORD='choose-a-password' \
  -e USE_HTTP=1 \
  -v "$(pwd)/example-data:/miRPipe/data" \
  docker.io/vivekruhela/mirpipe:latest
```

Open <http://localhost:8880/mirpipe> and run `mirpipe_pipeline.ipynb` in order. Enter the example adapter when prompted. For a published CLL comparison, use the historical option described in `Tutorial.md`; for a new hg38 analysis, use the matching hg38/miRBase choice and compatible references.

### A3. Check results

```bash
find example-data/output -maxdepth 1 -type f -print | sort
```

At minimum, check that all four sample columns occur exactly once in `miRNA_expression_counts.csv` and `piRNA_raw_counts.csv` before reading differential-expression results.

## B. Run the same example with HPC/Nextflow

### B1. Add the HPC sample sheet

Create `hpc_workflow/samplesheet.csv`:

```csv
sample,fastq,condition,batch
CLL_T01,../example-data/CLL_T01.fastq.gz,treated,batch1
CLL_T02,../example-data/CLL_T02.fastq.gz,treated,batch1
CLL_U01,../example-data/CLL_U01.fastq.gz,untreated,batch1
CLL_U02,../example-data/CLL_U02.fastq.gz,untreated,batch1
```

Relative FASTQ paths resolve from this sample sheet. On a cluster, absolute shared-storage paths are usually clearer.

### B2. Configure the run

From `hpc_workflow/`:

```bash
cp assets/params.example.yml params.yml
```

For an hg38 example, edit these fields while retaining the version-matched hg38 reference paths from the template:

```yaml
samplesheet: samplesheet.csv
outdir: results/cll_four_sample

adapter: TGGAATTCTCGGGTGCCAAGG
genome: hg38

run_mirna: true
run_pirna: true
run_deseq2: true
run_rfam: false

design: "~ condition"
contrast: condition,treated,untreated
min_count: 10
min_samples: 2
alpha: 0.05
```

Do not change only `genome`: all miRDeep, Bowtie, mature-miRNA, miRNA-GFF, and piRNA-GFF paths must change together when switching assemblies.

### B3. Submit locally or on SLURM

Local Docker:

```bash
nextflow run main.nf \
  -profile docker \
  -params-file params.yml \
  -resume
```

SLURM with the bundled wrapper:

```bash
sbatch submit_mirpipe.slurm
```

SLURM with site configuration:

```bash
nextflow run main.nf \
  -profile slurm \
  -c cluster.config \
  -params-file params.yml \
  -work-dir /shared/project/mirpipe/work \
  -resume
```

### B4. Review the run

```bash
test -s results/cll_four_sample/qc/multiqc/multiqc_report.html
test -s results/cll_four_sample/pipeline_info/execution_trace.tsv
find results/cll_four_sample/mirna/counts -maxdepth 1 -type f -print
find results/cll_four_sample/pirna/counts -maxdepth 1 -type f -print
```

Open MultiQC first. Then verify:

1. all four samples passed validation and preprocessing;
2. adapter trimming and length fractions are plausible for the library;
3. raw count matrices contain the four expected sample IDs;
4. the reported contrast is treated versus untreated;
5. filtered features and adjusted p-values are interpreted with the small sample size in mind;
6. candidate labels are treated as evidence categories, not official nomenclature.

## Comparing Docker and HPC outputs

Do not expect byte-identical output merely because the same FASTQs were used. The HPC implementation deliberately changes orchestration, sample matching, length boundaries, local annotation, provenance, and candidate identifiers. A fair comparison must hold adapter, genome, reference releases, thresholds, branches, and contrast constant, then compare QC, retained reads, count matrices, and calls stage by stage.
