# Docker notebook workflow

Use this workflow to run the original interactive `mirpipe_pipeline.ipynb` in the published miRPipe container. It is most suitable for a first run, a small experiment, or comparison with historical results.

## 1. Check the host

The legacy workflow was tested on 64-bit Linux and needs at least 8 GB RAM. Install Docker using the instructions for your distribution, then verify it:

```bash
docker version
docker run --rm hello-world
```

If Docker requires `sudo`, either use `sudo` for Docker commands or configure Docker's non-root access according to your institution's policy.

## 2. Prepare a data directory

Create one host directory containing all single-end FASTQ or FASTQ.GZ files and a comma-separated file named `sample_list.csv`:

```text
data/
├── CLL_T01.fastq.gz
├── CLL_T02.fastq.gz
├── CLL_U01.fastq.gz
├── CLL_U02.fastq.gz
└── sample_list.csv
```

The notebook expects these exact columns:

```csv
sample,file,condition
CLL_T01,CLL_T01.fastq.gz,treated
CLL_T02,CLL_T02.fastq.gz,treated
CLL_U01,CLL_U01.fastq.gz,untreated
CLL_U02,CLL_U02.fastq.gz,untreated
```

For the legacy notebook, keep sample IDs unique and sort rows alphanumerically by `sample`. The `file` value must exactly match a file in the mounted directory. Do not use the HPC column name `fastq` here.

Quick checks:

```bash
cd /path/to/miRPipe
find data -maxdepth 1 -type f -name '*.fastq*' -print
head -5 data/sample_list.csv
```

## 3. Pull and record the container image

```bash
docker pull docker.io/vivekruhela/mirpipe:latest
docker image inspect \
  --format '{{index .RepoDigests 0}}' \
  docker.io/vivekruhela/mirpipe:latest \
  | tee data/mirpipe-container-digest.txt
```

The digest makes the exact pulled image auditable. For a publication run, launch the recorded digest instead of relying on the mutable `latest` tag.

## 4. Start miRPipe

Run this command from the repository root so `$(pwd)/data` is the prepared directory:

```bash
docker run --rm --name mirpipe \
  -p 127.0.0.1:8880:8888 \
  -e PASSWORD='choose-a-password' \
  -e USE_HTTP=1 \
  -v "$(pwd)/data:/miRPipe/data" \
  docker.io/vivekruhela/mirpipe:latest
```

What the options do:

| Option | Purpose |
|---|---|
| `--rm` | Remove the stopped container; mounted results remain on the host |
| `--name mirpipe` | Give the running container a predictable name |
| `-p 127.0.0.1:8880:8888` | Expose Jupyter only on the local host |
| `PASSWORD=...` | Set the Jupyter login password |
| `USE_HTTP=1` | Enable the container's HTTP notebook service |
| `-v ...:/miRPipe/data` | Persist inputs and outputs in the host `data/` directory |

Do not expose this HTTP service directly to a public network. On a remote machine, bind it to localhost and use an SSH tunnel approved by your institution.

## 5. Open and run the notebook

1. Open <http://localhost:8880/mirpipe>.
2. Enter the password supplied to `docker run`.
3. Open `mirpipe_pipeline.ipynb`.
4. Confirm that `/miRPipe/data` contains every FASTQ and `sample_list.csv`.
5. Run cells in order. Do not skip cells that initialize paths or reference choices.
6. Enter the library's actual 3′ adapter when prompted.
7. Select a genome/miRBase option that matches every downstream reference.
8. Keep the terminal open until the notebook and background jobs finish.

For the historical CLL demonstration, the tutorial uses adapter `TGGAATTCTCGGGTGCCAAGG` and its historical option 1 reference combination. For a new hg38 analysis, the notebook's option 4 selects hg38 with miRBase 22. Treat this as an analysis decision, not a convenience default.

## 6. Verify the result

The CLL differential-expression example retains five principal result files:

```text
data/output/
├── final_diff_exp_miRNAs.csv
├── diff_exp_miRNAs_expression_counts.csv
├── miRNA_expression_counts.csv
├── piRNA_raw_counts.csv
└── significantly_DE_piRNA.csv
```

Before interpretation:

1. Check that every sample is present in the miRNA and piRNA count matrices.
2. Confirm count-matrix columns align with the sorted `sample_list.csv` rows.
3. Review adapter-trimming and FastQC outputs for unexpected read loss or contamination.
4. Confirm the intended untreated/control level was used in the differential-expression comparison.
5. Preserve `sample_list.csv`, the image digest, notebook selections, and unfiltered count tables.

Novel names emitted by the legacy notebook are putative. Recheck candidates against current, versioned reference databases before treating them as official identifiers.

## 7. Stop or rerun

Press `Ctrl+C` in the Docker terminal, or from another terminal run:

```bash
docker stop mirpipe
```

Because `data/` is mounted, files written there survive container removal. The notebook itself is not a resumable workflow engine: inspect partial outputs before rerunning, use a fresh output directory for a clean rerun, and never combine files from different parameter choices.

Next: follow the [worked example](Worked-Example.md) or review [Docker troubleshooting](Troubleshooting.md#docker-and-jupyter).
