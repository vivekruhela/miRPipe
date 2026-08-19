#!/usr/bin/env bash
#SBATCH --job-name=miRPipe_driver
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --output=miRPipe_driver.%j.out
#SBATCH --error=miRPipe_driver.%j.err

set -euo pipefail

workflow_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

: "${NXF_WORK:=${workflow_dir}/work}"
: "${NXF_APPTAINER_CACHEDIR:=${workflow_dir}/.apptainer}"
export NXF_WORK NXF_APPTAINER_CACHEDIR

command -v nextflow >/dev/null 2>&1 || {
    echo 'Nextflow is not on PATH. Load the site Nextflow module before submitting.' >&2
    exit 127
}
command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1 || {
    echo 'Apptainer/Singularity is not on PATH. Load the site container module.' >&2
    exit 127
}

cd "${workflow_dir}"
nextflow run main.nf \
    -profile slurm \
    -params-file params.yml \
    -work-dir "${NXF_WORK}" \
    -resume
