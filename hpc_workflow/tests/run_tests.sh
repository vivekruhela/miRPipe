#!/usr/bin/env bash
set -euo pipefail

workflow_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python3 -m unittest discover -s "${workflow_dir}/tests" -p 'test_*.py' -v
python3 -m compileall -q "${workflow_dir}/bin" "${workflow_dir}/tests"
bash -n "${workflow_dir}/submit_mirpipe.slurm"
