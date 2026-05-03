#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

chr_list=( {1..22} X Y )
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  index="$SLURM_ARRAY_TASK_ID"
  chr="chr${chr_list[$((index-1))]}"
  echo "Processing $chr"
  module load R/4.3.1
  Rscript --vanilla ./5b-gnomad_match_AF_chr.R "$chr"
else
  module load R/4.3.1
  for chr in "${chr_list[@]}"; do
    echo "Processing $chr"
    Rscript --vanilla ./5b-gnomad_match_AF_chr.R "$chr"
  done
fi
