#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
save_dir="${SAVE_DIR:-$repo_root/match_ture_tcga}"
summary_file="${SUMMARY_FILE:-$save_dir/matched_summary.txt}"
offset="${OFFSET:-0}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID is not set; run this script inside a Slurm array or set it manually." >&2
  exit 1
fi

index=$((SLURM_ARRAY_TASK_ID + offset))
tcga_id=$(sed -n "${index}p" "$summary_file" | awk '{print $1}')

if [[ -z "$tcga_id" ]]; then
  echo "No TCGA ID for index $index in $summary_file" >&2
  exit 1
fi

cd "$repo_root"
module load R/4.3.1
Rscript --vanilla ./3-process_features.R "$tcga_id"
