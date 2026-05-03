#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
data_dir="${DATA_DIR:-/gpfs/commons/groups/gursoy_lab/cwalker/projects/somatic_variant_prediction/data}"
save_dir="${SAVE_DIR:-$repo_root/match_ture_tcga}"

mkdir -p "$save_dir"
cd "$data_dir/tcga/true_somatic_mutations"

rm -f "$save_dir/tcga_id.txt"
for dir in */; do
  file_id="${dir%/}"
  tcga_id=$(zcat "${dir}"/*.maf.gz | grep 'TCGA-' | cut -f16 | cut -d'-' -f1-4 | uniq)
  if [[ -z "$tcga_id" ]]; then
    continue
  fi
  if [[ -f "${data_dir}/seven_bridges/processed_vcfs/raw/${tcga_id}.vcf.gz" ]]; then
    echo "${file_id} ${tcga_id}" >> "$save_dir/tcga_id.txt"
  fi
done

echo "Wrote $save_dir/tcga_id.txt"
