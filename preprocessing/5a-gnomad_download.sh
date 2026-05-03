#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$repo_root"

chr_list=( {1..22} X Y )
index="${SLURM_ARRAY_TASK_ID:-1}"
chr="chr${chr_list[$((index-1))]}"

wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.${chr}.vcf.bgz
