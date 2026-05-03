#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
oneKG_dir="${ONEKG_DIR:-/gpfs/commons/datasets/1000genomes/hg38}"
save_dir="${SAVE_DIR:-$repo_root/oneKG_germline}"

mkdir -p "$save_dir"
cd "$save_dir"

for i in $(seq 1 22); do
  zcat "$oneKG_dir/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz" \
    | grep -v '^##' | cut -f1-8 > "chr${i}.GRCh38_dbSNP.txt"
done

Rscript --vanilla <<'RSCRIPT'
variants_list <- list()
for (i in 1:22) {
  oneKG_chr <- read.table(paste0('chr', i, '.GRCh38_dbSNP.txt'), stringsAsFactors = FALSE)
  oneKG_chr <- oneKG_chr[(nchar(oneKG_chr$V4) == 1) & (nchar(oneKG_chr$V5) == 1) & (oneKG_chr$V6 == 100) & (oneKG_chr$V7 == 'PASS'), ]
  variants_list[[i]] <- paste0('chr', oneKG_chr$V1, '>', oneKG_chr$V2)
  print(i)
}
variants_all <- Reduce(c, variants_list)
saveRDS(variants_list, file = 'All_SNPs_filtered.rds')
write.table(variants_all, 'All_SNPs_filtered.txt', col.names = FALSE, quote = FALSE)
RSCRIPT

echo "OneKG processing completed in $save_dir"
