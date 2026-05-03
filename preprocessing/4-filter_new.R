#!/usr/bin/env Rscript
repo_root <- normalizePath(getwd())
library(vcfR)

features_dir <- file.path(repo_root, 'features')
filtered_dir <- file.path(repo_root, 'vcf_filtered')
output_dir <- file.path(repo_root, 'features_new')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_files <- read.table(file.path(features_dir, 'samples_id.txt'), stringsAsFactors = FALSE)$V1
n_all_ls <- list()
n_true_ls <- list()

for (sample_id in sample_files) {
  features_df <- readRDS(file.path(features_dir, paste0(sample_id, '_features_processed.rds')))
  n_all <- nrow(features_df)
  n_true <- sum(features_df$if_true_somatic)

  features_df <- features_df[features_df$gt == '0/1', ]
  n_all <- c(n_all, nrow(features_df))
  n_true <- c(n_true, sum(features_df$if_true_somatic))

  features_df <- features_df[features_df$`read_depth-11` >= 12, ]
  n_all <- c(n_all, nrow(features_df))
  n_true <- c(n_true, sum(features_df$if_true_somatic))

  features_df <- subset(features_df, (ref_mapping_qual_mean == 60) & (alt_mapping_qual_mean == 60))
  if (ncol(features_df) >= 74) {
    features_df <- features_df[, -c(71:74)]
  }
  n_all <- c(n_all, nrow(features_df))
  n_true <- c(n_true, sum(features_df$if_true_somatic))

  features_df <- subset(features_df, ref_freq <= 0.9)
  n_all <- c(n_all, nrow(features_df))
  n_true <- c(n_true, sum(features_df$if_true_somatic))

  filtered_vcf_path <- file.path(filtered_dir, paste0(sample_id, '.filtered.vcf.gz'))
  if (file.exists(filtered_vcf_path)) {
    filtered_vcf <- read.vcfR(filtered_vcf_path)
    passed <- filtered_vcf@fix[, 'FILTER'] == 'PASS'
    filtered_id <- apply(filtered_vcf@fix[, 1:2], 1, paste, collapse = '>')
    features_df <- features_df[rownames(features_df) %in% filtered_id[passed], ]
  }
  n_all <- c(n_all, nrow(features_df))
  n_true <- c(n_true, sum(features_df$if_true_somatic))

  n_all_ls[[sample_id]] <- n_all
  n_true_ls[[sample_id]] <- n_true

  write.csv(features_df, file = file.path(output_dir, paste0(sample_id, '_features_filtered_GT_RD_MQ_AF_vcf.csv')), row.names = TRUE, quote = FALSE)
}

all_num <- do.call(rbind, n_all_ls)
true_num <- do.call(rbind, n_true_ls)
rownames(all_num) <- rownames(true_num) <- sample_files
save(all_num, true_num, file = file.path(output_dir, 'summary_filtered_GT_RD_MQ_AF_vcf.rData'))

cat('All totals:', apply(all_num, 2, sum), '\n')
cat('True totals:', apply(true_num, 2, sum), '\n')
