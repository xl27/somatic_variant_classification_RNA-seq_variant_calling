#!/usr/bin/env Rscript
repo_root <- normalizePath(getwd())
sample_files <- read.table(file.path(repo_root, 'features', 'samples_id.txt'), stringsAsFactors = FALSE)$V1

for (sample_id in sample_files) {
  input_filtered <- file.path(repo_root, 'features_new', paste0(sample_id, '_features_filtered_GT_RD_MQ_AF_vcf_exons_jointAF1e-4.csv'))
  input_all <- file.path(repo_root, 'features_new', paste0(sample_id, '_features_filtered_GT_RD_MQ_AF_vcf_exons+jointAF.csv'))

  if (!file.exists(input_filtered) || !file.exists(input_all)) {
    message('Skipping ', sample_id, ' because required files are missing.')
    next
  }

  features_fil <- read.csv(input_filtered, row.names = 1)
  features_df <- read.csv(input_all, row.names = 1)
  features_germ <- features_df[features_df$AF >= 1e-2, ]
  features_germ$chr <- sapply(strsplit(rownames(features_germ), '_'), '[[', 1)
  features_germ$pos <- as.numeric(sapply(strsplit(rownames(features_germ), '_'), '[[', 2))
  features_germ <- features_germ[, c('chr', 'pos', 'ref_freq')]

  fil_vars <- data.frame(chr = sapply(strsplit(rownames(features_fil), '_'), '[[', 1),
                         pos = as.numeric(sapply(strsplit(rownames(features_fil), '_'), '[[', 2)))

  get_neighbor_vaf <- function(i) {
    germ_chr <- features_germ[features_germ$chr == fil_vars$chr[i], ]
    germ_chr_left <- germ_chr[(germ_chr$pos - fil_vars$pos[i]) < 0 & (germ_chr$pos - fil_vars$pos[i]) > -4e6, ]
    germ_chr_right <- germ_chr[(germ_chr$pos - fil_vars$pos[i]) > 0 & (germ_chr$pos - fil_vars$pos[i]) < 4e6, ]
    left_vaf <- ifelse(nrow(germ_chr_left) != 0, germ_chr_left$ref_freq[which.min(abs(germ_chr_left$pos - fil_vars$pos[i]))], 0.5)
    right_vaf <- ifelse(nrow(germ_chr_right) != 0, germ_chr_right$ref_freq[which.min(abs(germ_chr_right$pos - fil_vars$pos[i]))], 0.5)
    c(left_vaf, right_vaf)
  }

  vafs <- as.data.frame(t(sapply(seq_len(nrow(fil_vars)), get_neighbor_vaf)))
  features_fil$vaf_left <- vafs[, 1]
  features_fil$vaf_right <- vafs[, 2]

  write.csv(features_fil, file.path(repo_root, 'features_new', paste0(sample_id, '_features_filtered_GT_RD_MQ_AF_vcf_exons_jointAF_vafs.csv')), row.names = TRUE, quote = FALSE)
  message('Updated VAF for ', sample_id)
}
