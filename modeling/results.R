setwd("/gpfs/commons/groups/gursoy_lab/xli/somatic_variant_prediction/")
library(ggplot2)

model_ls <- c("xgboost-weights_noSampleFilter",
              "xgboost-SMOTE_noSampleFilter",
              "xgboost-weights_Filter-num20",
              "xgboost-SMOTE_Filter-num20",
              "xgboost-weights_Filter-tcga",
              "xgboost-SMOTE_Filter-tcga")


sample_files <- read.table("features/samples_id.txt")$V1
true_summary <- sapply(sample_files, function(id) nrow(read.table(paste0("match_ture_tcga/", id, "/tcga_true_variants.txt"))))
sample_id <- gsub("-01A|-01B|-01C", "", sample_files)
names(true_summary) <- sample_id

cases_info <- read.csv("tcga_info/cases_info.txt", sep='\t')
cases_info <- cases_info[cases_info$submitter_id %in% sample_id,]
cases_info$num_true_somatic <- true_summary[cases_info$submitter_id]


sample_sizes <- tapply(cases_info$submitter_id, cases_info$disease_type, length)
sample_sizes <- data.frame(disease_type = names(sample_sizes), sample_size = sample_sizes)
sample_sizes <- sample_sizes[order(sample_sizes$sample_size, decreasing = F),]
sample_sizes$disease_type <- factor(sample_sizes$disease_type, levels = sample_sizes$disease_type)



#metrics_1 <- read.csv(paste0("models/newfiltered_models/filter_models/", model_ls[1], "_metrics-by-disease.csv"))
#metrics_1$ratio <- metrics_1$num_class_1/metrics_1$total_cases
#metrics_1 <- metrics_1[match(sample_sizes$disease_type, metrics_1$disease_type, nomatch = 0),]
#metrics_1$disease_type <- factor(metrics_1$disease_type, levels = metrics_1$disease_type)



metrics_ls <- lapply(model_ls, function(m) {
  df <- read.csv(paste0("models/newfiltered_models/filter_models/", m, "_metrics-by-disease.csv"))
  df$ratio <- df$num_class_1/df$total_cases
  df$model_filter <- m
  return(df)
  })
metrics_df <- Reduce(rbind, metrics_ls)

metrics_df$filter <- sapply(strsplit(metrics_df$model_filter, "_"), "[[", 2)
metrics_df$filter <- factor(metrics_df$filter, levels = c("noSampleFilter", "Filter-num20", "Filter-tcga"))
levels(metrics_df$filter) <- c("No filter", "#somatic > 20", "#TCGA 0.1-0.9")
metrics_df$model <- sapply(strsplit(metrics_df$model_filter, "_"), "[[", 1)


pdf("models/newfiltered_models/filter_models/f1_ratio.pdf", 14, 8)
ggplot(metrics_df, aes(ratio, f1_score, color = disease_type)) +
  geom_point() +
  geom_text(aes(label = disease_type), size = 1.5, vjust = 2) +
  facet_wrap(~model+filter, scales = "free") +
  labs(y = "F1 score", x = "Ratio") +
  theme_bw()
dev.off()


pdf("models/newfiltered_models/filter_models/f1_samplesize.pdf", 14, 8)
ggplot(metrics_df, aes(sample_size, f1_score, color = disease_type)) +
  geom_point() +
  geom_text(aes(label = disease_type), size = 1.5, vjust = 2) +
  facet_wrap(~model+filter, scales = "free") +
  labs(y = "F1 score", x = "Sample size") +
  theme_bw()
dev.off()


pdf("models/newfiltered_models/filter_models/auprc_ratio.pdf", 14, 8)
ggplot(metrics_df, aes(ratio, auprc, color = disease_type)) +
  geom_point() +
  geom_text(aes(label = disease_type), size = 1.5, vjust = 2) +
  facet_wrap(~model+filter, scales = "free") +
  labs(y = "AUPRC", x = "Ratio") +
  theme_bw()
dev.off()



per_ls <- read.csv(paste0("models/newfiltered_models/filter_models/", model_ls[5], "_perDisease_metrics-by-disease.csv"))
per_df$ratio <- per_df$num_class_1/per_df$total_cases
per_df$model_filter <- model_ls[5]


per_ls <- lapply(model_ls[5:6], function(m) {
  df <- read.csv(paste0("models/newfiltered_models/filter_models/", m, "_perDisease_metrics-by-disease.csv"))
  df$ratio <- df$num_class_1/df$total_cases
  df$model_filter <- m
  metrics_sub <- metrics_df[metrics_df$model_filter == m & metrics_df$disease_type %in% df$disease_type,]
  df$previous_f1_score <- metrics_sub$f1_score[match(df$disease_type, metrics_sub$disease_type)]
  return(df)
})
per_df <- Reduce(rbind, per_ls)
per_df$model <- sapply(strsplit(per_df$model_filter, "_"), "[[", 1)


pdf("models/newfiltered_models/filter_models/f1_perDis.pdf", 8, 7.5)
ggplot(per_df, aes(f1_score, previous_f1_score, color = disease_type)) +
  geom_abline(slope = 1, intercept = 0, color = "grey") +
  geom_point() +
  geom_text(aes(label = disease_type), size = 1.5, vjust = 2) +
  facet_wrap(~model, ncol = 1) +
  labs(y = "Train: per disease", x = "Train: all", title = "F1 Score") +
  theme_bw()
dev.off()




