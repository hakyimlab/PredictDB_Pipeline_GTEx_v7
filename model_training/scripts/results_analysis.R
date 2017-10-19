setwd("/Volumes/im-lab/nas40t2/scott/gtex_v7_imputed_europeans/model_training/scripts/")
library(ggplot2)
library(dplyr)

model_performance <- read.table("../summary/Whole_Blood_nested_cv_chr22_model_summaries.txt", stringsAsFactors = F, header = T)

ggplot(model_performance, aes(x = nested_cv_fisher_pval)) + geom_histogram(bins = 20, na.rm = T) + ggtitle("Whole Blood GTEx chr 22: Distribution of p-values")

null_model_performance <- read.table("../summary/Whole_Blood_nested_cv_permuted_chr22_model_summaries.txt", stringsAsFactors = F, header = T)

ggplot(null_model_performance, aes(x = nested_cv_fisher_pval)) + geom_histogram(bins = 20, na.rm = T) + ggtitle(("Whole Blood GTEx chr 22 (Permuted): Distribution of p-values"))


ggplot(model_performance, aes(in_sample_R2, test_R2_avg)) + geom_point(aes(color = nested_cv_fisher_pval < 0.05))

ggplot(null_model_performance, aes(in_sample_R2, test_R2_avg)) + geom_point(aes(color = nested_cv_fisher_pval < 0.05))

model_performance2 <- read.table("../summary/Whole_Blood_nested_cv2_chr22_model_summaries.txt", stringsAsFactors = F, header = T)

ggplot(model_performance2, aes(cv_R2_avg, test_R2_avg)) + geom_point(aes(color = nested_cv_fisher_pval < 0.05))
