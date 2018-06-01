#setwd("/group/im-lab/nas40t2/scott/gtex_v7_imputed_europeans/model_training/scripts/")
source("gtex_v7_nested_cv_elnet.R")

snp_annot_file <- "../../prepare_data/genotype/gtex_v8_eur_shapeit2_phased_maf01_snp_annot.chr22.txt"
gene_annot_file <- "../../prepare_data/expression/gencode_v26_parsed.txt"
genotype_file <- "../../prepare_data/genotype/gtex_v8_eur_shapeit2_phased_maf01_qdimputed_maf0.01_chr22.txt.gz"
expression_file <- "../../prepare_data/expression/Whole_Blood_Analysis.expression.txt"
covariates_file <- "../../prepare_data/covariates/Whole_Blood_Analysis.combined_covariates.txt"
chrom <- 22
prefix1 <- "Whole_Blood_nested_cv"
prefix2 <- "Whole_Blood_nested_cv_permuted"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, chrom, prefix1)#prefix2, null_testing=TRUE)

