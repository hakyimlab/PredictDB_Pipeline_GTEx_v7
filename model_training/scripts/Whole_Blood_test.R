setwd("/group/im-lab/nas40t2/scott/gtex_v7_imputed_europeans/model_training/scripts/")
source("gtex_v7_nested_cv_elnet.R")

snp_annot_file <- "../../prepare_data/genotype/gtex_v7_hapmapceu_snp_annot.chr22.txt"
gene_annot_file <- "../../prepare_data/expression/gencode.v19.genes.patched_contigs.parsed.txt"
genotype_file <- "../../prepare_data/genotype/gtex_v7_eur_imputed_maf_0.01_R2_0.8_chr22.txt.gz"
expression_file <- "../../prepare_data/expression/Whole_Blood_Analysis.expression.txt"
covariates_file <- "../../prepare_data/covariates/Whole_Blood_Analysis.combined_covariates.txt"
chrom <- 22
prefix1 <- "Whole_Blood_nested_cv"
prefix2 <- "Whole_Blood_nested_cv_permuted"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, chrom, prefix2, null_testing=TRUE)

