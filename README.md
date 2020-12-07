# PredictDB_Pipeline_GTEx_v7
PredictDB Pipeline for GTEx v7 data

You may want to use this [tutorial](https://github.com/hakyimlab/PredictDB-Tutorial) which contains the full code 

--------------

Comparison of V6p and V7 models can be found [here](http://hakyimlab.org/post/2017/v7-v6p-analysis/)

---------------
Scott's write up

GTEx v7 PredictDB Models - README
=================================

Scott Dickinson
11/13/2017

## Overview

This directory contains gene expression prediction models for 48 different human
tissues, meant for use with the Im Lab's software PrediXcan and S-PrediXcan.
These models were trained using GTEx v7 data, subsampled to use European samples
only. Genotype data was also imputed using the University of Michigan Imputation
Server, with HRC (Version r1.1 2016) as a reference panel. For each model, there
is also a corresponding file with covariance data for the SNPs in each model.
These files are necessary to use S-PrediXcan.

## Methods

Genotype data - GTEx genotype data was downloaded from dbGaP in the form of a
vcf file with all original 635 samples. GTEx standard QC procedures had already
been applied at time of download.  We then filtered the samples to include
Europeans only, using race as the filtering criteria. The data was then uploaded
to the University of Michigan Imputation server for imputation, with the
specific goal of obtaining more complete coverage of HapMap SNPs. We selected
HRC (Version r1.1 2016) as the reference panel. SNPs were then filtered to those
with minor allele frequency >= 1%, R2 >= 0.8, biallelic only, and those which
would be unambiguously stranded, i.e. the polymorphism is not to the base's
complement. SNPs were then mapped to rsids in dbSNP version 150, using position
and alleles as matching criteria. Once genotype data was subsetted again when
doing analysis on a tissue, SNPs were again filtered to have minor allele
frequency >= 1% according to the samples that had RNA sequencing performed for
that tissue. For analysis, genotypes were encoded on a continuous range from 0
to 2, with the value denoting the estimated count of the second, or effect,
allele.

Expression data - Expression data for each tissue was quantified by the GTEx
consortium, and we made extensive use of their pipeline for eQTL analysis, which
can be found at https://github.com/broadinstitute/gtex-pipeline. We normalized
the RPKM values using the script normalize_expression.py. In short, this script
quantile normalizes the expression values to the average empirical distribution
observed across samples, and then for each gene, expression values are inverse
quantile normalized to a standard normal distribution across samples. Genes were
selected based on expression thresholds of >0.1 RPKM in >=10 samples and >5
reads in >=10 samples.

Covariates - Before training models to predict gene expression, expression
values were adjusted for the following covariates:

    PEER Factors - PEER Factors were calculated using the GTEx pipeline docker
    container.  The number of PEER factors was chosen according to GTEx v7
    protocol with the number of samples being the determining factor. If the
    number of samples was greater than or equal to 350, we used 60 PEER factors.
    If the number of samples was between 250 and 350, we used 45. Between 150
    and 250, we used 30, and less than 150 we used 15.

    Sex - As specified by GTEx Subject Phenotypes file

    Genetic Principal Components - We used the top 3 genetic principal
    compoonents as computed by GTEx.

    Sequencing Platform - Whole Genome Sequencing for samples was carried out on
    two different types of machines. This information can be found in the Sample
    Attributes file.

Expression was adjusted by performing a multivariate linear regression with all
covariates, pulling the residual values, and then assigning the residuals to be
the new expression values.

Annotations - Gene annotation was derived from gencode v19, using GTEx's
collapsed gene model. SNP annotation was derived from the vcf from GTEx, a
HapMap annotation file for the CEU population, and data from dbSNP version 150.

Nested Cross Validated Elastic-Net - In previous versions of PredictDB, we
employed 10-fold cross-validated elastic-net to tune the parameter lambda, and
then estimated the significance of the model. It recently became apparent that
this was biasing the significance measures because we were using the same data
to tune the parameter lambda and assess the performance. To correct for this
problem, we used the following "nested" cross validation procedure:

1. Randomly split the data into 5 folds.

2. For each fold:

    a. Remove the fold from the data.

    b. Use the remaining data to train an elastic-net model using 10-fold
    cross-validation to tune the lambda parameter.

    c. With the trained model, predict on the hold out fold, and get various
    test statistics for how the model performs.

3. Calculate the average and standard deviation of each of the significance
statistics, where applicable. This should provide a reasonable estimate for how
well the model will generalize to new data.

4. Train a new elastic-net model using all of the data. Again, use 10-fold cross
validation to tune the lambda parameter. The non-zero weights from this final
model are what are saved in the database, provided the model meets significance
criteria.

A model was determined to be "significant" if the average pearson correlation
between predicted and observed during nested cross validation was greater than
0.1 (equivalent to R2 > 0.01) and the estimated p-value for this statistic was
less than 0.05. See below for how the p-value was calculated.

## Schema

### Weights Table

gene - Ensembl ID of the gene

rsid - rsid of the SNP. rsids are from dbSNP version 150.

varID - String label of the format chromosome_position_allele1_allele2_build.
All varIDs are from build 37 of the Human Reference Genome.

ref_allele - Allele which appears on Reference Genome

eff_allele - Alternative/Effect allele.

weight - The weight for this SNP in predicting expression for the gene. In
predicting the expression for the gene, the weight is multiplied by the count,
or estimated count, of the effect allele in the individual. This value is added
to all other weighted SNPs in the model.

### Extra Table

This table contains summary statistics about each gene model for the tissue

gene - Ensembl ID of the gene

genename - HUGO symbol for the gene

gene_type - Type of the gene according to Gencode v19. We have only included
models for protein coding genes, pseudogenes, and lincRNAs.

alpha - The alpha parameter used with the R package glmnet to train the elastic
net model. Here, we used 0.5 for every gene.

n_snps_in_window - The number of cis-SNPs examined for model training.  This is
the number of SNPs found within 1 million base pairs upstream of the
transcription start site and 1 million base pairs downstream of the
transcription end site. Sites are determined from Gencode v19.

n.snps.in.model - The number of SNPs within the cis window that have
non-zero weights, as found by elastic-net.

lambda_min_mse - The lambda parameter for elastic-net found to have the minimum
average mean square error on the hold out folds during cross validation.

test_R2_avg - The average coefficient of determination when predicting values of
the hold out fold during nested cross validation.  R2 is defined as

  1 - sum((y_observed-y_predicted)**2) / sum((y_observed-mean(y_observed)**2)

NB: It is possible for this value to be negative and still achieve a positive
correlation between predicted and observed.

test_R2_sd - The standard deviation of the coefficient of determination for each
of the five hold out folds during nested cross- validation.

cv_R2_avg - The average coefficient of determination for each of the hold out
folds when cross-validation was performed on the entire data set.

cv_R2_sd - The standard deviation of the coefficient of determination for each
of the 10 hold out folds when doing cross- validation on the entire data set.

in_sample_R2 - After an optimal lambda parameter has been chosen, cv.glmnet
trains an elastic net model with all of the data provided, i.e. without any hold
out folds. Expression is then predicted using this model, and the coefficient of
determination is calculated.

nested_cv_fisher_pval - Our first attempt at calculating a p-value during nested
cross-validation. We calculated a p-value using a correlation test using on each
of the five hold out folds, and combined the five p-values using Fisher's
method. Since some models yielded predictions with a variance of 0, we were
unable to perform a correlation test, and the p-value we reported for that fold
was a random number from the uniform distribution on the interval (0,1), which
would be the correct distribution under the assumption of the null hypothesis.
Under simulated null results, we still found this method did not produce a
uniform p-value distribution, with values under 0.05 and over 0.95 being
under-represented. We decided not to use this statistic because it did not quite
produce the desired results and the random element involved when there was prior
knowledge that a model did not predict well.

rho_avg - Average correlation between predicted and observed on the hold out
folds when doing nested cross-validation

rho_se - Standard deviation of correlation between predicted and observed on the
hold out folds when doing nested cross-validation

rho_zscore - Transformation of rho_avg into a z-score using Stouffer's Method.

pred.perf.R2 - rho_avg squared

pred.perf.pval - p-value for rho_zscore

pred.perf.qval - Deprecated. Previously held q-values calculated from the
distribution of p-values, but it was later deemed this analysis was
inappropriate.

### Construction Table

This table contains info solely meant for reproducibility purposes.

chromosome - Chromosome being analyzed

cv_seed - seed that was set at start of script for randomly splitting samples
into different folds


### Sample Info Table

Basic info on samples used for training. Contains a single row

n_samples - Number of samples

population - The population studied

tissue - The tissue from which RNA was sequenced

## Citations

Das S, Forer L, Schönherr S, Sidore C, Locke AE, Kwong A, Vrieze S, Chew EY,
Levy S, McGue M, Schlessinger D, Stambolian D, Loh PR, Iacono WG, Swaroop A,
Scott LJ, Cucca F, Kronenberg F, Boehnke M, Abecasis GR, Fuchsberger C.
Next-generation genotype imputation service and methods. Nature Genetics 48,
1284-1287 (2016).

GTEx Pipeline https://github.com/broadinstitute/gtex-pipeline


