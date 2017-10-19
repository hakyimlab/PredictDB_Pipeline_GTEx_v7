#!/bin/bash

module load gcc/4.9.4
module load bcftools/1.4.0

bcftools concat -f vcf_files_by_chr.txt -Oz -o gtex_v7_eur_imputed.vcf.gz
