#!/bin/bash

for i in {1..22}
do
    echo $i
    out_file=gtex_v7_eur_imputed_maf_0.01_R2_0.8_chr${i}.txt
    echo -ne "varID\t" > ${out_file} 
    bcftools query -l ./imputation_results/chr${i}.dose.vcf.gz | tr '\n' '\t' >> $out_file
    echo "" >> ${out_file}
    bcftools query -e  'MAF[0]<0.01 | INFO/R2<0.8 |  TYPE!="snp" | N_ALT!=1' -f '%CHROM\_%POS\_%REF\_%ALT\_b37[\t%DS]\n' ./imputation_results/chr${i}.dose.vcf.gz >> ${out_file}
    gzip ${out_file}
done
