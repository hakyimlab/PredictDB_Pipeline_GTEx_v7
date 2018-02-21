# Chromosome 12 has a duplicate snp for some reason
zless gtex_v7_eur_imputed_maf_0.01_R2_0.8_chr12.txt.gz | head -n1 | gzip -c > new_chr12.txt.gz
zless gtex_v7_eur_imputed_maf_0.01_R2_0.8_chr12.txt.gz | tail -n+2 | sort | uniq -w 23 | gzip -c >> new_chr12.txt.gz
