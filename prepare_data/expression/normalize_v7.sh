#!/bin/bash

counts="../v7/54079/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2016-01-15_v7/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"
rpkm="../v7/54079/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2016-01-15_v7/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct"
gtf="gencode.v19.genes.v6p.patched_contigs.gtf"
vcf="../v7/56055/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2016-01-15_v7/genotypes/WGS/variant_calls/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz"

for tissue_donor_f in `ls *_donors.txt`;
do
    tiss=`echo $tissue_donor_f | sed 's/_donors.txt//' | sed 's/[()]//g'`
    echo $tiss
    python normalize_expression.py $rpkm $counts $gtf $tissue_donor_f ${tiss}_Analysis $vcf --min_samples 10 &
done
