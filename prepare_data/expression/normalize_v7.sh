#!/bin/bash

counts="/Volumes/gtex-group/v7/54079/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2016-01-15_v7/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"
rpkm="/Volumes/gtex-group/v7/54079/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2016-01-15_v7/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct"
gtf="gencode.v19.genes.v6p.patched_contigs.gtf"

for tissue_donor_f in `ls *_donors.txt`;
do
    tiss=`echo $tissue_donor_f | sed 's/_donors.txt//' | sed 's/[()]//g'`
    echo $tiss
    python normalize_expression.py $rpkm $counts $gtf $tissue_donor_f ${tiss}_Analysis --min_samples 10 &
done
