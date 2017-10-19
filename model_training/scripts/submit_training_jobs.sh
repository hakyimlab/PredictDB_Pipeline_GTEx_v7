#!/bin/bash

for tiss in `cat gtex_tissues.txt`
do
    echo $tiss
    for chrom in {1..22}
    do
        echo -e "    $chrom"
        qsub gtex_tiss_chr_elasticnet.pbs -N gtex_training_${tiss}_${chrom} -v tiss=${tiss},chrom=${chrom}
        sleep .5
    done
done
