#!/bin/bash

FAILED=not_ran.txt
#| awk '{ printf("%s", $0) }'
tissue=($(cat $FAILED | cut - -d " " -f 1))
chrom=($(cat $FAILED | cut - -d " " -f 2))

L=${#tissue[@]}
for i in $(eval echo "{0..$L}"); do
  ti_s=${tissue[$i]}
  ch_s=${chrom[$i]}
  echo "${ti_s} ${ch_s}"
  rm -f "../weights/${ti_s}_nested_cv_chr${ch_s}_weights.txt"
  rm -f "../summary/${ti_s}_nested_cv_chr${ch_s}_*"
  rm -f "../covariances/${ti_s}_nested_cv_chr${ch_s}_covariances.txt"
  qsub gtex_tiss_chr_elasticnet.pbs -N gtex_training_${ti_s}_${ch_s} -v tiss=${ti_s},chrom=${ch_s}
  sleep .5
done
