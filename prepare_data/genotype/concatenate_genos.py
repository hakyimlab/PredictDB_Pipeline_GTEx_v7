#!/usr/bin/env python
__author__ = "alvaro barbeira"

import os
import gzip

IR="gtex_v8_eur_shapeit2_phased_maf01_qdimputed.chr{}.txt.gz"

files = [IR.format(x) for x in range(1,23)]

with gzip.open("gtex_v8_eur_shapeit2_phased_maf01_qdimputed.txt.gz", "w") as o_:
    for file_ in files:
        print(file_)
        with gzip.open(file_) as i_:
            for i,line in enumerate(i_):
                if i==0:
                    if "chr1.txt.gz" in file_:
                        o_.write(line)
                    continue
                o_.write(line)
