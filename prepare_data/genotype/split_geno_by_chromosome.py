#!/usr/bin/env python
__author__ = "alvaro barbeira"
import os
import gzip
from convert_gtex_annotation_to_pred import read_snp_frequencies
INPUT="/scratch/abarbeira3/test/gtex_v8_eur_shapeit2_phased_maf01.txt.gz"
SNP_FREQUENCIES="/scratch/abarbeira3/test/gtex_v8_eur_shapeit2_phased_maf01_snps.txt.gz"
OUTPUT_FOLDER="."
OUTPUT_PREFIX="gtex_v8_eur_shapeit2_phased_maf01_qdimputed_f"

import numpy
def qdi(comps):
    d = comps[1:]
    mean = numpy.mean(numpy.array([x for x in d if x != "NA"], dtype=numpy.float32))/2
    mean=str(mean)
    comps = [comps[0]] + [str(x) if x != "NA" else mean for x in d]
    return "{}\n".format("\t".join(comps))

def run():
    if not os.path.exists(OUTPUT_FOLDER): os.makedirs(OUTPUT_FOLDER)

    print("reading gtex snp frequencies")
    snp_frequencies = read_snp_frequencies(SNP_FREQUENCIES)

    print("processing geno")
    last_chr = None
    chromosomes = {str(i) for i in range(1, 23)}

    o = None
    with gzip.open(INPUT) as f:
        header = f.readline().decode()
        for i,line in enumerate(f):
            line = line.decode()
            comps = line.strip().split()

            variant = comps[0]
            chr_ = variant.split("_")[0].split("chr")[1]
            if not chr_ in chromosomes: continue

            if not variant in snp_frequencies: continue

            if chr_ != last_chr:
                print(chr_)
                if o: o.close()
                last_chr = chr_
                o = gzip.open(os.path.join(OUTPUT_FOLDER, OUTPUT_PREFIX + ".chr{}.txt.gz".format(chr_)), "w")
                o.write(header.encode())

            line = qdi(comps)
            o.write(line.encode())
        o.close()

if __name__ == "__main__":
    run()