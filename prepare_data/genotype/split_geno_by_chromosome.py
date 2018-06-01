#!/usr/bin/env python
__author__ = "alvaro barbeira"
import os
import gzip
from convert_gtex_annotation_to_pred import read_snp_frequencies
from timeit import default_timer as timer

INPUT="/group/im-lab/nas40t2/abarbeira/projects/gtex_v8/data_formatting/vcf_process/gtex_v8_eur_filtered.txt.gz"
#SNP_FREQUENCIES="/scratch/abarbeira3/test/gtex_v8_eur_shapeit2_phased_maf01_snps.txt.gz"
OUTPUT_FOLDER="."
OUTPUT_PREFIX="gtex_v8_eur_shapeit2_phased_maf01_qdimputed"
MAF_FILTER=0.01

import numpy
def qdi(comps):
    d = comps[1:]
    _mean = float(numpy.mean(numpy.array([x for x in d if x != "NA"], dtype=numpy.float32)))
    mean=str(_mean)
    comps = [comps[0]] + [x if x != "NA" else mean for x in d]
    return "{}\n".format("\t".join(comps)), _mean

def run():
    if not os.path.exists(OUTPUT_FOLDER): os.makedirs(OUTPUT_FOLDER)

    #print("reading gtex snp frequencies")
    #snp_frequencies = read_snp_frequencies(SNP_FREQUENCIES)

    print("processing geno")
    last_chr = None
    chromosomes = {str(i) for i in range(1, 23)}

    o = None
    time = None
    with gzip.open(INPUT) as f:
        header = f.readline().decode()
        for i,line in enumerate(f):
            line = line.decode()
            comps = line.strip().split()

            variant = comps[0]
            chr_ = variant.split("_")[0].split("chr")[1]
            if not chr_ in chromosomes: continue

            #if not variant in snp_frequencies: continue

            if chr_ != last_chr:
                if not time:
                    time = timer()
                else:
                    _time = timer()
                    print(str(_time-time))
                    time=_time
                print(chr_)
                if o: o.close()
                last_chr = chr_
                o = gzip.open(os.path.join(OUTPUT_FOLDER, OUTPUT_PREFIX + "_maf{}_chr{}.txt.gz".format(MAF_FILTER,chr_)), "w")
                o.write(header.encode())

            line, mean = qdi(comps)

            if (mean/2<MAF_FILTER or 1-MAF_FILTER < mean/2): continue

            o.write(line.encode())
        o.close()

if __name__ == "__main__":
    run()