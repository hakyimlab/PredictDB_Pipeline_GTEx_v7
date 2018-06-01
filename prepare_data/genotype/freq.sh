#!/bin/bash
#PBS -N job_gtex_v8_frequency
#PBS -S /bin/bash
#PBS -l walltime=8:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err

INPUT="/group/im-lab/nas40t2/abarbeira/projects/gtex_v8/data_formatting/vcf_process/gtex_v8_eur_filtered.txt.gz"
OUTPUT="gtex_v8_freq.txt.gz"

module load gcc/6.2.0
module load python/3.5.3

cd $PBS_O_WORKDIR

python3 -c '
import sys
import gzip

with gzip.open(sys.argv[1]) as _i:
    _i.readline()
    with gzip.open(sys.argv[2], "w") as _o:
        _o.write("varID\tfrequency\n".encode())
        for i,l in enumerate(_i):
            comps = l.decode().strip().split()
            var = comps[0]
            d = [float(x) for x in comps[1:] if x != "NA"]
            freq = str(sum(d)/(2*len(d))) if len(d) > 0 else "NA"
            _o.write("{}\t{}\n".format(var,freq).encode())
' $INPUT $OUTPUT
