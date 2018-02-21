#!/bin/bash
#PBS -N job_gtex_v8_split_geno
#PBS -S /bin/bash
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -o ../joblogs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e ../logs/${PBS_JOBNAME}.e${PBS_JOBID}.err

cd $PBS_O_WORKDIR

module load gcc/6.2.0
module load python/3.5.3

python3 split_geno_by_chromosome.py