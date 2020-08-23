#!/bin/bash
#SBATCH -J GenomeCov # Job name
#SBATCH -o GenomeCov.o%j # Name of Output File
#SBATCH -e GenomeCov.e%j # Name of Error File
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for Job Info
#SBATCH --mail-type=all
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=2-20

file=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -n 1`
bash Run-Genomecov-d.sh ${file}