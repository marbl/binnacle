#!/bin/bash
#SBATCH -J Calculate_Coverages # Job name
#SBATCH -o Calculate_Coverages.o%j # Name of Output File
#SBATCH -e Calculate_Coverages.e%j # Name of Error File
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for Job Info
#SBATCH --mail-type=all
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=2-20

source activate hsmurali27
file=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -n 1`
python Calculate_Coverages.py ${file}