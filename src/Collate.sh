#!/bin/bash
#SBATCH -J Feature_Mat # Job name
#SBATCH -o Feature_Mat.o%j # Name of Output File
#SBATCH -e Feature_Mat.e%j # Name of Error File
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for Job Info
#SBATCH --mail-type=all
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=1-18

source /cbcbhomes/hsmurali/binnacle/bin/activate
op_dir_path=/fs/cbcb-scratch/hsmurali/binnacle/Sharon-Data/Summarize_Linear_Scaffold_Coverage/
ls -D ${op_dir_path} > samples.txt
s=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -n 1`
echo ${s}

output_path=${op_dir_path}${s}/

time python Collate.py -d ${output_path} -m metabat 