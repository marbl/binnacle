#!/bin/bash
#SBATCH -J Global-Coords # Job name
#SBATCH -o Global-Coord.o%j # Name of Output File
#SBATCH -e Global-Coord.e%j # Name of Error File
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for Job Info
#SBATCH --mail-type=all
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=1-1

source activaye hsmurali27
file=`head -n ${SLURM_ARRAY_TASK_ID} GenomeCov.txt | tail -n 1`
python Assign_Coords_BFS.py ${file}