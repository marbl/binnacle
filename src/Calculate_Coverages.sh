#!/bin/bash
#SBATCH -J Calculate_Coverages # Job name
#SBATCH -o Calculate_Coverages.o%j # Name of Output File
#SBATCH -e Calculate_Coverages.e%j # Name of Error File
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for Job Info
#SBATCH --mail-type=all
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=1-3

source /cbcbhomes/hsmurali/binnacle/bin/activate
abundance_path=/fs/cbcb-lab/mpop/projects/refining_bins_with_assembly_graphs/all_vs_all_alignments/Mapping_Bed_Files/
ls ${abundance_path} > samples.txt
s=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -n 1`
sample=${s:0:9}
echo ${s}
echo ${sample}
coverage_path=${abundance_path}${s}
output_path=/fs/cbcb-scratch/hsmurali/binnacle/Binnacle_Scaffold_Coverages/${sample}/
coords_path=${output_path}Coords_After_Delinking.txt

time python Estimate_Abundances.py -o ${coords_path} -a ${coverage_path} -d ${output_path} 