#!/bin/bash
#SBATCH -J Sharon-Data-Binnacle # Job name
#SBATCH -o Sharon-Data-Binnacle.o%j # Name of Output File
#SBATCH -e Sharon-Data-Binnacle.e%j # Name of Error File
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for Job Info
#SBATCH --mail-type=all
#SBATCH --time=0-18:00:00
#SBATCH --qos=large
#SBATCH --mem=128gb
#SBATCH --array=4-5

source /cbcbhomes/hsmurali/binnacle/bin/activate

abundance_path=/fs/cbcb-scratch/hsmurali/binnacle/Sharon-Data/all_vs_all_alignments/

ls ${abundance_path} > samples.txt
s=`head -n ${SLURM_ARRAY_TASK_ID} samples.txt | tail -n 1`
sample=${s:0:9}
echo ${s}
echo ${sample}

#graphpath=${abundance_path}${sample}/${sample}_scaffolds/oriented.gml
#contigpath=${abundance_path}${sample}/${sample}_assembly/${sample}.contigs.fa
coverage_path=${abundance_path}/${s}
output_path=/fs/cbcb-scratch/hsmurali/binnacle/Sharon-Data/Binnacle_Scaffold_Coverages/${sample}/
coords_path=${output_path}Coords_After_Delinking.txt

time python /fs/cbcb-scratch/hsmurali/binnacle/src/Estimate_Abundances.py -o ${coords_path} -a ${coverage_path}  -d ${output_path} 

#time python /fs/cbcb-scratch/hsmurali/binnacle/src/Estimate_Abundances.py -g ${graphpath} -a ${coverage_path} -c ${contigpath} -d ${output_path}