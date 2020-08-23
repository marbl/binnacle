#!/bin/bash
#SBATCH -J Calculate_Coverages # Job name
#SBATCH -o Calculate_Coverages.o%j # Name of Output File
#SBATCH -e Calculate_Coverages.e%j # Name of Error File
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for Job Info
#SBATCH --mail-type=all
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=2-5

source /cbcbhomes/hsmurali/binnacle/bin/activate
sample=`head -n ${SLURM_ARRAY_TASK_ID} samples_2.txt | tail -n 1`
echo ${sample}

graphpath=/fs/cbcb-lab/mpop/MetaCarvel_paper/hmp_scaffolds/stool/${sample}/${sample}_scaffolds/oriented.gml
coverage_path=/fs/cbcb-scratch/hsmurali/binnacle/all_vs_all_alignments/genomecov_d/${sample}_${sample}.txt
output_path=/fs/cbcb-scratch/hsmurali/binnacle/all_vs_all_alignments/Binnacle_Scaffold_Coverages/${sample}/
contigs_path=/fs/cbcb-lab/mpop/MetaCarvel_paper/hmp_scaffolds/stool/${sample}/${sample}.fna
alignpath=/fs/cbcb-scratch/hsmurali/binnacle/all_vs_all_alignments/genomecov_d/${sample}/

time python Calculate_Coverages.py -g ${graphpath} -a ${coverage_path} -d ${output_path} -c ${contigs_path} -b binnacle -A true -M ${alignpath}

