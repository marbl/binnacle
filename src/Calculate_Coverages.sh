#!/bin/bash
#SBATCH -J Sharon-Data-Binnacle # Job name
#SBATCH -o Sharon-Data-Binnacle.o%j # Name of Output File
#SBATCH -e Sharon-Data-Binnacle.e%j # Name of Error File
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for Job Info
#SBATCH --mail-type=all
#SBATCH --time=0-5:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=1-24

source /cbcbhomes/hsmurali/binnacle/bin/activate

abundance_path=/fs/cbcb-lab/mpop/projects/refining_bins_with_assembly_graphs/HMP_Samples_Analysis/all_vs_all_alignments/Mapping_Bed_Files/
metacarvel_path=/fs/cbcb-lab/mpop/MetaCarvel_paper/hmp_scaffolds/stool/

ls ${abundance_path} > samples.txt
s=`head -n ${SLURM_ARRAY_TASK_ID} Samples_Re_Run.txt | tail -n 1`
sample=${s:0:9}

echo ${s}
echo ${sample}

graphpath=${metacarvel_path}${sample}/${sample}_scaffolds/oriented.gml
contigpath=${metacarvel_path}${sample}/${sample}.fna
coverage_path=${abundance_path}/${s}
output_path=/fs/cbcb-scratch/hsmurali/binnacle/Binnacle_Scaffold_Coverages_Corrected_Span/${sample}/
coords_path=${output_path}Coords_After_Delinking.txt

time python /fs/cbcb-scratch/hsmurali/binnacle/src/Estimate_Abundances.py -o ${coords_path} -a ${coverage_path}  -d ${output_path} 

#time python /fs/cbcb-scratch/hsmurali/binnacle/src/Estimate_Abundances.py -g ${graphpath} -a ${coverage_path} -c ${contigpath} -d ${output_path}