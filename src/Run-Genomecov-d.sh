#!/bin/bash

module load bedtools

bed=${1}
contig_len=${1}
data_path=/fs/cbcb-lab/mpop/MetaCarvel_paper/hmp_scaffolds/stool/${contig_len}/${contig_len}_scaffolds/
contig_len_op=/fs/cbcb-lab/mpop/MetaCarvel_paper/hmp_scaffolds/stool/${contig_len}/${contig_len}_scaffolds/contig_length

echo ${bed}
#alignment_path=/fs/cbcb-scratch/hsmurali/binnacle/all_vs_all_alignments/Mapping_Between_Samples/Mapping_Bed_Files/${bed}.bed
alignment_path=${data_path}/alignment_sorted.bed
op_path=/fs/cbcb-scratch/hsmurali/genomecov_d/${contig_len}
#op_path=/fs/cbcb-scratch/hsmurali/binnacle/all_vs_all_alignments/genomecov_d/${bed}
#alignment_sorted=/fs/cbcb-scratch/hsmurali/binnacle/all_vs_all_alignments/Mapping_Between_Samples/Mapping_Bed_Files/${bed}_sorted.bed
#sort -k 1,1 ${alignment_path} > ${alignment_sorted}

genomeCoverageBed -d -i ${alignment_path} -g ${contig_len_op} > ${op_path}_temp.txt
LC_ALL=C sort -k1,1 ${op_path}_temp.txt > ${op_path}.txt
rm ${op_path}_temp.txt