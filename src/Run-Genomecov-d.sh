#!/bin/bash

module load bedtools
sample=${1}
alignment_path=/fs/cbcb-lab/mpop/MetaCarvel_paper/hmp_scaffolds/stool/${sample}/${sample}_scaffolds/alignmed_sorted.bed
genome_path=/fs/cbcb-lab/mpop/MetaCarvel_paper/hmp_scaffolds/${sample}/${sample}_scaffolds/contig_length
op_path=/fs/cbcb-lab/mpop/projects/refining_bins_with_assembly_graphs/genomecov_d/${sample}.txt

genomecov -d -i ${alignmed_path} -g ${genome_path} > op_path