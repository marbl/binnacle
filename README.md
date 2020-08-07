# Using scaffolds to improve the contiguity and quality of metagenomic bins

Metagenomics has revolutionized the field of microbiology, however, reconstructing complete genomes of organisms from metagenomic data is still challenging. Recovered genomes are often fragments, due to repeats within and across genomes, uneven abundance of organisms, sequencing errors, and strain-level variations within a single sample. To address the fragmented nature of metagenomic assemblies, scientists rely on a process called binning which clusters together contigs that are inferred to originate from the same organism. Existing binning algorithms use oligonucleotide frequencies and contig abundance (coverage) within and across samples to group together contigs from the same organism. However, these algorithms often miss short contigs and contigs from regions with unusual coverage or DNA composition characteristics, such as mobile elements. We propose that information from assembly graphs can assist current strategies for metagenomic binning. We use [MetaCarvel](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1791-3), a metagenomic scaffolding tool, to construct assembly graphs where contigs are nodes and edges are inferred based on mate pair or paired-end reads. We show that binning scaffolds, rather than contigs, improves the contiguity and quality of the resulting bins on twenty stool samples from the Human Microbiome Project.

## Results
For 20 complex stool samples from the human microbiome project, we ran MetaBat with contigs and MetaBat with our scaffolds. For all samples, there was an increase in the number of contigs and length of the overall bins when using scaffolds instead of contigs for binning. We also used Checkm to evaluate completeness and contamination of bins. 

<img src="img/cumulative_len.png" width=400 /> <img src="img/contamination2completeness.png" width=300 />


## Software Requirements:

<ol>
  <li> Python 3(3.7.3) </li>
  <li> Bedtools (>=2.2.3) </li>
</ol>

## Python Package Requirements
<ol>
  <li> Numpy (>=1.15.4) </li>
  <li> Pandas (>=0.23.4) </li>
  <li> Networkx (>=2.2) </li>
  <li> matplotlib (>=2.2.3) </li>
  <li> Bio-python (>=1.73) </li>
</ol>  

## Running binnacle

Before running binnacle, we run genomecov(-d) on the bed file obtained by mapping reads onto contigs. We also run metacarvel with repeat detection(-r) enabled
(https://github.com/marbl/MetaCarvel/wiki) and the assembly graph output is used to estimate coverages. 

To run binnacle, call the python wrapper script **src/Calculate_Coverages.py** with the following parameters, 
<ol>
  <li> Path to the metacarvel assembly graph</li>
  <li> Path to the output from running genomecov(-d option enabled) on the bed file obtained from mapping reads on to the contigs. </li>
  <li> Output directory. The location where you want the outputs to be written to. </li>
  <li> Path to the fasta file containing the contigs obtained as a result of assembly. </i>
</ol> 

The bash script **src/Calculate_Coverages.sh** gives a typical example of running an array job for 20 samples in the SLURM cluster. 
