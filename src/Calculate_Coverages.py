#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from Binnacle_IO_Utility import *


# In[2]:

graph_path = sys.argv[1]
coverage_path  =  sys.argv[2]
output_dir = sys.argv[3]
Contigs_Path = sys.argv[4]


if not isdir(output_dir): mkdir(output_dir)
df_coverage, G =Load_Read_Coverage_and_Assembly_Graph(graph_path, coverage_path) 
print('Loaded Coverage and Assembly Graph')
Write_Coverage_Outputs(G, df_coverage, output_dir)
Coords_Path = output_dir+'Coords_After_Delinking.txt'
op_path = output_dir+'Scaffolds.fasta'
Write_Scaffolds(Contigs_Path, Coords_Path, op_path)
print('Written Fasta Files')