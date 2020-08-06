#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Compute_Scaffold_Coverages_Utility import *
from os.path import isdir
from os import listdir, mkdir
import sys

# In[2]:

graphspath = '/Users/harihara/Mount/MetaCarvel_paper/hmp_scaffolds/stool/'
covpath = '/Users/harihara/Research-Activities/Data/Binnacle-Op/genomecov_d/'
op_path_before_delinking = '/Users/harihara/Research-Activities/Data/Binnacle-Op/Coverages_Before_Delinking/'
op_path_after_delinking = '/Users/harihara/Research-Activities/Data/Binnacle-Op/Coverages_After_Delinking/'

if not isdir(op_path_after_delinking): mkdir(op_path_after_delinking)
if not isdir(op_path_before_delinking): mkdir(op_path_before_delinking)

sample = 'SRS104311'
Graph_path= graphspath+sample+'/'+sample+'_scaffolds/oriented.gml'
df_coverage, G =Load_Read_Coverage_and_Assembly_Graph(Graph_path, covpath+sample+'.txt') 

print('Loaded Coverage and Assembly Graph')
Write_Coverage_Outputs(G, df_coverage, '/Users/harihara/Research-Activities/Data/Binnacle-Op/', sample)