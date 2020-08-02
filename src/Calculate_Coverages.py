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

s = sys.argv[1]
sample = s.replace(".txt","")

Graph_path= graphspath+sample+'/'+sample+'_scaffolds/oriented.gml'
df_coverage, G =Load_Read_Coverage_and_Assembly_Graph(Graph_path, covpath+s) 

print('Loaded Coverage and Assembly Graph')
(df_coords_before_delinking, 
 df_coords_after_delinking,
 df_coverage_scaffolds_before_delinking, 
 df_coverage_scaffolds_after_delinking) = Write_Coverage_Outputs(G, df_coverage)

df_cov_scaffolds_before_delinking_mean = df_coverage_scaffolds_before_delinking.groupby('Connected_Component_Before_Delinking').mean()[['Coverage']]
df_cov_scaffolds_before_delinking_mean = df_cov_scaffolds_before_delinking_mean.rename(columns = {'Coverage':'Cov_Mean'})
df_cov_scaffolds_before_delinking_dev = df_coverage_scaffolds_before_delinking.groupby('Connected_Component_Before_Delinking').std()[['Coverage']]
df_cov_scaffolds_before_delinking_dev = df_cov_scaffolds_before_delinking_dev.rename(columns = {'Coverage':'Cov_dev'})
df_cov_scaffolds_before_summ = df_cov_scaffolds_before_delinking_mean.join(df_cov_scaffolds_before_delinking_dev)

df_cov_scaffolds_after_delinking_mean = df_coverage_scaffolds_after_delinking.groupby('Connected_Component_After_Delinking').mean()[['Coverage']]
df_cov_scaffolds_after_delinking_mean = df_cov_scaffolds_after_delinking_mean.rename(columns = {'Coverage':'Cov_Mean'})
df_cov_scaffolds_after_delinking_dev = df_coverage_scaffolds_after_delinking.groupby('Connected_Component_After_Delinking').std()[['Coverage']]
df_cov_scaffolds_after_delinking_dev = df_cov_scaffolds_after_delinking_dev.rename(columns = {'Coverage':'Cov_dev'})
df_cov_scaffolds_after_summ = df_cov_scaffolds_after_delinking_mean.join(df_cov_scaffolds_after_delinking_dev)

df_coords_before_delinking.to_csv(op_path_before_delinking + sample+'_Coords.csv')
df_coverage_scaffolds_before_delinking.to_csv(op_path_before_delinking + sample+'_Coverages.csv')
df_cov_scaffolds_before_summ.to_csv(op_path_before_delinking+sample+'_Summary.csv')

df_coords_after_delinking.to_csv(op_path_after_delinking + sample+'_Coords.csv')
df_coverage_scaffolds_after_delinking.to_csv(op_path_after_delinking + sample+'_Coverages.csv')
df_cov_scaffolds_after_summ.to_csv(op_path_after_delinking+sample+'_Summary.csv')
