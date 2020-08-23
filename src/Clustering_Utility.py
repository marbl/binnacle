import numpy as np
import pandas as pd
from os import listdir
from Binnacle_IO_Utility import *
from Compute_Scaffold_Coverages_Utility import *

def Process_Scaffold_Coverages(df_coverages, df_coords):
	scaffolds = list(df_coords.index)
	counter = pd.DataFrame(data = {'Scaffold':scaffolds})
	counter['Counter'] = 1
	counter = counter.groupby('Scaffold').sum()
	counter = dict(zip(list(counter.index), list(counter['Counter'])))
	df_coords['Coords'] = list(zip(df_coords['Start'], df_coords['End']))
	mu_list, sigma_list = [], []

	for scaffold in np.unique(scaffolds):
		df_coords_filtered = df_coords.loc[scaffold]
		if counter[scaffold] > 1:
			coords = dict(zip(df_coords_filtered['Contig'], df_coords_filtered['Coords']))
		else: 
			coords = {df_coords_filtered['Contig']:df_coords_filtered['Coords']}
		contigs = list(coords.keys())
		df_coverages_scaffold = df_coverages.loc[contigs]
		coverages = Compute_Coverage(df_coverages_scaffold, coords)
		mean, std = round(np.mean(coverages),1), round(np.std(coverages), 1)
		mu_list.append(mean)
		sigma_list.append(std)

	df_summary = pd.DataFrame(data = {'Scaffold_id':np.unique(scaffolds), 'Mu':mu_list, 'Sigma':sigma_list})
	df_summary = df_summary.set_index('Scaffold_id')
	return df_summary

def Prepare_Feature_Matrix(parent_coords, parent_summary, align_flag, align_dir, output_dir):
	df_coords = pd.read_csv(parent_coords, names = ['cc_aft_dlink', 'cc_bef_dlink', 'Contig', 'Start', 'End'], 
							sep = '\t', index_col = ['cc_aft_dlink'])
	df_coords['Contig'] = df_coords['Contig'].astype(str)
	node_list = df_coords['Contig'].tolist()
	df_summary = pd.read_csv(parent_summary, names = ['Scaffold_id', 'Span', 'Mu_0', 'Sigma_0'], sep='\t', index_col = 'Scaffold_id')
	if align_flag.lower() == "false":
		return df_summary
	ctr = 1
	if len(align_dir) == 0:
		print('Alignment Files Not Found')
		return df_summary

	align_files = listdir(align_dir)
	for f in align_files:
		if (f[0].isalpha() or f[0].isnumeric()):
			filepath = align_dir+f
			df_coverage = Load_Read_Coverage(filepath, node_list, output_dir) 
			print('Loaded Coverage Summary....')
			df_summary_arg = Process_Scaffold_Coverages(df_coverage, df_coords)
			df_summary_arg = df_summary_arg.rename(columns={'Mu':'Mu_'+str(ctr), 'Sigma':'Sigma_'+str(ctr)})
			df_summary = df_summary.join(df_summary_arg)
			ctr += 1
	return df_summary

def Format_Outputs(binning_method, coords_path, summary_path, output_dir, align_flag, align_dir):
    df_summary = Prepare_Feature_Matrix(coords_path, summary_path, align_flag, align_dir, output_dir)
    if binning_method.lower().startswith("metabat"):
        df_summary['Mean'] = df_summary['Mu_0']
        num_samples = int((len(df_summary.columns.tolist())-2)/2)
        columns = ['Span','Mean']
        for n in range(num_samples):
            columns.append('Mu_'+str(n))
            columns.append('Sigma_'+str(n))
        df_summary = df_summary[columns]
        df_summary.to_csv(output_dir+'Feature_Matrix_'+binning_method.lower()+'.txt', sep = '\t') 
    elif binning_method.lower().startswith('maxbin') or binning_method.lower().startswith('concoct'):
        del df_summary['Span']
        num_samples = int(len(df_summary.columns.tolist())/2)
        columns = []
        for n in range(num_samples):
            columns.append('Mu_'+str(n))
        df_summary = df_summary[columns]
        df_summary.to_csv(output_dir+'Feature_Matrix_'+binning_method.lower()+'.txt', sep = '\t',header = False) 
    else:
        del df_summary['Span']
        num_samples = int(len(df_summary.columns.tolist())/2)
        columns = []
        for n in range(num_samples):
            columns.append('Mu_'+str(n))
            columns.append('Sigma_'+str(n))
        df_summary = df_summary[columns]
        df_summary.to_csv(output_dir+'Feature_Matrix_'+binning_method.lower()+'.txt',sep='\t',header=False) 
    