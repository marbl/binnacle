#!/usr/bin/env python
# coding: utf-8

'''
Program developed at Pop lab at the CBCB, University of Maryland by 
Harihara Subrahmaniam Muralidharan, Nidhi Shah, Jacquelyn S Meisel. 
'''

from Binnacle_IO_Utility import *

def Process_Scaffold_Coverages(df_coverages, df_coords):
    '''
    Function to generate the features for all vs all alignments. 
    Input:
        df_coverages: Coverage of the contigs obtained by mapping contigs to reads.
        df_coords: Coordinates for the contigs along the scaffolds along the global coorinates. 
    Output:
        df_summary: Dataframe object with the mean and deviation for the scaffold. 
    '''
    scaffolds = list(df_coords.index)
    counter = pd.DataFrame(data = {'Scaffold':scaffolds})
    counter['Counter'] = 1
    counter = counter.groupby('Scaffold').sum()
    counter = dict(zip(list(counter.index), list(counter['Counter'])))
    df_coords['Coords'] = list(zip(df_coords['Start'], df_coords['End']))
    mu_list, sigma_list, spanlist = [], [], []

    for scaffold in np.unique(scaffolds):
        df_coords_filtered = df_coords.loc[scaffold]
        if counter[scaffold] > 1:
            coords = dict(zip(df_coords_filtered['Contig'], df_coords_filtered['Coords']))
            max_coords = max(df_coords_filtered['Coords'].max())
        else: 
            coords = {df_coords_filtered['Contig']:df_coords_filtered['Coords']}
            max_coords  = df_coords_filtered[['Start','End']].max()

        contigs = list(coords.keys())
        df_coverages_scaffold = df_coverages.loc[contigs]
        coverages = Compute_Coverage(df_coverages_scaffold, coords)
        mean, std = round(np.mean(coverages),1), round(np.std(coverages), 1)
        mu_list.append(mean)
        sigma_list.append(std)
        spanlist.append(max_coords)

    df_summary = pd.DataFrame(data = {'Scaffold_id':np.unique(scaffolds), 'Span':spanlist, 
                                      'Mu':mu_list, 'Sigma':sigma_list})
    df_summary = df_summary.set_index('Scaffold_id')
    print(df_summary.head())
    return df_summary

'''def Prepare_Feature_Matrix(parent_coords, parent_summary, align_flag, align_dir, output_dir):
    Function to prepare all vs all alignments feature matrix for various binning algorithm
    Input:
        Parent_coords: Coordinate file for the scaffolds
        parent_summary: The path to the txt file containing the summary information by running binnacle on the original sample.
        align_flag: A flag that specifes whether to peform all vs all alignments
        align_dir: Directory that contains the output of running genomecov -d on the bed files obtaining mapping reads of all samples to contigs of all samples
        output_dir: Location where the putputs of binnacle is written to. 
    Output:
        df_summary: A dataframe containg the feature matrix 
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
    return df_summary'''

def Format_Outputs(summary_dir, binning_method):
    '''
    Function to format outputs to the binning method specified. 
    Input:
        binnning_method: Choice of binning method (metabat, concoct, maxbin)
        coords_path: Coordinate file for the scaffolds
        summary_path: The path to the txt file containing the summary information by running binnacle on the original sample.
        align_flag: A flag that specifes whether to peform all vs all alignments
        align_dir: Directory that contains the output of running genomecov -d on the bed files obtaining mapping reads of all samples to contigs of all samples
        output_dir: Location where the putputs of binnacle is written to.
    '''
    files = listdir(summary_dir)
    files.sort()
    df_summary = pd.DataFrame()
    ctr = 0
    for i in range(0, len(files)):
        if 'Summary.txt' in files[i]:
            print(files[i])
            col_prefix = files[i].replace("_Summary.txt","")
            df = pd.read_csv(summary_dir+files[i], names = ['Scaffold','Span', col_prefix+'_Mu', col_prefix+'_Sigma'],
                             sep='\t', index_col = 'Scaffold')
            print(df.head())
            if ctr > 0:
                del df['Span']
            if ctr == 0:
                df['Avg_Depth'] = 0
                df = df[['Span','Avg_Depth',col_prefix+'_Mu',col_prefix+'_Sigma']]
            ctr += 1
            df_summary = df_summary.join(df, how = 'outer')
    df_mu = df_summary.filter(regex='_Mu')
    if binning_method.lower().startswith('metabat'):
        df_summary['Avg_Depth'] = df_mu.mean(axis=1)
    elif binning_method.lower().startswith('maxbin') or binning_method.lower().startswith('concoct'):
        df_summary = df_mu
    else:
        del df_summary['Span'], df_summary['Avg_Depth']
    return df_summary