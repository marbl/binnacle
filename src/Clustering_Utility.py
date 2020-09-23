#!/usr/bin/env python
# coding: utf-8

'''
Program developed at Pop lab at the CBCB, University of Maryland by 
Harihara Subrahmaniam Muralidharan, Nidhi Shah, Jacquelyn S Meisel. 
'''

from Binnacle_IO_Utility import *

def Process_Scaffold_Coverages(df_coverages, df_coords, df_not_found):
    '''
    Function to generate the features for all vs all alignments. 
    Input:
        df_coverages: Coverage of the contigs obtained by mapping contigs to reads.
        df_coords: Coordinates for the contigs along the scaffolds along the global coorinates. 
    Output:
        df_summary: Dataframe object with the mean and deviation for the scaffold. 
    '''
    scaffolds = df_coords['cc_aft_dlink'].tolist()
    df_coords = df_coords.set_index('cc_aft_dlink')
    counter = pd.DataFrame(data = {'Scaffold':scaffolds})
    counter['Counter'] = 1
    counter = counter.groupby('Scaffold').sum()
    counter = dict(zip(list(counter.index), list(counter['Counter'])))
    df_coords['Coords'] = list(zip(df_coords['Start'].tolist(), df_coords['End'].tolist()))
    mu_list, sigma_list, spanlist, lengthlist = [], [], [], []

    for scaffold in np.unique(scaffolds):
        df_coords_filtered = df_coords.loc[scaffold]
        if counter[scaffold] > 1:
            coords = dict(zip(df_coords_filtered['Contig'].tolist(), df_coords_filtered['Coords'].tolist()))
        else: 
            coords = {df_coords_filtered['Contig']:df_coords_filtered['Coords']}

        contigs = list(coords.keys())
        df_coverages_scaffold = df_coverages.loc[contigs]
        coverages = Compute_Coverage(df_coverages_scaffold, coords)
        mean, std = round(np.mean(coverages),1), round(np.std(coverages), 1)
        length = np.sum(np.abs(df_coords_filtered['Start'] - df_coords_filtered['End']))

        mu_list.append(mean)
        sigma_list.append(std)
        spanlist.append(len(coverages))
        lengthlist.append(length)

    df_summary = pd.DataFrame(data = {'Scaffold_id':np.unique(scaffolds), 'Length':lengthlist, 
                                      'Span':spanlist, 'Mu':mu_list, 'Sigma':sigma_list})
    df_not_found['Scaffold_id'] = list(range(1, len(df_not_found)+1))
    df_not_found['Scaffold_id'] += len(mu_list)
    df_not_found = df_not_found.rename(columns = {'Mean':'Mu', 'Std':'Sigma'})
    df_not_found['Span'] = df_not_found['Length']
    df_summary = pd.concat([df_summary, df_not_found[['Scaffold_id','Length','Span','Mu', 'Sigma']]])
    df_summary = df_summary.set_index('Scaffold_id')
    print(df_summary.head())
    return df_summary

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
        if 'Summary.txt' in files[i] and files[i][0] != '.':
            print(files[i])
            col_prefix = files[i].replace("_Summary.txt","")
            df = pd.read_csv(summary_dir+files[i], names = ['Scaffold','Length','Span', col_prefix+'_Mu', col_prefix+'_Var'],
                             sep='\t', index_col = 'Scaffold')
            df[col_prefix+'_Var'] = df[col_prefix+'_Var']*df[col_prefix+'_Var']
            print(df.head())
            if ctr > 0:
                del df['Span'], df['Length']
            if ctr == 0:
                df['Avg_Depth'] = 0
                df = df[['Span','Avg_Depth',col_prefix+'_Mu',col_prefix+'_Var']]
            ctr += 1
            df_summary = df_summary.join(df, how = 'outer')
    df_mu = df_summary.filter(regex='_Mu')
    if binning_method.lower().startswith('metabat'):
        df_summary['Avg_Depth'] = df_mu.mean(axis=1)
    elif binning_method.lower().startswith('maxbin') or binning_method.lower().startswith('concoct'):
        df_summary = df_mu
    else:
        del df_summary['Span'], df_summary['Avg_Depth']
    print(df_summary.head())
    return df_summary