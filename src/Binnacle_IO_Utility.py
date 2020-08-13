import io
from Bio import SeqIO
from os.path import isfile
from Compute_Scaffold_Coverages_Utility import *

# <h2> Computing Coverages </h2>
#     
# To compute the read coverages we use the *genomecov* program, part of the *bedtools* suite. 
# We run *genomecov* with *-d* optional enabled so that we get the per base depth. 
# The output of the prgram is utilized to compute the depth along the scaffold and this function loads the output of genomecov as a dataframe.   



def Load_Read_Coverage_and_Assembly_Graph(graphpath, covpath):
    G = nx.read_gml(graphpath)
    nodes = list(G.nodes())
    df_coverage = pd.DataFrame()
    
    df_coverage_chunks = pd.read_csv(covpath,names = ['Contig','Loc','coverage'], 
                              sep = '\t', low_memory = False, memory_map = True, 
                              dtype = {'Contig': str, 'Loc': 'int32', 'coverage': 'int32'},
                              engine='c', chunksize = 5000000, index_col = ['Contig'])
    for chunk in df_coverage_chunks:
        chunk.index = chunk.index.astype('str')
        nodes_pres = list(set(chunk.index.tolist()).intersection(set(nodes)))
        temp = chunk.loc[nodes_pres]
        temp = temp.reset_index()
        df_coverage = pd.concat([df_coverage, temp])

    df_coverage['Loc'] = df_coverage['Loc']-1
    df_coverage = df_coverage.sort_values(by = ['Contig', 'Loc'])
    df_coverage = df_coverage.set_index('Contig')
    print(df_coverage.info())
    return df_coverage, G



def Write_Coverage_Outputs(graph,df_coverage, outdir):
    if not isdir(outdir):
        mkdir(outdir)
        
    weakly_connected_components = nx.weakly_connected_components(graph)
    
    coverage_before_delinking = io.FileIO(outdir +'Coverages_Before_Delinking.txt', 'w')
    coords_before_delinking = io.FileIO(outdir + 'Coords_Before_Delinking.txt', 'w')
    summary_before_delinking = io.FileIO(outdir + 'Summary_Before_Delinking.txt', 'w') 
    wb_cov_before_delinking = io.BufferedWriter(coverage_before_delinking)
    wb_coords_before_delinking = io.BufferedWriter(coords_before_delinking)
    wb_summary_before_delinking = io.BufferedWriter(summary_before_delinking)
    
    coverage_after_delinking = io.FileIO(outdir + 'Coverages_After_Delinking.txt', 'w')
    coords_after_delinking = io.FileIO(outdir + 'Coords_After_Delinking.txt', 'w')
    summary_after_delinking = io.FileIO(outdir + 'Summary_After_Delinking.txt', 'w') 
    wb_cov_after_delinking = io.BufferedWriter(coverage_after_delinking)
    wb_coords_after_delinking = io.BufferedWriter(coords_after_delinking)
    wb_summary_after_delinking = io.BufferedWriter(summary_after_delinking)
    
    cc_before_delinking, cc_after_delinking = 0, 0

    for conn in weakly_connected_components:
        test = nx.DiGraph(graph.subgraph(conn))
        nodes = list(test.nodes())
        
        if len(nodes) > 1:
            min_node, min_indegree = Return_Starting_Point(test)
            if min_indegree > 0: 
                print('Requires graph simplification')
                test = Random_Simplify(test, min_node)
                min_node, min_indegree = Return_Starting_Point(test)
        else: min_node = nodes[0]

        cc_before_delinking += 1
        df_coverage_cc = df_coverage.loc[nodes]
        coverage, coords = Compute_Coverage(test, df_coverage_cc, min_node)
        
            
        flag = False

        if len(nodes) == 1:
            for i in range(len(coverage)):
                d = bytes(str(cc_before_delinking)+'\t'+str(i)+'\t'+str(coverage[i])+'\t0\n', encoding = 'utf-8')
                wb_cov_before_delinking.write(d)  
            for c in coords:
                d = bytes(str(cc_before_delinking)+'\t'+c+'\t'+ str(coords[c][0]) + '\t' +  str(coords[c][1]) + '\n', encoding = 'utf-8')
                wb_coords_before_delinking.write(d)
            cc_after_delinking += 1
            flag =  True

        if len(nodes) > 1:
            mean_ratios = Helper_Changepoints_Z_Stat(deepcopy(coverage))
            for i in range(len(coverage)):
                d = bytes(str(cc_before_delinking)+'\t'+str(i)+'\t'+str(coverage[i])+'\t'+str(mean_ratios[i])+'\n', encoding = 'utf-8')
                wb_cov_before_delinking.write(d)  
            for c in coords:
                d = bytes(str(cc_before_delinking)+'\t'+c+'\t'+ str(coords[c][0]) + '\t' +  str(coords[c][1]) + '\n', encoding = 'utf-8')
                wb_coords_before_delinking.write(d)

            outliers = ID_outliers(mean_ratios, 99)
            outliers = Filter_Neighbors(outliers, mean_ratios)
            Pos_Dict = Return_Contig_Scaffold_Positions(coords)
            g_removed = Get_Outlier_Contigs(outliers, Pos_Dict, coords, test, 100)
            
            mu, dev = round(np.mean(coverage),1), round(np.std(coverage),1)
            d_before_dlink = bytes(str(cc_before_delinking) + '\t' + str(mu) + '\t' + str(dev) + '\n', encoding = 'utf-8')
            wb_summary_before_delinking.write(d_before_dlink)

            delinked_conn_comps = list(nx.weakly_connected_component_subgraphs(g_removed))
            print('Debug---->', cc_before_delinking, len(nodes), len(delinked_conn_comps), len(coverage), len(mean_ratios))

            if len(delinked_conn_comps) == 1:
                cc_after_delinking += 1
                flag = True  
            else:
                for cc in delinked_conn_comps:
                    cc_after_delinking += 1
                    nodes_cc = list(cc.nodes())

                    if len(nodes_cc) > 1:
                        min_node, min_indegree = Return_Starting_Point(cc)
                        if min_indegree > 0: 
                            cc = Random_Simplify(cc, min_node)
                            min_node, min_indegree = Return_Starting_Point(cc)
                    else: min_node = nodes_cc[0]
                    coverage_cc, coords_cc = Compute_Coverage(cc, df_coverage_cc, min_node)
                    mu, dev = round(np.mean(coverage_cc),1), round(np.std(coverage_cc),1)
                    d_after_dlink = bytes(str(cc_after_delinking)+'\t'+str(mu)+'\t'+str(dev)+'\n', encoding = 'utf-8')
                    wb_summary_after_delinking.write(d_after_dlink)
                    
                    for i in range(len(coverage_cc)):
                        d = bytes(str(cc_after_delinking)+'\t'+str(cc_before_delinking)+'\t'+str(i)+'\t'+str(coverage_cc[i])+'\n', encoding = 'utf-8')
                        wb_cov_after_delinking.write(d)      
                    for c in coords_cc:
                        d = bytes(str(cc_after_delinking)+'\t'+str(cc_before_delinking)+'\t'+c+'\t'+str(coords_cc[c][0])+'\t'+str(coords_cc[c][1])+'\n',encoding = 'utf-8')
                        wb_coords_after_delinking.write(d)

        if (flag):
            mu, dev = round(np.mean(coverage),1), round(np.std(coverage),1)
            d_before_dlink = bytes(str(cc_before_delinking) + '\t'+  str(mu) +'\t'+ str(dev) + '\n', encoding = 'utf-8')
            d_after_dlink = bytes(str(cc_after_delinking) + '\t'+  str(mu) +'\t'+ str(dev) + '\n', encoding = 'utf-8')
            wb_summary_before_delinking.write(d_before_dlink)
            wb_summary_after_delinking.write(d_after_dlink)
            
            for i in range(len(coverage)):
                d = bytes(str(cc_after_delinking)+'\t'+str(cc_before_delinking)+'\t'+str(i)+'\t'+str(coverage[i])+'\n', encoding = 'utf-8')
                wb_cov_after_delinking.write(d)
                    
            for c in coords:
                d = bytes(str(cc_after_delinking) + '\t' + str(cc_before_delinking) + '\t' +c+'\t'+str(coords[c][0]) + '\t' + str(coords[c][1]) + '\n', encoding = 'utf-8')
                wb_coords_after_delinking.write(d)

    del df_coverage
    wb_cov_before_delinking.flush()
    wb_coords_before_delinking.flush()
    wb_summary_before_delinking.flush()
    wb_cov_after_delinking.flush()
    wb_coords_after_delinking.flush()
    wb_summary_after_delinking.flush()
    
    print('Done.....')


def Load_FASTA_File(input_file):
    d = {}
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    for f in fasta_sequences:
        d[f.name] = str(f.seq)
    return d

def Get_Contigs_in_Scaffolds(input_file):
    df_coords = pd.read_csv(input_file, names = ['CC_after_dlnk', 'CC_before_dlnk', 
                                                  'Contig', 'Start', 'End'], sep = '\t')
    df_coords = df_coords[['CC_after_dlnk','Contig']]
    df_coords = df_coords.groupby('CC_after_dlnk')['Contig'].apply(list)
    return (df_coords.to_dict())

def Write_Scaffolds(Contigs_Path, Coords_Path, op_path):
    try:
        contigs = Load_FASTA_File(Contigs_Path)
        Scaffolds = Get_Contigs_in_Scaffolds(Coords_Path)
        connected_component_keys = list(Scaffolds.keys())
        f_op = io.FileIO(op_path, 'w')
        wb = io.BufferedWriter(f_op)

        for c in connected_component_keys:
            fasta_seq = '>'+str(c)+'\n'
            contigs_in_scaffold = list(Scaffolds[c])
            add_buff = 'N'*100
            for contig in contigs_in_scaffold[:-1]:
                fasta_seq += contigs[contig] + add_buff
            fasta_seq += contigs[contigs_in_scaffold[-1]]+'\n'
            wb.write(bytes(fasta_seq, encoding = 'utf-8'))
        wb.flush()
    except FileNotFoundError:
        print('Check Filepaths. File not Found')

def Format_Outputs(binning_method, coverage_path, outdir):
    df = pd.read_csv(coverage_path, sep='\t', names = ['Scaffold_id','Mean','Deviation'], index_col = ['Scaffold_id'])
    if (binning_method.lower().startswith("metabat")):
        df.to_csv(outdir+'Coverage_Metabat.txt', sep = '\t')
    elif (binning_method.lower().startwith("maxbin")):
        del df['Deviation']
        df.to_csv(outdir+'Coverage_Maxbin2.txt', sep = '\t', header = False)
    elif (binning_method.lower().startwith("concoct")) or (binning_method == 'CONCOCT'):
        del df['Deviation']
        df.to_csv(outdir+'Coverage_Concoct.txt', sep = '\t', header = False)
    else:
        ## Need to add Binnacle clustering
        pass
       

