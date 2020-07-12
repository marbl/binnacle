#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import networkx as nx
from os import listdir
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

rcParams = {'font.size': 17 , 'font.weight': 'normal', 'font.family': 'sans-serif',
            'axes.unicode_minus':False, 'axes.labelweight':'normal'}

graphspath = 'Mount-2/MetaCarvel_paper/hmp_scaffolds/stool/'
covpath = 'Research-Activities/Data/genomecov_d/'#'Mount-2/projects/refining_bins_with_assembly_graphs/genomecov_d/'
op_path = 'Scaffold_Coverage_After_Delinking/'


# <h2>Computing Global Coordinates</h2>
# <ol> 
#     <li>This function computes the global coordinate system for the scaffold(connected component) based on the length of the contigs in the connected component and the linking information as returned by MetaCarvel. We order the contigs based on the topological sort of the subgraph and assign coordinates in the global coorinate system in breadth first manner.</li>   
#     <li>We also consider the contig orientations and the edge orientations in accurate estimation of coordinates. </li>
#     <li>If there multiple possible assignments we pick the one that has the largest coordinate value. We choose the greedy approaach because the number of solutions grows exponentially and the problem is NP-Hard! </li>
#     <li>Finally, we normalize all other positions based on the least cooridnate value so that the coorinate system starts with 0.</li>
#     <li>The input to the program is the connected component subgraph of the assembly graph ``oriented.gml" which is a result of running MetaCarvel.</li>
# </ol>

# In[2]:


def Compute_Global_Coordinates(subgraph):
    if not nx.is_directed_acyclic_graph(subgraph):
        return []
    top_sort = list(nx.topological_sort(subgraph))
    v_0 = top_sort[0]
    v_0_orientation = subgraph.nodes[v_0]['orientation']
    v_0_length = int(subgraph.nodes[v_0]['length'])
    if v_0_orientation == 'FOW':
        start = 0
        end = start + v_0_length
    elif v_0_orientation == 'REV':
        start = 0
        end = start-v_0_length
    global_coords = {v_0 : (start,end)}
    visited = dict(zip(top_sort, [False]*len(top_sort)))
    Q = [v_0]
    g_undir = subgraph.to_undirected()
    
    while len(Q) > 0:
        src = Q.pop()
        visited[src] = True
        neighbors = g_undir.neighbors(src)
        for n in neighbors:
            if subgraph.has_edge(src,n):
                ###Estimate n's coordinates based on the conting ordering src,n
                edge = (src, n)
                edge_orientation = subgraph.edges[edge]['orientation']
                edge_overlap = int(float(g_undir.edges[edge]['mean']))
                c2_length = int(g_undir.nodes[n]['length'])
                s1, e1 = global_coords[src]
                if edge_orientation == 'EE':
                    e2 = e1 + edge_overlap
                    s2 = e2 + c2_length
                if edge_orientation == 'EB':
                    s2 = e1 + edge_overlap
                    e2 = s2 + c2_length
                if edge_orientation == 'BB':
                    s2 = s1 + edge_overlap
                    e2 = s2 + c2_length
                if edge_orientation == 'BE':
                    e2 = s1 + edge_overlap
                    s2 = e2 + c2_length
                start,end = (s2, e2)
            else:
                ###Estimate n's coordinates based on the conting ordering n,src
                edge = (n,src)
                edge_orientation = subgraph.edges[edge]['orientation']
                edge_overlap = int(float(g_undir.edges[edge]['mean']))
                c1_length = int(g_undir.nodes[n]['length'])
                s2, e2 = global_coords[src]
                if edge_orientation == 'EE':
                    e1 = e2 - edge_overlap
                    s1 = e1 - c1_length
                if edge_orientation == 'EB':
                    e1 = s2 - edge_overlap
                    s1 = e1 - c1_length
                if edge_orientation == 'BB':
                    s1 = s2 - edge_overlap
                    e1 = s1 - c1_length
                if edge_orientation == 'BE':
                    s1 = e2 - edge_overlap
                    e1 = s1 - c1_length
                start,end = (s1,e1)
            try:
                if global_coords[n][0] < start:
                    global_coords[n] = (start,end)
            except KeyError:
                global_coords[n] = (start,end)
            if visited[n] == False:
                visited[n] = True
                Q.append(n)
    ###Normalize for the global coordinate system to start from 0
    min_coord = np.inf
    for g in global_coords:
        start,end = global_coords[g]
        min_val = min(start,end)
        if min_val < min_coord:
            min_coord = min_val
    for g in global_coords:
        s,e = global_coords[g]
        global_coords[g] = s-min_coord, e-min_coord
    return global_coords


# <h2> Computing Coverages </h2>
#     
# To compute the read coverages we use the *genomecov* program, part of the *bedtools* suite. We run *genomecov* with *-d* optional enabled so that we get the per base depth. The output of the prgram is utilized to compute the depth along the scaffold and this function loads the out of genomecov as a dataframe.   

# In[3]:


def Load_Read_Coverage(sample):
    df_coverage = pd.read_csv(covpath+sample+'.txt',names = ['Contig','Loc','coverage'], 
                              sep = '\t', low_memory = False,  
                              dtype = {'Contig': str, 'Loc': np.int32, 'Coverage': np.int32},
                              index_col = 'Contig', engine='c')
    df_coverage['Loc'] = df_coverage['Loc']-1
    return df_coverage


# <h2> Computing the Depth of Coverage along the Scaffold's Global Coordinate System  </h2>
# 
# To compute the depth we run the *Compute_Global_Coordinates* and *Load_Read_Coverage* functions and pass this to the *Compute_Coverage* described below. This program estimates the per cooridnate depth by using the outputs of the said two information. The program returns the coverage along the entire scaffold and the coverage along the longest path and a dataframe contiaining the contig, its average depth and its coordinates in the scaffold.   

# In[4]:


def Compute_Coverage(connected_component, df_coverage):
    top_sort = list(nx.topological_sort(connected_component))
    coords = Compute_Global_Coordinates(connected_component)
    max_coord = -np.inf
    for c in top_sort:
        max_v = max(coords[c])
        if max_v >= max_coord: max_coord = max_v
    coverage = np.zeros(max_coord+1)
    d_arr = []
    for c in top_sort:
        contig_depth = np.array(df_coverage.loc[c]['coverage'])
        d = {'Contig':c, 'Average_Depth':np.median(contig_depth),'Coords':coords[c]}
        d_arr.append(d)
        s,e = coords[c]
        if (s > e): coverage[s:e:-1] += contig_depth[::-1]
        else: coverage[s:e] += contig_depth
    df_depths = pd.DataFrame(d_arr)
    df_depths = df_depths.set_index('Contig')
    df_depths = df_depths.loc[top_sort]
    longest_coverage = np.zeros(max_coord+1)
    longest_path = list(nx.dag_longest_path(connected_component))
    for c in longest_path:
        contig_depth = np.array(df_coverage.loc[c]['coverage'])
        s,e = coords[c]
        if (s > e): longest_coverage[s:e:-1] += contig_depth[::-1]
        else: longest_coverage[s:e] += contig_depth
    return coverage, longest_coverage, df_depths


# <h2> Change Point Detection on the Scaffold Coverage </h2>
# 
# To compute changepoints along the scaffold we slide a window of default size 1500 along the scaffold and take the ratios of the maximum of means of the window preceeding any point and the window suceeding the point and the minimum of the said quantities. Any abnromal spike magnitude indicates a change point. If the size of the contig is less than the size of the window we adaptively pick a smaller window. 

# In[5]:


def Helper_Changepoints(cov_vec, window_size = 1500):
    indices_non_zero = np.where(cov_vec > 0)
    sliced = cov_vec[indices_non_zero]
    while window_size >= len(sliced)/2:
        window_size = int(window_size/5)
        print(window_size)
    m = []
    for i in range(window_size, len(sliced)-window_size):
        pred = sliced[i-window_size:i].mean()
        succ = sliced[i:i+window_size].mean()
        if pred == 0 or succ == 0: 
            r = 0
        else: 
            r = max(pred,succ)/min(pred,succ)
        m.append(r)
        
    mean_ratios = [0]*window_size + m + [0]*window_size
    cpts = np.zeros(len(cov_vec))
    cpts[indices_non_zero] = mean_ratios
    return cpts


# <h2> Outlier Detection in Change Point Signals </h2>
# <ol>
#     <li>The following code segment is used to identify outliers in change points. To do so, we compute the peaks first.</li> 
#     <li>We don't directly identify outliers on the signal because doing so would pickup points near the peaks, these are indicators of sliding windows and not really outliers.</li> 
#     <li>To overcome the issue, we first identify peaks. A point is a peak if a point is larger than its predecessor and successor. This is performed by the function *ID_Peaks*.</li> 
#     <li>The output of this passed to *ID_outliers* which picks all those points that is gretaer than the point which is the *thresh*'s percentile. The default value of *thresh* is 98.5. </li>
#     <li>The filter outliers is still a work in progress. This is aimed at removing all the outliers points that is close to one another, but this is not super important. While the mthod described in the following block is data driven we are working on improving this method by examining the underlying graph structure.</li>
# </ol>

# In[6]:


def ID_outliers(change_point_vec, thresh):
    cutoff = np.percentile(change_point_vec, thresh)
    indices = np.where(change_point_vec >= cutoff)
    return indices

def ID_Peaks(change_point_vec, thresh = 98.5):
    indices = range(0, len(change_point_vec))
    Peak_Indices = []
    for i in range(1,len(indices)-1):
        prev_index, curr_index, next_index = indices[i-1], indices[i], indices[i+1]
        pt1, pt2, pt3 = change_point_vec[prev_index], change_point_vec[curr_index], change_point_vec[next_index]
        if max(pt1, pt2, pt3) == pt2 and pt2 > 0: Peak_Indices.append(curr_index)
    Peak_Indices = np.array(Peak_Indices)
    Peak_Values = change_point_vec[Peak_Indices]
    outlier_peaks = ID_outliers(Peak_Values, thresh)
    return Peak_Indices[outlier_peaks]#filtered_outliers

def Filter_Neighbors(outlier_list, changepoints, window_size = 100):
    if len(outlier_list) == 0:
        print(outlier_list)
        return outlier_list
    filtered_outlier_list = []
    i = 0
    cmp = outlier_list[0]
    outlier_set = set({cmp})
    for i in range(1, len(outlier_list)):
        curr = outlier_list[i]
        if (curr - cmp) <= window_size:
            if changepoints[curr] >= changepoints[cmp]:
                outlier_set.remove(cmp)
                outlier_set.update([curr])
                cmp = curr
            else:
                outlier_set.update([cmp])
        else:
            outlier_set.update([curr])
            cmp = curr
    outlier_set = list(outlier_set)
    outlier_set.sort()
    return outlier_set


# <h2> Rules for delinking contigs at change points. </h2> 
# 
# The rules for delinking the contigs are described below. The first row represent the orientation of the contig and the first colum represent the location of the contig relative to the contig at the change point. 
# 
# |#                   |  Forward            |  Reverse          |
# |:------------------:|---------------------|-------------------|
# |        Start       | Delink predecessor  | Delink successor  |
# |          End       | Delink successor    | Delink predecessor|
# 
# The input to the program is the set of outliers, the dictionary containing the contigs at each position, a vector of means for detecting change points and a *pos_cutoff* that specifies the number of basepairs from the contig ends to considered for delinking.

# In[7]:


def Get_Outlier_Contigs(outliers, positions, coordinates, graph, pos_cutoff = 100):
    g_ = deepcopy(graph)
    potential_contigs_removal = {}
    for o in outliers:
        contigs_intersecting = positions[o]
        closest_contig, closest_contig_val = '',np.inf
        forward, start = True, True
        for contig in contigs_intersecting:
            s,e = coordinates[contig]
            if s>e: f = False
            else: f = True
                
            if np.abs(s-o) < closest_contig_val:
                closest_contig_val= np.abs(s-o) 
                closest_contig = contig
                forward = f
            if np.abs(e-o) < closest_contig_val:
                closest_contig_val= np.abs(e-o) 
                closest_contig = contig
                forward = f
                start = False
        if closest_contig_val <= pos_cutoff:
            potential_contigs_removal[closest_contig] = (closest_contig_val, o, forward, start)
    
    for c in potential_contigs_removal:
        p = potential_contigs_removal[c]
        #print(c,p)
        fwd, start = p[2], p[3]       
        if(fwd == True) and (start == True): 
            print('Removing',c,'\'s Predecessors')
            predecessors = list(g_.predecessors(str(c)))
            if len(predecessors) > 0:
                edges = [(pred,c) for pred in predecessors]
                g_.remove_edges_from(edges)
        if(fwd == False) and (start == True): 
            print('Removing',c,'\'s Successors')
            successors = list(g_.successors(str(c)))
            if len(successors) > 0:
                edges = [(c,succ) for succ in successors]
                g_.remove_edges_from(edges)
        if(fwd == False) and (start == False): 
            print('Removing',c,'\'s Predecessors')
            predecessors = list(g_.predecessors(str(c)))
            if len(predecessors) > 0:
                edges = [(pred,c) for pred in predecessors]
                g_.remove_edges_from(edges)
        if(fwd == True) and (start == False): 
            print('Removing',c,'\'s Successors')
            successors = list(g_.successors(str(c)))
            if len(successors) > 0:
                edges = [(c,succ) for succ in successors]
                g_.remove_edges_from(edges)
    return g_
         


# <h2> Visualize the Results </h2>
# 
# In this code fragment we plot two figures, one describing 
# <ol>
#     <li>The Pileup of Contigs</li> 
#     <li>Coverages along the global Coordinate </li> 
#     <li>Change Point Metric </li> 
#     <li>The Coverages of the Various Contigs in the Scaffold</li> 
#     The program calls the methods above to perform the computations. 
# </ol>
# 
# The second figure plots the coverage of a contig with respect to its predecessors and successors. 

# In[31]:


def Return_Contig_Scaffold_Positions(coordinate_dictionary):
    pos_dict = {}
    for c in coordinate_dictionary:
        start, end = coordinate_dictionary[c]
        if start > end: a,b = end, start
        else: a,b = start,end
        for i in range(a, b+1):
            try: pos_dict[i].append(c)
            except KeyError: pos_dict[i] = [c]
    return pos_dict


def Write_Coverage_Outputs(graph,df_coverage):
    weakly_connected_components = list(nx.weakly_connected_component_subgraphs(graph))
    cc_ctr = 0
    list_coverages = []
    list_coords = []
    for i in range(len(weakly_connected_components)):
        test = weakly_connected_components[i]
        if nx.is_directed_acyclic_graph(test):
            print(i, len(test.nodes()))
            nodes = list(test.nodes())
            df_coverage_cc = df_coverage.loc[nodes]
            coverage, longest_coverage, df_depths = Compute_Coverage(test, df_coverage_cc)
            mean_ratios = Helper_Changepoints(deepcopy(coverage))
            outliers = ID_Peaks(mean_ratios)
            outliers = Filter_Neighbors(outliers, mean_ratios)
    
            top_sort = df_depths.index.tolist()
            coords = dict(zip(df_depths.index.tolist(), df_depths['Coords'].tolist()))
            Pos_Dict = Return_Contig_Scaffold_Positions(coords)
            df_depths = df_depths.loc[top_sort]

            g_removed = Get_Outlier_Contigs(outliers, Pos_Dict, coords, test, 100)
            delinked_conn_comps = list(nx.weakly_connected_component_subgraphs(g_removed))
            
            for cc in delinked_conn_comps:
                cc_ctr += 1
                coverage_cc, longest_coverage_cc, df_depths_cc = Compute_Coverage(cc, df_coverage_cc)
                coords_cc = dict(zip(df_depths_cc.index.tolist(), df_depths_cc['Coords'].tolist()))
                for i in range(len(coverage_cc)):
                    d = {'Connected_Component':cc_ctr, 'Coordinates':i,'Coverage':coverage_cc[i]}
                    list_coverages.append(d)
                for c in coords_cc:
                    d = {'Connected_Component':cc_ctr, 'Contig':c, 
                         'Start':coords_cc[c][0],'End':coords_cc[c][1]}
                    list_coords.append(d)
        print('\n')
    df_coords = pd.DataFrame(list_coords)
    df_coords = df_coords.sort_values(by = ['Connected_Component'])
    
    df_coverages_scaffolds = pd.DataFrame(list_coverages)
    df_coverages_scaffolds = df_coverages_scaffolds.sort_values(by = ['Connected_Component','Coordinates'])
    
    return df_coords, df_coverages_scaffolds

def Draw_Plots(connected_component, title, df_coverage):
    coverage, longest_coverage, df_depths = Compute_Coverage(connected_component, df_coverage)
    mean_ratios = Helper_Changepoints(deepcopy(coverage))
    outliers = ID_Peaks(mean_ratios)
    outliers = Filter_Neighbors(outliers, mean_ratios)
    
    top_sort = df_depths.index.tolist()
    coords = dict(zip(df_depths.index.tolist(), df_depths['Coords'].tolist()))
    Pos_Dict = Return_Contig_Scaffold_Positions(coords)
    df_depths = df_depths.loc[top_sort]
    
    contigs_set_outliers = set({})
    for o in outliers:
        contigs = Pos_Dict[o]
        contigs_set_outliers.update(contigs)
    contigs_set_outliers = list(contigs_set_outliers)
    
    fig_line_chart, ax_line_chart = plt.subplots(4, 1, figsize=(20, 40),
                    gridspec_kw={'height_ratios': [2.75 ,0.75, 0.75, 0.75]})
    ctr, ctr_inc  = 0.5, 0.5
    yticks, ytick_labels = [ctr], []
    for k in top_sort:
        start,end = coords[k]
        if start > end: color = 'red'
        else: color = 'blue'
        if k in contigs_set_outliers: color = 'green'
        ax_line_chart[0].plot([start, end],[ctr, ctr], linewidth = 3, color = color)
        ax_line_chart[0].axhline(ctr, color = 'black', linewidth=0.5)
        ctr += ctr_inc
        yticks.append(ctr)
        ytick_labels.append(k)
    ax_line_chart[0].set_yticks(yticks)
    colors = ['red', 'blue', 'green']
    lines = [Line2D([0], [0], color=c, linewidth=3) for c in colors]
    labels = ['Reverse', 'Forward','Outlier']
    ax_line_chart[0].legend(lines, labels)
    ax_line_chart[0].set_yticklabels(ytick_labels)
    #ax_line_chart[0].set_xlim([14000, 22000])
    
    ax_line_chart[1].plot(coverage, color = 'teal', label = 'Scaffold Coverage')
    ax_line_chart[1].plot(longest_coverage, color = 'orange', label = 'Longest Path Coverage')
    ax_line_chart[1].set_xlabel('Scaffold Coordinate')
    ax_line_chart[1].set_ylabel('Fold Coverage')
    ax_line_chart[1].legend()
    
    ax_line_chart[2].plot(mean_ratios, color = 'brown', linewidth = 2)
    ax_line_chart[2].set_ylabel('Change Points')
    for o in outliers:
        ax_line_chart[0].axvline(o, color = 'black', linewidth = 0.8, linestyle = '--')
        ax_line_chart[1].axvline(o, color = 'black', linewidth = 0.8, linestyle = '--')
        ax_line_chart[2].axvline(o, color = 'black', linewidth = 0.8, linestyle = '--')
    outliers.sort()
    temp = [0]+list(outliers)+[len(coverage)]
    for i in range(len(temp)-1): 
        if i%2 == 0: c = 'pink'
        else: c = 'black'
        ax_line_chart[0].axvspan(temp[i], temp[i+1], alpha=0.3, color=c)
        ax_line_chart[1].axvspan(temp[i], temp[i+1], alpha=0.3, color=c)
        ax_line_chart[2].axvspan(temp[i], temp[i+1], alpha=0.3, color=c)

    outlier_bars = []
    for i in range(0, len(top_sort)):
        if top_sort[i] in contigs_set_outliers:
            outlier_bars.append(i)
            
    bars = ax_line_chart[3].bar(df_depths.index, df_depths['Average_Depth'], color = 'green')
    ax_line_chart[3].set_xticklabels(df_depths.index, rotation=90)#, ha='right')
    for b in outlier_bars:
         bars[b].set_color('r')
    ax_line_chart[3].set_ylabel('Median Coverage')
    fig_line_chart.suptitle(title)
    fig_line_chart.tight_layout()
    fig_line_chart.subplots_adjust(top=0.97, left = 0.15)
    
    fig_depth_scatter, ax_depth_scatter = plt.subplots(2,1,figsize=(20,12), gridspec_kw={'height_ratios': [2,1]})
    pos, depth_ratio = 0, []
    for c in coords:
        depths_pred = list(df_depths.loc[list(connected_component.predecessors(c))]['Average_Depth'])
        depths_succ = list(df_depths.loc[list(connected_component.successors(c))]['Average_Depth'])
        x_pred = [pos for i in range(len(depths_pred))]
        x_succ = [pos for i in range(len(depths_succ))]
        ax_depth_scatter[0].scatter(x_pred, depths_pred, color = 'blue', marker = 'x')
        ax_depth_scatter[0].scatter(x_succ, depths_succ, color = 'olive', marker = '<')
        depth = df_depths.loc[c]['Average_Depth']
        depths = depths_pred+depths_succ
        depth_ratio.append(depth/max(depths))
        if max(depths) < depth and len(depths) > 2:
            ax_depth_scatter[0].scatter(pos, depth, color = 'black', s = 150, facecolors='none', linewidth = 3)
        ax_depth_scatter[0].scatter(pos, depth, color = 'red', s = 80, alpha = 0.6)
        ax_depth_scatter[0].axvline(pos, linewidth = 0.3, linestyle = '--', color = 'black')
        ax_depth_scatter[1].axvline(pos, linewidth = 0.3, linestyle = '--', color = 'black')
        pos += 1
    ax_depth_scatter[0].set_xticks([])
    
    ax_depth_scatter[1].plot(range(0, len(top_sort)), depth_ratio, marker = 'o', color = 'blue')
    ax_depth_scatter[1].set_xticks(range(0, len(top_sort)))
    ax_depth_scatter[1].set_xticklabels(top_sort, rotation = 90)
    ax_depth_scatter[1].axhline(1, color = 'black')
    
    fig_depth_scatter.suptitle(title)
    fig_depth_scatter.tight_layout()
    fig_depth_scatter.subplots_adjust(top=0.94, bottom = 0.25)
    
    g_removed = Get_Outlier_Contigs(outliers, Pos_Dict, coords, connected_component, 100)
    conns = list(nx.weakly_connected_component_subgraphs(g_removed))
    fig_removed, ax_removed = plt.subplots(len(conns),1, figsize = (20,40), 
                                          gridspec_kw={'height_ratios': [1.0/len(conns)]*len(conns)})

    for i in range(len(conns)):
        conn = conns[i]
        coverage, longest_coverage, df_depths = Compute_Coverage(conn, df_coverage)
        if len(conns) > 1:
            ax_removed[i].plot(coverage, color = 'red',label = 'Scaffold coverage')
            ax_removed[i].plot(longest_coverage, color = 'blue', label = 'Longest Path Coverage')
            ax_removed[i].legend()
            ax_removed[i].set_title('Number of Contigs:'+str(len(conn)))
        else:
            ax_removed.plot(coverage, color = 'red',label = 'Scaffold coverage')
            ax_removed.plot(longest_coverage, color = 'blue', label = 'Longest Path Coverage')
            ax_removed.legend()
            ax_removed.set_title('Number of Contigs:'+str(len(conn)))
    fig_removed.tight_layout()
    fig_removed.suptitle(title)
    fig_removed.subplots_adjust(top=0.94)
    
    return fig_line_chart, fig_depth_scatter, fig_removed


# In[32]:


samples = listdir(covpath)
samples.sort()
plt.rcParams.update(rcParams)
#samples = ['SRS019397']


# This part of the code just calls the "Draw_Plots" codes. This calls the various above methods. These are the various variables to look out for, 
# <ol>
#     <li>coverage = Scaffold level coverage</li>
#     <li>longest_coverage = Coverage along the longest path</li>
#     <li>mean_ratios = The ratios used to detect change points </li>
#     <li>outliers = Outliers detected on the change point signal </li>
#     <li>coords = The dictionary with the contigs as the keys and coordinates as the value</li>
#     <li>Pos_Dict = This dictionary gives the contig ids at every position along the scaffold</li>
# </ol>
# 
# Return whichever variable you want for your analysis.
# 
# If you want to run for a particular sample set that to the sample, you dont have to iterate through all the samples. 
# 
# Similarly just mention the connected component id you are interested in viewing the plot for the second loop. For an example the connected component we are looking at now is connected component id 16 in my machine. This could be different in various machines. We are looking at a scaffold with 73 contigs and the sample id is 'SRS012902.txt'. 
# 

# In[ ]:


for s in samples[16:]:
    sample = s.replace(".txt","")
    if sample[0] == '.':
        continue
    print(s)
    Graph_path= graphspath+sample+'/'+sample+'_scaffolds/oriented.gml'
    G = nx.read_gml(Graph_path)
    weakly_connected_components = list(nx.weakly_connected_component_subgraphs(G))
    df_coverage =Load_Read_Coverage(sample) 
    
    pdf_1 = PdfPages(op_path+sample+'.pdf')
    
    for i in range(len(weakly_connected_components)):
        test = weakly_connected_components[i]
        if len(test) >= 3 and nx.is_directed_acyclic_graph(test):
            nodes = list(test.nodes())
            coverage = df_coverage.loc[nodes]
            fig1, fig2, fig3 = Draw_Plots(test, 'Connected Component:'+str(i), coverage)
            pdf_1.savefig( fig1 )
            pdf_1.savefig( fig2 )
            pdf_1.savefig( fig3 )
            plt.close('all')
            print(i, len(weakly_connected_components))
            
    pdf_1.close()


# In[ ]:


'''for s in samples[4:]:
    sample = s.replace(".txt","")
    if sample[0] == '.':
        continue
    print(s)
    Graph_path= graphspath+sample+'/'+sample+'_scaffolds/oriented.gml'
    G = nx.read_gml(Graph_path)
    df_coverage =Load_Read_Coverage(sample) 
    df_coords, df_coverage_scaffolds = Write_Coverage_Outputs(G, df_coverage)
    df_coords.to_csv(op_path+sample+'_Gbl_Coords.csv')
    df_coverage_scaffolds.to_csv(op_path+sample+'_Coverages.csv')'''


# In[ ]:




