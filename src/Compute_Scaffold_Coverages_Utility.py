#!/usr/bin/env python
# coding: utf-8

'''
Program developed at Pop lab at the CBCB, University of Maryland by 
Harihara Subrahmaniam Muralidharan, Nidhi Shah, Jacquelyn S Meisel. 
'''

import numpy as np
import pandas as pd
import networkx as nx
from copy import deepcopy

def Mean(group):
    '''
    Function to estimate the mean coverage of a pandas group.
    Input:
        group: Pandas group grouped by scaffold id or contig
    Output:
        returns the mean of the group
    '''
    freq = group['Length'].tolist()
    val = group['Coverage'].tolist()
    vec = np.repeat(val, freq)
    return round(np.mean(vec),1)

def Deviation(group):
    '''
    Function to estimate the standard deviation of coverage of a pandas group.
    Input:
        group: Pandas group grouped by scaffold id or contig
    Output:
        returns the deviation of the group
    '''
    freq = group['Length'].tolist()
    val = group['Coverage'].tolist()
    vec = np.repeat(val, freq)
    return round(np.std(vec),1)

def Percentile(group, p):
    '''
    Function to estimate the pth percentile of coverage of a pandas group.
    Input:
        group: Pandas group grouped by scaffold id or contig
    Output:
        returns the pth percentile of the group
    '''
    freq = group['Length'].tolist()
    val = group['Coverage'].tolsit()
    vec = np.repeat(val,freq)
    return round(np.percentile(vec, p))
    
def Summarize_Coverages(df):
    '''
    Function to estimate the summary of coverage values of a pandas group.
    Input:
        group: Pandas dtafarame containing coverages
    Output:
        returns the summary of coverage of each contig in the dataframe.
    '''
    df['Length'] = df['End']-df['Start']
    df_length = df[['ContigID', 'End']].groupby('ContigID').max().rename(columns = {'End':'Length'})
    df_op = pd.DataFrame()
    df_op['Mean'] = df.groupby('ContigID').apply(Mean)
    df_op['Std'] = df.groupby('ContigID').apply(Deviation)
    df_op = df_op.join(df_length)
    df_op = df_op[['Length','Mean','Std']]
    return df_op

def Return_Starting_Point(subgraph):
    '''
    Function to return the starting point to assign coordinates on the global frame of reference. 
    The starting node has an indegree equal to 0
    Input:
        subgraph: The graph pertaining to the DAG.
    Output:
        min_node: Node with the lowest indegree
        min_indegree: min_node's indegree
    '''
    nodes = list(subgraph.nodes())
    min_node, min_indegree = '', np.inf
    for n in nodes:
        if subgraph.in_degree(n) < min_indegree and subgraph.out_degree(n) > 0:
            min_node, min_indegree = n, subgraph.in_degree(n)
    return min_node, min_indegree

def Random_Simplify(subgraph_, min_indegree_node):
    '''
    Function to randomly simplify a simple cycle. This routine is only called if there are no available 
    nodes in the graph to assign coordinates. 
    Input: 
        subgraph_: A graph to simplify and usually contains one or more simple cycles.
        min_indegree_node: Returns the node to simplify. We choose to simplify the node with the minimum 
        indegree to simplify.
    Output:
        subgraph: Returns the graph after simplification.
    '''
    subgraph = deepcopy(subgraph_)
    cycles = nx.simple_cycles(subgraph)
    edges_removal_set = []
    for c in cycles:
        edge = (c[-1], c[0])
        edges_removal_set += [edge]
    edges_removal_set = set(edges_removal_set)
    subgraph.remove_edges_from(list(edges_removal_set))
    return subgraph

def Compute_Global_Coordinates(subgraph, v_0):
    '''
    Function to assign coordinates to contigs in the Directed Acyclic Graph(DAG). 
    Input: 
        Subgraph: A graph pertaining to the DAG.
        v_0: Starting node to start assigning coordinate values. 
    Output:
        Returns a dictionary with contig ids as keys and coordinate values as values. 
    '''
    v_0_orientation = subgraph.nodes[v_0]['orientation']
    v_0_length = int(subgraph.nodes[v_0]['length'])
    if v_0_orientation == 'FOW':
        start = 0
        end = start + v_0_length
    elif v_0_orientation == 'REV':
        start = 0
        end = start-v_0_length
    global_coords = {v_0 : (start,end)}
    nodes = list(subgraph.nodes())
    visited = dict(zip(nodes, [False]*len(nodes)))
    Q = [v_0]
    g_undir = subgraph.to_undirected()
    
    while len(Q) > 0:
        src = Q.pop()
        visited[src] = True
        neighbors = list(g_undir.neighbors(src))
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

def Compute_Coverage(df_coverage, coords):
    '''
    Function to compute the coverages of the DAG based on the coordinates assigned by the previous functions. 
    Input:
        df_coverage: A pandas dataframe containing the perbase of coverage of the contigs in the DAG.
        coords: A dictionary containing the start and end points of the contig in the global frame of reference. 
    Output:
        A coverage vector for the entire DAG
    '''
    top_sort = list(coords.keys())
    max_coord = -np.inf
    for c in top_sort:
        max_v = max(coords[c])
        if max_v >= max_coord: max_coord = max_v
    coverage = np.zeros(max_coord+1)
    cov_ctr = df_coverage.reset_index().groupby(['Contig']).count()
    cov_ctr = dict(zip(cov_ctr.index.tolist(), cov_ctr['Start'].tolist()))
    try:
        for c in top_sort:
            s,e = coords[c]
            cov_coords = df_coverage.loc[c]
            if cov_ctr[c] > 1: loc, cov_contig = list(zip(cov_coords['Start'].tolist(), cov_coords['End'].tolist())),cov_coords['coverage'].tolist()
            else:  loc, cov_contig = list(zip([cov_coords['Start']], [cov_coords['End']])), [cov_coords['coverage']]
            contig_depth = np.zeros(np.abs(s-e))
            for i in range(len(loc)):
                l = loc[i]
                contig_depth[l[0]:l[1]] = cov_contig[i]
            if (s > e): coverage[s:e:-1] += contig_depth
            else: coverage[s:e] += contig_depth
            assert np.abs(s-e) == len(contig_depth), top_sort
    except Exception as e:
        print(e)
        print('Here')
        print('Did you sort the coverage files by contig ids before running binnacle???? Exiting with error...\n')
    return coverage

def Helper_Changepoints_Z_Stat(cov_vec, window_size):
    
    '''Function to compute outliers in coverage signals. We use the two sample two sample Z-tests. 
    The z-statistic is calculated by z = \frac{\mu_1-\mu_2}{\sqrt{\sigma_1^2+\sigma_2^2}}
    Input:
        cov_vec: The Coverage Vector estimated by the previous function. 
        window_size: Defaulting to 1500, this is the window size to identify change points. 
    Output:
        cpts: A vector of the change point statistic'''
    
    indices_non_zero = np.where(cov_vec > 0)
    sliced = cov_vec[indices_non_zero]
    while window_size >= len(sliced)/2:
        window_size = int(window_size/5)
    mean = np.array(pd.Series(sliced).rolling(window_size).mean())[window_size-1:]
    sd = np.array(pd.Series(sliced).rolling(window_size).std())[window_size-1:]
    mu_1, mu_2 = mean[0:len(mean)-window_size-1], mean[window_size+1:]
    sd_1, sd_2 = sd[0:len(sd)-window_size-1],sd[window_size+1:]
    test_stat = [0]*window_size + list((mu_1 - mu_2)/(np.sqrt(sd_1**2+sd_2**2))) + [0]*window_size
    cpts = np.zeros(len(cov_vec))
    cpts[indices_non_zero] = test_stat
    return cpts

def ID_outliers(change_point_vec, thresh):
    '''
    Function to eastimate outliers in the change point statistic vector. 
    Input:
        change_point_vec: A vector of changepoints. 
        thresh: Thresh to identify outliers. 
    Outpu:
        Indices: The indices on the change point statistic vector that are outliers. 
    '''
    cutoff_upper = np.percentile(change_point_vec, thresh)
    cutoff_lower = np.percentile(change_point_vec, 100-thresh)
    indices = np.where(((change_point_vec >= cutoff_upper) | (change_point_vec <= cutoff_lower))&(change_point_vec != 0))
    return indices[0]

def Filter_Neighbors(outlier_list, changepoints, window_size):
    '''
    Function to remove multiple outliers around a peak. 
    Input:
        outlier_list: The list of outlier indices as identified by the previous routines
        changepoints:  a vector of the change point statistic
        window_size: Defaulting to 100(typically the read length)
    Output:
        outlier_set: Filteres set of outliers after removing the neughbouring outliers.  
    '''
    if len(outlier_list) == 0:
        return outlier_list
    filtered_outlier_list = []
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

def Get_Outlier_Contigs(outliers, positions, coordinates, graph, pos_cutoff):
    '''
    Function to identify outlier contigs links based on change point using the following rules
    |#                   |  Forward            |  Reverse          |
    |:------------------:|---------------------|-------------------|
    |        Start       | Delink predecessor  | Delink successor  |
    |          End       | Delink successor    | Delink predecessor|
    Input:
        outliers: Set of outliers identified by running the previous routines. 
        positions: A dictionary that has keys equal to the length of the span of the scaffold 
                   and for each position we record the contigs intersecting at that coordinate
        coordinates: A dictionary contianing the startting and ending positions along a 
                     global frame reference for a scaffold for each contig. 
        graph: A graph of the scaffold.
        pos_cutoff: Outlier position on the scaffold to consider for delinking
    Output:
        g_ : delinked graph
    '''
    g_ = deepcopy(graph)
    potential_contigs_removal = {}
    counter_end_points = 0
    for o in outliers:
        try:
            contigs_intersecting = positions[o]
        except KeyError:
            print('Here')
            print('KeyError', o, np.max(list(positions.keys())))
            continue
        closest_contig, closest_contig_val = '',np.inf
        forward, start = True, True
        o_flag = False
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
            if np.abs(s-o) <= 3*pos_cutoff or np.abs(e-o) <= 3*pos_cutoff:
                o_flag = True
        if closest_contig_val <= pos_cutoff:
            potential_contigs_removal[closest_contig] = (closest_contig_val, o, forward, start)
        if o_flag:#closest_contig_val <= 2*pos_cutoff:
            counter_end_points += 1
            
    for c in potential_contigs_removal:
        p = potential_contigs_removal[c]
        fwd, start = p[2], p[3]       
        if(fwd == True) and (start == True): 
            predecessors = list(g_.predecessors(str(c)))
            if len(predecessors) > 0:
                edges = [(pred,c) for pred in predecessors]
                g_.remove_edges_from(edges)
        if(fwd == False) and (start == True): 
            successors = list(g_.successors(str(c)))
            if len(successors) > 0:
                edges = [(c,succ) for succ in successors]
                g_.remove_edges_from(edges)
        if(fwd == False) and (start == False): 
            predecessors = list(g_.predecessors(str(c)))
            if len(predecessors) > 0:
                edges = [(pred,c) for pred in predecessors]
                g_.remove_edges_from(edges)
        if(fwd == True) and (start == False): 
            successors = list(g_.successors(str(c)))
            if len(successors) > 0:
                edges = [(c,succ) for succ in successors]
                g_.remove_edges_from(edges)
    return g_

def Return_Contig_Scaffold_Positions(coordinate_dictionary):
    '''
    Function to return a dictionary that has keys equal to the span of the scaffold and 
    for each position we record the contigs intersecting at that coordinate
    Input: 
        coordinate_dictionary: A dictionary contianing the startting and ending positions along a global 
                               frame reference for a scaffold for each contig.
    Output:
        pos_dict: The forementioned dictionary that has keys that are positions and the values with
                  the contigs at that position. 
    '''
    pos_dict = {}
    for c in coordinate_dictionary.keys():
        start, end = coordinate_dictionary[c]
        if start > end: a,b = end, start
        else: a,b = start,end
        for i in range(a, b+1):
            try: 
                pos_dict[i].append(c)
            except KeyError: 
                pos_dict[i] = [c]
    return pos_dict