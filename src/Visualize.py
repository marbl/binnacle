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