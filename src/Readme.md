
<h1>Binnacle-Estimating Scaffold Coverages</h1>
 
 This repository is under active development. In this document we describe the various functions implemented in **Compute\_Scaffold\_Coverages\_Utility.py** and **Binnacle\_IO\_Utility.py**.
 
In the following section, we describe all the functions associated with **Compute\_Scaffold\_Coverages\_Utility.py**. 

<h2>Computing Global Coordinates</h2>

**Function**: Compute\_Global\_Coordinates <br/>
**Inputs**: subgraph of the weakly connected component, available vertex v<br/>
**Outputs**: dictionary describing the global coordinates for each vertex<br/>

<ol> 
    <li>This function computes the global coordinate system for the scaffold(connected component) based on the length of the contigs in the connected component and the linking information as returned by MetaCarvel. We assign (0, len(u)) where u is an "available vertex" of the subgraph. Available vertex is any vertex in the subgraph that has indegree 0. For all other vertices we assign coordinates in the global coorinate system in breadth first manner.</li>   
    <li>We also consider the contig orientations and the edge orientations in accurate estimation of coordinates. </li>
    <li>If there multiple possible assignments we pick the one that has the largest coordinate value. We choose the greedy approaach because the number of solutions grows exponentially and the problem can be shown to be NP-Complete! </li>
    <li>Finally, we normalize all other positions based on the least cooridnate value so that the coordinate system starts with 0.</li>
    <li>The input to the program is the connected component subgraph of the assembly graph "oriented.gml" which is a result of running MetaCarvel.</li>
</ol>



<h2>Computing the Depth of Coverage along the Scaffold's Global Coordinate System </h2>

**Function**: Compute_Coverage<br/>
**Input**: subgraph of the weakly connected component, available vertex u and the dataframe containing the coverage of each contig <br/>
**Output**: Coverage of the scaffold and per contigs depth<br/>

To compute the depth we run the Compute_Global_Coordinates and Load_Read_Coverage functions and pass this to the Compute_Coverage described below. This program estimates the per cooridnate depth by using the outputs of the said two information. The program returns the coverage along the entire scaffold and the coverage along the longest path and a dataframe contiaining the contig, its average depth and its coordinates in the scaffold.

<h2> Change Point Detection on the Scaffold Coverage </h2>

**Function**: Helper_Change_Points<br/>
**Input**: Coverage Vector, Window_Size(Default 1500bp)<br/>
**Output**: Change Point Statistic Vector<br/>

To compute changepoints along the scaffold we slide a window of default size 1500 along the scaffold and take the ratios of the maximum of means of the window preceeding any point and the window suceeding the point and the minimum of the said quantities. Any abnromal spike magnitude indicates a change point. If the size of the contig is less than the size of the window we adaptively pick a smaller window. 

<h2> Outlier Detection in Change Point Signals </h2>

**Function**: ID_Outliers<br/>
**Input**: Change Point Statistic Vector, threshold<br/>
**Output**: Outlier Indices<br/>

**Function**: ID_Peaks<br/>
**Input**: Change Point Statistic Vector, threshold <br/>
**Outputs**: Outlier Indices<br/>

**Function**: Filter_Neigbors<br/>
**Input**: Outlier List, Change Point Statistic Vector and Vicinity Window Size<br/>
**Output**: Filtered Outlier List<br/>

<ol>
    <li>The following code segment is used to identify outliers in change points. To do so, we compute the peaks first.</li> 
    <li>We don't directly identify outliers on the signal because doing so would pickup points near the peaks, these are indicators of sliding windows and not really outliers. 
    <li>To overcome the issue, we first identify peaks. A point is a peak if a point is larger than its predecessor and successor. This is performed by the function *ID_Peaks*.</li> 
    <li>The output of this passed to *ID_outliers* which picks all those points that is gretaer than the point which is the *thresh*'s percentile. The default value of *thresh* is 98.5. </li>
    <li>The filter outliers is still a work in progress. This is aimed at removing all the outliers points that is close to one another, but this is not super important. While the mthod described in the following block is data driven we are working on improving this method by examining the underlying graph structure.</li>
</ol>

<h2> Rules for delinking contigs at change points. </h2> 

**Function**: Get_Outlier_Contigs<br/>
**Input**: Outlier Positions, Subgraph of the Scaffold, Position cutoff to take action on the contig, and a dictionary describing the contig pileup at each position. <br/>
**Output**: subgraph after stratification. <br/>

The rules for delinking the contigs are described below. The first row represent the orientation of the contig and the first colum represent the location of the contig relative to the contig at the change point. 

|#                   |  Forward            |  Reverse          |
|:------------------:|---------------------|-------------------|
|        Start       | Delink predecessor  | Delink successor  |
|          End       | Delink successor    | Delink predecessor|

The input to the program is the set of outliers, the dictionary containing the contigs at each position, a vector of means for detecting change points and a *pos_cutoff* that specifies the number of basepairs from the contig ends to considered for delinking.

<h2> Contig Pileup Dictionary </h2>

**Function**: Return_Contig_Scaffold_Poistions<br/>
**Input**: Global Coordinate Dictionary<br/>
**Output**: Dictionary describing the pileup of contigs<br/>

We iterate through the global coordinate dictionary and return a dictionary that has the postion along the scaffold as the key and pile up of contig along each position as its values.

<h2> Computing an Available Vertex </h2>

**Function**: Return_Starting_Point<br/>
**Input**: Graph<br/>
**Output**: vertex with the lowest indegree and the indegree value <br/>

The above function returns the node with lowest indegree and the indegree value. If the minimum indegree is 0, then there is an available node and the program procceds by computing the global coordinates. If not we simplify the graph

<h2> Handling Non-DAG Elements </h2>

**Function**: Random_Simplify<br/>
**Input**: Graph<br/>
**Output**: Stratified Graph <br/>

Most of the scaffolds in the sample form a Directed Acyclic Graph. Some of the induce a cycle. Even if the scaffold has a cycle we look for an available vertex if not, the above code fragment converts the graph into a DAG by randonly snipping an edge in the simple cycle(s) in the scaffold. 

 
 
In the following section we describe the code fragments available in **Binnacle_IO_Utility.py** which deals with loading and processing the function described in this section.

<hr />

<h2> Loading Assembly Graph and Coverages </h2>

**Function**: Load_Read_Coverage_and_Assembly_Graph<br/>
**Input**: Path to the Assembly Graph, Path to the output og running genomecov<br/>
**Output**: Return a networkx graph object and dataframe containing the output of genomecov <br/>

To compute the read coverages we use the *genomecov* program, part of the *bedtools* suite. We run *genomecov* with *-d* optional enabled so that we get the per base depth. The output of the prgram is utilized to compute the depth along the scaffold and this function loads the out of genomecov as a dataframe.   


<h2> Computing Coverages and Writing the Results </h2> 

**Function**: Write_Coverage_Outputs<br/>
**Input** : networkx object of assembly graph, dataframe of coverage, output directory to dump outputs <br\>
**Output**: A tab seperated text file for each of the following, 
<ol>
    <li>Global Coordinates for the Scaffold Before Delinking with the Headers: 
        <ol><li>Scaffold id</li> 
            <li>Contig id</li> 
            <li>Start position in the global coordinate system</li> 
            <li>End position of the contig in the global coordinate system </li></ol>
    </li>
    <li>Global Coordinates for the Scaffold After Delinking with the Headers: 
        <ol><li>Scaffold id after delinking</li> 
            <li>Scaffold id before delinking</li> 
            <li>Contig id</li> 
            <li>Start position in the global coordinate system</li> 
            <li>End position of the contig in the global coordinate system</li></ol>
    </li>
    <li>Scaffold Coverages Before Delinking with the Headers: 
        <ol><li>Scaffold id before delinking</li> 
            <li>Location in the scaffold</li> 
            <li>coverage Value</li></ol>
    </li>
    <li>Scaffold Coverages After Delinking with the Headers: 
        <ol><li>Scaffold id after delinking</li> 
            <li>Scaffold id before delinking</li> 
            <li>Location in the scaffold</li> 
            <li>coverage Value</li></ol>
    </li>
    <li>Scaffold Level cover Summary Before Delinking with the Headers: 
        <ol><li>Scaffold id before delinking</li> 
            <li>Mean Coverage</li> 
            <li>Deviation in Coverage</li></ol>
    </li>
    <li>Scaffold Level cover Summary After Delinking with the Headers: 
        <ol><li>Scaffold id After delinking</li> 
            <li>Mean Coverage</li> 
            <li>Deviation in Coverage</li></ol>
    </li>
    
</ol>

The function imports the **Compute\_Scaffold\_Coverages\_Utility.py** and uses the function create the above output files. 
    

<h2> Creating Scaffold.fasta Files </h2>

**Function**: Load_FASTA_File<br/>
**Input**: Path to the fasta file containing the contig sequences. <br/>
**Output**: A dictionary containing the contig id as the keys and the sequences as values <br/>

**Function**: Get_Contigs_in_Scaffolds<br/>
**Input**: Path to the tsv file containing global coordinate information after delinking which is an output of running the Write_Coverage_Outputs<br/>
**Output**: A dictionary containing the scaffold id as the keys and the contigs in the scaffold as values <br/>

**Function**:  Write_Scaffolds<br/>
**Input**: Path to the fasta file of the contigs, Path to the tsv file containing global coordinate information after delinking which is an output of running the Write_Coverage_Outputs, output directory to write the scaffolds <br/>
**Output**: Writes the scaffolds to "Scaffolds.fasta" in the specified directory <br/>

The function  Write_Scaffolds calls the above tw methods to create the scaffold.fasta file. 
 * * *
 
<h2> Visualize the Results </h2>
This is under development....

For questions contact hsmurali@umd.edu


```python

```
