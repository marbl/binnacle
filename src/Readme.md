
<h1>Binnacle-Estimating Scaffold Coverages</h1>
 
 This repository is under active development. In this document we describe the various functions implemented in Compute_Scaffold_Coverages_Utility.py.
 
<h2>Computing Global Coordinates</h2>
Function: Compute_Global_Coordinates
<ol> 
    <li>This function computes the global coordinate system for the scaffold(connected component) based on the length of the contigs in the connected component and the linking information as returned by MetaCarvel. We order the contigs based on the topological sort of the subgraph and assign coordinates in the global coorinate system in breadth first manner.</li>   
    <li>We also consider the contig orientations and the edge orientations in accurate estimation of coordinates. </li>
    <li>If there multiple possible assignments we pick the one that has the largest coordinate value. We choose the greedy approaach because the number of solutions grows exponentially and the problem is NP-Hard! </li>
    <li>Finally, we normalize all other positions based on the least cooridnate value so that the coorinate system starts with 0.</li>
    <li>The input to the program is the connected component subgraph of the assembly graph ``oriented.gml" which is a result of running MetaCarvel.</li>
</ol>


<h2> Computing Coverages </h2>
Function: Load_Read_Coverage

To compute the read coverages we use the *genomecov* program, part of the *bedtools* suite. We run *genomecov* with *-d* optional enabled so that we get the per base depth. The output of the prgram is utilized to compute the depth along the scaffold and this function loads the out of genomecov as a dataframe.   

<h2>Computing the Depth of Coverage along the Scaffold's Global Coordinate System </h2>
Function: Compute_Coverage

To compute the depth we run the Compute_Global_Coordinates and Load_Read_Coverage functions and pass this to the Compute_Coverage described below. This program estimates the per cooridnate depth by using the outputs of the said two information. The program returns the coverage along the entire scaffold and the coverage along the longest path and a dataframe contiaining the contig, its average depth and its coordinates in the scaffold.

<h2> Change Point Detection on the Scaffold Coverage </h2>
Function: Helper_Change_Points

To compute changepoints along the scaffold we slide a window of default size 1500 along the scaffold and take the ratios of the maximum of means of the window preceeding any point and the window suceeding the point and the minimum of the said quantities. Any abnromal spike magnitude indicates a change point. If the size of the contig is less than the size of the window we adaptively pick a smaller window. 

<h2> Outlier Detection in Change Point Signals </h2>
Functions: ID_Outliers, ID_Peaks, Filter_Neigbors

<ol>
    <li>The following code segment is used to identify outliers in change points. To do so, we compute the peaks first.</li> 
    <li>We don't directly identify outliers on the signal because doing so would pickup points near the peaks, these are indicators of sliding windows and not really outliers.</li> 
    <li>To overcome the issue, we first identify peaks. A point is a peak if a point is larger than its predecessor and successor. This is performed by the function *ID_Peaks*.</li> 
    <li>The output of this passed to *ID_outliers* which picks all those points that is gretaer than the point which is the *thresh*'s percentile. The default value of *thresh* is 98.5. </li>
    <li>The filter outliers is still a work in progress. This is aimed at removing all the outliers points that is close to one another, but this is not super important. While the mthod described in the following block is data driven we are working on improving this method by examining the underlying graph structure.</li>
</ol>

<h2> Rules for delinking contigs at change points. </h2> 
Function: Get Outlier Contigs


The rules for delinking the contigs are described below. The first row represent the orientation of the contig and the first colum represent the location of the contig relative to the contig at the change point. 

|#                   |  Forward            |  Reverse          |
|:------------------:|---------------------|-------------------|
|        Start       | Delink predecessor  | Delink successor  |
|          End       | Delink successor    | Delink predecessor|

The input to the program is the set of outliers, the dictionary containing the contigs at each position, a vector of means for detecting change points and a *pos_cutoff* that specifies the number of basepairs from the contig ends to considered for delinking.

<h2> Visualize the Results </h2>

In this code fragment we plot two figures, one describing 
<ol>
    <li>The Pileup of Contigs</li> 
    <li>Coverages along the global Coordinate </li> 
    <li>Change Point Metric </li> 
    <li>The Coverages of the Various Contigs in the Scaffold</li> 
    The program calls the methods above to perform the computations. 
</ol>

The second figure plots the coverage of a contig with respect to its predecessors and successors.

For questions contact hsmurali@umd.edu
