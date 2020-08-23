from Clustering_Utility import *

align_flag = "true"
binning = "binnacle"
out_dir = '/Users/harihara/Research-Activities/Data/Binnacle-Op/all_vs_all_alignments/SRS062427/'
Coords_Path = out_dir+'Coords_After_Delinking.txt'
summary_path = out_dir+'Summary_After_Delinking.txt'
align_dir = '/Users/harihara/Research-Activities/Data/Binnacle-Op/all_vs_all_alignments/genomecov_d/SRS062427/'

Format_Outputs(binning, Coords_Path, summary_path, out_dir, align_flag, align_dir)