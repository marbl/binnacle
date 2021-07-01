import argparse as ap
import sys
from Binnacle_IO_Utility import *
from Clustering_Utility import *
from packaging import version

def cmd_exists(cmd):
    '''
    Function to check if a linux cmd exists
    Input:
        cmd: Linux Command
    Output:
        True if command exists, false otherwise. 
    '''
    return subprocess.call("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

def Check_Dependencies():
    '''
    Checks if all the dependencies are met! 
    '''
    pd_version = pd.__version__
    np_version = np.__version__
    nx_version = nx.__version__
    bio_version = bio.__version__

    if version.parse(pd_version) < version.parse("0.23.4"):
        print("Incompatible Pandas version. Install Pandas (>=0.23.4)")
        return False
    if version.parse(np_version) < version.parse("1.15.4"):
        print("Incompatible Numpy version. Install Numpy (>=1.15.4)")
        return False
    if version.parse(nx_version) < version.parse("2.2"):
        print("Incompatible Networkx version. Install Networkx (>=2.2)")
        return False
    if version.parse(bio_version) < version.parse("1.73"):
        print("Incompatible Biopython version. Install Biopython (>=1.73)")
        return False
    if not cmd_exists('samtools'):
        print('Samtools does not exist in PATH. Terminating....\n')#, file=sys.stderr)
        return False
    if not cmd_exists('genomeCoverageBed'):
        print('Bedtools does not exist in PATH. Terminating....\n')#, file=sys.stderr)
        return False
    return True

def Get_Coverage_Wrapper(bedpath, bampath, coveragepath, contigspath, coordspath, op_dir, prefix):
    '''
    Function to compute contig coverages given alignment information. 
    Input:
        bedpath: A bedfile describing coverages
        bampath: A bamfile describing coverages
        coveragepath: The output of running genomeCoverageBed with -bga -split flags
        contigspath: A fasta file describing contigs
        coordspath: The coordinate system as computed by binnacle
        op_dir: Location to write outputs to
        prefix: Prefix to attach to output files
    Output:
        coverage_path: The output of running genomeCoverageBed with -bga -split flags
    '''
    if len(bedpath) > 0:
        if(not isfile(bedpath)):
            print(bed_path +" file not found")
            sys.exit(1)
        if len(contigspath) > 0:
            if(not isfile(contigspath)):
                print(contigspath +" file not found")
                sys.exit(1)    
            faidx_cmd = 'samtools faidx '+contigspath 
            result = subprocess.getoutput(faidx_cmd)
            cut_cmd = 'cut -f 1,2 '+contigspath+'.fai > '+op_dir+prefix+'.length.txt'
            result = subprocess.getoutput(cut_cmd)
            length_path = op_dir+prefix+'.length.txt'
        elif len(coordspath) > 0:
            if(not isfile(coordspath)):
                print(coordspath + " file not found")
                sys.exit(1)
            df_coords = pd.read_csv(coordspath, names = ['cc_aft_dlink', 'cc_bef_dlink', 'Contig', 'Start', 'End', 'Ingraph', 'Length'], 
                                sep = '\t')
            df_coords['Contig'] = df_coords['Contig'].astype(str)
            df_length = df_coords[['Contig','Length']]
            df_length = df_length.set_index('Contig')
            df_length.to_csv(op_dir+prefix+'.length.txt', sep = "\t", header = False)
            length_path = op_dir+prefix+'.length.txt'
        else:
            print("Please provide either coords file or the contigs file. \
                   Failed to compute genome lengths for estimating perbase coverages from bedfile...\n")
            sys.exit(1)
        coverage_path = Get_Coverage_Bed(bedpath, length_path, op_dir, prefix)
    elif len(bampath) > 0:
        if(not isfile(bampath)):
            print(bampath +" file not found")
            sys.exit(1)
        coverage_path = Get_Coverage_BAM(bampath, op_dir, prefix)
    elif len(coveragepath) > 0:
        if(not isfile(coveragepath)):
            print(coveragepath +" file not found")
            sys.exit(1)
        coverage_path = coveragepath
    else:
        print("Please provide coverages as a BAM/BED/TXT file describing the perbase coverages...\n")
        sys.exit(1)
    return coverage_path

def Get_Coverage_BAM(bampath, opdir, prefix=""):
    '''
    Function to run genomeCoverageBed on bam file with -bga -split on ba file describing the alignments. 
    Inputs:
        bampath: A bamfile describing coverages.
        opdir: Location to write outputs to
        prefix: Prefix to attach to output files
    '''
    cov_cmd = 'samtools sort -@ 4 -T '+opdir+' '+bampath+' | genomeCoverageBed -ibam stdin -bga -split > ' + opdir+prefix+'.read.coverage.txt'
    result = subprocess.getoutput(cov_cmd)
    sort_cmd = 'LC_ALL=C sort -k1,1 '+opdir+prefix+'.read.coverage.txt > ' + opdir+prefix+'.read.coverage.sorted.txt'
    result = subprocess.getoutput(sort_cmd)
    remove(opdir+prefix+'.read.coverage.txt')
    return opdir+prefix+'.read.coverage.sorted.txt'

def  Get_Coverage_Bed(bedpath, lengthpath, opdir, prefix=""):
    '''
    Function to run genomeCoverageBed on bedfile with -bga -split on ba file describing the alignments. 
    Inputs:
        bedpath: A bedfile describing coverages.
        lengthpath: A genome file describing the lenghs of the contigs
        opdir: Location to write outputs to
        prefix: Prefix to attach to output files
    '''
    cov_cmd = 'sort -k 1,1 '+bedpath+'| genomeCoverageBed -i stdin -g '+lengthpath+' -bga -split > '+opdir+prefix+'.read.coverage.txt'
    result = subprocess.getoutput(cov_cmd)
    sort_cmd = 'LC_ALL=C sort -k1,1 '+opdir+prefix+'.read.coverage.txt > ' + opdir+prefix+'.read.coverage.sorted.txt'
    result = subprocess.getoutput(sort_cmd)
    remove(opdir+prefix+'.read.coverage.txt')
    return opdir+prefix+'.read.coverage.sorted.txt'