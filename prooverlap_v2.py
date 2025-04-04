#Import modules
from Bio import SeqIO
import pandas as pd
import statistics
from scipy import stats
import sys
import argparse
import os
import tempfile
import time
from pybedtools import BedTool
import pybedtools
import random
import warnings
from collections import Counter
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from multiprocessing import Pool


# Suppress all warnings
warnings.filterwarnings("ignore")
# Suppress specific warnings
warnings.filterwarnings("ignore", message="the index file is older than the FASTA file")

#Create parser
parser = argparse.ArgumentParser()
parser.add_argument("--mode", required = True, help = "Define mode: intersect or closest, intersect count the number of overlapping elements while closest test the distance. In closest if a feature overlap a target the distance is 0")
parser.add_argument("--input", required = True, help="Input bed file, accepted format 6 or 12 column BED files, name and score can be placeholder and are not used in standard analysis, score is used in Rank analysis, strand is used only if some strandess test are requested")
parser.add_argument("--targets", required = True, help="Target bed files (must contain 6-12 columns), to test enrichement against, if multiple files are supplied N independent test against each file are conducted, file names must be comma separated, the name of the file will be use as the name output")
parser.add_argument("--background", help="Background bed file (must contain 6-12 columns) to test enrichement aginst, should be a superset from wich input bed file is derived")
parser.add_argument("--randomization", help="Number of randomization, default 100", default = 100)
parser.add_argument("--genome", help="genome fasta file used to retrieve sequence from bed files, needed only for length or AT/GC content")
parser.add_argument("--tmp", default=".", help="Default is current working dir, Location of the directory to store temporary files, after running the sofware automatically clean up tmp files, if the software do not exit properly it may not clean up tmp file!")
parser.add_argument("--outfile", help="Full path to output file to store results, it will be created")
parser.add_argument("--orientation", required = True, help="Name of test/tests to be performed: concordant, discordant, strandless, or a combination of them, comma separated, no space allowed")
parser.add_argument("--ov_fraction", default="0.000000000000000000000000000000000000000001", help="Minimum overlap required as a fraction from input BED file to consider 2 features as overlapping. Default is 1E-9 (i.e. 1bp)")
parser.add_argument("--generate_bg", action = "store_true", help="Generates random bed intervals to test enrichment against, use this instead of background. Use only if background file cannot be used")
parser.add_argument("--exclude_intervals", default = None, help="Exclude those regions in both random background generation and feature testing")
parser.add_argument("--exclude_ov", action = "store_true", help="Do not count overlapping region in closest mode")
parser.add_argument("--exclude_upstream", action = "store_true", help="Do not count upstream region in closest mode, only for stranded files, not compatible with exclude_downstream")
parser.add_argument("--exclude_downstream", action = "store_true", help="Do not count downstream region in closest mode, only for stranded files, not compatible with exclude_upstream")
parser.add_argument("--test_AT_GC", action = "store_true", help="Test AT and GC content")
parser.add_argument("--test_length", action = "store_true", help="Test feature length")
parser.add_argument("--GenomicLocalization", action = "store_true", help= "Test the genomic localization and enrichment of founded overlaps, i.e TSS,Promoter,exons,introns,UTRs - Available only in intersect mode. Must provide a GTF file to extract genomic regions (--gtf), alternatively directly provide a bed file (--bed) with custom annotations")
parser.add_argument("--gtf", help="GTF file, only to test genomic localization of founded overlap, gtf file will be used to create genomic regions: promoter, tss, exons, intron, 3UTR and 5UTR")
parser.add_argument("--bed", help="BED file, only to test genomic localization of founded overlap, bed file will be used to test enrichment in different genomic regions, annotation must be stored as 4th column in bed file, i.e name field")
parser.add_argument("--RankTest", action = "store_true", help="Running Rank analyis")
parser.add_argument("--Ascending_RankOrder", action = "store_true", help="Sort Ascending")
parser.add_argument("--WeightRanking", action = "store_true", help="Weight the ranking test, this is done by multiply the score value in the BED file by the distance or number of overlapped nucleotides")
parser.add_argument("--alpha", default = "0.1", help="Strength of the Weight for the ranking test, only if --WeightRanking is active")
parser.add_argument("--thread", default = "1", help="Number of Threads")

#Parse argument
args = parser.parse_args()

################################## Functions ############################################

#Create directory if not exist, if exist print it and exit
def create_directory(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Directory created: {folder_path}")
    else:
        print("tmp dir exists, use:" + folder_path)

#Check input parameters
def check_input_parameters(mode):
    if mode == "intersect" or mode == "closest":
        print("Mode set to " + mode)
    else:
        print("Mode is not set, choose one of intersect or closest")
        print("Exit")
        sys.exit()

#Read bed file and if have more than 6 columns return only the first 6 columns
def get_bed_file(input_path):
    bed_df = pd.read_csv(input_path, sep="\t", header=None, comment='#')
    # Keep only the first 6 columns
    bed_df = bed_df.iloc[:, :6]
    return bed_df

#Get strandness parameters
def get_params_strand(strandness):
    if strandness == "concordant":
        return({"s": True, "S": False})
    if strandness == "discordant":
        return({"s": False, "S":True})
    if strandness == "strandless":
        return({"s": False, "S":False})
    else:
        print("Check strandness parameters: choose one or more from concordant, discordant or strandless")
        print("Exit")
        sys.exit()

#Test closeness
def test_closeness(bed,target,background,ov_fraction, randomization, frame, table,strandness="strandless", exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False):
    #get_strandness_parameters
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    if exclude_intervals is not None:
        print("Running closest in " + strandness + "mode and exclude intervals from ...." + exclude_intervals)
        real_count = statistics.mean([abs(x) for x in bed.sort().intersect(b = pybedtools.BedTool(exclude_intervals),s = s, S = S, v = True).closest(b = target.sort().intersect(b = pybedtools.BedTool(exclude_intervals),s = s, S = S, v = True), d=True,s=s,S=S,f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe().iloc[:,12].to_list()])
        print("Mean real distance is {0} using {1} elements".format(real_count, strandness))
        random_count = []
        background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals), s = s, S = S, v = True).to_dataframe() 
    else:    
        real_count = statistics.mean([abs(x) for x in bed.sort().closest(b = target.sort(), d=True,s=s,S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe().iloc[:,12].to_list()])
        print("Mean real distance is {0} for {1} elements".format(real_count, strandness))
        random_count = []
        background_df = background.to_dataframe()       
    #Begin randomizations
    for i in range(int(randomization)):
        random_features_index = random.sample(range(background_df.shape[0]), len(bed))
        background_df_random = background_df.iloc[random_features_index,:]
        random_count.append(statistics.mean([abs(x) for x in pybedtools.BedTool.from_dataframe(background_df_random).sort().closest(b = target.sort(), d=True, s=s, S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe().iloc[:,12].to_list()]))
    zscore = compute_z_score(real_count,random_count)
    pvalue = compute_pvalue(zscore)
    table = save_tables(real_count, random_count, strandness, os.path.basename(target.fn), table)
    frame = save_results(zscore,strandness,pvalue, os.path.basename(target.fn), real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table

#Test closeness against a random background  
def test_closeness_random_bg(bed, target, ov_fraction, randomization, frame, table,strandness="strandless", exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False):
    #get_strandness_parameters
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]  
    if exclude_intervals is not None:
        print("Excluding intervals")
        real_count = statistics.mean([abs(x) for x in bed.intersect(b = pybedtools.BedTool(exclude_intervals), v = True, s=s, S=S).sort().closest(b = target.sort().intersect(b = pybedtools.BedTool(exclude_intervals),s = s, S = S, v = True), d=True,s=s, S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe().iloc[:,12].to_list()])
        print("Mean real distance is: {0}, for {1} ".format(real_count, strandness))
        random_count = []
    else:
        print("Not excluding intervals")
        real_count = statistics.mean([abs(x) for x in bed.sort().closest(b = target.sort(), d=True,s=s,S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe().iloc[:,12].to_list()])
        print("Mean real distance is: {0}, for {1} ".format(real_count, strandness))
        random_count = []
    for i in range(int(randomization)):
        input_bed_file = bed.fn
        #Read the input BED file and calculate the mean length and chromosomal frequency to produce a correct random bg
        existing_intervals, mean_length, chrom_count = read_bed_file(input_bed_file)
        # Read the exclude regions
        if exclude_intervals is not None:
            print("Exclude intervals")
            exclude_intervals_bed = read_exclude_bed_file(exclude_intervals)
        else:
            print("Not excluding")
            exclude_intervals_bed = None

        # Generate the same number of random intervals as the input file, using the chromosomal frequencies
        num_random_intervals = len(existing_intervals)
        random_intervals = generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed)

        # Write the random intervals to a new BED file
        write_bed_file("Random_0123456789.bed", random_intervals)

        background_df_random = "Random_0123456789.bed"
        random_count.append(statistics.mean([abs(x) for x in pybedtools.BedTool(background_df_random).sort().closest(b = target.sort(), d=True, s= s, S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe().iloc[:,12].to_list()]))


    zscore = compute_z_score(real_count,random_count)
    pvalue = compute_pvalue(zscore)
    table = save_tables(real_count, random_count, strandness, os.path.basename(target.fn), table)
    frame = save_results(zscore,strandness, pvalue, os.path.basename(target.fn), real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table  

#Test length against a random background
def test_length_random_bg(bed,randomization, target, frame, table, genome, exclude_intervals = None):
    nucleotide_content = bed.nucleotide_content(fi = genome).to_dataframe()
    real_length = statistics.mean(nucleotide_content["15_seq_len"].to_list())
    random_count_length = []
    for i in range(int(randomization)):
        # Example usage:
        input_bed_file = args.input
        # Read the input BED file and calculate the mean length and chromosomal frequency
        existing_intervals, mean_length, chrom_count = read_bed_file(input_bed_file)
        # Read the exclude regions from the third BED file
        if exclude_intervals is not None:
            exclude_intervals_bed = read_exclude_bed_file(exclude_intervals)
        else:
            exclude_intervals_bed = None
        # Generate the same number of random intervals as the input file, using the chromosomal frequencies
        num_random_intervals = len(existing_intervals)
        random_intervals = generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed)
        # Write the random intervals to a new BED file
        write_bed_file(args.tmp + "/" + "Random_0123456789.bed", random_intervals)
        background_df_random = args.tmp + "/" + "Random_0123456789.bed"    
        random_count_length.append(statistics.mean(pybedtools.BedTool(background_df_random).nucleotide_content(fi = genome, s = True).to_dataframe()["15_seq_len"]))

    zscore = compute_z_score(real_length,random_count_length)
    pvalue = compute_pvalue(zscore)

    frame = save_results(zscore, "Length", pvalue, os.path.basename(target), real_length, statistics.mean(random_count_length), statistics.stdev(random_count_length), frame)
    table = save_tables(real_length, random_count_length, "Length", os.path.basename(target), table)
    return frame, table

#Test AT and GC content against a random background
def test_GC_AT_random_bg(bed,randomization, target, frame, table, genome, exclude_intervals = None):
    nucleotide_content = bed.nucleotide_content(fi = genome, s = True).to_dataframe()
    real_GC = statistics.mean(nucleotide_content["8_pct_gc"].to_list())
    real_AT = statistics.mean(nucleotide_content["7_pct_at"].to_list())
    random_count_GC = []
    random_count_AT = []
    for i in range(int(randomization)):
        # Example usage:
        input_bed_file = args.input
        # Read the input BED file and calculate the mean length and chromosomal frequency
        existing_intervals, mean_length, chrom_count = read_bed_file(input_bed_file)

        # Read the exclude regions from the third BED file
        if exclude_intervals is not None:
            exclude_intervals_bed = read_exclude_bed_file(exclude_intervals)
        else:
            exclude_intervals_bed = None

        # Generate the same number of random intervals as the input file, using the chromosomal frequencies
        num_random_intervals = len(existing_intervals)
        random_intervals = generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed)

        # Write the random intervals to a new BED file
        write_bed_file(args.tmp + "/" + "Random_0123456789.bed", random_intervals)

        background_df_random = args.tmp + "/" + "Random_0123456789.bed"

        random_count_GC.append(statistics.mean(pybedtools.BedTool(background_df_random).nucleotide_content(fi = genome, s = True).to_dataframe()["8_pct_gc"]))
        random_count_AT.append(statistics.mean(pybedtools.BedTool(background_df_random).nucleotide_content(fi = genome, s = True).to_dataframe()["7_pct_at"]))

    zscore_GC = compute_z_score(real_GC,random_count_GC)
    pvalue_GC = compute_pvalue(zscore_GC)
    zscore_AT = compute_z_score(real_AT,random_count_AT)
    pvalue_AT = compute_pvalue(zscore_AT)

    frame = save_results(zscore_GC, "GC", pvalue_GC, os.path.basename(target), real_GC, statistics.mean(random_count_GC), statistics.stdev(random_count_GC), frame)
    frame = save_results(zscore_AT, "AT", pvalue_AT, os.path.basename(target), real_AT, statistics.mean(random_count_AT), statistics.stdev(random_count_AT), frame)
    table = save_tables(real_GC, random_count_GC, "GC", os.path.basename(target), table)
    table = save_tables(real_AT, random_count_AT, "AT", os.path.basename(target), table)
    return frame, table

#Test overlap against a random background
def test_enrichement_random_bg(bed, target, ov_fraction, randomization, frame, table, strandness = "strandless", exclude_intervals = None):
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    real = bed.intersect(b = target, u = True, s=s, S=S, f = ov_fraction)
    real_count = len(real)
    print("There are {0} {1} overlapping elements".format(real_count, strandness))
    random_count = []
    
    for i in range(int(randomization)):
        input_bed_file = args.input
        # Read the input BED file and calculate the mean length and chromosomal frequency
        existing_intervals, mean_length, chrom_count = read_bed_file(input_bed_file)

        # Read the exclude regions from the third BED file
        if exclude_intervals is not None:
            exclude_intervals_bed = read_exclude_bed_file(exclude_intervals)
        else:
            exclude_intervals_bed = None


        # Generate the same number of random intervals as the input file, using the chromosomal frequencies
        num_random_intervals = len(existing_intervals)
        random_intervals = generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed)

        # Write the random intervals to a new BED file
        write_bed_file(args.tmp + "/" + "Random_0123456789.bed", random_intervals)

        background_df_random = args.tmp + "/" + "Random_0123456789.bed"
        random_count.append(len(pybedtools.BedTool(background_df_random).intersect(b = target, u = True, s=s, S=S, f = ov_fraction)))


    zscore = compute_z_score(real_count,random_count)
    pvalue = compute_pvalue(zscore)
    table = save_tables(real_count, random_count, strandness, os.path.basename(target.fn), table)
    frame = save_results(zscore,strandness, pvalue, os.path.basename(target.fn), real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table

def test_enrichement_regions_random_bg(bed, target, ov_fraction, randomization, frame, table, strandness = "strandless", exclude_intervals = None, regions_bed = "NA"):

    #read the bed file with annotations
    regions_bed = pd.read_table(args.tmp + "/" + "Reference_regions.bed")
    #Set negative regions to 0, this is bevcause some promoters may be calculated and negative start if the chromosome is too short
    regions_bed.iloc[:, 1] = regions_bed.iloc[:, 1].clip(lower=0)
    # Split the dataframe based on the values in column number 4, i.e name
    regions_dict = {category: group for category, group in regions_bed.groupby(regions_bed.columns[3])}

    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    
    for region in regions_dict.keys():
        real = bed.intersect(b = target, u = True, s=s, S=S, f = ov_fraction).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region]), u = True, s = s, S=S, f = ov_fraction)
        real_count = len(real)
        real_overlap = bed.intersect(b = target, u = True, s=s, S=S, f = ov_fraction)
    
        print("There are {0} {1} overlapping elements located in {2}".format(real_count, strandness, region))
        random_count = []

        for i in range(int(randomization)):
            input_bed_file = bed.intersect(b = target, u = True, s=s, S=S, f = ov_fraction).fn
            # Read the input BED file and calculate the mean length and chromosomal frequency
            existing_intervals, mean_length, chrom_count = read_bed_file(input_bed_file)

            # Read the exclude regions from the third BED file
            if exclude_intervals is not None:
                exclude_intervals_bed = read_exclude_bed_file(exclude_intervals)
            else:
                exclude_intervals_bed = None


            # Generate the same number of random intervals as the input file, using the chromosomal frequencies
            num_random_intervals = len(existing_intervals)
            random_intervals = generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed)

            # Write the random intervals to a new BED file
            write_bed_file(args.tmp + "/" + "Random_0123456789.bed", random_intervals)

            background_df_random = args.tmp + "/" + "Random_0123456789.bed"
            random_count.append(len(pybedtools.BedTool(background_df_random).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region]), u = True, s = s, S=S, f = ov_fraction)))


        zscore = compute_z_score(real_count,random_count)
        pvalue = compute_pvalue(zscore)
        table = save_tables(real_count, random_count, strandness, os.path.basename(target.fn) + "||" + region, table)
        frame = save_results(zscore,strandness, pvalue, os.path.basename(target.fn) + "||" + region, real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table

#Check if input path is a bed file
def is_bed_file(file_path):
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"Error: The file '{file_path}' does not exist.")
        return False
    
    # Open the file and check its contents
    with open(file_path, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            # Skip empty lines
            line = line.strip()
            if not line:
                continue
            
            # Split the line by tabs
            columns = line.split('\t')
            
            # Check if the line has exactly 6 or 12 columns
            if len(columns) != 6 and len(columns) != 12 :
                print(f"Error: Line {line_number} does not have exactly 6 or 12 columns.")
                return False
            
            # Validate the start and end positions (second and third columns)
            try:
                start = int(columns[1])
                end = int(columns[2])
            except ValueError:
                print(f"Error: Start and end positions are not integers at line {line_number}.")
                return False
            
            # Check that start < end
            if start >= end:
                print(f"Error: Start position is greater than or equal to end position at line {line_number}.")
                return False
            
            # Validate the name (fourth column) - it can be a string or "."
            name = columns[3]
            if name != '.' and not isinstance(name, str):
                print(f"Error: Invalid name format at line {line_number}. Name should be a non-empty string or '.' (dot).")
                return False
            
            # Validate the score (fifth column) - it can be a number or "."
            score = columns[4]
            if score != '.' and not is_number(score):
                print(f"Error: Score is not a valid number or '.' at line {line_number}.")
                return False
            
            # Validate the strand (sixth column) - it should be '+' or '-' or "."
            strand = columns[5]
            if strand not in ('+', '-', '.'):
                print(f"Error: Strand is not valid at line {line_number}. Expected '+', '-' or '.' (dot).")
                return False
    
    print(f"The file '{file_path}' appears to be a valid 6 or 12-column BED file.")
    return True

#Check if input is a number
def is_number(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

#Compute Zscore
def compute_z_score(real,random_counts):
    real = float(real)
    mean = statistics.mean(random_counts)
    stdev = statistics.stdev(random_counts)
    if stdev == 0:
        print("Standard deviation equal to zero, check files")
        zscore = "Not computable"
    else:
        zscore = (real - mean) / stdev
    return zscore

#Define function to compute p.value from Z-score
def compute_pvalue(zscore):
    if zscore == "Not computable":
        return("Not computable")
    else:
        pvalue = 2 * stats.norm.sf(abs(zscore))*2
        return(pvalue)

#Test length
def test_length(bed,background,randomization, target, frame, table, genome, exclude_intervals = None):
    if exclude_intervals is not None:
        nucleotide_content = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v = True).nucleotide_content(fi = genome).to_dataframe()
        real_length = statistics.mean(nucleotide_content["15_seq_len"].to_list())
        random_count_length = []
        background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals), v = True).to_dataframe()
    
    else:
        nucleotide_content = bed.nucleotide_content(fi = genome).to_dataframe()
        real_length = statistics.mean(nucleotide_content["15_seq_len"].to_list())
        random_count_length = []
        background_df = background.to_dataframe()
        
    for i in range(int(randomization)):
        random_features_index = random.sample(range(background_df.shape[0]), len(bed))
        background_df_random = background_df.iloc[random_features_index,:]
        random_count_length.append(statistics.mean(pybedtools.BedTool.from_dataframe(background_df_random).nucleotide_content(fi = genome, s = True).to_dataframe()["15_seq_len"]))

    zscore = compute_z_score(real_length,random_count_length)
    pvalue = compute_pvalue(zscore)

    frame = save_results(zscore, "Length", pvalue, os.path.basename(target), real_length, statistics.mean(random_count_length), statistics.stdev(random_count_length), frame)
    table = save_tables(real_length, random_count_length, "Length", os.path.basename(target), table)
    return frame, table

#Test AT and GC content
def test_GC_AT(bed,background,randomization, target, frame, table, genome, exclude_intervals = None):
    if exclude_intervals is not None:
        nucleotide_content = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v = True).nucleotide_content(fi = genome, s = True).to_dataframe()
        real_GC = statistics.mean(nucleotide_content["8_pct_gc"].to_list())
        real_AT = statistics.mean(nucleotide_content["7_pct_at"].to_list())
        random_count_GC = []
        random_count_AT = []
        background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals), v = True).to_dataframe()
      
    else:
        nucleotide_content = bed.nucleotide_content(fi = genome, s = True).to_dataframe()
        real_GC = statistics.mean(nucleotide_content["8_pct_gc"].to_list())
        real_AT = statistics.mean(nucleotide_content["7_pct_at"].to_list())
        random_count_GC = []
        random_count_AT = []
        background_df = background.to_dataframe()
    
    for i in range(int(randomization)):
        random_features_index = random.sample(range(background_df.shape[0]), len(bed))
        background_df_random = background_df.iloc[random_features_index,:]
        random_count_GC.append(statistics.mean(pybedtools.BedTool.from_dataframe(background_df_random).nucleotide_content(fi = genome, s = True).to_dataframe()["8_pct_gc"]))
        random_count_AT.append(statistics.mean(pybedtools.BedTool.from_dataframe(background_df_random).nucleotide_content(fi = genome, s = True).to_dataframe()["7_pct_at"]))
    
    zscore_GC = compute_z_score(real_GC,random_count_GC)
    pvalue_GC = compute_pvalue(zscore_GC)
    zscore_AT = compute_z_score(real_AT,random_count_AT)
    pvalue_AT = compute_pvalue(zscore_AT)
    
    frame = save_results(zscore_GC, "GC", pvalue_GC, os.path.basename(target), real_GC, statistics.mean(random_count_GC), statistics.stdev(random_count_GC), frame)
    frame = save_results(zscore_AT, "AT", pvalue_AT, os.path.basename(target), real_AT, statistics.mean(random_count_AT), statistics.stdev(random_count_AT), frame)
    table = save_tables(real_GC, random_count_GC, "GC", os.path.basename(target), table)
    table = save_tables(real_AT, random_count_AT, "AT", os.path.basename(target), table)
    return frame, table

#Test overlap
def test_enrichement(bed,target,background,ov_fraction, randomization, frame, table, strandness="strandless",exclude_intervals = None):
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    if exclude_intervals is not None:
        real = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b = target, u = True, s = s, S=S, f = ov_fraction)
        real_count = len(real)
        print("There are {0} {1} overlapping elements".format(real_count,strandness))
        random_count = []
        background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals),v =True, s=s, S=S).to_dataframe()
    
    else:
        real = bed.intersect(b = target, u = True, s = s, S=S, f = ov_fraction)
        real_count = len(real)
        print("There are {0} {1} overlapping elements".format(real_count, strandness))
        random_count = []
        background_df = background.to_dataframe()
        
    for i in range(int(randomization)):
        random_features_index = random.sample(range(background_df.shape[0]), len(bed))
        background_df_random = background_df.iloc[random_features_index,:]
        random_count.append(len(pybedtools.BedTool.from_dataframe(background_df_random).intersect(b = target, u = True, s = s,S=S, f = ov_fraction)))
    zscore = compute_z_score(real_count,random_count)
    pvalue = compute_pvalue(zscore)
    table = save_tables(real_count, random_count, strandness, os.path.basename(target.fn), table)
    frame = save_results(zscore,strandness,pvalue, os.path.basename(target.fn), real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table


def test_enrichement_regions(bed,target,background,ov_fraction, randomization, frame, table, strandness="strandless",exclude_intervals = None, regions_bed = "NA"):
    
    #read the bed file with annotations
    regions_bed = pd.read_table(args.tmp + "/" + "Reference_regions.bed")
    #Set negative regions to 0, this is bevcause some promoters may be calculated and negative start if the chromosome is too short
    regions_bed.iloc[:, 1] = regions_bed.iloc[:, 1].clip(lower=0)
    # Split the dataframe based on the values in column number 4, i.e name 
    regions_dict = {category: group for category, group in regions_bed.groupby(regions_bed.columns[3])}

    #get strand parameters
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    #Testing all regions
    for region in regions_dict.keys():
        if exclude_intervals is not None:
            real = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b = target, u = True, s = s, S=S, f = ov_fraction).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region]), u = True, s = s, S=S, f = ov_fraction)
            real_count = len(real)
            real_overlap = len(bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b = target, u = True, s = s, S=S, f = ov_fraction))
            print("There are {0} {1} overlapping elements located in {2}".format(real_count,strandness,region))
            random_count = []
            background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals),v =True, s=s, S=S).intersect(b = target, u = True, s = s, S=S, f = ov_fraction).to_dataframe()

        else:
            real = bed.intersect(b = target, u = True, s = s, S=S, f = ov_fraction).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region]), u = True, s = s, S=S, f = ov_fraction)
            real_count = len(real)
            real_overlap = len(bed.intersect(b = target, u = True, s = s, S=S, f = ov_fraction))
            print("There are {0} {1} overlapping elements located in {2}".format(real_count, strandness, region))
            random_count = []
            background_df = background.intersect(b = target, u = True, s = s, S=S, f = ov_fraction).to_dataframe()
            print("number of feature in bk is " + str(background_df.shape[0]) + "and bed is " + str(real_overlap))

        for i in range(int(randomization)):
            random_features_index = random.sample(range(background_df.shape[0]), real_overlap)
            background_df_random = background_df.iloc[random_features_index,:]
            random_count.append(len(pybedtools.BedTool.from_dataframe(background_df_random).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region]), u = True, s = s, S=S, f = ov_fraction) ))
        zscore = compute_z_score(real_count,random_count)
        pvalue = compute_pvalue(zscore)
        table = save_tables(real_count, random_count, strandness, os.path.basename(target.fn) + "||" + region , table)
        frame = save_results(zscore,strandness,pvalue, os.path.basename(target.fn) + "||" + region , real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
        print("tables saved")
    return frame, table



#Functions for Saving results and Tables for plotting
def save_tables(real,random, name, target, table):
    tmp = pd.concat([pd.DataFrame({"Name": name, "Target": target, "Type" : "Random" , "Count": random}, index=range(0,len(random))),pd.DataFrame({"Name": name, "Target": target,"Type" : "Real" , "Count": real}, index=[0])]).reset_index(drop=True)
    table = pd.concat([table, tmp])
    return table
def save_results(zscore, name, pvalue, target, real, random, stdev, frame):
    results = pd.concat([frame, pd.DataFrame({"Zscore" : zscore, "Type" : name, "P.value" : pvalue, "Target": target, "Real": real, "Random": random, "Sd" : stdev }, index=[0])])
    return results

#Function to read 6 columns BED file, name score and strand can be placeholder "."
def read_bed_file(file_path):
    intervals = []
    total_length = 0 
    chrom_count = Counter() 
    with open(file_path, 'r') as f:
        for line in f:
            # Skip comments and empty lines
            if line.startswith("#") or not line.strip():
                continue
            # Split line into chrom, start, end, name, score, strand
            parts = line.strip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) > 3 else "."
            score = parts[4] 
            strand = parts[5]

            # Calculate total length for mean length calculation
            length = end - start
            total_length += length

            intervals.append((chrom, start, end, name, score, strand, length))
            chrom_count[chrom] += 1  # Count occurrences of each chromosome

    # Calculate the mean length of intervals
    mean_length = total_length / len(intervals) if intervals else 0
    return intervals, mean_length, chrom_count

#Function to read the exclude BED file and store intervals to exclude
def read_exclude_bed_file(file_path):
    exclude_intervals = {}
    with open(file_path, 'r') as f:
        for line in f:
            # Skip comments and empty lines
            if line.startswith("#") or not line.strip():
                continue
            # Split line into chrom, start, end
            parts = line.strip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])

            if chrom not in exclude_intervals:
                exclude_intervals[chrom] = []
            exclude_intervals[chrom].append((start, end))
    return exclude_intervals

#Function to check if a random interval overlaps with any exclusion region
def is_overlapping(start, end, exclude_intervals):
    for ex_start, ex_end in exclude_intervals:
        if (start < ex_end) and (end > ex_start):  # Check for overlap
            return True
    return False

#Function to generate random intervals based on existing BED intervals
def generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals=None):
    random_intervals = []
    total_intervals = sum(chrom_count.values())  # Total number of intervals
    
    # Calculate min and max start and end positions for each chromosome
    chrom_ranges = {}
    for chrom, start, end, name, score, strand, length in existing_intervals:
        if chrom not in chrom_ranges:
            chrom_ranges[chrom] = {'min': start, 'max': end}
        else:
            chrom_ranges[chrom]['min'] = min(chrom_ranges[chrom]['min'], start)
            chrom_ranges[chrom]['max'] = max(chrom_ranges[chrom]['max'], end)
    
    # For each chromosome, generate intervals proportionally to its count
    for chrom, count in chrom_count.items():
        # Calculate the number of random intervals for this chromosome
        chrom_intervals = round((count / total_intervals) * num_random_intervals)
        
        # Define the range of possible positions on this chromosome
        chrom_range = chrom_ranges[chrom]
        chrom_min = chrom_range['min']
        chrom_max = chrom_range['max']
        
        for _ in range(chrom_intervals):
            # Keep trying to generate a valid interval
            while True:  
                # Generate a random start position within the chromosome's range
                start = random.randint(chrom_min, chrom_max)
                
                # Randomly select a length close to the mean length
                length = abs(random.gauss(mean_length, mean_length * 0.1))  # 10% variation from mean length
                length = max(1, int(length))  # Ensure the length is at least 1
                
                # Ensure the generated interval respects the chromosome's bounds
                end = min(start + length, chrom_max)
                
                # Make sure that start < end
                if end <= start:
                    continue  # Skip if the interval is invalid
                
                # Check if the generated interval overlaps with any exclusion region (if exclude_intervals is provided)
                if exclude_intervals and chrom in exclude_intervals:
                    # Exclude the interval if it overlaps with any of the exclude intervals
                    if any(start < ex_end and end > ex_start for ex_start, ex_end in exclude_intervals[chrom]):
                        continue  # Skip this interval if it overlaps
                
                # Generate a random name (could be something like 'random1', 'random2', etc.)
                name = f"random_{random.randint(1, 10000000)}"
                
                # Generate a random score between 0 and 1000 
                score = random.randint(0, 1000)
                
                # Randomly pick a strand ('+' or '-')
                strand = random.choice(["+", "-"])
                
                # Store the random interval with name, score, and strand
                random_intervals.append((chrom, start, end, name, score, strand))
                break  # Exit the loop once a valid interval is generated
    
    return random_intervals

#Function to write intervals to a 6 column BED file
def write_bed_file(file_path, intervals):
    with open(file_path, 'w') as f:
        for chrom, start, end, name, score, strand in intervals:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")


# Function to compute enrichment score (ES)
def compute_enrichment_score(ranked_list, gene_set):
    N = len(ranked_list)  # Total number of regions
    Nh = len(gene_set)  # Number of regions in the intersected set 
    Nm = N - Nh  # Number of regions NOT in the intersected set

    # Get indices of regions in intersected set
    hit_indices = [i for i, gene in enumerate(ranked_list) if gene in gene_set]

    if not hit_indices:
        return 0, [], []  # Return 0 if no overlap

    # Compute cumulative sums for hit and miss genes
    p_hit = np.cumsum([1/Nh if i in hit_indices else 0 for i in range(N)])
    p_miss = np.cumsum([1/Nm if i not in hit_indices else 0 for i in range(N)])

    # Compute running enrichment score
    ES_values = p_hit - p_miss
    ES = max(ES_values, key=abs)  # Maximum absolute deviation is ES
    peak_idx = np.argmax(np.abs(ES_values))  #Index where ES is max
    leading_regions = ranked_list[:peak_idx + 1]

    leading_regions_str = ",".join(leading_regions)

    return ES, ES_values, hit_indices, leading_regions_str


# Function to run permutation test
def permutation_test(ranked_list, gene_set, num_permutations=100):
    real_ES, _, _ , leading_regions = compute_enrichment_score(ranked_list, gene_set)

    # Generate null distribution by shuffling gene ranks
    null_ES = []
    for _ in range(num_permutations):
        shuffled_genes = np.random.permutation(ranked_list)
        shuffled_ES, _, _, _ = compute_enrichment_score(shuffled_genes, gene_set)
        null_ES.append(shuffled_ES)

    null_ES = np.array(null_ES)

    # Compute p-value
    p_value = np.sum(np.abs(null_ES) >= np.abs(real_ES)) / num_permutations

    # Compute Z-score
    null_mean = np.mean(null_ES)
    null_std = np.std(null_ES)
    Z_score = (real_ES - null_mean) / null_std if null_std != 0 else 0

    return real_ES, p_value, Z_score, null_ES, leading_regions

def scale_value(value, min_range, max_range, new_min, new_max):
    # Apply the Min-Max scaling formula
    scaled_value = ((value - min_range) / (max_range - min_range)) * (new_max - new_min) + new_min
    return scaled_value

#Weight the BED score field using other value, better if value are min max scaled between 0 and 1
def weighted_score(x, y, alpha=0.1):
                N = len(x) 
                rank_normalized = (np.argsort(x) + 1) / N 
                return [xi * (1 - alpha * (yi / rank)) for xi, yi, rank in zip(x, y, rank_normalized)]


def read_and_rank_bed(file_path, ascending=False):
    bed_df = pd.read_csv(file_path, sep="\t", header=None)
    
    bed_df = bed_df.iloc[:, :6]
    bed_df.columns = ["chrom", "start", "end", "name", "score", "strand"]
    
    # Convert score column to numeric and check for invalid values
    bed_df["score"] = pd.to_numeric(bed_df["score"], errors="coerce")
    if bed_df["score"].isna().any():
        print("The 'score' column contains non-numeric or NA or empty values, this analysis require all regions to have a numeric score")
        #sys.exit()
    
    # Generate names for regions where 'name' is "." or empty
    bed_df["name"] = bed_df["name"].replace({".": None})
    
    ranked_bed = bed_df.sort_values(by="score", ascending=ascending)
    bed_df["name"] = bed_df["name"].fillna(pd.Series([f"Region_{i+1}" for i in range(len(bed_df))], index=bed_df.index))
    #ranked_bed["score"] = [i for i in range(0,ranked_bed.shape[0])]
    return ranked_bed

#Test overlap
def test_enrichement_rank(bed,target,ov_fraction, randomization, strandness="strandless",exclude_intervals = None, ascending = False, WeightRanking = False, alpha = 0.1):
    #Get strand parameters
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]

    NAME = os.path.basename(bed)
    TARGET = os.path.basename(target)

    #Read Input and Target Bed files
    bed = pybedtools.BedTool.from_dataframe(read_and_rank_bed(bed, ascending=ascending ))
    target = pybedtools.BedTool.from_dataframe(get_bed_file(target))    
   
    if exclude_intervals is not None:
        #Perform intersection for overlapping regions, excluding some regions
        overlap_df = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b=target, s=s, S=S, f=ov_fraction, wo=True).to_dataframe(header=None)
        intersected_regions = overlap_df.iloc[:,3].tolist()
        #Perform intersection for non-overlapping regions
        non_overlap_df = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b=target, s=s, S=S, f=ov_fraction, v=True).to_dataframe(header=None)
        non_overlap_df.columns = range(0,non_overlap_df.shape[1])
    else:
        #Perform intersection for overlapping regions
        overlap_df = bed.intersect(b=target, s=s, S=S, f=ov_fraction, wo=True).to_dataframe(header=None)
        intersected_regions = overlap_df.iloc[:,3].tolist()
        #Perform intersection for non-overlapping regions
        non_overlap_df = bed.intersect(b=target, s=s, S=S, f=ov_fraction, v=True).to_dataframe(header=None)
        non_overlap_df.columns = range(0,non_overlap_df.shape[1])
            
            
    # Bind the two dataframes: overlap and non-overlap
    combined_df = pd.concat([overlap_df, non_overlap_df], ignore_index=True)
            
            
    # For non-overlapping regions, set overlap length (last column) to 0
    combined_df["overlap_length"] = combined_df.iloc[:, -1].fillna(0)
    combined_df["length"] = combined_df.iloc[:,2] - combined_df.iloc[:,1]
    combined_df["overlap_fraction"] = combined_df["overlap_length"] / combined_df["length"]
    combined_df["Zscore_overlap_fraction"] = (combined_df["overlap_fraction"] - np.mean(combined_df["overlap_fraction"])) / np.std(combined_df["overlap_length"])
    combined_df["Zscore_score"] = (combined_df.iloc[:, 4] - np.mean(combined_df.iloc[:, 4])) / np.std(combined_df.iloc[:, 4])
        
    #Set ranges for scaling
    min_range = combined_df["Zscore_score"].min()
    max_range = combined_df["Zscore_score"].max()
    new_min = 0
    new_max = 1
            
    combined_df["scaled_Zscore_score"] = combined_df["Zscore_score"].apply(lambda x: scale_value(x, min_range, max_range, new_min, new_max))
    combined_df["scaled_overlap_length"] = combined_df["overlap_length"].apply(lambda x: scale_value(x, min_range, max_range, new_min, new_max))
    
    if WeightRanking == True:
        # Sort by the weighted score
        combined_df["weighted_score"] = weighted_score(combined_df["scaled_Zscore_score"],combined_df["overlap_fraction"], alpha)
        regions_ranked = combined_df.sort_values(by="weighted_score", ascending=ascending)  
    if WeightRanking == False:
        # Sort by the score
        regions_ranked = combined_df.sort_values(by="Zscore_score", ascending=ascending)    
    
    # Extract the region names (from the 3rd column, which contains the region names)
    regions_ranked = regions_ranked.iloc[:, 3].tolist()

    
    ES, p_value, Zscore, null_ES, leading_regions = permutation_test(regions_ranked, intersected_regions, num_permutations=randomization)
    
    # Compute ES for plotting
    ES, ES_values, hit_indices, _ = compute_enrichment_score(regions_ranked, intersected_regions)

    res = pd.DataFrame({"Input" : NAME, "Target" : TARGET, "ES" : ES, "Pval" : p_value, "Type": "Intersect", "NRand" : str(randomization), "Orientation" : strandness, "Leading" : leading_regions}, index=[0])
    
    forplot = pd.DataFrame({"ES_values" : "|".join(map(str, ES_values)), "hit_indices" : "|".join(map(str, hit_indices)), "Null_ES": "|".join(map(str, null_ES)), "ES": str(ES)}, index=[0])
    
    return res, forplot


def RankPlot(forplot):
    
    ES_values = list(map(float, forplot["ES_values"][0].split("|")))
    hit_indices = list(map(float, forplot["hit_indices"][0].split("|")))
    null_ES = list(map(float, forplot["Null_ES"][0].split("|")))
    ES = float(forplot["ES"][0])

    plt.figure(figsize=(10, 5))
    plt.plot(range(len(ES_values)), ES_values, color="blue", linewidth=2, label="Running ES")
    plt.axhline(0, color="black", linestyle="--", linewidth=1)

    # Mark peak ES point
    peak_idx = np.argmax(np.abs(ES_values))
    
    plt.scatter(peak_idx, ES_values[peak_idx], color="red", marker="*", s=100, label=f"Peak ES = {ES:.3f}")

    # Add vertical lines where gene set genes are located
    for idx in hit_indices:
       plt.axvline(x=idx, color="gray", linestyle="dotted", alpha=0.5)

    plt.xlabel("Ranked Gene List")
    plt.ylabel("Enrichment Score (ES)")
    plt.title("Rank Enrichment Plot")
    plt.legend()
    plt.show() 
    
    plt.figure(figsize=(8, 5))
    plt.hist(null_ES, bins=10, alpha=0.7, color="gray", label="Null Distribution", density=True)
    plt.axvline(ES, color="red", linestyle="--", label=f"Real ES = {ES:.3f}")
    plt.xlabel("Enrichment Score (ES)")
    plt.ylabel("Frequency")
    plt.legend()
    plt.show()

#Test closest
def KS_test_closest(x, y, n_shuffles=100):
    
    # Calcolare la distribuzione cumulativa reale
    cumulative_real = np.cumsum(y) / np.sum(y)
    

    # Inizializzare liste per shuffle
    p_values = []
    enrichment_scores = []
    shuffle_cdfs = []

    # Effettuare n_shuffles per generare distribuzioni randomiche
    for _ in range(n_shuffles):
        y_shuffled = np.random.permutation(y)  # Shuffle delle etichette
        cumulative_shuffle = np.cumsum(y_shuffled) / np.sum(y_shuffled)
        
        
        # Test KS tra la distribuzione reale e quella shuffleata
        ks_stat, p_value = ks_2samp(cumulative_real, cumulative_shuffle)
        p_values.append(p_value)

        # Score di arricchimento (differenza massima tra le distribuzioni cumulative)
        max_index = np.argmax(abs(cumulative_real - cumulative_shuffle))
        enrichment_score =  (cumulative_real - cumulative_shuffle)[max_index]
        enrichment_scores.append(enrichment_score)

        # Salvare la distribuzione shuffleata
        shuffle_cdfs.append(cumulative_shuffle)

    # Calcolare la media delle distribuzioni shuffleate
    mean_shuffle_cdf = np.mean(shuffle_cdfs, axis=0)

    # Calcolare il p-value medio e lo score medio di arricchimento
    p_value_mean = np.mean(p_values)
    enrichment_score_max = enrichment_scores[np.argmax([abs(x) for x in enrichment_scores])]
    leading = x[0:np.argmax([abs(x) for x in enrichment_scores])]

    leading_str = ",".join(leading)

    return enrichment_score_max, p_value_mean, leading_str, cumulative_real, mean_shuffle_cdf


def RankPlot_closest(forplot):

    cumulative_real = list(map(float, forplot["cumulative_real"][0].split("|")))
    mean_shuffle_cdf = list(map(float, forplot["mean_shuffle_cdf"][0].split("|")))
    x = [x for x in range(len(cumulative_real))]

    # Plot delle distribuzioni cumulative
    plt.figure(figsize=(10,6))
    plt.plot(x, cumulative_real, label='Cumulativa Reale', color='blue', linewidth=2)
    plt.plot(x, mean_shuffle_cdf, label='Cumulativa Shuffleata Media', color='red', alpha=0.6)
    plt.fill_between(x, cumulative_real, mean_shuffle_cdf, color='gray', alpha=0.3)
    plt.title("Distribuzioni Cumulative: Reale vs Shuffle Media")
    plt.legend()
    plt.show()

def test_closest_rank(bed,target,ov_fraction, randomization, strandness="strandless", exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False, ascending = False, WeightRanking = False, alpha = 0.1):
    #Get strand parameters
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]

    NAME = os.path.basename(bed)
    TARGET = os.path.basename(target)

    #Read Input and Target Bed files
    bed = pybedtools.BedTool.from_dataframe(read_and_rank_bed(bed, ascending=ascending))
    target = pybedtools.BedTool.from_dataframe(get_bed_file(target))

    if exclude_intervals is not None:
        closest_df = bed.sort().intersect(b = pybedtools.BedTool(exclude_intervals),s = s, S = S, v = True).closest(b = target.sort().intersect(b = pybedtools.BedTool(exclude_intervals),s = s, S = S, v = True), d=True,s=s,S=S,f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe(header = None)
    else:
        closest_df = bed.sort().closest(b = target.sort(), d=True,s=s,S=S,f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe(header = None)
            
    closest_df = bed.sort().closest(b = target.sort(), d=True,s=s,S=S,f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe(header = None)
    closest_df["Distance"] = abs(closest_df.iloc[:,12])        
    closest_df["Zscore_Distance"] = (closest_df["Distance"] - np.mean(closest_df["Distance"])) / np.std(closest_df["Distance"])
    closest_df["Score"] = closest_df.iloc[:, 4]
    closest_df["Zscore_score"] = (closest_df.iloc[:, 4] - np.mean(closest_df.iloc[:, 4])) / np.std(closest_df.iloc[:, 4])

    min_range = closest_df["Zscore_score"].min()
    max_range = closest_df["Zscore_score"].max()
    new_min = 0
    new_max = 1
            

    closest_df["scaled_Zscore_score"] = closest_df["Zscore_score"].apply(lambda x: scale_value(x, min_range, max_range, new_min, new_max))
    closest_df["scaled_Distance"] = closest_df["Distance"].apply(lambda x: scale_value(x, min_range, max_range, new_min, new_max))
    
    if WeightRanking == False:
        #Sort by the  score in descending order
        regions_ranked = closest_df.sort_values(by="Score", ascending=ascending)
    if WeightRanking == True:
        closest_df["weighted_score"] = weighted_score(closest_df["scaled_Zscore_score"],closest_df["scaled_Distance"], alpha)
        # Sort by the weighted score in descending order
        regions_ranked = closest_df.sort_values(by="weighted_score", ascending=ascending)
    
    regions_ranked["Rank"] = [x for x in range(1,regions_ranked.shape[0] + 1 )] 
    ES, Pval, leading, cumulative_real, mean_shuffle_cdf = KS_test_closest(regions_ranked.iloc[:,3], regions_ranked["Distance"], randomization)
    
    res = pd.DataFrame({"Input": NAME, "Target": TARGET, "ES": ES, "Pval" : Pval, "Type" : "Closest", "NRand": str(randomization), "Orientation": strandness, "Leading": leading}, index = [0])
    forplot = pd.DataFrame({"cumulative_real": "|".join(map(str, cumulative_real)) , "mean_shuffle_cdf": "|".join(map(str, mean_shuffle_cdf))}, index = [0])
    
    return res, forplot


###### MAIN #########
def main(mode,input,targets,background,orientation,genome,ov_fraction,randomization,outfile, exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False, RankTest = False, WeightRanking = False, alpha= 0.1):
    #Check imput parameters
    check_input_parameters(mode)

    #Check and create tmp directory
    create_directory(args.tmp)

    #Set tmp dir for pybedtools
    pybedtools.set_tempdir(args.tmp)
    
    if args.GenomicLocalization:
        if args.gtf is not None:
            subprocess.run(["Rscript", "Create_bed_genomicRegions.R", args.gtf, args.tmp, args.genome])
            regions_bed = args.tmp + "/" + "Reference_regions.bed"
        if args.bed is not None:
            regions_bed = args.bed

        
    if args.generate_bg == False and mode == "intersect" and RankTest == False:
        print("Running intersect mode using provided background as " + background)
        #Check input files
        beds_targets = targets.split(",")
        beds = [input, background]
        beds.extend(beds_targets)
        for bed in beds:
            if is_bed_file(bed):
                print("")
            else:
                sys.exit()

        #Load bedfile and bakground file
        bedfile = pybedtools.BedTool.from_dataframe(get_bed_file(input))
        backgroundfile = pybedtools.BedTool.from_dataframe(get_bed_file(background))
        orientations = orientation.split(",")
    
        #Create empy df to store results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation  in orientations:
            for target in targets.split(","):
                targetfile = pybedtools.BedTool.from_dataframe(get_bed_file(target))
                res = test_enrichement(bedfile,targetfile,backgroundfile,ov_fraction,randomization, frame, table, strandness = orientation, exclude_intervals = exclude_intervals)
                frame = res[0]
                table = res[1]
                if args.GenomicLocalization:
                    res = test_enrichement_regions(bedfile,targetfile,backgroundfile,ov_fraction,randomization, frame, table, strandness = orientation, exclude_intervals = exclude_intervals, regions_bed = regions_bed)
                    frame = res[0]
                    table = res[1]

        if args.test_AT_GC:
            for target in targets.split(","):
                res = test_GC_AT(bedfile,backgroundfile,randomization,target, frame, table, genome, exclude_intervals = exclude_intervals)
                frame = res[0]
                table = res[1]
            
        if args.test_length:
            for target in targets.split(","):
                res = test_length(bedfile,backgroundfile,randomization, target, frame, table,genome, exclude_intervals =  exclude_intervals)
                frame = res[0]
                table = res[1]
            
        
            
        frame.to_csv(outfile, index=False, sep = "\t")
        table.to_csv("Tables.txt", index=False, sep = "\t")
    
    if args.generate_bg == False and mode == "closest" and RankTest == False:
        print("Running closest mode using provided background as " + background)
        #Check input files
        beds_targets = targets.split(",")
        beds = [input, background]
        beds.extend(beds_targets)
        for bed in beds:
            if is_bed_file(bed):
                print("")
            else:
                sys.exit()

        #Load bedfile and bakground file
        bedfile = pybedtools.BedTool.from_dataframe(get_bed_file(input))
        backgroundfile = pybedtools.BedTool.from_dataframe(get_bed_file(background))
        orientations = orientation.split(",")
    
        #Create empy df to stoire results
        frame = pd.DataFrame()
        table = pd.DataFrame()
        
        for orientation in orientations:
            for target in targets.split(","):
                targetfile = pybedtools.BedTool.from_dataframe(get_bed_file(target))
                res = test_closeness(bedfile, targetfile, backgroundfile, ov_fraction, randomization, frame, table,strandness=orientation, exclude_intervals = exclude_intervals, exclude_ov = exclude_ov, exclude_upstream = exclude_upstream, exclude_downstream=exclude_downstream)
                frame = res[0]
                table = res[1]
            
    
        if args.test_AT_GC:
            for target in targets.split(","):
                res = test_GC_AT(bedfile,backgroundfile,randomization,target, frame, table, genome, exclude_intervals = exclude_intervals)
                frame = res[0]
                table = res[1]
            
        if args.test_length in orientations:
            for target in targets.split(","):
                res = test_length(bedfile,backgroundfile,randomization, target, frame, table, genome, exclude_intervals = exclude_intervals)
                frame = res[0]
                table = res[1]
            
        
            
        frame.to_csv(outfile, index=False, sep = "\t")
        table.to_csv("Tables.txt", index=False, sep = "\t")

    
    if args.generate_bg == True and mode == "intersect" and RankTest == False:
        print("Running intersect mode using a random generated background")
        #Check input files
        beds_targets = targets.split(",")
        beds = [input]
        beds.extend(beds_targets)
        for bed in beds:
            if is_bed_file(bed):
                print("")
            else:
                sys.exit()

        #Load bedfile and bakground file
        bedfile = pybedtools.BedTool.from_dataframe(get_bed_file(input))
        orientations = orientation.split(",")
    
        #Create empy df to stoire results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation in orientations:
            for target in targets.split(","):
                targetfile = pybedtools.BedTool.from_dataframe(get_bed_file(target))
                res = test_enrichement_random_bg(bedfile,targetfile,ov_fraction,randomization, frame, table, strandness= orientation, exclude_intervals = exclude_intervals)
                frame = res[0]
                table = res[1]
                if args.GenomicLocalization:
                    res = test_enrichement_regions_random_bg(bedfile,targetfile,ov_fraction,randomization, frame, table, strandness= orientation, exclude_intervals = exclude_intervals, regions_bed = regions_bed)
                    frame = res[0]
                    table = res[1]
    
        if args.test_AT_GC:
            for target in targets.split(","):
                res = test_GC_AT_random_bg(bedfile,randomization,target, frame, table, genome, exclude_intervals = exclude_intervals)
                frame = res[0]
                table = res[1]
            
        if args.test_length:
            for target in targets.split(","):
                res = test_length_random_bg(bedfile,randomization, target, frame, table, genome, exclude_intervals = exclude_intervals)
                frame = res[0]
                table = res[1]
            
        
            
        frame.to_csv(outfile, index=False, sep = "\t")
        table.to_csv("Tables.txt", index=False, sep = "\t")
    
    if args.generate_bg == True and mode == "closest" and RankTest == False:
        print("Running closest mode using a random generated background")
        #Check input files
        beds_targets = targets.split(",")
        beds = [input]
        beds.extend(beds_targets)
        for bed in beds:
            if is_bed_file(bed):
                print("")
            else:
                sys.exit()

        #Load bedfile and bakground file
        bedfile = pybedtools.BedTool.from_dataframe(get_bed_file(input))
        orientations = orientation.split(",")
    
        #Create empy df to stoire results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation in orientations:
            for target in targets.split(","):
                targetfile = pybedtools.BedTool.from_dataframe(get_bed_file(target))
                res = test_closeness_random_bg(bedfile, targetfile, ov_fraction, randomization, frame, table,strandness=orientation, exclude_intervals = exclude_intervals, exclude_ov = exclude_ov, exclude_upstream = exclude_upstream, exclude_downstream=exclude_downstream)
                frame = res[0]
                table = res[1]
    
        if args.test_AT_GC:
            for target in targets.split(","):
                res = test_GC_AT_random_bg(bedfile,randomization,target, frame, table,genome, exclude_intervals = exclude_intervals)
                frame = res[0]
                table = res[1]
            
        if args.test_length:
            for target in targets.split(","):
                res = test_length_random_bg(bedfile,randomization, target, frame, table,genome,  exclude_intervals = exclude_intervals)
                frame = res[0]
                table = res[1]

        frame.to_csv(outfile, index=False, sep = "\t")
        table.to_csv("Tables.txt", index=False, sep = "\t")
            
    if RankTest == True and mode == "intersect":
        print("Perform Rank Mode Test using Intersect")
        #Check input files
        beds_targets = targets.split(",")
        beds = [input, background]
        beds.extend(beds_targets)
        for bed in beds:
            if is_bed_file(bed):
                print("")
            else:
                sys.exit()

        #Load bedfile and bakground file
        orientations = orientation.split(",")
    
        #Create empy df to store results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation  in orientations:
            for target in targets.split(","):
                print("Running test")
                res, forplot = test_enrichement_rank(input,target,ov_fraction, randomization, strandness=orientation,exclude_intervals = exclude_intervals, ascending = args.Ascending_RankOrder, WeightRanking = WeightRanking, alpha = alpha)
                RankPlot(forplot)
                frame = pd.concat([frame, res])
        
        frame.to_csv(outfile, index=False, sep = "\t")

    if RankTest == True and mode == "closest":
        print("Perform Rank Mode Test using Intersect")
        #Check input files
        beds_targets = targets.split(",")
        beds = [input, background]
        beds.extend(beds_targets)
        for bed in beds:
            if is_bed_file(bed):
                print("")
            else:
                sys.exit()

        #Load bedfile and bakground file
        orientations = orientation.split(",")
    
        #Create empy df to store results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation  in orientations:
            for target in targets.split(","):
                res, forplot = test_closest_rank(input,target,ov_fraction, randomization,strandness="strandless", exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False, ascending = args.Ascending_RankOrder, WeightRanking = WeightRanking, alpha = alpha)
                RankPlot_closest(forplot)
                frame = pd.concat([frame, res])
        
        frame.to_csv(outfile, index=False, sep = "\t")
        #table.to_csv("Tables_RankClosest.txt", index=False, sep = "\t")
        

if __name__ == "__main__":
    main(args.mode,args.input,args.targets,args.background,args.orientation,args.genome,args.ov_fraction,args.randomization, args.outfile, args.exclude_intervals, args.exclude_ov, args.exclude_upstream, args.exclude_downstream, args.RankTest, args.WeightRanking, args.alpha)