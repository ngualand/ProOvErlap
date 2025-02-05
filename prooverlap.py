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

# Suppress all warnings
warnings.filterwarnings("ignore")
# Suppress specific warnings
warnings.filterwarnings("ignore", message="the index file is older than the FASTA file")

#Create parser
parser = argparse.ArgumentParser()
parser.add_argument("--mode", required = True, help = "Define mode: intersect or closest, intersect count the number of overlapping elements while closest test the distance. In closest if a feature overlap a target the distance is 0")
parser.add_argument("--input", required = True, help="Input bed file, must contain 6 column, name and score can be placeholder and are not used, strand is used only if some strandess test are requested")
parser.add_argument("--targets", required = True, help="Target bed files (must contain 6 columns), to test enrichement against, if multiple files are supplied N independent test against each file are conducted, file names must be comma separated, the name of the file will be use as the name output")
parser.add_argument("--background", help="Background bed file (must contain 6 columns) to test enrichement aginst, should be a superset from wich input bed file is derived")
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

#Parse argument
args = parser.parse_args()

#Functions

#Create directory if not exist, if exist print and exit
def create_directory(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Directory created: {folder_path}")
    else:
        print("tmp dir exists, please choose a non-existing directory, it will be created")
        print("Exit")
        sys.exit()

#Check input parameters
def check_input_parameters(mode):
    if mode == "intersect" or mode == "closest":
        print("Mode set to " + mode)
    else:
        print("Mode is not set, choose one of intersect or closest")
        print("Exit")
        sys.exit()

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

#Test length againsta a random background
def test_length_random_bg(bed,randomization, target, frame, table, genome, exclude_intervals = None):
    nucleotide_content = bed.nucleotide_content(fi = genome).to_dataframe()
    real_length = statistics.mean(nucleotide_content["15_seq_len"].to_list())
    random_count_length = []
    for i in range(int(randomization)):
        # Example usage:
        input_bed_file = args.input
        #output_bed_file = "random_intervals_with_fields.bed"
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
        #output_bed_file = "random_intervals_with_fields.bed"
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
            
            # Check if the line has exactly 6 columns
            if len(columns) != 6:
                print(f"Error: Line {line_number} does not have exactly 6 columns.")
                return False
            
            # Validate the chromosome (first column)
            chrom = columns[0]
            if not (chrom.startswith('chr') or chrom.isdigit()):
                print(f"Error: Invalid chromosome format at line {line_number}. Expected 'chr' prefix or numeric chromosome name.")
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
    
    print(f"The file '{file_path}' appears to be a valid 6-column BED file.")
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
        print("tbles saved")
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

            # Handle score: check if it's a valid integer or if it's a '.' (missing value)
            score = parts[4] if len(parts) > 4 else "."
            if score == ".":
                score = 0  
            else:
                try:
                    score = int(score)
                except ValueError:
                    score = 0  

            # Handle strand: check if it's a valid strand or use "+" as default
            strand = parts[5] if len(parts) > 5 else "+"

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
            while True:  # Keep trying to generate a valid interval
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
                name = f"random_{random.randint(1, 100000)}"
                
                # Generate a random score between 0 and 1000 (common BED score range)
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


###### MAIN #########
def main(mode,input,targets,background,orientation,genome,ov_fraction,randomization,outfile, exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False):
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

        
    if args.generate_bg == False and mode == "intersect":
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
        bedfile = pybedtools.BedTool(input)
        backgroundfile = pybedtools.BedTool(background)
        orientations = orientation.split(",")
    
        #Create empy df to store results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation  in orientations:
            for target in targets.split(","):
                targetfile = pybedtools.BedTool(target)
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
    
    if args.generate_bg == False and mode == "closest":
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
        bedfile = pybedtools.BedTool(input)
        backgroundfile = pybedtools.BedTool(background)
        orientations = orientation.split(",")
    
        #Create empy df to stoire results
        frame = pd.DataFrame()
        table = pd.DataFrame()
        
        for orientation in orientations:
            for target in targets.split(","):
                targetfile = pybedtools.BedTool(target)
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

    
    if args.generate_bg == True and mode == "intersect":
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
        bedfile = pybedtools.BedTool(input)
        orientations = orientation.split(",")
    
        #Create empy df to stoire results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation in orientations:
            for target in targets.split(","):
                targetfile = pybedtools.BedTool(target)
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
    
    if args.generate_bg == True and mode == "closest":
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
        bedfile = pybedtools.BedTool(input)
        orientations = orientation.split(",")
    
        #Create empy df to stoire results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation in orientations:
            for target in targets.split(","):
                targetfile = pybedtools.BedTool(target)
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

if __name__ == "__main__":
    main(args.mode,args.input,args.targets,args.background,args.orientation,args.genome,args.ov_fraction,args.randomization, args.outfile, args.exclude_intervals, args.exclude_ov, args.exclude_upstream, args.exclude_downstream)




