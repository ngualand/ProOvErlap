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
from scipy.stats import ks_2samp
from multiprocessing import Pool
    
#Suppress some warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", message="the index file is older than the FASTA file")
################################## Define Functions ############################################

#Create directory if not exist, if exist print it and exit
def create_directory(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Directory created: {folder_path}")
    else:
        print("Directory exists, use: " + folder_path)

#Check input parameters
def check_input_parameters(mode,w,alpha):
    if mode == "intersect" or mode == "closest":
        print("Mode set to " + mode)
    else:
        print("Mode is not set, choose one of intersect or closest")
        print("Exit")
        sys.exit()
    
    if w > 1 or w < 0:
        print("w must be between 0 and 1, exit")
        print("Exit")
        sys.exit() 
    if alpha > 1 or alpha < 0:
        print("alpha must be between 0 and 1, exit")
        print("Exit")
        sys.exit() 

#Read bed file, if have more than 6 columns return only the first 6 columns
def get_bed_file(input_path):
    bed_df = pd.read_csv(input_path, sep="\t", header=None, comment='#')
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

#Function to compute mean distance of feature in Input with respect to target bed file
def compute_random_closest_mean_dist(background_df, bed, target, s, S, ov_fraction, exclude_ov, exclude_upstream, exclude_downstream):
    random_features_index = random.sample(range(background_df.shape[0]), len(bed))
    background_df_random = background_df.iloc[random_features_index, :]
    distances = pybedtools.BedTool.from_dataframe(background_df_random, header = False).sort().closest(b=target.sort(), d=True, s=s, S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe(header=None).iloc[:, 12].to_list()
    return statistics.mean([abs(x) for x in distances])

#Function to compute mean distance of feature in Input with respect to target bed file using a randomly generated background
def compute_random_closest_mean_dist_random_bg(target, s, S, ov_fraction, exclude_ov, exclude_upstream, exclude_downstream, existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed):
    #generate random intervals 
    background_df = pd.DataFrame(generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed))
    distances = pybedtools.BedTool.from_dataframe(background_df, header = False).sort().closest(b=target.sort(), d=True, s=s, S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe(header=None).iloc[:, 12].to_list()
    return statistics.mean([abs(x) for x in distances])

#Function to test mean distance
def test_closeness(bed,target,targetname,background,ov_fraction, randomization, frame, table,strandness="strandless", exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False, thread = 1):
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    if exclude_intervals is not None:
        real_count = statistics.mean([abs(x) for x in bed.sort().intersect(b = pybedtools.BedTool.from_dataframe(get_bed_file(exclude_intervals)),s = s, S = S, v = True).closest(b = target.sort().intersect(b = pybedtools.BedTool.from_dataframe(get_bed_file(exclude_intervals)),s = s, S = S, v = True), d=True,s=s,S=S,f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe(header=None).iloc[:,12].to_list()])
        random_count = []
        background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals), s = s, S = S, v = True).to_dataframe(header=None)
    else:    
        real_count = statistics.mean([abs(x) for x in bed.sort().closest(b = target.sort(), d=True,s=s,S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe().iloc[:,12].to_list()])
        random_count = []
        background_df = background.to_dataframe(header = None)        
    with Pool(processes=thread) as pool:  
        args = [(background_df, bed, target, s, S, ov_fraction, exclude_ov, exclude_upstream, exclude_downstream)] * randomization
        random_count = pool.starmap(compute_random_closest_mean_dist, args)
    zscore = compute_z_score(real_count,random_count)
    pvalue = compute_pvalue(zscore)
    table = save_tables(real_count, random_count, strandness, targetname, table)
    frame = save_results(zscore,strandness,pvalue, targetname, real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table

#Function to test mean distance compared to a randomly generated background  
def test_closeness_random_bg(bed, target, targetname, ov_fraction, randomization, frame, table,strandness="strandless", exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False, thread = 1):
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]  
    if exclude_intervals is not None:
        real_count = statistics.mean([abs(x) for x in bed.intersect(b = pybedtools.BedTool(exclude_intervals), v = True, s=s, S=S).sort().closest(b = target.sort().intersect(b = pybedtools.BedTool(exclude_intervals),s = s, S = S, v = True), d=True,s=s, S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe(header=None).iloc[:,12].to_list()])
    else:
        real_count = statistics.mean([abs(x) for x in bed.sort().closest(b = target.sort(), d=True,s=s,S=S, f=ov_fraction, io=exclude_ov, iu=exclude_upstream, id=exclude_downstream, D="a").to_dataframe(header=None).iloc[:,12].to_list()])  
    input_bed_file = bed.fn
    existing_intervals, mean_length, chrom_count = read_bed_file(input_bed_file)
    if exclude_intervals is not None:
        exclude_intervals_bed = read_exclude_bed_file(exclude_intervals)
    else:
        exclude_intervals_bed = None
    num_random_intervals = len(existing_intervals)
    with Pool(processes=thread) as pool:
        args = [(target, s, S, ov_fraction, exclude_ov, exclude_upstream, exclude_downstream, existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed)] * randomization
        random_count = pool.starmap(compute_random_closest_mean_dist_random_bg, args)  
    zscore = compute_z_score(real_count,random_count)
    pvalue = compute_pvalue(zscore)
    table = save_tables(real_count, random_count, strandness, targetname, table)
    frame = save_results(zscore,strandness, pvalue, targetname, real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table  

#Compute Length againts a randomly generate background, just a test, since backgrounf is length and chromosomal frequency matched
def compute_length_random_bg(input_bed_file, exclude_intervals, genome):
    existing_intervals, mean_length, chrom_count = read_bed_file(input_bed_file)
    exclude_intervals_bed = read_exclude_bed_file(exclude_intervals) if exclude_intervals else None
    num_random_intervals = len(existing_intervals)
    random_intervals = pd.DataFrame(generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed))
    nucleotide_content_df = pybedtools.BedTool.from_dataframe(random_intervals, header = False).nucleotide_content(fi=genome, s=True).to_dataframe()
    return statistics.mean(nucleotide_content_df["15_seq_len"])

#Function to test mean length against a randomly generated background
def test_length_random_bg(bed,randomization, target, frame, table, genome, exclude_intervals = None, thread = 1):
    nucleotide_content = bed.nucleotide_content(fi = genome).to_dataframe()
    real_length = statistics.mean(nucleotide_content["15_seq_len"].to_list())
    with Pool(processes=thread) as pool:  
        args_list = [(bed.fn, exclude_intervals, genome)] * randomization
        random_count_length = pool.starmap(compute_length_random_bg, args_list)
    zscore = compute_z_score(real_length,random_count_length)
    pvalue = compute_pvalue(zscore)
    frame = save_results(zscore, "Length", pvalue, os.path.basename(target), real_length, statistics.mean(random_count_length), statistics.stdev(random_count_length), frame)
    table = save_tables(real_length, random_count_length, "Length", os.path.basename(target), table)
    return frame, table

#Function to test AT and GC content
def compute_gc_at_content_random_bg(input_bed_file, exclude_intervals, genome):   
    existing_intervals, mean_length, chrom_count = read_bed_file(input_bed_file)   
    exclude_intervals_bed = read_exclude_bed_file(exclude_intervals) if exclude_intervals else None
    num_random_intervals = len(existing_intervals)
    random_intervals = pd.DataFrame(generate_random_intervals_frequency( existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed))
    nucleotide_content_df = pybedtools.BedTool.from_dataframe(random_intervals, header = False).nucleotide_content(fi=genome, s=True).to_dataframe()
    return ( statistics.mean(nucleotide_content_df["8_pct_gc"]), statistics.mean(nucleotide_content_df["7_pct_at"]))

##Function to test AT and GC content agains a randomly generated background
def test_GC_AT_random_bg(bed,randomization, target, frame, table, genome, exclude_intervals = None, thread = 1):
    nucleotide_content = bed.nucleotide_content(fi = genome, s = True).to_dataframe()
    real_GC = statistics.mean(nucleotide_content["8_pct_gc"].to_list())
    real_AT = statistics.mean(nucleotide_content["7_pct_at"].to_list())
    random_count_GC = []
    random_count_AT = []
    with Pool(processes=thread) as pool:  
        args_list = [(bed.fn, exclude_intervals, genome)] * randomization       
        results = pool.starmap(compute_gc_at_content_random_bg, args_list)
    random_count_GC, random_count_AT = zip(*results)
    zscore_GC = compute_z_score(real_GC,random_count_GC)
    pvalue_GC = compute_pvalue(zscore_GC)
    zscore_AT = compute_z_score(real_AT,random_count_AT)
    pvalue_AT = compute_pvalue(zscore_AT)
    frame = save_results(zscore_GC, "GC", pvalue_GC, os.path.basename(target), real_GC, statistics.mean(random_count_GC), statistics.stdev(random_count_GC), frame)
    frame = save_results(zscore_AT, "AT", pvalue_AT, os.path.basename(target), real_AT, statistics.mean(random_count_AT), statistics.stdev(random_count_AT), frame)
    table = save_tables(real_GC, random_count_GC, "GC", os.path.basename(target), table)
    table = save_tables(real_AT, random_count_AT, "AT", os.path.basename(target), table)
    return frame, table

#Function to test number of intersections against a randly generated background
def test_enrichement_random_bg(bed, target, targetname,ov_fraction, randomization, frame, table, strandness = "strandless", exclude_intervals = None, thread = 1):
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    real = bed.intersect(b = target, u = True, s=s, S=S, f = ov_fraction)
    real_count = len(real)
    existing_intervals, mean_length, chrom_count = read_bed_file(bed.fn)
    if exclude_intervals is not None:
        exclude_intervals_bed = read_exclude_bed_file(exclude_intervals)
    else:
        exclude_intervals_bed = None
    num_random_intervals = len(existing_intervals)     
    with Pool(processes=thread) as pool:  
        args = [(target, s, S, ov_fraction,existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed)] * randomization
        random_count = pool.starmap(compute_random_intersect_random_bg, args) 
    zscore = compute_z_score(real_count,random_count)
    pvalue = compute_pvalue(zscore)
    table = save_tables(real_count, random_count, strandness, targetname, table)
    frame = save_results(zscore,strandness, pvalue, targetname, real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table

#Function to test number of intersections against a randly generated background
def compute_random_intersect_random_bg(target, s, S, ov_fraction, existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed):
    background_df_random = pd.DataFrame(generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed))
    return len(pybedtools.BedTool.from_dataframe(background_df_random, header = False).intersect(b=target, u=True, s=s, S=S, f=ov_fraction))

#Function to test the genomic localization of intersections against a randly generated background
def test_enrichement_regions_random_bg(bed, target, targetname, ov_fraction, randomization, frame, table, tmp, strandness = "strandless", exclude_intervals = None, regions_bed = "NA"):
    regions_bed = pd.read_table(tmp + "/" + "Reference_regions.bed")
    regions_bed.iloc[:, 1] = regions_bed.iloc[:, 1].clip(lower=0)
    regions_dict = {category: group for category, group in regions_bed.groupby(regions_bed.columns[3])}
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]    
    for region in regions_dict.keys():
        real = bed.intersect(b = target, u = True, s=s, S=S, f = ov_fraction).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region], header = False), u = True, s = s, S=S, f = ov_fraction)
        real_count = len(real)
        real_overlap = bed.intersect(b = target, u = True, s=s, S=S, f = ov_fraction)
        random_count = []
        for i in range(int(randomization)):
            input_bed_file = bed.intersect(b = target, u = True, s=s, S=S, f = ov_fraction).fn
            existing_intervals, mean_length, chrom_count = read_bed_file(input_bed_file)
            if exclude_intervals is not None:
                exclude_intervals_bed = read_exclude_bed_file(exclude_intervals)
            else:
                exclude_intervals_bed = None

            num_random_intervals = len(existing_intervals)
            random_intervals = generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals_bed)
            write_bed_file(tmp + "/" + "Random_0123456789.bed", random_intervals)
            background_df_random = tmp + "/" + "Random_0123456789.bed"
            random_count.append(len(pybedtools.BedTool(background_df_random).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region], header = False), u = True, s = s, S=S, f = ov_fraction)))
        zscore = compute_z_score(real_count,random_count)
        pvalue = compute_pvalue(zscore)
        table = save_tables(real_count, random_count, strandness, targetname + "||" + region, table)
        frame = save_results(zscore,strandness, pvalue, targetname + "||" + region, real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table

#Check if input path is a bed file
def is_bed_file(file_path):
    if not os.path.isfile(file_path):
        print(f"Error: The file '{file_path}' does not exist.")
        return False
    with open(file_path, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            columns = line.split('\t')
            if len(columns) < 6 :
                print(f"Error: Line {line_number} have less than 6 column, please provide a BED file with at least 6 column, name, score and strand can also be placeholder '.'")
                return False
            try:
                start = int(columns[1])
                end = int(columns[2])
            except ValueError:
                print(f"Error: Start and end positions are not integers at line {line_number}.")
                return False
            if start >= end:
                print(f"Error: Start position is greater than or equal to end position at line {line_number}.")
                return False
            name = columns[3]
            if name != '.' and not isinstance(name, str):
                print(f"Error: Invalid name format at line {line_number}. Name should be a non-empty string or '.' (dot).")
                return False
            score = columns[4]
            if score != '.' and not is_number(score):
                print(f"Error: Score is not a valid number or '.' at line {line_number}.")
                return False
            strand = columns[5]
            if strand not in ('+', '-', '.'):
                print(f"Error: Strand is not valid at line {line_number}. Expected '+', '-' or '.' (dot).")
                return False
    
    print(f"The file '{file_path}' appears to be a valid BED file.")
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

def compute_length(background_df, bed, genome):
    random_features_index = random.sample(range(background_df.shape[0]), len(bed))
    background_df_random = background_df.iloc[random_features_index, :]
    nucleotide_content_df = pybedtools.BedTool.from_dataframe(background_df_random, header = False).nucleotide_content(fi=genome, s=True).to_dataframe()
    mean_length = statistics.mean(nucleotide_content_df["15_seq_len"])
    return mean_length

#Test length
def test_length(bed,background,randomization, target, frame, table, genome, exclude_intervals = None, thread = 1):
    if exclude_intervals is not None:
        nucleotide_content = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v = True).nucleotide_content(fi = genome).to_dataframe()
        real_length = statistics.mean(nucleotide_content["15_seq_len"].to_list())
        random_count_length = []
        background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals), v = True).to_dataframe(header=None)    
    else:
        nucleotide_content = bed.nucleotide_content(fi = genome).to_dataframe()
        real_length = statistics.mean(nucleotide_content["15_seq_len"].to_list())
        random_count_length = []
        background_df = background.to_dataframe(header=None)   
    with Pool(processes=thread) as pool:  # Adjust number of processes
        # Prepare arguments as list of tuples for starmap
        args = [(background_df, bed, genome)] * randomization
        random_count_length = pool.starmap(compute_length, args)  # Run in parallel
    zscore = compute_z_score(real_length,random_count_length)
    pvalue = compute_pvalue(zscore)
    frame = save_results(zscore, "Length", pvalue, os.path.basename(target), real_length, statistics.mean(random_count_length), statistics.stdev(random_count_length), frame)
    table = save_tables(real_length, random_count_length, "Length", os.path.basename(target), table)
    return frame, table

#Test AT and GC content
def compute_gc_at(background_df, bed, genome):
    random_features_index = random.sample(range(background_df.shape[0]), len(bed))
    background_df_random = background_df.iloc[random_features_index, :]
    nucleotide_content_df = pybedtools.BedTool.from_dataframe(background_df_random).nucleotide_content(fi=genome, s=True).to_dataframe()
    mean_gc = statistics.mean(nucleotide_content_df["8_pct_gc"])
    mean_at = statistics.mean(nucleotide_content_df["7_pct_at"])

    return mean_gc, mean_at

#Test AT and GC content
def test_GC_AT(bed,background,randomization, target, frame, table, genome, exclude_intervals = None, thread = 1):
    if exclude_intervals is not None:
        nucleotide_content = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v = True).nucleotide_content(fi = genome, s = True).to_dataframe()
        real_GC = statistics.mean(nucleotide_content["8_pct_gc"].to_list())
        real_AT = statistics.mean(nucleotide_content["7_pct_at"].to_list())
        random_count_GC = []
        random_count_AT = []
        background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals), v = True).to_dataframe(header=None)     
    else:
        nucleotide_content = bed.nucleotide_content(fi = genome, s = True).to_dataframe()
        real_GC = statistics.mean(nucleotide_content["8_pct_gc"].to_list())
        real_AT = statistics.mean(nucleotide_content["7_pct_at"].to_list())
        random_count_GC = []
        random_count_AT = []
        background_df = background.to_dataframe(header=None)   
    with Pool(processes=thread) as pool:
        args = [(background_df, bed, genome)] * randomization
        results = pool.starmap(compute_gc_at, args)  # Run in parallel
    random_count_GC, random_count_AT = zip(*results)
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
def compute_random_intersect(background_df, bed, target, s, S, ov_fraction):
    random_features_index = random.sample(range(background_df.shape[0]), len(bed))
    background_df_random = background_df.iloc[random_features_index, :]
    return len(pybedtools.BedTool.from_dataframe(background_df_random, header = False).intersect(b=target, u=True, s=s, S=S, f=ov_fraction))

def test_enrichement(bed,target,targetname, background,ov_fraction, randomization, frame, table, strandness="strandless",exclude_intervals = None, thread = 1):
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    if exclude_intervals is not None:
        real = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b = target, u = True, s = s, S=S, f = ov_fraction)
        real_count = len(real)
        background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals),v =True, s=s, S=S).to_dataframe(header=None)   
    else:
        real = bed.intersect(b = target, u = True, s = s, S=S, f = ov_fraction)
        real_count = len(real)
        background_df = background.to_dataframe(header=None)       
    with Pool(processes=thread) as pool:  
        args = [(background_df, bed, target, s, S, ov_fraction)] * randomization
        random_count = pool.starmap(compute_random_intersect, args)    
    zscore = compute_z_score(real_count,random_count)
    pvalue = compute_pvalue(zscore)
    table = save_tables(real_count, random_count, strandness, targetname, table)
    frame = save_results(zscore,strandness,pvalue, targetname, real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
    return frame, table


def test_enrichement_regions(bed,target,targetname,background,ov_fraction, randomization, frame, table, tmp, strandness="strandless",exclude_intervals = None, regions_bed = "NA"):   
    regions_bed = pd.read_table(tmp + "/" + "Reference_regions.bed")
    regions_bed.iloc[:, 1] = regions_bed.iloc[:, 1].clip(lower=0)
    regions_dict = {category: group for category, group in regions_bed.groupby(regions_bed.columns[3])}
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    for region in regions_dict.keys():
        if exclude_intervals is not None:
            real = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b = target, u = True, s = s, S=S, f = ov_fraction).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region]), u = True, s = s, S=S, f = ov_fraction)
            real_count = len(real)
            real_overlap = len(bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b = target, u = True, s = s, S=S, f = ov_fraction))
            random_count = []
            background_df = background.intersect(b = pybedtools.BedTool(exclude_intervals),v =True, s=s, S=S).intersect(b = target, u = True, s = s, S=S, f = ov_fraction).to_dataframe()
        else:
            real = bed.intersect(b = target, u = True, s = s, S=S, f = ov_fraction).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region], header = False), u = True, s = s, S=S, f = ov_fraction)
            real_count = len(real)
            real_overlap = len(bed.intersect(b = target, u = True, s = s, S=S, f = ov_fraction))
            random_count = []
            background_df = background.intersect(b = target, u = True, s = s, S=S, f = ov_fraction).to_dataframe(header=None)           
        for i in range(int(randomization)):
            random_features_index = random.sample(range(background_df.shape[0]), real_overlap)
            background_df_random = background_df.iloc[random_features_index,:]
            random_count.append(len(pybedtools.BedTool.from_dataframe(background_df_random, header = False).intersect(b = pybedtools.BedTool.from_dataframe(regions_dict[region]), u = True, s = s, S=S, f = ov_fraction) ))
        zscore = compute_z_score(real_count,random_count)
        pvalue = compute_pvalue(zscore)
        table = save_tables(real_count, random_count, strandness, targetname + "||" + region , table)
        frame = save_results(zscore,strandness,pvalue, targetname + "||" + region , real_count, statistics.mean(random_count), statistics.stdev(random_count), frame)
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
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) > 3 else "."
            score = parts[4] 
            strand = parts[5]
            length = end - start
            total_length += length
            intervals.append((chrom, start, end, name, score, strand, length))
            chrom_count[chrom] += 1 
    mean_length = total_length / len(intervals) if intervals else 0
    return intervals, mean_length, chrom_count

#Function to read the exclude BED file and store intervals to exclude
def read_exclude_bed_file(file_path):
    exclude_intervals = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
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
        if (start < ex_end) and (end > ex_start): 
            return True
    return False

#Function to generate random intervals based on existing BED intervals
def generate_random_intervals_frequency(existing_intervals, mean_length, chrom_count, num_random_intervals, exclude_intervals=None):
    random_intervals = []
    total_intervals = sum(chrom_count.values())  
    chrom_ranges = {}
    for chrom, start, end, name, score, strand, length in existing_intervals:
        if chrom not in chrom_ranges:
            chrom_ranges[chrom] = {'min': start, 'max': end}
        else:
            chrom_ranges[chrom]['min'] = min(chrom_ranges[chrom]['min'], start)
            chrom_ranges[chrom]['max'] = max(chrom_ranges[chrom]['max'], end)
    for chrom, count in chrom_count.items():
        chrom_intervals = round((count / total_intervals) * num_random_intervals)
        chrom_range = chrom_ranges[chrom]
        chrom_min = chrom_range['min']
        chrom_max = chrom_range['max']       
        for _ in range(chrom_intervals):
            while True:  
                start = random.randint(chrom_min, chrom_max)
                length = abs(random.gauss(mean_length, mean_length * 0.1)) 
                length = max(1, int(length))  
                end = min(start + length, chrom_max)
                if end <= start:
                    continue    
                if exclude_intervals and chrom in exclude_intervals:    
                    if any(start < ex_end and end > ex_start for ex_start, ex_end in exclude_intervals[chrom]):
                        continue  
                name = f"random_{random.randint(1, 10000000)}"
                score = random.randint(0, 1000)
                strand = random.choice(["+", "-"])
                random_intervals.append((chrom, start, end, name, score, strand))
                break    
    return random_intervals

#Function to write intervals to a 6 column BED file
def write_bed_file(file_path, intervals):
    with open(file_path, 'w') as f:
        for chrom, start, end, name, score, strand in intervals:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")


# Function to compute enrichment score (ES)
def compute_enrichment_score(ranked_list, gene_set):
    N = len(ranked_list)  
    Nh = len(gene_set)  
    Nm = N - Nh 
    hit_indices = [i for i, gene in enumerate(ranked_list) if gene in gene_set]
    if not hit_indices:
        return 0, [], []
    p_hit = np.cumsum([1/Nh if i in hit_indices else 0 for i in range(N)])
    p_miss = np.cumsum([1/Nm if i not in hit_indices else 0 for i in range(N)])
    ES_values = p_hit - p_miss
    ES = max(ES_values, key=abs) 
    peak_idx = np.argmax(np.abs(ES_values))  
    if ES >= 0:
        leading_regions = ranked_list[:peak_idx + 1]
    else:
        leading_regions = ranked_list[peak_idx + 1:]
    leading_regions_str = ",".join(leading_regions)
    return ES, ES_values, hit_indices, leading_regions_str


# Function to compute enrichment score for shuffled genes
def compute_null_ES(ranked_list, gene_set):
    shuffled_genes = np.random.permutation(ranked_list)
    shuffled_ES, _, _, _ = compute_enrichment_score(shuffled_genes, gene_set)
    return shuffled_ES

# Function to run permutation test
def permutation_test(ranked_list, gene_set, num_permutations=100, thread = 1):
    real_ES, _, _ , leading_regions = compute_enrichment_score(ranked_list, gene_set)
    args = [(ranked_list, gene_set)] * num_permutations 
    with Pool(processes=thread) as pool:  
        null_ES = pool.starmap(compute_null_ES, args)
    null_ES = np.array(null_ES)
    p_value = np.sum(np.abs(null_ES) >= np.abs(real_ES)) / num_permutations
    null_mean = np.mean(null_ES)
    null_std = np.std(null_ES)
    Z_score = (real_ES - null_mean) / null_std if null_std != 0 else 0
    return real_ES, p_value, Z_score, null_ES, leading_regions

def scale_value(value, min_range, max_range, new_min, new_max):
    scaled_value = ((value - min_range) / (max_range - min_range)) * (new_max - new_min) + new_min
    return scaled_value
                
def weighted_score(s, d, alpha=0.5, w=0.5, ascending=True):
    N = len(s)
    if ascending:
        r = (np.argsort(s) + 1) / N 
        W = d * alpha + (1 -r) * (1 - alpha)
        S_new = s * ((1 - w) + w * (1 - W))
    else:
        r = (np.argsort(-s) + 1) / N
        W = d * alpha + (1 - r) * (1 - alpha)
        S_new = s * ((1 - w) + w * W)   
    return S_new



def read_and_rank_bed(file_path, ascending=False):
    bed_df = pd.read_csv(file_path, sep="\t", header=None)   
    bed_df = bed_df.iloc[:, :6]
    bed_df.columns = ["chrom", "start", "end", "name", "score", "strand"]
    bed_df["score"] = pd.to_numeric(bed_df["score"], errors="coerce")
    if bed_df["score"].isna().any():
        print("The 'score' column contains non-numeric or NA or empty values, this analysis require all regions to have a numeric score")
        sys.exit() 
    bed_df["name"] = bed_df["name"].replace({".": None})  
    ranked_bed = bed_df.sort_values(by="score", ascending=ascending)
    bed_df["name"] = bed_df["name"].fillna(pd.Series([f"Region_{i+1}" for i in range(len(bed_df))], index=bed_df.index))
    return ranked_bed

#Test overlap
def test_enrichement_rank(bed,target,ov_fraction, randomization, strandness="strandless",exclude_intervals = None, ascending = False, WeightRanking = False, alpha = 0.5, w = 0.25, thread = 1):
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    NAME = os.path.basename(bed)
    TARGET = os.path.basename(target)
    bed = pybedtools.BedTool.from_dataframe(read_and_rank_bed(bed, ascending=ascending ), header = False)
    target = pybedtools.BedTool.from_dataframe(get_bed_file(target), header = False)      
    if exclude_intervals is not None:
        overlap_df = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b=target, s=s, S=S, f=ov_fraction, wo=True).to_dataframe(header=None)
        overlap_df = overlap_df.loc[overlap_df.groupby(overlap_df.columns[2])[overlap_df.columns[11]].idxmax()]
        intersected_regions = overlap_df.iloc[:,3].tolist()
        non_overlap_df = bed.intersect(b = pybedtools.BedTool(exclude_intervals), v =True, s=s, S=S).intersect(b=target, s=s, S=S, f=ov_fraction, v=True).to_dataframe(header=None)
        non_overlap_df.columns = range(0,non_overlap_df.shape[1])
    else:
        overlap_df = bed.intersect(b=target, s=s, S=S, f=ov_fraction, wo=True).to_dataframe(header=None)
        overlap_df = overlap_df.loc[overlap_df.groupby(overlap_df.columns[2])[overlap_df.columns[11]].idxmax()]
        intersected_regions = overlap_df.iloc[:,3].tolist()
        non_overlap_df = bed.to_dataframe(header = None)
        non_overlap_df = non_overlap_df.loc[~non_overlap_df.iloc[:,3].isin(intersected_regions),:]
        #non_overlap_df = bed.intersect(b=target, s=s, S=S, f=ov_fraction, v=True).to_dataframe(header=None)
        non_overlap_df.columns = range(0,non_overlap_df.shape[1])
    combined_df = pd.concat([overlap_df, non_overlap_df], ignore_index=True)
    combined_df["overlap_length"] = combined_df.iloc[:, -1].fillna(0)
    combined_df["length"] = combined_df.iloc[:,2] - combined_df.iloc[:,1]
    combined_df["overlap_fraction"] = combined_df["overlap_length"] / combined_df["length"]
    combined_df["Zscore_overlap_fraction"] = (combined_df["overlap_fraction"] - np.mean(combined_df["overlap_fraction"])) / np.std(combined_df["overlap_length"])
    combined_df["Zscore_score"] = (combined_df.iloc[:, 4] - np.mean(combined_df.iloc[:, 4])) / np.std(combined_df.iloc[:, 4])
    min_range = combined_df["Zscore_score"].min()
    max_range = combined_df["Zscore_score"].max()
    new_min = 0
    new_max = 1        
    combined_df["scaled_Zscore_score"] = combined_df["Zscore_score"].apply(lambda x: scale_value(x, min_range, max_range, new_min, new_max))
    combined_df["scaled_overlap_length"] = combined_df["overlap_length"].apply(lambda x: scale_value(x, min_range, max_range, new_min, new_max))
    if WeightRanking == True:
        combined_df["weighted_score"] = weighted_score(combined_df["scaled_Zscore_score"],combined_df["overlap_fraction"], alpha, w, ascending)
        regions_ranked = combined_df.sort_values(by="weighted_score", ascending=ascending)  
    if WeightRanking == False:
        regions_ranked = combined_df.sort_values(by="Zscore_score", ascending=ascending) 
    regions_ranked = regions_ranked.iloc[:, 3].tolist()
    ES, p_value, Zscore, null_ES, leading_regions = permutation_test(regions_ranked, intersected_regions, num_permutations=randomization, thread = thread)
    ES, ES_values, hit_indices, _ = compute_enrichment_score(regions_ranked, intersected_regions)
    res = pd.DataFrame({"Input" : NAME, "Target" : TARGET, "ES" : ES, "Pval" : p_value, "Type": "Intersect", "NRand" : str(randomization), "Orientation" : strandness, "Leading" : leading_regions}, index=[0])   
    #forplot = pd.DataFrame({"Pval" : p_value, "ES_values" : "|".join(map(str, ES_values)), "hit_indices" : "|".join(map(str, hit_indices)), "Null_ES": "|".join(map(str, null_ES)), "ES": str(ES)}, index=[0])  
    max_len = max(len(ES_values), len(hit_indices), len(null_ES))
    forplot = pd.DataFrame({"Pval": [p_value] * max_len,"ES": [ES] * max_len,"ES_values": list(ES_values) + [None] * (max_len - len(ES_values)) , "hit_indices": hit_indices + [None] * (max_len - len(hit_indices)), "Null_ES": list(null_ES) + [None] * (max_len - len(null_ES)) })
    return res, forplot

def compute_shuffle(i, y, cumulative_real):
    y_shuffled = np.random.permutation(y)  
    cumulative_shuffle = np.cumsum(y_shuffled) / np.sum(y_shuffled) 
    ks_stat, p_value = ks_2samp(cumulative_real, cumulative_shuffle)    
    max_index = np.argmax(abs(cumulative_real - cumulative_shuffle))
    enrichment_score = (cumulative_real - cumulative_shuffle)[max_index]
    return p_value, enrichment_score, cumulative_shuffle

def parallel_shuffling(n_shuffles, y, cumulative_real, thread):
    with Pool(processes=thread) as pool:
        inputs = [(i, y, cumulative_real) for i in range(n_shuffles)]
        results = pool.starmap(compute_shuffle, inputs)
    p_values, enrichment_scores, shuffle_cdfs = zip(*results)   
    return list(p_values), list(enrichment_scores), list(shuffle_cdfs)

#Test closest
def KS_test_closest(x, y, n_shuffles=100, thread = 1):
    cumulative_real = np.cumsum(y) / np.sum(y)  
    p_values, enrichment_scores, shuffle_cdfs = parallel_shuffling(n_shuffles, y, cumulative_real, thread)   
    mean_shuffle_cdf = np.mean(shuffle_cdfs, axis=0)
    diff_list_abs = list(np.abs(cumulative_real - mean_shuffle_cdf))
    diff_list = list(cumulative_real - mean_shuffle_cdf) 
    max_value = max(diff_list_abs)
    max_index = diff_list_abs.index(max_value)   
    enrichment_score_max = diff_list[max_index]
    enrichment_scores = cumulative_real - mean_shuffle_cdf 
    ks_stat, p_value = ks_2samp(cumulative_real, mean_shuffle_cdf)
    p_value_mean = np.mean(p_values)
    if enrichment_score_max >= 0:
        leading = x[max_index:]
    else:
        leading = x[:max_index]   
    leading_str = ",".join(leading)
    return enrichment_score_max, p_value, leading_str, cumulative_real, mean_shuffle_cdf, max_index

def test_closest_rank(bed,target,ov_fraction, randomization, strandness="strandless", exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False, ascending = False, WeightRanking = False, alpha = 0.5, w = 0.25, thread = 1):
    s = get_params_strand(strandness)["s"]
    S = get_params_strand(strandness)["S"]
    NAME = os.path.basename(bed)
    TARGET = os.path.basename(target)
    bed = pybedtools.BedTool.from_dataframe(read_and_rank_bed(bed, ascending=ascending), header = False)
    target = pybedtools.BedTool.from_dataframe(get_bed_file(target), header = False)
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
    closest_df["scaled_specular_Distance"] = 1 - closest_df["scaled_Distance"]
    if WeightRanking == False:
        regions_ranked = closest_df.sort_values(by="Score", ascending=ascending)
    if WeightRanking == True:
        closest_df["weighted_score"] = weighted_score(closest_df["scaled_Zscore_score"],closest_df["scaled_specular_Distance"], alpha, w, ascending)
        regions_ranked = closest_df.sort_values(by="weighted_score", ascending=ascending)   
    regions_ranked["Rank"] = [x for x in range(1,regions_ranked.shape[0] + 1 )] 
    ES, Pval, leading, cumulative_real, mean_shuffle_cdf, peak_idx = KS_test_closest(regions_ranked.iloc[:,3], regions_ranked["scaled_specular_Distance"], randomization, thread)    
    res = pd.DataFrame({"Input": NAME, "Target": TARGET, "ES": ES, "Pval" : Pval, "Type" : "Closest", "NRand": str(randomization), "Orientation": strandness, "Leading": leading}, index = [0])
    max_len = max(len(cumulative_real), len(mean_shuffle_cdf))
    forplot = pd.DataFrame({"Peak_idx" : [peak_idx] * max_len, "ES": [ES] * max_len, "Pval" : [Pval] * max_len, "cumulative_real": list(cumulative_real) + [None] * (max_len - len(cumulative_real))  , "mean_shuffle_cdf": list(mean_shuffle_cdf) + [None] * (max_len - len(mean_shuffle_cdf) )})
    
    return res, forplot


###### MAIN #########
def main(mode,input,targets,background,orientation,genome,ov_fraction,randomization,outfile,outdir, exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False, RankTest = False, WeightRanking = False, alpha = 0.5, w = 0.25, thread = 1, tmp = ".", GenomicLocalization = False, gtf = False, bed = None, generate_bg = False, test_AT_GC = False, test_lengths = False, Ascending_RankOrder = False):
    #Check input parameters
    check_input_parameters(mode,w,alpha)
    #Check and create tmp directory
    create_directory(tmp)
    #Create Table directory
    create_directory(outdir)
    #Set tmp dir for pybedtools
    pybedtools.set_tempdir(tmp)   
    if GenomicLocalization:
        if gtf is not None:
            script_r = os.path.join(os.path.dirname(__file__), "Create_bed_genomicRegions.R")
            subprocess.run(["Rscript", script_r, gtf, tmp, genome])
            regions_bed = tmp + "/" + "Reference_regions.bed"
        if bed is not None:
            regions_bed = bed
     
    if generate_bg == False and mode == "intersect" and RankTest == False:
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
        bedfile = pybedtools.BedTool.from_dataframe(get_bed_file(input), header = False)
        backgroundfile = pybedtools.BedTool.from_dataframe(get_bed_file(background), header = False)
        orientations = orientation.split(",")
    
        #Create empy df to store results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation  in orientations:
            for target in targets.split(","):
                targetname = os.path.basename(target)
                targetfile = pybedtools.BedTool.from_dataframe(get_bed_file(target), header = False)
                res = test_enrichement(bedfile,targetfile,targetname,backgroundfile,ov_fraction,randomization, frame, table, strandness = orientation, exclude_intervals = exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]
                if GenomicLocalization:
                    res = test_enrichement_regions(bedfile,targetfile,targetname,backgroundfile,ov_fraction,randomization, frame, table,tmp, strandness = orientation, exclude_intervals = exclude_intervals, regions_bed = regions_bed)
                    frame = res[0]
                    table = res[1]

        if test_AT_GC:
            for target in targets.split(","):
                res = test_GC_AT(bedfile,backgroundfile,randomization,target, frame, table, genome, exclude_intervals = exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]
            
        if test_lengths:
            for target in targets.split(","):
                res = test_length(bedfile,backgroundfile,randomization, target, frame, table,genome, exclude_intervals =  exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]
            
        
            
        frame.to_csv(outfile, index=False, sep = "\t")
        table.to_csv(outdir + "/Table_Intersect.txt", index=False, sep = "\t")
    
    if generate_bg == False and mode == "closest" and RankTest == False:
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
        bedfile = pybedtools.BedTool.from_dataframe(get_bed_file(input), header = False)
        backgroundfile = pybedtools.BedTool.from_dataframe(get_bed_file(background), header = False)
        orientations = orientation.split(",")
    
        #Create empy df to stoire results
        frame = pd.DataFrame()
        table = pd.DataFrame()
        
        for orientation in orientations:
            for target in targets.split(","):
                targetname = os.path.basename(target)
                targetfile = pybedtools.BedTool.from_dataframe(get_bed_file(target), header = False)
                res = test_closeness(bedfile, targetfile, targetname, backgroundfile, ov_fraction, randomization, frame, table,strandness=orientation, exclude_intervals = exclude_intervals, exclude_ov = exclude_ov, exclude_upstream = exclude_upstream, exclude_downstream=exclude_downstream, thread = thread)
                frame = res[0]
                table = res[1]
            
    
        if test_AT_GC:
            for target in targets.split(","):
                res = test_GC_AT(bedfile,backgroundfile,randomization,target, frame, table, genome, exclude_intervals = exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]
            
        if test_lengths:
            for target in targets.split(","):
                res = test_length(bedfile,backgroundfile,randomization, target, frame, table, genome, exclude_intervals = exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]
            
        
            
        frame.to_csv(outfile, index=False, sep = "\t")
        table.to_csv(outdir + "/Table_Closest.txt", index=False, sep = "\t")

    
    if generate_bg == True and mode == "intersect" and RankTest == False:
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
        bedfile = pybedtools.BedTool.from_dataframe(get_bed_file(input), header = False)
        orientations = orientation.split(",")
    
        #Create empy df to stoire results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation in orientations:
            for target in targets.split(","):
                targetname = os.path.basename(target)
                targetfile = pybedtools.BedTool.from_dataframe(get_bed_file(target), header = False)
                res = test_enrichement_random_bg(bedfile,targetfile,targetname,ov_fraction,randomization, frame, table, strandness= orientation, exclude_intervals = exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]
                if GenomicLocalization:
                    res = test_enrichement_regions_random_bg(bedfile,targetfile,targetname,ov_fraction,randomization, frame, table, tmp, strandness= orientation, exclude_intervals = exclude_intervals, regions_bed = regions_bed)
                    frame = res[0]
                    table = res[1]
    
        if test_AT_GC:
            for target in targets.split(","):
                res = test_GC_AT_random_bg(bedfile,randomization,target, frame, table, genome, exclude_intervals = exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]
            
        if test_lengths:
            for target in targets.split(","):
                res = test_length_random_bg(bedfile,randomization, target, frame, table, genome, exclude_intervals = exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]
            
        
            
        frame.to_csv(outfile, index=False, sep = "\t")
        table.to_csv(outdir + "/Table_Intersect.txt", index=False, sep = "\t")
    
    if generate_bg == True and mode == "closest" and RankTest == False:
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
        bedfile = pybedtools.BedTool.from_dataframe(get_bed_file(input), header = False)
        orientations = orientation.split(",")
    
        #Create empy df to stoire results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation in orientations:
            for target in targets.split(","):
                targetname = os.path.basename(target)
                targetfile = pybedtools.BedTool.from_dataframe(get_bed_file(target), header = False)
                res = test_closeness_random_bg(bedfile, targetfile, targetname, ov_fraction, randomization, frame, table,strandness=orientation, exclude_intervals = exclude_intervals, exclude_ov = exclude_ov, exclude_upstream = exclude_upstream, exclude_downstream=exclude_downstream, thread = thread)
                frame = res[0]
                table = res[1]
    
        if test_AT_GC:
            for target in targets.split(","):
                res = test_GC_AT_random_bg(bedfile,randomization,target, frame, table,genome, exclude_intervals = exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]
  
        if test_lengths:
            for target in targets.split(","):
                res = test_length_random_bg(bedfile,randomization, target, frame, table,genome,  exclude_intervals = exclude_intervals, thread = thread)
                frame = res[0]
                table = res[1]

        frame.to_csv(outfile, index=False, sep = "\t")
        table.to_csv(outdir + "/Table_Closest.txt", index=False, sep = "\t")
            
    if RankTest == True and mode == "intersect":
        print("Perform Rank Mode Test using Intersect")
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
        orientations = orientation.split(",")
    
        #Create empy df to store results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation  in orientations:
            for target in targets.split(","):
                res, forplot = test_enrichement_rank(input,target,ov_fraction, randomization, strandness=orientation,exclude_intervals = exclude_intervals, ascending = Ascending_RankOrder, WeightRanking = WeightRanking, alpha = alpha, w = w, thread = thread)
                frame = pd.concat([frame, res])
        
        frame.to_csv(outfile, index=False, sep = "\t")
        forplot.to_csv(outdir + "/Table_Rank_Intersect.txt", index=False, sep = "\t")

    if RankTest == True and mode == "closest":
        print("Perform Rank Mode Test using Closest")
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
        orientations = orientation.split(",")
    
        #Create empy df to store results
        frame = pd.DataFrame()
        table = pd.DataFrame()

        for orientation  in orientations:
            for target in targets.split(","):
                res, forplot = test_closest_rank(input,target,ov_fraction, randomization,strandness="strandless", exclude_intervals = None, exclude_ov = False, exclude_upstream = False, exclude_downstream=False, ascending = Ascending_RankOrder, WeightRanking = WeightRanking, alpha = alpha,w = w, thread = thread)
                frame = pd.concat([frame, res])
        
        frame.to_csv(outfile, index=False, sep = "\t")
        forplot.to_csv(outdir + "/Table_Rank_Closest.txt", index=False, sep = "\t")
        

def cli():
    #Create argument parser
    parser = argparse.ArgumentParser(description="ProOvErlap")
    parser.add_argument("--mode", required = True, help = "Define mode: intersect or closest: intersect count the number of overlapping elements while closest test the distance. In closest mode if a feature overlap a target the distance is 0, use --exclude_ov to test only for non-overlapping regions")
    parser.add_argument("--input", required = True, help="Input bed file, must contain 6 or more columns, name and score can be placeholder but score is required in --RankTest mode, strand is used only if some strandess test are requested")
    parser.add_argument("--targets", required = True, help="Target bed file(s) (must contain 6 or more columns) to test enrichement against, if multiple files are supplied N independent test against each file are conducted, file names must be comma separated, the name of the file will be use as the name output")
    parser.add_argument("--background", help="Background bed file (must contain 6 or more columns), should be a superset from wich input bed file is derived")
    parser.add_argument("--randomization", type=int, help="Number of randomization, default: 100", default = 100)
    parser.add_argument("--genome", help="Genome fasta file used to retrieve sequence features like AT or GC content and length, needed only for length or AT/GC content tests")
    parser.add_argument("--tmp", default=".", help="Temporary directory for storing intermediate files. Default is current working directory")
    parser.add_argument("--outfile", required=True, help="Full path to the output file to store final results in tab format")
    parser.add_argument("--outdir", required=True, help="Full path to output directory to store tables for plot, it is suggested to use a different directory for each analysis. It will be created")
    parser.add_argument("--orientation", default = "strandless", help="Name of test(s) to be performed: concordant, discordant, strandless, or a combination of them. If multiple tests are required tests names must be comma separated, no space allowed")
    parser.add_argument("--ov_fraction", default="0.0000000001", help="Minimum overlap required as a fraction from input BED file to consider 2 features as overlapping. Default is 1E-9 (i.e. 1bp)")
    parser.add_argument("--generate_bg", action = "store_true", help="This option activatates the generation of random bed intervals to test enrichment against, use this instead of background. Use only if background file cannot be used or is not available")
    parser.add_argument("--exclude_intervals", default = None, help="Exclude regions overlapping with regions in the supplied BED file")
    parser.add_argument("--exclude_ov", action = "store_true", help="Exclude overlapping regions between Input and Target file in closest mode")
    parser.add_argument("--exclude_upstream", action = "store_true", help="Exclude upstream region in closest mode, only for stranded files, not compatible with exclude_downstream")
    parser.add_argument("--exclude_downstream", action = "store_true", help="Exclude downstream region in closest mode, only for stranded files, not compatible with exclude_upstream")
    parser.add_argument("--test_AT_GC", action = "store_true", help="Test AT and GC content")
    parser.add_argument("--test_lengths", action = "store_true", help="Test feature length")
    parser.add_argument("--GenomicLocalization", action = "store_true", help= "Test also the genomic localization and enrichment of founded overlaps, i.e TSS,Promoter,exons,introns,UTRs - Available only in intersect mode. Must provide a GTF file to extract genomic regions (--gtf), alternatively directly provide a bed file (--bed) with custom annotations")
    parser.add_argument("--gtf", help="GTF file, only to test genomic localization of founded overlap, gtf file will be used to create genomic regions: promoter, tss, exons, intron, 3UTR and 5UTR")
    parser.add_argument("--bed", help="BED file, only to test genomic localization of founded overlap, bed file will be used to test enrichment in different genomic regions, annotation must be stored as 4th column in bed file, i.e name field")
    parser.add_argument("--RankTest", action = "store_true", help="Activates the Ranking analyis, require BED to contain numerical value in 4th column")
    parser.add_argument("--Ascending_RankOrder", action = "store_true", help="Activate the Sort Ascending in RankTest analysis")
    parser.add_argument("--WeightRanking", action = "store_true", help="Weight the ranking test, this is done by increase or decrease the score value in the BED file based on their relative rank and/or distance and/or fractional overlap")
    parser.add_argument("--alpha", default = "0.5", type=float, help="Relative Influence of the overlap fraction/distance (with respect to ranking) in weightRanked test, only if --WeightRanking is active, must be between 0 and 1")
    parser.add_argument("--w", default = "0.25", type=float, help="Strength of the Weight for the ranking test, only if --WeightRanking is active, must be between 0 and 1")
    parser.add_argument("--thread", default = "1",type = int, help="Number of Threads for parallel computation")

    #Parse argument
    args = parser.parse_args()

    #Run Main
    main(args.mode,args.input,args.targets,args.background,args.orientation,args.genome,args.ov_fraction,args.randomization, args.outfile, args.outdir, args.exclude_intervals, args.exclude_ov, args.exclude_upstream, args.exclude_downstream, args.RankTest, args.WeightRanking, args.alpha, args.w, args.thread, args.tmp, args.GenomicLocalization, args.gtf, args.bed, args.generate_bg, args.test_AT_GC, args.test_lengths, args.Ascending_RankOrder)



if __name__ == "__main__":
    cli()