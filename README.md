# ProOvErlap - Assessing feature proximity/overlap and testing statistical significance from genomic intervals
# Overview
Genomic feature overlap plays a crucial role in bioinformatics, occurring when two genomic intervals, often represented as BED files, are positioned within the same genomic regions. In contrast, feature proximity refers to the spatial closeness of genomic elements. For instance, gene promoters frequently overlap with or are located near the genes they regulate. Both overlap and proximity are particularly relevant in epigenetic studies, where regions enriched for specific epigenetic modifications or accessible chromatin can provide insights into complex molecular phenotypes. To facilitate the analysis of these genomic relationships, we introduce a computational tool designed to process BED-format data. This method quantitatively evaluates the extent of overlap or proximity between genomic features while assessing their statistical significance using a non-parametric randomization test. The goal is to determine whether the observed patterns deviate from what would be expected by chance. The tool is user-friendly, requiring only a single command-line execution for efficient analysis. Additionally, it generates clear visualizations and high-quality figures suitable for publication. Overall, this approach enhances the systematic assessment of feature overlap and proximity, offering a valuable resource for identifying meaningful genomic interactions in both normal and disease contexts.

![ProOvErlap Logo](Fig5.jpg)

# How to install:
ProOvErlap does not require installation; simply run it as a Python script using:  
python3 prooverlap.py --help  
Please note that certain Python and R libraries must be installed for the software to function properly. Additionally, ProOvErlap relies on an external R script for specific steps, so always ensure that you execute the code from within the main ProOvErlap directory.

# Needed Libraries
python Libraries:

- Biopython
- pandas
- statistics
- scipy
- sys
- argparse
- os
- tempfile
- time
- pybedtools
- random
- warnings
- collections
- subprocess

R Libraries:

- tidyverse
- argparse
- ggplot2
- AnnotationHub
- GenomicRanges
- rtracklayer
- GenomicFeatures
- Biostrings

# Input and Outputs:

ProOvErlap accepts three input files: two required BED files (input and target) and one optional BED file (background, optional but recommended). The software outputs a main table containing the results of the analysis. Additionally, it generates a second table that can be used as input for generating a density plot, which shows how far the real values deviate from what would be expected by chance. The density plot should be performed using the Density_plot.R script.

# Usage:

```
usage: python3 prooverlap.py --help

usage: prooverlap.py [-h] --mode MODE --input INPUT --targets TARGETS [--background BACKGROUND] [--randomization RANDOMIZATION] [--genome GENOME] [--tmp TMP]
                     [--outfile OUTFILE] --orientation ORIENTATION [--ov_fraction OV_FRACTION] [--generate_bg] [--exclude_intervals EXCLUDE_INTERVALS] [--exclude_ov]
                     [--exclude_upstream] [--exclude_downstream] [--test_AT_GC] [--test_length] [--GenomicLocalization] [--gtf GTF] [--bed BED]

options:
  -h, --help            show this help message and exit
  --mode MODE           Define mode: intersect or closest, intersect count the number of overlapping elements while closest test the distance. In closest if a feature
                        overlap a target the distance is 0
  --input INPUT         Input bed file, must contain 6 column, name and score can be placeholder and are not used, strand is used only if some strandess test are requested
  --targets TARGETS     Target bed files (must contain 6 columns), to test enrichement against, if multiple files are supplied N independent test against each file are
                        conducted, file names must be comma separated, the name of the file will be use as the name output
  --background BACKGROUND
                        Background bed file (must contain 6 columns) to test enrichement aginst, should be a superset from wich input bed file is derived
  --randomization RANDOMIZATION
                        Number of randomization, default 100
  --genome GENOME       genome fasta file used to retrieve sequence from bed files, needed only for length or AT/GC content
  --tmp TMP             Default is current working dir, Location of the directory to store temporary files, after running the sofware automatically clean up tmp files, if
                        the software do not exit properly it may not clean up tmp file!
  --outfile OUTFILE     Full path to output file to store results, it will be created
  --orientation ORIENTATION
                        Name of test/tests to be performed: concordant, discordant, strandless, or a combination of them, comma separated, no space allowed
  --ov_fraction OV_FRACTION
                        Minimum overlap required as a fraction from input BED file to consider 2 features as overlapping. Default is 1E-9 (i.e. 1bp)
  --generate_bg         Generates random bed intervals to test enrichment against, use this instead of background. Use only if background file cannot be used
  --exclude_intervals EXCLUDE_INTERVALS
                        Exclude those regions in both random background generation and feature testing
  --exclude_ov          Do not count overlapping region in closest mode
  --exclude_upstream    Do not count upstream region in closest mode, only for stranded files, not compatible with exclude_downstream
  --exclude_downstream  Do not count downstream region in closest mode, only for stranded files, not compatible with exclude_upstream
  --test_AT_GC          Test AT and GC content
  --test_length         Test feature length
  --GenomicLocalization
                        Test the genomic localization and enrichment of founded overlaps, i.e TSS,Promoter,exons,introns,UTRs - Available only in intersect mode. Must
                        provide a GTF file to extract genomic regions (--gtf), alternatively directly provide a bed file (--bed) with custom annotations
  --gtf GTF             GTF file, only to test genomic localization of founded overlap, gtf file will be used to create genomic regions: promoter, tss, exons, intron, 3UTR
                        and 5UTR
```

# How to plot results?
ProOvErlap supports the creation of two main types of graphical outputs (although you may also perform your own plots, as all data are saved to files). The first one is a density plot (generated by the Density_plot.R script), which shows how far the obtained results deviate from what would be expected by chance. Moreover, ProOvErlap also creates heatmaps of the Z-score for each target and, optionally, genomic regions or custom regions, using the Heatmap.R script.

```
Density_plot.R: Required arguments: 
input_table: the main output of Enrichard.py "ex: Results.txt",
randomizations: auto generated output of EnricharD.py containing the randomization table "ex: Tables.txt",
test: mode used in EnricharD.py, it must be intersect or closest (default: intersect)
outfile: name of the suffix of output file (default: Density_plot)
format: format used to save the output file, could be png, pdf or svg (default: png)

Heatmap.R: Required arguments:
input_table: main output of EnricharD.py when the option "GenomicLocalization" is set
outfile: name of output file (default = "Heatmap")
format: format used to save the output file, could be png, pdf or svg (default: png)
title: title of the plot (default: "")
```

# Development 
ProOvErlap was developed by Nicolò Gualandi (former post-doc in the Laboratory of Prof. Claudio Brancolini @ UniUd) and Alessio Bertozzo (PhD student in the Laboratory of Prof. Claudio Brancolini @ UniUd), under the supervision of Prof. Claudio Brancolini (Professor of Cell Biology, Department of Medicine, Università degli Studi di Udine, https://people.uniud.it/page/claudio.brancolini)  

ProOvErlap is actively being improved. If you would like to contribute, we welcome your comments and feedback.  
