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
- numpy
- scipy.stats
- multiprocessing

R Libraries:

- tidyverse
- argparse
- ggplot2
- AnnotationHub
- GenomicRanges
- rtracklayer
- GenomicFeatures
- Biostrings
- Argparse

# Input and Outputs:

ProOvErlap accepts three input files: two required BED files (input and target) and one optional BED file (background, optional but recommended). The software outputs a main table containing the results of the analysis. Additionally, it generates a second table that can be used as input for generating a density plot, which shows how far the real values deviate from what would be expected by chance. The density plot should be performed using the Density_plot.R script.

# Usage:

```
usage: prooverlap.py [-h] --mode MODE --input INPUT --targets TARGETS [--background BACKGROUND] [--randomization RANDOMIZATION] [--genome GENOME]
                     [--tmp TMP] --outfile OUTFILE --outdir OUTDIR [--orientation ORIENTATION] [--ov_fraction OV_FRACTION] [--generate_bg]
                     [--exclude_intervals EXCLUDE_INTERVALS] [--exclude_ov] [--exclude_upstream] [--exclude_downstream] [--test_AT_GC] [--test_length]
                     [--GenomicLocalization] [--gtf GTF] [--bed BED] [--RankTest] [--Ascending_RankOrder] [--WeightRanking] [--alpha ALPHA] [--w W]
                     [--thread THREAD]

options:
  -h, --help            show this help message and exit
  --mode MODE           Define mode: intersect or closest: intersect count the number of overlapping elements while closest test the distance. In closest
                        mode if a feature overlap a target the distance is 0, use --exclude_ov to test only for non-overlapping regions
  --input INPUT         Input bed file, must contain 6 or more columns, name and score can be placeholder but score is required in --RankTest mode,
                        strand is used only if some strandess test are requested
  --targets TARGETS     Target bed file(s) (must contain 6 or more columns) to test enrichement against, if multiple files are supplied N independent
                        test against each file are conducted, file names must be comma separated, the name of the file will be use as the name output
  --background BACKGROUND
                        Background bed file (must contain 6 or more columns), should be a superset from wich input bed file is derived
  --randomization RANDOMIZATION
                        Number of randomization, default: 100
  --genome GENOME       Genome fasta file used to retrieve sequence features like AT or GC content and length, needed only for length or AT/GC content
                        tests
  --tmp TMP             Temporary directory for storing intermediate files. Default is current working directory
  --outfile OUTFILE     Full path to the output file to store final results in tab format
  --outdir OUTDIR       Full path to output directory to store tables for plot, it is suggested to use a different directory for each analysis. It will
                        be created
  --orientation ORIENTATION
                        Name of test(s) to be performed: concordant, discordant, strandless, or a combination of them. If multiple tests are required
                        tests names must be comma separated, no space allowed
  --ov_fraction OV_FRACTION
                        Minimum overlap required as a fraction from input BED file to consider 2 features as overlapping. Default is 1E-9 (i.e. 1bp)
  --generate_bg         This option activatates the generation of random bed intervals to test enrichment against, use this instead of background. Use
                        only if background file cannot be used or is not available
  --exclude_intervals EXCLUDE_INTERVALS
                        Exclude regions overlapping with regions in the supplied BED file
  --exclude_ov          Exclude overlapping regions between Input and Target file in closest mode
  --exclude_upstream    Exclude upstream region in closest mode, only for stranded files, not compatible with exclude_downstream
  --exclude_downstream  Exclude downstream region in closest mode, only for stranded files, not compatible with exclude_upstream
  --test_AT_GC          Test AT and GC content
  --test_length         Test feature length
  --GenomicLocalization
                        Test also the genomic localization and enrichment of founded overlaps, i.e TSS,Promoter,exons,introns,UTRs - Available only in
                        intersect mode. Must provide a GTF file to extract genomic regions (--gtf), alternatively directly provide a bed file (--bed)
                        with custom annotations
  --gtf GTF             GTF file, only to test genomic localization of founded overlap, gtf file will be used to create genomic regions: promoter, tss,
                        exons, intron, 3UTR and 5UTR
  --bed BED             BED file, only to test genomic localization of founded overlap, bed file will be used to test enrichment in different genomic
                        regions, annotation must be stored as 4th column in bed file, i.e name field
  --RankTest            Activates the Ranking analyis, require BED to contain numerical value in 4th column
  --Ascending_RankOrder
                        Activate the Sort Ascending in RankTest analysis
  --WeightRanking       Weight the ranking test, this is done by increase or decrease the score value in the BED file based on their relative rank and/or
                        distance and/or fractional overlap
  --alpha ALPHA         Relative Influence of the overlap fraction/distance (with respect to ranking) in weightRanked test, only if --WeightRanking is
                        active, must be between 0 and 1
  --w W                 Strength of the Weight for the ranking test, only if --WeightRanking is active, must be between 0 and 1
  --thread THREAD       Number of Threads for parallel computation
```

# How to plot results?
ProOvErlap supports the creation of two main types of graphical outputs (although you may also perform your own plots, as all data are saved to files). The first one is a density plot (generated by the Density_plot.R script), which shows how far the obtained results deviate from what would be expected by chance. Moreover, ProOvErlap also creates heatmaps of the Z-score for each target and, optionally, genomic regions or custom regions, using the Heatmap.R script.
If RankTest is active the plots must be created using the RankPlot.R

```
Density_plot.R: Required arguments: 
input_table: the main output of prooverlap.py "ex: Results.txt",
randomizations: auto generated output of prooverlap.py containing the randomization table "ex: Tables.txt",
test: mode used in prooverlap.py, it must be intersect or closest (default: intersect)
outfile: name of the suffix of output file (default: Density_plot)
format: format used to save the output file, could be png, pdf or svg (default: png)

Heatmap.R: Required arguments:
input_table: main output of prooverlap.py when the option "GenomicLocalization" is set
outfile: name of output file (default = "Heatmap")
format: format used to save the output file, could be png, pdf or svg (default: png)
title: title of the plot (default: "")
```

# Development 
ProOvErlap was developed by Nicolò Gualandi (former post-doc in the Laboratory of Prof. Claudio Brancolini @ UniUd) and Alessio Bertozzo (PhD student in the Laboratory of Prof. Claudio Brancolini @ UniUd), under the supervision of Prof. Claudio Brancolini (Professor of Cell Biology, Department of Medicine, Università degli Studi di Udine, https://people.uniud.it/page/claudio.brancolini)  

ProOvErlap is actively being improved. If you would like to contribute, we welcome your comments and feedback.  
