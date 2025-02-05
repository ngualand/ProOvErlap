# Load required libraries
library(AnnotationHub)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(Biostrings)

#Read Arguments
args <- commandArgs(trailingOnly = TRUE)


### Similar to bedtools subtract.
setGeneric("subtract", signature=c("x", "y"),
           function(x, y, minoverlap=1L, ...) standardGeneric("subtract")
)

### Returns a GenomicRangesList object parallel to 'x'.
setMethod("subtract", c("GenomicRanges", "GenomicRanges"),
          function(x, y, minoverlap=1L, ignore.strand=FALSE)
          {
            y <- reduce(y, ignore.strand=ignore.strand)
            hits <- findOverlaps(x, y, minoverlap=minoverlap,
                                 ignore.strand=ignore.strand)
            ans <- psetdiff(x, extractList(y, as(hits, "IntegerList")))
            mcols(ans) <- mcols(x)
            setNames(ans, names(x))
          }
)

# Function to adjust invalid coordinates
adjust_coordinates <- function(gr) {
  # Ensure start positions are non-negative, gtf is 1 based the -1 is done after by the export funxtion
  start(gr) <- ifelse(start(gr) < 0, 1, start(gr))
  start(gr) <- ifelse(start(gr) == 0, 1, start(gr))
  
  # Ensure end positions are greater than or equal to start positions
  end(gr) <- pmax(end(gr), start(gr))
  
  return(gr)
}

#Read in fasta file and create bed file with whole genome coordinates
fasta_file <- args[3]
genome <- readDNAStringSet(fasta_file)
# Extract sequence names and lengths
seq_names <- unlist(lapply(names(genome), function(x) unlist(strsplit(x, " "))[1]))
seq_lengths <- width(genome)
#Create Bed file with whole genome
genome <- data.frame(
  chr = seq_names,
  start = 0,
  end = as.integer(seq_lengths))

#Load gtf and create TxDb
print("Creating reference annotation BED file... This may take some minutes based on GTF dimension")
txdb = makeTxDbFromGFF(args[1],
                       format = "gtf"
                       )

# Extract Transcripts (TSS is the first base of each transcript)
tss <- promoters(txdb, upstream = 1, downstream = 1)
mcols(tss)$name = "tss"
# Extract exons
exons <- exons(txdb)
mcols(exons)$name = "exon"
# Extract introns
introns <- unlist(intronsByTranscript(txdb))
mcols(introns)$name = "intron"
# Define promoter regions 
promoters <- promoters(txdb, upstream = 1000, downstream = 0)
mcols(promoters)$name = "promoter"
# Extract 5' UTRs
five_utr <- unlist(fiveUTRsByTranscript(txdb))
mcols(five_utr)$name = "5UTR"
# Extract 3' UTRs
three_utr <- unlist(threeUTRsByTranscript(txdb))
mcols(three_utr)$name = "3UTR"

# Combine all feature ranges into a single GRanges object (union of all genomic features)
all_features <- c(exons, introns, promoters, five_utr, three_utr, tss)
# Adjust negative starts
all_features <- adjust_coordinates(all_features)

#Create Genome Granges object for whole genome
genome = makeGRangesFromDataFrame(genome)

#Create intergenic regions as regions without any annotation
intergenic = subtract(genome,all_features)
intergenic = intergenic@unlistData
intergenic$name = "intergenic"

#Re-create all feature file by also adding intergenic
all_features = c(all_features, intergenic)

all_features <- adjust_coordinates(all_features)

#Export and save bed f
export(all_features, con = paste0(args[2],"/","Reference_regions.bed"), format = "bed")


