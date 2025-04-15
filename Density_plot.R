suppressMessages ({
library(tidyverse, verbose = F)
library(argparse, verbose = F)

#lettura argomenti
parser <- ArgumentParser()

parser$add_argument("--input_table", help = "Main output of ProOvErlap")
parser$add_argument("--randomizations", help = "Tables_Intersect|Closest.txt")
parser$add_argument("--test", help = "intersect or closest", default = "intersect")
parser$add_argument("--outfile", help = "Name for output file, default: Plot", default = "Density_plot")
parser$add_argument("--format", help = "Format of plot file, default: png", default = "png")
parser$add_argument("--width", help = "4", default = "4", type="numeric")
parser$add_argument("--heigth", help = "3", default = "3", type="numeric")


args <- parser$parse_args()

try(if (!(args$test %in% c("intersect", "closest"))) stop("test argument should be intersect or closest"))


if (args$test == "intersect"){
	Xlab = "Number of intersection"
}
if (args$test == "closest"){
        Xlab = "Mean distance"
}
#Legge la tabella con le randomizzazioni e seleziona solo la colonna Count
randomizations <- read_delim(args$randomizations, col_names = T)

#legge output principale di EnrichaRd.py
Results <- read_delim(args$input_table, col_names = T)


#Density plot
makeDensityplot <- function(Zscore, Type, Pvalue, Real, Random, Name, Outfile, Format, width, heigth) {

    #seleziona le righe del df con le randomizzazioni solo per il target della corrispondente riga di output
    random_df_name <- randomizations %>%
        dplyr::filter(Target == !!Name & Name == !!Type & Type != "Real") #qua secondo me sono da sistemare

    Random = mean(random_df_name$Count)

    NPerm = nrow(random_df_name)
    
    Plot <- random_df_name %>%
    dplyr::rename("Simulated_overlaps" = "Count") %>%
    ggplot() +
    geom_density(aes(Simulated_overlaps,
                     after_stat(count)),
                 adjust = 2,
                 linewidth = 1) +
    geom_vline(xintercept = Random,
               colour = "black",
               linewidth = 0.5,
	       linetype = "dashed") +
    geom_vline(xintercept = Real,
               colour = "#ef233c",
               linetype = "dashed",
               linewidth = 1) +
    labs(x = NULL,
         title = str_c(Name," ",Type),
         subtitle = str_c(
                        "Randomizations = ", NPerm,"\n",
                        "Z-score = ", round(Zscore, digits = 2),"\n",
                        "p-value = ", formatC(Pvalue, format = "e", digits = 2))) +
    ylab("Density") +
    xlab(Xlab) +
    theme_bw(base_size = 10)

ggsave(plot = Plot, filename = str_c(Outfile,"_",Name,"_",Type,".",Format), dpi = 300, height = heigth, width = width) 
}    

apply(Results, 1, function(x) makeDensityplot(as.numeric(x["Zscore"]), x["Type"], as.numeric(x["P.value"]), as.numeric(x["Real"]), as.numeric(x["Random"]), x["Target"], args$outfile, args$format, args$width, args$heigth))

})
