suppressMessages({

library(tidyverse)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--input_table", help = "output of EnrichaRd with Genomic features distrbution")
parser$add_argument("--outfile", help = "Name for output file directory, default: Heatmap", default = "Heatmap")
parser$add_argument("--format", help = "Format of plot file, default: png", default = "png")
parser$add_argument("--title", help = "title of the plot", default = "")
parser$add_argument("--width", help = "Plot Width", default = "4", type="numeric")
parser$add_argument("--heigth", help = "Plot Heigth", default = "3", type="numeric")

args <- parser$parse_args()

Tab = read_delim(args$input_table, col_names = T)
# Adjust table
Tab_for_plot = Tab |> 
  dplyr::filter(str_detect(Target, "\\|\\|")) |> 
  #dplyr::mutate(Zscore = ifelse(P.value > 0.05, "NA", Zscore)) |> 
  dplyr::mutate(Region = str_replace(Target, ".*\\|\\|", ""), 
                Sample = str_replace(Target, "\\|\\|.*", "")) |> 
  dplyr::select(Zscore, Region, Sample) |> 
  dplyr::mutate(Region = str_to_upper(Region)) |> 
  dplyr::mutate(Region = fct_relevel(Region, c("PROMOTER", "TSS", "5UTR", "EXON", "INTRON", "3UTR", "INTERGENIC")))

Tab_for_plot$Zscore = as.numeric(Tab_for_plot$Zscore)



# Heatmap
P <- ggplot(Tab_for_plot, aes(x = Sample, y = Region, fill = Zscore)) +
  geom_tile(color="black", linewidth = 0.5) +
  geom_text(aes(label = round(Zscore, digits = 2)), size = 5) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", midpoint = 0, na.value = "lightgrey") +
  theme_minimal() +
  labs(title = args$title) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),        
    panel.border = element_blank(),    
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90),
    panel.background = element_rect(
      fill = "white", color = NA
    )
  ) +
  coord_fixed()

if (args$format == "svg") {
  ggsave(P, filename = str_c(args$outfile,".",args$format), dpi = 300, width = args$width, height = args$heigth)
} else {
  ggsave(P, filename = str_c(args$outfile,".",args$format), bg = "white", dpi = 300, width = args$width, height = args$heigth)
}

  
})
