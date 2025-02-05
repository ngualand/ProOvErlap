suppressMessages({

library(tidyverse)
library(argparse)

parser <- ArgumentParser()

parser$add_argument("--input_table", help = "output of EnrichaRd with Genomic features distrbution")
parser$add_argument("--outfile", help = "Name for output file, default: Plot", default = "Plot")
parser$add_argument("--format", help = "Format of plot file, default: png", default = "png")

args <- parser$parse_args()


Tab = read_delim(args$input_table, col_names = T)

# Sistema la tabella
Tab_for_plot = Tab |> 
  dplyr::filter(str_detect(Target, "\\|\\|")) |> 
  #dplyr::mutate(Zscore = ifelse(P.value > 0.05, "NA", Zscore)) |> 
  dplyr::mutate(Region = str_replace(Target, ".*\\|\\|", ""), 
                Sample = str_replace(Target, "\\|\\|.*", "")) |> 
  dplyr::select(Zscore, Region, Sample) |> 
  dplyr::mutate(Region = str_to_upper(Region)) |> 
  dplyr::mutate(Region = fct_relevel(Region, c("PROMOTER", "TSS", "5UTR", "EXON", "INTRON", "3UTR", "INTERGENIC")))

#Converte valori dello Zscore in numerici
Tab_for_plot$Zscore = as.numeric(Tab_for_plot$Zscore)

# Heatmap
P <- ggplot(Tab_for_plot, aes(y = Sample, x = Region, fill = Zscore)) +
  geom_tile(color="black", linewidth = 0.5) +
  geom_text(aes(label = round(Zscore, digits = 2)), size = 4) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", midpoint = 0, na.value = "lightgrey") +
  theme_minimal() +
  theme(
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

ggsave(P, filename = paste(args$outfile,args$format, sep = "."), bg = "white", dpi = 300, width = 8, height = 3)
  
})
