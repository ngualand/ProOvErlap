library(ggplot2)
library(grid)

#arg1 = type intersect or closest
#arg2 = input table
#arg3 = filename
#arg4 = filetype
#arg5 = width
#arg6 = height

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript RankPlot.R intersect|closest table filename svg|png|jpeg|tif width height")
}

RankPlot <- function(forplot) {
  
  # Convert string inputs to numeric vectors
  ES_values <- as.numeric(forplot$ES_values[!is.na(forplot$ES_values)])
  hit_indices <- as.numeric(forplot$hit_indices)
  null_ES <- as.numeric(forplot$Null_ES)
  ES <- as.numeric(forplot$ES[1])
  Pval <- as.numeric(forplot$Pval[1])  
  
  # Create a data frame for ES values
  es_df <- data.frame(
    Rank = seq_along(ES_values),
    ES = ES_values
  )
  
  # Identify peak ES point
  peak_idx <- which.max(abs(ES_values))
  
  # Plot Running ES
  es_plot <- ggplot(es_df, aes(x = Rank, y = ES)) +
    geom_vline(xintercept = hit_indices, linetype = "dashed", color = "gray", alpha = 0.7) +
    geom_line(color = "blue", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_point(aes(x = peak_idx, y = ES_values[peak_idx]), color = "red", size = 2) +
    labs(
      x = "Ranked Region List",
      y = "Enrichment Score (ES)",
      title = paste("Rank Enrichment Plot\nES =", round(ES, 3), "P =", formatC(Pval, format = "e", digits = 2))
    ) +
    theme_minimal() +
    theme(#panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))  # No extra margin, occupy all space
  
  
  # Create a data frame for null distribution
  null_df <- data.frame(Null_ES = null_ES)
  
  # Plot Null Distribution
  null_plot <- ggplot(null_df, aes(x = Null_ES)) +
    geom_density(fill = "gray", alpha = 0.7, color = "black") +
    geom_vline(xintercept = ES, linetype = "dashed", color = "red") +
    labs(
      x = "Enrichment Score (ES)",
      y = "Frequency",
      title = paste("Null Distribution of ES\nES =", round(ES, 3), "P =", formatC(Pval, format = "e", digits = 2))
    ) +
    theme_minimal() +
    theme(plot.margin = margin(0, 0, 0, 0), 
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA))  # No extra margin, occupy all space
  

  
  return(list(es_plot, null_plot))
}
RankPlot_closest <- function(forplot) {
  
  # Convert string inputs to numeric vectors
  cumulative_real <- as.numeric(forplot$cumulative_real)
  mean_shuffle_cdf <- as.numeric(forplot$mean_shuffle_cdf)
  x <- seq_along(cumulative_real)
  
  # Extract ES and P-value for annotation
  ES <- as.numeric(forplot$ES[1])
  Pval <- as.numeric(forplot$Pval[1])
  
  # Create a data frame for plotting
  df <- data.frame(
    Peak_idx = forplot["Peak_idx"],
    x = x,
    Cumulative_Real = cumulative_real,
    Mean_Shuffle_CDF = mean_shuffle_cdf
  )

  # Plot the cumulative distributions
  p <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = Cumulative_Real, color = "Cumulative Real"), size = 1.2) +
    geom_line(aes(y = Mean_Shuffle_CDF, color = "Mean Cumulative Shuffled"), size = 1, alpha = 0.6) +
    geom_ribbon(aes(ymin = Mean_Shuffle_CDF, ymax = Cumulative_Real), fill = "gray", alpha = 0.3) +
    geom_point(aes(x = Peak_idx, y = cumulative_real[Peak_idx]), color = "red", size = 2) +
    scale_color_manual(values = c("Cumulative Real" = "blue", "Mean Cumulative Shuffled" = "red")) +
    labs(
      x = "Ranked Region list",
      y = "Cumulative Distribution (absolute distances)",
      title = paste("Cumulative Distributions: Real vs Shuffle\nES =", round(ES, 3), "| P =", formatC(Pval, format = "e", digits = 2)),
      color = "Legend"
    ) +
    theme_minimal() +
    theme(
      plot.margin = margin(5, 5, 5, 5, "pt"), 
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 6),   # Smaller legend text
      legend.title = element_text(size = 7),  # Smaller legend title
      legend.key.size = unit(0.5, "lines")    # Adjust size of legend keys
    )
  

}




if (args[1] == "intersect"){
  tab = read.delim(args[2])
  plot = RankPlot(tab)
  ggsave(plot = plot[[1]], filename = paste0(args[3], "_Enrichment", ".", args[4]), 
         width = as.numeric(args[5]), height = as.numeric(args[6]))
  ggsave(plot = plot[[2]], filename = paste0(args[3], "_RandomDist", ".", args[4]), 
         width = as.numeric(args[5]), height = as.numeric(args[6]))
}

if (args[1] == "closest"){
  tab = read.delim(args[2])
  plot = RankPlot_closest(tab)
  ggsave(plot = plot, filename = paste0(args[3], ".", args[4]), 
         width = as.numeric(args[5]), height = as.numeric(args[6]))
}

