
rm(list=ls())
library(dplyr)


DAA_saliva_taxa <- read.csv("DAA/DAA_taxa_saliva.csv")


DAA_plaque_taxa <- read.csv("DAA/DAA_taxa_plaque.csv")
DAA_saliva_taxa_subset <- DAA_saliva_taxa %>% filter(prevalence > 0.1)
DAA_plaque_taxa_subset <- DAA_plaque_taxa %>% filter(prevalence > 0.1)


DAA_saliva_ko <- read.csv("DAA/DAA_ko_saliva.csv")
DAA_plaque_ko <- read.csv("DAA/DAA_ko_plaque.csv")
DAA_saliva_ko_subset <- DAA_saliva_ko %>% filter(prevalence > 0.1)
DAA_plaque_ko_subset <- read.csv("DAA/DAA_ko_plaque.csv")
DAA_plaque_ko_subset <- DAA_plaque_ko %>% filter(prevalence > 0.1)


max(DAA_saliva_taxa_subset$log10foldchange)
min(DAA_saliva_ko_subset$log10foldchange)

plot_generation <- function(result, labeffsize=5, labpval=5, xmin=-5, xmax=5, ymax=9){
  
  result$neglog10pval <- -log10(result$pval)
  effsizes <- sort(result$log10foldchange)
  pvals <- sort(result$pval)
  
  # pick the features to color and label
  smallpval <- result$pval <= pvals[labpval]
  smalleffsize <- (result$log10foldchange <= effsizes[labeffsize] & result$pval < 0.05)
  largeeffsize <- (result$log10foldchange >= tail(effsizes, labeffsize)[1] & result$pval < 0.05)
  
  tolabel <- (smallpval| smalleffsize| largeeffsize)

  
  result$islabel <- tolabel
  result$Label  <- ""
  result$Label[tolabel] <- result$Taxa[tolabel]
  
  generated_plot <- ggplot(result, aes(x=.data$log10foldchange, y=.data$neglog10pval))+
    geom_point(alpha=0.8, aes(color=.data$islabel)) +
    ylab("-Log10 p-value") + theme_bw() + 
    scale_x_continuous(limits = c(xmin, xmax), breaks = seq(xmin, xmax))+
    scale_y_continuous(limits=c(0, ymax), breaks=seq(0, ymax))+
    theme(legend.position="none", axis.title=element_text(size=10), 
          axis.title.x=element_blank(),
          axis.text=element_text(size=10)) + 
    scale_color_manual(values=c("#616161", "#ff0066")) +
    geom_vline(xintercept=0, linetype="dashed", color = "blue") +
    geom_label_repel(aes(label = .data$Label),
                     size=2,
                     max.overlaps = 20,
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50')
  
  return(generated_plot)
  
}

DAA_saliva_taxa_plot <- plot_generation(DAA_saliva_taxa_subset, xmin=-5, xmax=5, ymax=9)
DAA_saliva_ko_plot <- plot_generation(DAA_saliva_ko_subset, xmin=-5, xmax=5, ymax=9)


DAA_plaque_taxa_plot <- plot_generation(DAA_plaque_taxa_subset, xmin=-7, xmax=7, ymax=6.5)
DAA_plaque_ko_plot <- plot_generation(DAA_plaque_ko_subset, xmin=-7, xmax=7, ymax=6.5)


