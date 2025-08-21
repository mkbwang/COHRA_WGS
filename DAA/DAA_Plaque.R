
library(ADAPT)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(dplyr)
rm(list=ls())



metadata_plaque <- read.table("metadata/metadata_yr1_imputed.tsv", sep="\t", header=1)
rownames(metadata_plaque) <- metadata_plaque$BabySubjectID
plaque_taxa_count <- read.csv("counts_cleaning/plaque_taxa_count_subset.csv",
                              row.names=1) |> t() |> as.data.frame()
plaque_taxa_count_corrected <- read.table("counts_cleaning/plaque_taxa_count_subset_corrected.tsv",
                                          sep='\t', header=1, row.names=1) 
rownames(plaque_taxa_count) <- rownames(plaque_taxa_count_corrected) <-
  gsub("-5", "", rownames(plaque_taxa_count_corrected))
colnames(plaque_taxa_count_corrected) <- colnames(plaque_taxa_count)
taxa_names <- sapply(colnames(plaque_taxa_count), function(longname){
  
  species_name <- strsplit(longname, split="[|]")[[1]][7]
  species_name <- gsub("s__", "", species_name)
  return(species_name)
  
})
colnames(plaque_taxa_count_corrected) <- colnames(plaque_taxa_count) <- unname(taxa_names)
metadata_plaque <- metadata_plaque[rownames(plaque_taxa_count_corrected), ]


# remove those that do not have genus names
ggbs <- grepl("GGB", colnames(plaque_taxa_count_corrected))
plaque_taxa_count <- plaque_taxa_count[, !ggbs]
plaque_taxa_count_corrected <- plaque_taxa_count_corrected[, !ggbs]

metadata_plaque$Case_status <- as.character(metadata_plaque$Case_status)

phyobj_taxa <- phyloseq(otu_table(plaque_taxa_count, taxa_are_rows = FALSE),
                        sample_data(metadata_plaque))

phyobj_taxa_corrected <- phyloseq(otu_table(plaque_taxa_count_corrected, taxa_are_rows = FALSE),
                                  sample_data(metadata_plaque))

plot_generation <- function(result, title_x, title_main, labeffsize=5, labpval=5){
  
  result$neglog10pval <- -log10(result$pval)
  effsizes <- sort(result$log10foldchange)
  pvals <- sort(result$pval)
  
  # pick the features to color and label
  tolabel <- (result$log10foldchange <= effsizes[labeffsize] | 
                result$log10foldchange >= tail(effsizes, labeffsize)[1] |
                result$pval <= pvals[labpval])
  
  result$islabel <- tolabel
  result$Label  <- ""
  result$Label[tolabel] <- result$Taxa[tolabel]
  
  generated_plot <- ggplot(result, aes(x=.data$log10foldchange, y=.data$neglog10pval))+
    geom_point(alpha=0.8, aes(color=.data$islabel)) +
    xlab(title_x) + ylab("-Log10 p-value") + theme_bw() + 
    theme(legend.position="none", axis.title=element_text(size=10), 
          axis.text=element_text(size=10)) + 
    scale_color_manual(values=c("#616161", "#ff0066")) +
    geom_vline(xintercept=0, linetype="dashed", color = "blue") +
    geom_label_repel(aes(label = .data$Label),
                     size=2,
                     max.overlaps = 20,
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50')+
    ggtitle(title_main)
  
  return(generated_plot)
  
}


DAA_taxa <- adapt(input_data=phyobj_taxa,
                  cond.var="Case_status",base.cond="0",
                  censor=1, prev.filter=0.05)
DAA_taxa_table <- DAA_taxa@details 
DAA_taxa_plot <- plot_generation(result=DAA_taxa_table,
                                 title_x="Log10 fold change (Case vs Control)",
                                 title_main="Taxa Differential Abundance in Plaque",
                                 labeffsize=3, labpval=3)


# DAA_taxa_corrected <- adapt(input_data=phyobj_taxa_corrected,
#                   cond.var="Case_status",base.cond="0",
#                   censor=1, prev.filter=0.05)
# DAA_taxa_corrected_table <- DAA_taxa_corrected@details 
# DAA_taxa_corrected_plot <- plot_generation(result=DAA_taxa_corrected_table,
#                                  title_x="Log10 fold change (Case vs Control)",
#                                  title_main="Taxa Differential Abundance in plaque",
#                                  labeffsize=5, labpval=5)



plaque_ko_count <- read.csv("counts_cleaning/plaque_ko_abundance_subset.csv",
                            row.names=1) |> t() |> as.data.frame()

rownames(plaque_ko_count) <- rownames(plaque_taxa_count)

phyobj_ko <- phyloseq(otu_table(plaque_ko_count, taxa_are_rows = FALSE),
                      sample_data(metadata_plaque))

DAA_ko <- adapt(input_data=phyobj_ko,
                cond.var="Case_status",base.cond="0",
                censor=min(plaque_ko_count[plaque_ko_count > 0]), 
                prev.filter=0.05)

DAA_ko_table <- DAA_ko@details 
DAA_ko_plot <- plot_generation(result=DAA_ko_table,
                               title_x="Log10 fold change (Case vs Control)",
                               title_main="KEGG Differential Abundance in Plaque",
                               labeffsize=5, labpval=5)

library(patchwork)
combined_plots <- wrap_plots(DAA_taxa_plot, DAA_ko_plot,
                             ncol=2)

ggsave(filename="DAA/plaque_plot.pdf",
       width=12, height=5,
       plot=combined_plots)

write.csv(DAA_taxa_table, "DAA/DAA_taxa_plaque.csv",
          row.names=FALSE)
write.csv(DAA_ko_table, "DAA/DAA_ko_plaque.csv",
          row.names=FALSE)


