
library(ADAPT)
library(phyloseq)
library(dplyr)
rm(list=ls())

folder <- "counts_cleaning/strata"

metadata_plaque_yr1 <- read.table(file.path(folder, "metadata_plaque_yr1.tsv"), sep='\t')


## DAA of taxa
plaque_taxa_count_yr1 <- read.table(file.path(folder, "plaque_taxa_counts_yr1.tsv"), sep='\t')
taxa_names <- colnames(plaque_taxa_count_yr1)
noname <- grepl("GGB", taxa_names) | grepl("phylum", taxa_names)
subset_taxa_names <- taxa_names[!noname]
plaque_taxa_count_yr1 <- plaque_taxa_count_yr1[, subset_taxa_names]

min_pos_value <- min(plaque_taxa_count_yr1[plaque_taxa_count_yr1 > 0])
plaque_taxa_phyloseq <- phyloseq(otu_table(plaque_taxa_count_yr1, taxa_are_rows=FALSE),
                                 sample_data(metadata_plaque_yr1))
plaque_DA_taxa_result <- adapt(plaque_taxa_phyloseq, cond.var="Case_status",
                               censor=min_pos_value, prev.filter=0.1)
plot(plaque_DA_taxa_result, n.label=8)

plaque_DA_taxa_table <- summary(plaque_DA_taxa_result) |> arrange(pval)
write.table(plaque_DA_taxa_table, "DAA/yr1/plaque/taxa/plaque_DA_taxa_table.tsv",
            sep='\t', quote=FALSE, row.names=FALSE)



## DAA of KEGG
metabolism_KOs <- read.csv("counts_cleaning/metabolism_KOs.csv")
metabolism_KOs <- metabolism_KOs$KEGG
plaque_ko_count_yr1 <- read.table(file.path(folder, "plaque_ko_counts_yr1.tsv"), sep='\t')
subset_kegg <- intersect(metabolism_KOs, colnames(plaque_ko_count_yr1))
plaque_ko_count_yr1 <- plaque_ko_count_yr1[, subset_kegg]
min_pos_value <- min(plaque_ko_count_yr1[plaque_ko_count_yr1 > 0])
plaque_ko_phyloseq <- phyloseq(otu_table(plaque_ko_count_yr1, taxa_are_rows = F),
                               sample_data(metadata_plaque_yr1))
plaque_DA_ko_result <- adapt(plaque_ko_phyloseq, cond.var="Case_status",
                               censor=min_pos_value, prev.filter=0.1)
plot(plaque_DA_ko_result, n.label=20)
plaque_DA_ko_table <- summary(plaque_DA_ko_result) |> arrange(pval)
write.table(plaque_DA_ko_table, "DAA/yr1/plaque/ko/plaque_DA_ko_table.tsv",
            sep='\t', quote=FALSE, row.names=FALSE)
