
library(ADAPT)
library(phyloseq)
library(dplyr)
rm(list=ls())

folder <- "counts_cleaning/strata"

metadata_saliva_yr1 <- read.table(file.path(folder, "metadata_saliva_yr1.tsv"), sep='\t')


## DAA of taxa
saliva_taxa_count_yr1 <- read.table(file.path(folder, "saliva_taxa_counts_yr1.tsv"), sep='\t')
taxa_names <- colnames(saliva_taxa_count_yr1)
noname <- grepl("GGB", taxa_names) | grepl("phylum", taxa_names)
subset_taxa_names <- taxa_names[!noname]
saliva_taxa_count_yr1 <- saliva_taxa_count_yr1[, subset_taxa_names]

min_pos_value <- min(saliva_taxa_count_yr1[saliva_taxa_count_yr1 > 0])
saliva_taxa_phyloseq <- phyloseq(otu_table(saliva_taxa_count_yr1, taxa_are_rows=FALSE),
                                 sample_data(metadata_saliva_yr1))
saliva_DA_taxa_result <- adapt(saliva_taxa_phyloseq, cond.var="Case_status",
                               censor=1, prev.filter=0.1)
plot(saliva_DA_taxa_result, n.label=10)

saliva_DA_taxa_table <- summary(saliva_DA_taxa_result) |> arrange(pval)
write.table(saliva_DA_taxa_table, "DAA/yr1/saliva/taxa/saliva_DA_taxa_table.tsv",
            sep='\t', quote=FALSE, row.names=FALSE)



## DAA of KEGG
metabolism_KOs <- read.csv("counts_cleaning/metabolism_KOs.csv")
metabolism_KOs <- metabolism_KOs$KEGG
saliva_ko_count_yr1 <- read.table(file.path(folder, "saliva_ko_counts_yr1.tsv"), sep='\t')
subset_kegg <- intersect(metabolism_KOs, colnames(saliva_ko_count_yr1))
saliva_ko_count_yr1 <- saliva_ko_count_yr1[, subset_kegg]
min_pos_value <- min(saliva_ko_count_yr1[saliva_ko_count_yr1 > 0])
saliva_ko_phyloseq <- phyloseq(otu_table(saliva_ko_count_yr1, taxa_are_rows = F),
                               sample_data(metadata_saliva_yr1))
saliva_DA_ko_result <- adapt(saliva_ko_phyloseq, cond.var="Case_status",
                               censor=min_pos_value, prev.filter=0.1)
plot(saliva_DA_ko_result, n.label=20)
saliva_DA_ko_table <- summary(saliva_DA_ko_result) |> arrange(pval)
write.table(saliva_DA_ko_table, "DAA/yr1/saliva/ko/saliva_DA_ko_table.tsv",
            sep='\t', quote=FALSE, row.names=FALSE)


## DAA of uniref
saliva_uniref_count_yr1 <- read.table(file.path(folder, "saliva_uniref_counts_yr1.tsv"), sep='\t')
min_pos_value <- min(saliva_uniref_count_yr1[saliva_uniref_count_yr1 > 0])
saliva_uniref_phyloseq <- phyloseq(otu_table(saliva_uniref_count_yr1, taxa_are_rows = F),
                                   sample_data(metadata_saliva_yr1))
saliva_DA_uniref_result <- adapt(saliva_uniref_phyloseq, cond.var="Case_status",
                                 censor=min_pos_value, prev.filter=0.1)
plot(saliva_DA_uniref_result, n.label=20)
saliva_DA_uniref_table <- summary(saliva_DA_uniref_result) |> arrange(pval)
write.table(saliva_DA_uniref_table, "DAA/yr1/saliva/saliva_DA_uniref_table.tsv",
            sep='\t', quote=FALSE, row.names=FALSE)

