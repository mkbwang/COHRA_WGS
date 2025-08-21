
library(dplyr)
library(ggplot2)

rm(list=ls())
# remove samples with low library sizes and remove rare features

saliva_taxa_count <- read.table("saliva_preprocessing/metaphlan_output/joint_taxonomic_counts.tsv",
                                sep='\t')
saliva_ko_abundance <- read.table("saliva_preprocessing/humann_output/joint_ko_table_concise.tsv",
                                  sep='\t')

colnames(saliva_taxa_count) <- gsub("[.]", "-", colnames(saliva_taxa_count))
colnames(saliva_taxa_count) <- gsub("X", "", colnames(saliva_taxa_count))
colnames(saliva_ko_abundance) <- gsub("[.]", "-", colnames(saliva_ko_abundance))
colnames(saliva_ko_abundance) <- gsub("X", "", colnames(saliva_ko_abundance))


yr1_filter <- grepl("-5", colnames(saliva_taxa_count))
saliva_taxa_count <- saliva_taxa_count[, yr1_filter]
yr1_filter <- grepl("-5", colnames(saliva_ko_abundance))
saliva_ko_abundance <- saliva_ko_abundance[, yr1_filter]


## filter samples based on library size
saliva_libsizes <- colSums(saliva_taxa_count)
saliva_taxa_count_subset <- saliva_taxa_count[, saliva_libsizes > 5e5]
saliva_ko_abundance_subset <- saliva_ko_abundance[, saliva_libsizes > 5e5]
saliva_ko_abundance_subset <- saliva_ko_abundance_subset[-1, ]


## filter features based on prevalence
prevalence_taxa <- rowMeans(saliva_taxa_count_subset > 0)
saliva_taxa_count_subset <- saliva_taxa_count_subset[prevalence_taxa > 0.05,]
prevalence_ko <- rowMeans(saliva_ko_abundance_subset > 0)
saliva_ko_abundance_subset <- saliva_ko_abundance_subset[prevalence_ko > 0.05, ]


write.csv(saliva_taxa_count_subset, "counts_cleaning/saliva_taxa_count_subset.csv",
          quote=FALSE)

write.csv(saliva_ko_abundance_subset, "counts_cleaning/saliva_ko_abundance_subset.csv",
          quote=FALSE)


# Do the same thing for uniref90s

saliva_uniref90_abundance <- read.table("saliva_preprocessing/humann_output/joint_uniref90_subset.tsv",
                                        sep='\t')

colnames(saliva_uniref90_abundance) <- gsub("[.]", "-", colnames(saliva_uniref90_abundance))
colnames(saliva_uniref90_abundance) <- gsub("X", "", colnames(saliva_uniref90_abundance))


yr1_filter <- grepl("-5", colnames(saliva_uniref90_abundance))
saliva_uniref90_abundance <- saliva_uniref90_abundance[, yr1_filter]

saliva_uniref90_abundance_subset <- saliva_uniref90_abundance[, saliva_libsizes > 5e5]
prevalence_uniref90 <- rowMeans(saliva_uniref90_abundance_subset > 0)
saliva_uniref90_abundance_subset <- saliva_uniref90_abundance_subset[prevalence_uniref90 > 0.05, ]

write.csv(saliva_uniref90_abundance_subset, "counts_cleaning/saliva_uniref90_abundance_subset.csv",
          quote=FALSE)


