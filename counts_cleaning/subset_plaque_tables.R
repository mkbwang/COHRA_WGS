
library(dplyr)
library(ggplot2)

rm(list=ls())
# remove samples with low library sizes and remove rare features

plaque_taxa_count <- read.table("plaque_preprocessing/metaphlan_output/joint_taxonomic_counts.tsv",
                                sep='\t')
plaque_ko_abundance <- read.table("plaque_preprocessing/humann_output/joint_ko_table_concise.tsv",
                                  sep='\t')

metadata <- read.csv(file.path("metadata", "metadata_yr1.csv"))

colnames(plaque_taxa_count) <- gsub("[.]", "-", colnames(plaque_taxa_count))
colnames(plaque_taxa_count) <- gsub("X", "", colnames(plaque_taxa_count))
colnames(plaque_ko_abundance) <- gsub("[.]", "-", colnames(plaque_ko_abundance))
colnames(plaque_ko_abundance) <- gsub("X", "", colnames(plaque_ko_abundance))


yr1_filter <- grepl("-5", colnames(plaque_taxa_count))
plaque_taxa_count <- plaque_taxa_count[, yr1_filter]
yr1_filter <- grepl("-5", colnames(plaque_ko_abundance))
plaque_ko_abundance <- plaque_ko_abundance[, yr1_filter]

sample_ids <- gsub("-5", "", colnames(plaque_taxa_count))
metadata <- metadata %>% filter(BabySubjectID %in% sample_ids)


## filter samples based on library size
plaque_libsizes <- colSums(plaque_taxa_count)
plaque_taxa_count_subset <- plaque_taxa_count[, plaque_libsizes > 5e5]
plaque_ko_abundance_subset <- plaque_ko_abundance[, plaque_libsizes > 5e5]
plaque_ko_abundance_subset <- plaque_ko_abundance_subset[-1, ]
sample_ids <- sample_ids[plaque_libsizes > 5e5]
metadata <- metadata %>% filter(BabySubjectID %in% sample_ids)


## filter features based on prevalence
prevalence_taxa <- rowMeans(plaque_taxa_count_subset > 0)
plaque_taxa_count_subset <- plaque_taxa_count_subset[prevalence_taxa > 0.05,]
prevalence_ko <- rowMeans(plaque_ko_abundance_subset > 0)
plaque_ko_abundance_subset <- plaque_ko_abundance_subset[prevalence_ko > 0.05, ]


write.csv(plaque_taxa_count_subset, "counts_cleaning/plaque_taxa_count_subset.csv",
          quote=FALSE)
write.csv(plaque_ko_abundance_subset, "counts_cleaning/plaque_ko_abundance_subset.csv",
          quote=FALSE)



# Do the same thing for uniref90s

plaque_uniref90_abundance <- read.table("plaque_preprocessing/humann_output/joint_uniref90_subset.tsv",
                                  sep='\t')

colnames(plaque_uniref90_abundance) <- gsub("[.]", "-", colnames(plaque_uniref90_abundance))
colnames(plaque_uniref90_abundance) <- gsub("X", "", colnames(plaque_uniref90_abundance))


yr1_filter <- grepl("-5", colnames(plaque_uniref90_abundance))
plaque_uniref90_abundance <- plaque_uniref90_abundance[, yr1_filter]
plaque_uniref90_abundance_subset <- plaque_uniref90_abundance[, plaque_libsizes > 5e5]
prevalence_uniref90 <- rowMeans(plaque_uniref90_abundance_subset > 0)
plaque_uniref90_abundance_subset <- plaque_uniref90_abundance_subset[prevalence_uniref90 > 0.05, ]
write.csv(plaque_uniref90_abundance_subset, "counts_cleaning/plaque_uniref90_abundance_subset.csv",
          quote=FALSE)



