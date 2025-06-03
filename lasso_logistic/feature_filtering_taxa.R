

rm(list=ls())
saliva_taxa_counts_yr1 <- read.table("counts_cleaning/strata/saliva_taxa_counts_yr1.tsv",
                                     header=TRUE, sep="\t", row.names=1)
taxa_saliva <- colnames(saliva_taxa_counts_yr1)
plaque_taxa_counts_yr1 <- read.table("counts_cleaning/strata/plaque_taxa_counts_yr1.tsv",
                                   header=TRUE, sep="\t", row.names=1)
taxa_plaque <- colnames(plaque_taxa_counts_yr1)


# remove taxa that does not have genus names
saliva_filter <- grepl("GGB", taxa_saliva)
plaque_filter <- grepl("GGB", taxa_plaque)

subset_saliva_taxa <-  taxa_saliva[!saliva_filter]
subset_plaque_taxa <- taxa_plaque[!plaque_filter]

metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
metadata_plaque_yr1 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                                  header=T, sep='\t')

library(phyloseq)
library(ADAPT)
library(dplyr)

saliva_taxa_counts_subset <- saliva_taxa_counts_yr1[, subset_saliva_taxa]
prevalence_saliva <- colMeans(saliva_taxa_counts_subset > 0)
saliva_taxa_counts_subset <- saliva_taxa_counts_subset[, prevalence_saliva > 0.1]
write.csv(saliva_taxa_counts_subset,"lasso_logistic/feature_filter_taxa/saliva_taxa_subset.csv",
          row.names=T)


plaque_taxa_counts_subset <- plaque_taxa_counts_yr1[, subset_plaque_taxa]
prevalence_plaque <- colMeans(plaque_taxa_counts_subset > 0)
plaque_taxa_counts_subset <- plaque_taxa_counts_subset[, prevalence_plaque > 0.1]
write.csv(plaque_taxa_counts_subset, "lasso_logistic/feature_filter_taxa/plaque_taxa_subset.csv",
          row.names=T)






