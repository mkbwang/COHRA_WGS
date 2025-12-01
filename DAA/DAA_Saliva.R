
library(ADAPT)
library(phyloseq)
library(ggplot2)
library(ggrepel)
library(dplyr)
rm(list=ls())



metadata_saliva <- read.table("metadata/metadata_yr1_imputed.tsv", sep="\t", header=1)
rownames(metadata_saliva) <- metadata_saliva$BabySubjectID
saliva_taxa_count <- read.csv("counts_cleaning/saliva_taxa_count_subset.csv",
                              row.names=1) |> t() |> as.data.frame()
saliva_taxa_count_corrected <- read.table("counts_cleaning/saliva_taxa_count_subset_corrected.tsv",
                                          sep='\t', header=1, row.names=1) 
rownames(saliva_taxa_count) <- rownames(saliva_taxa_count_corrected) <-
  gsub("-5", "", rownames(saliva_taxa_count_corrected))
colnames(saliva_taxa_count_corrected) <- colnames(saliva_taxa_count)
taxa_names <- sapply(colnames(saliva_taxa_count), function(longname){
  
  species_name <- strsplit(longname, split="[|]")[[1]][7]
  species_name <- gsub("s__", "", species_name)
  return(species_name)
  
})
colnames(saliva_taxa_count_corrected) <- colnames(saliva_taxa_count) <- unname(taxa_names)
metadata_saliva <- metadata_saliva[rownames(saliva_taxa_count_corrected), ]


# remove those that do not have genus names
ggbs <- grepl("GGB", colnames(saliva_taxa_count_corrected))
saliva_taxa_count <- saliva_taxa_count[, !ggbs]
saliva_taxa_count_corrected <- saliva_taxa_count_corrected[, !ggbs]

metadata_saliva$Case_status <- as.character(metadata_saliva$Case_status)

phyobj_taxa <- phyloseq(otu_table(saliva_taxa_count, taxa_are_rows = FALSE),
                        sample_data(metadata_saliva))

phyobj_taxa_corrected <- phyloseq(otu_table(saliva_taxa_count_corrected, taxa_are_rows = FALSE),
                        sample_data(metadata_saliva))

## first test presence absence
presence_absence <- 1*(saliva_taxa_count_corrected > 0)
fisher_pvals <- rep(0, ncol(saliva_taxa_count_corrected))
for (j in 1:ncol(saliva_taxa_count_corrected)){
  ctable <- table(metadata_saliva$Case_status, presence_absence[, j])
  if (dim(ctable)[2] == 2){
    test_result <- fisher.test(ctable)  
    fisher_pvals[j] <- test_result$p.value
  }else{
    fisher_pvals[j] <- 1
  }
}
fisher_pvals_adjusted <- p.adjust(fisher_pvals, method="BH")


censor_value <- min(saliva_taxa_count_corrected[saliva_taxa_count_corrected > 0])
DAA_taxa_corrected <- adapt(input_data=phyobj_taxa_corrected,
                  cond.var="Case_status",base.cond="0",
                  censor=censor_value, prev.filter=0.05)
DAA_taxa_corrected_table <- DAA_taxa_corrected@details
DAA_taxa_corrected_table$fisher_pvals <- fisher_pvals
DAA_taxa_corrected_table$fisher_pvals_adjusted <- fisher_pvals_adjusted
write.csv(DAA_taxa_corrected_table, "DAA/DAA_taxa_saliva.csv",
          row.names=FALSE)


####### KEGG ##########

saliva_ko_count <- read.csv("counts_cleaning/saliva_ko_abundance_subset.csv",
                            row.names=1) |> t() |> as.data.frame()
saliva_ko_count_corrected <- read.table("counts_cleaning/saliva_ko_abundance_subset_corrected.tsv",
                                        sep='\t', header=1, row.names=1)

rownames(saliva_ko_count) <- rownames(saliva_ko_count_corrected) <- rownames(saliva_taxa_count)


presence_absence <- 1*(saliva_ko_count_corrected > 0)
fisher_pvals <- rep(0, ncol(saliva_ko_count_corrected))
for (j in 1:ncol(saliva_ko_count_corrected)){
  ctable <- table(metadata_saliva$Case_status, presence_absence[, j])
  if (dim(ctable)[2] == 2){
    test_result <- fisher.test(ctable)  
    fisher_pvals[j] <- test_result$p.value
  }else{
    fisher_pvals[j] <- 1
  }
}
fisher_pvals_adjusted <- p.adjust(fisher_pvals, method="BH")

phyobj_ko_corrected <- phyloseq(otu_table(saliva_ko_count_corrected, taxa_are_rows = FALSE),
                                sample_data(metadata_saliva))

censor_value <- min(saliva_ko_count_corrected[saliva_ko_count_corrected > 0])
DAA_ko_corrected <- adapt(input_data=phyobj_ko_corrected,
                            cond.var="Case_status",base.cond="0",
                            censor=censor_value, prev.filter=0.05)
DAA_ko_corrected_table <- DAA_ko_corrected@details
DAA_ko_corrected_table$fisher_pvals <- fisher_pvals
DAA_ko_corrected_table$fisher_pvals_adjusted <- fisher_pvals_adjusted
write.csv(DAA_ko_corrected_table, "DAA/DAA_ko_saliva.csv",
          row.names=FALSE)


