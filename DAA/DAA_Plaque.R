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

## first test presence absence
presence_absence <- 1*(plaque_taxa_count_corrected > 0)
fisher_pvals <- rep(0, ncol(plaque_taxa_count_corrected))
for (j in 1:ncol(plaque_taxa_count_corrected)){
  ctable <- table(metadata_plaque$Case_status, presence_absence[, j])
  if (dim(ctable)[2] == 2){
    test_result <- fisher.test(ctable)  
    fisher_pvals[j] <- test_result$p.value
  }else{
    fisher_pvals[j] <- 1
  }
}
fisher_pvals_adjusted <- p.adjust(fisher_pvals, method="BH")


censor_value <- min(plaque_taxa_count_corrected[plaque_taxa_count_corrected > 0])
DAA_taxa_corrected <- adapt(input_data=phyobj_taxa_corrected,
                            cond.var="Case_status",base.cond="0",
                            censor=censor_value, prev.filter=0.05)
DAA_taxa_corrected_table <- DAA_taxa_corrected@details
DAA_taxa_corrected_table$fisher_pvals <- fisher_pvals
DAA_taxa_corrected_table$fisher_pvals_adjusted <- fisher_pvals_adjusted
write.csv(DAA_taxa_corrected_table, "DAA/DAA_taxa_plaque.csv",
          row.names=FALSE)


####### KEGG ##########

plaque_ko_count <- read.csv("counts_cleaning/plaque_ko_abundance_subset.csv",
                            row.names=1) |> t() |> as.data.frame()
plaque_ko_count_corrected <- read.table("counts_cleaning/plaque_ko_abundance_subset_corrected.tsv",
                                        sep='\t', header=1, row.names=1)

rownames(plaque_ko_count) <- rownames(plaque_ko_count_corrected) <- rownames(plaque_taxa_count)


presence_absence <- 1*(plaque_ko_count_corrected > 0)
fisher_pvals <- rep(0, ncol(plaque_ko_count_corrected))
for (j in 1:ncol(plaque_ko_count_corrected)){
  ctable <- table(metadata_plaque$Case_status, presence_absence[, j])
  if (dim(ctable)[2] == 2){
    test_result <- fisher.test(ctable)  
    fisher_pvals[j] <- test_result$p.value
  }else{
    fisher_pvals[j] <- 1
  }
}
fisher_pvals_adjusted <- p.adjust(fisher_pvals, method="BH")

phyobj_ko_corrected <- phyloseq(otu_table(plaque_ko_count_corrected, taxa_are_rows = FALSE),
                                sample_data(metadata_plaque))

censor_value <- min(plaque_ko_count_corrected[plaque_ko_count_corrected > 0])
DAA_ko_corrected <- adapt(input_data=phyobj_ko_corrected,
                          cond.var="Case_status",base.cond="0",
                          censor=censor_value, prev.filter=0.05)
DAA_ko_corrected_table <- DAA_ko_corrected@details
DAA_ko_corrected_table$fisher_pvals <- fisher_pvals
DAA_ko_corrected_table$fisher_pvals_adjusted <- fisher_pvals_adjusted
write.csv(DAA_ko_corrected_table, "DAA/DAA_ko_plaque.csv",
          row.names=FALSE)


