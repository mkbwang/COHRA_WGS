rm(list=ls())

library(dplyr)
library(purrr)

ID_matching <- read.csv(file.path("metadata", "synonym_IDs.csv"))
metadata_yr1 <- read.csv(file.path("metadata", "metadata_yr1.csv"))
metadata_yr2 <- read.csv(file.path("metadata", "metadata_yr2.csv"))

# taxa counts

load_counts <- function(fname, folder){
  
  sampleID <- strsplit(fname, split="_")[[1]][1]
  count_table_indv <- read.table(file.path(folder, "metaphlan_output", fname),
                                 sep="", header=FALSE)
  colnames(count_table_indv) <- c("Taxonomy", "ID", "Relabd", "RPK", "Count")
  levels <- sapply(count_table_indv$Taxonomy, function(fulltax) strsplit(fulltax, split="[|]")[[1]] |> length())
  
  count_species <- count_table_indv[levels == 7, ] |> select(Taxonomy, Count)
  colnames(count_species)[2] <- sampleID
  
  return(count_species)
}


combine_counts <- function(folder){
  
  all_files <- list.files(file.path(folder, "metaphlan_output"))
  indv_file_filter <- grepl("cat", all_files)
  indv_files <- all_files[indv_file_filter]
  
  
  taxa_counts <- lapply(indv_files, load_counts, folder=folder)
  combined_taxa_counts <- reduce(taxa_counts, full_join, by = "Taxonomy")
  combined_taxa_counts[is.na(combined_taxa_counts)] <- 0
  
  return(combine_taxa_counts)
  write.table(combined_taxa_counts, 
              file=file.path(folder, "metaphlan_output", "joint_taxonomic_counts.tsv"),
              sep='\t', row.names = FALSE, quote=FALSE)
}

combine_counts("saliva_preprocessing")
combine_counts("plaque_preprocessing")


# pathway abundances

load_pathabundance <- function(fname, folder){
  
  sampleID <- strsplit(fname, split="_")[[1]][1]
  abundance_table_indv <- read.table(file.path(folder, "humann_output", fname),
                                 sep="\t", header=FALSE)
  colnames(abundance_table_indv) <- c("Pathway", sampleID)
  pathway_names <- abundance_table_indv$Pathway
  pathway_filter <- !grepl("[|]", pathway_names)
  
  abundance_table_filtered <- abundance_table_indv[pathway_filter, ]
  
  return(abundance_table_filtered)
  
}

combine_abundance <- function(folder){

  all_files <- list.files(file.path(folder, "humann_output"))
  indv_file_filter <- grepl("cat", all_files) & grepl("abundance", all_files)
  indv_files <- all_files[indv_file_filter]
  
  
  pathway_abundances <- lapply(indv_files, load_pathabundance, folder=folder)
  combined_pathway <- reduce(pathway_abundances, full_join, by = "Pathway")
  combined_pathway[is.na(combined_pathway)] <- 0
  
  return(combined_pathway)
  
}

combined_pathway_abundances <- combine_abundance("plaque_preprocessing")
write.table(combined_pathway_abundances, 
            file=file.path("plaque_preprocessing", "humann_output", "joint_pathway_abundances.tsv"),
            sep='\t', row.names = FALSE, quote=FALSE)


combined_pathway_abundances <- combine_abundance("saliva_preprocessing")
write.table(combined_pathway_abundances, 
            file=file.path("saliva_preprocessing", "humann_output", "joint_pathway_abundances.tsv"),
            sep='\t', row.names = FALSE, quote=FALSE)






