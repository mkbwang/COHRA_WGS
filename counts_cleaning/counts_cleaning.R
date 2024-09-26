rm(list=ls())

library(dplyr)
library(purrr)

ID_matching <- read.csv(file.path("metadata", "synonym_IDs.csv"))

ID_matching_yr1 <- ID_matching
ID_matching_yr1$MotherSubjectID <- sprintf("%d-5", ID_matching_yr1$MotherSubjectID)
ID_matching_yr1$BabySubjectID <- sprintf("%d-5", ID_matching_yr1$BabySubjectID)
ID_matching_yr1$AltSubjectID <- sprintf("%s-5", ID_matching_yr1$AltSubjectID)

ID_matching_yr2 <- ID_matching
ID_matching_yr2$MotherSubjectID <- sprintf("%d-7", ID_matching_yr2$MotherSubjectID)
ID_matching_yr2$BabySubjectID <- sprintf("%d-7", ID_matching_yr2$BabySubjectID)
ID_matching_yr2$AltSubjectID <- sprintf("%s-7", ID_matching_yr2$AltSubjectID)
ID_matching <- rbind(ID_matching_yr1, ID_matching_yr2)

metadata_yr1 <- read.csv(file.path("metadata", "metadata_yr1.csv"))
metadata_yr2 <- read.csv(file.path("metadata", "metadata_yr2.csv"))

# combine taxa counts

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



# change the sample names of KO tables


plaque_ko_table <- read.table(file.path("plaque_preprocessing", "humann_output", "joint_ko_table.tsv"),
                              sep='\t', header=T)

columns <- colnames(plaque_ko_table)

## filter out the marginal KO names
all_kos <- plaque_ko_table$Gene.Family
ko_filter <- !grepl("[|]", all_kos)
plaque_ko_table <- plaque_ko_table[ko_filter, ]
rownames(plaque_ko_table) <- all_kos[ko_filter]
plaque_ko_table$Gene.Family <- NULL

columns <- colnames(plaque_ko_table)
columns <- gsub("[A-Za-z]", "", columns)
columns <- gsub("[[:punct:]]+$", "", columns)
columns <- gsub("[.]", "-", columns)
for (j in 1:length(columns)){
  originalID <- columns[j]
  location <- which(ID_matching$AltSubjectID == originalID)
  if (length(location > 0)){
    columns[j] <- ID_matching$BabySubjectID[location]
  }
}

colnames(plaque_ko_table) <- columns
write.table(plaque_ko_table, 
            file=file.path("plaque_preprocessing", "humann_output", "joint_ko_table_concise.tsv"),
            sep='\t', quote=FALSE)








saliva_ko_table <- read.table(file.path("saliva_preprocessing", "humann_output", "joint_ko_table.tsv"),
                              sep='\t', header=T)

columns <- colnames(saliva_ko_table)

# filter out the marginal KO names
all_kos <- saliva_ko_table$Gene.Family
ko_filter <- !grepl("[|]", all_kos)
saliva_ko_table <- saliva_ko_table[ko_filter, ]
rownames(saliva_ko_table) <- all_kos[ko_filter]
saliva_ko_table$Gene.Family <- NULL

columns <- colnames(saliva_ko_table)
columns <- gsub("[A-Za-z]", "", columns)
columns <- gsub("[[:punct:]]+$", "", columns)
columns <- gsub("[.]", "-", columns)
for (j in 1:length(columns)){
  originalID <- columns[j]
  location <- which(ID_matching$AltSubjectID == originalID)
  if (length(location > 0)){
    columns[j] <- ID_matching$BabySubjectID[location]
  }
}

colnames(saliva_ko_table) <- columns
write.table(saliva_ko_table, 
            file=file.path("saliva_preprocessing", "humann_output", "joint_ko_table_concise.tsv"),
            sep='\t', quote=FALSE)




# change the sample names of taxa tables


plaque_taxa_table <- read.table(file.path("plaque_preprocessing", "metaphlan_output", "joint_taxonomic_counts.tsv"),
                                sep='\t', header=T)

rownames(plaque_taxa_table) <- plaque_taxa_table$Taxonomy
plaque_taxa_table$Taxonomy <- NULL
columns <- colnames(plaque_taxa_table)
columns <- gsub("[A-Za-z]", "", columns)
columns <- gsub("[[:punct:]]+$", "", columns)
columns <- gsub("[.]", "-", columns)
for (j in 1:length(columns)){
  originalID <- columns[j]
  location <- which(ID_matching$AltSubjectID == originalID)
  if (length(location > 0)){
    columns[j] <- ID_matching$BabySubjectID[location]
  }
}
colnames(plaque_taxa_table) <- columns
write.table(plaque_taxa_table, 
            file=file.path("plaque_preprocessing", "metaphlan_output", "joint_taxonomic_counts.tsv"),
            sep='\t', quote=FALSE)






saliva_taxa_table <- read.table(file.path("saliva_preprocessing", "metaphlan_output", "joint_taxonomic_counts.tsv"),
                                sep='\t', header=T)

rownames(saliva_taxa_table) <- saliva_taxa_table$Taxonomy
saliva_taxa_table$Taxonomy <- NULL
columns <- colnames(saliva_taxa_table)
columns <- gsub("[A-Za-z]", "", columns)
columns <- gsub("[[:punct:]]+$", "", columns)
columns <- gsub("[.]", "-", columns)
for (j in 1:length(columns)){
  originalID <- columns[j]
  location <- which(ID_matching$AltSubjectID == originalID)
  if (length(location > 0)){
    columns[j] <- ID_matching$BabySubjectID[location]
  }
}
colnames(saliva_taxa_table) <- columns
write.table(saliva_taxa_table, 
            file=file.path("saliva_preprocessing", "metaphlan_output", "joint_taxonomic_counts.tsv"),
            sep='\t', quote=FALSE)



