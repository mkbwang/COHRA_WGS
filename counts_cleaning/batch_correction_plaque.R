rm(list=ls())
library(dplyr)
library(openxlsx)
library(ConQuR)
library(doParallel)

plaque_batchinfo <- read.xlsx("metadata/Plaque_batchinfo.xlsx", sheet="Sheet1")


synonym_IDs <- read.csv("metadata/synonym_IDs.csv")

plaque_batchinfo$findmatch <- 0


clean_ID <- function(sampleid){
  all_components <- strsplit(sampleid, split="-")[[1]]
  cleaned_ID <- paste(all_components[c(1,2)], collapse="-")
  cleaned_ID
}


for (ID in seq(1, length(synonym_IDs$AltSubjectID))){
  
  original_ID <- synonym_IDs$AltSubjectID[ID]
  replace_ID <- as.character(synonym_IDs$BabySubjectID[ID])
  exist_in_plaque <- grepl(original_ID, plaque_batchinfo$Sample) + 
    grepl(replace_ID, plaque_batchinfo$Sample)
  plaque_batchinfo$findmatch <- plaque_batchinfo$findmatch + exist_in_plaque
  
  if (sum(exist_in_plaque) > 0){ # the ID being searched exist
    replace_loc <- which(exist_in_plaque > 0)
    plaque_batchinfo$Sample[replace_loc] <- gsub(original_ID, replace_ID, 
                                                 plaque_batchinfo$Sample[replace_loc])
  } 
  
}

plaque_batchinfo <- plaque_batchinfo %>% filter(findmatch > 0) %>%
  select(Sample, Batch)
plaque_batchinfo$Sample <- sapply(plaque_batchinfo$Sample, clean_ID)



plaque_taxa_count <- read.table("plaque_preprocessing/metaphlan_output/joint_taxonomic_counts.tsv",
                                sep='\t')
colnames(plaque_taxa_count) <- gsub("[.]", "-", colnames(plaque_taxa_count))
colnames(plaque_taxa_count) <- gsub("X", "", colnames(plaque_taxa_count))
plaque_ko_abundance <- read.table("plaque_preprocessing/humann_output/joint_ko_table_concise.tsv", 
                                  header=T, sep='\t')
colnames(plaque_ko_abundance) <- gsub("[.]", "-", colnames(plaque_ko_abundance))
colnames(plaque_ko_abundance) <- gsub("X", "", colnames(plaque_ko_abundance))
plaque_uniref_abundance <- read.table("plaque_preprocessing/humann_output/joint_uniref90_subset.tsv", 
                                      header=T, sep='\t')
colnames(plaque_uniref_abundance) <- gsub("[.]", "-", colnames(plaque_uniref_abundance))
colnames(plaque_uniref_abundance) <- gsub("X", "", colnames(plaque_uniref_abundance))



## batch correction
## plaque
batch_numbers <- factor(plaque_batchinfo$Batch)
names(batch_numbers) <- plaque_batchinfo$Sample
plaque_taxa_count <- t(plaque_taxa_count)
plaque_ko_abundance <- t(plaque_ko_abundance)
plaque_uniref_abundance <- t(plaque_uniref_abundance)
plaque_ko_abundance <- plaque_ko_abundance[, -1] # removed ungrouped and 
subset_batch_numbers <- batch_numbers[rownames(plaque_taxa_count)]
# retain samples whose library size bigger than 5e5
libsizes <- rowSums(plaque_taxa_count) 
subset_plaque_taxa_count <- plaque_taxa_count[libsizes > 5e5, ]
subset_plaque_ko_abundance <- plaque_ko_abundance[libsizes > 5e5, ]
subset_plaque_uniref_abundance <- plaque_uniref_abundance[libsizes > 5e5, ]
subset_batch_numbers <- subset_batch_numbers[libsizes > 5e5]
prevalences <- colMeans(subset_plaque_taxa_count > 0)
# retain taxa whose taxa prevalence bigger than 0.05
subset_plaque_taxa_count <- subset_plaque_taxa_count[, prevalences > 0.1] 
prevalences <- colMeans(subset_plaque_ko_abundance > 0)
subset_plaque_ko_abundance <- subset_plaque_ko_abundance[, prevalences > 0.1]
prevalences <- colMeans(subset_plaque_uniref_abundance > 0)
subset_plaque_uniref_abundance <- subset_plaque_uniref_abundance[, prevalences > 0.1]


# retrieve metadata
metadata <- read.table("metadata/metadata_yr1_imputed.tsv",
                       sep='\t', header=1)
rownames(metadata) <- as.character(metadata$BabySubjectID)
indvonly <- function(sampleID) {
  strsplit(sampleID, split="-")[[1]][1]
}
indvs <- sapply(rownames(subset_plaque_taxa_count), indvonly)
metadata_allplaque <- metadata[indvs, c("Cigarettes", "region")]
# metadata_allplaque$Education_HS <- factor(metadata_allplaque$Education_HS)
# metadata_allplaque$HouseholdIncome_cat2 <- factor(metadata_allplaque$HouseholdIncome_cat2)
metadata_allplaque$Cigarettes <- factor(metadata_allplaque$Cigarettes)
metadata_allplaque$region <- factor(metadata_allplaque$region)




## batch correct plaque taxa counts
start.time <- Sys.time()
tune_plaque_taxa_count <- Tune_ConQuR(tax_tab=subset_plaque_taxa_count,
                                      batchid=subset_batch_numbers,
                                      covariates = metadata_allplaque,
                                      batch_ref_pool = c("3", "5"),
                                      logistic_lasso_pool = T,
                                      quantile_type_pool = c("standard", "lasso"),
                                      simple_match_pool = F,
                                      lambda_quantile_pool = "2p/n",
                                      interplt_pool=T,
                                      frequencyL = 0.1,
                                      frequencyU = 1,
                                      taus=seq(0.01, 0.99, by=0.01),
                                      num_core=12)
end.time <- Sys.time()

tune_plaque_taxa_count$method_final
subset_plaque_taxa_count_corrected <- tune_plaque_taxa_count$tax_final
subset_plaque_taxa_count <- data.frame(subset_plaque_taxa_count)
subset_plaque_taxa_count_corrected <- data.frame(subset_plaque_taxa_count_corrected)

write.table(subset_plaque_taxa_count, "subset_plaque_taxa_count.tsv", sep="\t",
            quote=FALSE)
write.table(subset_plaque_taxa_count_corrected, "subset_plaque_taxa_count_corrected.tsv", sep="\t",
            quote=FALSE)

# batch correction on plaque KEGG counts
start.time <- Sys.time()
tune_plaque_ko_abundance <- Tune_ConQuR(tax_tab=subset_plaque_ko_abundance,
                                        batchid=subset_batch_numbers,
                                        covariates = metadata_allplaque,
                                        batch_ref_pool = c("3", "5"),
                                        logistic_lasso_pool = T,
                                        quantile_type_pool = c("standard", "lasso"),
                                        simple_match_pool = F,
                                        lambda_quantile_pool = "2p/n",
                                        interplt_pool=T,
                                        frequencyL = 0.1,
                                        frequencyU = 1,
                                        taus=seq(0.01, 0.99, by=0.01),
                                        num_core=12)
end.time <- Sys.time()
time.taken <- end.time - start.time

# subset_plaque_ko_abundance <- data.frame(subset_plaque_ko_abundance)
subset_plaque_ko_abundance_corrected <- tune_plaque_ko_abundance$tax_final
# subset_plaque_ko_abundance_corrected <- data.frame(subset_plaque_ko_abundance_corrected)

write.table(subset_plaque_ko_abundance, "subset_plaque_ko_abundance.tsv", sep="\t",
            quote=FALSE)
write.table(subset_plaque_ko_abundance_corrected, "subset_plaque_ko_abundance_corrected.tsv", sep="\t",
            quote=FALSE)





# batch correction on plaque uniref counts
start.time <- Sys.time()
tune_plaque_uniref_abundance <- Tune_ConQuR(tax_tab=subset_plaque_uniref_abundance,
                                        batchid=subset_batch_numbers,
                                        covariates = metadata_allplaque,
                                        batch_ref_pool = c("3", "5"),
                                        logistic_lasso_pool = T,
                                        quantile_type_pool = c("standard", "lasso"),
                                        simple_match_pool = F,
                                        lambda_quantile_pool = "2p/n",
                                        interplt_pool=T,
                                        frequencyL = 0.1,
                                        frequencyU = 1,
                                        taus=seq(0.01, 0.99, by=0.01),
                                        num_core=12)
end.time <- Sys.time()
time.taken <- end.time - start.time

# subset_plaque_ko_abundance <- data.frame(subset_plaque_ko_abundance)
subset_plaque_uniref_abundance_corrected <- tune_plaque_uniref_abundance$tax_final
# subset_plaque_ko_abundance_corrected <- data.frame(subset_plaque_ko_abundance_corrected)

write.table(subset_plaque_uniref_abundance, "subset_plaque_uniref_abundance.tsv", sep="\t",
            quote=FALSE)
write.table(subset_plaque_uniref_abundance_corrected, "subset_plaque_uniref_abundance_corrected.tsv", sep="\t",
            quote=FALSE)





