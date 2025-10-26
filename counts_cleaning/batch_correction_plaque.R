rm(list=ls())
library(dplyr)
library(openxlsx)
library(ConQuR)
library(sva)
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



plaque_taxa_count <- read.csv("counts_cleaning/plaque_taxa_count_subset.csv", row.names=1)
colnames(plaque_taxa_count) <- gsub("[.]", "-", colnames(plaque_taxa_count))
colnames(plaque_taxa_count) <- gsub("X", "", colnames(plaque_taxa_count))
plaque_ko_abundance <- read.csv("counts_cleaning/plaque_ko_abundance_subset.csv", row.names=1)
colnames(plaque_ko_abundance) <- gsub("[.]", "-", colnames(plaque_ko_abundance))
colnames(plaque_ko_abundance) <- gsub("X", "", colnames(plaque_ko_abundance))



## batch correction
## plaque
batch_numbers <- factor(plaque_batchinfo$Batch)
names(batch_numbers) <- plaque_batchinfo$Sample
plaque_taxa_count <- t(plaque_taxa_count)
plaque_ko_abundance <- t(plaque_ko_abundance)

subset_batch_numbers <- batch_numbers[rownames(plaque_taxa_count)]
subset_batch_numbers <- as.integer(subset_batch_numbers)

# normalization
normalize <- function(mymat, libsize=1e5){
  output_mat <- mymat
  for (j in 1:nrow(mymat)){
    output_mat[j, ] <- mymat[j, ] / sum(mymat[j, ]) * libsize
  }
  return(output_mat)
}

med_libsize <- median(rowSums(plaque_taxa_count))
plaque_taxa_count_normalized <- normalize(plaque_taxa_count, libsize=med_libsize)
ps_taxa <- min(plaque_taxa_count_normalized[plaque_taxa_count_normalized > 0])
plaque_taxa_count_normalized_log <- log(plaque_taxa_count_normalized + ps_taxa/2)

med_libsize <- median(rowSums(plaque_ko_abundance))
plaque_ko_abundance_normalized <- normalize(plaque_ko_abundance, libsize=med_libsize)
ps_ko <- min(plaque_ko_abundance_normalized[plaque_ko_abundance_normalized > 0])
plaque_ko_abundance_normalized_log <- log(plaque_ko_abundance_normalized + ps_ko/2)


# retrieve metadata
metadata <- read.table("metadata/metadata_yr1_imputed.tsv",
                       sep='\t', header=1)
rownames(metadata) <- as.character(metadata$BabySubjectID)
indvonly <- function(sampleID) {
  strsplit(sampleID, split="-")[[1]][1]
}
indvs <- sapply(rownames(plaque_taxa_count), indvonly)
metadata_allplaque <- metadata[indvs, c("Cigarettes", "region")]
covar_matrix <- cbind(metadata_allplaque$Cigarettes == "Yes", 
                      metadata_allplaque$region == "WV") |> as.matrix()




## batch correct plaque taxa counts
plaque_taxa_count_log_corrected <- ComBat(dat=t(plaque_taxa_count_normalized_log),
                                          batch=subset_batch_numbers,
                                          mod=covar_matrix,
                                          par.prior=FALSE)
plaque_taxa_count_corrected <- exp(plaque_taxa_count_log_corrected) - ps_taxa/2
plaque_taxa_count_corrected[plaque_taxa_count_corrected < ps_taxa] <- 0

plaque_taxa_count_corrected <- data.frame(t(plaque_taxa_count_corrected))
write.table(plaque_taxa_count_corrected, "counts_cleaning/plaque_taxa_count_subset_corrected.tsv", sep="\t",
            quote=FALSE)



## batch correct plaque KO abundance
plaque_ko_abundance_log_corrected <- ComBat(dat=t(plaque_ko_abundance_normalized_log),
                                            batch=subset_batch_numbers,
                                            mod=covar_matrix,
                                            par.prior=FALSE)
plaque_ko_abundance_corrected <- exp(plaque_ko_abundance_log_corrected) - ps_ko/2
plaque_ko_abundance_corrected[plaque_ko_abundance_corrected < ps_ko] <- 0

plaque_ko_abundance_corrected <- as.data.frame(t(plaque_ko_abundance_corrected))
write.table(plaque_ko_abundance_corrected, "counts_cleaning/plaque_ko_abundance_subset_corrected.tsv", sep="\t",
            quote=FALSE)





# batch correction on plaque uniref counts
# start.time <- Sys.time()
# tune_plaque_uniref_abundance <- Tune_ConQuR(tax_tab=subset_plaque_uniref_abundance,
#                                         batchid=subset_batch_numbers,
#                                         covariates = metadata_allplaque,
#                                         batch_ref_pool = c("3", "5"),
#                                         logistic_lasso_pool = T,
#                                         quantile_type_pool = c("standard", "lasso"),
#                                         simple_match_pool = F,
#                                         lambda_quantile_pool = "2p/n",
#                                         interplt_pool=T,
#                                         frequencyL = 0.1,
#                                         frequencyU = 1,
#                                         taus=seq(0.01, 0.99, by=0.01),
#                                         num_core=12)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# 
# # subset_plaque_ko_abundance <- data.frame(subset_plaque_ko_abundance)
# subset_plaque_uniref_abundance_corrected <- tune_plaque_uniref_abundance$tax_final
# # subset_plaque_ko_abundance_corrected <- data.frame(subset_plaque_ko_abundance_corrected)
# 
# write.table(subset_plaque_uniref_abundance, "subset_plaque_uniref_abundance.tsv", sep="\t",
#             quote=FALSE)
# write.table(subset_plaque_uniref_abundance_corrected, "subset_plaque_uniref_abundance_corrected.tsv", sep="\t",
#             quote=FALSE)
# 
# 
# 
# 
# 
