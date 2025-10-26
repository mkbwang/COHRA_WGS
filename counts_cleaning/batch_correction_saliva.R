rm(list=ls())
library(dplyr)
library(openxlsx)
library(ConQuR)
library(sva)
library(doParallel)


clean_ID <- function(sampleid){
  all_components <- strsplit(sampleid, split="-")[[1]]
  cleaned_ID <- paste(all_components[c(1,2)], collapse="-")
  cleaned_ID
}


saliva_batchinfo <- read.xlsx("metadata/Saliva_batchinfo.xlsx", sheet="QC and metadata")
synonym_IDs <- read.csv("metadata/synonym_IDs.csv")
saliva_batchinfo$findmatch <- 0


# saliva
for (ID in seq(1, length(synonym_IDs$AltSubjectID))){
  
  original_ID <- synonym_IDs$AltSubjectID[ID]
  replace_ID <- as.character(synonym_IDs$BabySubjectID[ID])
  exist_in_saliva <- grepl(original_ID, saliva_batchinfo$Sample_ID) + 
    grepl(replace_ID, saliva_batchinfo$Sample_ID)
  saliva_batchinfo$findmatch <- saliva_batchinfo$findmatch + exist_in_saliva
  
  if (sum(exist_in_saliva) > 0){ # the ID being searched exist
    replace_loc <- which(exist_in_saliva > 0)
    saliva_batchinfo$Sample_ID[replace_loc] <- gsub(original_ID, replace_ID, 
                                                    saliva_batchinfo$Sample_ID[replace_loc])
  } 
  
}
saliva_batchinfo <- saliva_batchinfo %>% filter(findmatch > 0) %>%
  select(Sample_ID, Sample_type)
saliva_batchinfo$Sample_ID <- sapply(saliva_batchinfo$Sample_ID, clean_ID)



saliva_taxa_count <- read.csv("counts_cleaning/saliva_taxa_count_subset.csv", row.names=1)
saliva_ko_abundance <- read.csv("counts_cleaning/saliva_ko_abundance_subset.csv", row.names=1)



colnames(saliva_taxa_count) <- gsub("[.]", "-", colnames(saliva_taxa_count))
colnames(saliva_taxa_count) <- gsub("X", "", colnames(saliva_taxa_count))
colnames(saliva_ko_abundance) <- gsub("[.]", "-", colnames(saliva_ko_abundance))
colnames(saliva_ko_abundance) <- gsub("X", "", colnames(saliva_ko_abundance))



batch_numbers <- factor(saliva_batchinfo$Sample_type)
names(batch_numbers) <- saliva_batchinfo$Sample_ID
saliva_taxa_count <- t(saliva_taxa_count)
saliva_ko_abundance <- t(saliva_ko_abundance)
subset_batch_numbers <- batch_numbers[rownames(saliva_taxa_count)]
subset_batch_numbers <- as.integer(subset_batch_numbers)

# normalization
normalize <- function(mymat, libsize=1e5){
  output_mat <- mymat
  for (j in 1:nrow(mymat)){
    output_mat[j, ] <- mymat[j, ] / sum(mymat[j, ]) * libsize
  }
  return(output_mat)
}


med_libsize <- median(rowSums(saliva_taxa_count))
saliva_taxa_count_normalized <- normalize(saliva_taxa_count, libsize=med_libsize)
ps_taxa <- min(saliva_taxa_count_normalized[saliva_taxa_count_normalized > 0])
saliva_taxa_count_normalized_log <- log(saliva_taxa_count_normalized + ps_taxa/2)

med_libsize <- median(rowSums(saliva_ko_abundance))
saliva_ko_abundance_normalized <- normalize(saliva_ko_abundance, libsize=med_libsize)
ps_ko <- min(saliva_ko_abundance_normalized[saliva_ko_abundance_normalized > 0])
saliva_ko_abundance_normalized_log <- log(saliva_ko_abundance_normalized + ps_ko/2)



metadata <- read.table("metadata/metadata_yr1_imputed.tsv",
                       sep='\t', header=1)
rownames(metadata) <- as.character(metadata$BabySubjectID)
indvonly <- function(sampleID) {
  strsplit(sampleID, split="-")[[1]][1]
}
indvs <- sapply(rownames(saliva_taxa_count), indvonly)
metadata_allsaliva <- metadata[indvs, ]
casestatus <- metadata_allsaliva$Case_status
covar_matrix <- cbind(metadata_allsaliva$Cigarettes == "Yes", 
                      metadata_allsaliva$region == "WV") |> as.matrix()


## batch correct saliva taxa counts

saliva_taxa_count_log_corrected <- ComBat(dat=t(saliva_taxa_count_normalized_log),
                                          batch=subset_batch_numbers,
                                          mod=covar_matrix,
                                          par.prior=FALSE)
saliva_taxa_count_corrected <- exp(saliva_taxa_count_log_corrected) - ps_taxa/2
saliva_taxa_count_corrected[saliva_taxa_count_corrected < ps_taxa] <- 0

saliva_taxa_count_corrected <- data.frame(t(saliva_taxa_count_corrected))
write.table(saliva_taxa_count_corrected, "counts_cleaning/saliva_taxa_count_subset_corrected.tsv", sep="\t",
            quote=FALSE)



## batch correct saliva KO abundance
saliva_ko_abundance_log_corrected <- ComBat(dat=t(saliva_ko_abundance_normalized_log),
                                        batch=subset_batch_numbers,
                                        mod=covar_matrix,
                                        par.prior=FALSE)
saliva_ko_abundance_corrected <- exp(saliva_ko_abundance_log_corrected) - ps_ko/2
saliva_ko_abundance_corrected[saliva_ko_abundance_corrected < ps_ko] <- 0

saliva_ko_abundance_corrected <- as.data.frame(t(saliva_ko_abundance_corrected))
write.table(saliva_ko_abundance_corrected, "counts_cleaning/saliva_ko_abundance_subset_corrected.tsv", sep="\t",
            quote=FALSE)




## batch correct saliva uniref abundance
# saliva_uniref_abundance_log_corrected <- ComBat(dat=t(saliva_uniref_abundance_normalized_log),
#                                             batch=subset_batch_numbers,
#                                             mod=covar_matrix,
#                                             par.prior=FALSE)
# saliva_uniref_abundance_corrected <- exp(saliva_uniref_abundance_log_corrected) - ps_uniref/2
# saliva_uniref_abundance_corrected[saliva_uniref_abundance_corrected < 0] <- 0

# start.time <- Sys.time()
# tune_saliva_uniref_abundance <- Tune_ConQuR(tax_tab=saliva_uniref_abundance,
#                                         batchid=subset_batch_numbers,
#                                         covariates = metadata_allsaliva,
#                                         batch_ref_pool = "human saliva",
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

# saliva_uniref_abundance_corrected <- as.data.frame(t(saliva_uniref_abundance_corrected))
# write.table(saliva_uniref_abundance_corrected, "saliva_uniref_abundance_subset_corrected.tsv", sep="\t",
#             quote=FALSE)


