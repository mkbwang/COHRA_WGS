rm(list=ls())
library(dplyr)
library(openxlsx)
library(ConQuR)
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



saliva_taxa_count <- read.table("saliva_preprocessing/metaphlan_output/joint_taxonomic_counts.tsv",
                                sep='\t')
saliva_ko_abundance <- read.table("saliva_preprocessing/humann_output/joint_ko_table_concise.tsv",
                                  sep='\t')
colnames(saliva_taxa_count) <- gsub("[.]", "-", colnames(saliva_taxa_count))
colnames(saliva_taxa_count) <- gsub("X", "", colnames(saliva_taxa_count))
colnames(saliva_ko_abundance) <- gsub("[.]", "-", colnames(saliva_ko_abundance))
colnames(saliva_ko_abundance) <- gsub("X", "", colnames(saliva_ko_abundance))



batch_numbers <- factor(saliva_batchinfo$Sample_type)
names(batch_numbers) <- saliva_batchinfo$Sample_ID
saliva_taxa_count <- t(saliva_taxa_count)
saliva_ko_abundance <- t(saliva_ko_abundance)
saliva_ko_abundance <- saliva_ko_abundance[, -1] # removed ungrouped and 
subset_batch_numbers <- batch_numbers[rownames(saliva_taxa_count)]
# retain samples whose library size bigger than 5e5
libsizes <- rowSums(saliva_taxa_count) 
subset_saliva_taxa_count <- saliva_taxa_count[libsizes > 5e5, ]
subset_saliva_ko_abundance <- saliva_ko_abundance[libsizes > 5e5, ]
subset_batch_numbers <- subset_batch_numbers[libsizes > 5e5]
prevalences <- colMeans(subset_saliva_taxa_count > 0)
# retain taxa whose taxa prevalence bigger than 0.1
subset_saliva_taxa_count <- subset_saliva_taxa_count[, prevalences > 0.1] 
prevalences <- colMeans(subset_saliva_ko_abundance > 0)
subset_saliva_ko_abundance <- subset_saliva_ko_abundance[, prevalences > 0.1]


metadata <- read.table("metadata/metadata_yr1_imputed.tsv",
                       sep='\t', header=1)
rownames(metadata) <- as.character(metadata$BabySubjectID)
indvonly <- function(sampleID) {
  strsplit(sampleID, split="-")[[1]][1]
}
indvs <- sapply(rownames(subset_saliva_taxa_count), indvonly)
metadata_allsaliva <- metadata[indvs, c("Cigarettes", "region")]
#metadata_allsaliva$Education_HS <- factor(metadata_allsaliva$Education_HS)
#metadata_allsaliva$HouseholdIncome_cat2 <- factor(metadata_allsaliva$HouseholdIncome_cat2)
metadata_allsaliva$Cigarettes <- factor(metadata_allsaliva$Cigarettes)
metadata_allsaliva$region <- factor(metadata_allsaliva$region)


## batch correct saliva taxa counts
start.time <- Sys.time()
tune_saliva_taxa_count <- Tune_ConQuR(tax_tab=subset_saliva_taxa_count,
                                      batchid=subset_batch_numbers,
                                      covariates = metadata_allsaliva,
                                      batch_ref_pool = "human saliva",
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


subset_saliva_taxa_count_corrected <- tune_saliva_taxa_count$tax_final


subset_saliva_taxa_count <- data.frame(subset_saliva_taxa_count)
subset_saliva_taxa_count_corrected <- data.frame(subset_saliva_taxa_count_corrected)


write.table(subset_saliva_taxa_count, "counts_cleaning/subset_saliva_taxa_count.tsv", sep="\t",
            quote=FALSE)
write.table(subset_saliva_taxa_count_corrected, "counts_cleaning/subset_saliva_taxa_count_corrected.tsv", sep="\t",
            quote=FALSE)


## batch correct saliva KO abundance
start.time <- Sys.time()
tune_saliva_ko_abundance <- Tune_ConQuR(tax_tab=subset_saliva_ko_abundance,
                                        batchid=subset_batch_numbers,
                                        covariates = metadata_allsaliva,
                                        batch_ref_pool = "human saliva",
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

subset_saliva_ko_abundance_corrected <- tune_saliva_ko_abundance$tax_final



write.table(subset_saliva_ko_abundance, "counts_cleaning/subset_saliva_ko_abundance.tsv", sep="\t",
            quote=FALSE)
write.table(subset_saliva_ko_abundance_corrected, "counts_cleaning/subset_saliva_ko_abundance_corrected.tsv", sep="\t",
            quote=FALSE)
