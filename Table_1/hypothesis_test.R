
rm(list=ls())
# calculate statistics for metadata
library(dplyr)

metadata_full <- read.table("metadata/metadata_yr1_imputed.tsv", sep="\t", header=1)
rownames(metadata_full) <- metadata_full$BabySubjectID

plaque_taxa_count <- read.csv("counts_cleaning/plaque_taxa_count_subset_corrected.tsv",
                              sep='\t', header=1, row.names=1) |> t() |> as.data.frame()
plaque_ids <- gsub("-5", "", colnames(plaque_taxa_count))


saliva_taxa_count <- read.csv("counts_cleaning/saliva_taxa_count_subset_corrected.tsv",
                              sep='\t', header=1, row.names=1) |> t() |> as.data.frame()
saliva_ids <- gsub("-5", "", colnames(saliva_taxa_count))


intersection_ids <- intersect(plaque_ids, saliva_ids) |> as.integer()
union_ids <- unique(c(plaque_ids, saliva_ids)) |> as.integer()


metadata_union_subset <- metadata_full %>% filter(BabySubjectID %in% union_ids)
metadata_intersection_subset <- metadata_full %>% filter(BabySubjectID %in% intersection_ids)
metadata_plaque_subset <- metadata_full %>% filter(BabySubjectID %in% plaque_ids)
metadata_saliva_subset <- metadata_full %>% filter(BabySubjectID %in% saliva_ids)



demographics_summary <- function(mymetadata){
  
  output <- list()
  # mother's age
  mom_age <- mymetadata %>% group_by(Case_status) %>%
    summarise(mean_age = mean(MotherAgeAtExam), sd_age = sd(MotherAgeAtExam))
  lm_result <- lm(MotherAgeAtExam ~ Case_status, data=mymetadata) |> summary()
  output$mom_age <- mom_age
  output$mom_age_test <- lm_result
  
  # region
  region_summary <- mymetadata %>% group_by(Case_status, region) %>%
    summarise(count=n())
  
  output$region <- region_summary
  count_matrix <- matrix(region_summary$count, nrow=2)
  output$region_test <- chisq.test(count_matrix)
  
  # education
  education_summary <- mymetadata %>% group_by(Case_status, Education_HS) %>%
    summarise(count=n())
  output$education <- education_summary
  
  count_matrix <- matrix(education_summary$count, nrow=2)
  output$education_test <- chisq.test(count_matrix)
  
  # DMFT scores of moms
  DMFT_mom_summary <- mymetadata %>% group_by(Case_status) %>%
    summarise(mean_DMFT=mean(PERM_D2MFT), sd_DMFT=sd(PERM_D2MFT))
  output$DMFT <- DMFT_mom_summary
  dmft_model <- wilcox.test(PERM_D2MFT ~ Case_status,
                   data=mymetadata)
  output$DMFT_test <- dmft_model
  
  
  # household income
  income_summary <- mymetadata %>% group_by(Case_status, HouseholdIncome_cat2) %>%
    summarise(count=n())
  output$income <- income_summary
  count_matrix <- matrix(income_summary$count, ncol=2)
  output$income_test <- chisq.test(count_matrix)
  
  # smoking
  smoking_summary <- mymetadata %>% group_by(Case_status, Cigarettes) %>%
    summarise(count=n())
  output$smoking <- smoking_summary
  count_matrix <- matrix(smoking_summary$count, nrow=2)
  output$smoking_test <- chisq.test(count_matrix)
  
  # delivery
  delivery_summary <- mymetadata %>% group_by(Case_status, Delivery) %>%
    summarise(count=n())
  output$delivery <- delivery_summary
  count_matrix <- matrix(delivery_summary$count, nrow=2)
  output$delivery_test <- chisq.test(count_matrix)
  
  
  # baby sex
  babysex_summary <- mymetadata %>% group_by(Case_status, BabySex) %>%
    summarise(count=n())
  output$babysex <- babysex_summary
  
  count_matrix <- matrix(babysex_summary$count, nrow=2)
  output$babysex_test <- chisq.test(count_matrix)
  
  
  # primary teeth count
  primtooth_summary <- mymetadata %>% group_by(Case_status) %>%
    summarise(mean_tooth = mean(Prim_Tot_Teeth_Present),
              sd_tooth = sd(Prim_Tot_Teeth_Present))
  output$primtooth <- primtooth_summary
  
  primtooth_model <- lm(Prim_Tot_Teeth_Present ~ Case_status,
                        data=mymetadata) |> summary()
  output$primtooth_test <- primtooth_model
  
  
  # breastfeeding
  breastfed_summary <- mymetadata %>% group_by(Case_status, Breastfed) %>%
    summarise(count=n())
  output$breastfed <- breastfed_summary
  
  count_matrix <- matrix(breastfed_summary$count, nrow=2)
  output$breastfed_test <- chisq.test(count_matrix)
  
  return(output)

}


union_output <- demographics_summary(metadata_union_subset)


# diet
categorical_testing <- function(covariate) {
  
  counts_summary <- metadata_yr1_imputed %>% group_by(Case_status, .data[[covariate]]) %>%
    summarise(count=n())
  unique_categories <- unique(metadata_yr1[, covariate])
  count_matrix <- counts_summary %>%
    pivot_wider(names_from = covariate, values_from = count)
  count_matrix[is.na(count_matrix)] <- 0
  test_result <- chisq.test(count_matrix[, 2:(length(unique_categories) + 1)])
  test_pval <- test_result$p.value
  return(list(variable=covariate, contingency_table=count_matrix, test_pval=test_pval))
  
}

diet_variables <- colnames(metadata_yr1_imputed)[12: 29]

diet_tests <- list()

for (j in 1:length(diet_variables)){
  
  diet_tests[[length(diet_tests) + 1]] <- categorical_testing(diet_variables[j])
  
}

