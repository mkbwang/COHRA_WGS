
rm(list=ls())
library(mice)
metadata_yr1 <- read.csv("metadata/metadata_yr1.csv")
# add region columns
region <- metadata_yr1$BabySubjectID %/% 1e6
metadata_yr1$region <- "Pitt"
metadata_yr1$region[metadata_yr1$BabySubjectID %/% 1e6 > 11] <- "WV"

# change breastfed variable to binary
metadata_yr1$Breastfed <- metadata_yr1$Breastfed == 1
metadata_yr1$BFCurrent <- NULL
metadata_yr1$BFStopMnth <- NULL


# missing_count_rows <- rowSums(is.na(metadata_yr1))

# remove four individuals with a large number of missing data
# indv_missing_a_lot <- which(missing_count_rows > 20)
# metadata_yr1 <- metadata_yr1[-indv_missing_a_lot, ]

missing_count_columns <- colSums(is.na(metadata_yr1))
# remove EatPizza column: too many missing
metadata_yr1$EatPizza <- NULL
# D1FT scores are zero for all infants, remove this column
metadata_yr1$Prim_d1ft <- NULL



# impute missing data of household income, delivery method and breastfeed information
## first change set up all the factor variables

missing_impute <- function(metadata){
  
  metadata$HouseholdIncome_cat2 <- factor(metadata$HouseholdIncome_cat2, 
                                              levels = c("<$25,000", "$25,000-$74,999", ">=$75,000"))
  metadata$Education_HS <- factor(metadata$Education_HS,
                                      levels=c("High school degree or less", "Associates degree or higher")) 
  metadata$Cigarettes <- factor(metadata$Cigarettes,
                                    levels=c("No", "Yes"))
  metadata$Delivery <- factor(metadata$Delivery,
                                  levels=c("Vaginal", "C-section"))
  metadata$BabySex <- factor(metadata$BabySex,
                                 levels=c("Male", "Female"))
  metadata$region <- factor(metadata$region, levels=c("Pitt", "WV"))
  
  # use mice to impute the metadata
  column_missingness <- colSums(is.na(metadata))
  ## select the major metadata for imputations
  meta_major <- metadata[, c("MotherAgeAtExam", "HouseholdIncome_cat2", 
                                 "Education_HS", "PERM_D2MFT", "Cigarettes", "Delivery", "Prim_Tot_Teeth_Present",
                                 "Breastfed", "region")]
  
  
  imputed_meta_major <- mice(meta_major) |> complete()
  metadata[, c("MotherAgeAtExam", "HouseholdIncome_cat2", 
                   "Education_HS", "PERM_D2MFT", "Cigarettes", "Delivery", "Prim_Tot_Teeth_Present",
                   "Breastfed", "region")] <- imputed_meta_major
  
  return(metadata)
  
}
set.seed(2024)
metadata_yr1_imputed <- missing_impute(metadata_yr1)

# use mice to impute the metadata
column_missingness <- colSums(is.na(metadata_yr1_imputed))

write.table(metadata_yr1_imputed, file.path("metadata/metadata_yr1_imputed.tsv"),
            sep='\t', row.names=F, quote=F)


# impute the missing values in metadata at yr 2 as well
metadata_yr2 <- read.csv("metadata/metadata_yr2.csv")
metadata_yr2$region <- "Pitt"
metadata_yr2$region[metadata_yr2$BabySubjectID %/% 1e6 > 11] <- "WV"
metadata_yr2 <- metadata_yr2[, colnames(metadata_yr1)]
subject_incomemissing <- which(is.na(metadata_yr2$HouseholdIncome_cat2))
for (id in subject_incomemissing){
  metadata_yr2$HouseholdIncome_cat2[id] <- 
    metadata_yr1_imputed$HouseholdIncome_cat2[metadata_yr1_imputed$MotherSubjectID == metadata_yr2$MotherSubjectID[id]]
}
subject_cigarettemissing <- which(is.na(metadata_yr2$Cigarettes))
for (id in subject_cigarettemissing){
  metadata_yr2$Cigarettes[id] <- 
    metadata_yr1_imputed$Cigarettes[metadata_yr1_imputed$MotherSubjectID == metadata_yr2$MotherSubjectID[id]]
}
subject_deliverymissing <- which(is.na(metadata_yr2$Delivery))
for (id in subject_deliverymissing){
  metadata_yr2$Delivery[id] <- 
    metadata_yr1_imputed$Delivery[metadata_yr1_imputed$MotherSubjectID == metadata_yr2$MotherSubjectID[id]]
}
## breast feeding has stopped at 24 months old for all kids
metadata_yr2$Breastfed <- metadata_yr2$Breastfed == 1
set.seed(2024)
metadata_yr2_imputed <- missing_impute(metadata_yr2)
metadata_yr2_imputed$Breastfed <- NULL
write.table(metadata_yr2_imputed, file.path("metadata/metadata_yr2_imputed.tsv"),
            sep='\t', row.names=F, quote=F)


# calculate statistics for metadata
library(dplyr)
# mother's age
mom_age <- metadata_yr1_imputed %>% group_by(Case_status) %>%
  summarise(mean_age = mean(MotherAgeAtExam), sd_age = sd(MotherAgeAtExam))
lm_result <- lm(MotherAgeAtExam ~ Case_status, data=metadata_yr1)


# region
region_summary <- metadata_yr1_imputed %>% group_by(Case_status, region) %>%
  summarise(count=n())

count_matrix <- matrix(region_summary$count, nrow=2)
chisq.test(count_matrix)

# education
education_summary <- metadata_yr1_imputed %>% group_by(Case_status, Education_HS) %>%
  summarise(count=n())
count_matrix <- matrix(education_summary$count, nrow=2)
chisq.test(count_matrix)

# DMFT scores of moms
DMFT_mom_summary <- metadata_yr1_imputed %>% group_by(Case_status) %>%
  summarise(mean_DMFT=mean(PERM_D2MFT), sd_DMFT=sd(PERM_D2MFT))
dmft_model <- lm(PERM_D2MFT ~ Case_status,
                  data=metadata_yr1)


# household income
income_summary <- metadata_yr1_imputed %>% group_by(Case_status, HouseholdIncome_cat2) %>%
  summarise(count=n())
count_matrix <- matrix(income_summary$count, ncol=2)
chisq.test(count_matrix)

# smoking
smoking_summary <- metadata_yr1_imputed %>% group_by(Case_status, Cigarettes) %>%
  summarise(count=n())
count_matrix <- matrix(smoking_summary$count, nrow=2)
chisq.test(count_matrix)

# delivery
delivery_summary <- metadata_yr1_imputed %>% group_by(Case_status, Delivery) %>%
  summarise(count=n())
count_matrix <- matrix(delivery_summary$count, nrow=2)
chisq.test(count_matrix)


# baby sex
babysex_summary <- metadata_yr1_imputed %>% group_by(Case_status, BabySex) %>%
  summarise(count=n())
count_matrix <- matrix(babysex_summary$count, nrow=2)
chisq.test(count_matrix)


# primary teeth count
primtooth_summary <- metadata_yr1_imputed %>% group_by(Case_status) %>%
  summarise(mean_tooth = mean(Prim_Tot_Teeth_Present),
            sd_tooth = sd(Prim_Tot_Teeth_Present))
primtooth_model <- lm(Prim_Tot_Teeth_Present ~ Case_status,
                 data=metadata_yr1)

# breastfeeding
breastfed_summary <- metadata_yr1_imputed %>% group_by(Case_status, Breastfed) %>%
  summarise(count=n())
count_matrix <- matrix(breastfed_summary$count, nrow=2)
chisq.test(count_matrix)

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

