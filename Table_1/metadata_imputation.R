
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


