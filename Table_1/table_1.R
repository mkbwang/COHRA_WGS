

metadata_yr1 <- read.csv("metadata/metadata_yr1.csv")
# add region columns
region <- metadata_yr1$BabySubjectID %/% 1e6
metadata_yr1$region <- "Pitt"
metadata_yr1$region[metadata_yr1$BabySubjectID %/% 1e6 > 11] <- "WV"

# change breastfed variable to binary
metadata_yr1$Breastfed <- metadata_yr1$Breastfed == 1
metadata_yr1$BFCurrent <- NULL
metadata_yr1$BFStopMnth <- NULL


missing_count_rows <- rowSums(is.na(metadata_yr1))

# remove four individuals with a large number of missing data
indv_missing_a_lot <- which(missing_count_rows > 20)
metadata_yr1 <- metadata_yr1[-indv_missing_a_lot, ]

missing_count_columns <- colSums(is.na(metadata_yr1))
# remove EatPizza column: too many missing
metadata_yr1$EatPizza <- NULL
# D1FT scores are zero for all infants, remove this column
metadata_yr1$Prim_d1ft <- NULL



# impute missing data of household income, delivery method and breastfeed information
## first change set up all the factor variables
metadata_yr1$HouseholdIncome_cat2 <- factor(metadata_yr1$HouseholdIncome_cat2, 
                                               levels = c("<$25,000", "$25,000-$74,999", ">=$75,000"))
metadata_yr1$Education_HS <- factor(metadata_yr1$Education_HS,
                                    levels=c("High school degree or less", "Associates degree or higher")) 
metadata_yr1$Cigarettes <- factor(metadata_yr1$Cigarettes,
                                  levels=c("No", "Yes"))
metadata_yr1$Delivery <- factor(metadata_yr1$Delivery,
                                levels=c("Vaginal", "C-section"))
metadata_yr1$BabySex <- factor(metadata_yr1$BabySex,
                               levels=c("Male", "Female"))
metadata_yr1$region <- factor(metadata_yr1$region, levels=c("Pitt", "WV"))


# use mice to impute the metadata
library(mice)
column_missingness <- colSums(is.na(metadata_yr1))
## select the major metadata for imputations
meta_major <- metadata_yr1[, c("MotherAgeAtExam", "HouseholdIncome_cat2", 
                                  "Education_HS", "PERM_D2MFT", "Cigarettes", "Delivery", "Prim_Tot_Teeth_Present",
                                  "Breastfed", "region")]

imputed_meta_major <- mice(meta_major) |> complete()
metadata_yr1[, c("MotherAgeAtExam", "HouseholdIncome_cat2", 
                 "Education_HS", "PERM_D2MFT", "Cigarettes", "Delivery", "Prim_Tot_Teeth_Present",
                 "Breastfed", "region")] <- imputed_meta_major

write.table(metadata_yr1, file.path("metadata/metadata_yr1_imputed.tsv"),
            sep='\t', row.names=F, quote=F)

# calculate statistics for metadata
library(dplyr)
# mother's age
mom_age <- metadata_yr1 %>% group_by(Case_status) %>%
  summarise(mean_age = mean(MotherAgeAtExam), sd_age = sd(MotherAgeAtExam))
lm_result <- lm(MotherAgeAtExam ~ Case_status, data=metadata_yr1)


# region
region_summary <- metadata_yr1 %>% group_by(Case_status, region) %>%
  summarise(count=n())

count_matrix <- matrix(region_summary$count, nrow=2)
chisq.test(count_matrix)

# education
education_summary <- metadata_yr1 %>% group_by(Case_status, Education_HS) %>%
  summarise(count=n())
count_matrix <- matrix(education_summary$count, nrow=2)
chisq.test(count_matrix)

# DMFT scores of moms
DMFT_mom_summary <- metadata_yr1 %>% group_by(Case_status) %>%
  summarise(mean_DMFT=mean(PERM_D2MFT), sd_DMFT=sd(PERM_D2MFT))
dmft_model <- lm(PERM_D2MFT ~ Case_status,
                  data=metadata_yr1)


# smoking
smoking_summary <- metadata_yr1 %>% group_by(Case_status, Cigarettes) %>%
  summarise(count=n())
count_matrix <- matrix(smoking_summary$count, nrow=2)
chisq.test(count_matrix)

# delivery
delivery_summary <- metadata_yr1 %>% group_by(Case_status, Delivery) %>%
  summarise(count=n())
count_matrix <- matrix(delivery_summary$count, nrow=2)
chisq.test(count_matrix)


# baby sex
babysex_summary <- metadata_yr1 %>% group_by(Case_status, BabySex) %>%
  summarise(count=n())
count_matrix <- matrix(babysex_summary$count, nrow=2)
chisq.test(count_matrix)


# primary teeth count
primtooth_summary <- metadata_yr1 %>% group_by(Case_status) %>%
  summarise(mean_tooth = mean(Prim_Tot_Teeth_Present),
            sd_tooth = sd(Prim_Tot_Teeth_Present))
primtooth_model <- lm(Prim_Tot_Teeth_Present ~ Case_status,
                 data=metadata_yr1)


# diet
allcols <- colnames(metadata_yr1)
diet_items <- allcols[9:26]
diet_items <- diet_items[diet_items != "EatPizza"]
category_counts <- matrix(0, nrow=length(diet_items), ncol=4)
na_count <- rep(0, length(diet_items))
for (j in 1:length(diet_items)){
  selected_col <- diet_items[j]
  na_count[j] <- sum(is.na(metadata_yr1[, selected_col]))
  category_counts[j, ] <- table(metadata_yr1[, selected_col])
}

diet_summary <- as.data.frame(category_counts, row.names=diet_items)
colnames(diet_summary) <- c("Every_few_days", "Never_or_once", "Once_a_day", "Several_times_a_day")
diet_summary$NAs <- na_count
write.csv(diet_summary, "metadata/diet_survey.csv")

# breast fed
table(metadata_yr1$Breastfed)



