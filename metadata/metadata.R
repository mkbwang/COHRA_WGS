# clean up metadata
library(openxlsx)
library(dplyr)

# load basic metadata (babyID, motherID and eventual case status)
kids_diagnosis  <- read.xlsx("metadata/COHRA_2_metagenomic_samples_list.xlsx", sheet="Project_2")
kids_diagnosis <- kids_diagnosis %>% filter(Visit==5) %>% select(BabysubjectID, Case_status, Yr_2) %>%
  rename(BabySubjectID = BabysubjectID)

kids_diagnosis$MotherSubjectID <- kids_diagnosis$BabySubjectID - 1

metadata_yr1 <- kids_diagnosis %>% select(MotherSubjectID, BabySubjectID, Case_status)
metadata_yr2 <- kids_diagnosis %>% filter(Yr_2 == 1) %>% select(MotherSubjectID, BabySubjectID, Case_status)



# identify all the alternative IDs
alternative_IDs <- read.xlsx("metadata/COHRA_2_Readable.xlsx")
alternative_IDs <- na.omit(alternative_IDs)
subset_name <- function(myname, num_items=2){
  selected_items <- strsplit(myname, split="-")[[1]][1:num_items]
  subset_ID <- paste(selected_items, collapse="-")
  return(subset_ID)
}

alternative_IDs$readable <- sapply(alternative_IDs$readable, subset_name, num_items=2)
alternative_IDs$ID <- sapply(alternative_IDs$ID, subset_name, num_items=1)

synonym_IDs <- alternative_IDs %>% select(ID, readable) %>% unique()
colnames(synonym_IDs) <- c("BabySubjectID", "AltSubjectID")
synonym_IDs$BabySubjectID <- as.integer(synonym_IDs$BabySubjectID)

subset_synonym <- metadata_yr1 %>% left_join(synonym_IDs, by="BabySubjectID")

write.csv(subset_synonym, "metadata/synonym_IDs.csv", row.names=FALSE)




# load extra metadata
load("metadata/meta_sampling.Rdata")

survey_cols <- colnames(kid_cross12)
usefulcols <- c("BabySubjectID", "MotherAgeAtExam", "HouseholdIncome_cat2", "Education_HS", "PERM_D2MFT",
                "Cigarettes", "Delivery", "BabySex", "Prim_d1ft", "Prim_Tot_Teeth_Present",
                "EatCereals", "EatCrackers", "EatFruits", "EatPotatoes", "EatOthVeg", "EatCheese",
                "EatMeat", "EatPizza", "EatChips", "EatDesserts", "EatCandies", "CowAnimalMilk", "PowderMix",
                "PlantMilk", "FlavWater", "SportsDrink", "Juice", "Soda", "MealDrink")

selected_metadata_yr1 <- kid_cross12_sample %>% select(usefulcols)

selected_metadata_yr2 <- kid_cross24_sample %>% select(usefulcols)


othermetadata <- read.csv("metadata/AdditionalCOHRAMetaData.csv")
othermetadata_yr1 <- othermetadata %>% filter(PhCall == 5) %>% select(BabysubjectID, Breastfed, 
                                                                      BFCurrent, BFStopMnth)
colnames(othermetadata_yr1)[1] <- "BabySubjectID"
othermetadata_yr1[othermetadata_yr1 < 0] <- NA

othermetadata_yr2 <- othermetadata %>% filter(PhCall == 7) %>% select(BabysubjectID, Breastfed, 
                                                                      BFCurrent, BFStopMnth)
colnames(othermetadata_yr2)[1] <- "BabySubjectID"
othermetadata_yr2[othermetadata_yr2 < 0] <- NA


metadata_yr1 <- metadata_yr1 %>% left_join(selected_metadata_yr1, by="BabySubjectID") %>%
  left_join(othermetadata_yr1, by="BabySubjectID")

metadata_yr2 <- metadata_yr2 %>% left_join(selected_metadata_yr2, by="BabySubjectID") %>%
  left_join(othermetadata_yr2, by="BabySubjectID")


write.csv(metadata_yr1, "metadata/metadata_yr1.csv", row.names=FALSE)
write.csv(metadata_yr2, "metadata/metadata_yr2.csv", row.names=FALSE)


metadata_yr1 <- read.csv("metadata/metadata_yr1.csv")
region <- metadata_yr1$BabySubjectID %/% 1e6
table(region)

# age
hist(metadata_yr1$MotherAgeAtExam, nclass=20, xlab="Mother Age", main="Mother Age")


# delivery
table(metadata_yr1$Delivery)


# baby sex
table(metadata_yr1$BabySex)

# primary teeth count
hist(metadata_yr1$Prim_Tot_Teeth_Present, nclass=20, main="Primary teeth count", xlab="")


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

