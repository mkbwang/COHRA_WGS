
rm(list=ls())
# load metadata
metadata_yr1 <- read.table(file="metadata/metadata_yr1_imputed.tsv",
                           sep='\t', header=T)
metadata_yr2 <- read.table(file="metadata/metadata_yr2_imputed.tsv",
                           sep='\t', header=T)


# load taxa counts
plaque_taxa_counts <- read.table(file="counts_cleaning/subset_plaque_taxa_count_corrected.tsv",
                                 header=T)
saliva_taxa_counts <- read.table(file="counts_cleaning/subset_saliva_taxa_count_corrected.tsv",
                                 header=T)
plaque_ko_counts <- read.table(file="counts_cleaning/subset_plaque_ko_abundance_corrected.tsv",
                               header=T)
saliva_ko_counts <- read.table(file="counts_cleaning/subset_saliva_ko_abundance_corrected.tsv",
                               header=T)



split_by_year <- function(count_table){
  
  # split the samples into 12-month and 24-month samples
  samples <- rownames(count_table)
  yr1_samples <- grepl("-5", samples)
  yr2_samples <- grepl("-7", samples)
  count_yr1 <- count_table[yr1_samples, ]
  rownames(count_yr1) <- gsub("-5", "", rownames(count_yr1))
  
  count_yr2 <- count_table[yr2_samples, ]
  rownames(count_yr2) <- gsub("-7", "", rownames(count_yr2))
  
  return(list(count_yr1=count_yr1, count_yr2=count_yr2))
  
}


# simplify the column name in the count table and separate out a taxonomy table
gen_tax_table <- function(count_table){
  
  taxonomies <- colnames(count_table)
  tax_table <- sapply(taxonomies, function(fulltaxo) strsplit(fulltaxo, split="[.]")[[1]]) |> t() |>
    as.data.frame()
  colnames(tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Group", "Species")
  species_names <- sapply(tax_table[, 7], function(name) gsub("s__", "", name))
  colnames(count_table) <- species_names
  rownames(tax_table) <- species_names
  
  return(list(count_table=count_table, taxonomy=tax_table))
  
}




## saliva taxa

saliva_taxa_counts_split <- split_by_year(saliva_taxa_counts)
saliva_taxa_yr1 <- gen_tax_table(saliva_taxa_counts_split$count_yr1)
saliva_taxa_counts_yr1 <- saliva_taxa_yr1$count_table
write.table(saliva_taxa_counts_yr1, "counts_cleaning/strata/saliva_taxa_counts_yr1.tsv",
            sep='\t', quote=F)
saliva_taxonomy <- saliva_taxa_yr1$taxonomy
write.table(saliva_taxonomy, "counts_cleaning/strata/saliva_taxonomy.tsv",
            sep='\t', quote=F)

saliva_taxa_yr2 <- gen_tax_table(saliva_taxa_counts_split$count_yr2)
saliva_taxa_counts_yr2 <- saliva_taxa_yr2$count_table
write.table(saliva_taxa_counts_yr2, "counts_cleaning/strata/saliva_taxa_counts_yr2.tsv",
            sep='\t', quote=F)


## plaque taxa

plaque_taxa_counts_split <- split_by_year(plaque_taxa_counts)
plaque_taxa_yr1 <- gen_tax_table(plaque_taxa_counts_split$count_yr1)
plaque_taxa_counts_yr1 <- plaque_taxa_yr1$count_table
write.table(plaque_taxa_counts_yr1, "counts_cleaning/strata/plaque_taxa_counts_yr1.tsv",
            sep='\t', quote=F)
plaque_taxonomy <- plaque_taxa_yr1$taxonomy
write.table(plaque_taxonomy, "counts_cleaning/strata/plaque_taxonomy.tsv",
            sep='\t', quote=F)

plaque_taxa_yr2 <- gen_tax_table(plaque_taxa_counts_split$count_yr2)
plaque_taxa_counts_yr2 <- plaque_taxa_yr2$count_table
write.table(plaque_taxa_counts_yr2, "counts_cleaning/strata/plaque_taxa_counts_yr2.tsv",
            sep='\t', quote=F)


## saliva KEGG ortholog

saliva_ko_counts_split <- split_by_year(saliva_ko_counts)
write.table(saliva_ko_counts_split$count_yr1, "counts_cleaning/strata/saliva_ko_counts_yr1.tsv",
            sep='\t', quote=F)
write.table(saliva_ko_counts_split$count_yr2, "counts_cleaning/strata/saliva_ko_counts_yr2.tsv",
            sep='\t', quote=F)


## plaque KEGG ortholog

plaque_ko_counts_split <- split_by_year(plaque_ko_counts)
write.table(plaque_ko_counts_split$count_yr1, "counts_cleaning/strata/plaque_ko_counts_yr1.tsv",
            sep='\t', quote=F)
write.table(plaque_ko_counts_split$count_yr2, "counts_cleaning/strata/plaque_ko_counts_yr2.tsv",
            sep='\t', quote=F)

## metadata at year 1


rownames(metadata_yr1) <- as.character(metadata_yr1$BabySubjectID)
metadata_yr1 <- metadata_yr1[, c("BabySubjectID", "Case_status", "HouseholdIncome_cat2",
                                 "Education_HS", "PERM_D2MFT", "Cigarettes", "Delivery", "BabySex",
                                 "Breastfed", "region")]
metadata_yr1_saliva <- metadata_yr1[rownames(saliva_taxa_counts_yr1), ]
metadata_yr1_plaque <- metadata_yr1[rownames(plaque_taxa_counts_yr1), ]
write.table(metadata_yr1_saliva, "counts_cleaning/strata/metadata_saliva_yr1.tsv",
            sep='\t', quote=F)
write.table(metadata_yr1_plaque, "counts_cleaning/strata/metadata_plaque_yr1.tsv",
            sep='\t', quote=F)

## metadata at year 2


rownames(metadata_yr2) <- as.character(metadata_yr2$BabySubjectID)
metadata_yr2_saliva <- metadata_yr2[rownames(saliva_taxa_counts_yr2), ]
metadata_yr2_plaque <- metadata_yr2[rownames(plaque_taxa_counts_yr2), ]
write.table(metadata_yr2_saliva, "counts_cleaning/strata/metadata_saliva_yr2.tsv",
            sep='\t', quote=F)
write.table(metadata_yr2_plaque, "counts_cleaning/strata/metadata_plaque_yr2.tsv",
            sep='\t', quote=F)



