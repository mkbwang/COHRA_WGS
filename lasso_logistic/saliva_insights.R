
rm(list=ls())
library(openxlsx)
# further insight into saliva prediction results (Taxa)
# load observed counts and metadata
saliva_counts <- read.table("counts_cleaning/saliva_taxa_count_subset_corrected.tsv",
                            header=TRUE, sep="\t", row.names=1) |> as.matrix()
metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
diagnoses <- (metadata_saliva_yr1$Case_status)
saliva_counts_imputed <- simple_impute(saliva_counts, scale=0.5) |> t()


# load DAA results
DAA_taxa_results <- read.csv("DAA/DAA_taxa_saliva.csv")
marker_taxa <- DAA_taxa_results$Taxa[DAA_taxa_results$pval < 0.1]


# start with DAA markers
shorten_names <- function(longname){
  strsplit(longname, split="__")[[1]] |> tail(1)
}
species_names <- sapply(colnames(saliva_counts), shorten_names) |> unname()
colnames(saliva_counts) <- colnames(saliva_counts_imputed) <- species_names

saliva_counts <- saliva_counts[, marker_taxa]
saliva_counts_imputed <- saliva_counts_imputed[, marker_taxa]

clr_saliva_counts <- clr_transform(saliva_counts_imputed)
colnames(clr_saliva_counts) <- colnames(saliva_counts)


# load the useful predictive markers
saliva_taxa_biomarkers <- read.csv("lasso_logistic/taxa/saliva_taxa_biomarkers.csv")

predictive_taxa <- saliva_taxa_biomarkers$Taxa
correlated_taxa_results <- list()
for (j in 1:length(predictive_taxa)){
  
  marker_values <- clr_saliva_counts[, predictive_taxa[j]]
  
  correlations_df <- data.frame(Taxa=colnames(clr_saliva_counts),
                                Correlation=0)
  for (k in 1:ncol(clr_saliva_counts)){
    correlations_df$Correlation[k] <- cor(marker_values, clr_saliva_counts[, k],
                                          method="spearman")
  }
  correlations_df <- correlations_df %>% arrange(-Correlation)
  correlated_taxa_results[[predictive_taxa[j]]] <- correlations_df[-1, ]
  
}


write.xlsx(correlated_taxa_results, 
           file = "lasso_logistic/taxa/saliva_correlated_marker_taxa.xlsx")


# further insight into saliva prediction results (KEGG)
# load observed counts and metadata
saliva_counts <- read.table("counts_cleaning/saliva_ko_abundance_subset_corrected.tsv",
                            header=TRUE, sep="\t", row.names=1) |> as.matrix()
metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
diagnoses <- (metadata_saliva_yr1$Case_status)
saliva_counts_imputed <- simple_impute(saliva_counts, scale=0.5) |> t()


# load DAA results
DAA_ko_results <- read.csv("DAA/DAA_ko_saliva.csv")
marker_ko <- DAA_ko_results$Taxa[DAA_ko_results$pval < 0.05]


# start with DAA markers
saliva_counts <- saliva_counts[, marker_ko]
saliva_counts_imputed <- saliva_counts_imputed[, marker_ko]

clr_saliva_counts <- clr_transform(saliva_counts)
colnames(clr_saliva_counts) <- colnames(saliva_counts)


saliva_KEGG_biomarkers <- read.csv("lasso_logistic/KEGG/saliva_ko_biomarkers.csv")

kegg_uniref_map_file <- "/home/wangmk/UM/Research/COHRA_WGS/lasso_logistic/map_ko_uniref90.txt"
all_uniref90s_file <- "/home/wangmk/UM/Research/COHRA_WGS/saliva_preprocessing/humann_output/unique_uniref90.txt"
all_uniref90s <- read.table(all_uniref90s_file, sep='\t', header=TRUE)
all_uniref90s <- all_uniref90s$Gene.Family[-1]
uniref_species_file <- "/home/wangmk/UM/Research/COHRA_WGS/saliva_preprocessing/humann_output/uniref90_species_source.txt"
uniref_species <- read.table(uniref_species_file, sep='\t', header=FALSE)

KEGG_details <- list()
for (j in 1:nrow(saliva_KEGG_biomarkers)){
  
  selected_kegg <- saliva_KEGG_biomarkers$KEGG[j]
  command <- sprintf("grep %s %s", selected_kegg, kegg_uniref_map_file)
  output <- system(command, intern=TRUE)
  related_unirefs <- strsplit(output, split="\t")[[1]][-1]
  
  subset_unirefs <- intersect(related_unirefs, all_uniref90s)
  KEGG_species_df <- data.frame(Protein=character(0), 
                                Species=character(0))
  for (k in 1:length(subset_unirefs)){
    
    command <- sprintf("grep %s %s", subset_unirefs[k], uniref_species_file)
    output <- system(command, intern=TRUE)[-1]
    relevant_species <- sapply(output, function(longname){
      strsplit(longname, split="[|]")[[1]][2]
    })  
    KEGG_species_df <- rbind(KEGG_species_df, data.frame(Protein=subset_unirefs[k],
                                                         Species=relevant_species))  
  }
  KEGG_details[[selected_kegg]] <- KEGG_species_df
  
}


write.xlsx(KEGG_details, 
           file = "lasso_logistic/KEGG/saliva_KEGG_marker_details.xlsx")


