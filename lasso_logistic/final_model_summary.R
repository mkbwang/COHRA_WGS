
library(BenchmarkDenoise)
library(ggplot2)
library(pROC)
rm(list=ls())
phenotype <- c("Control", "Case")
metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                                  header=T, sep='\t')
diagnoses_saliva_yr1 <- metadata_saliva_yr1$Case_status
metadata_plaque_yr1 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                           header=T, sep='\t')
diagnoses_plaque_yr1 <- metadata_plaque_yr1$Case_status


summarize_result <- function(outcome, count_table, marker,
                             input, type){
  
  count_table_imputed <- simple_impute(count_table) |> t()
  colnames(count_table_imputed) <- colnames(count_table)
  count_table_subset <- count_table_imputed[, c(marker$negative_features, marker$positive_features)]
  count_ratios <- rowSums(count_table_subset[, marker$positive_features, drop=FALSE]) / 
    rowSums(count_table_subset[, marker$negative_features, drop=FALSE])
  log_count_ratios <- log10(count_ratios)
  result_df <- data.frame(Diagnoses=phenotype[outcome+1],
                          diagnoses_number=outcome,
                          LogRatio=log_count_ratios) 
  biomarker_boxplot <- ggplot(result_df, aes(x=Diagnoses, y=LogRatio)) + 
    geom_boxplot() + xlab("Diagnoses") + ylab("Biomarker Value") +
    ggtitle(sprintf("%s %s", input, type))
  
  logistic_regression <- glm(diagnoses_number ~ LogRatio, family="binomial",
                             data=result_df)
  predicted_value <- logistic_regression$fitted.values
  roc_curve <- roc(response=outcome,
                   predictor=predicted_value) |> suppressMessages()
  AUC <- auc(roc_curve)
  roc_coords <- roc_curve$thresholds
  true_positive_rates <- roc_curve$sensitivities
  false_positive_rates <- 1 - roc_curve$specificities
  roc_coordinates <- data.frame(
    Source=input,
    Type=type,
    TPR = sort(true_positive_rates),
    FPR = sort(false_positive_rates)
  )
  
  
  output <- list(biomarkers=marker, biomarker_plot=biomarker_boxplot,
                 roc_coordinates=roc_coordinates, AUC=AUC)
  
  return(output)
  
}


map_uniref_ko <- function(uniref, mapping_file){
  
  command <- sprintf("grep %s %s", uniref, mapping_file)
  output <- system(command, intern=TRUE)
  kegg <- strsplit(output, split="\t")[[1]][1]
  return(kegg)
  
}


map_uniref_species <- function(uniref, mapping_file){
  
  command <- sprintf("grep %s %s", uniref, mapping_file)
  output <- system(command, intern=TRUE)
  species <- sapply(output, function(y) strsplit(y, "[|]")[[1]][2])
  species <- species[!is.na(species)]
  all_species <- paste(species, collapse=";")
  return(all_species)
  
}



# load taxa counts
saliva_taxa_counts <- read.table("counts_cleaning/strata/saliva_taxa_counts_yr1.tsv",
                                 header=TRUE, sep="\t", row.names=1) |> as.matrix()
saliva_predictive_taxa <- readRDS("lasso_logistic/taxa/saliva_predictive_features_yr1.rds")
saliva_taxa_results <- summarize_result(outcome=diagnoses_saliva_yr1,
                                        count_table=saliva_taxa_counts,
                                        marker=saliva_predictive_taxa,
                                        input="Saliva",
                                        type="Taxa")



plaque_taxa_counts <- read.table("counts_cleaning/strata/plaque_taxa_counts_yr1.tsv",
                                 header=TRUE, sep="\t", row.names=1) |> as.matrix()
plaque_predictive_taxa <- readRDS("lasso_logistic/taxa/plaque_predictive_features_yr1.rds")
plaque_taxa_results <- summarize_result(outcome=diagnoses_plaque_yr1,
                                        count_table=plaque_taxa_counts,
                                        marker=plaque_predictive_taxa,
                                        input="Plaque",
                                        type="Taxa")


# load KEGG counts
saliva_ko_counts <- read.table("counts_cleaning/strata/saliva_ko_counts_yr1.tsv",
                                 header=TRUE, sep="\t", row.names=1) |> as.matrix()
saliva_predictive_ko <- readRDS("lasso_logistic/KEGG/saliva_predictive_features_yr1.rds")
saliva_ko_results <- summarize_result(outcome=diagnoses_saliva_yr1,
                                        count_table=saliva_ko_counts,
                                        marker=saliva_predictive_ko,
                                        input="Saliva",
                                        type="KEGG")


plaque_ko_counts <- read.table("counts_cleaning/strata/plaque_ko_counts_yr1.tsv",
                               header=TRUE, sep="\t", row.names=1) |> as.matrix()
plaque_predictive_ko <- readRDS("lasso_logistic/KEGG/plaque_predictive_features_yr1.rds")
plaque_ko_results <- summarize_result(outcome=diagnoses_plaque_yr1,
                                      count_table=plaque_ko_counts,
                                      marker=plaque_predictive_ko,
                                      input="Plaque",
                                      type="KEGG")



# load uniref counts
saliva_uniref_counts <- read.table("counts_cleaning/strata/saliva_uniref_counts_yr1.tsv",
                               header=TRUE, sep="\t", row.names=1) |> as.matrix()
saliva_predictive_uniref <- readRDS("lasso_logistic/uniref/saliva_predictive_features_yr1.rds")
saliva_uniref_results <- summarize_result(outcome=diagnoses_saliva_yr1,
                                      count_table=saliva_uniref_counts,
                                      marker=saliva_predictive_uniref,
                                      input="Saliva",
                                      type="Uniref")


ko_ref_file <- "/home/wangmk/UM/Research/COHRA_WGS/lasso_logistic/map_ko_uniref90.txt"
species_ref_file <- "/home/wangmk/UM/Research/COHRA_WGS/lasso_logistic/uniref/saliva_uniref90_details.txt"
saliva_positive_unirefs <- saliva_uniref_results$biomarkers$positive_features
saliva_positive_keggs <- rep("", length(saliva_positive_unirefs))
saliva_positive_species <- rep("", length(saliva_positive_unirefs))
for (j in 1:length(saliva_positive_unirefs)){
  saliva_positive_keggs[j] <- map_uniref_ko(saliva_positive_unirefs[j],
                                            ko_ref_file)
  saliva_positive_species[j] <- map_uniref_species(saliva_positive_unirefs[j],
                                                species_ref_file)
}


saliva_positive_df <- data.frame(Unirefs=saliva_positive_unirefs,
                                 KEGGs=saliva_positive_keggs,
                                 Species=saliva_positive_species)
write.csv(saliva_positive_df, "lasso_logistic/uniref/saliva_positive.csv",
          row.names=FALSE)


saliva_negative_unirefs <- saliva_uniref_results$biomarkers$negative_features
saliva_negative_keggs <- rep("", length(saliva_negative_unirefs))
saliva_negative_species <- rep("", length(saliva_negative_unirefs))
for (j in 1:length(saliva_negative_unirefs)){
  saliva_negative_keggs[j] <- map_uniref_ko(saliva_negative_unirefs[j],
                                            ko_ref_file)
  saliva_negative_species[j] <- map_uniref_species(saliva_negative_unirefs[j],
                                                   species_ref_file)
}


saliva_negative_df <- data.frame(Unirefs=saliva_negative_unirefs,
                                 KEGGs=saliva_negative_keggs,
                                 Species=saliva_negative_species)
write.csv(saliva_negative_df, "lasso_logistic/uniref/saliva_negative.csv",
          row.names=FALSE)



plaque_uniref_counts <- read.table("counts_cleaning/strata/plaque_uniref_counts_yr1.tsv",
                               header=TRUE, sep="\t", row.names=1) |> as.matrix()
plaque_predictive_uniref <- readRDS("lasso_logistic/uniref/plaque_predictive_features_yr1.rds")
plaque_uniref_results <- summarize_result(outcome=diagnoses_plaque_yr1,
                                      count_table=plaque_uniref_counts,
                                      marker=plaque_predictive_uniref,
                                      input="Plaque",
                                      type="Uniref")


ko_ref_file <- "/home/wangmk/UM/Research/COHRA_WGS/lasso_logistic/map_ko_uniref90.txt"
species_ref_file <- "/home/wangmk/UM/Research/COHRA_WGS/lasso_logistic/uniref/plaque_uniref90_details.txt"
plaque_positive_unirefs <- plaque_uniref_results$biomarkers$positive_features
plaque_positive_keggs <- rep("", length(plaque_positive_unirefs))
plaque_positive_species <- rep("", length(plaque_positive_unirefs))
for (j in 1:length(plaque_positive_unirefs)){
  plaque_positive_keggs[j] <- map_uniref_ko(plaque_positive_unirefs[j],
                                            ko_ref_file)
  plaque_positive_species[j] <- map_uniref_species(plaque_positive_unirefs[j],
                                                   species_ref_file)
}


plaque_positive_df <- data.frame(Unirefs=plaque_positive_unirefs,
                                 KEGGs=plaque_positive_keggs,
                                 Species=plaque_positive_species)
write.csv(plaque_positive_df, "lasso_logistic/uniref/plaque_positive.csv",
          row.names=FALSE)


plaque_negative_unirefs <- plaque_uniref_results$biomarkers$negative_features
plaque_negative_keggs <- rep("", length(plaque_negative_unirefs))
plaque_negative_species <- rep("", length(plaque_negative_unirefs))
for (j in 1:length(plaque_negative_unirefs)){
  plaque_negative_keggs[j] <- map_uniref_ko(plaque_negative_unirefs[j],
                                            ko_ref_file)
  plaque_negative_species[j] <- map_uniref_species(plaque_negative_unirefs[j],
                                                   species_ref_file)
}


plaque_negative_df <- data.frame(Unirefs=plaque_negative_unirefs,
                                 KEGGs=plaque_negative_keggs,
                                 Species=plaque_negative_species)
write.csv(plaque_negative_df, "lasso_logistic/uniref/plaque_negative.csv",
          row.names=FALSE)



roc_curves <- rbind(saliva_taxa_results$roc_coordinates,
                    plaque_taxa_results$roc_coordinates,
                    saliva_ko_results$roc_coordinates,
                    plaque_ko_results$roc_coordinates,
                    saliva_uniref_results$roc_coordinates,
                    plaque_uniref_results$roc_coordinates)


roc_curves_plot <- ggplot(roc_curves, aes(x=FPR, y=TPR)) +
  geom_line(aes(color=Source, linetype=Type))

library(patchwork)

marker_plots <- (saliva_taxa_results$biomarker_plot + saliva_ko_results$biomarker_plot + saliva_uniref_results$biomarker_plot) /
  (plaque_taxa_results$biomarker_plot + plaque_ko_results$biomarker_plot + plaque_uniref_results$biomarker_plot) +
  plot_layout(axis_titles = "collect")



