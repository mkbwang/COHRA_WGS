
library(readxl)
library(dplyr)
library(ggraph)
library(igraph)
library(ggplot2)

# based on the initial lasso log contrast model, select an even smaller set of KEGGs
rm(list=ls())

# load data
saliva_ko_counts_yr1 <- read.table("counts_cleaning/strata/saliva_ko_counts_yr1.tsv",
                                   header=TRUE, sep="\t", row.names=1)
plaque_ko_counts_yr1 <- read.table("counts_cleaning/strata/plaque_ko_counts_yr1.tsv",
                                   header=TRUE, sep="\t", row.names=1)

saliva_relabd <- readRDS("lasso_logistic/feature_filter/saliva_relabd_imputed.rds")
plaque_relabd <- readRDS("lasso_logistic/feature_filter/plaque_relabd_imputed.rds")
metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
diagnoses_1 <- metadata_saliva_yr1$Case_status
metadata_plaque_yr1 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                                  header=T, sep='\t')
diagnoses_2 <- metadata_plaque_yr1$Case_status

# load lasso log contrast results
codalasso_saliva <- readRDS("lasso_logistic/feature_filter/codalasso_saliva.rds")
codalasso_plaque <- readRDS("lasso_logistic/feature_filter/codalasso_plaque.rds")
saliva_KEGG_probs <- codalasso_saliva$choice_probs
colnames(saliva_KEGG_probs) <- colnames(saliva_relabd)
saliva_KEGG_importance <- saliva_KEGG_probs[7, ]
plaque_KEGG_probs <- codalasso_plaque$choice_probs
colnames(plaque_KEGG_probs) <- colnames(plaque_relabd)
plaque_KEGG_importance <- plaque_KEGG_probs[6, ]


# KEGG detailed information and DAA result
saliva_KEGG_info <- read_excel("lasso_logistic/feature_filter/KEGG_Information.xlsx",
                               sheet="Saliva")
saliva_DAA_result <- read.csv("lasso_logistic/feature_filter/saliva_DAA_details.csv")

plaque_KEGG_info <- read_excel("lasso_logistic/feature_filter/KEGG_Information.xlsx",
                               sheet="Plaque")
plaque_DAA_result <- read.csv("lasso_logistic/feature_filter/plaque_DAA_details.csv")




predictive_feature_selection <- function(KEGG_info, KEGG_importance, 
                                         KEGG_counts, KEGG_relabd, KEGG_DAA_result){
  
  ## KEGG functionality is important, should have good classification
  KEGG_info <- KEGG_info |> filter(Function !="Unclassified")
  
  ## KEGG with high prevalence are worth considering
  KEGG_count_subset <- KEGG_counts[, KEGG_info$KEGG]
  KEGG_relabd_subset <- KEGG_relabd[, KEGG_info$KEGG]
  KEGG_info$prevalence <- colMeans(KEGG_count_subset != 0)
  
  ## KEGG with large DAA effect size are worth considering
  KEGG_DAA_result <- KEGG_DAA_result |> filter(Taxa %in% KEGG_info$KEGG)
  KEGG_info$logfold <- KEGG_DAA_result$log10foldchange
  
  ## KEGGs that are frequently selected into lasso log contrast are worth considering
  KEGG_info$importance <- KEGG_importance[KEGG_info$KEGG]
  
  ## network analysis to find highly correlated pairs
  library(SpiecEasi)
  KEGG_relabd_se <- spiec.easi(as.matrix(KEGG_relabd_subset), method="mb",
                        lambda.min.ratio=1e-2, nlambda=40)
  KEGG_relabd_network <- getRefit(KEGG_relabd_se) |> as.matrix()
  rownames(KEGG_relabd_network) <- colnames(KEGG_relabd_network) <- colnames(KEGG_relabd_subset)
  KEGG_relabd_graph <- graph_from_adjacency_matrix(KEGG_relabd_network, mode="undirected")
  KEGG_relabd_graph_plot <- ggraph(KEGG_relabd_graph, layout="fr")+
    geom_edge_link(aes(edge_alpha = 0.5), show.legend = FALSE) +
    geom_node_point(size = 5, color = "steelblue") +
    geom_node_text(aes(label = name), vjust = -0.5, size=3.5) +
    theme_void()
  
  
  KEGG_count_se <- spiec.easi(as.matrix(KEGG_count_subset), method="mb",
                               lambda.min.ratio=1e-2, nlambda=40)
  KEGG_count_network <- getRefit(KEGG_count_se) |> as.matrix()
  rownames(KEGG_count_network) <- colnames(KEGG_count_network) <- colnames(KEGG_count_subset)
  KEGG_count_graph <- graph_from_adjacency_matrix(KEGG_count_network, mode="undirected")
  KEGG_count_graph_plot <- ggraph(KEGG_count_graph, layout="fr")+
    geom_edge_link(aes(edge_alpha = 0.5), show.legend = FALSE) +
    geom_node_point(size = 5, color = "steelblue") +
    geom_node_text(aes(label = name), vjust = -0.5, size=3.5) +
    theme_void()
  
  
  return(list(KEGG_info=KEGG_info, KEGG_relabd_graph=KEGG_relabd_graph_plot,
              KEGG_count_graph=KEGG_count_graph_plot))

}


saliva_outcome <- predictive_feature_selection(KEGG_info = saliva_KEGG_info,
                                               KEGG_importance = saliva_KEGG_importance,
                                               KEGG_counts = saliva_ko_counts_yr1,
                                               KEGG_relabd = saliva_relabd,
                                               KEGG_DAA_result = saliva_DAA_result)

library(dplyr)
saliva_KEGG_info <- saliva_outcome$KEGG_info
saliva_KEGG_info <- saliva_KEGG_info |> select(KEGG, Function, Type, importance, logfold, prevalence, Name, Notes) |>
  arrange(desc(importance))

write.csv(saliva_KEGG_info, 
          "lasso_logistic/feature_filter/saliva_KEGG_info.csv", row.names=FALSE)


plaque_outcome <- predictive_feature_selection(KEGG_info = plaque_KEGG_info,
                                               KEGG_importance = plaque_KEGG_importance,
                                               KEGG_counts = plaque_ko_counts_yr1,
                                               KEGG_relabd = plaque_relabd,
                                               KEGG_DAA_result = plaque_DAA_result)

plaque_KEGG_info <- plaque_outcome$KEGG_info
plaque_KEGG_info <- plaque_KEGG_info |> select(KEGG, Function, Type, importance, logfold, prevalence, Name, Notes) |>
  arrange(desc(importance))


write.csv(plaque_KEGG_info, 
          "lasso_logistic/feature_filter/plaque_KEGG_info.csv", row.names=FALSE)



View(saliva_outcome$KEGG_info[, c("KEGG", "Function", "Type", "prevalence", "logfold", "importance")])
View(plaque_outcome$KEGG_info[, c("KEGG", "Function", "Type", "prevalence", "logfold", "importance")])

