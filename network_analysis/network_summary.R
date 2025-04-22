
library(ggplot2)

# plot the number of times selected against node degree
network_summary <- function(source="saliva", type="taxa", year="yr1", title="Saliva Taxa"){
  
  regression_selection <- read.csv(sprintf("lasso_logistic/%s_%s_%s.csv", source, type, year),
                                   row.names = 1)
  adjmat <- read.table(sprintf("network_analysis/%s_%s_adjmat.tsv", source, type),sep='\t', 
                       row.names = 1)
  regression_selection_subset <- regression_selection[rownames(adjmat), ]
  selected_times <- rowSums(regression_selection_subset != 0)
  graph_degrees <- rowSums(adjmat)
  
  summary_df <- data.frame(Names=rownames(adjmat),
                           lasso_selection = selected_times,
                           graph_degree = graph_degrees)
  
  generated_plot <- ggplot(summary_df, aes(x=lasso_selection, y=graph_degree)) + 
    geom_point() + 
    scale_x_continuous(name="Number of times selected in LASSO logistic regression", 
                       limits=c(0, 100), breaks=seq(0, 100, 10))+
    scale_y_continuous(name="Node degree in correlation rraph",
                       limits=c(0, max(graph_degrees)), breaks=seq(0, max(graph_degrees)))+
    ggtitle(title)
  
}


saliva_taxa_plot <- network_summary(source="saliva", 
                                     type="taxa", 
                                     year="yr1", 
                                     title="Saliva Taxa")

saliva_gene_plot <- network_summary(source="saliva",
                                    type="ko",
                                    year="yr1",
                                    title="Saliva KEGG")


plaque_taxa_plot <- network_summary(source="plaque",
                                    type="taxa",
                                    year="yr1",
                                    title="Plaque Taxa")


plaque_gene_plot <- network_summary(source="plaque",
                                    type="ko",
                                    year="yr1",
                                    title="Plaque KEGG")


