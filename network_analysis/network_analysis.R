library(SpiecEasi)

run_network <- function(source="saliva", type="taxa", year="yr1", pseudocount=0.5,
                        feature_filter=5){
  
  counts <- read.table(sprintf("counts_cleaning/strata/%s_%s_counts_%s.tsv", source, type, year), 
                       header=T,
                       sep='\t') |> as.matrix()
  
  performance <- read.csv(sprintf("lasso_logistic/%s_%s_%s.csv", source, type, year),
                          row.names = 1)
  
  performance <- performance[-1, ]
  selection_times <- rowSums(performance != 0)
  # only select the ones that are selected more than five times
  subset_features <- names(selection_times[selection_times >= feature_filter])
  
  counts_subset <- counts[, subset_features]
  counts_subset[counts_subset == 0] <- pseudocount
  
  se <- spiec.easi(counts_subset, 
                   method='mb', 
                   lambda.min.ratio=1e-2, nlambda=10)
  
  adjmat <- getRefit(se)
  rownames(adjmat) <- colnames(counts_subset)
  colnames(adjmat) <- colnames(counts_subset)
  
  result <- list(spieceasi=se, adjmat=adjmat)
  
  return(result)
}

saliva_taxa_network <- run_network(source="saliva", type="taxa", year="yr1")
saliva_ko_network <- run_network(source="saliva", type="ko", year="yr1")
plaque_taxa_network <- run_network(source="plaque", type="taxa", year="yr1")
plaque_ko_network <- run_network(source="plaque", type="ko", year="yr1")
  

saliva_taxa_adjmat_dense <- as.matrix(saliva_taxa_network$adjmat)
saliva_ko_adjmat_dense <- as.matrix(saliva_ko_network$adjmat)
plaque_taxa_adjmat_dense <- as.matrix(plaque_taxa_network$adjmat)
plaque_ko_adjmat_dense <- as.matrix(plaque_ko_network$adjmat)


write.table(as.data.frame(saliva_taxa_adjmat_dense), "network_analysis/saliva_taxa_adjmat.tsv",
            sep='\t', quote=F)

write.table(as.data.frame(saliva_ko_adjmat_dense), "network_analysis/saliva_ko_adjmat.tsv",
            sep='\t', quote=F)

write.table(as.data.frame(plaque_ko_adjmat_dense), "network_analysis/plaque_ko_adjmat.tsv",
            sep='\t', quote=F)

write.table(as.data.frame(plaque_taxa_adjmat_dense), "network_analysis/plaque_taxa_adjmat.tsv",
            sep='\t', quote=F)


library(igraph)

taxa_names <- colnames(taxa_adjmat_dense)
taxa_degrees <- rowSums(taxa_adjmat_dense)
subset_taxa <- names(taxa_degrees[taxa_degrees > 3])
graph_taxa <- graph_from_adjacency_matrix(taxa_adjmat_dense[subset_taxa, subset_taxa], 
                                          mode="undirected")
# taxa_labels <- taxa_names
# taxa_labels[!(taxa_labels %in% subset_taxa)] <- NA
par(mar = c(1, 1, 1, 1))
plot(graph_taxa, # vertex.label = taxa_labels,
     vertex.label.cex = 0.8,  vertex.size = 1,
     layout = layout_with_kk)


ko_names <- colnames(ko_adjmat_dense)
ko_degrees <- rowSums(ko_adjmat_dense)
subset_ko <- names(ko_degrees[ko_degrees > 4])
graph_ko <- graph_from_adjacency_matrix(ko_adjmat_dense[subset_ko, subset_ko], 
                                          mode="undirected")
# taxa_labels <- taxa_names
# taxa_labels[!(taxa_labels %in% subset_taxa)] <- NA
par(mar = c(1, 1, 1, 1))
plot(graph_ko, # vertex.label = taxa_labels,
     vertex.label.cex = 0.8,  vertex.size = 1,
     layout = layout_with_kk)


