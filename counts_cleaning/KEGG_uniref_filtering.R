
rm(list=ls())
library(KEGGREST)
library(stringr)

saliva_ko_counts_yr1 <- read.csv("counts_cleaning/saliva_ko_abundance_subset.csv", 
                                 row.names=1)
KEGG_saliva <- rownames(saliva_ko_counts_yr1)
plaque_ko_counts_yr1 <- read.csv("counts_cleaning/plaque_ko_abundance_subset.csv", 
                                 row.names=1)
KEGG_plaque <- rownames(plaque_ko_counts_yr1)


# limit scope to KEGGs that are related to metabolism
KEGG_union <- union(KEGG_saliva, KEGG_plaque)
retain_indicator <- rep(FALSE, length(KEGG_union))


begin <- proc.time()
for(k in 1:length(KEGG_union)){
  if(k %% 10 == 0) print(k)
  retain_indicator[k] = tryCatch({
    query <-keggGet(sprintf("ko:%s", KEGG_union[k]))
    any(str_detect(query[[1]]$BRITE, "(?i)metabolism"))
  }, warning=function(w){
    message(sprintf("Warning for %d", k), conditionMessage(w))
    FALSE
  }, error=function(e){
    message(sprintf("Error for %d", k), conditionMessage(e))
    FALSE
  }, finally={}
  )
  Sys.sleep(0.2)
}
end <- proc.time()
KEGG_selected <- KEGG_union[retain_indicator]
KEGG_df <- data.frame(KEGG=KEGG_selected)
write.csv(KEGG_df, file="counts_cleaning/metabolism_KOs.csv", row.names=FALSE, quote=FALSE)


metabolism_KOs <- read.csv("counts_cleaning/metabolism_KOs.csv")
subset_saliva_ko <- intersect(metabolism_KOs$KEGG, rownames(saliva_ko_counts_yr1))
saliva_ko_counts_yr1 <- saliva_ko_counts_yr1[subset_saliva_ko, ]

subset_plaque_ko <- intersect(metabolism_KOs$KEGG, rownames(plaque_ko_counts_yr1))
plaque_ko_counts_yr1 <- plaque_ko_counts_yr1[subset_plaque_ko, ]

write.csv(saliva_ko_counts_yr1, "counts_cleaning/saliva_ko_abundance_subset.csv",
          quote=FALSE)

write.csv(plaque_ko_counts_yr1, "counts_cleaning/plaque_ko_abundance_subset.csv",
          quote=FALSE)



# load the unirefs that map to each KEGGs
saliva_unirefs <- read.csv("counts_cleaning/saliva_valid_unirefs.txt", header=FALSE)
saliva_unirefs <- sprintf("UniRef90_%s", saliva_unirefs$V1)


saliva_ko_uniref_match <- data.frame(KEGGs=subset_saliva_ko,
                                     unirefs="")
for (j in 1:length(subset_saliva_ko)){
  if (j %%10 == 0){print(j)}
  command <- sprintf("grep -i '%s' /home/wangmk/UM/Research/COHRA_WGS/plaque_preprocessing/humann_output/map_ko_uniref90.txt",
                     subset_saliva_ko[j])
  result <- system(command, intern=T)
  unirefs <- strsplit(result, split="\t")[[1]]
  existing_unirefs <- intersect(unirefs, saliva_unirefs)
  saliva_ko_uniref_match$unirefs[j] <- paste(existing_unirefs, collapse="|")
  
}
write.csv(saliva_ko_uniref_match,
          "counts_cleaning/saliva_ko_uniref_match.csv", quote=FALSE, row.names=FALSE)

saliva_ko_uniref_match <- read.csv("counts_cleaning/saliva_ko_uniref_match.csv")
saliva_subset_unirefs <- paste(saliva_ko_uniref_match$unirefs, collapse="|")
saliva_subset_unirefs <- strsplit(saliva_subset_unirefs, split="[|]")[[1]]
saliva_subset_unirefs <- unique(saliva_subset_unirefs)
writeLines(saliva_subset_unirefs, "saliva_preprocessing/humann_output/KEGG_matched_unirefs.txt")


plaque_unirefs <- read.csv("counts_cleaning/plaque_valid_unirefs.txt", header=FALSE)
plaque_unirefs <- sprintf("UniRef90_%s", plaque_unirefs$V1)
plaque_ko_uniref_match <- data.frame(KEGGs=subset_plaque_ko,
                                     unirefs="")
for (j in 1:length(subset_plaque_ko)){
  if (j %%10 == 0){print(j)}
  command <- sprintf("grep -i '%s' /home/wangmk/UM/Research/COHRA_WGS/plaque_preprocessing/humann_output/map_ko_uniref90.txt",
                     subset_plaque_ko[j])
  result <- system(command, intern=T)
  unirefs <- strsplit(result, split="\t")[[1]]
  existing_unirefs <- intersect(unirefs, plaque_unirefs)
  plaque_ko_uniref_match$unirefs[j] <- paste(existing_unirefs, collapse="|")
  
}
write.csv(plaque_ko_uniref_match,
          "counts_cleaning/plaque_ko_uniref_match.csv", quote=FALSE, row.names=FALSE)

plaque_ko_uniref_match <- read.csv("counts_cleaning/plaque_ko_uniref_match.csv")
plaque_subset_unirefs <- paste(plaque_ko_uniref_match$unirefs, collapse="|")
plaque_subset_unirefs <- strsplit(plaque_subset_unirefs, split="[|]")[[1]]
plaque_subset_unirefs <- unique(plaque_subset_unirefs)
writeLines(plaque_subset_unirefs, "plaque_preprocessing/humann_output/KEGG_matched_unirefs.txt")



## based on the filtered unirefs, recalibrate the KEGG orthologs
rm(list=ls())
saliva_ko_counts <- read.csv("counts_cleaning/saliva_ko_abundance_subset.csv", 
                                 row.names=1)
KEGG_saliva <- rownames(saliva_ko_counts)
plaque_ko_counts <- read.csv("counts_cleaning/plaque_ko_abundance_subset.csv", 
                                 row.names=1)
KEGG_plaque <- rownames(plaque_ko_counts)


plaque_ko_uniref_match <- read.csv("counts_cleaning/plaque_ko_uniref_match.csv")
plaque_ko_uniref_match_long <- list()
for (j in 1:nrow(plaque_ko_uniref_match)){
  
  unirefs <- strsplit(plaque_ko_uniref_match$unirefs[j], "[|]")[[1]]
  kegg <- plaque_ko_uniref_match$KEGGs[j]
  plaque_ko_uniref_match_long[[j]] <- data.frame(Protein=unirefs, KEGG=kegg)
    
}
plaque_ko_uniref_match_long <- do.call(rbind, plaque_ko_uniref_match_long)
plaque_uniref_counts <- read.csv("counts_cleaning/plaque_uniref90_abundance_subset.csv",
                                 row.names=1)
unirefs <- rownames(plaque_uniref_counts)
plaque_ko_uniref_match_long <- plaque_ko_uniref_match_long %>% filter(Protein %in% unirefs)
write.csv(plaque_ko_uniref_match_long,
          "counts_cleaning/plaque_ko_uniref_match_long.csv", row.names=FALSE)

plaque_new_ko_table <- matrix(0, nrow=length(unique(plaque_ko_uniref_match_long$KEGG)),
                              ncol=ncol(plaque_uniref_counts))
rownames(plaque_new_ko_table) <- unique(plaque_ko_uniref_match_long$KEGG)
colnames(plaque_new_ko_table) <- colnames(plaque_uniref_counts)
for (j in 1:nrow(plaque_new_ko_table)){
  
  kegg <- rownames(plaque_new_ko_table)[j]
  subset_proteins <- plaque_ko_uniref_match_long %>% filter(KEGG==kegg)
  plaque_new_ko_table[j, ] <- colSums(plaque_uniref_counts[subset_proteins$Protein, ])
  
}
plaque_new_ko_table <- as.data.frame(plaque_new_ko_table)
write.csv(plaque_new_ko_table, "counts_cleaning/plaque_ko_abundance_subset.csv")



## repeat for saliva
saliva_ko_uniref_match <- read.csv("counts_cleaning/saliva_ko_uniref_match.csv")
saliva_ko_uniref_match_long <- list()
for (j in 1:nrow(saliva_ko_uniref_match)){
  
  unirefs <- strsplit(saliva_ko_uniref_match$unirefs[j], "[|]")[[1]]
  kegg <- saliva_ko_uniref_match$KEGGs[j]
  saliva_ko_uniref_match_long[[j]] <- data.frame(Protein=unirefs, KEGG=kegg)
  
}
saliva_ko_uniref_match_long <- do.call(rbind, saliva_ko_uniref_match_long)
saliva_uniref_counts <- read.csv("counts_cleaning/saliva_uniref90_abundance_subset.csv",
                                 row.names=1)
unirefs <- rownames(saliva_uniref_counts)
saliva_ko_uniref_match_long <- saliva_ko_uniref_match_long %>% filter(Protein %in% unirefs)
write.csv(saliva_ko_uniref_match_long,
          "counts_cleaning/saliva_ko_uniref_match_long.csv", row.names=FALSE)

saliva_new_ko_table <- matrix(0, nrow=length(unique(saliva_ko_uniref_match_long$KEGG)),
                              ncol=ncol(saliva_uniref_counts))
rownames(saliva_new_ko_table) <- unique(saliva_ko_uniref_match_long$KEGG)
colnames(saliva_new_ko_table) <- colnames(saliva_uniref_counts)
for (j in 1:nrow(saliva_new_ko_table)){
  
  kegg <- rownames(saliva_new_ko_table)[j]
  subset_proteins <- saliva_ko_uniref_match_long %>% filter(KEGG==kegg)
  saliva_new_ko_table[j, ] <- colSums(saliva_uniref_counts[subset_proteins$Protein, ])
  
}
saliva_new_ko_table <- as.data.frame(saliva_new_ko_table)
write.csv(saliva_new_ko_table, "counts_cleaning/saliva_ko_abundance_subset.csv")


