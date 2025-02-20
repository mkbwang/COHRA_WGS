


saliva_ko_counts_yr1 <- read.table("counts_cleaning/strata/saliva_ko_counts_yr1.tsv",
                                     header=TRUE, sep="\t", row.names=1)
KEGG_saliva <- colnames(saliva_ko_counts_yr1)
plaque_ko_counts_yr1 <- read.table("counts_cleaning/strata/plaque_ko_counts_yr1.tsv",
                                   header=TRUE, sep="\t", row.names=1)
KEGG_plaque <- colnames(plaque_ko_counts_yr1)


# limit scope to KEGGs that are related to metabolism
KEGG_intersect <- intersect(KEGG_saliva, KEGG_plaque)
retain_indicator <- rep(FALSE, length(KEGG_intersect))

# kegg_example <- sprintf("ko:%s", KEGG_intersect[1:10])
begin <- proc.time()
for(k in 1:length(KEGG_intersect)){
  if(k %% 10 == 0) print(k)
  retain_indicator[k] = tryCatch({
    query <-keggGet(sprintf("ko:%s", KEGG_intersect[k]))
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
KEGG_selected <- KEGG_intersect[retain_indicator]
saveRDS(KEGG_selected, file="lasso_logistic/feature_filter/metabolism_KOs.rds")

metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
metadata_plaque_yr1 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                                  header=T, sep='\t')

library(phyloseq)
library(ADAPT)
library(dplyr)

# saliva DA KEGG
saliva_ko_counts_subset <- saliva_ko_counts_yr1[, KEGG_selected]
phy_saliva <- phyloseq(otu_table(saliva_ko_counts_subset, taxa_are_rows=FALSE),
                       sample_data(metadata_saliva_yr1))
censor_val <- min(saliva_ko_counts_subset[saliva_ko_counts_subset > 0])
DAA_saliva <- adapt(phy_saliva, cond.var="Case_status", censor=censor_val,
                    prev.filter=0.2, alpha=0.1)
saliva_DA_KEGG <- DAA_saliva@signal
saveRDS(saliva_DA_KEGG, file="lasso_logistic/feature_filter/saliva_DA_KEGG.rds")


# plaque DA KEGG
plaque_ko_counts_subset <- plaque_ko_counts_yr1[, KEGG_selected]
phy_plaque <- phyloseq(otu_table(plaque_ko_counts_subset, taxa_are_rows=FALSE),
                       sample_data(metadata_plaque_yr1))
censor_val <- min(plaque_ko_counts_subset[plaque_ko_counts_subset > 0])
DAA_plaque <- adapt(phy_plaque, cond.var="Case_status", censor=censor_val,
                    prev.filter=0.2, alpha=0.1)
plaque_DA_KEGG <- DAA_plaque@signal
saveRDS(plaque_DA_KEGG, file="lasso_logistic/feature_filter/plaque_DA_KEGG.rds")




