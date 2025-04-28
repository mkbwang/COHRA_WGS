
rm(list=ls())
library(KEGGREST)
library(stringr)

saliva_ko_counts_yr1 <- read.table("counts_cleaning/strata/saliva_ko_counts_yr1.tsv",
                                   header=TRUE, sep="\t", row.names=1)
KEGG_saliva <- colnames(saliva_ko_counts_yr1)
plaque_ko_counts_yr1 <- read.table("counts_cleaning/strata/plaque_ko_counts_yr1.tsv",
                                   header=TRUE, sep="\t", row.names=1)
KEGG_plaque <- colnames(plaque_ko_counts_yr1)


# limit scope to KEGGs that are related to metabolism
KEGG_union <- union(KEGG_saliva, KEGG_plaque)
retain_indicator <- rep(FALSE, length(KEGG_union))

# kegg_example <- sprintf("ko:%s", KEGG_intersect[1:10])
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
