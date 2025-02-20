

rm(list=ls())
# load observed counts and metadata
saliva_ko_counts_yr1 <- read.table("counts_cleaning/strata/saliva_ko_counts_yr1.tsv",
                                   header=TRUE, sep="\t", row.names=1)
plaque_ko_counts_yr1 <- read.table("counts_cleaning/strata/plaque_ko_counts_yr1.tsv",
                                   header=TRUE, sep="\t", row.names=1)

saliva_DA_KEGG <- readRDS("lasso_logistic/feature_filter/saliva_DA_KEGG.rds")
plaque_DA_KEGG <- readRDS("lasso_logistic/feature_filter/plaque_DA_KEGG.rds")

saliva_counts <- saliva_ko_counts_yr1[, saliva_DA_KEGG]
plaque_counts <- plaque_ko_counts_yr1[, plaque_DA_KEGG]


library(zCompositions)
saliva_relabd_imputed <- cmultRepl(as.matrix(saliva_counts))
plaque_relabd_imputed <- cmultRepl(as.matrix(plaque_counts),
                                   z.warning = 0.9)


saveRDS(saliva_relabd_imputed, 
        file="lasso_logistic/feature_filter/saliva_relabd_imputed.rds")


saveRDS(plaque_relabd_imputed, 
        file="lasso_logistic/feature_filter/plaque_relabd_imputed.rds")










