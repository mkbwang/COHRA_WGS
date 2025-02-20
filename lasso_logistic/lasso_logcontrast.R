
rm(list=ls())
# load observed counts and metadata

saliva_relabd <- readRDS("lasso_logistic/feature_filter/saliva_relabd_imputed.rds")
plaque_relabd <- readRDS("lasso_logistic/feature_filter/plaque_relabd_imputed.rds")


metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
diagnoses_1 <- metadata_saliva_yr1$Case_status

metadata_plaque_yr1 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                                  header=T, sep='\t')
diagnoses_2 <- metadata_plaque_yr1$Case_status


source("lasso_logistic/functions_coda_penalized_regression.R")


codalasso_saliva <-  coda_logistic_lasso(y=diagnoses_1, X=saliva_relabd, lambda=0.2)
predicted_values <- codalasso_saliva$betas[1] + 
  log(as.matrix(saliva_relabd)) %*% codalasso_saliva$betas[-1]
predicted_prob <- exp(predicted_values) / (1+exp(predicted_values))

auc(roc(diagnoses_1, predicted_prob))


codalasso_plaque <-  coda_logistic_lasso(y=diagnoses_2, X=plaque_relabd, lambda=0.2)
predicted_values <- codalasso_plaque$betas[1] + 
  log(as.matrix(plaque_relabd)) %*% codalasso_plaque$betas[-1]
predicted_prob <- exp(predicted_values) / (1+exp(predicted_values))
auc(roc(diagnoses_2, predicted_prob))

tune_saliva <-  lambdaRange_codalasso(y=diagnoses_1, X=saliva_relabd, 
                                         lambdaSeq = seq(0.1, 0.5, 0.05))
tune_plaque <- lambdaRange_codalasso(y=diagnoses_2, X=plaque_relabd, 
                                     lambdaSeq = seq(0.05, 0.6, 0.05))



