
rm(list=ls())
# load observed counts and metadata

saliva_relabd <- readRDS("lasso_logistic/feature_filter/saliva_relabd_imputed.rds")
plaque_relabd <- readRDS("lasso_logistic/feature_filter/plaque_relabd_imputed.rds")

saliva_counts <- read.table("counts_cleaning/strata/saliva_ko_counts_yr1.tsv",
                            header=TRUE, sep="\t", row.names=1)
plaque_counts <- read.table("counts_cleaning/strata/plaque_ko_counts_yr1.tsv",
                            header=TRUE, sep='\t', row.names=1)


metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
diagnoses_1 <- metadata_saliva_yr1$Case_status

metadata_plaque_yr1 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                                  header=T, sep='\t')
diagnoses_2 <- metadata_plaque_yr1$Case_status


source("lasso_logistic/functions_coda_penalized_regression.R")
source("lasso_logistic/functions_lasso_tune.R")

library(ggplot2)
library(patchwork)


codalasso_saliva <- lambda_tune(y=diagnoses_1, lambdas=seq(0.05, 0.50, 0.05), X=saliva_relabd)
saveRDS(codalasso_saliva, file="lasso_logistic/feature_filter/codalasso_saliva.rds")

codalasso_plaque <- lambda_tune(y=diagnoses_2, lambdas=seq(0.05, 0.50, 0.05), X=plaque_relabd)
saveRDS(codalasso_plaque, file="lasso_logistic/feature_filter/codalasso_plaque.rds")




saliva_plot <- performance_plot(codalasso_saliva)
plaque_plot <- performance_plot(codalasso_plaque)

library(dplyr)
saliva_choice <- data.frame(Features=colnames(saliva_relabd), 
                            Choice_prob=codalasso_saliva$choice_probs[7, ]) |> 
  arrange(desc(Choice_prob)) |> head(10)

saliva_choice$Features <- factor(saliva_choice$Features, 
                                 levels=rev(saliva_choice$Features))

saliva_choice_plot <- ggplot(saliva_choice, aes(x=Features, y=Choice_prob)) + 
  geom_bar(stat="identity") + coord_flip() + 
  labs(x="Features", y="Selection probability") + 
  theme_bw()


plaque_choice <- data.frame(Features=colnames(plaque_relabd), 
                                        Choice_prob=codalasso_plaque$choice_probs[3, ]) |> 
  arrange(desc(Choice_prob)) |> head(10)
plaque_choice$Features <- factor(plaque_choice$Features, 
                                 levels=rev(plaque_choice$Features))

plaque_choice_plot <- ggplot(plaque_choice, aes(x=Features, y=Choice_prob)) + 
  geom_bar(stat="identity") + coord_flip() + 
  labs(x="Features", y="Selection probability") + 
  theme_bw()
