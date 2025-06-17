
rm(list=ls())
library(BenchmarkDenoise)
library(dplyr)


# load observed counts and metadata
saliva_counts <- read.table("counts_cleaning/strata/saliva_uniref_counts_yr1.tsv",
                            header=TRUE, sep="\t", row.names=1) |> as.matrix()
metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
diagnoses <- (metadata_saliva_yr1$Case_status)
saliva_counts_imputed <- simple_impute(saliva_counts, scale=0.5) |> t()


# load DAA results
DAA_uniref_results <- read.table("DAA/yr1/saliva/saliva_DA_uniref_table.tsv",
                                 sep='\t', header=1)
DAA_uniref_results_positive <- DAA_uniref_results %>% filter(log10foldchange > 0) %>% arrange(pval)
DAA_uniref_results_negative <- DAA_uniref_results %>% filter(log10foldchange < 0) %>% arrange(pval)

marker_uniref <- c(DAA_uniref_results_positive$Taxa[DAA_uniref_results_positive$pval < 0.01],
                   DAA_uniref_results_negative$Taxa[DAA_uniref_results_negative$pval < 0.01])


# start with DAA markers to fit clr-lasso model
saliva_counts <- saliva_counts_imputed[, marker_uniref]

clr_saliva_counts <- clr_transform(saliva_counts)
colnames(clr_saliva_counts) <- colnames(saliva_counts)



coefficients <- matrix(0, nrow=100, ncol=ncol(clr_saliva_counts))
train_auc <- rep(0, 100)
test_auc <- rep(0, 100)
for (j in 1:100){
  print(j)
  model <- fit_logistic_lasso(input=clr_saliva_counts,
                                output=diagnoses, train_prop=0.7, seed=j)
  coefficients[j, ] <- model$coefs[-1]
  train_auc[j] <- model$train_auc
  test_auc[j] <- model$test_auc
  
}

# check variables that are most frequently selected, split into positive and negative
selection_frequency <- colMeans(coefficients != 0)
names(selection_frequency) <- colnames(coefficients) <- colnames(clr_saliva_counts)
hist(selection_frequency, nclass=20)

library(pROC)
AUCs <- matrix(0, nrow=100, ncol=6)
frequency_cutoffs <- seq(0.1, 0.6, 0.1)
variables_selected <- rep("", 6)

for (j in 1:6){
  
  marker_uniref <- names(which(selection_frequency > frequency_cutoffs[j]))
  subset_coefficient_mean <- colMeans(coefficients[, marker_uniref])
  
  
  positive_features <- which(subset_coefficient_mean > 0) |> names()
  negative_features <- which(subset_coefficient_mean < 0) |> names()
  variables_selected[j] <- sprintf("%d+%d", 
                                   length(positive_features), 
                                   length(negative_features))
  
  counts_positive_features <- rowSums(saliva_counts_imputed[, positive_features])
  counts_negative_features <- rowSums(saliva_counts_imputed[, negative_features])
  
  count_ratio <- counts_positive_features / counts_negative_features
  
  ## check the relative abundance of these unirefs
  
  df <- data.frame(cbind(diagnoses, log(count_ratio)))
  colnames(df) <- c("diagnoses", "log_count_ratio")
  
  
  for (k in 1:100){
    
    df_bootstrap <- df[sample(nrow(df), nrow(df)*0.8, replace=T), ]
    
    logistic_regression <- glm(diagnoses ~ log_count_ratio,
                               family=binomial, data=df_bootstrap)
    prediction <- logistic_regression$fitted.values
    AUCs[k, j] <- auc(roc(response=df_bootstrap$diagnoses, 
                          predictor=prediction)) |> suppressMessages()
    
  }

}


performance <- data.frame(Threshold=rep(frequency_cutoffs, each=100),
                          VarSelect=rep(variables_selected, each=100),
                          AUC=as.vector(AUCs))
performance$VarSelect <- factor(performance$VarSelect,
                                levels=variables_selected)

avg_performance <- performance %>% group_by(VarSelect) %>%
  summarise(AUC = mean(AUC),
            Threshold = mean(Threshold))


library(ggplot2)
plot_performance <- ggplot(performance, aes(x=VarSelect, y=AUC)) + 
  geom_boxplot() + 
  xlab("Number of Features") + ylab("AUROC")

best_cutoff <- avg_performance$Threshold[which.max(avg_performance$AUC)]
marker_uniref <- names(which(selection_frequency > best_cutoff))
subset_coefficient_mean <- colMeans(coefficients[, marker_uniref])


positive_features <- which(subset_coefficient_mean > 0) |> names()
negative_features <- which(subset_coefficient_mean < 0) |> names()

predictive_features <- list(positive_features=positive_features,
                            negative_features=negative_features,
                            AUCs = performance$AUC[performance$Threshold == best_cutoff])

saveRDS(predictive_features,
        "lasso_logistic/uniref/saliva_predictive_features_yr1.rds")


# archived: train/test split 

# train_aucs <- rep(0, 100)
# test_aucs <- rep(0, 100)
# coefficients_final <- rep(0, 100)
# # colnames(coefficients_final) <- colnames(log10relabd)
# 
# for(j in 1:100){
#   
#   train_id <- createDataPartition(y=factor(diagnoses), p=0.7, list=FALSE)
#   df_train <- df[train_id, ]
#   df_test <- df[-train_id, ]
#   
#   logistic_regression <- glm(diagnoses ~ log_count_ratio,
#                              family=binomial, data=df_train)
#   
#   coefficients_final[j] <- logistic_regression$coefficients[-1]
#   
#   probabilities_train <- predict(logistic_regression, newdata=df_train,
#                                  type="response")
#   probabilities_test <- predict(logistic_regression, newdata=df_test,
#                                 type="response")
#   
#   roc_obj_train <- roc(response=diagnoses[train_id], predictor=probabilities_train) |> suppressMessages()
#   train_aucs[j] <- auc(roc_obj_train)
#   roc_obj_test <- roc(response=diagnoses[-train_id], predictor=probabilities_test) |> suppressMessages()
#   test_aucs[j] <- auc(roc_obj_test)
#   
# }

