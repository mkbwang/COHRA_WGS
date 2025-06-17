
rm(list=ls())
library(BenchmarkDenoise)

# load observed counts and metadata
plaque_counts <- read.table("counts_cleaning/strata/plaque_taxa_counts_yr1.tsv",
                            header=TRUE, sep="\t", row.names=1) |> as.matrix()
metadata_plaque_yr1 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                                  header=T, sep='\t')
diagnoses <- (metadata_plaque_yr1$Case_status)
plaque_counts_imputed <- simple_impute(plaque_counts, scale=0.5) |> t()


# load DAA results
DAA_taxa_results <- read.table("DAA/yr1/plaque/taxa/plaque_DA_taxa_table.tsv",
                             sep='\t', header=1)
marker_taxa <- DAA_taxa_results$Taxa[DAA_taxa_results$pval < 0.05]




# start with DAA markers
plaque_counts <- plaque_counts[, marker_taxa]
clr_plaque_counts <- clr_transform(plaque_counts)
colnames(clr_plaque_counts) <- colnames(plaque_counts)
coefficients <- matrix(0, nrow=100, ncol=ncol(clr_plaque_counts))
colnames(coefficients) <- colnames(clr_plaque_counts)
train_auc <- rep(0, 100)
test_auc <- rep(0, 100)
for (j in 1:100){
  print(j)
  model <- fit_logistic_lasso(input=clr_plaque_counts,
                                output=diagnoses, train_prop=0.7, seed=j)
  coefficients[j, ] <- model$coefs[-1]
  train_auc[j] <- model$train_auc
  test_auc[j] <- model$test_auc
  
}

# check variables that are most frequently selected, split into positive and negative
selection_frequency <- colMeans(coefficients != 0)
names(selection_frequency) <- colnames(coefficients) <- colnames(clr_plaque_counts)
hist(selection_frequency, nclass=20)

library(pROC)
AUCs <- matrix(0, nrow=100, ncol=5)
frequency_cutoffs <- c(0.1, 0.2, 0.3, 0.4, 0.5)
variables_selected <- rep("", 5)


for (j in 1:5){
  
  marker_taxa <- names(which(selection_frequency > frequency_cutoffs[j]))
  subset_coefficient_mean <- colMeans(coefficients[, marker_taxa])
  
  
  positive_features <- which(subset_coefficient_mean > 0) |> names()
  negative_features <- which(subset_coefficient_mean < 0) |> names()
  variables_selected[j] <- sprintf("%d+%d", 
                                   length(positive_features), 
                                   length(negative_features))
  
  counts_positive_features <- rowSums(plaque_counts_imputed[, positive_features, drop=FALSE])
  counts_negative_features <- rowSums(plaque_counts_imputed[, negative_features, drop=FALSE])
  
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
marker_taxa <- names(which(selection_frequency > best_cutoff))
subset_coefficient_mean <- colMeans(coefficients[, marker_taxa])


positive_features <- which(subset_coefficient_mean > 0) |> names()
negative_features <- which(subset_coefficient_mean < 0) |> names()

predictive_features <- list(positive_features=positive_features,
                            negative_features=negative_features,
                            AUCs = performance$AUC[performance$Threshold == best_cutoff])

saveRDS(predictive_features,
        "lasso_logistic/taxa/plaque_predictive_features_yr1.rds")




