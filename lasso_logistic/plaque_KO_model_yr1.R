
rm(list=ls())
library(BenchmarkDenoise)
library(dplyr)

# load observed counts and metadata
plaque_counts <- read.table("counts_cleaning/strata/plaque_ko_counts_yr1.tsv",
                            header=TRUE, sep="\t", row.names=1) |> as.matrix()
metadata_plaque_yr1 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                                  header=T, sep='\t')
diagnoses <- (metadata_plaque_yr1$Case_status)
plaque_counts_imputed <- simple_impute(plaque_counts, scale=0.5) |> t()


# load DAA results
DAA_ko_results <- read.table("DAA/yr1/plaque/ko/plaque_DA_ko_table.tsv",
                             sep='\t', header=1)
DAA_ko_results_pos <- DAA_ko_results %>% filter(log10foldchange > 0) %>% arrange(pval)
DAA_ko_results_neg <- DAA_ko_results %>% filter(log10foldchange < 0) %>% arrange(pval)

marker_KO <- DAA_ko_results$Taxa[DAA_ko_results$pval < 0.01]




# start with DAA markers
plaque_counts_imputed_0 <- plaque_counts_imputed[, marker_KO]


source("lasso_logistic/codalasso_utils.R")
lambdas_0 <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
result_0 <- repeat_codalasso(y=diagnoses, X=plaque_counts_imputed_0,
                             lambdas=lambdas_0, train_prop=0.6, times=100, ncores=6)

result_0$features <- marker_KO

train_auc_0 <- result_0$train_auc
rowMeans(train_auc_0)
test_auc_0 <- result_0$test_auc
rowMeans(test_auc_0)
coefs_list_0 <- result_0$coefs_list
best_choice_0 <- which.max(rowMeans(test_auc_0))
best_coefs_0 <- coefs_list_0[[best_choice_0]]
selection_frequency_0 <- colMeans(best_coefs_0 != 0)


# pick subset of markers which are selected over half of the times

subset_features <- which(selection_frequency_0 >= 0.5)
plaque_counts_imputed_1 <- plaque_counts_imputed_0[, subset_features]
lambdas_1 <- c(0, 0.05, 0.1, 0.15, 0.2)
result_1 <- repeat_codalasso(y=diagnoses, X=plaque_counts_imputed_1,
                             lambdas=lambdas_1, train_prop=0.6, times=100, ncores=6)
result_1$features <- colnames(plaque_counts_imputed_1)

train_auc_1 <- result_1$train_auc
rowMeans(train_auc_1)
test_auc_1 <- result_1$test_auc
rowMeans(test_auc_1)
coefs_list_1 <- result_1$coefs_list
best_choice_1 <- which.max(rowMeans(test_auc_1))
best_coefs_1 <- coefs_list_1[[best_choice_1]]
selection_frequency_1 <- colMeans(best_coefs_1 != 0)




# pick subset of markers which are selected over 80% of times

subset_features <- which(selection_frequency_1 > 0.8)
plaque_counts_imputed_2 <- plaque_counts_imputed_1[, subset_features]
lambdas_2 <- c(0, 0.05, 0.1, 0.15, 0.2)
result_2 <- repeat_codalasso(y=diagnoses, X=plaque_counts_imputed_2,
                             lambdas=lambdas_2, train_prop=0.6, times=100, ncores=6)
result_2$features <- colnames(plaque_counts_imputed_2)

train_auc_2 <- result_2$train_auc
rowMeans(train_auc_2)
test_auc_2 <- result_2$test_auc
rowMeans(test_auc_2)
coefs_list_2 <- result_2$coefs_list
best_choice_2 <- which.max(rowMeans(test_auc_2))
best_coefs_2 <- coefs_list_2[[best_choice_2]]
selection_frequency_2 <- colMeans(best_coefs_2 != 0)


output <- list(result_0, result_1, result_2)

saveRDS(output, "lasso_logistic/KEGG/plaque_yr1.rds")


# visualize the results
combined_results <- readRDS("lasso_logistic/KEGG/plaque_yr1.rds")
result_1 <- combined_results[[1]]
selected_index <- which.max(rowMeans(result_1$test_auc))
auc_df_1 <- data.frame(Train=result_1$train_auc[selected_index, ],
                       Test=result_1$test_auc[selected_index, ],
                       Iteration="First Iteration")

result_2 <- combined_results[[2]]
selected_index <- which.max(rowMeans(result_2$test_auc))
auc_df_2 <- data.frame(Train=result_2$train_auc[selected_index, ],
                       Test=result_2$test_auc[selected_index, ],
                       Iteration="Second Iteration")

result_3 <- combined_results[[3]]
selected_index <- which.max(rowMeans(result_3$test_auc))
auc_df_3 <- data.frame(Train=result_3$train_auc[selected_index, ],
                       Test=result_3$test_auc[selected_index, ],
                       Iteration="Third Iteration")

auc_df <- rbind(auc_df_1, auc_df_2, auc_df_3)
library(reshape2)

auc_df_long <- melt(auc_df, id.vars="Iteration", 
                    variable.name="Type", value.name="AUC")
  
auc_df_long$Type <- factor(auc_df_long$Type, levels=c("Train", "Test"))

ggplot(auc_df_long, aes(x=Type, y=AUC)) + 
  geom_boxplot() + xlab("") + ylab("AUC") + 
  ylim(0.5, 1)+
  facet_wrap(vars(Iteration), ncol=3) 



features_3 <- result_3$features
DAA_subset_result <- DAA_ko_results %>% filter(Taxa %in% features_3)
pos_features <- DAA_subset_result$Taxa[DAA_subset_result$log10foldchange < 0]
neg_features <- DAA_subset_result$Taxa[DAA_subset_result$log10foldchange > 0]

feature_pairs <- expand.grid(pos_features, neg_features) |> as.data.frame()
feature_pairs$Var1 <- as.character(feature_pairs$Var1)
feature_pairs$Var2 <- as.character(feature_pairs$Var2)

scatter_plots <- list()
for (j in 1:nrow(feature_pairs)){
  
  feature_x <- feature_pairs$Var1[j]
  feature_y <- feature_pairs$Var2[j]
  counts_subset <- plaque_counts_imputed[, c(feature_x, feature_y)] |>
    as.data.frame()
  counts_subset$Caries <- as.factor(diagnoses)
  
  scatter_plots[[j]] <- ggplot(counts_subset, aes(x=.data[[feature_x]], 
                            y=.data[[feature_y]], color=Caries)) + 
    geom_point(size=1, alpha=0.7) + 
    scale_x_log10() + scale_y_log10()+
    xlab(feature_x) +
    ylab(feature_y)
  
}

library(patchwork)

all_plots <- wrap_plots(scatter_plots, nrow=length(neg_features)) +
  plot_layout(guides="collect")



