rm(list=ls())
library(glmnet)
library(compositions)
library(pROC)
library(doParallel)
library(foreach)
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# clr_transformation <- function(count_vec, pseudocount=0.5){
#   count_vec[count_vec == 0] <- pseudocount
#   mean_logcount <- mean(log(count_vec))
#   clr_counts <- log(count_vec) - mean_logcount
#   return(clr_counts)
# }


fit_logistic_lasso <- function(input, output, seed=1){
  
  set.seed(seed)
  index_positive <- which(output == TRUE)
  split_mask <- rbinom(n=length(index_positive), 1, prob=0.8)
  train_positive <- index_positive[split_mask == 1]
  test_positive <- index_positive[split_mask == 0]
  
  index_negative <- which(output == FALSE)
  split_mask <- rbinom(n=length(index_negative), 1, prob=0.8)
  train_negative <- index_negative[split_mask == 1]
  test_negative <- index_negative[split_mask == 0]
  
  train_labels <- c(output[train_positive], output[train_negative]) |> as.integer()
  test_labels <- c(output[test_positive], output[test_negative]) |> as.integer()
  
  train_data <- input[c(train_positive, train_negative), ]
  test_data <- input[c(test_positive, test_negative), ]
  
  cv_fit <- cv.glmnet(x=train_data, y=train_labels, family="binomial", 
                      alpha=1, type.measure = "class")
  
  optimal_coefs <- coef(cv_fit, s = "lambda.min") |> as.vector()
  
  predicted_probabilities_train <-
    predict(cv_fit, s = "lambda.min", newx = train_data, type = "response") |> as.vector()
  
  roc_obj <- roc(train_labels, predicted_probabilities_train)
  train_auc_value <- auc(roc_obj) |> as.numeric()
  
  predicted_probabilities_test <- 
    predict(cv_fit, s = "lambda.min", newx = test_data, type = "response") |> as.vector()
  
  roc_obj <- roc(test_labels, predicted_probabilities_test)
  test_auc_value <- auc(roc_obj) |> as.numeric()
  
  output <- list(coefs=optimal_coefs, train_auc=train_auc_value, test_auc=test_auc_value)
  return(output)
  
}


lasso_logistic <- function(source="saliva", type="taxa", year="yr1",pseudocount=0.5,  prev=0.1){
  
  metadata <- read.table(sprintf("counts_cleaning/strata/metadata_%s_%s.tsv", source, year),
                         header=T, sep='\t')
  
  counts <- read.table(sprintf("counts_cleaning/strata/%s_%s_counts_%s.tsv", source, type, year), 
                       header=T,
                       sep='\t') |> as.matrix()
  
  prevalences <- colMeans(counts > 0)
  counts_subset <- counts[, prevalences > prev]
  counts_subset[counts_subset == 0] <- pseudocount
  
  clr_count_subset <- matrix(0, nrow=nrow(counts_subset), ncol=ncol(counts_subset))
  for (j in 1:nrow(counts_subset)){
    clr_count_subset[j, ] <- clr(counts_subset[j, ])
  }
  
  caries <- as.integer(metadata$Case_status == 1)
  
  logitlasso_results <- foreach(j=1:100, .packages=c("glmnet", "pROC"), .export="fit_logistic_lasso") %dopar% {
    result <- fit_logistic_lasso(input=clr_count_subset, 
                                        output=caries, seed=j)
    result
  }
  
  coefficients <- do.call(cbind, lapply(logitlasso_results, function(result) result$coefs ))
  rownames(coefficients) <- c("Intercept", colnames(counts_subset))
  train_aucs <- do.call(c, lapply(logitlasso_results, function(result) result$train_auc))
  test_aucs <- do.call(c, lapply(logitlasso_results, function(result) result$test_auc))
  
  return(list(coefs=coefficients, train_aucs=train_aucs, test_aucs=test_aucs))
  
}


saliva_taxa_yr1 <- lasso_logistic(source="saliva",
                                  type="taxa",
                                  year="yr1",
                                  prev=0.1)


library(dplyr)
# sort by descending test AUC
auc_df <- data.frame(ID = seq(1, 100),
                              Source="Saliva Taxa", 
                              Train_AUC = saliva_taxa_yr1$train_aucs,
                              Test_AUC = saliva_taxa_yr1$test_aucs)
auc_df <- auc_df %>% arrange(desc(Test_AUC))
coef_df <- saliva_taxa_yr1$coefs[, auc_df$ID] |> as.data.frame()

# only keep the train/test split where the test AUC is larger than 0.8
run_mask <- auc_df$Test_AUC >= 0.8
coef_df_subset <- as.data.frame(coef_df[-1, run_mask])

## sort the features based on the number of times they are selected
coef_df_subset$selection <- rowSums(coef_df_subset != 0)
coef_df_subset <- coef_df_subset %>% arrange(desc(selection))

## keep the features that are selected in at least half of all the train/test splits
nruns <- ncol(coef_df_subset) - 1
feature_mask <- coef_df_subset$selection >= nruns/2

coef_viz <- coef_df_subset[feature_mask, 1:nruns]
coef_viz$Feature = rownames(coef_viz)
library(reshape2)
coef_viz_long <- melt(coef_viz, id.vars="Feature", variable.name="Run", 
                      value.name="Coefficient")

coef_viz_long$Feature <- factor(coef_viz_long$Feature,
                                levels=rev(rownames(coef_viz)))
coef_viz_long$Run <- factor(coef_viz_long$Run,
                            levels=colnames(coef_viz))
coef_viz_long$Direction <- 0
coef_viz_long$Direction[coef_viz_long$Coefficient > 0] <- 1
coef_viz_long$Direction[coef_viz_long$Coefficient < 0] <- -1
coef_viz_long$Direction <- factor(coef_viz_long$Direction,
                                  levels=c(-1, 0, 1))

## binary visualization
library(ggplot2)
library(scales) 
logistic_lasso_binary <- ggplot(coef_viz_long, aes(Run, Feature, fill = Direction))+
  geom_tile(color="gray") + 
  scale_fill_manual(values=c("#077DE8", "white", '#F04520')) +
  labs(y="Taxon")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none",
        axis.title = element_blank())

logistic_lasso_continuous <- ggplot(coef_viz_long, aes(Run, Feature, fill = Coefficient))+
  geom_tile(color="gray") + 
  scale_fill_gradient2(low="#077DE8", mid="white", high='#F04520', 
                       limits=c(-0.2, 0.2), midpoint=0, oob = squish)+
  labs(y="Taxon", fill="Coefficient")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

coef_viz_long_nonzeros <- coef_viz_long %>% filter(Coefficient != 0)
nonzero_coef_violin_plot <- ggplot(coef_viz_long, aes(x=Feature, y=Coefficient)) + 
  geom_violin() + geom_boxplot(width=0.1, fill="white")+
  labs(x="Taxon", y="Coefficient") + geom_hline(yintercept=0, col="blue", linetype="dashed")+
  theme_bw() + ylim(-0.3, 0.3) + coord_flip()



occurrences <- rowSums(saliva_taxa_coefficients_subset != 0)



write.csv(saliva_taxa_yr1$coefs, "lasso_logistic/saliva_taxa_yr1.csv",
          quote=F)



saliva_ko_yr1 <- lasso_logistic(source="saliva",
                                  type="ko",
                                  year="yr1",
                                  prev=0.1)


saliva_ko_auc <- saliva_ko_auc %>% arrange(desc(Test_AUC))
saliva_ko_coefficients <- saliva_ko_yr1$coefs[, saliva_taxa_auc$ID]

good_performance_mask <- saliva_ko_auc$Test_AUC >= 0.8
saliva_ko_coefficients_subset <- as.data.frame(saliva_ko_coefficients[-1, good_performance_mask])
saliva_ko_coefficients_subset$selection <- rowSums(saliva_ko_coefficients_subset != 0)
saliva_ko_coefficients_subset <- saliva_ko_coefficients_subset %>% arrange(desc(selection))

top_coefficients <- as.matrix(saliva_ko_coefficients_subset[1:20, 1:64])




write.csv(saliva_ko_yr1$coefs, "lasso_logistic/saliva_ko_yr1.csv",
          quote=F)



saliva_ko_auc <- data.frame(Source="Saliva KEGG",
                            Train_AUC = saliva_ko_yr1$train_aucs,
                            Test_AUC = saliva_ko_yr1$test_aucs)


plaque_taxa_yr1 <- lasso_logistic(source="plaque",
                                type="taxa",
                                year="yr1",
                                prev=0.1)


write.csv(plaque_taxa_yr1$coefs, "lasso_logistic/plaque_taxa_yr1.csv",
          quote=F)


plaque_taxa_auc <- data.frame(Source="Plaque Taxa",
                            AUC=plaque_taxa_yr1$aucs)




plaque_ko_yr1 <- lasso_logistic(source="plaque",
                                  type="ko",
                                  year="yr1",
                                  prev=0.1)


write.csv(plaque_ko_yr1$coefs, "lasso_logistic/plaque_ko_yr1.csv",
          quote=F)

plaque_ko_auc <- data.frame(Source="Plaque KEGG",
                              AUC=plaque_ko_yr1$aucs)


AUC_combined <- rbind(saliva_taxa_auc, saliva_ko_auc,
                      plaque_taxa_auc, plaque_ko_auc)


library(ggplot2)
ggplot(AUC_combined, aes(x=Source, y=AUC)) + 
  geom_boxplot() +xlab("Data") + ylab("AUC")


