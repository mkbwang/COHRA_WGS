
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
library(caret)
library(pROC)

expit <- function(xval) {
  exp(xval)/(1+exp(xval))
}


lambda_tune <- function(y, X, lambdas=seq(0.05, 0.55, 0.05), rep.times=50){

  
  train_auc_mean <- rep(0, length(lambdas))
  train_auc_sd <- rep(0, length(lambdas))
  test_auc_mean <- rep(0, length(lambdas))
  test_auc_sd <- rep(0, length(lambdas))
  num_features_mean <- rep(0, length(lambdas))
  num_features_sd <- rep(0, length(lambdas))
  instability <- rep(0, length(lambdas))
  choice_probs <- matrix(0, nrow=length(lambdas), 
                         ncol=ncol(X))
  
  for (j in 1:length(lambdas)){
    
    train_aucs <- rep(0, rep.times)
    test_aucs <- rep(0, rep.times)
    selected_mat <- matrix(0, nrow=rep.times, ncol=ncol(X))
    lambda <- lambdas[j]
    print(lambda)
    
    for(k in 1:rep.times){
      train_index <- createDataPartition(y, p=0.8, list=F)
      
      train_labels <- y[train_index]
      test_labels <- y[-train_index]
      train_features <- X[train_index, ]
      test_features <- X[-train_index, ]
      
      model <-  coda_logistic_lasso(y=train_labels, X=train_features, lambda=lambda)
      selected_mat[k, ] <- model$betas[-1] != 0
      
      predicted_links_train <- model$betas[1] + 
        log(as.matrix(train_features)) %*% model$betas[-1]
      predicted_probs_train <- expit(predicted_links_train)
      train_aucs[k] <- auc(roc(train_labels, as.vector(predicted_probs_train))) |> suppressMessages()
      
      predicted_links_test <- model$betas[1] + 
        log(as.matrix(test_features)) %*% model$betas[-1]
      predicted_probs_test <- expit(predicted_links_test)
      test_aucs[k] <- auc(roc(test_labels, as.vector(predicted_probs_test))) |> suppressMessages()
      
    }
    train_auc_mean[j] <- mean(train_aucs)
    train_auc_sd[j] <- sd(train_aucs)
    test_auc_mean[j] <- mean(test_aucs)
    test_auc_sd[j] <- sd(test_aucs)
    
    selection_probability <- colMeans(selected_mat)
    num_features <- rowSums(selected_mat)
    num_features_mean[j] <- mean(num_features)
    num_features_sd[j] <- sd(num_features)
    
    choice_probs[j, ] <- selection_probability
    instability[j] <- mean(selection_probability * (1-selection_probability))
    
  }
  output <- list(lambdas=lambdas,train_auc_mean=train_auc_mean, train_auc_sd=train_auc_sd,
                 test_auc_mean=test_auc_mean, test_auc_sd=test_auc_sd,
                 num_features_mean=num_features_mean, num_features_sd=num_features_sd,
                 instability=instability, choice_probs=choice_probs)
  
  return(output)
}

library(ggplot2)
library(patchwork)


codalasso_saliva <- lambda_tune(y=diagnoses_1, X=saliva_relabd)


performance_plot <- function(model){
  
  saliva_pred_performance <- data.frame(Lambdas=rep(model$lambdas, 2),
                                        AUC_mean=c(model$train_auc_mean, model$test_auc_mean),
                                        AUC_sd=c(model$train_auc_sd, model$test_auc_sd),
                                        Type=rep(c("Train", "Test"), each=length(model$lambdas)))
  
  pred_performance_plot <- ggplot(saliva_pred_performance, aes(x=Lambdas, y=AUC_mean, color=Type)) + 
    geom_point() + geom_line() +
    scale_x_continuous(breaks=model$lambdas)+ylim(0.5, 1)+
    labs(x="Lambda", y="AUC") + 
    theme_bw() + theme(legend.position="top")
  
  var_selection_df <- data.frame(Lambdas=model$lambdas,
                                 Instability = model$instability,
                                 num_features = model$num_features_mean)
  
  instability_plot <- ggplot(var_selection_df, aes(x=Lambdas, y=Instability)) +
    geom_point() + geom_line() + 
    scale_x_continuous(breaks=model$lambdas) +
    labs(x="Lambda", y="Instability") + 
    theme_bw()
  
  num_feature_plot <- ggplot(var_selection_df, aes(x=Lambdas, y=num_features))+
    geom_point() + geom_line() + 
    scale_x_continuous(breaks=model$lambdas) +
    labs(x="Lambda", y="Number of features") + 
    theme_bw()
  
  
  combined_plot <- pred_performance_plot +(instability_plot / num_feature_plot)
  
  return(combined_plot)
  
}

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
