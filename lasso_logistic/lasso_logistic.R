rm(list=ls())
library(glmnet)
library(compositions)
library(pROC)
library(doParallel)
library(foreach)
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)



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



lasso_logistic <- function(source="saliva", type="taxa", year="yr1",pseudocount=0.5,  prev=0.1,
                           feature_subset = NULL){
  
  metadata <- read.table(sprintf("counts_cleaning/strata/metadata_%s_%s.tsv", source, year),
                         header=T, sep='\t')
  
  counts <- read.table(sprintf("counts_cleaning/strata/%s_%s_counts_%s.tsv", source, type, year), 
                       header=T,
                       sep='\t') |> as.matrix()
  
  prevalences <- colMeans(counts > 0)
  counts_subset <- counts[, prevalences > prev]
  if (!is.null(feature_subset)){
    counts_subset <- counts_subset[, feature_subset]
  }
  
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



library(dplyr)
## binary visualization
library(ggplot2)
library(scales) 

performance_summary <- function(results, auc_bound=0.8, axis_title="Taxon"){
  # sort by descending test AUC
  auc_df <- data.frame(ID = seq(1, 100),
                       Train_AUC = results$train_aucs,
                       Test_AUC = results$test_aucs)
  auc_df <- auc_df %>% arrange(desc(Test_AUC))
  
  # lasso logistic regression coefficients
  coef_df <- results$coefs[, auc_df$ID] |> as.data.frame()
  
  # only keep the train/test split where the test AUC is larger than 0.8
  run_mask <- auc_df$Test_AUC >= auc_bound
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
  
  
  logistic_lasso_binary <- ggplot(coef_viz_long, aes(Run, Feature, fill = Direction))+
    geom_tile(color="gray") + 
    scale_fill_manual(values=c("#077DE8", "white", '#F04520')) +
    labs(y=axis_title)+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none",
          axis.title = element_blank())
  
  logistic_lasso_continuous <- ggplot(coef_viz_long, aes(Run, Feature, fill = Coefficient))+
    geom_tile(color="gray") + 
    scale_fill_gradient2(low="#077DE8", mid="white", high='#F04520', 
                         limits=c(-0.3, 0.3), midpoint=0, oob = squish)+
    labs(y=axis_title, fill="Coefficient")+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank())
  
  coef_viz_long_nonzeros <- coef_viz_long %>% filter(Coefficient != 0)
  nonzero_coef_violin_plot <- ggplot(coef_viz_long, aes(x=Feature, y=Coefficient)) + 
    geom_violin() + geom_boxplot(width=0.1, fill="white")+
    labs(x=axis_title, y="Coefficient") + geom_hline(yintercept=0, col="blue", linetype="dashed")+
    theme_bw() + ylim(-0.3, 0.3) + coord_flip()
  
  
  output <- list(AUCdf = auc_df, coef = coef_df,
                 
                 binary_heatmap=logistic_lasso_binary,
                 coef_heatmap=logistic_lasso_continuous,
                 nonzero_coef = nonzero_coef_violin_plot)
  
  return(output)
  
}

saliva_subset_taxa <- c('Neisseria_subflava', 'Abiotrophia_defectiva', 'Kingella_bonacorsii', 'Lachnoanaerobaculum_sp_ICM7',
                 'Actinomyces_graevenitzii', 'Candidatus_Nanosynbacter_sp_HMT_352', 'Oribacterium_asaccharolyticum',
                 'Veillonella_rogosae', 'Catonella_massiliensis', 'Prevotella_melaninogenica', 'Campylobacter_concisus',
                 'Neisseria_meningitidis', 'Ottowia_sp_Marseille_P4747', 'Veillonella_dispar', 'Streptococcus_lactarius',
                 'Prevotella_sp_HJM029')

saliva_subset_ko <- c('K00799', 'K03676', 'K03856', 'K07101', 'K09774', 'K12410',
                                   'K22292', 'K00287','K02914', 'K07248', 'K07718', 'K19265', 'K05795', 
                                   'K10843', 'K13057','K03286', 'K07662' ,'K02237', 'K00247', 'K01639',
                                   'K04760', 'K13640', 'K02392', 'K03150')


plaque_subset_taxa <- c('Neisseria_subflava', 'Neisseria_sicca', 'Streptococcus_mutans',
                        'Streptococcus_lactarius', 'Selenomonas_noxia')

plaque_subset_ko <- c('K00639', 'K14058', 'K01619', 'K02492', 'K01207',
                      'K01991', 'K03442', 'K18011', 'K01297', 'K05593')


fit_and_evaluate <- function(source="saliva", type="taxa", subset_feature=NULL){
  
  result_all <- lasso_logistic(source=source,
                               type=type,
                               year="yr1",
                               prev=0.1)
  auc_all <- data.frame(AUC=result_all$test_aucs,
                       Predictor=sprintf("%s %s (all)", source, type))
  
  performance_all <- performance_summary(results=result_all,
                                         auc_bound=0.8,
                                         axis_title="Taxon")
  
  # refit the models but start with a subset of features
  result_subset <- lasso_logistic(source=source,
                                  type=type,
                                  year="yr1",
                                  prev=0.1,
                                  feature_subset = subset_feature)
  
  auc_subset <- data.frame(AUC=result_subset$test_aucs,
                           Predictor=sprintf("%s %s (subset)", source, type))
  
  performance_subset <- performance_summary(results=result_subset,
                                            auc_bound=0.8,
                                            axis_title=source)
  
  return(list(auc_all=auc_all, auc_subset=auc_subset,
              coef_all=result_all$coef, coef_subset=result_subset$coef))
  
}


saliva_taxa <- fit_and_evaluate(source="saliva", type="taxa",
                                subset_feature=saliva_subset_taxa)

saliva_ko <- fit_and_evaluate(source="saliva", type="ko",
                              subset_feature=saliva_subset_ko)

plaque_taxa <- fit_and_evaluate(source="plaque", type="taxa",
                                subset_feature = plaque_subset_taxa)

plaque_ko <- fit_and_evaluate(source="plaque", type="ko",
                              subset_feature = plaque_subset_ko)


auc_df <- data.frame(AUC=c(saliva_taxa$auc_all$AUC, saliva_taxa$auc_subset$AUC,
                           saliva_ko$auc_all$AUC, saliva_ko$auc_subset$AUC,
                           plaque_taxa$auc_all$AUC, plaque_taxa$auc_subset$AUC,
                           plaque_ko$auc_all$AUC, plaque_ko$auc_subset$AUC),
                     Predictor = rep(c("Saliva Taxa (all)", "Saliva Taxa (subset)",
                               "Saliva KEGG (all)", "Saliva KEGG (subset)",
                               "Plaque Taxa (all)", "Plaque Taxa (subset)",
                               "Plaque KEGG (all)", "Plaque KEGG (subset)"), each=100))



# AUC_combined <- rbind(saliva_taxa_auc, saliva_ko_auc,
#                       plaque_taxa_auc, plaque_ko_auc)
# 
# 
library(ggplot2)
ggplot(auc_df, aes(x=Predictor, y=AUC)) +
  geom_boxplot() +xlab("Predictor") + ylab("Test AUC")+
  ylim(0.5, 1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


