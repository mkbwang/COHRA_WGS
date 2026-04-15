
rm(list=ls())
library(BenchmarkDenoise)
library(ggplot2)
library(dplyr)
library(openxlsx)

# load observed counts and metadata
saliva_counts <- read.table("counts_cleaning/saliva_ko_abundance_subset_corrected.tsv",
                            header=TRUE, sep="\t", row.names=1) |> as.matrix()
metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
saliva_batchinfo <- read.csv("metadata/saliva_batchinfo.csv")
rownames(saliva_batchinfo) <- saliva_batchinfo$Sample_ID
saliva_batchinfo <- saliva_batchinfo[rownames(saliva_counts), ]
diagnoses <- (metadata_saliva_yr1$Case_status)
saliva_counts_imputed <- simple_impute(saliva_counts, scale=0.5) |> t()


# load DAA results
DAA_ko_results <- read.csv("DAA/DAA_ko_saliva.csv")
subset_ko <- DAA_ko_results$Taxa
marker_ko <- DAA_ko_results$Taxa[DAA_ko_results$pval < 0.05]


# start with DAA markers
saliva_counts <- saliva_counts[, subset_ko]
saliva_counts_imputed <- saliva_counts_imputed[, subset_ko]

# clr transformation before selecting marker KEGGs
clr_saliva_counts <- clr_transform(saliva_counts)
colnames(clr_saliva_counts) <- colnames(saliva_counts)
clr_saliva_counts <- clr_saliva_counts[, marker_ko]
coefficients <- matrix(0, nrow=100, ncol=ncol(clr_saliva_counts))
colnames(coefficients) <- colnames(clr_saliva_counts)
train_auc <- rep(0, 100)
test_auc <- rep(0, 100)
set.seed(2025)
for (j in 1:100){
  print(j)
  model <- fit_logistic_lasso(input=clr_saliva_counts,
                              output=diagnoses, train_prop=0.7, seed=j)
  coefficients[j, ] <- model$coefs[-1]
  train_auc[j] <- model$train_auc
  test_auc[j] <- model$test_auc
  
}

auc_df <- data.frame(AUC=c(train_auc, test_auc),
                     Sample=c(rep("Train", 100), rep("Test", 100)))
auc_df$Sample <- factor(auc_df$Sample, levels=c("Train", "Test"))
AUC_plot_1 <- ggplot(auc_df, aes(x=Sample, y=AUC)) + 
  geom_boxplot() + scale_y_continuous(breaks=seq(0.5, 1, 0.05), limits=c(0.5, 1)) +
  xlab("Sample") + ylab("AUROC")

ggsave(filename="lasso_logistic/KEGG/saliva_logisticlasso_raw.svg", AUC_plot_1,
       width=5, height=4)


# check variables that are most frequently selected, split into positive and negative
selection_frequency <- colMeans(coefficients != 0)
names(selection_frequency) <- colnames(coefficients) <- colnames(clr_saliva_counts)
selection_frequency_df <- data.frame(Taxa=names(selection_frequency),
                                     Frequency=selection_frequency)
selection_frequency_df <- DAA_ko_results %>% select(Taxa, prevalence, log10foldchange) %>% 
  right_join(selection_frequency_df, by="Taxa")
colnames(selection_frequency_df)[1] <- "KEGG"
selection_frequency_df$Direction <- ifelse(selection_frequency_df$log10foldchange > 0, "Enriched in Cases", "Enriched in Controls")

# look at ko prevalence and 
prev_freq_biplot <-  ggplot(selection_frequency_df, aes(x=prevalence, y=Frequency, color=Direction)) + 
  geom_point()+ xlab("Prevalence") + ylab("LASSO Selection Probability") + 
  scale_x_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1))+
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) + 
  scale_color_manual(values=c("#8B0000", "#00008B"))+
  geom_hline(yintercept = 0.4, color = "black", linetype = "dashed", linewidth = 1)

ggsave(filename="lasso_logistic/KEGG/saliva_predictive_feature.svg", prev_freq_biplot,
       width=6, height=4)


subset_ko_df <- selection_frequency_df %>% filter(Frequency > 0.4)
subset_ko_df_numerator <- subset_ko_df %>% filter(Direction == "Enriched in Cases")
subset_ko_df_denominator <- subset_ko_df %>% filter(Direction == "Enriched in Controls")

ko_numerator <- subset_ko_df_numerator$KEGG
ko_denominator <- subset_ko_df_denominator$KEGG

# confirm that it is not the batch effects that affect the presence absence issue of KEGGs
predictors_df <- data.frame(CaseStatus=ifelse(metadata_saliva_yr1$Case_status == 1, "Case", "Control"),
                   Batch=saliva_batchinfo$Sample_type)
predictors_df <- cbind(predictors_df, saliva_counts[, subset_ko_df$KEGG])
predictors_df_batch1 <- predictors_df %>% filter(Batch == "human saliva")
predictors_df_batch2 <- predictors_df %>% filter(Batch == "DNA")
ctable_list <- list()
for (KEGG in c(ko_numerator, ko_denominator)){
  ctable_list[[sprintf("Batch1_%s", KEGG)]] <- as.data.frame(table(predictors_df_batch1$CaseStatus, 
                                                                   predictors_df_batch1[, KEGG] > 0))
  ctable_list[[sprintf("Batch2_%s", KEGG)]] <- as.data.frame(table(predictors_df_batch2$CaseStatus, 
                                                                   predictors_df_batch2[, KEGG] > 0))
  
}
write.xlsx(ctable_list, file="lasso_logistic/KEGG/saliva_ko_presence_absence.xlsx")


# fit logistic regression
ko_count_numerator <- rowSums(saliva_counts_imputed[, ko_numerator, drop=FALSE])
ko_count_denominator <- rowSums(saliva_counts_imputed[, ko_denominator, drop=FALSE])

final_df <- data.frame(CaseStatus=ifelse(metadata_saliva_yr1$Case_status == 1, "Case", "Control"),
                       Batch=saliva_batchinfo$Sample_type,
                       LogRatio=log(ko_count_numerator) - log(ko_count_denominator))

final_df$Outcome <- 1*(final_df$CaseStatus == "Case")

marker_boxplot <- ggplot(final_df, aes(x=CaseStatus, y=LogRatio)) + 
  geom_boxplot() + xlab("Status") + ylab("Biomarker Value")


ggsave(filename="lasso_logistic/KEGG/saliva_biomarker_boxplot.svg", marker_boxplot,
       width=4, height=4)

logit_regression <-  glm(Outcome~LogRatio, data=final_df, family="binomial")

predicted_logit <- predict(logit_regression)
predicted_probability <- exp(predicted_logit) / (1 + exp(predicted_logit))

library(pROC)
roc_obj <- roc(final_df$Outcome, predicted_probability)
auc_value <- auc(roc_obj)

roc_df <- data.frame(
  Specificity = rev(roc_obj$specificities),
  Sensitivity = rev(roc_obj$sensitivities)
)

roc_curve <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity)) +
  geom_line(color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  labs(
    title = sprintf("Saliva KEGG (AUC %.3f)", auc_value),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_bw()


ggsave(filename="lasso_logistic/KEGG/saliva_AUROC.svg", roc_curve,
       width=5, height=4)

# wilcoxon test of relative abundance
libsizes <- rowSums(saliva_counts)
pvals <- rep(0, nrow(subset_ko_df))
for(j in 1:nrow(subset_ko_df)){
  
  feature_name <- subset_ko_df$KEGG[j]
  relabd <- saliva_counts[, feature_name] / libsizes
  test_result <- wilcox.test(relabd, metadata_saliva_yr1$Case_status)
  pvals[j] <- test_result$p.value
  
}

subset_ko_df$Pval <- pvals

write.csv(subset_ko_df, "lasso_logistic/KEGG/saliva_ko_biomarkers.csv", row.names=FALSE,
          quote=FALSE)


