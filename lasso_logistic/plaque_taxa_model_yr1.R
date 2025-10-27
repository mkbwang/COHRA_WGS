
rm(list=ls())
library(BenchmarkDenoise)
library(ggplot2)
library(dplyr)

# load observed counts and metadata
plaque_counts <- read.table("counts_cleaning/plaque_taxa_count_subset_corrected.tsv",
                            header=TRUE, sep="\t", row.names=1) |> as.matrix()
metadata_plaque_yr1 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                                  header=T, sep='\t')
diagnoses <- (metadata_plaque_yr1$Case_status)
plaque_counts_imputed <- simple_impute(plaque_counts, scale=0.5) |> t()


# load DAA results
DAA_taxa_results <- read.csv("DAA/DAA_taxa_plaque.csv")
marker_taxa <- DAA_taxa_results$Taxa[DAA_taxa_results$pval < 0.1]


# start with DAA markers
shorten_names <- function(longname){
  strsplit(longname, split="__")[[1]] |> tail(1)
}
species_names <- sapply(colnames(plaque_counts), shorten_names) |> unname()
colnames(plaque_counts) <- colnames(plaque_counts_imputed) <- species_names

plaque_counts <- plaque_counts[, marker_taxa]
plaque_counts_imputed <- plaque_counts_imputed[, marker_taxa]

clr_plaque_counts <- clr_transform(plaque_counts_imputed)
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

auc_df <- data.frame(AUC=c(train_auc, test_auc),
                     Sample=c(rep("Train", 100), rep("Test", 100)))
auc_df$Sample <- factor(auc_df$Sample, levels=c("Train", "Test"))
AUC_plot_1 <- ggplot(auc_df, aes(x=Sample, y=AUC)) + 
  geom_boxplot() + scale_y_continuous(breaks=seq(0.5, 1, 0.05), limits=c(0.5, 1)) +
  xlab("Sample") + ylab("AUROC")

ggsave(filename="lasso_logistic/taxa/plaque_logisticlasso_raw.svg", AUC_plot_1,
       width=4, height=4)


# check variables that are most frequently selected, split into positive and negative
selection_frequency <- colMeans(coefficients != 0)
names(selection_frequency) <- colnames(coefficients) <- colnames(clr_plaque_counts)
selection_frequency_df <- data.frame(Taxa=names(selection_frequency),
                                     Frequency=selection_frequency)
selection_frequency_df <- DAA_taxa_results %>% select(Taxa, prevalence, log10foldchange) %>% 
  right_join(selection_frequency_df, by="Taxa")
selection_frequency_df$Direction <- ifelse(selection_frequency_df$log10foldchange > 0, "Enriched in Cases", "Enriched in Controls")

# look at taxa prevalence and 
prev_freq_biplot <-  ggplot(selection_frequency_df, aes(x=prevalence, y=Frequency, color=Direction)) + 
  geom_point()+ xlab("Prevalence") + ylab("LASSO Selection Probability") + 
  scale_x_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1))+
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) + 
  scale_color_manual(values=c("#8B0000", "#00008B"))+
  geom_hline(yintercept = 0.6, color = "black", linetype = "dashed", linewidth = 1)

ggsave(filename="lasso_logistic/taxa/plaque_predictive_feature.svg", prev_freq_biplot,
       width=6, height=4)


subset_taxa_df <- selection_frequency_df %>% filter(Frequency > 0.6)
subset_taxa_df_numerator <- subset_taxa_df %>% filter(Direction == "Enriched in Cases")
subset_taxa_df_denominator <- subset_taxa_df %>% filter(Direction == "Enriched in Controls")

taxa_numerator <- subset_taxa_df_numerator$Taxa
taxa_denominator <- subset_taxa_df_denominator$Taxa

taxa_count_numerator <- rowSums(plaque_counts_imputed[, taxa_numerator])
taxa_count_denominator <- rowSums(plaque_counts_imputed[, taxa_denominator])


final_df <- data.frame(CaseStatus=ifelse(metadata_plaque_yr1$Case_status == 1, "Case", "Control"),
                       LogRatio=log(taxa_count_numerator) - log(taxa_count_denominator))
final_df$Outcome <- 1*(final_df$CaseStatus == "Case")

marker_boxplot <- ggplot(final_df, aes(x=CaseStatus, y=LogRatio)) + 
  geom_boxplot() + xlab("Status") + ylab("Biomarker Value")

ggsave(filename="lasso_logistic/taxa/plaque_biomarker_boxplot.svg", marker_boxplot,
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
    title = sprintf("Plaque Taxa (AUC %.3f)", auc_value),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_bw()


ggsave(filename="lasso_logistic/taxa/plaque_AUROC.svg", roc_curve,
       width=5, height=4)

write.csv(subset_taxa_df, "lasso_logistic/taxa/plaque_taxa_biomarkers.csv", row.names=FALSE,
          quote=FALSE)

