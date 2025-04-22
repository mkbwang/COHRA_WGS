# 

rm(list=ls())

saliva_relabd <- readRDS("lasso_logistic/feature_filter/saliva_relabd_imputed.rds")
plaque_relabd <- readRDS("lasso_logistic/feature_filter/plaque_relabd_imputed.rds")

saliva_counts <- read.table("counts_cleaning/strata/saliva_ko_counts_yr1.tsv",
                            header=TRUE, sep="\t", row.names=1)
plaque_counts <- read.table("counts_cleaning/strata/plaque_ko_counts_yr1.tsv",
                            header=TRUE, sep='\t', row.names=1)


saliva_KEGG_subset <- c("K00799", "K13057", "K03856", "K01625", "K07248", "K01664", "K12410")
plaque_KEGG_subset <- c("K00156", "K08659", "K15371", "K00333", "K01682", "K03342")

# adding counts
saliva_counts <- saliva_counts[, saliva_KEGG_subset]
plaque_counts <- plaque_counts[, plaque_KEGG_subset]


metadata_saliva_yr1 <- read.table("counts_cleaning/strata/metadata_saliva_yr1.tsv",
                                  header=T, sep='\t')
diagnoses_1 <- metadata_saliva_yr1$Case_status
metadata_plaque_yr2 <- read.table("counts_cleaning/strata/metadata_plaque_yr1.tsv",
                                  header=T, sep='\t')
diagnoses_2 <- metadata_plaque_yr2$Case_status

library(patchwork)

zero_counts_plot <- list()
for(j in 1:length(saliva_KEGG_subset)){
  zero_counts <- table(diagnoses_1 != 0, saliva_counts[, saliva_KEGG_subset[j]] != 0) |> as.vector()
  cases_control <- c(0,1,0,1)
  nonzero <- c(0,0,1,1)
  zero_count_df <- data.frame(Counts=zero_counts, Diagnoses = cases_control,
                              Nonzero = nonzero)
  zero_count_df$Nonzero <- factor(zero_count_df$Nonzero)
  zero_count_df$Diagnoses <- factor(zero_count_df$Diagnoses, labels=c("Control", "Case"))
  zero_counts_plot[[j]] <- ggplot(zero_count_df, aes(x=Diagnoses, y=Counts, fill=Nonzero)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values=c("#3c7de6", "#d13d78")) +
    theme_classic() + 
    xlab("Diagnoses") +
    ylab("Counts") +
    ggtitle(saliva_KEGG_subset[j])
  
}
saliva_zero_plot <- wrap_plots(zero_counts_plot, nrow=2, guides='collect')


zero_counts_plot <- list()
for(j in 1:length(plaque_KEGG_subset)){
  zero_counts <- table(diagnoses_2 != 0, plaque_counts[, plaque_KEGG_subset[j]] != 0) |> as.vector()
  cases_control <- c(0,1,0,1)
  nonzero <- c(0,0,1,1)
  zero_count_df <- data.frame(Counts=zero_counts, Diagnoses = cases_control,
                              Nonzero = nonzero)
  zero_count_df$Nonzero <- factor(zero_count_df$Nonzero)
  zero_count_df$Diagnoses <- factor(zero_count_df$Diagnoses, labels=c("Control", "Case"))
  zero_counts_plot[[j]] <- ggplot(zero_count_df, aes(x=Diagnoses, y=Counts, fill=Nonzero)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values=c("#3c7de6", "#d13d78")) +
    theme_classic() + 
    xlab("Diagnoses") +
    ylab("Counts") +
    ggtitle(plaque_KEGG_subset[j])
  
}
plaque_zero_plot <- wrap_plots(zero_counts_plot, nrow=2, guides='collect')


relabd_plots <- list()
for(j in 1:length(saliva_KEGG_subset)){
  relabd_df <- data.frame(Diagnoses = diagnoses_1,
                          logrelabd = log10(saliva_relabd[, saliva_KEGG_subset[j]]))
  
  relabd_df$Diagnoses <- factor(relabd_df$Diagnoses, labels=c("Control", "Case"))
  relabd_plots[[j]] <- ggplot(relabd_df, aes(x=Diagnoses, y=logrelabd)) +
    geom_violin() +
    theme_classic() +
    xlab("Diagnoses") +
    ylab("log10 Relative Abundance") +
    ggtitle(saliva_KEGG_subset[j])
}
saliva_relabd_plot <- wrap_plots(relabd_plots, nrow=2)




relabd_plots <- list()
for(j in 1:length(plaque_KEGG_subset)){
  relabd_df <- data.frame(Diagnoses = diagnoses_2,
                          logrelabd = log10(plaque_relabd[, plaque_KEGG_subset[j]]))
  
  relabd_df$Diagnoses <- factor(relabd_df$Diagnoses, labels=c("Control", "Case"))
  relabd_plots[[j]] <- ggplot(relabd_df, aes(x=Diagnoses, y=logrelabd)) +
    geom_violin() +
    theme_classic() +
    xlab("Diagnoses") +
    ylab("log10 Relative Abundance") +
    ggtitle(plaque_KEGG_subset[j])
}
plaque_relabd_plot <- wrap_plots(relabd_plots, nrow=2)




