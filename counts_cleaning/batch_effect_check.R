library(vegan)
library(ConQuR)
# check whether batch correction is effective

rm(list=ls())

## plaque
plaque_batchinfo <- read.csv("metadata/plaque_batchinfo.csv")
plaque_batchinfo$Batch <- as.factor(plaque_batchinfo$Batch)
plaque_taxa_counts <- read.table("counts_cleaning/subset_plaque_taxa_count.tsv",
                                           sep='\t', header=T)
plaque_taxa_counts_corrected <- read.table("counts_cleaning/subset_plaque_taxa_count_corrected.tsv",
                                 sep='\t', header=T)
plaque_ko_abundance <- read.table("counts_cleaning/subset_plaque_ko_abundance.tsv",
                                 sep='\t', header=T)
plaque_ko_abundance_corrected <- read.table("counts_cleaning/subset_plaque_ko_abundance_corrected.tsv",
                                           sep='\t', header=T)
batch_plaque <- plaque_batchinfo$Batch
names(batch_plaque) <- plaque_batchinfo$Sample
batch_plaque <- batch_plaque[rownames(plaque_taxa_counts)]


Plot_PCoA(plaque_taxa_counts, factor=batch_plaque, 
          dissimilarity = "Bray")
Plot_PCoA(plaque_taxa_counts_corrected, factor=batch_plaque, 
          dissimilarity = "Bray")
Plot_PCoA(plaque_ko_abundance, factor=batch_plaque, 
          dissimilarity = "Bray")
Plot_PCoA(plaque_ko_abundance_corrected, factor=batch_plaque, 
          dissimilarity = "Bray")
# Plot_PCoA(plaque_ko_abundance, factor=batch_plaque, 
#           dissimilarity = "Aitch")
# Plot_PCoA(plaque_ko_abundance_corrected, factor=batch_plaque, 
#           dissimilarity = "Aitch")


adonis2(coda.base::dist(plaque_taxa_counts+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
adonis2(coda.base::dist(plaque_taxa_counts_corrected+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
adonis2(plaque_taxa_counts ~ batch_plaque, method="hellinger")
adonis2(plaque_taxa_counts_corrected ~ batch_plaque, method="hellinger")
adonis2(plaque_taxa_counts ~ batch_plaque, method="jaccard")
adonis2(plaque_taxa_counts_corrected ~ batch_plaque, method="jaccard")
adonis2(coda.base::dist(plaque_ko_abundance+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
adonis2(coda.base::dist(plaque_ko_abundance_corrected+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
adonis2(plaque_ko_abundance ~ batch_plaque, method="hellinger")
adonis2(plaque_ko_abundance_corrected ~ batch_plaque, method="hellinger")
adonis2(plaque_ko_abundance ~ batch_plaque, method="jaccard")
adonis2(plaque_ko_abundance_corrected ~ batch_plaque, method="jaccard")


## saliva
saliva_batchinfo <- read.csv("metadata/saliva_batchinfo.csv")
saliva_batchinfo$Sample_type <- as.factor(saliva_batchinfo$Sample_type)
saliva_taxa_counts <- read.table("counts_cleaning/subset_saliva_taxa_count.tsv",
                                 sep='\t', header=T)
saliva_taxa_counts_corrected <- read.table("counts_cleaning/subset_saliva_taxa_count_corrected.tsv",
                                           sep='\t', header=T)
saliva_ko_abundance <- read.table("counts_cleaning/subset_saliva_ko_abundance.tsv",
                                  sep='\t', header=T)
saliva_ko_abundance_corrected <- read.table("counts_cleaning/subset_saliva_ko_abundance_corrected.tsv",
                                            sep='\t', header=T)
batch_saliva <- saliva_batchinfo$Sample_type
names(batch_saliva) <- saliva_batchinfo$Sample_ID
batch_saliva <- batch_saliva[rownames(saliva_taxa_counts)]

Plot_PCoA(saliva_taxa_counts, factor=batch_saliva, 
          dissimilarity = "Bray")
Plot_PCoA(saliva_taxa_counts_corrected, factor=batch_saliva, 
          dissimilarity = "Bray")
Plot_PCoA(saliva_ko_abundance, factor=batch_saliva, 
          dissimilarity = "Bray")
Plot_PCoA(saliva_ko_abundance_corrected, factor=batch_saliva, 
          dissimilarity = "Bray")
# Plot_PCoA(saliva_ko_abundance, factor=batch_saliva, 
#           dissimilarity = "Aitch")
# Plot_PCoA(saliva_ko_abundance_corrected, factor=batch_saliva, 
#           dissimilarity = "Aitch")


adonis2(coda.base::dist(saliva_taxa_counts+0.5, method="aitchison") ~ batch_saliva, method="euclidean")
adonis2(coda.base::dist(saliva_taxa_counts_corrected+0.5, method="aitchison") ~ batch_saliva, method="euclidean")
adonis2(saliva_taxa_counts ~ batch_saliva, method="hellinger")
adonis2(saliva_taxa_counts_corrected ~ batch_saliva, method="hellinger")
adonis2(saliva_taxa_counts ~ batch_saliva, method="jaccard")
adonis2(saliva_taxa_counts_corrected ~ batch_saliva, method="jaccard")
adonis2(coda.base::dist(saliva_ko_abundance+0.5, method="aitchison") ~ batch_saliva, method="euclidean")
adonis2(coda.base::dist(saliva_ko_abundance_corrected+0.5, method="aitchison") ~ batch_saliva, method="euclidean")
adonis2(saliva_ko_abundance ~ batch_saliva, method="hellinger")
adonis2(saliva_ko_abundance_corrected ~ batch_saliva, method="hellinger")
adonis2(saliva_ko_abundance ~ batch_saliva, method="jaccard")
adonis2(saliva_ko_abundance_corrected ~ batch_saliva, method="jaccard")

