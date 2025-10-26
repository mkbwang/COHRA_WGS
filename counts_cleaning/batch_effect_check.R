library(vegan)
library(ConQuR)
library(sva)
# check whether batch correction is effective

rm(list=ls())

# normalization
normalize <- function(mymat, libsize=1e5){
  output_mat <- mymat
  for (j in 1:nrow(mymat)){
    output_mat[j, ] <- mymat[j, ] / sum(mymat[j, ]) * libsize
  }
  return(output_mat)
}


## plaque
plaque_batchinfo <- read.csv("metadata/plaque_batchinfo.csv")
plaque_batchinfo$Batch <- as.factor(plaque_batchinfo$Batch)
plaque_taxa_counts <- read.csv("counts_cleaning/plaque_taxa_count_subset.csv",
                               row.names=1) |> t()
med_libsize <- median(rowSums(plaque_taxa_counts))
plaque_taxa_count_normalized <- normalize(plaque_taxa_counts, libsize=med_libsize)
plaque_taxa_counts_corrected <- read.table("counts_cleaning/plaque_taxa_count_subset_corrected.tsv",
                                 sep='\t', header=T)
plaque_ko_abundance <- read.csv("counts_cleaning/plaque_ko_abundance_subset.csv", row.names=1) |> t()
plaque_ko_abundance_normalized <- normalize(plaque_ko_abundance, libsize=5e5)
plaque_ko_abundance_corrected <- read.table("counts_cleaning/plaque_ko_abundance_subset_corrected.tsv",
                                           sep='\t', header=T)

batch_plaque <- plaque_batchinfo$Batch
names(batch_plaque) <- plaque_batchinfo$Sample
batch_plaque <- batch_plaque[rownames(plaque_taxa_counts_corrected)]


Plot_PCoA(plaque_taxa_count_normalized, factor=batch_plaque, 
          dissimilarity = "Bray")
Plot_PCoA(plaque_taxa_counts_corrected, factor=batch_plaque, 
          dissimilarity = "Bray")
Plot_PCoA(plaque_ko_abundance_normalized, factor=batch_plaque, 
          dissimilarity = "Bray")
Plot_PCoA(plaque_ko_abundance_corrected, factor=batch_plaque, 
          dissimilarity = "Bray")

# Plot_PCoA(plaque_ko_abundance, factor=batch_plaque, 
#           dissimilarity = "Aitch")
# Plot_PCoA(plaque_ko_abundance_corrected, factor=batch_plaque, 
#           dissimilarity = "Aitch")


adonis2(coda.base::dist(plaque_taxa_count_normalized+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
adonis2(coda.base::dist(plaque_taxa_counts_corrected+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
adonis2(plaque_taxa_count_normalized ~ batch_plaque, method="hellinger")
adonis2(plaque_taxa_counts_corrected ~ batch_plaque, method="hellinger")
adonis2(plaque_taxa_count_normalized ~ batch_plaque, method="jaccard")
adonis2(plaque_taxa_counts_corrected ~ batch_plaque, method="jaccard")


adonis2(coda.base::dist(plaque_ko_abundance_normalized+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
adonis2(coda.base::dist(plaque_ko_abundance_corrected+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
adonis2(plaque_ko_abundance_normalized ~ batch_plaque, method="hellinger")
adonis2(plaque_ko_abundance_corrected ~ batch_plaque, method="hellinger")
# adonis2(plaque_ko_abundance ~ batch_plaque, method="jaccard")
# adonis2(plaque_ko_abundance_corrected ~ batch_plaque, method="jaccard")

# adonis2(coda.base::dist(plaque_uniref_abundance+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
# adonis2(coda.base::dist(plaque_uniref_abundance_corrected+0.5, method="aitchison") ~ batch_plaque, method="euclidean")
# adonis2(plaque_uniref_abundance ~ batch_plaque, method="hellinger")
# adonis2(plaque_uniref_abundance_corrected ~ batch_plaque, method="hellinger")
# adonis2(plaque_uniref_abundance ~ batch_plaque, method="jaccard")
# adonis2(plaque_uniref_abundance_corrected ~ batch_plaque, method="jaccard")



## saliva
saliva_batchinfo <- read.csv("metadata/saliva_batchinfo.csv")
saliva_batchinfo$Sample_type <- as.factor(saliva_batchinfo$Sample_type)
saliva_taxa_counts <- read.csv("counts_cleaning/saliva_taxa_count_subset.csv", row.names=1) |> t()
med_libsize <- median(rowSums(saliva_taxa_counts))
saliva_taxa_count_normalized <- normalize(saliva_taxa_counts, libsize=med_libsize)
saliva_taxa_counts_corrected <- read.table("counts_cleaning/saliva_taxa_count_subset_corrected.tsv",
                                           sep='\t', header=T)

saliva_ko_abundance <- read.csv("counts_cleaning/saliva_ko_abundance_subset.csv", row.names=1) |> t()
saliva_ko_abundance_normalized <- normalize(saliva_ko_abundance, libsize=5e5)
saliva_ko_abundance_corrected <- read.table("counts_cleaning/saliva_ko_abundance_subset_corrected.tsv",
                                            sep='\t', header=T)

# saliva_uniref_abundance <- read.csv("counts_cleaning/saliva_uniref90_abundance_subset.csv",
#                                       row.names=1) |> t()
# saliva_uniref_abundance_normalized <- normalize(saliva_uniref_abundance, libsize=5e5)
# saliva_uniref_abundance_corrected <- read.table("counts_cleaning/saliva_uniref_abundance_corrected.tsv",
#                                             sep='\t', header=T)


batch_saliva <- saliva_batchinfo$Sample_type
names(batch_saliva) <- saliva_batchinfo$Sample_ID
batch_saliva <- batch_saliva[rownames(saliva_taxa_counts_corrected)]

Plot_PCoA(saliva_taxa_count_normalized, factor=batch_saliva, 
          dissimilarity = "Bray")
Plot_PCoA(saliva_taxa_counts_corrected, factor=batch_saliva, 
          dissimilarity = "Bray")
Plot_PCoA(saliva_ko_abundance_normalized, factor=batch_saliva, 
          dissimilarity = "Bray")
Plot_PCoA(saliva_ko_abundance_corrected, factor=batch_saliva, 
          dissimilarity = "Bray")
# Plot_PCoA(saliva_uniref_abundance, factor=batch_saliva, 
#           dissimilarity = "Bray")
# Plot_PCoA(saliva_uniref_abundance_corrected, factor=batch_saliva, 
#           dissimilarity = "Bray")


# Plot_PCoA(saliva_ko_abundance, factor=batch_saliva, 
#           dissimilarity = "Aitch")
# Plot_PCoA(saliva_ko_abundance_corrected, factor=batch_saliva, 
#           dissimilarity = "Aitch")


adonis2(coda.base::dist(saliva_taxa_count_normalized+0.5, method="aitchison") ~ batch_saliva, method="euclidean")
adonis2(coda.base::dist(saliva_taxa_counts_corrected+0.5, method="aitchison") ~ batch_saliva, method="euclidean")
adonis2(saliva_taxa_count_normalized ~ batch_saliva, method="hellinger")
adonis2(saliva_taxa_counts_corrected ~ batch_saliva, method="hellinger")
adonis2(saliva_taxa_count_normalized ~ batch_saliva, method="jaccard")
adonis2(saliva_taxa_counts_corrected ~ batch_saliva, method="jaccard")


adonis2(coda.base::dist(saliva_ko_abundance_normalized+0.1, method="aitchison") ~ batch_saliva, method="euclidean")
adonis2(coda.base::dist(saliva_ko_abundance_corrected+0.1, method="aitchison") ~ batch_saliva, method="euclidean")
adonis2(saliva_ko_abundance_normalized ~ batch_saliva, method="hellinger")
adonis2(saliva_ko_abundance_corrected ~ batch_saliva, method="hellinger")
adonis2(saliva_ko_abundance ~ batch_saliva, method="jaccard")
adonis2(saliva_ko_abundance_corrected ~ batch_saliva, method="jaccard")


# adonis2(coda.base::dist(saliva_uniref_abundance+0.5, method="aitchison") ~ batch_saliva, method="euclidean")
# adonis2(coda.base::dist(saliva_uniref_abundance_corrected+0.5, method="aitchison") ~ batch_saliva, method="euclidean")
# adonis2(saliva_uniref_abundance ~ batch_saliva, method="hellinger")
# adonis2(saliva_uniref_abundance_corrected ~ batch_saliva, method="hellinger")
# adonis2(saliva_uniref_abundance ~ batch_saliva, method="jaccard")
# adonis2(saliva_uniref_abundance_corrected ~ batch_saliva, method="jaccard")


