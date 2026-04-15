# clean up metadata
library(openxlsx)
library(dplyr)


gene_counts <- read.table("transcriptome/all_samples.featureCounts.txt",
                          header=T)
gene_counts <- gene_counts[, -seq(2, 6)]
gene_counts_colnames <- colnames(gene_counts)
long_snames <- gene_counts_colnames[-1]
short_snames <- sapply(long_snames, function(longname){
  
  short_name <- strsplit(longname, split="_")[[1]][5]
  short_name <- gsub("mapped.", "", short_name)
  short_name <- gsub("[.]", "-", short_name)
  
})
colnames(gene_counts)[-1] <- short_snames

# load RNAseq metadata
metadata <- read.csv("transcriptome/DemuxStats_7792-FB.csv")
metadata$Sample_ID <- sapply(metadata$Sample_ID, function(name){
  strsplit(name, "[_]")[[1]][1]
})
metadata$IndvID <- sapply(metadata$Description, function(name){
  elements <- strsplit(name, split="-")[[1]]
  paste(elements[c(1,2)], collapse="-")
})
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[short_snames, ]

# load original ID excel sheet
alternative_IDs <- read.xlsx("metadata/COHRA_2_Readable.xlsx")
alternative_IDs <- na.omit(alternative_IDs)
subset_name <- function(myname, num_items=2){
  selected_items <- strsplit(myname, split="-")[[1]][1:num_items]
  subset_ID <- paste(selected_items, collapse="-")
  return(subset_ID)
}

alternative_IDs$readable <- sapply(alternative_IDs$readable, subset_name, num_items=2)
alternative_IDs$ID <- sapply(alternative_IDs$ID, subset_name, num_items=1)
alternative_IDs_subset <- alternative_IDs %>% filter(readable %in% metadata$IndvID) %>% 
  select(readable, ID) %>% distinct(readable, .keep_all = TRUE)
rownames(alternative_IDs_subset) <- alternative_IDs_subset$readable
alternative_IDs_subset <- alternative_IDs_subset[metadata$IndvID, ]

ids_match <- as.character(alternative_IDs_subset$ID)
names(ids_match) <- alternative_IDs_subset$readable

colnames(gene_counts)[-1] <- alternative_IDs_subset$ID

write.csv(gene_counts, file="transcriptome/gene_counts.csv", quote=FALSE)

#load microbiome metadata
metadata_microbiome <- read.xlsx("metadata/COHRA_2_metagenomic_samples_list.xlsx",
                                sheet="Project_2")

ids_intersect <- intersect(ids_match, metadata_microbiome$BabysubjectID)

rownames(gene_counts) <- gene_counts$Geneid



gene_counts$Geneid <- NULL


metadata_microbiome_subset <- metadata_microbiome %>% filter(BabysubjectID %in% ids_intersect) %>%
  filter(Visit == 5) %>%
  arrange(Case_status) %>% select(BabysubjectID, Case_status)

# median of ratios normalization


library(DESeq2)

conditions <- c(rep("Control", 3), "Case")
sample_info <- data.frame(row.names=colnames(gene_counts_subset),
                          Sample=
                          condition=conditions)
dds <- DESeqDataSetFromMatrix(countData=gene_counts_subset,
                              colData=sample_info,
                              design=~condition)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)


prevalences <- rowSums(normalized_counts > 0)

# pick genes that are present in at least two samples
normalized_counts_subset <- normalized_counts[prevalences == 4, ]
group_comparison <- rowSums(normalized_counts_subset[, c(1,2,3)] < normalized_counts_subset[, 4])



normalized_counts_enriched <- normalized_counts_subset[group_comparison == 3, ]
enriched_ratios <- normalized_counts_enriched[, 4] / rowMeans(normalized_counts_enriched[, c(1,2,3)])
enriched_df <- data.frame(Gene=rownames(normalized_counts_enriched),
                          Ratio=enriched_ratios)
enriched_df <- enriched_df %>% arrange(desc(Ratio))


normalized_counts_depleted <- normalized_counts_subset[group_comparison == 0, ]
depleted_ratios <- rowMeans(normalized_counts_depleted[, c(1,2,3)]) /  normalized_counts_depleted[, 4]
depleted_df <- data.frame(Gene=rownames(normalized_counts_depleted),
                          Ratio=depleted_ratios)
depleted_df <- depleted_df %>% arrange(desc(Ratio))

output <- list(samples=metadata_microbiome_subset,
               Case_enrich=enriched_df,
               Case_deplete=depleted_df)

write.xlsx(output, file = "transcriptome/gene_expression_summary.xlsx")

