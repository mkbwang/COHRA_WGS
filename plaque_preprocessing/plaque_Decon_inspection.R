library(dplyr)

stats_with_spikein = read.table("plaque_preprocessing/stats_with_spikein.tsv", 
                                header=T)
stats_with_spikein$file <- basename(stats_with_spikein$file)
stats_with_spikein$file <- gsub("UM1798_", "", stats_with_spikein$file)
stats_with_spikein$num_seqs <- gsub(",", "", stats_with_spikein$num_seqs) |> 
  as.integer()
fname_filter <- grepl("R1", stats_with_spikein$file) # | grepl("R2", stats_with_spikein$file)
stats_with_spikein <- stats_with_spikein[fname_filter, ]
seq_count_spikein <- stats_with_spikein %>% select(file, num_seqs) %>% rename(num_seq_with_s=num_seqs)


stats_without_spikein = read.table("plaque_preprocessing/stats_without_spikein.tsv",
                                   header=T)
stats_without_spikein$file <- basename(stats_without_spikein$file)
stats_without_spikein$num_seqs <- gsub(",", "", stats_without_spikein$num_seqs) |>
  as.integer()
fname_filter <- grepl("R1", stats_without_spikein$file) # | grepl("R2", stats_without_spikein$file)
stats_without_spikein <- stats_without_spikein[fname_filter, ]
seq_count_nospikein <- stats_without_spikein %>% select(file, num_seqs) %>% 
  rename(num_seq_no_s=num_seqs)


seq_count_combined <- seq_count_nospikein %>% left_join(seq_count_spikein, by="file")
seq_count_combined$proportion_rm <- (seq_count_combined$num_seq_with_s-seq_count_combined$num_seq_no_s) / 
  seq_count_combined$num_seq_with_s

write.table(seq_count_combined, "plaque_preprocessing/plaque_spikein_removal_details.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
