# README

* `counts_cleaning.R`: Shorten the taxonomy names and replace the sample names with standard naming convention (eight digit indvID-visit)
* `batch_correction_plaque.R` and `batch_correction_saliva.R`: Retain samples whose library sizes based on taxa table are at least 5e5, then apply batch correction to both taxa tables and KEGG ortholog tables. `batch_effect_check.R` confirms whether batch effect is corrected.
* `stratify.R`: Split the data into 12 month samples and 24 month samples.
* `KEGG_filtering.R`: Select KEGGs that are related to metabolism.
