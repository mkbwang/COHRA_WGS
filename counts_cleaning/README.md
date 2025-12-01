# README

* `counts_cleaning.R`: Shorten the taxonomy names and replace the sample names with standard naming convention (eight digit indvID-visit)
* `subset_saliva_tables.R` and `subset_plaque_tables.R` subsample the saliva and plaque count tables based on library sizes and select subset of features (taxa and KEGG) based on prevalence.
* `KEGG_uniref_filtering.R` retains KEGGs that are related to metabolism. Uniref90s are subsetted by verifying their existence in the uniprot database (5_uniref_filter.py in the plaque_preprocessing folder). Match uniref90 proteins to the KEGGs based on the dictionary provided by Chocophlan and the observed uniref90s. 
* `batch_correction_plaque.R` and `batch_correction_saliva.R`: Apply batch correction to both taxa tables and KEGG ortholog tables. `batch_effect_check.R` confirms whether batch effect is corrected.
* `stratify.R`: Split the data into 12 month samples and 24 month samples.
* `KEGG_filtering.R`: Select KEGGs that are related to metabolism.
