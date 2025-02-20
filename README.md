# COHRA_WGS
Preprocessing and Analyzing whole genome sequencing of saliva and plaque samples of children in COHRA2 study.

* `plaque_preprocessing` and `saliva_preprocessing`: Use biobakery tools ([metaphlan4](https://huttenhower.sph.harvard.edu/metaphlan/) and [humann3](https://huttenhower.sph.harvard.edu/humann/)) to preprocess the fastq files and generate the microbiome taxa count and functional gene abundance (counts per million).

* `metadata`: Aggregate the demographics and batch number of samples for both plaque and saliva

* `counts_cleaning`: Aggregate the metaphlan and humann output counts into a count table. Apply batch effect correction.

* `Table_1`: Evaluate the correlation between demographics (e.g., income, mothers' smoking habit, education, region) and dental caries rate

* `DAA`: differential abundance analysis for both saliva and plaque with MaasLin3

* `lasso_logistic`: Lasso logistic regression on CLR transformed counts/abundances. Repeat train/test split to find out frequently selected features and average AUC.

* `network_analysis`: Among the taxa/KEGG orthologs that were selected at least once among the lasso logistic regressions, I use [SPIEC-EASI](https://github.com/zdk123/SpiecEasi) to identify highly correlated taxa/gene pairs, then use [leiden algorithm](https://www.nature.com/articles/s41598-019-41695-z) to cluster the nodes in the graph into clusters. Nodes with high degrees in each cluster and are frequently selected in lasso logistic regressions are considered "predictive" features.
