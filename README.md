# COHRA_WGS
Preprocessing and Analyzing whole genome sequencing of saliva and plaque samples of children in COHRA2 study.

* `plaque_preprocessing` and `saliva_preprocessing`: Use biobakery tools ([metaphlan4](https://huttenhower.sph.harvard.edu/metaphlan/) and [humann3](https://huttenhower.sph.harvard.edu/humann/)) to preprocess the fastq files and generate the microbiome taxa count and functional gene abundance (counts per million).

* `metadata`: Aggregate the demographics and batch number of samples for both plaque and saliva

* `counts_cleaning`: Aggregate the metaphlan and humann output counts into a count table. Apply batch effect correction.

* `DAA`: differential abundance analysis for both saliva and plaque with [ADAPT](https://doi.org/10.1093/bioinformatics/btae661).

* `lasso_logistic`: Lasso logistic regression on CLR transformed counts/abundances. 


