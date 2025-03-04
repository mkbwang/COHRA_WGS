# README

* `feature_filtering.R`: only select KEGGs that are related to metabolism and are differentially abundant

* `zero_imputation.R`: Impute the zeros with `zComposition` package

* `lasso_logcontrast.R`: LASSO log contrast model for the filtered features

* `predictive_features.R`: Further select a subset of predictive features

* `lasso_logcontrast_final.R`: Fit a final model with the selected predictive features.

* `functions_coda_penalized_regression.R`: Codes for log contrast model based on this [publication](https://doi.org/10.1093/nargab/lqaa029).

* `functions_lasso_tune.R`: Tuning penalty parameters.
