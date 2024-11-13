library(maaslin3)
library(dplyr)
library(ggplot2)

source("DAA/utils.R")



plaque_income_taxa_yr1 <- DAA(title="Income over $75,000 vs lower than $75,000",
                              source="plaque",
                              timestamp="yr1",
                              type="taxa",
                              covariate="inc75",
                              raw_pval = 0.05)
write.csv(plaque_income_taxa_yr1$signals, "DAA/yr1/plaque/taxa/income.csv",
          row.names=F, quote=F)
plaque_income_taxa_yr1$visual




plaque_smoker_taxa_yr1 <- DAA(title="Mother Smoker vs Nonsmoker",
                              source="plaque",
                              timestamp="yr1",
                              type="taxa",
                              covariate="Cigarettes",
                              raw_pval = 0.05)
write.csv(plaque_smoker_taxa_yr1$signals, "DAA/yr1/plaque/taxa/smoker.csv",
          row.names=F, quote=F)
plaque_smoker_taxa_yr1$visual



plaque_region_taxa_yr1 <- DAA(title="West Virginia vs Pittsburgh",
                              source="plaque",
                              timestamp="yr1",
                              type="taxa",
                              covariate="region",
                              raw_pval = 0.05)
write.csv(plaque_region_taxa_yr1$signals, "DAA/yr1/plaque/taxa/region.csv",
          row.names=F, quote=F)
plaque_region_taxa_yr1$visual



plaque_education_taxa_yr1 <- DAA(title="College educated vs non college educated mother",
                                 source="plaque",
                                 timestamp="yr1",
                                 type="taxa",
                                 covariate="higherED",
                                 raw_pval = 0.05)
write.csv(plaque_education_taxa_yr1$signals, "DAA/yr1/plaque/taxa/education.csv",
          row.names=F, quote=F)
plaque_education_taxa_yr1$visual


