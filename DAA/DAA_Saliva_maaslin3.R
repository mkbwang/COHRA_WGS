library(maaslin3)
library(dplyr)
library(ggplot2)

source("DAA/utils.R")

saliva_case_taxa_yr1 <- DAA(title="Case vs control",
                            source="saliva",
                            timestamp="yr1",
                            type="taxa",
                            covariate="Case_status",
                            raw_pval = 0.1)


saliva_income_taxa_yr1 <- DAA(title="Income over $75,000 vs lower than $75,000",
               source="saliva",
               timestamp="yr1",
               type="taxa",
               covariate="inc75",
               raw_pval = 0.05)
write.csv(saliva_income_taxa_yr1$signals, "DAA/yr1/saliva/taxa/income.csv",
          row.names=F, quote=F)
saliva_income_taxa_yr1$visual




saliva_smoker_taxa_yr1 <- DAA(title="Mother Smoker vs Nonsmoker",
                              source="saliva",
                              timestamp="yr1",
                              type="taxa",
                              covariate="Cigarettes",
                              raw_pval = 0.05)
write.csv(saliva_smoker_taxa_yr1$signals, "DAA/yr1/saliva/taxa/smoker.csv",
          row.names=F, quote=F)
saliva_smoker_taxa_yr1$visual



saliva_region_taxa_yr1 <- DAA(title="West Virginia vs Pittsburgh",
                              source="saliva",
                              timestamp="yr1",
                              type="taxa",
                              covariate="region",
                              raw_pval = 0.05)
write.csv(saliva_region_taxa_yr1$signals, "DAA/yr1/saliva/taxa/region.csv",
          row.names=F, quote=F)
saliva_region_taxa_yr1$visual



saliva_education_taxa_yr1 <- DAA(title="College educated vs non college educated mother",
                              source="saliva",
                              timestamp="yr1",
                              type="taxa",
                              covariate="higherED",
                              raw_pval = 0.05)
write.csv(saliva_education_taxa_yr1$signals, "DAA/yr1/saliva/taxa/education.csv",
          row.names=F, quote=F)
saliva_education_taxa_yr1$visual


## KEGG ortholog


saliva_income_ko_yr1 <- DAA(title="Income over $75,000 vs lower than $75,000",
                              source="saliva",
                              timestamp="yr1",
                              type="ko",
                              covariate="inc75",
                              raw_pval = 0.01)
write.csv(saliva_income_ko_yr1$signals, "DAA/yr1/saliva/ko/income.csv",
          row.names=F, quote=F)
saliva_income_ko_yr1$visual




saliva_smoker_ko_yr1 <- DAA(title="Mother Smoker vs Nonsmoker",
                              source="saliva",
                              timestamp="yr1",
                              type="ko",
                              covariate="Cigarettes",
                              raw_pval = 0.01)
write.csv(saliva_smoker_ko_yr1$signals, "DAA/yr1/saliva/ko/smoker.csv",
          row.names=F, quote=F)
saliva_smoker_taxa_yr1$visual



saliva_region_ko_yr1 <- DAA(title="West Virginia vs Pittsburgh",
                              source="saliva",
                              timestamp="yr1",
                              type="taxa",
                              covariate="region",
                              raw_pval = 0.05)
write.csv(saliva_region_taxa_yr1$signals, "DAA/yr1/saliva/taxa/region.csv",
          row.names=F, quote=F)
saliva_region_taxa_yr1$visual



saliva_education_ko_yr1 <- DAA(title="College educated vs non college educated mother",
                                 source="saliva",
                                 timestamp="yr1",
                                 type="taxa",
                                 covariate="higherED",
                                 raw_pval = 0.05)
write.csv(saliva_education_ko_yr1$signals, "DAA/yr1/saliva/taxa/education.csv",
          row.names=F, quote=F)
saliva_education_taxa_yr1$visual


