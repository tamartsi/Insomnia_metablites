``` r
library(tidyverse)
library(survey)
library(readxl)
```

# Reading Data

``` r
exposures <- c("Insomnia", "WHIIRS")

metab_file_b1 <- "Processed_metab_b1.RDS"
metab_file_b2 <- "Processed_metab_b2.RDS"
pheno_file <- "Processed_pheno.RDS"

metab_b1 <- readRDS(metab_file_b1)
metab_b2 <- readRDS(metab_file_b2)
pheno <- readRDS(pheno_file)

col_chems <- setdiff(1:ncol(metab_b1), grep("SOL_ID", colnames(metab_b1)))
colnames(metab_b1)[col_chems] <- colnames(metab_b2)[col_chems] <- paste0("chem_", colnames(metab_b1)[col_chems])

covariates <- c("AGE", "GENDER", "CENTER", "Background")

design_vars <- c("ID", "STRAT", "PSU_ID", "WEIGHT_FINAL_NORM_OVERALL")

pheno <- pheno[,c(design_vars, exposures, covariates)]
pheno <- pheno[!duplicated(pheno$ID), ]

pheno$SOL_ID <- as.character(pheno$ID)

dat_b1 <- merge(pheno, metab_b1, by = "SOL_ID", all = TRUE)
dat_b2 <- merge(pheno, metab_b2, by = "SOL_ID", all = TRUE)

dat_b1$keep <- ifelse(!is.na(dat_b1[[exposure]]) & !is.na(dat_b1$chem_35), 1, 0)
dat_b2$keep <- ifelse(!is.na(dat_b2[[exposure]]) & !is.na(dat_b2$chem_35), 1, 0)
```

# Association Analysis Batch 1 Model 1

``` r
survey_base <- svydesign(
  id = ~PSU_ID, strata = ~STRAT,
  weights = ~WEIGHT_FINAL_NORM_OVERALL,
  data = dat_b1
)
survey_base <- subset(survey_base, keep == 1)



res_list <- vector(mode = "list", length = length(col_chems))

for (i in 1:length(col_chems)) {
  chem_i <- col_chems[i]

  model_formula_i <- paste0(chem_i, " ~ ", exposure, "+", paste(covariates, collapse = "+"))

  mod_i <- svyglm(as.formula(model_formula_i),
    design = survey_base,
    na.action = na.omit,
    family = "gaussian"
  )

  assoc_i <- summary(mod_i)$coef[exposure, ]

  res_list[[i]] <- data.frame(
    chem = chem_i,
    n = nrow(model.matrix(mod_i)),
    Est = assoc_i["Estimate"],
    SE = assoc_i["Std. Error"],
    pval = assoc_i["Pr(>|t|)"]
  )
}

res_exposures <- do.call(rbind, res_list)
res_exposures$FDR_p <- p.adjust(res_exposures$pval, "BH")
# write.csv(res_exposures, file = "association_result_m1_b1")
```

# Replication Analysis Batch 2
