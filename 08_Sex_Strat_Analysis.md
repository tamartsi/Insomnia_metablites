``` r
library(tidyverse)
library(readxl)
library(survey)
```

# Read data

``` r
metab_file_b1 <- "Processed_metab_b1_updated.RDS"
metab_file_b2 <- "Processed_metab_b2_updated.RDS"
pheno_file <- "Processed_pheno.RDS"

metab_b1 <- readRDS(metab_file_b1)
metab_b2 <- readRDS(metab_file_b2)
pheno <- readRDS(pheno_file)

all(colnames(metab_b1) == colnames(metab_b2))
```

    ## [1] TRUE

``` r
col_chems <- setdiff(1:ncol(metab_b1), grep("SOL_ID", colnames(metab_b1)))
colnames(metab_b1)[col_chems] <- colnames(metab_b2)[col_chems] <- paste0("chem_", colnames(metab_b1)[col_chems])
```

# rank-normalize metabolite values and combine batches

``` r
rankNorm <- function(x) {
  qnorm((rank(x, ties.method = "random", na.last = "keep") - 0.5) / length(!is.na(x)))
}

chem_columns <- grep("chem", colnames(metab_b1), value = TRUE)

set.seed(89)
for (i in 1:length(chem_columns)) {
  chem_i <- chem_columns[i]
  metab_b1[, chem_i] <- rankNorm(metab_b1[, chem_i])

  metab_b2[, chem_i] <- rankNorm(metab_b2[, chem_i])
}

metab_comb <- rbind(metab_b1, metab_b2)
dat_comb <- merge(pheno, metab_comb, by = "SOL_ID", all = TRUE)
```

# create survey design object and perform regression

``` r
ass_result <- list()
design_vars <- c("ID", "STRAT", "PSU_ID", "WEIGHT_FINAL_NORM_OVERALL")
exposures <- c("WHIIRS")

for (sex in c("comb", "F", "M")) {
  for (batch in c("comb")) {
    for (exposure in exposures) {
      if (sex == "comb") {
        covariates <- c("AGE", "GENDER", "CENTER", "Background")
        dat_comb$keep <- ifelse(!is.na(dat_comb[[exposure]]) & !is.na(dat_comb$chem_35), 1, 0)
        strat <- "Sex combined"
      } else if (sex == "F") {
        covariates <- c("AGE", "CENTER", "Background")
        dat_comb$keep <- ifelse(!is.na(dat_comb[[exposure]]) & !is.na(dat_comb$chem_35) & dat_comb$GENDER == "F", 1, 0)
        strat <- "Female"
      } else {
        covariates <- c("AGE", "CENTER", "Background")
        dat_comb$keep <- ifelse(!is.na(dat_comb[[exposure]]) & !is.na(dat_comb$chem_35) & dat_comb$GENDER == "M", 1, 0)
        strat <- "Male"
      }

      survey_base <- svydesign(
        id = ~PSU_ID, strata = ~STRAT,
        weights = ~WEIGHT_FINAL_NORM_OVERALL,
        data = dat_comb
      )

      survey_base <- subset(survey_base, keep == 1)

      res_list <- vector(mode = "list", length = length(chem_columns))

      for (i in 1:length(chem_columns)) {
        chem_i <- chem_columns[i]

        model_formula_i <- paste0(chem_i, " ~ ", exposure, "+", paste(covariates, collapse = "+"))

        mod_i <- svyglm(as.formula(model_formula_i),
          design = survey_base,
          na.action = na.omit,
          family = "gaussian"
        )

        assoc_i <- summary(mod_i)$coef[exposure, ]

        res_list[[i]] <- data.frame(
          Chem_ID = chem_i,
          n = nrow(model.matrix(mod_i)),
          Est = assoc_i["Estimate"],
          SE = assoc_i["Std. Error"],
          pval = assoc_i["Pr(>|t|)"],
          Stratum = strat
        )

        res_exposures <- do.call(rbind, res_list)
        res_exposures$FDR_p <- p.adjust(res_exposures$pval, "BH")

        name <- paste0(exposure, "_", sex)

        ass_result[[name]] <- res_exposures
      }
    }
  }
}

result_tab <- do.call(rbind, ass_result)
row.names(result_tab) <- NULL
head(result_tab)
```

    ##    Chem_ID    n          Est          SE        pval      Stratum      FDR_p
    ## 1  chem_35 6016 -0.003816647 0.003285277 0.245795820 Sex combined 0.50744943
    ## 2  chem_55 6016 -0.009642438 0.003050533 0.001651025 Sex combined 0.04528525
    ## 3  chem_93 6016  0.001861773 0.003225951 0.564069461 Sex combined 0.79270261
    ## 4  chem_98 6016 -0.008309554 0.003111176 0.007767931 Sex combined 0.08904136
    ## 5 chem_111 6016 -0.007409725 0.003058292 0.015691181 Sex combined 0.11814536
    ## 6 chem_112 6016  0.004176499 0.003193098 0.191375400 Sex combined 0.46298391
