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

outcome <- "Insomnia"

pheno <- pheno[,c(design_vars, outcome, covariates)]
pheno <- pheno[!duplicated(pheno$ID), ]

pheno$SOL_ID <- as.character(pheno$ID)

dat_b1 <- merge(pheno, metab_b1, by = "SOL_ID", all = TRUE)
dat_b2 <- merge(pheno, metab_b2, by = "SOL_ID", all = TRUE)

dat_b1$keep <- ifelse(!is.na(dat_b1[[exposure]]) & !is.na(dat_b1$chem_35), 1, 0)
dat_b2$keep <- ifelse(!is.na(dat_b2[[exposure]]) & !is.na(dat_b2$chem_35), 1, 0)
```

# MRS Development Using Batch 1

``` r
dat_lasso <- subset(dat_b1, keep == 1)
x <- model.matrix(as.formula(paste(outcome, "~", paste(setdiff(names(dat_lasso), c("SOL_ID", outcome)), collapse = "+"))), data = dat_lasso)

y <- data_imputation[,outcome]

n_met <- length(col_chems)
p.fac <- c(rep(0, ncol(x) - n_met), rep(1, n_met))

set.seed(123)
kcvlasso <- cv.glmnet(x, y, nfolds = 10,
                          family = "binomial",
                          type.measure = "class",
                          penalty.factor = p.fac)

ind_min <- which.min(kcvlasso$cvm)
selected_lambda <- kcvlasso$lambda[ind_min]

set.seed(123)
lasso <- cv.glmnet(x, y, lambda = selected_lambda,
                          family = "binomial",
                          type.measure = "class",
                          penalty.factor = p.fac)

lasso_beta <- lasso$beta
mrs_coefs <- lasso_beta[grep("chem", rownames(lasso_beta)), "s0"]
mrs_coefs <- mrs_coefs[which(mrs_coefs != 0)]
```

# MRS Construction in Batch 2

``` r
row.names(metab_b2) <- metab_b2$SOL_ID
tmp <- as.matrix(metab_b2[, names(mrs_coefs)]) %*% mrs_coefs
tmp <- scale(tmp)
colnames(tmp) <- paste0(outcome, "_MRS")
tmp <- data.frame(tmp)
tmp$SOL_ID <- rownames(tmp)
mrs_score <- tmp
row.names(mrs_score) <- NULL
```

# Association Analysis

\`\`\`{, eval=FALSEr} survey_base \<- svydesign( id = ~PSU_ID, strata =
~STRAT, weights = ~WEIGHT_FINAL_NORM_OVERALL, data = dat_b2 )
survey_base \<- subset(survey_base, keep == 1)

mod \<- svyglm(as.formula(paste(outcome, “~”, paste0(outcome, “\_MRS”),
“+”, paste(covariates, collapse = ” + “))), design = survey_base,
na.action = na.omit, family = family_type) \`\`\`
