``` r
library(tidyverse)
library(readxl)
library(glmnet)
library(survey)
```

# Functions

## SOL Association Analysis FUnction

``` r
sol_association_analysis <- function(pheno, mrs_df, id_col, outcome, cov_names, family_type = "gaussian") {
  if (family_type == "binomial") {
    family_type <- "quasibinomial"
  }
  mrsdat <- merge(pheno, mrs_df, by = id_col, all = TRUE)
  mrsdat <- na.omit(mrsdat)
  survey <- svydesign(
    id = ~PSU_ID, strata = ~STRAT,
    weights = ~WEIGHT_FINAL_NORM_OVERALL,
    data = mrsdat
  )
  mod <- svyglm(as.formula(paste(outcome, "~", paste0(outcome, "_", model, "_MRS"), "+", paste(cov_names, collapse = " + "))),
    design = survey,
    na.action = na.omit,
    family = family_type
  )
  return(mod)
}
```

## LASSO function

``` r
lasso_selected <- function(data, chem_columns, outcome, family_var, cov_names) {
  if (family_var == "binomial") {
    type_var <- "class"
  } else {
    type_var <- "mse"
  }

  x <- model.matrix(as.formula(paste(outcome, "~", paste(c(cov_names, chem_columns), collapse = "+"))), data = data)
  y <- data[, outcome]

  # apply penalty only on metabolites and not on other variables
  n_chem <- length(chem_columns)
  p.fac <- c(rep(0, ncol(x) - n_chem), rep(1, n_chem))

  # CV
  set.seed(1997)
  kcvlasso <- cv.glmnet(x, y,
    nfolds = 5,
    family = family_var,
    type.measure = type_var,
    penalty.factor = p.fac
  )

  ind_min <- which.min(kcvlasso$cvm)
  selected_lambda <- kcvlasso$lambda[ind_min]

  # lasso without cross validation
  set.seed(1997)
  lasso <- glmnet(x, y,
    lambda = selected_lambda,
    family = family_var,
    type.measure = type_var,
    penalty.factor = p.fac
  )

  lasso_beta <- lasso$beta
  mrs_coefs <- lasso_beta[grep("chem", rownames(lasso_beta)), "s0"]
  mrs_coefs <- mrs_coefs[which(mrs_coefs != 0)]

  return(mrs_coefs)
}
```

## MRS Construction Function

``` r
compute_mrs <- function(dataset, coef_set) {
  row.names(dataset) <- dataset$SOL_ID
  tmp <- as.matrix(dataset[, names(coef_set)]) %*% coef_set
  colnames(tmp) <- outcome
  tmp <- data.frame(tmp)
  tmp$SOL_ID <- rownames(tmp)
  mrs <- tmp
  row.names(mrs) <- NULL
  return(mrs)
}
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

# rank-normalize metabolite values

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
```

# Define variables and merge phenotype with metabolites

``` r
pheno$AGE_nomed <- pheno$AGE
pheno$AGE_nomed[which(pheno$SLEA8 != 1)] <- NA

covariates <- list(
  M1 = c("AGE", "GENDER", "CENTER", "Background"),
  M2 = c("AGE", "GENDER", "CENTER", "Background", "HYPERTENSION", "DIABETES2_INDICATOR"),
  M3 = c("AGE", "GENDER", "CENTER", "Background", "HYPERTENSION", "DIABETES2_INDICATOR", "AHEI2010", "CIGARETTE_USE", "ALCOHOL_USE"),
  M1_nomed = c("AGE_nomed", "GENDER", "CENTER", "Background")
)
design_vars <- c("ID", "STRAT", "PSU_ID", "WEIGHT_FINAL_NORM_OVERALL")

dat_b1 <- merge(pheno, metab_b1, by = "SOL_ID", all.y = TRUE)

##############################
####### Outcome Define #######

outcomes <- c("Insomnia", "WHIIRS", "SLEA4", "SLEA5", "SLEA6", "SLEA7")
families <- c("binomial", rep("gaussian", 5))

##############################
##############################
```

# MRS Coefficients Generation

``` r
coef_list <- list()
b1_sample_size <- list()
metabolite_num <- list()

for (model in names(covariates)) {
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[i]
    family_var <- families[i]
    dat_b1_cur <- dat_b1[, c("SOL_ID", outcome, covariates[[model]], chem_columns)]
    dat_b1_cur <- dat_b1_cur[complete.cases(dat_b1_cur), ]
    mrs_coefs <- lasso_selected(dat_b1_cur, chem_columns, outcome, family_var, cov_names = covariates[[model]])
    coef_list[[paste(outcome, model, sep = "_")]] <- mrs_coefs
    b1_sample_size[[paste(outcome, model, sep = "_")]] <- nrow(dat_b1_cur)
    cat("There are ", length(mrs_coefs), " metablites selected into the MRS for", paste(outcome, model, sep = "_"), "\n")
    metabolite_num[[paste(outcome, model, sep = "_")]] <- length(mrs_coefs)
  }
}
```

    ## There are  76  metablites selected into the MRS for Insomnia_M1 
    ## There are  68  metablites selected into the MRS for WHIIRS_M1 
    ## There are  55  metablites selected into the MRS for SLEA4_M1 
    ## There are  51  metablites selected into the MRS for SLEA5_M1 
    ## There are  47  metablites selected into the MRS for SLEA6_M1 
    ## There are  48  metablites selected into the MRS for SLEA7_M1 
    ## There are  90  metablites selected into the MRS for Insomnia_M2 
    ## There are  69  metablites selected into the MRS for WHIIRS_M2 
    ## There are  58  metablites selected into the MRS for SLEA4_M2 
    ## There are  35  metablites selected into the MRS for SLEA5_M2 
    ## There are  40  metablites selected into the MRS for SLEA6_M2 
    ## There are  47  metablites selected into the MRS for SLEA7_M2 
    ## There are  81  metablites selected into the MRS for Insomnia_M3 
    ## There are  65  metablites selected into the MRS for WHIIRS_M3 
    ## There are  52  metablites selected into the MRS for SLEA4_M3 
    ## There are  19  metablites selected into the MRS for SLEA5_M3 
    ## There are  37  metablites selected into the MRS for SLEA6_M3 
    ## There are  32  metablites selected into the MRS for SLEA7_M3 
    ## There are  42  metablites selected into the MRS for Insomnia_M1_nomed 
    ## There are  36  metablites selected into the MRS for WHIIRS_M1_nomed 
    ## There are  44  metablites selected into the MRS for SLEA4_M1_nomed 
    ## There are  31  metablites selected into the MRS for SLEA5_M1_nomed 
    ## There are  23  metablites selected into the MRS for SLEA6_M1_nomed 
    ## There are  32  metablites selected into the MRS for SLEA7_M1_nomed

# MRS Construction

``` r
mrs_list <- list()

for (model in names(covariates)) {
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[i]
    mrs_coefs <- coef_list[[paste(outcome, model, sep = "_")]]
    mrs <- compute_mrs(metab_b2, mrs_coefs)
    mrs[, 1] <- as.vector(scale(mrs[, 1]))
    colnames(mrs)[1] <- paste0(outcome, "_", model, "_MRS")
    mrs_list[[paste(outcome, model, sep = "_")]] <- mrs
  }
}
```

# Association Test

``` r
result_list <- list()

for (model in names(covariates)) {
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[i]
    family_var <- families[i]
    mrs <- mrs_list[[paste(outcome, model, sep = "_")]]
    mod <- sol_association_analysis(
      pheno = pheno, mrs_df = mrs, id_col = "SOL_ID", outcome = outcome,
      cov_names = covariates[[model]], family_type = family_var
    )

    summary_mod <- summary(mod)
    estimate_first_covariate <- summary_mod$coefficients[2, "Estimate"]
    stderr_first_covariate <- summary_mod$coefficients[2, "Std. Error"]
    pvalue_first_covariate <- summary_mod$coefficients[2, "Pr(>|t|)"]
    ci_low <- estimate_first_covariate - (1.96 * stderr_first_covariate)
    ci_high <- estimate_first_covariate + (1.96 * stderr_first_covariate)
    dat_b2 <- merge(metab_b2, pheno, by = "SOL_ID")
    if (model == "M1_nomed") {
      dat_b2 <- dat_b2[!is.na(dat_b2[[outcome]]) & !is.na(dat_b2$AGE_nomed), ]
    } else {
      dat_b2 <- dat_b2[!is.na(dat_b2[[outcome]]), ]
    }
    size_b2 <- nrow(dat_b2)

    if (family_var == "binomial") {
      mod_df <- data.frame(
        Exposures = outcome, Model = model, 
        Sample_Size = b1_sample_size[[paste(outcome, model, sep = "_")]], 
        Sample_Size_b2 = size_b2, Est = exp(estimate_first_covariate),
        CI_Lower = exp(ci_low), CI_Upper = exp(ci_high),
        P_Value = pvalue_first_covariate, Metabolites = metabolite_num[[paste(outcome, model, sep = "_")]]
      )
    } else {
      mod_df <- data.frame(
        Exposures = outcome, Model = model, 
        Sample_Size = b1_sample_size[[paste(outcome, model, sep = "_")]], 
        Sample_Size_b2 = size_b2, Est = estimate_first_covariate,
        CI_Lower = ci_low, CI_Upper = ci_high,
        P_Value = pvalue_first_covariate, Metabolites = metabolite_num[[paste(outcome, model, sep = "_")]]
      )
    }
    result_list[[paste(outcome, model, sep = "_")]] <- mod_df
  }
}
```

# Table

``` r
tab_temp <- do.call(rbind, result_list)
tab_report <- tab_temp %>%
  mutate(
    Est = round(Est, 4),
    CI = paste0("[", round(CI_Lower, 4), ", ", round(CI_Upper, 4), "]"),
    P_Value = sprintf("%.4e", P_Value)
  ) %>%
  select(Exposures, Model, Sample_Size, Sample_Size_b2, Est, CI, P_Value, Metabolites)
print(tab_report)
```

    ##                   Exposures    Model Sample_Size Sample_Size_b2    Est
    ## Insomnia_M1        Insomnia       M1        3895           2121 1.2787
    ## WHIIRS_M1            WHIIRS       M1        3895           2121 0.8114
    ## SLEA4_M1              SLEA4       M1        3954           2150 0.1677
    ## SLEA5_M1              SLEA5       M1        3953           2150 0.1641
    ## SLEA6_M1              SLEA6       M1        3954           2149 0.1570
    ## SLEA7_M1              SLEA7       M1        3905           2128 0.1156
    ## Insomnia_M2        Insomnia       M2        3895           2121 1.2302
    ## WHIIRS_M2            WHIIRS       M2        3895           2121 0.7420
    ## SLEA4_M2              SLEA4       M2        3954           2150 0.1852
    ## SLEA5_M2              SLEA5       M2        3953           2150 0.0369
    ## SLEA6_M2              SLEA6       M2        3954           2149 0.1498
    ## SLEA7_M2              SLEA7       M2        3905           2128 0.1119
    ## Insomnia_M3        Insomnia       M3        3857           2121 1.1509
    ## WHIIRS_M3            WHIIRS       M3        3857           2121 0.4949
    ## SLEA4_M3              SLEA4       M3        3915           2150 0.1009
    ## SLEA5_M3              SLEA5       M3        3914           2150 0.0071
    ## SLEA6_M3              SLEA6       M3        3916           2149 0.1108
    ## SLEA7_M3              SLEA7       M3        3867           2128 0.1015
    ## Insomnia_M1_nomed  Insomnia M1_nomed        3361           1757 1.1744
    ## WHIIRS_M1_nomed      WHIIRS M1_nomed        3361           1757 0.5939
    ## SLEA4_M1_nomed        SLEA4 M1_nomed        3413           1781 0.1064
    ## SLEA5_M1_nomed        SLEA5 M1_nomed        3412           1781 0.1396
    ## SLEA6_M1_nomed        SLEA6 M1_nomed        3412           1781 0.1479
    ## SLEA7_M1_nomed        SLEA7 M1_nomed        3369           1764 0.0736
    ##                                  CI    P_Value Metabolites
    ## Insomnia_M1        [1.0494, 1.5582] 1.5268e-02          76
    ## WHIIRS_M1          [0.2377, 1.3851] 5.8616e-03          68
    ## SLEA4_M1           [0.0108, 0.3247] 3.6927e-02          55
    ## SLEA5_M1           [0.0015, 0.3267] 4.8670e-02          51
    ## SLEA6_M1            [0.0259, 0.288] 1.9449e-02          47
    ## SLEA7_M1          [-0.0124, 0.2436] 7.7668e-02          48
    ## Insomnia_M2        [1.0084, 1.5009] 4.1878e-02          90
    ## WHIIRS_M2          [0.1163, 1.3676] 2.0676e-02          69
    ## SLEA4_M2           [0.0166, 0.3538] 3.1992e-02          58
    ## SLEA5_M2          [-0.1346, 0.2084] 6.7348e-01          35
    ## SLEA6_M2            [0.0166, 0.283] 2.8131e-02          40
    ## SLEA7_M2          [-0.0164, 0.2402] 8.8305e-02          47
    ## Insomnia_M3        [0.9527, 1.3902] 1.4580e-01          81
    ## WHIIRS_M3         [-0.0864, 1.0763] 9.6073e-02          65
    ## SLEA4_M3          [-0.0454, 0.2473] 1.7726e-01          52
    ## SLEA5_M3          [-0.1547, 0.1689] 9.3153e-01          19
    ## SLEA6_M3          [-0.0209, 0.2424] 9.9964e-02          37
    ## SLEA7_M3          [-0.0228, 0.2258] 1.1029e-01          32
    ## Insomnia_M1_nomed   [0.977, 1.4117] 8.7703e-02          42
    ## WHIIRS_M1_nomed    [0.1304, 1.0575] 1.2479e-02          36
    ## SLEA4_M1_nomed    [-0.0189, 0.2317] 9.7023e-02          44
    ## SLEA5_M1_nomed    [-0.0142, 0.2933] 7.6018e-02          31
    ## SLEA6_M1_nomed     [0.0245, 0.2713] 1.9337e-02          23
    ## SLEA7_M1_nomed    [-0.0487, 0.1958] 2.3913e-01          32
