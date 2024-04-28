``` r
library(tidyverse)
library(readxl)
library(glmnet)
library(survey)
```

# Functions

## SOL Association Analysis Function

``` r
sol_association_analysis <- function(pheno, mrs_df, id_col, outcome, cov_names, family_type, exposure) {
  if (family_type == "binomial") {
    family_type <- "quasibinomial"
  }
  mrsdat <- merge(pheno, mrs_df, by = id_col, all.y = TRUE)
  mrsdat$keep <- ifelse(!is.na(mrsdat[[exposure]]), 1, 0)
  survey <- svydesign(
    id = ~PSU_ID, strata = ~STRAT,
    weights = ~WEIGHT_FINAL_NORM_OVERALL,
    data = mrsdat
  )
  survey <- subset(survey, keep == 1)
  mod <- svyglm(as.formula(paste(exposure, "~", paste0(outcome, "_", model, "_MRS"), "+", paste(cov_names, collapse = " + "))),
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

# MRS for Insomnia and WHIIRS

## Read data

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

## rank-normalize metabolite values and combine batches

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
```

## Define variables and merge phenotype with metabolites

``` r
pheno$AGE_nomed <- pheno$AGE
pheno$AGE_nomed[which(pheno$SLEA8 != 1)] <- NA

covariates <- list(
  M1 = c("AGE", "GENDER", "CENTER", "Background"),
  M1_nomed = c("AGE_nomed", "GENDER", "CENTER", "Background")
)
design_vars <- c("ID", "STRAT", "PSU_ID", "WEIGHT_FINAL_NORM_OVERALL")

dat_b1 <- merge(pheno, metab_b1, by = "SOL_ID", all.y = TRUE)

##############################
####### Outcome Define #######

outcomes <- c("Insomnia", "WHIIRS")
families <- c("binomial", "gaussian")

##############################
##############################
```

## MRS Coefficients Generation

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
    ## There are  42  metablites selected into the MRS for Insomnia_M1_nomed 
    ## There are  36  metablites selected into the MRS for WHIIRS_M1_nomed

## MRS Construction

``` r
mrs_list <- list()

for (model in names(covariates)) {
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[i]
    mrs_coefs <- coef_list[[paste(outcome, model, sep = "_")]]
    mrs <- compute_mrs(metab_comb, mrs_coefs)
    mrs[, 1] <- as.vector(scale(mrs[, 1]))
    colnames(mrs)[1] <- paste0(outcome, "_", model, "_MRS")
    mrs_list[[paste(outcome, model, sep = "_")]] <- mrs
  }
}
```

# Association with Other Sleep and Cognitive Phenotypes

``` r
pheno <- readRDS(pheno_file)
pheno$AGE_nomed <- pheno$AGE
pheno$AGE_nomed[which(pheno$SLEA8 != 1)] <- NA
pheno$Sleep_Duration <- pheno$SLPDUR
pheno$OSA <- as.factor(pheno$slpa54gt5)
pheno$SLEA8 <- as.numeric(pheno$SLEA8)
pheno$Sleep_Med_Use <- pheno$SLEA8
pheno$shift_work <- ifelse(pheno$OCEA13 %in% c("2", "3", "4", "5", "6"), 1, 0)
cog <- read.csv("20240205_global_cog_base_and_inca.csv", row.names = 1)
pheno <- merge(pheno, cog, by = "ID")
cog <- read.csv("sol_cognitive_cov_20210408.csv")
cog <- cog[c("ID", "MCI")]
pheno <- merge(pheno, cog, by = "ID")
pheno$COG <- pheno$global_cog_score
pheno$COG_Change <- pheno$global_cog_score_change

exposures <- c(
  "Sleep_Duration", "OSA", "Sleep_Med_Use", "shift_work",
  "COG", "COG_Change", "MCI"
)
families <- c(
  "gaussian", "binomial", "gaussian", "binomial",
  "gaussian", "gaussian", "binomial"
)

result_list <- list()

for (model in names(covariates)) {
  for (i in seq_along(exposures)) {
    exposure <- exposures[i]
    family_var <- families[i]
    for (outcome in outcomes) {
      mrs <- mrs_list[[paste(outcome, model, sep = "_")]]
      mod <- sol_association_analysis(
        pheno = pheno, mrs_df = mrs, id_col = "SOL_ID", outcome = outcome,
        cov_names = covariates[["M1"]], family_type = family_var, exposure = exposure
      )

      summary_mod <- summary(mod)
      estimate_first_covariate <- summary_mod$coefficients[2, "Estimate"]
      stderr_first_covariate <- summary_mod$coefficients[2, "Std. Error"]
      pvalue_first_covariate <- summary_mod$coefficients[2, "Pr(>|t|)"]
      ci_low <- estimate_first_covariate - (1.96 * stderr_first_covariate)
      ci_high <- estimate_first_covariate + (1.96 * stderr_first_covariate)
      dat_comb <- merge(metab_comb, pheno, by = "SOL_ID")
      dat_comb <- dat_comb[!is.na(dat_comb[[exposure]]), ]
      size_b2 <- nrow(dat_comb)

      if (family_var == "binomial") {
        mod_df <- data.frame(
          Exposures = exposure, MRS = outcome, Model = model,
          Sample_Size = size_b2, Est = exp(estimate_first_covariate),
          CI_Lower = exp(ci_low), CI_Upper = exp(ci_high),
          P_Value = pvalue_first_covariate, Metabolites = metabolite_num[[paste(outcome, model, sep = "_")]]
        )
      } else {
        mod_df <- data.frame(
          Exposures = exposure, MRS = outcome, Model = model,
          Sample_Size = size_b2, Est = estimate_first_covariate,
          CI_Lower = ci_low, CI_Upper = ci_high,
          P_Value = pvalue_first_covariate, Metabolites = metabolite_num[[paste(outcome, model, sep = "_")]]
        )
      }
      result_list[[paste(exposure, outcome, "MRS", model, sep = "_")]] <- mod_df
    }
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
    P_Value = sprintf("%.4e", P_Value),
    MRS = paste0(MRS, "_MRS")
  ) %>%
  select(Exposures, MRS, Model, Sample_Size, Est, CI, P_Value)
print(tab_report)
```

    ##                                           Exposures          MRS    Model
    ## Sleep_Duration_Insomnia_MRS_M1       Sleep_Duration Insomnia_MRS       M1
    ## Sleep_Duration_WHIIRS_MRS_M1         Sleep_Duration   WHIIRS_MRS       M1
    ## OSA_Insomnia_MRS_M1                             OSA Insomnia_MRS       M1
    ## OSA_WHIIRS_MRS_M1                               OSA   WHIIRS_MRS       M1
    ## Sleep_Med_Use_Insomnia_MRS_M1         Sleep_Med_Use Insomnia_MRS       M1
    ## Sleep_Med_Use_WHIIRS_MRS_M1           Sleep_Med_Use   WHIIRS_MRS       M1
    ## shift_work_Insomnia_MRS_M1               shift_work Insomnia_MRS       M1
    ## shift_work_WHIIRS_MRS_M1                 shift_work   WHIIRS_MRS       M1
    ## COG_Insomnia_MRS_M1                             COG Insomnia_MRS       M1
    ## COG_WHIIRS_MRS_M1                               COG   WHIIRS_MRS       M1
    ## COG_Change_Insomnia_MRS_M1               COG_Change Insomnia_MRS       M1
    ## COG_Change_WHIIRS_MRS_M1                 COG_Change   WHIIRS_MRS       M1
    ## MCI_Insomnia_MRS_M1                             MCI Insomnia_MRS       M1
    ## MCI_WHIIRS_MRS_M1                               MCI   WHIIRS_MRS       M1
    ## Sleep_Duration_Insomnia_MRS_M1_nomed Sleep_Duration Insomnia_MRS M1_nomed
    ## Sleep_Duration_WHIIRS_MRS_M1_nomed   Sleep_Duration   WHIIRS_MRS M1_nomed
    ## OSA_Insomnia_MRS_M1_nomed                       OSA Insomnia_MRS M1_nomed
    ## OSA_WHIIRS_MRS_M1_nomed                         OSA   WHIIRS_MRS M1_nomed
    ## Sleep_Med_Use_Insomnia_MRS_M1_nomed   Sleep_Med_Use Insomnia_MRS M1_nomed
    ## Sleep_Med_Use_WHIIRS_MRS_M1_nomed     Sleep_Med_Use   WHIIRS_MRS M1_nomed
    ## shift_work_Insomnia_MRS_M1_nomed         shift_work Insomnia_MRS M1_nomed
    ## shift_work_WHIIRS_MRS_M1_nomed           shift_work   WHIIRS_MRS M1_nomed
    ## COG_Insomnia_MRS_M1_nomed                       COG Insomnia_MRS M1_nomed
    ## COG_WHIIRS_MRS_M1_nomed                         COG   WHIIRS_MRS M1_nomed
    ## COG_Change_Insomnia_MRS_M1_nomed         COG_Change Insomnia_MRS M1_nomed
    ## COG_Change_WHIIRS_MRS_M1_nomed           COG_Change   WHIIRS_MRS M1_nomed
    ## MCI_Insomnia_MRS_M1_nomed                       MCI Insomnia_MRS M1_nomed
    ## MCI_WHIIRS_MRS_M1_nomed                         MCI   WHIIRS_MRS M1_nomed
    ##                                      Sample_Size     Est                 CI
    ## Sleep_Duration_Insomnia_MRS_M1              5893  0.0536   [0.0014, 0.1059]
    ## Sleep_Duration_WHIIRS_MRS_M1                5893  0.0865   [0.0305, 0.1425]
    ## OSA_Insomnia_MRS_M1                         5513  1.1508   [1.0579, 1.2519]
    ## OSA_WHIIRS_MRS_M1                           5513  1.1374   [1.0408, 1.2429]
    ## Sleep_Med_Use_Insomnia_MRS_M1               6107  0.1586   [0.1104, 0.2067]
    ## Sleep_Med_Use_WHIIRS_MRS_M1                 6107  0.1840   [0.1396, 0.2283]
    ## shift_work_Insomnia_MRS_M1                  6180  0.8230   [0.7474, 0.9061]
    ## shift_work_WHIIRS_MRS_M1                    6180  0.8230   [0.7486, 0.9048]
    ## COG_Insomnia_MRS_M1                         4070 -0.0867 [-0.1177, -0.0558]
    ## COG_WHIIRS_MRS_M1                           4070 -0.1149 [-0.1453, -0.0845]
    ## COG_Change_Insomnia_MRS_M1                  2936 -0.0078  [-0.0335, 0.0179]
    ## COG_Change_WHIIRS_MRS_M1                    2936 -0.0037  [-0.0292, 0.0219]
    ## MCI_Insomnia_MRS_M1                         3169  1.3370   [1.1478, 1.5574]
    ## MCI_WHIIRS_MRS_M1                           3169  1.4539   [1.2361, 1.7101]
    ## Sleep_Duration_Insomnia_MRS_M1_nomed        5893  0.0554   [0.0015, 0.1093]
    ## Sleep_Duration_WHIIRS_MRS_M1_nomed          5893  0.0830   [0.0284, 0.1376]
    ## OSA_Insomnia_MRS_M1_nomed                   5513  1.1176   [1.0267, 1.2165]
    ## OSA_WHIIRS_MRS_M1_nomed                     5513  1.1128     [1.02, 1.2141]
    ## Sleep_Med_Use_Insomnia_MRS_M1_nomed         6107  0.0787   [0.0302, 0.1272]
    ## Sleep_Med_Use_WHIIRS_MRS_M1_nomed           6107  0.1189   [0.0745, 0.1634]
    ## shift_work_Insomnia_MRS_M1_nomed            6180  0.8772   [0.8043, 0.9568]
    ## shift_work_WHIIRS_MRS_M1_nomed              6180  0.8497   [0.7707, 0.9369]
    ## COG_Insomnia_MRS_M1_nomed                   4070 -0.0837 [-0.1127, -0.0548]
    ## COG_WHIIRS_MRS_M1_nomed                     4070 -0.1066 [-0.1359, -0.0773]
    ## COG_Change_Insomnia_MRS_M1_nomed            2936 -0.0091   [-0.0332, 0.015]
    ## COG_Change_WHIIRS_MRS_M1_nomed              2936 -0.0165  [-0.0416, 0.0087]
    ## MCI_Insomnia_MRS_M1_nomed                   3169  1.3284    [1.1422, 1.545]
    ## MCI_WHIIRS_MRS_M1_nomed                     3169  1.4486   [1.2316, 1.7039]
    ##                                         P_Value
    ## Sleep_Duration_Insomnia_MRS_M1       4.4613e-02
    ## Sleep_Duration_WHIIRS_MRS_M1         2.5516e-03
    ## OSA_Insomnia_MRS_M1                  1.1353e-03
    ## OSA_WHIIRS_MRS_M1                    4.6308e-03
    ## Sleep_Med_Use_Insomnia_MRS_M1        2.2105e-10
    ## Sleep_Med_Use_WHIIRS_MRS_M1          2.4097e-15
    ## shift_work_Insomnia_MRS_M1           8.1717e-05
    ## shift_work_WHIIRS_MRS_M1             6.3170e-05
    ## COG_Insomnia_MRS_M1                  5.7582e-08
    ## COG_WHIIRS_MRS_M1                    4.2071e-13
    ## COG_Change_Insomnia_MRS_M1           5.5303e-01
    ## COG_Change_WHIIRS_MRS_M1             7.7894e-01
    ## MCI_Insomnia_MRS_M1                  2.1043e-04
    ## MCI_WHIIRS_MRS_M1                    7.5795e-06
    ## Sleep_Duration_Insomnia_MRS_M1_nomed 4.4316e-02
    ## Sleep_Duration_WHIIRS_MRS_M1_nomed   3.0003e-03
    ## OSA_Insomnia_MRS_M1_nomed            1.0428e-02
    ## OSA_WHIIRS_MRS_M1_nomed              1.6403e-02
    ## Sleep_Med_Use_Insomnia_MRS_M1_nomed  1.5412e-03
    ## Sleep_Med_Use_WHIIRS_MRS_M1_nomed    2.1446e-07
    ## shift_work_Insomnia_MRS_M1_nomed     3.2225e-03
    ## shift_work_WHIIRS_MRS_M1_nomed       1.1397e-03
    ## COG_Insomnia_MRS_M1_nomed            2.2020e-08
    ## COG_WHIIRS_MRS_M1_nomed              2.8605e-12
    ## COG_Change_Insomnia_MRS_M1_nomed     4.6071e-01
    ## COG_Change_WHIIRS_MRS_M1_nomed       1.9965e-01
    ## MCI_Insomnia_MRS_M1_nomed            2.5010e-04
    ## MCI_WHIIRS_MRS_M1_nomed              9.2741e-06
