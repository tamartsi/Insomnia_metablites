``` r
library(tidyverse)
library(survey)
library(readxl)
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
chem_columns <- grep("chem", colnames(metab_b1), value = TRUE)

rankNorm <- function(x) {
  qnorm((rank(x, ties.method = "random", na.last = "keep") - 0.5) / length(!is.na(x)))
}

set.seed(89)
for (i in 1:length(chem_columns)) {
  chem_i <- chem_columns[i]
  metab_b1[, chem_i] <- rankNorm(metab_b1[, chem_i])
  metab_b2[, chem_i] <- rankNorm(metab_b2[, chem_i])
}

metab_comb <- rbind(metab_b1, metab_b2)
```

# Merge phenotypes with metabolites

``` r
pheno$Sleep_Duration <- pheno$SLPDUR
pheno$OSA <- ifelse(pheno$slpa54gt5 == 1, 1, 0)
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

exposures <- c("Insomnia", "WHIIRS", "SLEA4", "SLEA5", "SLEA6", "SLEA7",
               "Sleep_Duration", "OSA", "Sleep_Med_Use", "shift_work",
               "COG", "COG_Change", "MCI")

pheno$AGE_nomed <- pheno$AGE
pheno$AGE_nomed[which(pheno$SLEA8 != 1)] <- NA

covariates <- list(
  M1 = c("AGE", "GENDER", "CENTER", "Background"),
  M2 = c("AGE", "GENDER", "CENTER", "Background", "HYPERTENSION", "DIABETES2_INDICATOR"),
  M3 = c("AGE", "GENDER", "CENTER", "Background", "HYPERTENSION", "DIABETES2_INDICATOR", "AHEI2010", "CIGARETTE_USE", "ALCOHOL_USE"),
  M1_nomed = c("AGE_nomed", "GENDER", "CENTER", "Background")
)

design_vars <- c("ID", "STRAT", "PSU_ID", "WEIGHT_FINAL_NORM_OVERALL")

dat_comb <- merge(pheno, metab_comb, by = "SOL_ID", all = TRUE)
```

# Association Analysis

``` r
ass_result <- list()
reped_chem <- c("chem_100000010", "chem_100001083")

for (batch in c("comb")) {
  for (exposure in exposures) {
    dat_comb$keep <- ifelse(!is.na(dat_comb[[exposure]]) & !is.na(dat_comb$chem_35), 1, 0)
    if (batch == "comb") {
      survey_base <- svydesign(
        id = ~PSU_ID, strata = ~STRAT,
        weights = ~WEIGHT_FINAL_NORM_OVERALL,
        data = dat_comb
      )
    }

    survey_base <- subset(survey_base, keep == 1)

    for (model in c("M1", "M2", "M3", "M1_nomed")) {
      if (exposure == "Sleep_Med_Use" & model == "M1_nomed") {
        next
      }
      res_list <- vector(mode = "list", length = length(chem_columns))

      for (i in 1:length(reped_chem)) {
        chem_i <- reped_chem[i]

        model_formula_i <- paste0(chem_i, " ~ ", exposure, "+", paste(covariates[[model]], collapse = "+"))

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
          pval = assoc_i["Pr(>|t|)"]
        )
      }

      res_exposures <- do.call(rbind, res_list)
      res_exposures$FDR_p <- p.adjust(res_exposures$pval, "BH")

      name <- paste0(exposure, "_", batch, "_", model)

      ass_result[[name]] <- res_exposures
    }
  }
}
```

# Table

``` r
combined_df <- do.call(rbind, lapply(names(ass_result), function(name) {
  df <- ass_result[[name]]
  df$Model <- name
  return(df)
}))

rownames(combined_df) <- NULL

metab_batch2_file <- file.path("batch2_data.xlsx")
metabolites_info_sheet_name <- "metabolites.info"

metab_info <- as.data.frame(read_excel(metab_batch2_file, sheet = metabolites_info_sheet_name))
metab_info$CHEM_ID <- paste0("chem_", metab_info$CHEM_ID)
metab_info <- metab_info[,c("CHEM_ID","CHEMICAL_NAME", "HMDB", "KEGG")]
names(metab_info)[1] <- "Chem_ID"

merged_results <- lapply(ass_result, function(df) merge(df, metab_info, by = "Chem_ID"))

combined_df <- do.call(rbind, lapply(names(merged_results), function(name) {
  df <- merged_results[[name]]
  df$Model <- name
  return(df)
}))


combined_df <- subset(combined_df, select = -Chem_ID)
combined_df <- combined_df[c(6:9,1:5)]
print(combined_df)
```

    ##                           CHEMICAL_NAME        HMDB   KEGG
    ## 1   3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 2                      indolepropionate HMDB0002302   <NA>
    ## 3   3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 4                      indolepropionate HMDB0002302   <NA>
    ## 5   3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 6                      indolepropionate HMDB0002302   <NA>
    ## 7   3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 8                      indolepropionate HMDB0002302   <NA>
    ## 9   3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 10                     indolepropionate HMDB0002302   <NA>
    ## 11  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 12                     indolepropionate HMDB0002302   <NA>
    ## 13  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 14                     indolepropionate HMDB0002302   <NA>
    ## 15  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 16                     indolepropionate HMDB0002302   <NA>
    ## 17  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 18                     indolepropionate HMDB0002302   <NA>
    ## 19  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 20                     indolepropionate HMDB0002302   <NA>
    ## 21  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 22                     indolepropionate HMDB0002302   <NA>
    ## 23  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 24                     indolepropionate HMDB0002302   <NA>
    ## 25  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 26                     indolepropionate HMDB0002302   <NA>
    ## 27  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 28                     indolepropionate HMDB0002302   <NA>
    ## 29  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 30                     indolepropionate HMDB0002302   <NA>
    ## 31  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 32                     indolepropionate HMDB0002302   <NA>
    ## 33  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 34                     indolepropionate HMDB0002302   <NA>
    ## 35  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 36                     indolepropionate HMDB0002302   <NA>
    ## 37  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 38                     indolepropionate HMDB0002302   <NA>
    ## 39  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 40                     indolepropionate HMDB0002302   <NA>
    ## 41  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 42                     indolepropionate HMDB0002302   <NA>
    ## 43  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 44                     indolepropionate HMDB0002302   <NA>
    ## 45  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 46                     indolepropionate HMDB0002302   <NA>
    ## 47  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 48                     indolepropionate HMDB0002302   <NA>
    ## 49  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 50                     indolepropionate HMDB0002302   <NA>
    ## 51  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 52                     indolepropionate HMDB0002302   <NA>
    ## 53  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 54                     indolepropionate HMDB0002302   <NA>
    ## 55  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 56                     indolepropionate HMDB0002302   <NA>
    ## 57  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 58                     indolepropionate HMDB0002302   <NA>
    ## 59  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 60                     indolepropionate HMDB0002302   <NA>
    ## 61  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 62                     indolepropionate HMDB0002302   <NA>
    ## 63  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 64                     indolepropionate HMDB0002302   <NA>
    ## 65  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 66                     indolepropionate HMDB0002302   <NA>
    ## 67  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 68                     indolepropionate HMDB0002302   <NA>
    ## 69  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 70                     indolepropionate HMDB0002302   <NA>
    ## 71  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 72                     indolepropionate HMDB0002302   <NA>
    ## 73  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 74                     indolepropionate HMDB0002302   <NA>
    ## 75  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 76                     indolepropionate HMDB0002302   <NA>
    ## 77  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 78                     indolepropionate HMDB0002302   <NA>
    ## 79  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 80                     indolepropionate HMDB0002302   <NA>
    ## 81  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 82                     indolepropionate HMDB0002302   <NA>
    ## 83  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 84                     indolepropionate HMDB0002302   <NA>
    ## 85  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 86                     indolepropionate HMDB0002302   <NA>
    ## 87  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 88                     indolepropionate HMDB0002302   <NA>
    ## 89  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 90                     indolepropionate HMDB0002302   <NA>
    ## 91  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 92                     indolepropionate HMDB0002302   <NA>
    ## 93  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 94                     indolepropionate HMDB0002302   <NA>
    ## 95  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 96                     indolepropionate HMDB0002302   <NA>
    ## 97  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 98                     indolepropionate HMDB0002302   <NA>
    ## 99  3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 100                    indolepropionate HMDB0002302   <NA>
    ## 101 3-phenylpropionate (hydrocinnamate) HMDB0000764 C05629
    ## 102                    indolepropionate HMDB0002302   <NA>
    ##                            Model    n           Est          SE         pval
    ## 1               Insomnia_comb_M1 6016 -0.1235040780 0.035451341 5.299543e-04
    ## 2               Insomnia_comb_M1 6016 -0.1401314613 0.042306608 9.802613e-04
    ## 3               Insomnia_comb_M2 6016 -0.1183079174 0.035170785 8.172425e-04
    ## 4               Insomnia_comb_M2 6016 -0.1345046662 0.041462575 1.243461e-03
    ## 5               Insomnia_comb_M3 5949 -0.1006999716 0.035354406 4.545212e-03
    ## 6               Insomnia_comb_M3 5949 -0.1118505152 0.043388024 1.017611e-02
    ## 7         Insomnia_comb_M1_nomed 5118 -0.1089210463 0.039457434 5.950669e-03
    ## 8         Insomnia_comb_M1_nomed 5118 -0.1142717910 0.045322630 1.195376e-02
    ## 9                 WHIIRS_comb_M1 6016 -0.0150135899 0.003242125 4.457271e-06
    ## 10                WHIIRS_comb_M1 6016 -0.0157319328 0.003335146 2.971537e-06
    ## 11                WHIIRS_comb_M2 6016 -0.0145739501 0.003202620 6.464303e-06
    ## 12                WHIIRS_comb_M2 6016 -0.0152363317 0.003308785 5.033196e-06
    ## 13                WHIIRS_comb_M3 5949 -0.0124012447 0.003241979 1.442703e-04
    ## 14                WHIIRS_comb_M3 5949 -0.0123954415 0.003345260 2.304926e-04
    ## 15          WHIIRS_comb_M1_nomed 5118 -0.0138180339 0.003563587 1.173233e-04
    ## 16          WHIIRS_comb_M1_nomed 5118 -0.0142774383 0.003722277 1.386569e-04
    ## 17                 SLEA4_comb_M1 6104 -0.0460991166 0.011165267 4.156464e-05
    ## 18                 SLEA4_comb_M1 6104 -0.0460782302 0.011796073 1.042584e-04
    ## 19                 SLEA4_comb_M2 6104 -0.0450699740 0.011023256 4.924705e-05
    ## 20                 SLEA4_comb_M2 6104 -0.0449744784 0.011694226 1.328386e-04
    ## 21                 SLEA4_comb_M3 6036 -0.0393819757 0.011085535 4.113200e-04
    ## 22                 SLEA4_comb_M3 6036 -0.0384025996 0.011598565 9.852246e-04
    ## 23           SLEA4_comb_M1_nomed 5194 -0.0383983533 0.012860372 2.945353e-03
    ## 24           SLEA4_comb_M1_nomed 5194 -0.0350578400 0.013901710 1.193562e-02
    ## 25                 SLEA5_comb_M1 6103 -0.0355843442 0.011431900 1.940731e-03
    ## 26                 SLEA5_comb_M1 6103 -0.0524298502 0.011148922 3.180084e-06
    ## 27                 SLEA5_comb_M2 6103 -0.0333648597 0.011376826 3.487159e-03
    ## 28                 SLEA5_comb_M2 6103 -0.0495646245 0.011090308 9.371986e-06
    ## 29                 SLEA5_comb_M3 6035 -0.0273134614 0.011318487 1.611152e-02
    ## 30                 SLEA5_comb_M3 6035 -0.0422109546 0.011004298 1.383006e-04
    ## 31           SLEA5_comb_M1_nomed 5193 -0.0270878846 0.012334314 2.846901e-02
    ## 32           SLEA5_comb_M1_nomed 5193 -0.0434645623 0.011847176 2.657002e-04
    ## 33                 SLEA6_comb_M1 6103 -0.0310782887 0.012558993 1.361083e-02
    ## 34                 SLEA6_comb_M1 6103 -0.0338288515 0.012089008 5.299827e-03
    ## 35                 SLEA6_comb_M2 6103 -0.0305585205 0.012513556 1.488959e-02
    ## 36                 SLEA6_comb_M2 6103 -0.0331713488 0.012028975 5.997921e-03
    ## 37                 SLEA6_comb_M3 6036 -0.0241834393 0.012636409 5.612004e-02
    ## 38                 SLEA6_comb_M3 6036 -0.0241652246 0.011894027 4.262086e-02
    ## 39           SLEA6_comb_M1_nomed 5193 -0.0288723900 0.013319692 3.058339e-02
    ## 40           SLEA6_comb_M1_nomed 5193 -0.0354788970 0.012811988 5.795620e-03
    ## 41                 SLEA7_comb_M1 6033 -0.0420764506 0.012670169 9.509965e-04
    ## 42                 SLEA7_comb_M1 6033 -0.0342689155 0.013837234 1.353602e-02
    ## 43                 SLEA7_comb_M2 6033 -0.0416063196 0.012509307 9.341782e-04
    ## 44                 SLEA7_comb_M2 6033 -0.0339439048 0.013610612 1.289919e-02
    ## 45                 SLEA7_comb_M3 5966 -0.0367216847 0.012323713 3.000890e-03
    ## 46                 SLEA7_comb_M3 5966 -0.0275021990 0.013890053 4.815856e-02
    ## 47           SLEA7_comb_M1_nomed 5133 -0.0412820214 0.013773084 2.838019e-03
    ## 48           SLEA7_comb_M1_nomed 5133 -0.0368118822 0.015307570 1.648661e-02
    ## 49        Sleep_Duration_comb_M1 5893  0.0100166725 0.012279616 4.149834e-01
    ## 50        Sleep_Duration_comb_M1 5893  0.0100188099 0.012804995 4.342768e-01
    ## 51        Sleep_Duration_comb_M2 5893  0.0125280444 0.012286082 3.082818e-01
    ## 52        Sleep_Duration_comb_M2 5893  0.0131368057 0.012789538 3.047598e-01
    ## 53        Sleep_Duration_comb_M3 5828  0.0117990813 0.011697324 3.135255e-01
    ## 54        Sleep_Duration_comb_M3 5828  0.0127856906 0.011908678 2.834135e-01
    ## 55  Sleep_Duration_comb_M1_nomed 5031  0.0118148450 0.014164214 4.045439e-01
    ## 56  Sleep_Duration_comb_M1_nomed 5031  0.0129414176 0.014055388 3.575609e-01
    ## 57                   OSA_comb_M1 5513 -0.2592453976 0.037505697 1.223703e-11
    ## 58                   OSA_comb_M1 5513 -0.2082206361 0.040130381 2.902899e-07
    ## 59                   OSA_comb_M2 5513 -0.2367952828 0.036462419 1.756503e-10
    ## 60                   OSA_comb_M2 5513 -0.1852934506 0.039570382 3.507716e-06
    ## 61                   OSA_comb_M3 5453 -0.2209554534 0.036795870 3.332574e-09
    ## 62                   OSA_comb_M3 5453 -0.1705006257 0.038459341 1.105392e-05
    ## 63             OSA_comb_M1_nomed 4731 -0.2977367746 0.041244482 1.630662e-12
    ## 64             OSA_comb_M1_nomed 4731 -0.2297768967 0.042314106 8.249234e-08
    ## 65         Sleep_Med_Use_comb_M1 6107 -0.0638082063 0.017686852 3.343706e-04
    ## 66         Sleep_Med_Use_comb_M1 6107 -0.0520618894 0.016701962 1.912523e-03
    ## 67         Sleep_Med_Use_comb_M2 6107 -0.0627511023 0.017334144 3.191259e-04
    ## 68         Sleep_Med_Use_comb_M2 6107 -0.0500067629 0.016443660 2.459003e-03
    ## 69         Sleep_Med_Use_comb_M3 6039 -0.0607744269 0.017387564 5.083266e-04
    ## 70         Sleep_Med_Use_comb_M3 6039 -0.0434377796 0.016228347 7.638246e-03
    ## 71            shift_work_comb_M1 6180  0.1430568055 0.047014859 2.445004e-03
    ## 72            shift_work_comb_M1 6180  0.1128487397 0.043253643 9.304115e-03
    ## 73            shift_work_comb_M2 6180  0.1365602240 0.046875710 3.708816e-03
    ## 74            shift_work_comb_M2 6180  0.1064408942 0.043133118 1.387244e-02
    ## 75            shift_work_comb_M3 6104  0.1367333936 0.046370087 3.314285e-03
    ## 76            shift_work_comb_M3 6104  0.1020140049 0.041734860 1.479727e-02
    ## 77      shift_work_comb_M1_nomed 5269  0.1331152451 0.047440474 5.181386e-03
    ## 78      shift_work_comb_M1_nomed 5269  0.1032647437 0.046297955 2.609107e-02
    ## 79                   COG_comb_M1 4070  0.0469086581 0.031357587 1.352101e-01
    ## 80                   COG_comb_M1 4070  0.0661047929 0.026920833 1.435669e-02
    ## 81                   COG_comb_M2 4070  0.0351480813 0.031132712 2.593707e-01
    ## 82                   COG_comb_M2 4070  0.0482214095 0.026855081 7.307158e-02
    ## 83                   COG_comb_M3 4027  0.0232801765 0.031207434 4.559807e-01
    ## 84                   COG_comb_M3 4027  0.0346846770 0.027270436 2.039261e-01
    ## 85             COG_comb_M1_nomed 3377  0.0320999196 0.034715183 3.555342e-01
    ## 86             COG_comb_M1_nomed 3377  0.0604601139 0.032117394 6.028458e-02
    ## 87            COG_Change_comb_M1 2936 -0.0168036303 0.052589969 7.494539e-01
    ## 88            COG_Change_comb_M1 2936 -0.0564786930 0.040626506 1.650398e-01
    ## 89            COG_Change_comb_M2 2936 -0.0137234652 0.052142946 7.925049e-01
    ## 90            COG_Change_comb_M2 2936 -0.0572263743 0.040066861 1.537906e-01
    ## 91            COG_Change_comb_M3 2910 -0.0168795917 0.052251700 7.467878e-01
    ## 92            COG_Change_comb_M3 2910 -0.0559015639 0.040477215 1.678342e-01
    ## 93      COG_Change_comb_M1_nomed 2457  0.0009488888 0.054478272 9.861100e-01
    ## 94      COG_Change_comb_M1_nomed 2457 -0.0579473132 0.044745472 1.958801e-01
    ## 95                   MCI_comb_M1 3169 -0.0264310351 0.078756734 7.372972e-01
    ## 96                   MCI_comb_M1 3169 -0.0860089311 0.064364363 1.820067e-01
    ## 97                   MCI_comb_M2 3169 -0.0008998430 0.078460971 9.908537e-01
    ## 98                   MCI_comb_M2 3169 -0.0455251313 0.062515852 4.667899e-01
    ## 99                   MCI_comb_M3 3135  0.0231218326 0.078799769 7.693085e-01
    ## 100                  MCI_comb_M3 3135 -0.0167672916 0.064673515 7.955317e-01
    ## 101            MCI_comb_M1_nomed 2636  0.0129011919 0.090419225 8.865959e-01
    ## 102            MCI_comb_M1_nomed 2636 -0.0465810161 0.068114867 4.943643e-01
    ##            FDR_p
    ## 1   9.802613e-04
    ## 2   9.802613e-04
    ## 3   1.243461e-03
    ## 4   1.243461e-03
    ## 5   9.090424e-03
    ## 6   1.017611e-02
    ## 7   1.190134e-02
    ## 8   1.195376e-02
    ## 9   4.457271e-06
    ## 10  4.457271e-06
    ## 11  6.464303e-06
    ## 12  6.464303e-06
    ## 13  2.304926e-04
    ## 14  2.304926e-04
    ## 15  1.386569e-04
    ## 16  1.386569e-04
    ## 17  8.312927e-05
    ## 18  1.042584e-04
    ## 19  9.849410e-05
    ## 20  1.328386e-04
    ## 21  8.226400e-04
    ## 22  9.852246e-04
    ## 23  5.890706e-03
    ## 24  1.193562e-02
    ## 25  1.940731e-03
    ## 26  6.360167e-06
    ## 27  3.487159e-03
    ## 28  1.874397e-05
    ## 29  1.611152e-02
    ## 30  2.766011e-04
    ## 31  2.846901e-02
    ## 32  5.314004e-04
    ## 33  1.361083e-02
    ## 34  1.059965e-02
    ## 35  1.488959e-02
    ## 36  1.199584e-02
    ## 37  5.612004e-02
    ## 38  5.612004e-02
    ## 39  3.058339e-02
    ## 40  1.159124e-02
    ## 41  1.901993e-03
    ## 42  1.353602e-02
    ## 43  1.868356e-03
    ## 44  1.289919e-02
    ## 45  6.001780e-03
    ## 46  4.815856e-02
    ## 47  5.676038e-03
    ## 48  1.648661e-02
    ## 49  4.342768e-01
    ## 50  4.342768e-01
    ## 51  3.082818e-01
    ## 52  3.082818e-01
    ## 53  3.135255e-01
    ## 54  3.135255e-01
    ## 55  4.045439e-01
    ## 56  4.045439e-01
    ## 57  2.447405e-11
    ## 58  2.902899e-07
    ## 59  3.513006e-10
    ## 60  3.507716e-06
    ## 61  6.665147e-09
    ## 62  1.105392e-05
    ## 63  3.261323e-12
    ## 64  8.249234e-08
    ## 65  6.687412e-04
    ## 66  1.912523e-03
    ## 67  6.382517e-04
    ## 68  2.459003e-03
    ## 69  1.016653e-03
    ## 70  7.638246e-03
    ## 71  4.890009e-03
    ## 72  9.304115e-03
    ## 73  7.417632e-03
    ## 74  1.387244e-02
    ## 75  6.628571e-03
    ## 76  1.479727e-02
    ## 77  1.036277e-02
    ## 78  2.609107e-02
    ## 79  1.352101e-01
    ## 80  2.871338e-02
    ## 81  2.593707e-01
    ## 82  1.461432e-01
    ## 83  4.559807e-01
    ## 84  4.078522e-01
    ## 85  3.555342e-01
    ## 86  1.205692e-01
    ## 87  7.494539e-01
    ## 88  3.300796e-01
    ## 89  7.925049e-01
    ## 90  3.075813e-01
    ## 91  7.467878e-01
    ## 92  3.356684e-01
    ## 93  9.861100e-01
    ## 94  3.917601e-01
    ## 95  7.372972e-01
    ## 96  3.640135e-01
    ## 97  9.908537e-01
    ## 98  9.335797e-01
    ## 99  7.955317e-01
    ## 100 7.955317e-01
    ## 101 8.865959e-01
    ## 102 8.865959e-01
