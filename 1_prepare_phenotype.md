``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(tableone)
library(survey)
```

    ## Loading required package: grid
    ## Loading required package: Matrix

    ## Warning: package 'Matrix' was built under R version 4.3.1

    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack
    ## 
    ## Loading required package: survival
    ## 
    ## Attaching package: 'survey'
    ## 
    ## The following object is masked from 'package:graphics':
    ## 
    ##     dotchart

``` r
library(readxl)
```

``` r
base_path <- "~/OneDrive-BethIsraelLaheyHealth"
project_folder <- "2022_Insomnia_metabolites_project"
# phenotype_file <- file.path(base_path, project_folder, "Phenotypes/metsleep_covariates_20200715.csv")
phenotype_file <- "metsleep_covariates_20200715.csv"

result_folder <- "Results"

# metab_data_file <- file.path(base_path, "metabolites info/SOL_with_Batch2/2batch_combined_data_V1_only.xlsx")
```

``` r
dat <- read.csv(phenotype_file)
cols_to_keep <- c(
  "ID", "STRAT", "PSU_ID", "WEIGHT_FINAL_NORM_OVERALL",
  "CENTER", "AGE", "BMI", "GENDER", "ALCOHOL_USE",
  "CIGARETTE_USE", "HYPERTENSION", "BKGRD1_C7", "DIABETES2_INDICATOR",
  "ESS", "WHIIRS", "slpa54gt5", "SLPDUR", "GPAQ_TOTAL_MET"
)

dat <- dat[, cols_to_keep]
```

``` r
dat <- dat %>%
  mutate(
    Gender = factor(GENDER, levels = c("F", "M"), labels = c("Female", "Male")),
    Center = factor(CENTER, levels = c("B", "C", "M", "S"), labels = c("Bronx", "Chicago", "Miami", "San Diego")),
    Background = factor(BKGRD1_C7,
      levels = c("0", "1", "2", "3", "4", "5", "6"),
      labels = c(
        "Domician", "Central American", "Cuban", "Mexican",
        "Puerto Rican", "South American", "More than one/Other heritage"
      )
    ),
    Survey_weight = WEIGHT_FINAL_NORM_OVERALL,
    Alcohol_use = factor(ALCOHOL_USE, levels = c("1", "2", "3"), labels = c("Never", "Former", "current")),
    Smoking = factor(CIGARETTE_USE, levels = c("1", "2", "3"), labels = c("Never", "Former", "current")),
    Diabetes = factor(DIABETES2_INDICATOR, c("0", "1"), labels = c("No", "Yes")),
    EDS = factor(ifelse(ESS > 10, "Yes", "No")),
    Insomnia = factor(ifelse(WHIIRS > 9, "Yes", "No")),
    Mild_OSA = factor(slpa54gt5, c("0", "1"), labels = c("No", "Yes")),
    Sleep_duration = SLPDUR,
    Total_phys = GPAQ_TOTAL_MET
  )
```

``` r
#saveRDS(dat, "Results/Processed_pheno.RDS")
#saveRDS(dat, file.path(base_path, project_folder, result_folder, "Processed_pheno.RDS"))
```
