``` r
library(tidyverse)
```

# Set path

``` r
base_path <- "~/OneDrive-BethIsraelLaheyHealth"
project_folder <- "2022_Insomnia_metabolites_project"
# phenotype_file <- file.path(base_path, project_folder, "Phenotypes/metsleep_covariates_20200715.csv")
phenotype_file <- "metsleep_covariates_20200715.csv"

result_folder <- "Results"

# metab_data_file <- file.path(base_path, "metabolites info/SOL_with_Batch2/2batch_combined_data_V1_only.xlsx")
```

# Read and process phenotype data

``` r
dat <- read.csv(phenotype_file)

dat$Insomnia <- ifelse(dat$WHIIRS >=9, 1, 0)
dat$CENTER <- factor(dat$CENTER)
dat$GENDER <- factor(dat$GENDER)
dat$Background <- factor(dat$BKGRD1_C7)

dat$SLEA4 <- as.numeric(as.character(dat$SLEA4))
dat$SLEA5 <- as.numeric(as.character(dat$SLEA5))
dat$SLEA6 <- as.numeric(as.character(dat$SLEA6))
dat$SLEA7 <- as.numeric(as.character(dat$SLEA7))

dat <- dat[!duplicated(dat$ID),]
dat$SOL_ID  <- as.character(dat$ID)
```

# Save phenotype data

``` r
#saveRDS(dat, "Processed_pheno.RDS")
#saveRDS(dat, file.path(base_path, project_folder, result_folder, "Processed_pheno.RDS"))
```
