``` r
library(tidyverse)
library(tableone)
library(survey)
library(readxl)
```

# Batch strata

# Read data

``` r
dat <- read.csv("metsleep_covariates_20200715.csv")

cols_to_keep <- c(
  "ID", "STRAT", "PSU_ID", "WEIGHT_FINAL_NORM_OVERALL",
  "CENTER", "AGE", "BMI", "GENDER", "ALCOHOL_USE",
  "CIGARETTE_USE", "HYPERTENSION", "BKGRD1_C7", "DIABETES2_INDICATOR",
  "ESS", "WHIIRS", "slpa54gt5", "SLPDUR", "GPAQ_TOTAL_MET",
  "OCEA13"
)

dat <- dat[, cols_to_keep]

metab_file_b1 <- "Processed_metab_b1_updated.RDS"
metab_file_b2 <- "Processed_metab_b2_updated.RDS"
metab_b1 <- readRDS(metab_file_b1)
metab_b2 <- readRDS(metab_file_b2)

col_chems <- setdiff(1:ncol(metab_b1), grep("SOL_ID", colnames(metab_b1)))
colnames(metab_b1)[col_chems] <- colnames(metab_b2)[col_chems] <- paste0("chem_", colnames(metab_b1)[col_chems])

dat$SOL_ID <- as.character(dat$ID)
dat <- dat[!duplicated(dat$SOL_ID), ]

dat_b1 <- merge(dat, metab_b1, by = "SOL_ID", all.y = TRUE)
dat_b1 <- dat_b1[!is.na(dat_b1$WHIIRS), ]
dat_b2 <- merge(dat, metab_b2, by = "SOL_ID", all.y = TRUE)
dat_b2 <- dat_b2[!is.na(dat_b2$WHIIRS), ]
dat_b1$Batch <- "Batch 1"
dat_b2$Batch <- "Batch 2"

dat_comb <- rbind(dat_b1, dat_b2)
```

# Recoding variables

``` r
dat_comb <- dat_comb %>%
  mutate(
    Gender = factor(GENDER, levels = c("F", "M"), labels = c("Female", "Male")),
    Center = factor(CENTER, levels = c("B", "C", "M", "S"), labels = c("Bronx", "Chicago", "Miami", "San Diego")),
    Background = factor(BKGRD1_C7, levels = c("0", "1", "2", "3", "4", "5", "6"), labels = c("Domician", "Central American", "Cuban", "Mexican", "Puerto Rican", "South American", "More than one/Other heritage")),
    Survey_weight = WEIGHT_FINAL_NORM_OVERALL,
    Alcohol_use = factor(ALCOHOL_USE, levels = c("1", "2", "3"), labels = c("Never", "Former", "current")),
    Smoking = factor(CIGARETTE_USE, levels = c("1", "2", "3"), labels = c("Never", "Former", "current")),
    Diabetes = factor(DIABETES2_INDICATOR, c("0", "1"), labels = c("No", "Yes")),
    EDS = factor(ifelse(ESS > 10, "Yes", "No")),
    Insomnia = factor(ifelse(WHIIRS >= 9, "Yes", "No")),
    Mild_to_severe_OSA = factor(slpa54gt5, c("0", "1"), labels = c("No", "Yes")),
    Sleep_duration = SLPDUR,
    Total_phys = GPAQ_TOTAL_MET,
    Shift_work = factor(ifelse(OCEA13 %in% c("2", "3", "4", "5", "6"), "Yes", "No"))
  )
```

# Batch stratified (Table 1)

``` r
survey_obj <- svydesign(id = ~PSU_ID, strata = ~STRAT, weights = ~Survey_weight, data = dat_comb)

tbl1_vars <- c(
  "Gender", "Batch", "Center", "BMI", "AGE",
  "Background", "Alcohol_use",
  "Smoking", "Diabetes", "EDS", "Insomnia",
  "Mild_to_severe_OSA", "Sleep_duration", "Total_phys",
  "Shift_work"
)


tbl1_noweight <- print(CreateTableOne(vars = tbl1_vars, data = dat_comb, strata = "Batch"), missing = TRUE, varLabels = TRUE, digits = 3, pDigits = 3, showAllLevels = TRUE)
```

    ##                             Stratified by Batch
    ##                              level                        Batch 1        
    ##   n                                                         3895         
    ##   Gender (%)                 Female                         2235 ( 57.4) 
    ##                              Male                           1660 ( 42.6) 
    ##   Batch (%)                  Batch 1                        3895 (100.0) 
    ##                              Batch 2                           0 (  0.0) 
    ##   Center (%)                 Bronx                          1018 ( 26.1) 
    ##                              Chicago                         903 ( 23.2) 
    ##                              Miami                          1088 ( 27.9) 
    ##                              San Diego                       886 ( 22.7) 
    ##   BMI (mean (SD))                                          29.80 (6.08)  
    ##   AGE (mean (SD))                                          45.90 (13.80) 
    ##   Background (%)             Domician                        377 (  9.7) 
    ##                              Central American                398 ( 10.2) 
    ##                              Cuban                           653 ( 16.8) 
    ##                              Mexican                        1417 ( 36.5) 
    ##                              Puerto Rican                    693 ( 17.8) 
    ##                              South American                  227 (  5.8) 
    ##                              More than one/Other heritage    122 (  3.1) 
    ##   Alcohol_use (%)            Never                           718 ( 18.4) 
    ##                              Former                         1253 ( 32.2) 
    ##                              current                        1923 ( 49.4) 
    ##   Smoking (%)                Never                          2281 ( 58.6) 
    ##                              Former                          778 ( 20.0) 
    ##                              current                         833 ( 21.4) 
    ##   Diabetes (%)               No                             3160 ( 81.1) 
    ##                              Yes                             735 ( 18.9) 
    ##   EDS (%)                    No                             3278 ( 84.4) 
    ##                              Yes                             606 ( 15.6) 
    ##   Insomnia (%)               No                             2489 ( 63.9) 
    ##                              Yes                            1406 ( 36.1) 
    ##   Mild_to_severe_OSA (%)     No                             2409 ( 69.4) 
    ##                              Yes                            1062 ( 30.6) 
    ##   Sleep_duration (mean (SD))                                7.93 (1.43)  
    ##   Total_phys (mean (SD))                                  627.37 (974.05)
    ##   Shift_work (%)             No                             3154 ( 81.0) 
    ##                              Yes                             741 ( 19.0) 
    ##                             Stratified by Batch
    ##                              Batch 2         p      test Missing
    ##   n                            2121                             
    ##   Gender (%)                   1373 ( 64.7)  <0.001       0.0   
    ##                                 748 ( 35.3)                     
    ##   Batch (%)                       0 (  0.0)  <0.001       0.0   
    ##                                2121 (100.0)                     
    ##   Center (%)                    547 ( 25.8)   0.143       0.0   
    ##                                 538 ( 25.4)                     
    ##                                 597 ( 28.1)                     
    ##                                 439 ( 20.7)                     
    ##   BMI (mean (SD))             30.20 (6.11)    0.015       0.3   
    ##   AGE (mean (SD))             52.73 (10.57)  <0.001       0.0   
    ##   Background (%)                269 ( 12.7)  <0.001       0.2   
    ##                                 227 ( 10.7)                     
    ##                                 393 ( 18.5)                     
    ##                                 653 ( 30.8)                     
    ##                                 380 ( 17.9)                     
    ##                                 170 (  8.0)                     
    ##                                  27 (  1.3)                     
    ##   Alcohol_use (%)               526 ( 24.8)  <0.001       0.0   
    ##                                 696 ( 32.8)                     
    ##                                 897 ( 42.3)                     
    ##   Smoking (%)                  1244 ( 58.7)   0.003       0.1   
    ##                                 486 ( 22.9)                     
    ##                                 389 ( 18.4)                     
    ##   Diabetes (%)                 1522 ( 71.8)  <0.001       0.0   
    ##                                 599 ( 28.2)                     
    ##   EDS (%)                      1773 ( 83.9)   0.676       0.3   
    ##                                 339 ( 16.1)                     
    ##   Insomnia (%)                 1242 ( 58.6)  <0.001       0.0   
    ##                                 879 ( 41.4)                     
    ##   Mild_to_severe_OSA (%)       1192 ( 62.3)  <0.001      10.5   
    ##                                 722 ( 37.7)                     
    ##   Sleep_duration (mean (SD))   7.88 (1.47)    0.250       3.5   
    ##   Total_phys (mean (SD))     487.87 (846.86) <0.001       0.2   
    ##   Shift_work (%)               1804 ( 85.1)  <0.001       0.0   
    ##                                 317 ( 14.9)

``` r
tbl1_weighted <- print(svyCreateTableOne(vars = tbl1_vars, data = survey_obj, strata = "Batch"), missing = TRUE, varLabels = TRUE, digits = 3, pDigits = 3, showAllLevels = TRUE)
```

    ##                             Stratified by Batch
    ##                              level                        Batch 1          
    ##   n                                                        3988.9          
    ##   Gender (%)                 Female                        2037.4 ( 51.1)  
    ##                              Male                          1951.5 ( 48.9)  
    ##   Batch (%)                  Batch 1                       3988.9 (100.0)  
    ##                              Batch 2                          0.0 (  0.0)  
    ##   Center (%)                 Bronx                         1148.5 ( 28.8)  
    ##                              Chicago                        586.9 ( 14.7)  
    ##                              Miami                         1265.4 ( 31.7)  
    ##                              San Diego                      988.2 ( 24.8)  
    ##   BMI (mean (SD))                                           29.42 (6.31)   
    ##   AGE (mean (SD))                                           41.57 (15.14)  
    ##   Background (%)             Domician                       410.0 ( 10.3)  
    ##                              Central American               265.8 (  6.7)  
    ##                              Cuban                          892.1 ( 22.4)  
    ##                              Mexican                       1383.1 ( 34.8)  
    ##                              Puerto Rican                   670.5 ( 16.9)  
    ##                              South American                 178.2 (  4.5)  
    ##                              More than one/Other heritage   177.9 (  4.5)  
    ##   Alcohol_use (%)            Never                          665.9 ( 16.7)  
    ##                              Former                        1153.0 ( 28.9)  
    ##                              current                       2169.7 ( 54.4)  
    ##   Smoking (%)                Never                         2375.5 ( 59.6)  
    ##                              Former                         632.7 ( 15.9)  
    ##                              current                        978.8 ( 24.5)  
    ##   Diabetes (%)               No                            3383.5 ( 84.8)  
    ##                              Yes                            605.3 ( 15.2)  
    ##   EDS (%)                    No                            3390.1 ( 85.2)  
    ##                              Yes                            590.0 ( 14.8)  
    ##   Insomnia (%)               No                            2656.6 ( 66.6)  
    ##                              Yes                           1332.3 ( 33.4)  
    ##   Mild_to_severe_OSA (%)     No                            2566.4 ( 73.2)  
    ##                              Yes                            940.2 ( 26.8)  
    ##   Sleep_duration (mean (SD))                                 7.98 (1.45)   
    ##   Total_phys (mean (SD))                                   698.12 (1025.17)
    ##   Shift_work (%)             No                            3206.2 ( 80.4)  
    ##                              Yes                            782.7 ( 19.6)  
    ##                             Stratified by Batch
    ##                              Batch 2          p      test Missing
    ##   n                           1625.6                             
    ##   Gender (%)                   923.6 ( 56.8)   0.004       0.0   
    ##                                702.0 ( 43.2)                     
    ##   Batch (%)                      0.0 (  0.0)  <0.001       0.0   
    ##                               1625.6 (100.0)                     
    ##   Center (%)                   450.0 ( 27.7)   0.009       0.0   
    ##                                236.1 ( 14.5)                     
    ##                                623.1 ( 38.3)                     
    ##                                316.5 ( 19.5)                     
    ##   BMI (mean (SD))              29.73 (5.82)    0.217       0.3   
    ##   AGE (mean (SD))              51.17 (13.19)  <0.001       0.0   
    ##   Background (%)               209.3 ( 12.9)  <0.001       0.2   
    ##                                109.8 (  6.8)                     
    ##                                494.9 ( 30.5)                     
    ##                                436.2 ( 26.9)                     
    ##                                263.6 ( 16.2)                     
    ##                                 90.6 (  5.6)                     
    ##                                 19.1 (  1.2)                     
    ##   Alcohol_use (%)              425.0 ( 26.2)  <0.001       0.0   
    ##                                513.3 ( 31.6)                     
    ##                                686.4 ( 42.2)                     
    ##   Smoking (%)                  964.0 ( 59.3)  <0.001       0.1   
    ##                                367.5 ( 22.6)                     
    ##                                293.5 ( 18.1)                     
    ##   Diabetes (%)                1197.7 ( 73.7)  <0.001       0.0   
    ##                                428.0 ( 26.3)                     
    ##   EDS (%)                     1369.2 ( 84.5)   0.650       0.3   
    ##                                250.2 ( 15.5)                     
    ##   Insomnia (%)                1014.6 ( 62.4)   0.021       0.0   
    ##                                611.0 ( 37.6)                     
    ##   Mild_to_severe_OSA (%)       921.0 ( 64.2)  <0.001      10.5   
    ##                                513.8 ( 35.8)                     
    ##   Sleep_duration (mean (SD))    7.99 (1.49)    0.817       3.5   
    ##   Total_phys (mean (SD))      495.10 (880.54) <0.001       0.2   
    ##   Shift_work (%)              1370.5 ( 84.3)   0.006       0.0   
    ##                                255.1 ( 15.7)

Combine the count values (from the unweighted table) and percentage
(from the weighted table) together

``` r
tbl1_w <- tbl1_weighted[, -which(colnames(tbl1_noweight) %in% c("p", "test"))]
tbl1 <- tbl1_noweight[, -which(colnames(tbl1_noweight) %in% c("p", "test"))]

tbl1_comb <- tbl1
col_inds_to_update <- which(colnames(tbl1_w) %in% c("Batch 1", "Batch 2"))

# update tbl1_comb with the percentages from the weighted table
for (i in col_inds_to_update) {
  counts <- sapply(tbl1[, i], function(x) {
    strsplit(x, split = "(", fixed = TRUE)[[1]][1]
  })
  percent <- sapply(tbl1_w[, i], function(x) {
    paste0("(", strsplit(x, split = "(", fixed = TRUE)[[1]][2])
  })
  tbl1_comb[, i] <- paste0(counts, percent)
}


#write.csv(tbl1_comb, "Results/20240408_table1_Batch_strat_V3.csv")
tbl1_comb
```

    ##                             Stratified by Batch
    ##                              level                          Batch 1           
    ##   n                          ""                             "  3895(NA"       
    ##   Gender (%)                 "Female"                       "  2235 ( 51.1) " 
    ##                              "Male"                         "  1660 ( 48.9) " 
    ##   Batch (%)                  "Batch 1"                      "  3895 (100.0) " 
    ##                              "Batch 2"                      "     0 (  0.0) " 
    ##   Center (%)                 "Bronx"                        "  1018 ( 28.8) " 
    ##                              "Chicago"                      "   903 ( 14.7) " 
    ##                              "Miami"                        "  1088 ( 31.7) " 
    ##                              "San Diego"                    "   886 ( 24.8) " 
    ##   BMI (mean (SD))            ""                             " 29.80 (6.31)"   
    ##   AGE (mean (SD))            ""                             " 45.90 (15.14)"  
    ##   Background (%)             "Domician"                     "   377 ( 10.3) " 
    ##                              "Central American"             "   398 (  6.7) " 
    ##                              "Cuban"                        "   653 ( 22.4) " 
    ##                              "Mexican"                      "  1417 ( 34.8) " 
    ##                              "Puerto Rican"                 "   693 ( 16.9) " 
    ##                              "South American"               "   227 (  4.5) " 
    ##                              "More than one/Other heritage" "   122 (  4.5) " 
    ##   Alcohol_use (%)            "Never"                        "   718 ( 16.7) " 
    ##                              "Former"                       "  1253 ( 28.9) " 
    ##                              "current"                      "  1923 ( 54.4) " 
    ##   Smoking (%)                "Never"                        "  2281 ( 59.6) " 
    ##                              "Former"                       "   778 ( 15.9) " 
    ##                              "current"                      "   833 ( 24.5) " 
    ##   Diabetes (%)               "No"                           "  3160 ( 84.8) " 
    ##                              "Yes"                          "   735 ( 15.2) " 
    ##   EDS (%)                    "No"                           "  3278 ( 85.2) " 
    ##                              "Yes"                          "   606 ( 14.8) " 
    ##   Insomnia (%)               "No"                           "  2489 ( 66.6) " 
    ##                              "Yes"                          "  1406 ( 33.4) " 
    ##   Mild_to_severe_OSA (%)     "No"                           "  2409 ( 73.2) " 
    ##                              "Yes"                          "  1062 ( 26.8) " 
    ##   Sleep_duration (mean (SD)) ""                             "  7.93 (1.45)"   
    ##   Total_phys (mean (SD))     ""                             "627.37 (1025.17)"
    ##   Shift_work (%)             "No"                           "  3154 ( 80.4) " 
    ##                              "Yes"                          "   741 ( 19.6) " 
    ##                             Stratified by Batch
    ##                              Batch 2           Missing
    ##   n                          "  2121(NA"       "    " 
    ##   Gender (%)                 "  1373 ( 56.8) " " 0.0" 
    ##                              "   748 ( 43.2) " "    " 
    ##   Batch (%)                  "     0 (  0.0) " " 0.0" 
    ##                              "  2121 (100.0) " "    " 
    ##   Center (%)                 "   547 ( 27.7) " " 0.0" 
    ##                              "   538 ( 14.5) " "    " 
    ##                              "   597 ( 38.3) " "    " 
    ##                              "   439 ( 19.5) " "    " 
    ##   BMI (mean (SD))            " 30.20 (5.82)"   " 0.3" 
    ##   AGE (mean (SD))            " 52.73 (13.19)"  " 0.0" 
    ##   Background (%)             "   269 ( 12.9) " " 0.2" 
    ##                              "   227 (  6.8) " "    " 
    ##                              "   393 ( 30.5) " "    " 
    ##                              "   653 ( 26.9) " "    " 
    ##                              "   380 ( 16.2) " "    " 
    ##                              "   170 (  5.6) " "    " 
    ##                              "    27 (  1.2) " "    " 
    ##   Alcohol_use (%)            "   526 ( 26.2) " " 0.0" 
    ##                              "   696 ( 31.6) " "    " 
    ##                              "   897 ( 42.2) " "    " 
    ##   Smoking (%)                "  1244 ( 59.3) " " 0.1" 
    ##                              "   486 ( 22.6) " "    " 
    ##                              "   389 ( 18.1) " "    " 
    ##   Diabetes (%)               "  1522 ( 73.7) " " 0.0" 
    ##                              "   599 ( 26.3) " "    " 
    ##   EDS (%)                    "  1773 ( 84.5) " " 0.3" 
    ##                              "   339 ( 15.5) " "    " 
    ##   Insomnia (%)               "  1242 ( 62.4) " " 0.0" 
    ##                              "   879 ( 37.6) " "    " 
    ##   Mild_to_severe_OSA (%)     "  1192 ( 64.2) " "10.5" 
    ##                              "   722 ( 35.8) " "    " 
    ##   Sleep_duration (mean (SD)) "  7.88 (1.49)"   " 3.5" 
    ##   Total_phys (mean (SD))     "487.87 (880.54)" " 0.2" 
    ##   Shift_work (%)             "  1804 ( 84.3) " " 0.0" 
    ##                              "   317 ( 15.7) " "    "

# Sex stratified (Supplementary Table 1)

``` r
survey_obj <- svydesign(id = ~PSU_ID, strata = ~STRAT, weights = ~Survey_weight, data = dat_comb)

tbl1_vars <- c(
  "Gender", "Batch", "Center", "BMI", "AGE",
  "Background", "Alcohol_use",
  "Smoking", "Diabetes", "EDS", "Insomnia",
  "Mild_to_severe_OSA", "Sleep_duration", "Total_phys",
  "Shift_work"
)


tbl1_noweight <- print(CreateTableOne(vars = tbl1_vars, data = dat_comb, strata = "Gender"), missing = TRUE, varLabels = TRUE, digits = 3, pDigits = 3, showAllLevels = TRUE)
```

    ##                             Stratified by Gender
    ##                              level                        Female         
    ##   n                                                         3608         
    ##   Gender (%)                 Female                         3608 (100.0) 
    ##                              Male                              0 (  0.0) 
    ##   Batch (%)                  Batch 1                        2235 ( 61.9) 
    ##                              Batch 2                        1373 ( 38.1) 
    ##   Center (%)                 Bronx                           978 ( 27.1) 
    ##                              Chicago                         809 ( 22.4) 
    ##                              Miami                           987 ( 27.4) 
    ##                              San Diego                       834 ( 23.1) 
    ##   BMI (mean (SD))                                          30.53 (6.47)  
    ##   AGE (mean (SD))                                          48.83 (12.76) 
    ##   Background (%)             Domician                        434 ( 12.0) 
    ##                              Central American                400 ( 11.1) 
    ##                              Cuban                           561 ( 15.6) 
    ##                              Mexican                        1275 ( 35.4) 
    ##                              Puerto Rican                    607 ( 16.9) 
    ##                              South American                  241 (  6.7) 
    ##                              More than one/Other heritage     84 (  2.3) 
    ##   Alcohol_use (%)            Never                          1023 ( 28.4) 
    ##                              Former                         1219 ( 33.8) 
    ##                              current                        1363 ( 37.8) 
    ##   Smoking (%)                Never                          2422 ( 67.1) 
    ##                              Former                          607 ( 16.8) 
    ##                              current                         578 ( 16.0) 
    ##   Diabetes (%)               No                             2813 ( 78.0) 
    ##                              Yes                             795 ( 22.0) 
    ##   EDS (%)                    No                             3044 ( 84.7) 
    ##                              Yes                             550 ( 15.3) 
    ##   Insomnia (%)               No                             2066 ( 57.3) 
    ##                              Yes                            1542 ( 42.7) 
    ##   Mild_to_severe_OSA (%)     No                             2357 ( 72.8) 
    ##                              Yes                             882 ( 27.2) 
    ##   Sleep_duration (mean (SD))                                7.96 (1.44)  
    ##   Total_phys (mean (SD))                                  408.11 (703.01)
    ##   Shift_work (%)             No                             3070 ( 85.1) 
    ##                              Yes                             538 ( 14.9) 
    ##                             Stratified by Gender
    ##                              Male             p      test Missing
    ##   n                            2408                              
    ##   Gender (%)                      0 (  0.0)   <0.001       0.0   
    ##                                2408 (100.0)                      
    ##   Batch (%)                    1660 ( 68.9)   <0.001       0.0   
    ##                                 748 ( 31.1)                      
    ##   Center (%)                    587 ( 24.4)   <0.001       0.0   
    ##                                 632 ( 26.2)                      
    ##                                 698 ( 29.0)                      
    ##                                 491 ( 20.4)                      
    ##   BMI (mean (SD))             29.05 (5.34)    <0.001       0.3   
    ##   AGE (mean (SD))             47.52 (13.72)   <0.001       0.0   
    ##   Background (%)                212 (  8.8)   <0.001       0.2   
    ##                                 225 (  9.4)                      
    ##                                 485 ( 20.2)                      
    ##                                 795 ( 33.1)                      
    ##                                 466 ( 19.4)                      
    ##                                 156 (  6.5)                      
    ##                                  65 (  2.7)                      
    ##   Alcohol_use (%)               221 (  9.2)   <0.001       0.0   
    ##                                 730 ( 30.3)                      
    ##                                1457 ( 60.5)                      
    ##   Smoking (%)                  1103 ( 45.9)   <0.001       0.1   
    ##                                 657 ( 27.3)                      
    ##                                 644 ( 26.8)                      
    ##   Diabetes (%)                 1869 ( 77.6)    0.773       0.0   
    ##                                 539 ( 22.4)                      
    ##   EDS (%)                      2007 ( 83.6)    0.249       0.3   
    ##                                 395 ( 16.4)                      
    ##   Insomnia (%)                 1665 ( 69.1)   <0.001       0.0   
    ##                                 743 ( 30.9)                      
    ##   Mild_to_severe_OSA (%)       1244 ( 58.0)   <0.001      10.5   
    ##                                 902 ( 42.0)                      
    ##   Sleep_duration (mean (SD))   7.84 (1.44)     0.003       3.5   
    ##   Total_phys (mean (SD))     833.13 (1152.82) <0.001       0.2   
    ##   Shift_work (%)               1888 ( 78.4)   <0.001       0.0   
    ##                                 520 ( 21.6)

``` r
tbl1_weighted <- print(svyCreateTableOne(vars = tbl1_vars, data = survey_obj, strata = "Gender"), missing = TRUE, varLabels = TRUE, digits = 3, pDigits = 3, showAllLevels = TRUE)
```

    ##                             Stratified by Gender
    ##                              level                        Female          
    ##   n                                                        2961.0         
    ##   Gender (%)                 Female                        2961.0 (100.0) 
    ##                              Male                             0.0 (  0.0) 
    ##   Batch (%)                  Batch 1                       2037.4 ( 68.8) 
    ##                              Batch 2                        923.6 ( 31.2) 
    ##   Center (%)                 Bronx                          914.1 ( 30.9) 
    ##                              Chicago                        390.0 ( 13.2) 
    ##                              Miami                          976.0 ( 33.0) 
    ##                              San Diego                      680.9 ( 23.0) 
    ##   BMI (mean (SD))                                           30.17 (6.75)  
    ##   AGE (mean (SD))                                           44.92 (15.14) 
    ##   Background (%)             Domician                       391.2 ( 13.3) 
    ##                              Central American               214.2 (  7.3) 
    ##                              Cuban                          666.2 ( 22.6) 
    ##                              Mexican                        970.0 ( 32.9) 
    ##                              Puerto Rican                   451.2 ( 15.3) 
    ##                              South American                 159.7 (  5.4) 
    ##                              More than one/Other heritage    97.3 (  3.3) 
    ##   Alcohol_use (%)            Never                          818.3 ( 27.6) 
    ##                              Former                         926.2 ( 31.3) 
    ##                              current                       1215.4 ( 41.1) 
    ##   Smoking (%)                Never                         2034.5 ( 68.7) 
    ##                              Former                         394.2 ( 13.3) 
    ##                              current                        532.0 ( 18.0) 
    ##   Diabetes (%)               No                            2420.5 ( 81.7) 
    ##                              Yes                            540.5 ( 18.3) 
    ##   EDS (%)                    No                            2536.6 ( 85.9) 
    ##                              Yes                            416.9 ( 14.1) 
    ##   Insomnia (%)               No                            1776.7 ( 60.0) 
    ##                              Yes                           1184.3 ( 40.0) 
    ##   Mild_to_severe_OSA (%)     No                            2037.2 ( 77.8) 
    ##                              Yes                            580.2 ( 22.2) 
    ##   Sleep_duration (mean (SD))                                 8.05 (1.49)  
    ##   Total_phys (mean (SD))                                   426.10 (717.67)
    ##   Shift_work (%)             No                            2545.6 ( 86.0) 
    ##                              Yes                            415.4 ( 14.0) 
    ##                             Stratified by Gender
    ##                              Male              p      test Missing
    ##   n                           2653.5                              
    ##   Gender (%)                     0.0 (  0.0)   <0.001       0.0   
    ##                               2653.5 (100.0)                      
    ##   Batch (%)                   1951.5 ( 73.5)    0.004       0.0   
    ##                                702.0 ( 26.5)                      
    ##   Center (%)                   684.4 ( 25.8)    0.013       0.0   
    ##                                432.9 ( 16.3)                      
    ##                                912.5 ( 34.4)                      
    ##                                623.7 ( 23.5)                      
    ##   BMI (mean (SD))              28.78 (5.37)    <0.001       0.3   
    ##   AGE (mean (SD))              43.71 (15.33)    0.046       0.0   
    ##   Background (%)               228.1 (  8.6)   <0.001       0.2   
    ##                                161.3 (  6.1)                      
    ##                                720.8 ( 27.2)                      
    ##                                849.4 ( 32.0)                      
    ##                                482.9 ( 18.2)                      
    ##                                109.1 (  4.1)                      
    ##                                 99.7 (  3.8)                      
    ##   Alcohol_use (%)              272.7 ( 10.3)   <0.001       0.0   
    ##                                740.1 ( 27.9)                      
    ##                               1640.7 ( 61.8)                      
    ##   Smoking (%)                 1305.0 ( 49.2)   <0.001       0.1   
    ##                                605.9 ( 22.9)                      
    ##                                740.2 ( 27.9)                      
    ##   Diabetes (%)                2160.7 ( 81.4)    0.805       0.0   
    ##                                492.8 ( 18.6)                      
    ##   EDS (%)                     2222.7 ( 84.0)    0.158       0.3   
    ##                                423.3 ( 16.0)                      
    ##   Insomnia (%)                1894.5 ( 71.4)   <0.001       0.0   
    ##                                759.0 ( 28.6)                      
    ##   Mild_to_severe_OSA (%)      1450.2 ( 62.4)   <0.001      10.5   
    ##                                873.8 ( 37.6)                      
    ##   Sleep_duration (mean (SD))    7.91 (1.42)     0.023       3.5   
    ##   Total_phys (mean (SD))      877.31 (1179.39) <0.001       0.2   
    ##   Shift_work (%)              2031.1 ( 76.5)   <0.001       0.0   
    ##                                622.4 ( 23.5)

Combine the count values (from the unweighted table) and percentage
(from the weighted table) together

``` r
tbl1_w <- tbl1_weighted[, -which(colnames(tbl1_noweight) %in% c("p", "test"))]
tbl1 <- tbl1_noweight[, -which(colnames(tbl1_noweight) %in% c("p", "test"))]

tbl1_comb <- tbl1
col_inds_to_update <- which(colnames(tbl1_w) %in% c("Female", "Male"))

# update tbl1_comb with the percentages from the weighted table
for (i in col_inds_to_update) {
  counts <- sapply(tbl1[, i], function(x) {
    strsplit(x, split = "(", fixed = TRUE)[[1]][1]
  })
  percent <- sapply(tbl1_w[, i], function(x) {
    paste0("(", strsplit(x, split = "(", fixed = TRUE)[[1]][2])
  })
  tbl1_comb[, i] <- paste0(counts, percent)
}

tbl1_comb
```

    ##                             Stratified by Gender
    ##                              level                          Female           
    ##   n                          ""                             "  3608(NA"      
    ##   Gender (%)                 "Female"                       "  3608 (100.0) "
    ##                              "Male"                         "     0 (  0.0) "
    ##   Batch (%)                  "Batch 1"                      "  2235 ( 68.8) "
    ##                              "Batch 2"                      "  1373 ( 31.2) "
    ##   Center (%)                 "Bronx"                        "   978 ( 30.9) "
    ##                              "Chicago"                      "   809 ( 13.2) "
    ##                              "Miami"                        "   987 ( 33.0) "
    ##                              "San Diego"                    "   834 ( 23.0) "
    ##   BMI (mean (SD))            ""                             " 30.53 (6.75)"  
    ##   AGE (mean (SD))            ""                             " 48.83 (15.14)" 
    ##   Background (%)             "Domician"                     "   434 ( 13.3) "
    ##                              "Central American"             "   400 (  7.3) "
    ##                              "Cuban"                        "   561 ( 22.6) "
    ##                              "Mexican"                      "  1275 ( 32.9) "
    ##                              "Puerto Rican"                 "   607 ( 15.3) "
    ##                              "South American"               "   241 (  5.4) "
    ##                              "More than one/Other heritage" "    84 (  3.3) "
    ##   Alcohol_use (%)            "Never"                        "  1023 ( 27.6) "
    ##                              "Former"                       "  1219 ( 31.3) "
    ##                              "current"                      "  1363 ( 41.1) "
    ##   Smoking (%)                "Never"                        "  2422 ( 68.7) "
    ##                              "Former"                       "   607 ( 13.3) "
    ##                              "current"                      "   578 ( 18.0) "
    ##   Diabetes (%)               "No"                           "  2813 ( 81.7) "
    ##                              "Yes"                          "   795 ( 18.3) "
    ##   EDS (%)                    "No"                           "  3044 ( 85.9) "
    ##                              "Yes"                          "   550 ( 14.1) "
    ##   Insomnia (%)               "No"                           "  2066 ( 60.0) "
    ##                              "Yes"                          "  1542 ( 40.0) "
    ##   Mild_to_severe_OSA (%)     "No"                           "  2357 ( 77.8) "
    ##                              "Yes"                          "   882 ( 22.2) "
    ##   Sleep_duration (mean (SD)) ""                             "  7.96 (1.49)"  
    ##   Total_phys (mean (SD))     ""                             "408.11 (717.67)"
    ##   Shift_work (%)             "No"                           "  3070 ( 86.0) "
    ##                              "Yes"                          "   538 ( 14.0) "
    ##                             Stratified by Gender
    ##                              Male               Missing
    ##   n                          "  2408(NA"        "    " 
    ##   Gender (%)                 "     0 (  0.0) "  " 0.0" 
    ##                              "  2408 (100.0) "  "    " 
    ##   Batch (%)                  "  1660 ( 73.5) "  " 0.0" 
    ##                              "   748 ( 26.5) "  "    " 
    ##   Center (%)                 "   587 ( 25.8) "  " 0.0" 
    ##                              "   632 ( 16.3) "  "    " 
    ##                              "   698 ( 34.4) "  "    " 
    ##                              "   491 ( 23.5) "  "    " 
    ##   BMI (mean (SD))            " 29.05 (5.37)"    " 0.3" 
    ##   AGE (mean (SD))            " 47.52 (15.33)"   " 0.0" 
    ##   Background (%)             "   212 (  8.6) "  " 0.2" 
    ##                              "   225 (  6.1) "  "    " 
    ##                              "   485 ( 27.2) "  "    " 
    ##                              "   795 ( 32.0) "  "    " 
    ##                              "   466 ( 18.2) "  "    " 
    ##                              "   156 (  4.1) "  "    " 
    ##                              "    65 (  3.8) "  "    " 
    ##   Alcohol_use (%)            "   221 ( 10.3) "  " 0.0" 
    ##                              "   730 ( 27.9) "  "    " 
    ##                              "  1457 ( 61.8) "  "    " 
    ##   Smoking (%)                "  1103 ( 49.2) "  " 0.1" 
    ##                              "   657 ( 22.9) "  "    " 
    ##                              "   644 ( 27.9) "  "    " 
    ##   Diabetes (%)               "  1869 ( 81.4) "  " 0.0" 
    ##                              "   539 ( 18.6) "  "    " 
    ##   EDS (%)                    "  2007 ( 84.0) "  " 0.3" 
    ##                              "   395 ( 16.0) "  "    " 
    ##   Insomnia (%)               "  1665 ( 71.4) "  " 0.0" 
    ##                              "   743 ( 28.6) "  "    " 
    ##   Mild_to_severe_OSA (%)     "  1244 ( 62.4) "  "10.5" 
    ##                              "   902 ( 37.6) "  "    " 
    ##   Sleep_duration (mean (SD)) "  7.84 (1.42)"    " 3.5" 
    ##   Total_phys (mean (SD))     "833.13 (1179.39)" " 0.2" 
    ##   Shift_work (%)             "  1888 ( 76.5) "  " 0.0" 
    ##                              "   520 ( 23.5) "  "    "
