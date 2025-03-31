``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(readxl)
library(glmnet)
```

    ## Loading required package: Matrix
    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack
    ## 
    ## Loaded glmnet 4.1-8

``` r
library(survey)
```

    ## Loading required package: grid
    ## Loading required package: survival
    ## 
    ## Attaching package: 'survey'
    ## 
    ## The following object is masked from 'package:graphics':
    ## 
    ##     dotchart

``` r
library(patchwork)
library(ggpubr)
```

# Read Data

``` r
m1_dat <- readRDS("Results/20240616_MRS_forPlot_M1_Ranked_V2.RDS")
m1_nomed_dat <- readRDS("Results/20240616_MRS_forPlot_M1_nomed_Ranked_V2.RDS")
```

# Process Data

``` r
m1_dat$Model <- "M1"
m1_nomed_dat$Model <- "M1_No_Med"
sum_dat <- rbind(m1_dat, m1_nomed_dat)
sum_dat$Exposures <- as.factor(sum_dat$Exposures)
new_order <- c("Insomnia", "WHIIRS", "SLEA4", "SLEA5", "SLEA6", "SLEA7", "SLEA11",
               "Insomnia1", "WHIIRS1", "SLEA41", "SLEA51", "SLEA61", "SLEA71", "SLEA111")
sum_dat$Exposures <- factor(sum_dat$Exposures, levels = new_order)
sum_dat$MODEL <- ifelse(sum_dat$Exposures == "Insomnia" & sum_dat$Model == "M1", "M1_OR", 
                                 ifelse(sum_dat$Exposures == "Insomnia" & sum_dat$Model == "M1_No_Med", "M1_No_Med_OR", 
                                        ifelse(sum_dat$Exposures != "Insomnia" & sum_dat$Model == "M1", "M1_Est", "M1_No_Med_Est")))
model_order <- c("M1_Est", "M1_No_Med_Est", "M1_OR", "M1_No_Med_OR")
sum_dat$MODEL <- factor(sum_dat$MODEL, levels = model_order)
```

# seperate table

``` r
tab_1 <- sum_dat[sum_dat$Exposures == "Insomnia",]
tab_2 <- sum_dat[sum_dat$Exposures != "Insomnia",]
```

# Plot 1

``` r
my_colors <- c("darkblue", "lightblue")

p1 <- ggplot(tab_1, aes(x = Exposures, y = Est, color = Model)) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), position = position_dodge(0.7), width = 0.2, linewidth = 1.5) +
  geom_point(position = position_dodge(0.7), size = 3) +
  geom_text(aes(label = Metabolites), hjust = -0.2, vjust = -0.3, size = 7, position = position_dodge(0.7)) +
  scale_color_manual(values = my_colors) +
  theme_bw(base_size = 20) +
  ylab("Odds Ratio") +
  theme(plot.caption = element_text(hjust = 0.5, size = 22)) +
  theme(legend.position = "none")+
  #theme(legend.position = "bottom")+
  theme(axis.title.x = element_blank()) +
  geom_hline(aes(yintercept = 1), color = "red", linetype = "dashed", linewidth = 1.5)
```

# Plot 2

``` r
my_colors <- c("darkblue", "lightblue")

p2 <- ggplot(tab_2, aes(x = Exposures, y = Est, color = Model)) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), position = position_dodge(0.7), width = 0.2, linewidth = 1.5) +
  geom_point(position = position_dodge(0.7), size = 3) +
  geom_text(aes(label = Metabolites), hjust = -0.2, vjust = -0.3, size = 7, position = position_dodge(0.7)) +
  scale_color_manual(values = my_colors) +
  theme_bw(base_size = 20) +
  ylab("Estimate") +
  theme(plot.caption = element_text(hjust = 0.5, size = 22)) +
  #theme(legend.position = "bottom")+
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  geom_hline(aes(yintercept = 0), color = "red", linetype = "dashed", linewidth = 1.5)
```

# Combine Plots

``` r
combined_plot <- ggarrange(
  p1, p2, 
  widths = c(1, 3), 
  ncol = 2, 
  common.legend = TRUE, 
  legend = "bottom"
)
```
