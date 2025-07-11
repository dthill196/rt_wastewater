Estimating the effective reproduction number from wastewater (Rt)
================

## Introduction

This repository is a collection of methods for estimating the effective
reprodouction number from wastewater (the Rt). The Rt is a number that
indicates the growth or decline of disease spread. Eight distinct
methods were identified from the literature and subsequently compared to
one another for accuracy, computation time, and efficiency. Each method
has been wrapped in a function with data prep provided. On this page,
you will find some descriptions of the methods and examples. Full
explanations and documentation of each Rt method are available in the
**functions** folder.

## Citation information

When using code or other content from this repository, please cite the
following:

Hill DT, Zhu Y, Dunham C, Moran J, Zhou Y, Collins MB, Kmush BL, Larsen DA. 
Estimating the effective reproduction number from wastewater (Rt): 
A methods comparison. Epidemics. 2025 Jun 18:100839.
[https://doi.org/10.1016/j.epidem.2025.100839](https://doi.org/10.1016/j.epidem.2025.100839)

## EpiEstim and Rt functions

Much of the work that has gone into the calculation of the Rt originates
with models and methods proposed by Cori et al (2013) and their
`EpiEstim` R package.

``` r
library(dplyr)
library(tidyr)
library(tidyverse)
library(EpiEstim)
library(ggplot2)
library(zoo)
library(scales)
library(broom)

# 
# DATA
dat <- readRDS("data/Rt_data_state.rds")

# Reference model - case based Rt
source("functions/Rt-EpiEstim-functions.R")

# example county df
cases_df <- dat %>%
  ungroup() %>%
  select(date, state_new_cases_7avg) %>%
  rename(Date = date) 

# Rt from  case data with SI 4 (not much difference in the unknown si results)
rt_cases_weekly <- rt_function(cases_df, mean_si = 4, std_si = 1, weekly = "Yes")

# start and end dates for time series
start_date <- as.Date("2022-09-01")
end_date <- as.Date("2024-02-20")

# plot
rt_case_plot <- 
  ggplot(data = rt_cases_weekly, aes(x = date, y = mean_rt))+
  geom_ribbon(aes(ymin = ll_95_rt, ymax = ul_95_rt), alpha = 0.4, fill = "orange") +
  geom_line(color = "orange", linewidth = 1) +
  geom_hline(yintercept = 1, lty = "dashed")+
  geom_vline(xintercept = as.Date("2022-12-21"), lty = "longdash")+
  geom_vline(xintercept = as.Date("2023-09-21"), lty = "longdash")+
  geom_vline(xintercept = as.Date("2024-01-05"), lty = "longdash")+
  theme_bw() +
  labs(title = expression("Statewide effective reproduction number (R"["t"]~") from case data"),
       x = "",
       y = expression("Mean R"["t"]~"cases"))+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month", limits = c(start_date, end_date))+
  theme(axis.text.x = element_text(angle = 90))
rt_case_plot
```

![](README_files/figure-gfm/EpiEstim%20example-1.png)<!-- -->

## Example wastewater-based Rt method

Eight wastewater methods were identified and are explained in detail in
the manuscript. Examples for running each method can be found in the
`functions` folder.

``` r
# Method 1 - fit line
source("functions/Rt-fit-line-function.R")

#extract ww data and join it to rt data
ww_data <- dat %>%
  select(date, state_sars2)%>%
  arrange(date) 

# one dataframe for the county with ww data and rt from case data
data_fit_line <- inner_join(ww_data, rt_cases_weekly, c("date")) %>%
  rename(rt = mean_rt)

# run function
rt_fit_weekly <- rt_fit_line_function(dataframe = data_fit_line, weekly = "Yes", range = 45, predictor = "state_sars2")
rt_cases_weekly$model <- "Rt cases"
rt_fit_weekly$model <- "Fit line"
values <- c("orange", "darkblue")

# plot
plot_fit_line <- 
  ggplot()+
  geom_line(data = rt_cases_weekly, aes(x = date, y = mean_rt, color = model), linewidth = 1) +
  geom_ribbon(data = rt_cases_weekly, aes(x = date, y = mean_rt, ymin = ll_95_rt, ymax = ul_95_rt, fill = model), alpha = 0.4) +
  geom_line(data = rt_fit_weekly, aes(x= week, y = mean_rt, color = model), linewidth = 1)+
  geom_ribbon(data = rt_fit_weekly, aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt, fill = model), alpha = 0.2) +
  geom_hline(yintercept = 1, lty = "dashed") +
  labs(title = expression("R"["t"]~"fit line"),
       x = "Date",
       y = "Rt")+
  scale_colour_manual(name = "",
                      values =values)+
  scale_fill_manual(name = "", values = values)+
  theme_bw()+
  guides(color = guide_legend(override.aes = list(linewidth= 2)))+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month", 
               #limits = c(as.Date("2023-10-01"), as.Date("2024-02-20"))
               limits = c(start_date, end_date)
               )+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom")
plot_fit_line
```

![](README_files/figure-gfm/ww%20example-1.png)<!-- -->

## Repository outline

### Data folder

This folder contains data for generating the figures and tables
published in the paper. Data are from New York State’s Wastewater
Surveillance Network [Open Data
website](https://coronavirus.health.ny.gov/covid-19-data-new-york).

### Functions folder

This folder contains a function for each Rt method as well as additional
functions used in the calculation of each Rt method. The functions are
wrappers for the methods and packages that we review in the paper.

### Figure and table scripts folder

This folder contains the R scripts for generating the figures and tables
included in the Hill et al. publication linked to this repository.

### High-throughput folder

This folder contains example scripts for setting up a parallel process
on separate workers with HT Condor and Conda. Consult your IT department
for assistance. These scripts are intended as an example and will not
work without adaptation to your specific computing environment.

#### Note on high-throughput computing

Six of the methods that we evaluated can be run on most machines,
however, two methods (`EpiSewer` and the Goldstein et al. method)
require additional computing resources. We used a high-throughput
computing approach through our institution’s computing cluster that was
based in HT Condor and used Conda to set up R and other environments. We
provide an example of that process in the **high-throughput** folder.
Depending on your institution, you may need to set up your own process.
This is only needed if you intend to iterate Rt estimates over multiple
separate sets of data. As an example, New York has 62 counties and we
wanted to compare the accuracy of the models across counties based on
population size. Thus, we needed 62 separate Rt estimates for each
method. A `for` loop and a function for the six other methods, however,
for `EpiSewer` and the Goldstein method, we needed to distribute each
county calculation as separate processes. Thus, the need for
high-throughput computing to set up 62 separate calculations.

If you are estimating the Rt for one jurisdiction, you will likely not
need high-througput computing. Note that each method varies in the
length of time it takes to run, and speed will be impacted by the data
structure. Please consult the individual package documentation and
papers associated with each method for more details and to troubleshoot
any issues.

## References and resources

### EpiEstim

Cori A, Ferguson NM, Fraser C, Cauchemez S. A New Framework and Software
to Estimate Time-Varying Reproduction Numbers During Epidemics. Am J
Epidemiol. 2013 Nov 1;178(9):1505–12.

Cori A. EpiEstim: Estimate Time Varying Reproduction Numbers from
Epidemic Curves \[Internet\]. 2021. Available from:
<https://CRAN.R-project.org/package=EpiEstim>

GitHub page: <https://github.com/mrc-ide/EpiEstim>

### Huisman et al. method

Huisman JS, Scire J, Caduff L, Fernandez -Cassi Xavier,
Ganesanandamoorthy P, Kull A, et al. Wastewater-Based Estimation of the
Effective Reproductive Number of SARS-CoV-2. Environ Health Perspect.
2022;130(5):057011.

GitHub page: <https://github.com/JSHuisman/wastewaterRe>

### Goldstein et al. method

Goldstein IH, Parker DM, Jiang S, Minin VM. Semiparametric inference of
effective reproduction number dynamics from wastewater pathogen
surveillance data. ArXiv. 2023 Aug 31;arXiv:2308.15770v2.

Isaac H Goldstein, Daniel M Parker, Sunny Jiang, Volodymyr M Minin,
Semiparametric inference of effective reproduction number dynamics from
wastewater pathogen surveillance data, Biometrics, Volume 80, Issue 3,
September 2024, ujae074, <https://doi.org/10.1093/biomtc/ujae074>

GitHub page: <https://github.com/igoldsteinh/concRt>

### EpiSewer method

GitHub page: <https://github.com/adrian-lison/EpiSewer>

### ERN method

Champredon D, Papst I, Yusuf W. ern: An R package to estimate the
effective reproduction number using clinical and wastewater surveillance
data. PLOS ONE. 2024 Jun 21;19(6):e0305550.

GitHub page: <https://github.com/phac-nml-phrsd/ern>
