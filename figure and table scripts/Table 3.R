# PACKAGES
# data wrangling
library(dplyr)
library(zoo)
library(tidyr)
library(lubridate)

# analysis
library(EpiEstim)
library(quantmod) # find peaks

# plotting
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(cowplot)
library(ggpubr)

# DATA

# case rt (EpiEstim)
rt_cases_weekly <- readRDS("data/rt_cases_weekly.rds") %>%
  mutate(method = "Rt cases")

# ww substitution
rt_ww_sub <- readRDS("data/rt_ww_sub.rds")
rt_ww_sub <- rt_ww_sub %>%
  arrange(date) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(method = "EpiEstim Substitution")

# EpiSewer method - statewide
episewer_state <- readRDS("data/episewer_state.rds") %>%
  arrange(week) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(method = "EpiSewer method")

# exp change rate
rt_exp_weekly <- readRDS("data/rt_exp_weekly.rds") %>%
  mutate(method = "Exp change rate")

# fit line
rt_fit_weekly <- readRDS("data/rt_fit_weekly.rds") %>%
  mutate(method = "Fit line")

# goldstein method statewide (note - need to rerun this with 200 resamples)
eirr_state <- readRDS("data/eirr_state.rds") %>%
  arrange(week) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(method = "Goldstein method")

# huisman method
rt_huisman_weekly <- readRDS("data/rt_huisman_weekly.rds") %>%
  arrange(week) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(method = "Huisman method")

# ern state
ern_state <- readRDS("data/ern_state_weekly.rds") %>%
  arrange(week) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(method = "ERN method")

# rolling GLM
rt_glm_weekly <- readRDS("data/rt_glm_weekly.rds") %>%
  mutate(method = "Rolling GLM")
# 5 week average
rt_glm_weekly <- rt_glm_weekly %>%
  arrange(date) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  )

# one large dataframe
rt_state_all <- bind_rows(
  rt_cases_weekly, rt_exp_weekly, rt_fit_weekly, rt_glm_weekly, eirr_state, ern_state, episewer_state, rt_huisman_weekly, rt_ww_sub
)
# TABLE - EVALUATION OF EACH METHOD AGAINST THE REFERENCE CASE RT METHOD

source("functions/Rt-evaluation-function.R")

# fit line
eval_fit_line <- rt_evaluation_function(method_name = "Fit line method", reference_dataframe = rt_cases_weekly, method_dataframe = rt_fit_weekly)

# rolling glm
eval_glm <- rt_evaluation_function(method_name = "Rolling GLM", reference_dataframe = rt_cases_weekly, method_dataframe = rt_glm_weekly)

# epiestim sub
eval_epi_sub <- rt_evaluation_function(method_name = "EpiEstim Substitution", reference_dataframe = rt_cases_weekly, method_dataframe = rt_ww_sub)

# exp change rate
eval_exp <- rt_evaluation_function(method_name = "Exp change rate", reference_dataframe = rt_cases_weekly, method_dataframe = rt_exp_weekly)

# huisman
eval_huisman <- rt_evaluation_function(method_name = "Huisman method", reference_dataframe = rt_cases_weekly, method_dataframe = rt_huisman_weekly)

# goldstein
eval_goldstein <- rt_evaluation_function(method_name = "Goldstein method", reference_dataframe = rt_cases_weekly, method_dataframe = eirr_state)

# episewer
eval_episewer <- rt_evaluation_function(method_name = "EpiSewer method", reference_dataframe = rt_cases_weekly, method_dataframe = episewer_state)

# ern
eval_ern <- rt_evaluation_function(method_name = "ERN method", reference_dataframe = rt_cases_weekly, method_dataframe = ern_state)

# merge together into one table and round digits
state_eval_table <- left_join(eval_fit_line, eval_glm, by = c("metric_names")) %>%
  left_join(eval_epi_sub, by = c("metric_names")) %>%
  left_join(eval_exp, by = c("metric_names")) %>%
  left_join(eval_huisman, by = c("metric_names")) %>%
  left_join(eval_goldstein, by = c("metric_names")) %>%
  left_join(eval_episewer, by = c("metric_names")) %>%
  left_join(eval_ern, by = c("metric_names"))
str(state_eval_table)

# change to numeric so we can round in our table
state_eval_table[2:9] <- lapply(state_eval_table[2:9], as.numeric)
str(state_eval_table)

# round column results
state_eval_table[,-1] <- round(state_eval_table[,-1], 4)

# TABLE - EVAL METRICS TABLE COUNTY

# Using median of county values ?
eval_county <- metrics_long %>%
  group_by(model, metric_names) %>%
  summarize(median_value = round(median(value, na.rm = TRUE),4)
  ) %>%
  ungroup()

# factor the values so we can order them
eval_county$metric_names <- factor(eval_county$metric_names,
                                   levels = c("RMSE",
                                              "Pearson Correlation",
                                              "Percent of peaks that coincide",
                                              "Above or below 1 percent agreement",
                                              "Mean absolute difference",
                                              "Sharpness"))
table(eval_county$metric_names)
eval_county$model <- factor(eval_county$model,
                            levels = c("Fit line method",
                                       "Rolling GLM",
                                       "EpiEstim Substitution", 
                                       "Exp change rate",
                                       "Huisman method",
                                       "Goldstein - EIRR",
                                       "EpiSewer",
                                       "ERN method"))

# sewershed level evaluation table and boxplot (once episewer finishes running) ...
# one large evaluation table
eval_sewershed <- readRDS("data/eval_sewershed.rds") %>%
  rename(value = median_value)
eval_sewershed$scale <- "Sewershed"

# county
eval_county$scale <- "County"
eval_county <- eval_county %>%
  rename(value = median_value)
eval_county$metric_names[eval_county$metric_names == "Above or below 1 percent agreement "] <- "Above or below 1 agreement"

# state
eval_state <- state_eval_table %>%
  pivot_longer(cols = c(`Fit line method`, `Rolling GLM`, `EpiEstim Substitution`, `Exp change rate`, `Huisman method`,
                        `Goldstein method`, `EpiSewer method`, `ERN method`)) %>%
  rename(model = name) %>%
  mutate(scale = "State")

# combine
eval_table_final <- bind_rows(eval_state, eval_county, eval_sewershed)
eval_table_final$metric_names[eval_table_final$metric_names == "Above or below 1 percent agreement"] <- "Above or below 1 agreement"

# try making a plot
ggplot(data = eval_table_final)+
  geom_point(aes(x = metric_names, y = value, shape = model))+
  facet_grid(~model)

# make it wider
eval_table_final_wide <- eval_table_final %>%
  pivot_wider(names_from = metric_names, values_from = value)

# save
tab_df(eval_table_final_wide %>% arrange(model),
       file = "E:Table - evaluation - all scales.doc",
       title = "Table : Evaluation of Rt methods.",
       digits = 4)