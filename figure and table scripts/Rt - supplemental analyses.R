########
######## REPRODUCTION NUMBER FROM WASTEWATER - METHODS COMPARISON
######## SUPPLEMENTARY ANALYSES
########


# PACKAGES
# data wrangling
library(dplyr)
library(zoo)
library(tidyr)
library(lubridate)

# analysis
library(EpiEstim)
library(rstatix) # outlier detection
library(quantmod) # find peaks
library(broom)

# plotting
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(cowplot)

# export tables
library(sjPlot)

# DATA
dat <- readRDS("data/Rt_data_county.rds")
dat_state <- readRDS("data/Rt_data_state.rds")

# color palette
pal <- MetBrewer::met.brewer(n = 8, name = "Cross", type = "discrete")

# match to methods
values = c(pal, c("EpiEstim Substitution", 
                  "EpiSewer method",
                  "Exp. change rate",
                  "Fit line", 
                  "Goldstein method",
                  "Huisman method",
                  "Nourbakhsh method",
                  "Rolling GLM"))

# run comparisons with the fit line method since it is currently achieving the greatest accuracy relative to the reference model

# Reference model - case based Rt
source("functions/Rt-EpiEstim-functions.R")

county_select <- "Orange"

# example county df
cases_df <- dat %>%
  ungroup() %>%
  filter(county == county_select)%>%
  select(date, cases_new.7avg) %>%
  rename(Date = date) 

# Rt from  case data with SI 4 (not much difference in the unknown si results)
rt_cases_weekly <- rt_function(cases_df, mean_si = 4, std_si = 1, weekly = "Yes")

# --------------------------------------------------------------------------------------------------------------------------------------

# # # # # # # # # # # # # #
# lagged v. unlagged data #
# # # # # # # # # # # # # #

# Method 1 - fit line
source("functions/Rt-fit-line-function.R")

#extract ww data and join it to rt data
ww_data <- dat %>%
  select(date, county, sars2.7avg, start_date, intensity.7avg)%>%
  group_by(county)%>%
  arrange(date) %>%
  mutate(lag_sars2.7avg = lag(sars2.7avg, 4, default = NA),
         lag_intensity.7avg = lag(intensity.7avg, 4, default = NA)
  )%>%
  ungroup() %>%
  filter(county == county_select)  %>%
  filter(date >= start_date) %>%
  filter(!is.na(lag_sars2.7avg)) %>%
  select(-start_date)

# one dataframe for the county with ww data and rt from case data
data_fit_line <- inner_join(ww_data, rt_cases_weekly, c("date")) %>%
  rename(rt = mean_rt)

# run function
rt_fit_weekly_lag <- rt_fit_line_function(dataframe = data_fit_line, weekly = "Yes", range = 45, predictor = "lag_sars2.7avg")

# unlagged data
rt_fit_weekly_nolag <- rt_fit_line_function(dataframe = data_fit_line, weekly = "Yes", range = 45, predictor = "sars2.7avg")

# plot side by side
supplement_fig_1 <- 
  ggplot()+
  geom_line(data = rt_cases_weekly, aes(x = week, y = mean_rt, color = "#004f63"), linewidth = 1)+
  geom_line(data = rt_fit_weekly_lag, aes(x = week, y = mean_rt, color = "#ce4441"), linewidth = 1)+
  geom_line(data = rt_fit_weekly_nolag, aes(x = week, y = mean_rt, color = "#eb7926"), linewidth = 1)+
  theme_bw() +
  labs(title = expression("R"[t]~" fit line method: Lagged v. unlagged data"),
       x = "",
       y = expression("R"[t]))+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")+
  scale_colour_manual(name = "",
                      values =c(
                                '#ce4441'='#ce4441',
                                '#eb7926' = '#eb7926',
                                '#004f63' = '#004f63'), labels = c(expression("Case R"[t]), "Lagged data","No lag"
                                ))+
  geom_hline(yintercept = 1, linetype = "dashed")
supplement_fig_1

# save
png("Supplemental Figure 1 - Lag v. no lag.png",
     units = "in",
     width = 7, height =5,
     res = 600)
supplement_fig_1
dev.off()

# compare evaluation metrics for each of the fit line models to the reference method
source("functions/Rt-evaluation-function.R")

lag_eval <- rt_evaluation_function(method_name = "Fit line - lagged", reference_dataframe = rt_cases_weekly, method_dataframe = rt_fit_weekly_lag)
nolag_eval <- rt_evaluation_function(method_name = "Fit line - no lag", reference_dataframe = rt_cases_weekly, method_dataframe = rt_fit_weekly_nolag)

# merge into one df
metrics_table <- lag_eval %>%
  left_join(nolag_eval, by = c("metric_names")) 

metrics_table$county <- county_select

# not much difference, by the non-lagged data is the winner here. See what Dave thinks about this one
tab_df(metrics_table,
       file = "Supplemental Table 1 - Lag v. no lag.doc",
       title = "Supplemental Table 1: Lagged data v. no lagged data for Rt calculation. There is a small difference in the values when the wasteater data are 
       lagged 4 days. The unlagged data is closest to the reference model for the example county.")


# lag evaluation - should we lag our data (regardless of other people's data, not a method specific decision)

# lag script
## are the lags influenced by low pop counties ?

data.county <- dat

# empty list
datalist_lags <- list()

# number of lags to test
c <- 1:10

# for loop for daily lags
for(i in c){
  # create one lag for daily data and test correlation
  data.county <- data.county %>%
    group_by(county) %>%
    arrange(-desc(date)) %>%
    mutate(sars2.7avg_lag = lag(sars2.7avg, i, default = NA)) %>%
    ungroup()
  
  # correlation
  cor_result <- cor(log(data.county$sars2.7avg_lag), log1p(data.county$cases_new.7avg), use = "na.or.complete")
  
  # group
  group <- "daily"
  
  # store lag
  lag_n <- i
  
  # dataframe of values
  df <- as.data.frame(cbind(group, lag_n, cor_result), col.names = c("group", "lag_n", "cor_result"))
  
  datalist_lags[[i]] <- df
}
lags <- do.call(rbind, datalist_lags)

# plot them
lag_plot <- 
  ggplot(data = lags, aes(x = as.numeric(lag_n), y = as.numeric(cor_result), color = group, fill = group))+
  geom_point(color = "black", size = 1)+
  geom_line(color = "black", linewidth = 0.5)+
  theme_bw() +
  labs(title = "Pearson correlation between lagged wastewater data and new cases",
       x = "Lag",
       y = "Pearson's r")+
  theme(legend.position = "none")

lag_plot

# save
png("Supplemental Figure 2 - Lag correlation.png",
     units = "in",
     width = 7, height =5,
     res = 600)
lag_plot
dev.off()

# --------------------------------------------------------------------------------------------------------------------------------------

# # # # # # # # # # # # # # # #
# intensity v. raw gene copies#
# # # # # # # # # # # # # # # #

# Method 1 - fit line
source("functions/Rt-fit-line-function.R")

#extract ww data and join it to rt data
ww_data <- dat %>%
  select(date, county, sars2.7avg, start_date, intensity.7avg)%>%
  group_by(county)%>%
  arrange(date) %>%
  mutate(lag_sars2.7avg = lag(sars2.7avg, 4, default = NA),
         lag_intensity.7avg = lag(intensity.7avg, 4, default = NA)
  )%>%
  ungroup() %>%
  filter(county == county_select)  %>%
  filter(date >= start_date) %>%
  filter(!is.na(lag_sars2.7avg)) %>%
  select(-start_date)

# one dataframe for the county with ww data and rt from case data
data_fit_line <- inner_join(ww_data, rt_cases_weekly, c("date")) %>%
  rename(rt = mean_rt)

# run function
rt_fit_weekly_rawcopies <- rt_fit_line_function(dataframe = data_fit_line, weekly = "Yes", range = 45, predictor = "sars2.7avg")

# unlagged data
rt_fit_weekly_intensity <- rt_fit_line_function(dataframe = data_fit_line, weekly = "Yes", range = 45, predictor = "intensity.7avg")

# plot side by side
supplement_fig_3 <- 
  ggplot()+
  geom_line(data = rt_cases_weekly, aes(x = week, y = mean_rt, color = "#004f63"), linewidth = 1)+
  geom_line(data = rt_fit_weekly_rawcopies, aes(x = week, y = mean_rt, color = "#eb7926"), linewidth = 1)+
  geom_line(data = rt_fit_weekly_intensity, aes(x = week, y = mean_rt, color = "#ce4441"), linewidth = 1)+
  theme_bw() +
  labs(title = expression("R"[t]~" fit line method: Normalized v. not normalized data"),
       x = "",
       y = expression("R"[t]))+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")+
  scale_colour_manual(name = "",
                      values =c(
                        '#ce4441'='#ce4441',
                        '#eb7926' = '#eb7926',
                        '#004f63' = '#004f63'), labels = c( expression("Case R"[t]),  "Normalized data","Raw copies"
                        ))+
  geom_hline(yintercept = 1, linetype = "dashed")
supplement_fig_3

# save
png("Supplemental Figure 3 - Normalized v. not normalized.png",
     units = "in",
     width = 7, height =5,
     res = 600)
supplement_fig_3
dev.off()

# eval metrics
rawcopies_eval <- rt_evaluation_function(method_name = "Fit line - raw copies", reference_dataframe = rt_cases_weekly, method_dataframe = rt_fit_weekly_rawcopies)
intensity_eval <- rt_evaluation_function(method_name = "Fit line - normalized", reference_dataframe = rt_cases_weekly, method_dataframe = rt_fit_weekly_intensity)

# merge into one df
metrics_table <- rawcopies_eval %>%
  left_join(intensity_eval, by = c("metric_names")) 

metrics_table$county <- county_select

# not much difference, by the non-lagged data is the winner here. See what Dave thinks about this one
tab_df(metrics_table,
       file = "Supplemental Table 2 - normalized v. not normalized.doc",
       title = "Supplemental Table 2: Normalized data v. not normalized data for Rt calculation. There is a small difference in the values when the wasteater data are 
       lagged 4 days. Neither method is necessarily superior.")

# --------------------------------------------------------------------------------------------------------------------------------------

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# nyc data in the state model v. not in the state model #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# --------------------------------------------------------------------------------------------------------------------------------------

dat_state <- readRDS("data/Rt_data_state.rds")

dat_state_nonyc <- readRDS("data/Rt_data_state_noNYC.rds")


# reference method for state data

# Reference model - case based Rt
source("functions/Rt-EpiEstim-functions.R")

# including nyc reference data
cases_df <- dat_state %>%
  ungroup() %>%
  select(date, state_new_cases_7avg) %>%
  rename(Date = date) 

# Rt from  case data with SI 4 (not much difference in the unknown si results)
rt_cases_weekly_nyc <- rt_function(cases_df, mean_si = 4, std_si = 1, weekly = "Yes")

# example county df
cases_df <- dat_state_nonyc %>%
  ungroup() %>%
  select(date, state_new_cases_7avg) %>%
  rename(Date = date) 

# Rt from  case data with SI 4 (not much difference in the unknown si results)
rt_cases_weekly_nonyc <- rt_function(cases_df, mean_si = 4, std_si = 1, weekly = "Yes")


# Use the fit line method, for now, for the comparison
source("functions/Rt-fit-line-function.R")

#extract ww data and join it to rt data
ww_data <- dat_state %>%
  select(date, state_sars2)%>%
  arrange(date) %>%
  mutate(lag_sars2.7avg = lag(state_sars2, 4, default = NA)
  )%>%
  ungroup() %>%
  filter(!is.na(lag_sars2.7avg)) 

# one dataframe for the county with ww data and rt from case data
data_fit_line <- inner_join(ww_data, rt_cases_weekly_nyc, c("date")) %>%
  rename(rt = mean_rt)

# run function
rt_fit_weekly_nyc <- rt_fit_line_function(dataframe = data_fit_line, weekly = "Yes", range = 45, predictor = "lag_sars2.7avg")

# no nyc

#extract ww data and join it to rt data
ww_data <- dat_state_nonyc %>%
  select(date, state_sars2)%>%
  arrange(date) %>%
  mutate(lag_sars2.7avg = lag(state_sars2, 4, default = NA)
  )%>%
  ungroup() %>%
  filter(!is.na(lag_sars2.7avg)) 

# one dataframe for the county with ww data and rt from case data
data_fit_line <- inner_join(ww_data, rt_cases_weekly_nonyc, c("date")) %>%
  rename(rt = mean_rt)

# run function
rt_fit_weekly_nonyc <- rt_fit_line_function(dataframe = data_fit_line, weekly = "Yes", range = 45, predictor = "lag_sars2.7avg")

# comparison Figures

# Reference data figure first
# plot side by side
supplement_fig_4 <- 
  ggplot()+
  geom_line(data = rt_cases_weekly_nyc, aes(x = week, y = mean_rt, color = "#ce4441"), linewidth = 1)+
  geom_line(data = rt_cases_weekly_nonyc, aes(x = week, y = mean_rt, color = "#004f63"), linewidth = 1)+
  theme_bw() +
  labs(title = expression("R"[t]~" NYC v. no NYC data (case data)"),
       x = "",
       y = expression("R"[t]))+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")+
  scale_colour_manual(name = "",
                      values =c(
                        '#ce4441'='#ce4441',
                        '#004f63' = '#004f63'), labels = c("No NYC data", "NYC"
                        ))+
  geom_hline(yintercept = 1, linetype = "dashed")
supplement_fig_4


# model comparison
supplement_fig_5 <- 
  ggplot()+
  geom_line(data = rt_fit_weekly_nyc, aes(x = week, y = mean_rt, color = "#ce4441"), linewidth = 1)+
  geom_line(data = rt_fit_weekly_nonyc, aes(x = week, y = mean_rt, color = "#004f63"), linewidth = 1)+
  theme_bw() +
  labs(title = expression("R"[t]~" NYC v. no NYC data (Wastewater fit line method)"),
       x = "",
       y = expression("R"[t]))+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")+
  scale_colour_manual(name = "",
                      values =c(
                        '#ce4441'='#ce4441',
                        '#004f63' = '#004f63'), labels = c("No NYC data", "NYC"
                        ))+
  geom_hline(yintercept = 1, linetype = "dashed")

supplement_fig_5

grid.arrange(supplement_fig_4, supplement_fig_5, nrow = 2)

# save
png("Supplemental Figure 4 - NYC v. no NYC.png",
     units = "in",
     width = 7, height =10,
     res = 600)
grid.arrange(supplement_fig_4, supplement_fig_5, nrow = 2)
dev.off()


# -----------------------------------------------------------------------------------------------------------------------------

# ERN METHOD 
# Omicron or no omicron comparison
# ern data
ern_table <- readRDS("data/ern_county_metrics.rds")
ern_table_no_omicron <- readRDS("data/ern_county_metrics_no_omicron.rds")

# compare if omicron influenced ern results or no
ern_table_comparison <- ern_table %>%
  rename(ern_omicron = `ERN method`)%>%
  left_join(ern_table_no_omicron, by = c("metric_names", "county")
  ) %>%
  rename(ern_no_omicron = `ERN method`)

# plot to compare values across metrics
ggplot()+
  geom_point(data = ern_table_comparison, aes(x = as.numeric(ern_omicron), y = as.numeric(ern_no_omicron)))+
  facet_wrap(~metric_names, scales = "free")

# calculate difference in values, see what counties had the highest differences
ern_table_comparison$difference <- abs(as.numeric(ern_table_comparison$ern_omicron) - as.numeric(ern_table_comparison$ern_no_omicron))
summary(ern_table_comparison$difference)

# the inclusion of omicron seems to have a higher influence in more highly populated counties (NYC, Long Island) but also in some
# small counties (Schuyler). This is only for some metrics, not all for each county. For example, the pearson correlation was
# highly impacted by the Omicron wave for Bronx but the other values are not very different

# see how many metrics are 1 sd higher than the mean difference for that metric
metric_sd <- ern_table_comparison %>%
  group_by(metric_names) %>%
  summarize(mean= mean(as.numeric(difference)),
            sd = sd(as.numeric(difference))
  ) %>%
  ungroup() %>%
  mutate(mean_1sd = mean + sd)

ern_table_comparison <- left_join(ern_table_comparison, metric_sd, by = c("metric_names"))
ern_table_comparison$difference_high <- ifelse(ern_table_comparison$difference > ern_table_comparison$mean_1sd,
                                               "Yes", "No")

# any metrics that are off the most?
ern_comparison_summary <- ern_table_comparison %>%
  group_by(metric_names, difference_high) %>%
  summarize(n = n()) %>%
  ungroup()

# is there a county that appears on multiple metric lists
ern_comparison_yes <- ern_table_comparison %>%
  filter(difference_high == "Yes")
# Madison, Orange, New York, Bronx all appear at least 3 times for 3 metrics
# all of them have omicron data, the ones that are 0 are counties with no omicron data
# so, inclusion of omicron does alter the final results quite a bit

# ---------------------------------------------------------------------------------------------------------------------------------------------

# EPIESTIM - unknown v. known SI

dat_county <- readRDS("data/Rt_data_county.rds")
dat_state <- readRDS("data/Rt_data_state.rds")

source("functions/Rt-EpiEstim-functions.R")

county_select <- "Albany"

# example county df
cases_county_df <- dat_county%>%
  ungroup() %>%
  filter(county == county_select)%>%
  select(date, cases_new.7avg) %>%
  rename(Date = date) 

# Rt from  case data with SI 4 (not much difference in the unknown si results)
rt_cases_weekly_county <- rt_function(cases_county_df, mean_si = 4, std_si = 1, weekly = "Yes")
rt_cases_daily_county <- rt_function(cases_county_df, mean_si = 4, std_si = 1, weekly = "No")

# unknown serial interval
rt_cases_weekly_USI_county <- rt_function_unknown_si(cases_county_df,  weekly = "Yes")
rt_cases_daily_USI_county <- rt_function_unknown_si(cases_county_df,  weekly = "No")

# plot next to each other to comare results -> are they substantially different? 

# new package
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

daily_rt_plot_county <- ggplot()+
  geom_line(data = rt_cases_daily_county, aes(x = date, y = mean_rt, color = "#ce4441"), linewidth = 1) +
  geom_ribbon(data = rt_cases_daily_county, aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt,fill = "#ce4441"), alpha = 0.5)+
  geom_line(data = rt_cases_daily_USI_county, aes(x = date, y = mean_rt, color = "#004f63"), linewidth = 1)+
  #geom_ribbon(data = rt_cases_daily_USI, aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt, fill = "#004f63"), alpha = 0.5)+
  scale_colour_manual(name = "",
                      values =c(
                        '#ce4441'='#ce4441',
                        '#004f63' = '#004f63'), labels = c("Data-derived SI", "Mean SI of 4, SD of 1"
                        ))+
  scale_fill_manual(name = "",
                    values =c(
                      '#ce4441'='#ce4441',
                      '#004f63' = '#004f63'), labels = c("Data-derived SI", "Mean SI of 4, SD of 1"
                      ))+
  theme_bw()+
  labs(x = "",
       y = expression("R"[t]),
       title = expression("Albany County Daily R"[t]))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme(legend.position = "bottom")+
  geom_ribbon_pattern(data = rt_cases_daily_USI_county %>% filter(!is.na(mean_rt)), aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt, fill = "#004f63"),
                      color = "#004f63", 
                      pattern = "stripe",
                      pattern_fill = "gray80",
                      pattern_angle = 45,
                      pattern_density = 0.2,
                      pattern_spacing = 0.05,
                      pattern_alpha = 0.2,
                      alpha = 0.2) 

daily_rt_plot_county

# which one follows the case curve more closely?
weekly_rt_plot_county <- ggplot()+
  geom_line(data = rt_cases_weekly_county, aes(x = date, y = mean_rt, color = "#ce4441"), linewidth = 1) +
  geom_ribbon(data = rt_cases_weekly_county, aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt,fill = "#ce4441"), alpha = 0.5)+
  geom_line(data = rt_cases_weekly_USI_county, aes(x = date, y = mean_rt, color = "#004f63"), linewidth = 1)+
  #geom_ribbon(data = rt_cases_weekly_USI, aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt, fill = "#004f63"), alpha = 0.5)+
  scale_colour_manual(name = "",
                      values =c(
                        '#ce4441'='#ce4441',
                        '#004f63' = '#004f63'), labels = c("Data-derived SI", "Mean SI of 4, SD of 1"
                        ))+
  scale_fill_manual(name = "",
                    values =c(
                      '#ce4441'='#ce4441',
                      '#004f63' = '#004f63'), labels = c("Data-derived SI", "Mean SI of 4, SD of 1"
                      ))+
  theme_bw()+
  labs(x = "",
       y = expression("R"[t]),
       title = expression("Albany County Weekly R"[t]))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme(legend.position = "bottom")+
  geom_ribbon_pattern(data = rt_cases_weekly_USI_county %>% filter(!is.na(mean_rt)), aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt, fill = "#004f63"),
                      color = "#004f63", 
                      pattern = "stripe",
                      pattern_fill = "gray80",
                      pattern_angle = 45,
                      pattern_density = 0.2,
                      pattern_spacing = 0.05,
                      pattern_alpha = 0.2,
                      alpha = 0.2) 

weekly_rt_plot_county

grid.arrange(daily_rt_plot_county, weekly_rt_plot_county, nrow = 2)

# statewide
cases_state_df <- dat_state%>%
  ungroup() %>%
  select(date, state_new_cases_7avg) %>%
  rename(Date = date) 

# Rt from  case data with SI 4 (not much difference in the unknown si results)
rt_cases_weekly_state <- rt_function(cases_state_df, mean_si = 4, std_si = 1, weekly = "Yes")
rt_cases_daily_state <- rt_function(cases_state_df, mean_si = 4, std_si = 1, weekly = "No")

# unknown serial interval
rt_cases_weekly_USI_state <- rt_function_unknown_si(cases_state_df,  weekly = "Yes")
rt_cases_daily_USI_state <- rt_function_unknown_si(cases_state_df,  weekly = "No")

# plot next to each other to comare results -> are they substantially different? 
daily_rt_plot_state <- ggplot()+
  geom_line(data = rt_cases_daily_state, aes(x = date, y = mean_rt, color = "#ce4441"), linewidth = 1) +
  geom_ribbon(data = rt_cases_daily_state, aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt,fill = "#ce4441"), alpha = 0.5)+
  geom_line(data = rt_cases_daily_USI_state, aes(x = date, y = mean_rt, color = "#004f63"), linewidth = 1)+
  scale_colour_manual(name = "",
                      values =c(
                        '#ce4441'='#ce4441',
                        '#004f63' = '#004f63'), labels = c("Data-derived SI", "Mean SI of 4, SD of 1"
                        ))+
  scale_fill_manual(name = "",
                    values =c(
                      '#ce4441'='#ce4441',
                      '#004f63' = '#004f63'), labels = c("Data-derived SI", "Mean SI of 4, SD of 1"
                      ))+
  theme_bw()+
  labs(x = "",
       y = expression("R"[t]),
       title = expression("NY Statewide Daily R"[t]))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme(legend.position = "bottom")+
  geom_ribbon_pattern(data = rt_cases_daily_USI_state %>% filter(!is.na(mean_rt)), aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt, fill = "#004f63"),
                      color = "#004f63", 
                      pattern = "stripe",
                      pattern_fill = "gray80",
                      pattern_angle = 45,
                      pattern_density = 0.2,
                      pattern_spacing = 0.05,
                      pattern_alpha = 0.2,
                      alpha = 0.2) 

daily_rt_plot_state

# which one follows the case curve more closely?
weekly_rt_plot_state <- ggplot()+
  geom_line(data = rt_cases_weekly_state, aes(x = date, y = mean_rt, color = "#ce4441"), linewidth = 1) +
  geom_ribbon(data = rt_cases_weekly_state, aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt,fill = "#ce4441"), alpha = 0.5)+
  geom_line(data = rt_cases_weekly_USI_state, aes(x = date, y = mean_rt, color = "#004f63"), linewidth = 1)+
  scale_colour_manual(name = "",
                      values =c(
                        '#ce4441'='#ce4441',
                        '#004f63' = '#004f63'), labels = c("Data-derived SI", "Mean SI of 4, SD of 1"
                        ))+
  scale_fill_manual(name = "",
                    values =c(
                      '#ce4441'='#ce4441',
                      '#004f63' = '#004f63'), labels = c("Data-derived SI", "Mean SI of 4, SD of 1"
                      ))+
  theme_bw()+
  labs(x = "",
       y = expression("R"[t]),
       title = expression("NY Statewide Weekly R"[t]))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme(legend.position = "bottom")+
  geom_ribbon_pattern(data = rt_cases_weekly_USI_state %>% filter(!is.na(mean_rt)), aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt, fill = "#004f63"),
                      color = "#004f63", 
                      pattern = "stripe",
                      pattern_fill = "gray80",
                      pattern_angle = 45,
                      pattern_density = 0.2,
                      pattern_spacing = 0.05,
                      pattern_alpha = 0.2,
                      alpha = 0.2) 

weekly_rt_plot_state

grid.arrange(daily_rt_plot_state, weekly_rt_plot_state, nrow = 2)

grid.arrange( daily_rt_plot_county, weekly_rt_plot_county, 
  daily_rt_plot_state, weekly_rt_plot_state, nrow = 2)

# save
png("Supplemental Figure 5 - Known SI v. unknown.png",
    units = "in",
    width = 11, height =10,
    res = 600)
grid.arrange( daily_rt_plot_county, weekly_rt_plot_county, 
              daily_rt_plot_state, weekly_rt_plot_state, nrow = 2)
dev.off()

# --------------------------------------------------------------------------------------------------------------------------------------

# Extra smoothing

ern_state <- readRDS("data/ern_state_weekly.rds")

ern_state<- ern_state %>%
  mutate(center_3_mean = rollmean(mean_rt, 3, align = "center", na.pad = TRUE),
         center_5_mean = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         center_7_mean = rollmean(mean_rt, 7, align = "center", na.pad = TRUE))

ggplot()+
  geom_line(data = ern_state, aes(x = week, y = center_3_mean), color = "seagreen", linewidth = 1)+
  geom_line(data = ern_state, aes(x = week, y = mean_rt), color = "dodgerblue", linewidth = 1)+
  geom_line(data = ern_state, aes(x = week, y = center_5_mean), color = "gold", linewidth = 1)+
  labs(title = "3 week mean v. standard output",
       x = "",
       y = expression("R"[t]))+
  theme_bw()

p1 <- ggplot()+
  geom_line(data = rt_cases_weekly, aes(x = week, y = mean_rt), color = "orange", linewidth = 1)+
  geom_line(data = ern_state, aes(x = week, y = mean_rt), color = "dodgerblue", linewidth = 1)+
  labs(title = "standard output v. cases",
       x = "",
       y = expression("R"[t]))+
  theme_bw()

p2 <- ggplot()+
  geom_line(data = ern_state, aes(x = week, y = center_3_mean), color = "seagreen", linewidth = 1)+
  geom_line(data = rt_cases_weekly, aes(x = week, y = mean_rt), color = "orange", linewidth = 1)+
  labs(title = "3 week mean",
       x = "",
       y = expression("R"[t]))+
  theme_bw()

p3 <- ggplot()+
  geom_line(data = ern_state, aes(x = week, y = center_5_mean), color = "seagreen", linewidth = 1)+
  geom_line(data = rt_cases_weekly, aes(x = week, y = mean_rt), color = "orange", linewidth = 1)+
  labs(title = "5 week mean",
       x = "",
       y = expression("R"[t]))+
  theme_bw()

p4 <- ggplot()+
  geom_line(data = ern_state, aes(x = week, y = center_7_mean), color = "gold", linewidth = 1)+
  geom_line(data = rt_cases_weekly, aes(x = week, y = mean_rt), color = "orange", linewidth = 1)+
  labs(title = "7 week mean",
       x = "",
       y = expression("R"[t]))+
  theme_bw()
grid.arrange(p1, p2, p3, p4)

# save
png("Supplemental Figure 6 - Additional smoothing.png",
    units = "in",
    width = 10, height =7,
    res = 600)
grid.arrange(p1, p2, p3, p4)
dev.off()


# load the data
huisman_rt_smooth <- readRDS("data/huisman_rt_smooth.rds")
rt_huisman_county <- readRDS("data/huisman_rt_county.rds")
rt_exp_county_smooth <- readRDS("data/rt_exp_county_smooth.rds")
rt_exp_county <- readRDS("data/rt_exp_county.rds")
rt_sub_county <- readRDS("data/rt_sub_county.rds")
rt_sub_county_smooth <- readRDS("data/rt_sub_county_smooth.rds")
rt_glm_county <- readRDS("data/rt_glm_county.rds")
rt_glm_county_smooth <- readRDS("data/rt_glm_county_smooth.rds")
rt_fit_county <- readRDS("data/rt_fit_county.rds")
rt_cases_weekly_county <- readRDS("data/rt_cases_weekly_county.rds")

# evaluate each one by county
source("functions/Rt-evaluation-function.R")


# ---------------------------------------------------------------------------------------------------------------------

# compare evaluations of smooth v. not smooth -> is one better? 
# huisman smooth v. not smooth

# empty list for loop
eval_list <- list()

for(i in unique(huisman_rt_smooth$county)){
  
  # filter dfs for the county 
  eval_df <- huisman_rt_smooth %>%
    filter(county == i)
  reference_df <- rt_cases_weekly_county %>%
    filter(county == i)
  
  eval_table <- rt_evaluation_function(method_name = "Huisman - smooth", reference_dataframe = reference_df, method_dataframe = eval_df)
  eval_table$county <- i
  # store in list
  eval_list[[i]] <- eval_table
  
}

# make the list a df
huisman_smooth_eval <- do.call(rbind, eval_list)

# not smoothed
# empty list for loop
eval_list <- list()

for(i in unique(rt_huisman_county$county)){
  
  # filter dfs for the county 
  eval_df <- rt_huisman_county%>%
    filter(county == i)
  reference_df <- rt_cases_weekly_county %>%
    filter(county == i)
  
  eval_table <- rt_evaluation_function(method_name = "Huisman - no smooth", reference_dataframe = reference_df, method_dataframe = eval_df)
  eval_table$county <- i
  # store in list
  eval_list[[i]] <- eval_table
  
}

# make the list a df
huisman_eval <- do.call(rbind, eval_list)
huisman_eval$group <- "Huisman no smooth"
huisman_smooth_eval$group <- "Huismain smooth"
huisman_eval <- huisman_eval %>%
  rename(value = `Huisman - no smooth`)
huisman_smooth_eval <- huisman_smooth_eval %>%
  rename(value = `Huisman - smooth`)
t <- bind_rows(huisman_eval, huisman_smooth_eval)



density_plot <- 
  ggplot(t )+
  geom_density(aes( y = as.numeric(value), group = group, color = group, fill = group), alpha = 0.5, 
               stat = "density",
               position = "identity")+
  coord_flip()+
  theme_bw()+
  facet_wrap(~metric_names, scales = "free")+
  theme(legend.position = "bottom")
density_plot

# save for supplement


# -------------------------------------------------------------------------------------------------------------------------------------

# Epi Estim - different waves v. continuous

dat <- readRDS("data/Rt_data_county.rds")

source("functions/Rt-EpiEstim-functions.R")

county_select <- "Albany"

# example county df
cases_df <- dat %>%
  ungroup() %>%
  filter(county == county_select)%>%
  select(date, cases_new.7avg) %>%
  rename(Date = date) 

# Rt from  case data with SI 4 (not much difference in the unknown si results)
rt_cases_weekly <- rt_function(cases_df, mean_si = 4, std_si = 1, weekly = "Yes")

# waves -> chunk out each season separately
wave_1 <- dat %>%
  ungroup() %>%
  filter(county == county_select)%>%
  select(date, cases_new.7avg) %>%
  rename(Date = date) %>%
  filter(Date <= as.Date("2023-06-01"))


wave_2 <- dat %>%
  ungroup() %>%
  filter(county == county_select)%>%
  select(date, cases_new.7avg) %>%
  rename(Date = date) %>%
  filter(Date > as.Date("2023-06-01"))

rt_wave_1 <- rt_function(wave_1, mean_si = 4, std_si = 1, weekly = "Yes")
rt_wave_2 <- rt_function(wave_2, mean_si = 4, std_si = 1, weekly = "Yes")

plot1 <- 
  ggplot()+
  geom_line(data = rt_cases_weekly, aes(x = date, y = mean_rt, color = "#ce4441"), linewidth = 1) +
  geom_ribbon(data = rt_cases_weekly, aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt), alpha = 0.5,fill = "#ce4441")+
  geom_line(data =rt_wave_1, aes(x = date, y = mean_rt, color = "#004f63"), linewidth = 1)+
  geom_line(data = rt_wave_2, aes(x = date, y = mean_rt, color = "gold"), linewidth = 1)+
  scale_color_manual(name = "",
                     values =c(
                       '#ce4441'='#ce4441',
                       '#004f63' = '#004f63',
                       'gold' = 'gold'), labels = c( "2022-2023","full time series", "2023-2024"
                       ))+
  theme_bw()+
  labs(title = "EpiEstim - broken time series v. continuous data",
       x = "",
       y = expression("R"[t]))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme(legend.position = "bottom")

# save
png("Supplemental Figure 7 - epiestim broken time series.png",
    units = "in",
    width = 7, height =5,
    res = 600)
plot1
dev.off()
  

##########################
# LOG V. NOT LOG GENE COPIES
##########################

dat_state <- readRDS("data/Rt_data_state.rds")

# statewide weekly
log_copies <- dat_state %>%
  ungroup() %>%
  select(date, state_new_cases_7avg) %>%
  mutate(log_copies = log(state_new_cases_7avg)) %>%
  select(-state_new_cases_7avg)

rt_state_log <- rt_function(dataframe = log_copies, mean_si = 4, std_si = 1, weekly = "Yes")

# not logged
not_log_copies <- dat_state %>%
  ungroup() %>%
  select(date, state_new_cases_7avg)

rt_state_not_log <- rt_function(dataframe = not_log_copies, mean_si = 4, std_si = 1, weekly = "Yes")

# plot them both
# which one follows the case curve more closely?
log_plot <- ggplot()+
  geom_line(data = rt_state_log, aes(x = date, y = mean_rt, color = "#ce4441"), linewidth = 1) +
  geom_ribbon(data = rt_state_log, aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt,fill = "#ce4441"), alpha = 0.5)+
  geom_line(data = rt_state_not_log, aes(x = date, y = mean_rt, color = "#004f63"), linewidth = 1)+
  scale_colour_manual(name = "",
                      values =c(
                        '#ce4441'='#ce4441',
                        '#004f63' = '#004f63'), labels = c("Not log transformed","Log transformed"
                        ))+
  scale_fill_manual(name = "",
                    values =c(
                      '#ce4441'='#ce4441',
                      '#004f63' = '#004f63'), labels = c("Not log transformed", "Log transformed"
                      ))+
  theme_bw()+
  labs(x = "",
       y = expression("R"[t]),
       title = expression("NY Statewide Weekly R"[t]))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  theme(legend.position = "bottom")+
  geom_ribbon_pattern(data = rt_state_not_log %>% filter(!is.na(mean_rt)), aes(x = date, ymin = ll_95_rt, ymax = ul_95_rt, fill = "#004f63"),
                      color = "#004f63", 
                      pattern = "stripe",
                      pattern_fill = "gray80",
                      pattern_angle = 45,
                      pattern_density = 0.2,
                      pattern_spacing = 0.05,
                      pattern_alpha = 0.2,
                      alpha = 0.2) 

log_plot

png("Supplemental Figure 8 - log v. not log transformed.png",
    units = "in",
    width = 7, height =5,
    res = 600)
log_plot
dev.off()


# # # # # # # #
# ADDITIONAL SMOOTHING OF WW DATA GOING INTO MODELS
# # # # # # # # 

# try 14 day and 30 day lags of ww data. does this improve the smoothness and the model output
dat_state_smooth <- dat_state %>%
  mutate(sars2_14 = rollmean(state_sars2, 14, default = NA, na.pad = TRUE, align = "right"),
         sars2_30 = rollmean(state_sars2, 30, default = NA, na.pad = TRUE, align = "right")
         )

# plot to compare
ggplot()+
  geom_line(data = dat_state_smooth, aes(x = date,y = state_sars2), linewidth = 1, color = "black")+
  geom_line(data = dat_state_smooth, aes(x = date, y = sars2_14), linewidth = 1, color = "gold")+
  geom_line(data = dat_state_smooth, aes(x = date, y = sars2_30), linewidth = 1, color = "dodgerblue")



# now use each in the rolling glm function 
source("functions/Rt-glm-roll.R")
dat_14 <- dat_state_smooth %>%
  select(date, state_new_cases_7avg, sars2_14) %>%
  filter(!is.na(sars2_14))
rt_ww_14 <- rt_glm_function(dataframe = dat_14, predictor = "sars2_14", case_data = "state_new_cases_7avg", range = 45, unknown_si = "No",
                            mean_si = 4, std_si = 1, weekly = "Yes")

# 30 day model
dat_30 <- dat_state_smooth %>%
  select(date, state_new_cases_7avg, sars2_30) %>%
  filter(!is.na(sars2_30))
rt_ww_30 <- rt_glm_function(dataframe = dat_30, predictor = "sars2_30", case_data = "state_new_cases_7avg", range = 45, unknown_si = "No",
                            mean_si = 4, std_si = 1, weekly = "Yes")

# original, 7 day average
rt_ww_7 <- rt_glm_function(dataframe = dat_state, predictor = "state_sars2", case_data = "state_new_cases_7avg", range = 45, unknown_si = "No",
                           mean_si = 4, std_si = 1, weekly = "Yes")

# evalute them against the case data from epiestim method
source("functions/Rt-EpiEstim-functions.R")
case_data <- dat_state %>%
  select(date, state_new_cases_7avg)
rt_cases <- rt_function(dataframe = case_data, weekly = "Yes", mean_si = 4, std_si = 1 )

# does it improve the evaluation metrics
source("functions/Rt-evaluation-function.R")

eval_14 <- rt_evaluation_function(method_name = "14 day smoothed average", reference_dataframe = rt_cases, method_dataframe = rt_ww_14)
eval_30 <- rt_evaluation_function(method_name = "30 day smoothed average", reference_dataframe = rt_cases, method_dataframe = rt_ww_30)
eval_7 <- rt_evaluation_function(method_name = "7 day smoothed average", reference_dataframe = rt_cases, method_dataframe = rt_ww_7)

# make a table
table_smooth <- left_join(eval_7, eval_14, by = c("metric_names"))
table_smooth <- left_join(table_smooth, eval_30, by = c("metric_names"))

# not definitive benefit based on this table, so leave it as 7 day average for now
# save table for supplement
tab_df(table_smooth, file = "supplemental table - smoothed average comparison.doc",
       digits = 3)
