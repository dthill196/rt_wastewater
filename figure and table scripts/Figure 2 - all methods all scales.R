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

# FIGURE 2 - RT METHODS DURING A SPECIFIC SURGE, PANEL WITH METHODS AND CASES / CASE RT


plot_all_models <- 
  ggplot()+
  geom_line(data = rt_state_all %>% filter(method != "EpiSewer method"), aes(x= week, y = mean_rt, color = method), linewidth = 1)+
  geom_ribbon(data = rt_state_all %>% filter(method != "EpiSewer method"), aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt, fill = method), alpha = 0.2) +
  geom_hline(yintercept = 1, lty = "dashed") +
  labs(title = expression("R"["t"], "methods comparison"),
       x = "Date",
       y = "Rt")+
  scale_colour_manual(
    values =values)+
  scale_fill_manual(values = values)+
  theme_bw()+
  guides(color = guide_legend(override.aes = list(linewidth= 2)))+
  # scale_x_date(labels = date_format("%b %y"),
  #              date_breaks = "1 month", limits = c(as.Date("2023-10-01"), as.Date("2024-02-20")))+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month", limits = c(start_date, end_date))+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "bottom",
        legend.title = element_blank())

plot_all_models

# panel for all of them
figure_rt_all <- 
  cowplot::plot_grid(incidence_plot,  ww_plot, plot_all_models,rel_heights = c(1, 1, 1.5),
                     ncol = 1, align = "v",
                     labels = c("A)", "B)", "C)"), hjust = -2)
figure_rt_all

# export
png("Figure - all methods all weeks.png",
    units = "in",
    width = 8.5, height =10,
    res = 600)
figure_rt_all
dev.off()