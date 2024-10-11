
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

# statewide data
dat_state <- readRDS("data/Rt_data_state.rds")

# case rt (EpiEstim)
rt_cases_weekly <- readRDS("data/rt_cases_weekly.rds") %>%
  mutate(method = "Rt cases")

# FIGURE 1 - PANEL WITH CASES, HOSP, WW, AND Rt FROM CASE DATA

incidence_plot <- 
  dat_state %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  mutate(cases_new.7avg_week = mean(state_new_cases_7avg, na.rm  = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = week, y = cases_new.7avg_week))+
  geom_bar(position = "dodge", stat = "identity", fill = "darkgoldenrod2")+
  geom_vline(xintercept = as.Date("2022-12-21"), lty = "longdash")+
  geom_vline(xintercept = as.Date("2023-09-21"), lty = "longdash")+
  geom_vline(xintercept = as.Date("2024-01-05"), lty = "longdash")+
  theme_bw() +
  labs(title = "New cases",
       x = "",
       y = "Cases\n(weekly total)")+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month", limits = c(start_date, end_date))+
  theme(axis.text.x = element_text(angle = 90))
incidence_plot

hosp_plot <-
  dat_hosp %>%
  group_by(date) %>%
  # daily statewide admissions (7 day average)
  mutate(hos.7avg_week = sum(hos.7avg, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  # weekly total admissions
  mutate(hos.7avg_week = mean(hos.7avg_week, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(week>= min(dat_state$date)) %>%
  ggplot(aes(x = week, y = hos.7avg_week)) +
  geom_bar(position = "dodge", stat = "identity", fill = "tomato")+
  geom_vline(xintercept = as.Date("2022-12-21"), lty = "longdash")+
  geom_vline(xintercept = as.Date("2023-09-21"), lty = "longdash")+
  geom_vline(xintercept = as.Date("2024-01-05"), lty = "longdash")+
  theme_bw() +
  labs(title = "New hospital admissions",
       x = "",
       y = "Hospital admissions\n(weekly total)")+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month", limits = c(start_date, end_date))+
  theme(axis.text.x = element_text(angle = 90))
hosp_plot

ww_plot <- 
  dat_state %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  mutate(sars2.7avg_week = mean(state_sars2, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = week, y = sars2.7avg_week))+
  geom_line(color = "darkblue", linewidth = 1)+
  geom_vline(xintercept = as.Date("2022-12-21"), lty = "longdash")+
  geom_vline(xintercept = as.Date("2023-09-21"), lty = "longdash")+
  geom_vline(xintercept = as.Date("2024-01-05"), lty = "longdash")+
  theme_bw() +
  labs(title = "SARS-CoV-2 wastewater concentration",
       x = "",
       y = "Concentration\n(copies/mL)")+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month", limits = c(start_date, end_date))+
  theme(axis.text.x = element_text(angle = 90))
ww_plot

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

# panel
#grid.arrange(incidence_plot, hosp_plot, ww_plot, rt_case_plot, nrow = 4)
figure_1 <- 
  cowplot::plot_grid(incidence_plot, hosp_plot, ww_plot, rt_case_plot, ncol = 1, align = "v",
                     labels = c("A)", "B)", "C)", "D)"), hjust = -2)



cowplot::plot_grid(incidence_plot,  rt_case_plot, ncol = 1, align = "v",
                   labels = c("A)", "B)"), hjust = -2)

# no hosp admissions
figure_1 <- 
  cowplot::plot_grid(incidence_plot, ww_plot, rt_case_plot, ncol = 1, align = "v",
                     labels = c("A)", "B)", "C)"), hjust = -2)

# export
png("/Figure 1 - main panel_no hosp.png",
    units = "in",
    width = 8.5, height =9,
    res = 600)
figure_1
dev.off()