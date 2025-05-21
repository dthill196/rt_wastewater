# 
# Rt paper - z score analysis for evaluation
#

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

# load data from Yiquan
metrics_table <- readRDS("data/metrics_table.rds")

# waiting for updated data with 200 resmaples from orangegrid computing run
goldstein_table <- readRDS("data/goldstein_eval_table.rds")

# ern data
ern_table <- readRDS("data/ern_county_metrics.rds")
ern_table_no_omicron <- readRDS("data/ern_county_metrics_no_omicron.rds")

# episewer data
episewer_eval <- readRDS("data/episewer_eval_table.rds")

# because we include omicron for all of the other methods, we can either 1) keep omicron in and compare or 2) leave out omicron
# for all methods, or 3) include both ERN summaries in the comparison to see if that changes the final outcome of which method is better
# in the density plot, overall, inclusion of omicron dara slightly helps out ERN visually
# will let the state comparison help see if that matters most

# merge to main table
metrics_table <- left_join(metrics_table, goldstein_table, by = c("metric_names", "county"))
metrics_table <- left_join(metrics_table, ern_table, by = c("metric_names", "county"))
metrics_table <- left_join(metrics_table, episewer_eval, by = c("metric_names", "county"))

# make columns numeric
metrics_table[4:10] <- lapply(metrics_table[4:10], as.numeric)
summary(metrics_table)
metrics_table$`Rolling GLM` <- as.numeric(metrics_table$`Rolling GLM`)
metrics_table$metric_names[metrics_table$metric_names == "Above or below 1 agreement"] <- "Above or below 1 percent agreement"

# pivot dataset
metrics_long <- metrics_table %>%
  pivot_longer(cols = c(`Rolling GLM`, `Fit line method`, `EpiEstim Substitution`, `Exp change rate`, `Huisman method`, `Goldstein - EIRR`,
                        `ERN method`, `EpiSewer`),
               names_to = c("model"))

# z score analysis

# z score = (x - mean) /sd

# calculate for each metric

z_scores <- metrics_long %>%
  group_by(metric_names) %>%
  mutate(mean_value = mean(value, na.rm = TRUE),
         sd_value = sd(value, na.rm = TRUE)
  ) %>%
  mutate(
    z_score = (value - mean_value) / sd_value
  ) %>%
  mutate(mean_z = mean(z_score, na.rm = TRUE)
  ) %>%
  ungroup()

# mean z score by model
mean_z <- z_scores %>%
  group_by(model) %>%
  summarize(
    mean_z = mean(z_score, na.rm = TRUE),
    sd_z = sd(z_score, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(dist_to_zero = abs(0 - abs(mean_z)),
         se_z = sd_z/sqrt(n)
  ) %>%
  mutate(
    ll = dist_to_zero - 1.96*se_z,
    ul = dist_to_zero + 1.96*se_z
  )

# boxplot of z scores
ggplot()+
  geom_boxplot(data = z_scores,
               aes(x = model, y = z_score, fill = model))


# factor values
mean_z$model <- factor(mean_z$model, c("Rt cases",
                                       "Fit line method",
                                       "Rolling GLM",
                                       "EpiEstim Substitution", 
                                       "Exp change rate",
                                       "Huisman method",
                                       
                                       "Goldstein - EIRR",
                                       "EpiSewer",
                                       
                                       "ERN method"))
# line plot
county_z_plot <- 
  ggplot()+
  annotate(geom = "rect",
           xmin = 0.5, xmax = 1.5,
           ymin = -0.079,
           ymax = 0.12,
           fill = NA,
           color = "black",
           lwd = 1)+
  geom_pointrange(
    data = mean_z,
    aes(x = reorder(model, dist_to_zero),
        y = dist_to_zero,
        ymin = ll,
        ymax = ul,
        color = model),
    size = 1, lwd = 1
  )+
  scale_color_manual(values = values)+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        legend.position = "none")+
  labs(x = "",
       y = "Z-score distance from zero",
       color = "",
       title = "B) County level best method")+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

county_z_plot

# z score for state data
z_scores_state <- state_eval_table %>%
  pivot_longer(
    cols = c(`Fit line method`, `Rolling GLM`, `EpiEstim Substitution`, `Exp change rate`, `Huisman method`,
             `Goldstein method`, `EpiSewer method`, `ERN method`)
  ) %>%
  group_by(metric_names) %>%
  mutate(mean_value = mean(value, na.rm = TRUE),
         sd_value = sd(value, na.rm = TRUE)
  ) %>%
  mutate(
    z_score = (value - mean_value) / sd_value
  ) %>%
  mutate(mean_z = mean(z_score, na.rm = TRUE)
  ) %>%
  ungroup()

mean_z_state <- z_scores_state %>%
  group_by(name) %>%
  summarize(
    mean_z = mean(z_score, na.rm = TRUE),
    sd_z = sd(z_score, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(dist_to_zero = abs(0 - abs(mean_z)),
         se_z = sd_z/sqrt(n)
  ) %>%
  mutate(
    ll = dist_to_zero - 1.96*se_z,
    ul = dist_to_zero + 1.96*se_z
  )

# factor values
mean_z_state$name[mean_z_state$name == "Goldstein method"] <- "Goldstein - EIRR"
mean_z_state$name[mean_z_state$name == "EpiSewer method"] <- "EpiSewer"

mean_z_state$model <- factor(mean_z_state$name, c("Rt cases",
                                                  "Fit line method",
                                                  "Rolling GLM",
                                                  "EpiEstim Substitution", 
                                                  "Exp change rate",
                                                  "Huisman method",
                                                  
                                                  "Goldstein - EIRR",
                                                  "EpiSewer",
                                                  
                                                  "ERN method"))

# line plot
state_z_plot <- 
  ggplot(data = mean_z_state)+
  annotate(geom = "rect",
           xmin = 0.5, xmax = 1.5,
           ymin = -0.5,
           ymax = 0.54,
           fill = NA,
           color = "black",
           lwd = 1)+
  geom_pointrange(
    data = mean_z_state,
    aes(x = reorder(model, dist_to_zero),
        y = dist_to_zero,
        ymin = ll,
        ymax = ul,
        color =model),
    size = 1, lwd = 1
  )+
  scale_color_manual(values = values)+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        legend.position = "none")+
  labs(x = "",
       y = "Z-score distance from zero",
       color = "",
       title = "A) State level best method")+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

state_z_plot

# z score sewershed data
metrics_sewershed <- readRDS("data/sewershed_evaluation_metrics.rds")

# z score sewershed
z_scores_sw <- metrics_sewershed %>%
  group_by(metric_names) %>%
  mutate(mean_value = mean(as.numeric(value), na.rm = TRUE),
         sd_value = sd(as.numeric(value), na.rm = TRUE)
  ) %>%
  mutate(
    z_score = (as.numeric(value) - mean_value) / sd_value
  ) %>%
  mutate(mean_z = mean(z_score, na.rm = TRUE)
  ) %>%
  ungroup()

# mean by sewershed
mean_z_sw <- z_scores_sw %>%
  group_by(model) %>%
  summarize(
    mean_z = mean(z_score, na.rm = TRUE),
    sd_z = sd(z_score, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup() %>%
  mutate(dist_to_zero = abs(0 - abs(mean_z)),
         se_z = sd_z/sqrt(n)
  ) %>%
  mutate(
    ll = dist_to_zero - 1.96*se_z,
    ul = dist_to_zero + 1.96*se_z
  )

# factor values
mean_z_sw$model <- factor(mean_z_sw$model, c("Rt cases",
                                             "Fit line method",
                                             "Rolling GLM",
                                             "EpiEstim Substitution", 
                                             "Exp change rate",
                                             "Huisman method",
                                             
                                             "Goldstein - EIRR",
                                             "EpiSewer",
                                             
                                             "ERN method"))

# line plot
sw_z_plot <- 
  ggplot(data = mean_z_sw)+
  annotate(geom = "rect",
           xmin = 0.5, xmax = 1.5,
           ymin = -0.13,
           ymax = 0.31,
           fill = NA,
           color = "black",
           lwd = 1)+
  geom_pointrange(
    data = mean_z_sw,
    aes(x = reorder(model, dist_to_zero),
        y = dist_to_zero,
        ymin = ll,
        ymax = ul,
        color = model),
    size = 1, lwd = 1
  )+
  scale_color_manual(values = values)+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        legend.position = "none")+
  labs(x = "",
       y = "Z-score distance from zero",
       color = "",
       title = "C) Sewershed level best method")+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

sw_z_plot

grid.arrange(state_z_plot, county_z_plot, sw_z_plot, nrow = 1)

# save plot
png("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/Figures and tables/Figure - z scores.png",
    units = "in",
    width = 10.5, height =3.5,
    res = 300)
grid.arrange(state_z_plot, county_z_plot, sw_z_plot, nrow = 1)
dev.off()

# jpeg
jpeg("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/paper/Figure - z scores.jpg",
     units = "in",
     width = 10.5, height =3.5,
     res = 300)
grid.arrange(state_z_plot, county_z_plot, sw_z_plot, nrow = 1)
dev.off()

# tiff
tiff("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/paper/Figure - z scores.tiff",
     units = "in",
     width = 10.5, height =3.5,
     res = 300)
grid.arrange(state_z_plot, county_z_plot, sw_z_plot, nrow = 1)
dev.off()

# one plot, side by side, ranking them by closest value to 0 for best score
# closest to 0 = closest to the average, best relative to the other scores
