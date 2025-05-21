##
# REVISED FIGURE 2

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

# panel 1 - state data
# # #
# plot with models at state level (statewide weighted mean)

# case rt (EpiEstim)
rt_cases_weekly <- readRDS("data/rt_cases_weekly.rds") %>%
  mutate(model = "Rt cases")

# ww substitution
rt_ww_sub <- readRDS("data/rt_ww_sub.rds")
rt_ww_sub <- rt_ww_sub %>%
  arrange(date) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(model = "EpiEstim substitution")

# EpiSewer method - statewide
episewer_state <- readRDS("data/episewer_state.rds") %>%
  arrange(week) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(model= "EpiSewer")

# exp change rate
rt_exp_weekly <- readRDS("data/rt_exp_weekly.rds") %>%
  mutate(model = "Exp change rate")

# fit line
rt_fit_weekly <- readRDS("data/rt_fit_weekly.rds") %>%
  mutate(model = "Fit line method")

# goldstein method statewide (note - need to rerun this with 200 resamples)
eirr_state <- readRDS("data/eirr_state.rds") %>%
  arrange(week) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(model = "Goldstein method")

# huisman method
rt_huisman_weekly <- readRDS("data/rt_huisman_weekly.rds") %>%
  arrange(week) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(model = "Huisman method")

# ern state
ern_state <- readRDS("data/ern_state_weekly.rds") %>%
  arrange(week) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  ) %>%
  mutate(model = "ERN method")

# rolling GLM
rt_glm_weekly <- readRDS("data/rt_glm_weekly.rds") %>%
  mutate(model = "Rolling GLM")
# 5 week average
rt_glm_weekly <- rt_glm_weekly %>%
  arrange(date) %>%
  mutate(mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
         se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
         ll_95_rt = rollmean(ll_95_rt, 5, align = "center", na.pad = TRUE),
         ul_95_rt = rollmean(ul_95_rt, 5, alidng = "center", na.pad = TRUE)
  )

# combine
# one large dataframe
rt_state_all <- bind_rows(
  rt_cases_weekly, rt_exp_weekly, rt_fit_weekly, rt_glm_weekly, eirr_state, ern_state, episewer_state, rt_huisman_weekly, rt_ww_sub
)

# factor the data
rt_state_all$method[rt_state_all$method == "Goldstein method"] <- "Goldstein - EIRR"
rt_state_all$method[rt_state_all$method == "EpiSewer method"] <- "EpiSewer"
rt_state_all$method[rt_state_all$method == "Rt cases"] <- "Case Rt"
rt_state_all$method[rt_state_all$method == "Fit line"] <- "Fit line method"

rt_state_all$model <- factor(rt_state_all$method, levels = c("Case Rt",
                                                            "Fit line method",
                                                            "Rolling GLM", 
                                                            "EpiEstim Substitution", 
                                                            "Exp change rate",
                                                            "Huisman method",
                                                            "Goldstein - EIRR",
                                                            "EpiSewer",
                                                            "ERN method"))

# color palette
pal <- MetBrewer::met.brewer(n = 9, name = "Cross", type = "discrete")

# match to methods
values = c(c("black", "firebrick", pal), 
           c("Rt cases",
             "EpiEstim Substitution", 
             "EpiSewer method",
             "Exp change rate",
             "Fit line", 
             "Goldstein method",
             "Huisman method",
             "ERN method",
             "Rolling GLM"))

# add case counts and ww levels as a shared plot to the facets*
state_add_data <-
  dat_state %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  mutate(sars2.7avg_week = mean(state_sars2, na.rm = TRUE),
         cases_new.7avg_week = mean(state_new_cases_7avg, na.rm  = TRUE)) %>%
  ungroup() %>%
  mutate(mean_rt = sars2.7avg_week) %>%
  mutate(model = "Observed data")

# adding another panel
# https://stackoverflow.com/questions/10673481/add-a-geom-layer-for-a-single-panel-in-a-faceted-plot

rt_state_all2 <- bind_rows(rt_state_all, state_add_data)

# factor the values
rt_state_all2$model <- factor(rt_state_all2$model, levels = c("Observed data",
                                                             "Case Rt",
                                                             "Fit line method",
                                                             "Rolling GLM", 
                                                             "EpiEstim Substitution", 
                                                             "Exp change rate",
                                                             "Huisman method",
                                                             "Goldstein - EIRR",
                                                             "EpiSewer",
                                                             "ERN method"))

# statewide plot 
plot_all_models <- 
  ggplot(data = rt_state_all2 %>%
           filter(week <= as.Date("2024-01-23"))
           )+ 
  layer(data = rt_state_all2,  geom = c( "bar"),
                    aes(x = week, y = cases_new.7avg_week/50,
                        fill = "orange"),
                    position = "dodge",
                    stat = "identity",
        #params = list(fill = "orange"),
        inherit.aes = FALSE
        )+
  geom_line(aes(x= week, y = mean_rt, color = model), linewidth = 1)+
  geom_ribbon( aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt, fill = model), alpha = 0.2) +
  geom_hline(yintercept = 1, lty = "dashed") +
  labs(title = "State level Rt comparison",
       x = "",
       y = "Rt")+
  scale_colour_manual(
    values =values,
    guide = "none")+
  scale_fill_manual(values = values)+
  theme_bw()+
  #guides(color = guide_legend(override.aes = list(linewidth= 2)))+http://127.0.0.1:26643/graphics/plot_zoom_png?width=990&height=856
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "2 months"
               #, 
               #limits = c(as.Date("2023-10-23"), as.Date("2024-01-23"))
               )+
  
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(.7,.05),
        legend.title = element_blank(),#http://127.0.0.1:26643/graphics/plot_zoom_png?width=990&height=856
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color = "darkgrey"),
        legend.text = element_text(size = 15)
        )+
  facet_wrap(~model, scales = "free",
             nrow = 4)

p <- 
  plot_all_models + 
  layer(data = rt_state_all2,  geom = c( "point"),
        aes(x = week, y = sars2.7avg_week,
            shape = "16"),
        position = "dodge",
        stat = "identity",
        inherit.aes = FALSE
        )+
  scale_fill_manual(values = c("orange" = "orange"),
                    labels = c("orange" = "Reported\ncases"))+
  scale_shape_manual(values = c("16" = 16),
                                labels = c("Wastewater\nconcentration"))+
  guides(nrow = 1)
p

# save
png("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/Figures and tables/Figure - all methods statewide panel.png",
    units = "in",
    width = 11, height =8.5,
    res = 600)
p
dev.off()

### 
# state level, no smooth
# # #
# plot with models at state level (statewide weighted mean)

# case rt (EpiEstim)
rt_cases_weekly <- readRDS("data/rt_cases_weekly.rds") %>%
  mutate(model = "Rt cases")

# ww substitution
rt_ww_sub <- readRDS("data/rt_ww_sub.rds")
rt_ww_sub <- rt_ww_sub %>%
  mutate(model = "EpiEstim Substitution")

# EpiSewer method - statewide
episewer_state <- readRDS("data/episewer_state.rds") %>%
  mutate(model= "EpiSewer")

# exp change rate
rt_exp_weekly <- readRDS("data/rt_exp_weekly.rds") %>%
  mutate(model = "Exp change rate")

# fit line
rt_fit_weekly <- readRDS("data/rt_fit_weekly.rds") %>%
  mutate(model = "Fit line method")

# goldstein method statewide (note - need to rerun this with 200 resamples)
eirr_state <- readRDS("data/eirr_state.rds") %>%
  mutate(model = "Goldstein method")

# huisman method
rt_huisman_weekly <- readRDS("data/rt_huisman_weekly.rds") %>%
  mutate(model = "Huisman method")

# ern state
ern_state <- readRDS("data/ern_state_weekly.rds") %>%
  mutate(model = "ERN method")

# rolling GLM
rt_glm_weekly <- readRDS("data/rt_glm_weekly.rds") %>%
  mutate(model = "Rolling GLM")

# combine and make the plot
# one large dataframe
rt_state_all <- bind_rows(
  rt_cases_weekly, rt_exp_weekly, rt_fit_weekly, rt_glm_weekly, eirr_state, ern_state, episewer_state, rt_huisman_weekly, rt_ww_sub
)
rt_state_all$method <- rt_state_all$model
# factor the data
rt_state_all$method[rt_state_all$method == "Goldstein method"] <- "Goldstein - EIRR"
rt_state_all$method[rt_state_all$method == "EpiSewer method"] <- "EpiSewer"
rt_state_all$method[rt_state_all$method == "Rt cases"] <- "Case Rt"
rt_state_all$method[rt_state_all$method == "Fit line"] <- "Fit line method"

rt_state_all$model <- factor(rt_state_all$method, levels = c("Case Rt",
                                                             "Fit line method",
                                                             "Rolling GLM", 
                                                             "EpiEstim Substitution", 
                                                             "Exp change rate",
                                                             "Huisman method",
                                                             "Goldstein - EIRR",
                                                             "EpiSewer",
                                                             "ERN method"))

# color palette
pal <- MetBrewer::met.brewer(n = 9, name = "Cross", type = "discrete")

# match to methods
values = c(c("black", "firebrick", pal), 
           c("Rt cases",
             "EpiEstim Substitution", 
             "EpiSewer method",
             "Exp change rate",
             "Fit line", 
             "Goldstein method",
             "Huisman method",
             "ERN method",
             "Rolling GLM"))

# add case counts and ww levels as a shared plot to the facets*
state_add_data <-
  dat_state %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  mutate(sars2.7avg_week = mean(state_sars2, na.rm = TRUE),
         cases_new.7avg_week = mean(state_new_cases_7avg, na.rm  = TRUE)) %>%
  ungroup() %>%
  mutate(mean_rt = sars2.7avg_week) %>%
  mutate(model = "Observed data")

# adding another panel
# https://stackoverflow.com/questions/10673481/add-a-geom-layer-for-a-single-panel-in-a-faceted-plot

rt_state_all2 <- bind_rows(rt_state_all, state_add_data)

# factor the values
rt_state_all2$model <- factor(rt_state_all2$model, levels = c("Observed data",
                                                              "Case Rt",
                                                              "Fit line method",
                                                              "Rolling GLM", 
                                                              "EpiEstim Substitution", 
                                                              "Exp change rate",
                                                              "Huisman method",
                                                              "Goldstein - EIRR",
                                                              "EpiSewer",
                                                              "ERN method"))

# statewide plot 
plot_all_models <- 
  ggplot(data = rt_state_all2 %>%
           filter(week <= as.Date("2024-01-23"))
  )+ 
  layer(data = rt_state_all2,  geom = c( "bar"),
        aes(x = week, y = cases_new.7avg_week/50,
            fill = "orange"),
        position = "dodge",
        stat = "identity",
        #params = list(fill = "orange"),
        inherit.aes = FALSE
  )+
  geom_line(aes(x= week, y = mean_rt, color = model), linewidth = 1)+
  geom_ribbon( aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt, fill = model), alpha = 0.2) +
  geom_hline(yintercept = 1, lty = "dashed") +
  labs(title = "State level Rt comparison",
       x = "",
       y = "Rt")+
  scale_colour_manual(
    values =values,
    guide = "none")+
  scale_fill_manual(values = values)+
  theme_bw()+
  #guides(color = guide_legend(override.aes = list(linewidth= 2)))+http://127.0.0.1:26643/graphics/plot_zoom_png?width=990&height=856
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "2 months"
               #, 
               #limits = c(as.Date("2023-10-23"), as.Date("2024-01-23"))
  )+
  
  theme(axis.text.x = element_text(angle = 90),
        legend.position = c(.7,.05),
        legend.title = element_blank(),#http://127.0.0.1:26643/graphics/plot_zoom_png?width=990&height=856
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color = "darkgrey"),
        legend.text = element_text(size = 15)
  )+
  facet_wrap(~model, scales = "free",
             nrow = 4)

p <- 
  plot_all_models + 
  layer(data = rt_state_all2,  geom = c( "point"),
        aes(x = week, y = sars2.7avg_week,
            shape = "16"),
        position = "dodge",
        stat = "identity",
        inherit.aes = FALSE
  )+
  scale_fill_manual(values = c("orange" = "orange"),
                    labels = c("orange" = "Reported\ncases"))+
  scale_shape_manual(values = c("16" = 16),
                     labels = c("Wastewater\nconcentration"))+
  guides(nrow = 1)
p

# save
png("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/Figures and tables/Figure - all methods statewide no smooth panel.png",
    units = "in",
    width = 11, height =8.5,
    res = 600)
p
dev.off()

# ---------------------------------------------------------------------------------------------------------------------------


### Rt 90 day window - sewershed level data


# start and end dates for time series
start_date <- as.Date("2022-09-01")
end_date <- as.Date("2024-02-20")

# color palette
pal <- MetBrewer::met.brewer(n = 9, name = "Cross", type = "discrete")

# match to methods
values = c(pal, c("Rt cases",
                  "EpiEstim substitution", 
                  "EpiSewer",
                  "Exp change rate",
                  "Fit line method", 
                  "Goldstein method",
                  "Huisman method",
                  "ERN method",
                  "Rolling GLM"))
# plot with models at sewershed level (one sewershed)

# data load
dat <- readRDS("data/rt_sewershed_data_fill.rds")
cases <- readRDS("data/cases_county_sewershed.rds")
dat <- dat %>%
  group_by(sw_id) %>%
  fill(county) %>%
  left_join(cases, by = c("county", "date")) %>%
  ungroup()

episewer_sewershed <- readRDS("data/rt_episewer_sw.rds")
eirr_sewershed <- readRDS("data/eirr_sewershed_9.rds")

ggplot()+
  geom_line(data = episewer_sewershed, aes(x = week, y = mean_rt))+
  geom_ribbon(data = episewer_sewershed, aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt), alpha = 0.5)+
  facet_wrap(~sw_id)

# compare to case Rt to pick our exemplary example
s_list <- readRDS("data/s_list.rds") # list of sewersheds to iterate over
s_list <- head(unique(s_list$sw_id), 9)

# Reference model - case based Rt
source("functions/Rt-EpiEstim-functions.R")

# loop through counties
reference_rt_list <- list()

for(i in unique(s_list)){
  
  # example county df
  cases_df <- dat%>%
    ungroup() %>%
    filter(sw_id == i)%>%
    select(date, cases_new.7avg) %>%
    rename(Date = date) 
  
  # Rt from  case data with SI 4 (not much difference in the unknown si results)
  orange_rt_cases_weekly <- rt_function(cases_df, mean_si = 4, std_si = 1, weekly = "Yes")
  
  # add county
  orange_rt_cases_weekly$sw_id <- i
  
  # store in list
  reference_rt_list[[i]] <- orange_rt_cases_weekly
}

# make list object a dataframe
rt_cases_weekly_sewer <- do.call(rbind, reference_rt_list)

ggplot()+
  geom_line(data = rt_cases_weekly_sewer, aes(x = week, y = mean_rt))+
  geom_ribbon(data = rt_cases_weekly_sewer, aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt), alpha = 0.5)+
  facet_wrap(~sw_id)

# sewershed we selected
sewer_select <- "36085NY0026107NYCWWN"

# fit line
source("functions/Rt-fit-line-function.R")

# case rt
cases_df <- dat%>%
  ungroup() %>%
  filter(sw_id == sewer_select)%>%
  select(date, cases_new.7avg) %>%
  rename(Date = date) 

# Rt from  case data with SI 4 (not much difference in the unknown si results)
rt_cases_weekly <- rt_function(cases_df, mean_si = 4, std_si = 1, weekly = "No")

#extract ww data and join it to rt data
ww_data <- dat %>%
  select(date, sw_id, sars2.7avg) %>%
  filter(sw_id == sewer_select)

# one dataframe for the county with ww data and rt from case data
data_fit_line <- inner_join(ww_data, rt_cases_weekly, c("date")) %>%
  rename(rt = mean_rt)

# run function
rt_fit_weekly <- rt_fit_line_function(dataframe = data_fit_line, weekly = "Yes", range = 45, predictor = "sars2.7avg")
rt_fit_weekly$model <- "Fit line method"

# rolling glm
source("functions/Rt-glm-roll.R")
# wastewater data and case data
ww_data <- dat%>%
  select(date, sw_id, sars2.7avg, cases_new.7avg)%>%
  group_by(sw_id)%>%
  arrange(date) %>%
  mutate(
    log_sars2.7avg = (log(sars2.7avg))
  )%>%
  ungroup() %>%
  filter(sw_id == sewer_select)

# run the function
rt_glm_weekly <- rt_glm_function(dataframe = ww_data,
                                 predictor = "log_sars2.7avg",
                                 case_data = "cases_new.7avg",
                                 range = 45,
                                 unknown_si = "Yes",
                                 mean_si = 4,
                                 std_si = 1,
                                 weekly = "Yes")
# 5 week rolling window
rt_glm_sewershed_smooth <- rt_glm_weekly %>%
  arrange(week) %>%
  mutate(
    mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
    se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
    ll_95_rt = rollmean(ll_95_rt,5, align = "center", na.pad = TRUE),
    ul_95_rt = rollmean(ul_95_rt,5, align = "center", na.pad = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model = "Rolling GLM")

# epiestim sub
source("functions/Rt-EpiEstim-substitution.R")
# wastewater data
ww_data <- dat %>%
  select(date, sw_id, sars2.7avg) %>%
  filter(sw_id == sewer_select) %>%
  ungroup() %>%
  select(-sw_id)

# run function
rt_ww_sub <- rt_ww_sub_function(dataframe = ww_data, mean_si = 4, std_si = 1,  weekly = "Yes", unknown_si = "Yes")
# smooth
rt_sub_sewershed_smooth <- rt_ww_sub %>%
  arrange(week) %>%
  mutate(
    mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
    se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
    ll_95_rt = rollmean(ll_95_rt,5, align = "center", na.pad = TRUE),
    ul_95_rt = rollmean(ul_95_rt,5, align = "center", na.pad = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model = "EpiEstim substitution")

# exp change rate
source("functions/Rt-exp-change-rate.R")
# select county level data
ww_data <- dat %>%
  select(date, sw_id, sars2.7avg) %>%
  filter(sw_id == sewer_select) %>%
  select(-sw_id)

# log sars 2
ww_data$sars2.7avg <- log(ww_data$sars2.7avg)

# test the function
rt_exp_weekly <- rt_change_rate_function(dataframe = ww_data, change_window = 1, weekly = "Yes")
rt_exp_weekly$model <- "Exp change rate"

# huisman
source("functions/Rt-huisman-functions.r")
# filter for county we want to analyze
ww_data <- dat %>%
  select(date, sw_id, sars2.7avg) %>%
  filter(sw_id == sewer_select) %>%
  select(-sw_id)

# run function
rt_huisman_weekly <- rt_huisman_function(dataframe = ww_data, region = i, weekly = "Yes")

# smooth step
# 5 week rolling window by county
rt_huisman_weekly2  <- rt_huisman_weekly %>%
  arrange(week) %>%
  mutate(
    mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
    se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
    ll_95_rt = rollmean(ll_95_rt,5, align = "center", na.pad = TRUE),
    ul_95_rt = rollmean(ul_95_rt,5, align = "center", na.pad = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model = "Huisman method")

# ern
library(ern)
source("functions/Rt-ERN-function.r")

# filter for county we want to analyze
ww_data <- dat%>%
  select(date, sw_id, sars2.7avg, cases_new.7avg) %>%
  filter(sw_id == sewer_select) %>%
  select(-sw_id)

# run function
rt_ern_weekly <- Rt_ERN_function(dataframe = ww_data, weekly = "Yes")

# smooth step
# 5 week rolling window by county
ern_rt_smooth <- rt_ern_weekly %>%
  arrange(week) %>%
  mutate(
    mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
    se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
    ll_95_rt = rollmean(ll_95_rt,5, align = "center", na.pad = TRUE),
    ul_95_rt = rollmean(ul_95_rt,5, align = "center", na.pad = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model = "ERN method")


# plotting object
episewer_sewershed_select <- episewer_sewershed %>%
  filter(sw_id == sewer_select) %>%
  mutate(model = "EpiSewer")
rt_case_select <- rt_cases_weekly_sewer %>%
  filter(sw_id == sewer_select) %>%
  mutate(model = "Rt cases")


eirr_sewer_select <- eirr_sewer[2:7]
eirr_sewer_select <- eirr_sewer_select %>%
  filter(sw_id %in% sewer_select) %>%
  mutate(model = "Goldstein - EIRR")
rt_huisman_weekly2 <- rt_huisman_weekly2 %>%
  filter(week >= as.Date("2023-11-01") & week <= as.Date("2024-02-01")) 
  

sewer_results <- bind_rows(episewer_sewershed_select, rt_case_select, rt_fit_weekly, rt_glm_sewershed_smooth, rt_sub_sewershed_smooth,
                           rt_exp_weekly, rt_huisman_weekly2 , ern_rt_smooth,
                           eirr_sewer_select) 
sewer_results$model <- factor(sewer_results$model, levels = c("Rt cases",
                                                              "EpiEstim substitution", 
                                                              "EpiSewer",
                                                              "Exp change rate",
                                                              "Fit line method", 
                                                              "Goldstein - EIRR",
                                                              "Huisman method",
                                                              "ERN method",
                                                              "Rolling GLM"))

# add case and ww data
# SEWERSHED DATA
sewer_select <- "36085NY0026107NYCWWN"
sewer_raw_data <- readRDS("data/rt_sewershed_data_fill.rds") %>%
  filter(sw_id == sewer_select) %>%
  filter(date >= as.Date("2023-11-01") & date <= as.Date("2024-02-01")) %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  summarize(
    mean_sars2 = mean(sars2.7avg, na.rm = TRUE)
  ) %>%
  mutate(mean_rt = mean_sars2) %>%
  mutate(model = "Observed data")
sewer_cases <- readRDS("data/cases_county_sewershed.rds") %>%
  filter(county == "Richmond") %>%
  filter(date >= as.Date("2023-11-01") & date <= as.Date("2024-02-01")) %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  summarize(
    cases = sum(cases_new.7avg, na.rm = TRUE)
  )
sewer_raw_data <- left_join(sewer_raw_data, sewer_cases, by = c("week"))

# add to sewer data
sewer_results_2 <- bind_rows(sewer_results, sewer_raw_data) %>%
  filter(week >= as.Date("2023-11-01") & week <= as.Date("2024-02-01"))

# add to factor
sewer_results_2$model[sewer_results_2$model == "Rt cases"] <- "Case Rt"
sewer_results_2$model[sewer_results_2$model == "EpiEstim substitution"] <- "EpiEstim Substitution"

sewer_results_2$model <- factor(sewer_results_2$model, levels = c("Observed data",
                                                              "Case Rt",
                                                              "Fit line method",
                                                              "Rolling GLM", 
                                                              "EpiEstim Substitution", 
                                                              "Exp change rate",
                                                              "Huisman method",
                                                              "Goldstein - EIRR",
                                                              "EpiSewer",
                                                              "ERN method"))

# summarize by week and bin to the dates we are using
# match to methods
values = c(c("black", "firebrick", pal), 
           c("Rt cases",
             "EpiEstim Substitution", 
             "EpiSewer method",
             "Exp change rate",
             "Fit line", 
             "Goldstein method",
             "Huisman method",
             "ERN method",
             "Rolling GLM"))

rt_sewer_plot <- 
  ggplot(data = sewer_results_2 %>%
           filter(week <= "2024-01-27")
         )+
  layer(data = sewer_results_2%>%
          filter(week <= "2024-01-27"),  geom = c( "bar"),
        aes(x = week, y = cases/5,
            fill = "orange"),
        position = "dodge",
        stat = "identity",
        #params = list(fill = "orange"),
        inherit.aes = FALSE
  ) +
  geom_line(aes(x = week, y = mean_rt, color = model), linewidth = 1)+
  geom_ribbon(aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt, fill = model), alpha = 0.2)+
  scale_color_manual(values = values,
                     guide = "none")+
  scale_fill_manual(values = values,
                    guide = "none")+
  theme_bw()+
  #guides(color = guide_legend(override.aes = list(linewidth= 2)))+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month", limits = c(as.Date("2023-11-01"), as.Date("2024-02-01")))+
  theme(#axis.text.x = element_text(angle = 90),
    legend.position = c(.7,.1),
    legend.title = element_blank(),
    #plot.margin=unit(c(0.5,0,-0.5,0), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "darkgrey"),
    legend.text = element_text(size = 15)
    
    )+
  #ylim(0.2,3)+
  labs(x = "",
       y = expression("R"[t]),
       title = "Sewershed level Rt comparison (Port Richmond WWTP, Richmond County, NY)"
       )+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_wrap(~model, scales = "free")

p <- 
  rt_sewer_plot+ 
  layer(data =  sewer_results_2%>%
          filter(week <= "2024-01-27"),  geom = c( "point"),
        aes(x = week, y = mean_sars2,
            shape = "16"),
        position = "dodge",
        stat = "identity",
        inherit.aes = FALSE
  )+
  scale_fill_manual(values = c("orange" = "orange"),
                    labels = c("orange" = "Reported\ncases"))+
  scale_shape_manual(values = c("16" = 16),
                     labels = c("Wastewater\nconcentration"))+
  guides(nrow = 1)

p

# save
png("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/Figures and tables/Figure - all methods sewershed panel.png",
    units = "in",
    width = 11.5, height =8.5,
    res = 600)
p
dev.off()

# ### 90 day sewershed comparison - no smoothing


# start and end dates for time series
start_date <- as.Date("2022-09-01")
end_date <- as.Date("2024-02-20")

# color palette
pal <- MetBrewer::met.brewer(n = 9, name = "Cross", type = "discrete")

# match to methods
values = c(pal, c("Rt cases",
                  "EpiEstim substitution", 
                  "EpiSewer",
                  "Exp change rate",
                  "Fit line method", 
                  "Goldstein method",
                  "Huisman method",
                  "ERN method",
                  "Rolling GLM"))
# plot with models at sewershed level (one sewershed)

# data load
dat <- readRDS("data/rt_sewershed_data_fill.rds")
cases <- readRDS("data/cases_county_sewershed.rds")
dat <- dat %>%
  group_by(sw_id) %>%
  fill(county) %>%
  left_join(cases, by = c("county", "date")) %>%
  ungroup()

episewer_sewershed <- readRDS("data/rt_episewer_sw.rds")
eirr_sewershed <- readRDS("data/eirr_sewershed_9.rds")

ggplot()+
  geom_line(data = episewer_sewershed, aes(x = week, y = mean_rt))+
  geom_ribbon(data = episewer_sewershed, aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt), alpha = 0.5)+
  facet_wrap(~sw_id)

# compare to case Rt to pick our exemplary example
s_list <- readRDS("data/s_list.rds") # list of sewersheds to iterate over
s_list <- head(unique(s_list$sw_id), 9)

# Reference model - case based Rt
source("functions/Rt-EpiEstim-functions.R")

# loop through counties
reference_rt_list <- list()

for(i in unique(s_list)){
  
  # example county df
  cases_df <- dat%>%
    ungroup() %>%
    filter(sw_id == i)%>%
    select(date, cases_new.7avg) %>%
    rename(Date = date) 
  
  # Rt from  case data with SI 4 (not much difference in the unknown si results)
  orange_rt_cases_weekly <- rt_function(cases_df, mean_si = 4, std_si = 1, weekly = "Yes")
  
  # add county
  orange_rt_cases_weekly$sw_id <- i
  
  # store in list
  reference_rt_list[[i]] <- orange_rt_cases_weekly
}

# make list object a dataframe
rt_cases_weekly_sewer <- do.call(rbind, reference_rt_list)

ggplot()+
  geom_line(data = rt_cases_weekly_sewer, aes(x = week, y = mean_rt))+
  geom_ribbon(data = rt_cases_weekly_sewer, aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt), alpha = 0.5)+
  facet_wrap(~sw_id)

# sewershed we selected
sewer_select <- "36085NY0026107NYCWWN"

# fit line
source("functions/Rt-fit-line-function.R")

# case rt
cases_df <- dat%>%
  ungroup() %>%
  filter(sw_id == sewer_select)%>%
  select(date, cases_new.7avg) %>%
  rename(Date = date) 

# Rt from  case data with SI 4 (not much difference in the unknown si results)
rt_cases_weekly <- rt_function(cases_df, mean_si = 4, std_si = 1, weekly = "No")

#extract ww data and join it to rt data
ww_data <- dat %>%
  select(date, sw_id, sars2.7avg) %>%
  filter(sw_id == sewer_select)

# one dataframe for the county with ww data and rt from case data
data_fit_line <- inner_join(ww_data, rt_cases_weekly, c("date")) %>%
  rename(rt = mean_rt)

# run function
rt_fit_weekly <- rt_fit_line_function(dataframe = data_fit_line, weekly = "Yes", range = 45, predictor = "sars2.7avg")
rt_fit_weekly$model <- "Fit line method"

# rolling glm
source("functions/Rt-glm-roll.R")
# wastewater data and case data
ww_data <- dat%>%
  select(date, sw_id, sars2.7avg, cases_new.7avg)%>%
  group_by(sw_id)%>%
  arrange(date) %>%
  mutate(
    log_sars2.7avg = (log(sars2.7avg))
  )%>%
  ungroup() %>%
  filter(sw_id == sewer_select)

# run the function
rt_glm_weekly <- rt_glm_function(dataframe = ww_data,
                                 predictor = "log_sars2.7avg",
                                 case_data = "cases_new.7avg",
                                 range = 45,
                                 unknown_si = "Yes",
                                 mean_si = 4,
                                 std_si = 1,
                                 weekly = "Yes")
# 5 week rolling window
rt_glm_sewershed_smooth <- rt_glm_weekly %>%
  mutate(model = "Rolling GLM")

# epiestim sub
source("functions/Rt-EpiEstim-substitution.R")
# wastewater data
ww_data <- dat %>%
  select(date, sw_id, sars2.7avg) %>%
  filter(sw_id == sewer_select) %>%
  ungroup() %>%
  select(-sw_id)

# run function
rt_ww_sub <- rt_ww_sub_function(dataframe = ww_data, mean_si = 4, std_si = 1,  weekly = "Yes", unknown_si = "Yes")
# smooth
rt_sub_sewershed_smooth <- rt_ww_sub %>%
  mutate(model = "EpiEstim substitution")

# exp change rate
source("functions/Rt-exp-change-rate.R")
# select county level data
ww_data <- dat %>%
  select(date, sw_id, sars2.7avg) %>%
  filter(sw_id == sewer_select) %>%
  select(-sw_id)

# log sars 2
ww_data$sars2.7avg <- log(ww_data$sars2.7avg)

# test the function
rt_exp_weekly <- rt_change_rate_function(dataframe = ww_data, change_window = 1, weekly = "Yes")
rt_exp_weekly$model <- "Exp change rate"

# huisman
source("functions/Rt-huisman-functions.r")
# filter for county we want to analyze
ww_data <- dat %>%
  select(date, sw_id, sars2.7avg) %>%
  filter(sw_id == sewer_select) %>%
  select(-sw_id)

# run function
rt_huisman_weekly <- rt_huisman_function(dataframe = ww_data, region = i, weekly = "Yes")

# smooth step
# 5 week rolling window by county
rt_huisman_weekly2  <- rt_huisman_weekly %>%
  mutate(model = "Huisman method")

# ern
library(ern)
source("functions/Rt-ERN-function.r")

# filter for county we want to analyze
ww_data <- dat%>%
  select(date, sw_id, sars2.7avg, cases_new.7avg) %>%
  filter(sw_id == sewer_select) %>%
  select(-sw_id)

# run function
rt_ern_weekly <- Rt_ERN_function(dataframe = ww_data, weekly = "Yes")

# smooth step
# 5 week rolling window by county
ern_rt_smooth <- rt_ern_weekly %>%
  mutate(model = "ERN method")


# plotting object
episewer_sewershed_select <- episewer_sewershed %>%
  filter(sw_id == sewer_select) %>%
  mutate(model = "EpiSewer")
rt_case_select <- rt_cases_weekly_sewer %>%
  filter(sw_id == sewer_select) %>%
  mutate(model = "Rt cases")


eirr_sewer_select <- eirr_sewer[2:7]
eirr_sewer_select <- eirr_sewer_select %>%
  filter(sw_id %in% sewer_select) %>%
  mutate(model = "Goldstein - EIRR")
rt_huisman_weekly2 <- rt_huisman_weekly2 %>%
  filter(week >= as.Date("2023-11-01") & week <= as.Date("2024-02-01")) 


sewer_results <- bind_rows(episewer_sewershed_select, rt_case_select, rt_fit_weekly, rt_glm_sewershed_smooth, rt_sub_sewershed_smooth,
                           rt_exp_weekly, rt_huisman_weekly2 , ern_rt_smooth,
                           eirr_sewer_select) 
sewer_results$model <- factor(sewer_results$model, levels = c("Rt cases",
                                                              "EpiEstim substitution", 
                                                              "EpiSewer",
                                                              "Exp change rate",
                                                              "Fit line method", 
                                                              "Goldstein - EIRR",
                                                              "Huisman method",
                                                              "ERN method",
                                                              "Rolling GLM"))

# add case and ww data
# SEWERSHED DATA
sewer_select <- "36085NY0026107NYCWWN"
sewer_raw_data <- readRDS("data/rt_sewershed_data_fill.rds") %>%
  filter(sw_id == sewer_select) %>%
  filter(date >= as.Date("2023-11-01") & date <= as.Date("2024-02-01")) %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  summarize(
    mean_sars2 = mean(sars2.7avg, na.rm = TRUE)
  ) %>%
  mutate(mean_rt = mean_sars2) %>%
  mutate(model = "Observed data")
sewer_cases <- readRDS("data/cases_county_sewershed.rds") %>%
  filter(county == "Richmond") %>%
  filter(date >= as.Date("2023-11-01") & date <= as.Date("2024-02-01")) %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  summarize(
    cases = sum(cases_new.7avg, na.rm = TRUE)
  )
sewer_raw_data <- left_join(sewer_raw_data, sewer_cases, by = c("week"))

# add to sewer data
sewer_results_2 <- bind_rows(sewer_results, sewer_raw_data) %>%
  filter(week >= as.Date("2023-11-01") & week <= as.Date("2024-02-01"))

# add to factor
sewer_results_2$model[sewer_results_2$model == "Rt cases"] <- "Case Rt"
sewer_results_2$model[sewer_results_2$model == "EpiEstim substitution"] <- "EpiEstim Substitution"

sewer_results_2$model <- factor(sewer_results_2$model, levels = c("Observed data",
                                                                  "Case Rt",
                                                                  "Fit line method",
                                                                  "Rolling GLM", 
                                                                  "EpiEstim Substitution", 
                                                                  "Exp change rate",
                                                                  "Huisman method",
                                                                  "Goldstein - EIRR",
                                                                  "EpiSewer",
                                                                  "ERN method"))

# summarize by week and bin to the dates we are using
# match to methods
values = c(c("black", "firebrick", pal), 
           c("Rt cases",
             "EpiEstim Substitution", 
             "EpiSewer method",
             "Exp change rate",
             "Fit line", 
             "Goldstein method",
             "Huisman method",
             "ERN method",
             "Rolling GLM"))

rt_sewer_plot <- 
  ggplot(data = sewer_results_2 %>%
           filter(week <= "2024-01-27")
  )+
  layer(data = sewer_results_2%>%
          filter(week <= "2024-01-27"),  geom = c( "bar"),
        aes(x = week, y = cases/5,
            fill = "orange"),
        position = "dodge",
        stat = "identity",
        #params = list(fill = "orange"),
        inherit.aes = FALSE
  ) +
  geom_line(aes(x = week, y = mean_rt, color = model), linewidth = 1)+
  geom_ribbon(aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt, fill = model), alpha = 0.2)+
  scale_color_manual(values = values,
                     guide = "none")+
  scale_fill_manual(values = values,
                    guide = "none")+
  theme_bw()+
  #guides(color = guide_legend(override.aes = list(linewidth= 2)))+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month", limits = c(as.Date("2023-11-01"), as.Date("2024-02-01")))+
  theme(#axis.text.x = element_text(angle = 90),
    legend.position = c(.7,.1),
    legend.title = element_blank(),
    #plot.margin=unit(c(0.5,0,-0.5,0), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "darkgrey"),
    legend.text = element_text(size = 15)
    
  )+
  #ylim(0.2,3)+
  labs(x = "",
       y = expression("R"[t]),
       title = "Sewershed level Rt comparison without additional smoothing\n(Port Richmond WWTP, Richmond County, NY)"
  )+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_wrap(~model, scales = "free")

p <- 
  rt_sewer_plot+ 
  layer(data =  sewer_results_2%>%
          filter(week <= "2024-01-27"),  geom = c( "point"),
        aes(x = week, y = mean_sars2,
            shape = "16"),
        position = "dodge",
        stat = "identity",
        inherit.aes = FALSE
  )+
  scale_fill_manual(values = c("orange" = "orange"),
                    labels = c("orange" = "Reported\ncases"))+
  scale_shape_manual(values = c("16" = 16),
                     labels = c("Wastewater\nconcentration"))+
  guides(nrow = 1)

p

# save
png("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/Figures and tables/Figure - all methods sewershed no smooth panel.png",
    units = "in",
    width = 11.5, height =8.5,
    res = 600)
p
dev.off()

# -------------------------------------------------------------------------------------

# county Rt plot - all dates
# # # 

# plot with models at county level (one county)
huisman_rt_smooth <- readRDS("data/huisman_rt_smooth.rds") %>%
  mutate(model = "Huisman method")
rt_exp_county <- readRDS("data/rt_exp_county.rds") %>%
  mutate(model = "Exp change rate")
rt_sub_county_smooth <- readRDS("data/rt_sub_county_smooth.rds") %>%
  mutate(model = "EpiEstim substitution")
rt_glm_county_smooth <- readRDS("data/rt_glm_county_smooth.rds") %>%
  mutate(model = "Rolling GLM")
rt_fit_county <- readRDS("data/rt_fit_county.rds") %>%
  mutate(model = "Fit line method")
rt_cases_weekly_county <- readRDS("data/rt_cases_weekly_county.rds") %>%
  mutate(model = "Rt cases")
episewer_county <- readRDS("data/episewer_rt_counties.rds") %>%
  mutate(model = "EpiSewer") %>%
  mutate(
    mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
    se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
    ll_95_rt = rollmean(ll_95_rt,5, align = "center", na.pad = TRUE),
    ul_95_rt = rollmean(ul_95_rt,5, align = "center", na.pad = TRUE)
  )
goldstein_rt_counties <- readRDS("data/goldstein_rt_counties.R") %>%
  group_by(county) %>%
  arrange(week) %>%
  mutate(
    mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
    se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
    ll_95_rt = rollmean(ll_95_rt,5, align = "center", na.pad = TRUE),
    ul_95_rt = rollmean(ul_95_rt,5, align = "center", na.pad = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model = "Goldstein method")

ern_county <- readRDS("data/ern_county_raw_data.rds") %>%
  group_by(county) %>%
  arrange(week) %>%
  mutate(
    mean_rt = rollmean(mean_rt, 5, align = "center", na.pad = TRUE),
    se_rt = rollmean(se_rt, 5, align = "center", na.pad = TRUE),
    ll_95_rt = rollmean(ll_95_rt,5, align = "center", na.pad = TRUE),
    ul_95_rt = rollmean(ul_95_rt,5, align = "center", na.pad = TRUE)
  ) %>%
  ungroup() %>%
  mutate(model = "ERN method")

rt_county <- bind_rows(huisman_rt_smooth,  rt_glm_county_smooth, rt_sub_county_smooth, rt_fit_county, rt_cases_weekly_county, 
                       episewer_county, goldstein_rt_counties, ern_county,
                       rt_exp_county) %>%
  filter(county == "Richmond") %>% filter(ll_95_rt > 0.525) %>% filter(ul_95_rt < 2.2)



# county incidence and ww conc
county_case <- readRDS("data/cases_county.rds") %>%
  filter(county == "Richmond") %>%
  filter(date >= start_date & date <= end_date) %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  summarize(
    cases = sum(cases_new.7avg, na.rm = TRUE)
  )

county_ww <- readRDS("data/ww.county.rds") %>%
  filter(county == "Richmond") %>%
  filter(date >= start_date & date <= end_date) %>%
  
  group_by(week = floor_date(date, unit = "weeks")) %>%
  summarize(
    mean_sars2 = mean(sars2.7avg, na.rm = TRUE)
  ) %>%
  mutate(mean_rt = mean_sars2) 



d <- left_join(county_case, county_ww, by = c("week")) %>%
  mutate(model = "Observed data")

rt_county2 <- bind_rows(rt_county, d)

rt_county2$model[rt_county2$model == "EpiEstim substitution"] <- "EpiEstim Substitution"
rt_county2$model[rt_county2$model == "Goldstein method"] <- "Goldstein - EIRR"

rt_county2$model <- factor(rt_county2$model, levels = c("Observed data",
                                                        "Rt cases",
                                                      "EpiEstim Substitution", 
                                                      "EpiSewer",
                                                      "Exp change rate",
                                                      "Fit line method", 
                                                      "Goldstein - EIRR",
                                                      "Huisman method",
                                                      "ERN method",
                                                      "Rolling GLM"))
rt_county_plot <- 
  ggplot(data = rt_county2 )+
  layer(data = rt_county2,
        geom = c( "bar"),
        aes(x = week, y = cases/20,
            fill = "orange"),
        position = "dodge",
        stat = "identity",
        #params = list(fill = "orange"),
        inherit.aes = FALSE
  )+
  geom_line(aes(x = week, y = mean_rt, color = model), linewidth = 1)+
  geom_ribbon(aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt, fill = model), alpha = 0.2)+
  scale_color_manual(values = values,
                     guide = "none")+
  scale_fill_manual(values = values,
                    guide = "none")+
  theme_bw()+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "3 months", limits = c(start_date, end_date))+
  theme(#axis.text.x = element_text(angle = 90),
    legend.position = c(.7,.05),
    legend.title = element_blank(),#http://127.0.0.1:26643/graphics/plot_zoom_png?width=990&height=856
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "darkgrey"),
    legend.text = element_text(size = 15)
    )+
  labs(x = "",
       y = expression("R"[t]),
       title = "County level Rt comparison (Richmond County/Staten Island)")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_wrap(~model, scales = "free")

p <- 
  rt_county_plot +
  layer(data =  rt_county2,
        geom = c( "point"),
        aes(x = week, y = mean_sars2,
            shape = "16"),
        position = "dodge",
        stat = "identity",
        inherit.aes = FALSE
  )+
  scale_fill_manual(values = c("orange" = "orange"),
                    labels = c("orange" = "Reported\ncases"))+
  scale_shape_manual(values = c("16" = 16),
                     labels = c("Wastewater\nconcentration"))+
  guides(nrow = 1)
p

# save
png("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/Figures and tables/Figure - all methods county panel.png",
    units = "in",
    width = 13, height =8.5,
    res = 600)
p
dev.off()


# ## Rt county - 90 day version
rt_county_plot <- 
  ggplot(data = rt_county2 %>%
           filter(week >= as.Date("2023-10-23") & week <= as.Date("2024-01-23"))
  )+
  layer(data = rt_county2%>%
          filter(week >= as.Date("2023-10-23") & week <= as.Date("2024-01-23")),
        geom = c( "bar"),
        aes(x = week, y = cases/20,
            fill = "orange"),
        position = "dodge",
        stat = "identity",
        #params = list(fill = "orange"),
        inherit.aes = FALSE
  )+
  geom_line(aes(x = week, y = mean_rt, color = model), linewidth = 1)+
  geom_ribbon(aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt, fill = model), alpha = 0.2)+
  scale_color_manual(values = values,
                     guide = "none")+
  scale_fill_manual(values = values,
                    guide = "none")+
  theme_bw()+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "1 month",limits = c(as.Date("2023-10-23"), as.Date("2024-01-23"))
               )+
  theme(#axis.text.x = element_text(angle = 90),
    legend.position = c(.7,.05),
    legend.title = element_blank(),#http://127.0.0.1:26643/graphics/plot_zoom_png?width=990&height=856
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "darkgrey"),
    legend.text = element_text(size = 15)
  )+
  labs(x = "",
       y = expression("R"[t]),
       title = "County level Rt comparison (Richmond County/Staten Island)")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_wrap(~model, scales = "free")

p <- 
  rt_county_plot +
  layer(data =  rt_county2%>%
          filter(week >= as.Date("2023-10-23") & week <= as.Date("2024-01-23")),
        geom = c( "point"),
        aes(x = week, y = mean_sars2,
            shape = "16"),
        position = "dodge",
        stat = "identity",
        inherit.aes = FALSE
  )+
  scale_fill_manual(values = c("orange" = "orange"),
                    labels = c("orange" = "Reported\ncases"))+
  scale_shape_manual(values = c("16" = 16),
                     labels = c("Wastewater\nconcentration"))+
  guides(nrow = 1)
p

# save
# save
png("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/Figures and tables/Figure - all methods county 90 days panel.png",
    units = "in",
    width = 13, height =8.5,
    res = 600)
p
dev.off()

#
# county data - unadjusted, no smoothing of the rt values
#

# plot with models at county level (one county)
huisman_rt_smooth <- readRDS("data/huisman_rt_county.rds") %>%
  mutate(model = "Huisman method")
rt_exp_county <- readRDS("data/rt_exp_county.rds") %>%
  mutate(model = "Exp change rate")
rt_sub_county_smooth <- readRDS("data/rt_sub_county.rds") %>%
  mutate(model = "EpiEstim substitution")
rt_glm_county_smooth <- readRDS("data/rt_glm_county.rds") %>%
  mutate(model = "Rolling GLM")
rt_fit_county <- readRDS("data/rt_fit_county.rds") %>%
  mutate(model = "Fit line method")
rt_cases_weekly_county <- readRDS("data/rt_cases_weekly_county.rds") %>%
  mutate(model = "Rt cases")
episewer_county <- readRDS("data/episewer_rt_counties.rds") %>%
  mutate(model = "EpiSewer") 

goldstein_rt_counties <- readRDS("data/goldstein_rt_counties.R") %>%
  mutate(model = "Goldstein method")

ern_county <- readRDS("data/ern_county_raw_data.rds") %>%
  mutate(model = "ERN method")

rt_county <- bind_rows(huisman_rt_smooth,  rt_glm_county_smooth, rt_sub_county_smooth, rt_fit_county, rt_cases_weekly_county, 
                       episewer_county, goldstein_rt_counties, ern_county,
                       rt_exp_county) %>%
  filter(county == "Richmond") %>% filter(ll_95_rt > 0.525) %>% filter(ul_95_rt < 2.2)



# county incidence and ww conc
county_case <- readRDS("data/cases_county.rds") %>%
  filter(county == "Richmond") %>%
  filter(date >= start_date & date <= end_date) %>%
  group_by(week = floor_date(date, unit = "weeks")) %>%
  summarize(
    cases = sum(cases_new.7avg, na.rm = TRUE)
  )

county_ww <- readRDS("data/ww.county.rds") %>%
  filter(county == "Richmond") %>%
  filter(date >= start_date & date <= end_date) %>%
  
  group_by(week = floor_date(date, unit = "weeks")) %>%
  summarize(
    mean_sars2 = mean(sars2.7avg, na.rm = TRUE)
  ) %>%
  mutate(mean_rt = mean_sars2) 



d <- left_join(county_case, county_ww, by = c("week")) %>%
  mutate(model = "Observed data")

rt_county2 <- bind_rows(rt_county, d)

rt_county2$model[rt_county2$model == "EpiEstim substitution"] <- "EpiEstim Substitution"
rt_county2$model[rt_county2$model == "Goldstein method"] <- "Goldstein - EIRR"

rt_county2$model <- factor(rt_county2$model, levels = c("Observed data",
                                                        "Rt cases",
                                                        "EpiEstim Substitution", 
                                                        "EpiSewer",
                                                        "Exp change rate",
                                                        "Fit line method", 
                                                        "Goldstein - EIRR",
                                                        "Huisman method",
                                                        "ERN method",
                                                        "Rolling GLM"))
rt_county_plot <- 
  ggplot(data = rt_county2 )+
  layer(data = rt_county2,
        geom = c( "bar"),
        aes(x = week, y = cases/20,
            fill = "orange"),
        position = "dodge",
        stat = "identity",
        #params = list(fill = "orange"),
        inherit.aes = FALSE
  )+
  geom_line(aes(x = week, y = mean_rt, color = model), linewidth = 1)+
  geom_ribbon(aes(x = week, ymin = ll_95_rt, ymax = ul_95_rt, fill = model), alpha = 0.2)+
  scale_color_manual(values = values,
                     guide = "none")+
  scale_fill_manual(values = values,
                    guide = "none")+
  theme_bw()+
  scale_x_date(labels = date_format("%b %y"),
               date_breaks = "3 months", limits = c(start_date, end_date))+
  theme(#axis.text.x = element_text(angle = 90),
    legend.position = c(.7,.05),
    legend.title = element_blank(),#http://127.0.0.1:26643/graphics/plot_zoom_png?width=990&height=856
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "darkgrey"),
    legend.text = element_text(size = 15)
  )+
  labs(x = "",
       y = expression("R"[t]),
       title = "County level Rt comparison without smoothing select methods\n(Richmond County/Staten Island)")+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_wrap(~model, scales = "free")

p <- 
  rt_county_plot +
  layer(data =  rt_county2,
        geom = c( "point"),
        aes(x = week, y = mean_sars2,
            shape = "16"),
        position = "dodge",
        stat = "identity",
        inherit.aes = FALSE
  )+
  scale_fill_manual(values = c("orange" = "orange"),
                    labels = c("orange" = "Reported\ncases"))+
  scale_shape_manual(values = c("16" = 16),
                     labels = c("Wastewater\nconcentration"))+
  guides(nrow = 1)
p

# save
png("E:/Dropbox/CEMI/Wastewater/Papers/Reproductive number/Figures and tables/Figure - all methods county no smooth panel.png",
    units = "in",
    width = 13, height =8.5,
    res = 600)
p
dev.off()