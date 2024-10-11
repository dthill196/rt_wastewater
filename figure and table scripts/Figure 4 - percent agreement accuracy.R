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


# load data
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

county_covariates <- readRDS("data/county_covariates.rds")

metrics_long <- metrics_table %>%
  pivot_longer(cols = c(`Rolling GLM`, `Fit line method`, `EpiEstim Substitution`, `Exp change rate`, `Huisman method`, `Goldstein - EIRR`,
                        `ERN method`, `EpiSewer`),
               names_to = c("model"))

# remove infinities and genesee / orleans and greene counties
metrics_long <- metrics_long %>%
  filter(county != "Genesee" & county != "Orleans" & county != "Greene") %>%
  filter(value < 100)

# join to metrics data
metrics_long <- left_join(metrics_long, county_covariates, by = c("county"))

metrics_long$model <- factor(metrics_long$model, c("Rt cases",
                                                   "Fit line method",
                                                   "Rolling GLM",
                                                   "EpiEstim Substitution", 
                                                   "Exp change rate",
                                                   "Huisman method",
                                                   
                                                   "Goldstein - EIRR",
                                                   "EpiSewer",
                                                   
                                                   "ERN method"))

# above and below 1 agreement scatterplots for the main text
accuracy_scatter2 <-
  ggplot(data = metrics_long %>% filter(metric_names == "Above or below 1 percent agreement"), aes(x = log(sewer_pop), y = value))+
  geom_point(alpha = 0.5, size = 2)+
  facet_wrap(~ model,
             labeller = label_wrap_gen(20),
             nrow = 2)+
  geom_smooth(method = "lm", color = "darkblue", linewidth = 1)+
  stat_cor(aes(label = ..r.label..))+
  theme_bw()+
  labs(title = expression("R"["t"]~"evaluation metric for above or below 1 percent agreement and correlation with sewershed population coverage"),
       x = "log(sewershed population)",
       y = "Above or below 1 percent agreement")
accuracy_scatter2

# save
png("Figure - percent agreement accuracy scatter.png",
    units = "in",
    width = 10, height =8,
    res = 300)
accuracy_scatter2
dev.off()