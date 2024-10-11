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

# join to metrics data
metrics_long <- left_join(metrics_long, county_covariates, by = c("county"))
# TABLE - REGRESSION RESULTS

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

# try a simple regression
model_fit_line <- lm(value ~ 
                       #Prop_sewer_pop 
                       + log(sewer_pop)
                     + number_plants 
                     + mean_frequency
                     #+ frequency_factor
                     #+ frequency_factor_2
                     ,
                     data = metrics_long %>% 
                       filter(model == "Fit line method" & metric_names == "Above or below 1 percent agreement")
)
summary(model_fit_line)

model_glm <- lm(value ~ 
                  #Prop_sewer_pop 
                  + log(sewer_pop)
                + number_plants 
                + mean_frequency
                #+ frequency_factor
                #+ frequency_factor_2
                #+ log(county_pop) 
                ,
                data = metrics_long %>% 
                  filter(model == "Rolling GLM" & metric_names == "Above or below 1 percent agreement")
)
summary(model_glm)

model_episewer <- lm(value ~ 
                       #Prop_sewer_pop 
                       + log(sewer_pop)
                     + number_plants 
                     + mean_frequency
                     #+ frequency_factor
                     #+ frequency_factor_2
                     ,
                     data = metrics_long %>% 
                       filter(model == "EpiSewer" & metric_names == "Above or below 1 percent agreement")
)

model_wwsub <- lm(value ~ 
                    #Prop_sewer_pop 
                    + log(sewer_pop)
                  + number_plants 
                  + mean_frequency
                  #+ frequency_factor
                  #+ frequency_factor_2
                  ,
                  data = metrics_long %>% 
                    filter(model == "EpiEstim Substitution" & metric_names == "Above or below 1 percent agreement")
)

model_ern <- lm(value ~ 
                  #Prop_sewer_pop 
                  + log(sewer_pop)
                + number_plants 
                + mean_frequency
                #+ frequency_factor
                #+ frequency_factor_2
                ,
                data = metrics_long %>% 
                  filter(model == "ERN method" & metric_names == "Above or below 1 percent agreement")
)

model_exp <- lm(value ~ 
                  #Prop_sewer_pop 
                  + log(sewer_pop)
                + number_plants 
                + mean_frequency
                #+ frequency_factor
                #+ frequency_factor_2
                ,
                data = metrics_long %>% 
                  filter(model == "Exp change rate" & metric_names == "Above or below 1 percent agreement")
)

model_goldstein <- lm(value ~ 
                        #Prop_sewer_pop 
                        + log(sewer_pop)
                      + number_plants 
                      + mean_frequency
                      #+ frequency_factor
                      #+ frequency_factor_2
                      ,
                      data = metrics_long %>% 
                        filter(model == "Goldstein - EIRR" & metric_names == "Above or below 1 percent agreement")
)

model_huisman <- lm(value ~ 
                      #Prop_sewer_pop 
                      + log(sewer_pop)
                    + number_plants 
                    + mean_frequency
                    #+ frequency_factor
                    #+ frequency_factor_2
                    ,
                    data = metrics_long %>% 
                      filter(model == "Huisman method" & metric_names == "Above or below 1 percent agreement")
)

# coef table function
coef_table_function_lm <- function(model){
  cc <- coef(summary(model))
  cc <- within(as.data.frame(cc),
               {   `Std. Error` <- `Std. Error`
               `z value` <- Estimate/`Std. Error`
               `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
               })
  t <- printCoefmat(cc,digits=3)
  final_table <- as.data.frame(cbind(t$Estimate,t$`Std. Error`,t$`Pr(>|z|)`))
  colnames(final_table) <- c("Estimate", "SE", "p value")
  r2 <- summary(model)$r.squared
  n <- nobs(model)
  # round data to 3 digits
  A <- function(x){
    round(x, digits = 3)
  }
  final_table <- final_table %>%
    mutate(across(Estimate:`p value`, A)
    )
  final_table <- final_table %>%
    mutate(p_value = case_when(
      `p value` < 0.01 ~ paste("<0.01***"),
      `p value` >= 0.01 & `p value` <=0.05 ~ paste(`p value`, "**", sep = ""),
      `p value` > 0.05 & `p value` <= 0.01 ~ paste(`p value`, "*", sep = ""),
      `p value` > 0.01 & `p value` <= 0.2 ~ paste(`p value`, ".", sep = ""),
      `p value` > 0.2 ~ paste(`p value`)
    )) %>%
    select(-`p value`) %>%
    rename(`p value` = p_value)  
  final_table$`Estimate (P value)` <- paste(final_table$Estimate, " (", final_table$`p value`, ")", sep = "")
  final_table <- final_table %>%
    select(-Estimate, - SE) %>%
    select(`Estimate (P value)`)
  vars <- rownames(cc)
  final_table <- cbind(vars, final_table)
  
  # add r2
  r2 <- as.character(round(summary(model)$r.squared, 2))
  r_table <- as.data.frame(r2)
  rownames(r_table) <- "R2"
  r_table <- cbind(rownames(r_table), r_table)
  colnames(r_table) <- c("vars", "Estimate (P value)")
  # number of observations
  n <- as.data.frame(as.character(nobs(model)))
  rownames(n) <- "n"
  n <- cbind(rownames(n), n)
  colnames(n) <- c("vars", "Estimate (P value)")
  
  
  # combine table
  final_table <- bind_rows( final_table, r_table, n)
  final_table <- final_table %>%
    rename(t = `Estimate (P value)`)
  return(final_table)
}

table_episewer <- coef_table_function_lm(model = model_episewer) 
table_episewer <- table_episewer %>%
  rename(EpiSewer = t)
table_ern <- coef_table_function_lm(model = model_ern) %>%
  rename(`ERN method` = t)
table_exp <- coef_table_function_lm(model = model_exp) %>%
  rename(`Exp change rate` = t)
table_fit_line <- coef_table_function_lm(model = model_fit_line) %>%
  rename(`Fit line method` = t)
table_glm <- coef_table_function_lm(model = model_glm) %>%
  rename(`Rolling GLM` = t)
table_huisman <- coef_table_function_lm(model= model_huisman) %>%
  rename(`Huisman method` = t)
table_wwsub <- coef_table_function_lm(model = model_wwsub) %>%
  rename(`EpiEstim substituion` = t)
table_goldstein <- coef_table_function_lm(model = model_goldstein) %>%
  rename(`Goldstein method` = t)

# one large table
table_final <- left_join(table_episewer, table_ern, by = c("vars"))
table_final <- left_join(table_final, table_exp, by = c("vars"))
table_final <- left_join(table_final, table_fit_line, by = c("vars"))
table_final <- left_join(table_final, table_glm, by = c("vars"))
table_final <- left_join(table_final, table_huisman, by = c("vars"))
table_final <- left_join(table_final, table_wwsub, by = c("vars"))
table_final <- left_join(table_final, table_goldstein, by = c("vars"))

library(sjPlot)
tab_df(table_final,
       digits = 4,
       title = "Table: Regression results for Rt accuracy for percent agreement with above or below one.",
       file = "Table - regression Rt results.doc")
