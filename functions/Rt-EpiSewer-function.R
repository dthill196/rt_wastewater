###
### EPISEWER function
###

# LOAD PACKAGES
# link to github https://github.com/adrian-lison/EpiSewer?tab=readme-ov-file

# pre set up

# Install CMDSTANR

# follow these steps for a windows machine if you run into trouble
#https://www.maxmantei.com/2020/05/16/cmdstanr-windows/

#devtools::install_github("stan-dev/cmdstanr", force = TRUE)
library(cmdstanr)
#install_cmdstan()

# Steps
# 1 setup packages and stan connection

#remotes::install_github("adrian-lison/EpiSewer", dependencies = TRUE)
library(EpiSewer)
library(data.table)
library(ggplot2)
library(Matrix)
library(dplyr)
library(lubridate)

#cmdstanr::check_cmdstan_toolchain()
#cmdstanr::install_cmdstan()

rt_episewer_function <- function(dataframe, weekly, cores_use, chains, iter_warmup_n, iter_sampling_n){
  
  if(weekly == "No"){
    
    # fill missing flow data with the mean of the county
    dat <- dataframe %>%
      group_by(county) %>%
      mutate(flow_mean_for_missing = mean(mean_flow_w, na.rm = TRUE)
      ) %>%
      mutate(mean_flow_w = ifelse(is.na(mean_flow_w), flow_mean_for_missing, mean_flow_w)
      ) %>%
      ungroup()
    
    # data for our test county -> use winter surge: october to january 15, 2023
    cases <- dataframe %>%
      filter(county == county_name) %>%
      select(date, cases_new.7avg) %>%
      rename(cases = cases_new.7avg) %>%
      arrange(date)
    
    # concentration data
    conc_data <- dataframe %>%
      filter(county == county_name) %>%
      select(date, sars2.7avg) %>%
      mutate(weekday = weekdays(date)
      ) %>%
      rename(concentration = sars2.7avg)%>%
      filter(!is.na(concentration)) %>%
      arrange(date)
    
    # flow data
    flow <- dataframe %>%
      filter(county == county_name) %>%
      select(date, mean_flow_w) %>%
      rename(flow = mean_flow_w)%>%
      arrange(date)
    
    # combine into episewer format
    ww_data <- sewer_data(measurements = conc_data, flows = flow)
    
    # 2 generate assumptions, using defaults from EpiSewer main page
    
    # assumptions
    generation_dist <- get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4, maxX = 12)
    incubation_dist <- get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4, maxX = 12)
    shedding_dist <- get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397, maxX = 30)
    
    # shedding load per case (from Huisman method)
    load_per_case  <- suggest_load_per_case(
      conc_data,
      cases,
      #flow,
      # assume constant flow for our analysis
      flow_constant = mean(flow$flow),
      ascertainment_prop = 1
    )
    
    # assumptions object
    # combine assumptions
    ww_assumptions <- sewer_assumptions(
      generation_dist = generation_dist,
      incubation_dist = incubation_dist,
      shedding_dist = shedding_dist,
      load_per_case = load_per_case
    )
    
    # 4 make the estimate
    
    # option 1- decrease iterations per chain, increase chains and cores
    
    cores_use <- cores_use
    chains <- chains
    iter_warmup_n <- iter_warmup_n
    iter_sampling_n <- iter_warmup_n
    
    options(mc.cores = cores_use) # allow stan to use 4 cores, i.e. one for each chain
    
    ww_result <- EpiSewer(
      data = ww_data,
      assumptions = ww_assumptions,
      fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(iter_warmup = iter_warmup_n, iter_sampling = iter_sampling_n, chains = chains))
      ,
      # add overdispersion parameter for infections. better fit for covid data suggested by Adrian
      infections = model_infections(
        infection_noise = infection_noise_estimate(overdispersion = TRUE, overdispersion_prior_mu = 0.05, overdispersion_fixed = TRUE)
      ),
      # # optional variation in shedding by 10%
      shedding = model_shedding(
        load_variation = load_variation_estimate(cv_prior_mu = 0.1, cv_prior_sigma = 0.025)
      ),
      sewage = model_sewage(
        flows = flows_assume(flow_constant = mean(flow$flow, na.rm = TRUE))
      )
      
    )
    
    # extract to df
    episewer_rt <- as.data.frame(ww_result$summary$R)
    
    # filter to match input start date
    episewer_rt <- episewer_rt %>%
      filter(date >= min(dat$date))
    
    # calculate SE from the 95 % error bars
    # E = (upper limit – lower limit) / 3.92.
    episewer_rt$se_rt <- (episewer_rt$upper_0.95 - episewer_rt$lower_0.95) / 3.92
    
    # rename values and keep what we want
    episewer_rt <- episewer_rt %>%
      rename(mean_rt = mean,
             ll_95_rt = lower_0.95, 
             ul_95_rt = upper_0.95) %>%
      select(date, mean_rt, se_rt, ll_95_rt, ul_95_rt) %>%
      mutate(week = floor_date(date, unit = "week")
      )
    
    return(episewer_rt)
  } else if(weekly == "Yes"){
    # fill missing flow data with the mean of the county
    dat <- dataframe %>%
      group_by(county) %>%
      mutate(flow_mean_for_missing = mean(mean_flow_w, na.rm = TRUE)
      ) %>%
      mutate(mean_flow_w = ifelse(is.na(mean_flow_w), flow_mean_for_missing, mean_flow_w)
      ) %>%
      ungroup()
    
    # data for our test county -> use winter surge: october to january 15, 2023
    cases <- dataframe %>%
      filter(county == county_name) %>%
      select(date, cases_new.7avg) %>%
      rename(cases = cases_new.7avg) %>%
      arrange(date)
    
    # concentration data
    conc_data <- dataframe %>%
      filter(county == county_name) %>%
      select(date, sars2.7avg) %>%
      mutate(weekday = weekdays(date)
      ) %>%
      rename(concentration = sars2.7avg)%>%
      filter(!is.na(concentration)) %>%
      arrange(date)
    
    # flow data
    flow <- dataframe %>%
      filter(county == county_name) %>%
      select(date, mean_flow_w) %>%
      rename(flow = mean_flow_w)%>%
      arrange(date)
    
    # combine into episewer format
    ww_data <- sewer_data(measurements = conc_data, flows = flow)
    
    # 2 generate assumptions, using defaults from EpiSewer main page
    
    # assumptions
    generation_dist <- get_discrete_gamma_shifted(gamma_mean = 3, gamma_sd = 2.4, maxX = 12)
    incubation_dist <- get_discrete_gamma(gamma_shape = 8.5, gamma_scale = 0.4, maxX = 12)
    shedding_dist <- get_discrete_gamma(gamma_shape = 0.929639, gamma_scale = 7.241397, maxX = 30)
    
    # shedding load per case (from Huisman method)
    load_per_case  <- suggest_load_per_case(
      conc_data,
      cases,
      #flow,
      # assume constant flow for our analysis
      flow_constant = mean(flow$flow),
      ascertainment_prop = 1
    )
    
    # assumptions object
    # combine assumptions
    ww_assumptions <- sewer_assumptions(
      generation_dist = generation_dist,
      incubation_dist = incubation_dist,
      shedding_dist = shedding_dist,
      load_per_case = load_per_case
    )
    
    # 4 make the estimate
    
    # option 1- decrease iterations per chain, increase chains and cores
    
    cores_use <- cores_use
    chains <- chains
    iter_warmup_n <- iter_warmup_n
    iter_sampling_n <- iter_warmup_n
    
    options(mc.cores = cores_use) # allow stan to use 4 cores, i.e. one for each chain
    
    ww_result <- EpiSewer(
      data = ww_data,
      assumptions = ww_assumptions,
      fit_opts = set_fit_opts(sampler = sampler_stan_mcmc(iter_warmup = iter_warmup_n, iter_sampling = iter_sampling_n, chains = chains))
      ,
      # add overdispersion parameter for infections. better fit for covid data suggested by Adrian
      infections = model_infections(
        infection_noise = infection_noise_estimate(overdispersion = TRUE, overdispersion_prior_mu = 0.05, overdispersion_fixed = TRUE)
      ),
      # # optional variation in shedding by 10%
      shedding = model_shedding(
        load_variation = load_variation_estimate(cv_prior_mu = 0.1, cv_prior_sigma = 0.025)
      ),
      sewage = model_sewage(
        flows = flows_assume(flow_constant = mean(flow$flow, na.rm = TRUE))
      )
      
    )
    
    # extract to df
    episewer_rt <- as.data.frame(ww_result$summary$R)
    
    # filter to match input start date
    episewer_rt <- episewer_rt %>%
      filter(date >= min(dat$date))
    
    # calculate SE from the 95 % error bars
    # E = (upper limit – lower limit) / 3.92.
    episewer_rt$se_rt <- (episewer_rt$upper_0.95 - episewer_rt$lower_0.95) / 3.92
    
    # rename values and keep what we want
    episewer_rt <- episewer_rt %>%
      rename(mean_rt = mean,
             ll_95_rt = lower_0.95, 
             ul_95_rt = upper_0.95) %>%
      select(date, mean_rt, se_rt, ll_95_rt, ul_95_rt) %>%
      mutate(week = floor_date(date, unit = "week")
      )
    
    # weekly means for the final data
    final_weekly <- episewer_rt %>%
      group_by(week) %>%
      summarize(mean_rt = mean(mean_rt, na.rm = TRUE),
                ll_95_rt = mean(ll_95_rt, na.rm = TRUE),
                ul_95_rt = mean(ul_95_rt, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      mutate(se_rt = (ul_95_rt - ul_95_rt / 3.92) ) %>%
      select(week, mean_rt, se_rt, ll_95_rt, ul_95_rt)
    
    return(final_weekly)
  }
  
}