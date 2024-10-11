# ERN method

Rt_ERN_function <- function(dataframe, weekly){
  
  # assign span based on size of df
  if(length(dataframe$date) < 200){
    span = 0.2
  } else if(length(dataframe$date) >= 200 & length(dataframe$date) < 500){
    span = 0.15
  } else if(length(dataframe$date) >= 500){
    span = 0.1
  }
  
  # change name of case data
  colnames(dataframe) <- c("date", "value", "cases_new.7avg")
  
  # compute weekly mean for weekly output
  if(weekly =="Yes"){
    
    # Define SARS-CoV-2 fecal shedding and generation interval distributions
    dist.fec = ern::def_dist(
      dist = "gamma",
      mean = 12.90215,
      mean_sd = 1.136829,
      shape = 1.759937,
      shape_sd = 0.2665988,
      max = 33
    )
    
    dist.gi  = ern::def_dist(
      dist     = "gamma",
      mean     = 6.84,
      mean_sd  = 0.7486,
      shape    = 2.39,
      shape_sd = 0.3573,
      max      = 15
    )
    
    # scaling factor
    # ratio of ww conc to new cases
    model <- glm(cases_new.7avg ~ log(value), data = dataframe, family = "poisson")
    scaling.factor = exp(model$coefficients[2])
    
    # remove case data from df
    ww.conc <- dataframe %>%
      dplyr::select(-cases_new.7avg)
    
    # Initializing smoothing parameters
    prm.smooth = list(
      align  = 'center', # smoothing alignment
      method = 'loess',  # smoothing method
      span   = span,     # smoothing span (used for loess smoothing only)
      floor  = 5        # minimum smoothed concentration value (optional, loess smoothing only)
    )
    
    
    # Initialzing Rt settings
    prm.R = list(
      iter   = 20,   # number of iterations in Rt ensemble
      CI     = 0.95, # confidence interval
      window = 10,   # Time window for Rt calculations
      config.EpiEstim = NULL # optional EpiEstim configuration for Rt calculations
    )
    
    # estimate rt
    r.estim = ern::estimate_R_ww(
      ww.conc        = ww.conc,
      dist.fec       = dist.fec,
      dist.gi        = dist.gi,
      scaling.factor = scaling.factor,
      prm.smooth     = prm.smooth,
      prm.R = prm.R,
      silent = TRUE # suppress output messages
    )
    
    # return rt dataframe with necessary columns
    rt_df <- r.estim$R
    
    # add se
    rt_df$se_rt <- (rt_df$upr - rt_df$lwr) / 3.92
    
    # rename fields
    rt_df <- rt_df %>%
      rename(mean_rt = mean,
             ll_95_rt = lwr,
             ul_95_rt = upr) %>%
      mutate(week = lubridate::floor_date(date, unit = "week"))
    
    # weekly means for the final data
    final_weekly <- rt_df %>%
      group_by(week) %>%
      summarize(mean_rt = mean(mean_rt, na.rm = TRUE),
                ll_95_rt = mean(ll_95_rt, na.rm = TRUE),
                ul_95_rt = mean(ul_95_rt, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      mutate(se_rt = (ul_95_rt - ul_95_rt / 3.92) ) %>%
      select(week, mean_rt, se_rt, ll_95_rt, ul_95_rt)
    
    return(final_weekly)
    
  } else if(weekly == "No"){
    
    
    # Define SARS-CoV-2 fecal shedding and generation interval distributions
    dist.fec = ern::def_dist(
      dist = "gamma",
      mean = 12.90215,
      mean_sd = 1.136829,
      shape = 1.759937,
      shape_sd = 0.2665988,
      max = 33
    )
    
    dist.gi  = ern::def_dist(
      dist     = "gamma",
      mean     = 6.84,
      mean_sd  = 0.7486,
      shape    = 2.39,
      shape_sd = 0.3573,
      max      = 15
    )
    
    # scaling factor
    # ratio of ww conc to new cases
    model <- glm(cases_new.7avg ~ log(value), data = dataframe, family = "poisson")
    scaling.factor = exp(model$coefficients[2])
    
    # remove case data from df
    ww.conc <- dataframe %>%
      dplyr::select(-cases_new.7avg)
    
    # Initializing smoothing parameters
    prm.smooth = list(
      align  = 'center', # smoothing alignment
      method = 'loess',  # smoothing method
      span   = span,     # smoothing span (used for loess smoothing only)
      floor  = 5        # minimum smoothed concentration value (optional, loess smoothing only)
    )
    
    
    # Initialzing Rt settings
    prm.R = list(
      iter   = 20,   # number of iterations in Rt ensemble
      CI     = 0.95, # confidence interval
      window = 10,   # Time window for Rt calculations
      config.EpiEstim = NULL # optional EpiEstim configuration for Rt calculations
    )
    
    # estimate rt
    r.estim = ern::estimate_R_ww(
      ww.conc        = ww.conc,
      dist.fec       = dist.fec,
      dist.gi        = dist.gi,
      scaling.factor = scaling.factor,
      prm.smooth     = prm.smooth,
      prm.R = prm.R,
      silent = TRUE # suppress output messages
    )
    
    # return rt dataframe with necessary columns
    rt_df <- r.estim$R
    
    # add se
    rt_df$se_rt <- (rt_df$upr - rt_df$lwr) / 3.92
    
    # rename fields
    rt_df <- rt_df %>%
      rename(mean_rt = mean,
             ll_95_rt = lwr,
             ul_95_rt = upr) %>%
      mutate(week = lubridate::floor_date(date, unit = "week"))
    
   return(rt_df)
    
  }
  
  
}
