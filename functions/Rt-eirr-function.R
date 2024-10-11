eirr_function <- function(dataframe, log, E_prior, I_prior, R1_prior, Rt_prior){
  
  colnames(dataframe) <- c("Date", "predictor")
  
  # add time change
  dataframe$week <- floor_date(dataframe$Date, unit = "week")
  dataframe$param_change_times <- seq(1:nrow(dataframe))
  param_change_times <- dataframe %>%
    arrange(desc(Date)) %>%
    filter(!duplicated(week))%>%
    select(param_change_times)
  param_change_times <- rev(param_change_times$param_change_times)
  
  # new observations
  dataframe$count <- seq(1:nrow(dataframe))
  obstimes <- as.numeric(dataframe$count) # numeric, triplicate
  
  # data (note if it needs to be log transformed)
  if(log == "Yes"){
    data <- log(dataframe$predictor)
    
  } else if(log == "No"){
    data <- dataframe$predictor
    
  }
  
  # choose to sample from prior or posterior
  priors_only <- FALSE
  
  # choose number of samples, number of chains, and seed
  n_samples <- 25L
  n_chains <- 4L
  seed <- 1L
  
  start_time <- Sys.time()
  
  posterior_samples_eirr <- fit_eirrc(data, 
                                      obstimes, 
                                      param_change_times, 
                                      priors_only, 
                                      n_samples, 
                                      n_chains, 
                                      seed,
                                      E_init_mean = E_prior,
                                      I_init_mean = I_prior,
                                      R1_init_mean = R1_prior,
                                      rt_init_mean = log(Rt_prior))
  # next steps
  # create dataframes of posterior/posterior predictive draws
  posterior_output_eirr <- generate_eirrc(posterior_samples_eirr,
                                          data,
                                          obstimes, 
                                          param_change_times,
                                          seed = seed,
                                          E_init_mean = E_prior,
                                          I_init_mean = I_prior,
                                          R1_init_mean = R1_prior,
                                          rt_init_mean = log(Rt_prior))
  
  # create quantiles of time-varying parameters
  eirr_quantiles <- make_timevarying_quantiles(posterior_output_eirr[[2]])
  
  eirr_rt_quantiles <- eirr_quantiles %>% dplyr::filter(name == "rt_t_values")
  
  end_time <- Sys.time()
  
  end_time - start_time
  
  # check for extra row*
  
  # seems to be an extra observation, maybe because it had an incomplete week
  eirr_rt_quantiles <- eirr_rt_quantiles %>%
    filter(time <= length(unique(dataframe$week)))
  
  # merge to time series
  week <- c(param_change_times, param_change_times, param_change_times)
  eirr_rt_quantiles$param_change_times <- week
  
  # select the 95% quantile and merge in the weeks
  eirr_rt_95 <- eirr_rt_quantiles %>%
    filter(.width == 0.95)
  
  weeks <- dataframe %>%
    select(week, param_change_times)
  
  eirr_rt_95 <- left_join(eirr_rt_95, weeks, by = c("param_change_times")) %>%
    mutate(se_rt = (.upper - .lower) / 3.02) %>%
    select(week, value, .lower, .upper, se_rt) %>%
    rename(mean_rt = value, 
           ll_95_rt = .lower,
           ul_95_rt = .upper)
  
  return(eirr_rt_95)
}
