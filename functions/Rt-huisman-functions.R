###########################################################
# Function to estimate Re from WW
# Author: JS Huisman
###########################################################
# These are mostly wrapper functions around the 
# deconvolution and Re estimation functions from 
# https://github.com/covid-19-Re/shiny-dailyRe

###########################################################
# Delay/Shedding Load Distributions #####

# find gamma parameters from mean/sd of distribution
getGammaParams <- function(meanParam, sdParam){
  shapeParam <- meanParam^2 / (sdParam^2)
  scaleParam <- (sdParam^2) / meanParam
  return(list(shape = shapeParam, scale = scaleParam))
}

#find mean/sd from gamma shape/scale parameters
getInvGammaParams <- function(shapeParam, scaleParam){
  meanParam <- scaleParam * shapeParam
  sdParam <- sqrt(scaleParam^2 * shapeParam)
  return(list(mean = meanParam, sd = sdParam))
}

# specific distributions we use in the paper
getCountParams <- function(obs_type){
  switch(obs_type,
         incubation = getGammaParams(5.3, 3.2),
         zero = list(shape = 0, scale = 0),
         #confirmed = getGammaParams(5.5, 3.8),
         confirmed_zh = getGammaParams(4, 4),#confirmed_zh = getGammaParams(2.83, 2.96), # they are using zurich and california values, we may want to calculate our own
         confirmed_cali = getGammaParams(4.51, 3.16),
         death = getGammaParams(15.0, 6.9),
         han = getGammaParams(4.7, 1.7),
         wolfel = getGammaParams(8.6, 0.9),
         benefield = list(shape = 0.929639, scale = 7.241397))
}

###########################################################
# Deconvolve #####

# To achieve an input format required by the 
# get_infection_incidence_by_deconvolution function
addUselessColumns <- function(df, inc_var = 'n1'){
  if (!'region' %in% colnames(df)){
    df <- df %>%
      mutate(region = 'CHE')
  }
  
  observation_df <- df %>%
    dplyr::select(date, region, value = inc_var) %>%
    mutate(data_type = inc_var,
           source = 'ETH',
           variable = 'incidence',
           country = 'Switzerland',
           date_type = 'report',
           local_infection = TRUE)
  
  observation_df <- observation_df %>%
    filter(!is.na(value))
  
  return(observation_df)
}

# wrapper around the actual deconvolution function
deconvolveIncidence <- function(df, incidence_var = 'n1',
                                IncubationParams, OnsetToCountParams,
                                smooth_param = FALSE, n_boot = 50){
  infection_df <- addUselessColumns(df, inc_var = incidence_var)
  
 constant_delay_distributions <- list("Simulated" = get_vector_constant_waiting_time_distr(
  IncubationParams$shape, IncubationParams$scale,
   OnsetToCountParams$shape, OnsetToCountParams$scale),
    "Symptoms" = get_vector_constant_waiting_time_distr(
      IncubationParams$shape, IncubationParams$scale,
    0, 0))
  
  estimatedInfections <- get_infection_incidence_by_deconvolution(
    infection_df,
    is_local_cases = T,
    constant_delay_distribution = constant_delay_distributions[['Simulated']],
    constant_delay_distribution_incubation = constant_delay_distributions[["Symptoms"]],
    max_iterations = 100,
    smooth_incidence = smooth_param,
    empirical_delays = tibble(),
    n_bootstrap = n_boot,
    verbose = FALSE)
  
  return(estimatedInfections)
}

###########################################################
# Estimate Re ####

# wrapper around the Re estimation function
getReBootstrap <- function(deconvoluted_data){
  
  all_delays <- lapply(unique(deconvoluted_data$data_type), function(x){ c(Cori = 0)})
  names(all_delays) <- unique(deconvoluted_data$data_type)
  
  truncations <- list(left = c(Cori = 5),
                      right = c(Cori = 0))
  
  rawReEstimates <- suppressWarnings(doAllReEstimations(
    deconvoluted_data,
    slidingWindow = 3,
    methods = c("Cori"),
    variationTypes = c("slidingWindow"),
    all_delays,
    truncations,
    interval_ends = list()) )
  
  cleanEstimates <- cleanCountryReEstimate(rawReEstimates, method = 'bootstrap',
                                           rename_types = F, alpha=0.95)
  
  return(cleanEstimates)
}

###########################################################
compareTraces <- function(Re_i, Re_j){
  compare_df = Re_i %>%
    left_join(Re_j, by = 'date', suffix = c(".i", ".j")) %>%
    mutate(se = (median_R_mean.i - median_R_mean.j)^2,
           rele = abs((median_R_mean.j - median_R_mean.i)/median_R_mean.j),
           coverage = (median_R_mean.i > median_R_lowHPD.j) & (median_R_mean.i < median_R_highHPD.j) ) 
  
  se = compare_df %>% pull(se)
  rele = compare_df %>% pull(rele)
  coverage = compare_df %>% pull(coverage) %>% sum(na.rm = T) /length(Re_i$date)
  
  rmse = sqrt(sum(se, na.rm = T)/length(Re_i$date))
  mape = sum(rele, na.rm = T)/length(Re_i$date)
  
  return(c(rmse, coverage, mape))
}

############ additional functions
get_vector_constant_waiting_time_distr <- function(shape_incubation,
                                                   scale_incubation,
                                                   shape_onset_to_report,
                                                   scale_onset_to_report,
                                                   length_out = 200,
                                                   n_random_samples = 1E6) {
  
  F_h <- make_ecdf_from_gammas(shape = c(shape_incubation, shape_onset_to_report), scale = c(scale_incubation, scale_onset_to_report))
  
  f <- Vectorize(function(x){
    if(x < 0) {
      return(0)
    } else if(x < 0.5) {
      return(F_h(0.5))
    } else {
      return(F_h(round(x + 1E-8) + 0.5) - F_h(round(x + 1E-8) - 0.5))
    }
  })
  
  x <- 0:(length_out - 1)
  
  return(f(x))
}
get_infection_incidence_by_deconvolution <- function(
  data_subset,
  constant_delay_distribution,
  constant_delay_distribution_incubation = c(),
  is_onset_data = F,
  is_local_cases = T,
  smooth_incidence = T,
  days_incl = 21,
  empirical_delays  = tibble(),
  n_bootstrap = 5,
  days_further_in_the_past = 30,
  days_further_in_the_past_incubation = 5,
  max_iterations = 100,
  verbose = FALSE) {
  
  #TODO make the days_further_in_the_past type specific
  
  if(nrow(data_subset) == 0) {
    return(tibble())
  }
  
  data_type_subset <- unique(data_subset$data_type)[1]
  
  # exclude leading zeroes
  data_subset <- data_subset %>%
    arrange(date) %>%
    filter(cumsum(value) > 0)
  
  if(nrow(data_subset) == 0) {
    return(tibble())
  }
  
  minimal_date <- min(data_subset$date) - days_further_in_the_past
  maximal_date <- max(data_subset$date)
  all_dates <- seq(minimal_date, maximal_date, by = "days")
  
  is_empirical = (nrow(empirical_delays) > 0)
  
  if(verbose && is_empirical) {
    cat("\tEmpirical delay distribution available\n")
  }
  
  if( is_onset_data ) {
    delay_distribution_matrix_incubation <- get_matrix_constant_waiting_time_distr(
      constant_delay_distribution_incubation,
      all_dates)
    
    initial_delta_incubation <- min(which(cumsum(constant_delay_distribution_incubation) > 0.5)) - 1 # take median value (-1 because index 1 corresponds to zero days)
    
    
    # account for additional right-truncation of onset data (needs to be reported first)
    if(is_empirical) {
      delay_distribution_matrix_onset_to_report <- get_matrix_empirical_waiting_time_distr(
        empirical_delays,
        seq.Date(min(data_subset$date), max(data_subset$date), by = "days"))
    } else {
      delay_distribution_matrix_onset_to_report <- get_matrix_constant_waiting_time_distr(
        constant_delay_distribution,
        seq.Date(min(data_subset$date), max(data_subset$date), by = "days"))
    }
    
    data_subset <- data_subset %>%
      complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0))
    
    Q_vector_onset_to_report <- apply(delay_distribution_matrix_onset_to_report, MARGIN = 2, sum)
    
    if(unique(data_subset$region)[1] == "ESP") { # hack to work around spanish data between symptom onset dates only
      right_truncation <- 3
      # need to offset the Q vector by how many days were truncated off originally
      Q_vector_onset_to_report <- c(rep(1, right_truncation), Q_vector_onset_to_report[1:(length(Q_vector_onset_to_report) - right_truncation)] )
    }
    
    data_subset <- data_subset %>%
      mutate(value = value / Q_vector_onset_to_report) %>% 
      mutate(value = if_else(value == Inf, 0, value))
    
  } else {
    if(is_empirical) {
      delay_distribution_matrix_onset_to_report <- get_matrix_empirical_waiting_time_distr(
        empirical_delays,
        all_dates[(days_further_in_the_past_incubation + 1):length(all_dates)])
      
      delay_distribution_matrix_incubation <- get_matrix_constant_waiting_time_distr(
        constant_delay_distribution_incubation,
        all_dates)
      
      initial_delta_incubation <- min(which(cumsum(constant_delay_distribution_incubation) > 0.5)) - 1 # take median value (-1 because index 1 corresponds to zero days)
      initial_delta_report <-  median(empirical_delays$delay, na.rm = T)
    } else {
      delay_distribution_matrix <- get_matrix_constant_waiting_time_distr(
        constant_delay_distribution,
        all_dates)
      
      initial_delta <- min(which(cumsum(constant_delay_distribution) > 0.5)) - 1 # take median value (-1 because index 1 corresponds to zero days)
    }
  }
  
  
  
  results <- list(tibble())
  
  for (bootstrap_replicate_i in 0:n_bootstrap) {
    
    if (verbose == T) {
      cat("    Bootstrap replicate: ", bootstrap_replicate_i, "\n")
    }
    
    if (bootstrap_replicate_i == 0) {
      time_series <- data_subset
    } else {
      time_series <- get_bootstrap_replicate(data_subset)
    }
    
    if (smooth_incidence == T) {
      smoothed_incidence_data <- time_series %>%
        complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>% 
        mutate(value = getLOESSCases(dates = date, count_data = value, days_incl))
      
      raw_total_incidence <- sum(time_series$value, na.rm = TRUE)
      smoothed_total_incidence <- sum(smoothed_incidence_data$value, na.rm = T)
      
      if (smoothed_total_incidence > 0) {
        smoothed_incidence_data <- smoothed_incidence_data %>%
          mutate(value = value * raw_total_incidence / smoothed_total_incidence)
      }
      
    } else {
      smoothed_incidence_data <- time_series  %>%
        complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0))
    }
    
    
    if (is_onset_data) {
      deconvolved_infections <-  do_deconvolution(smoothed_incidence_data,
                                                  delay_distribution_matrix = delay_distribution_matrix_incubation,
                                                  days_further_in_the_past = days_further_in_the_past,
                                                  initial_delta = initial_delta_incubation,
                                                  max_iterations = max_iterations,
                                                  verbose = verbose)
    } else {
      if(is_empirical) {
        # perform the deconvolution in two steps
        deconvolved_symptom_onsets <- do_deconvolution(smoothed_incidence_data,
                                                       delay_distribution_matrix = delay_distribution_matrix_onset_to_report,
                                                       days_further_in_the_past = days_further_in_the_past - days_further_in_the_past_incubation,
                                                       initial_delta = initial_delta_report,
                                                       max_iterations = max_iterations,
                                                       verbose = verbose)
        
        deconvolved_infections <- do_deconvolution(deconvolved_symptom_onsets,
                                                   delay_distribution_matrix = delay_distribution_matrix_incubation,
                                                   days_further_in_the_past = days_further_in_the_past_incubation,
                                                   initial_delta = initial_delta_incubation,
                                                   max_iterations = max_iterations,
                                                   verbose = verbose)
      } else {
        deconvolved_infections <-  do_deconvolution(smoothed_incidence_data,
                                                    delay_distribution_matrix = delay_distribution_matrix,
                                                    days_further_in_the_past = days_further_in_the_past,
                                                    initial_delta = initial_delta,
                                                    max_iterations = max_iterations,
                                                    verbose = verbose)
      }
    }
    
    
    deconvolved_infections <- deconvolved_infections %>% slice((days_further_in_the_past -5 + 1):n())
    
    data_type_name <- paste0("infection_", data_type_subset)
    
    ## dataframe containing results
    deconvolved_infections <- tibble(
      date = deconvolved_infections$date,
      region = unique(time_series$region)[1],
      country = unique(time_series$country)[1],
      source = unique(time_series$source)[1],
      local_infection = is_local_cases,
      data_type = data_type_name,
      replicate = bootstrap_replicate_i,
      value = deconvolved_infections$value
    )
    
    results <- c(results, list(deconvolved_infections))
  }
  
  return(bind_rows(results))
}
make_ecdf_from_gammas <- function(shape, scale, numberOfSamples = 1E6) {
  draws <-
    rgamma(numberOfSamples, shape = shape[1], scale = scale[1]) +
    rgamma(numberOfSamples, shape = shape[2], scale = scale[2])
  return(Vectorize(ecdf(draws)))
}

get_matrix_constant_waiting_time_distr <- function(waiting_time_distr,
                                                   all_dates) {
  N <- length(all_dates)
  
  if(length(all_dates) >= length(waiting_time_distr)) {
    waiting_time_distr <- c(waiting_time_distr, rep(0, times = N - length(waiting_time_distr)))
  }
  
  delay_distribution_matrix <- matrix(0, nrow = N, ncol = N)
  for(i in 1:N) {
    delay_distribution_matrix[, i ] <-  c(rep(0, times = i - 1 ), waiting_time_distr[1:(N - i + 1)])
  }
  
  return(delay_distribution_matrix)
}
getLOESSCases <- function(dates, count_data, days_incl = 21, degree = 1, truncation = 0) {
  
  if (truncation != 0) {
    dates <- dates[1:(length(dates) - truncation)]
    count_data <- count_data[1:(length(count_data) - truncation)]
  }
  
  n_points <- length(unique(dates))
  sel_span <- days_incl / n_points
  
  n_pad <- round(length(count_data) * sel_span * 0.5)
  
  c_data <- data.frame(value = c(rep(0, n_pad), count_data),
                       date_num = c(seq(as.numeric(dates[1]) - n_pad, as.numeric(dates[1]) - 1),
                                    as.numeric(dates)))
  c_data.lo <- loess(value ~ date_num, data = c_data, span = sel_span, degree = degree)
  smoothed <- predict(c_data.lo)
  smoothed[smoothed < 0] <- 0
  raw_smoothed_counts <- smoothed[(n_pad + 1):length(smoothed)]
  normalized_smoothed_counts <-
    raw_smoothed_counts * sum(count_data, na.rm = T) / sum(raw_smoothed_counts, na.rm = T)
  
  if (truncation != 0) {
    normalized_smoothed_counts <- append(normalized_smoothed_counts, rep(NA, truncation))
  }
  return(normalized_smoothed_counts)
}
do_deconvolution <- function(
  incidence_data,
  days_further_in_the_past = 30,
  verbose = FALSE,
  delay_distribution_matrix,
  initial_delta,
  max_iterations = 100
) {
  
  # use mode of 'constant_delay_distribution'. -1 because indices are offset by one as the delay can be 0.
  
  first_guess_delay <- ceiling(initial_delta)
  
  if (verbose) {
    cat("\tDelay on first guess: ", first_guess_delay, "\n")
  }
  
  first_recorded_incidence <-  with(filter(incidence_data, cumsum(value) > 0), value[which.min(date)])
  last_recorded_incidence <- with(incidence_data, value[which.max(date)])
  
  minimal_date <- min(incidence_data$date) - days_further_in_the_past
  maximal_date <- max(incidence_data$date)
  
  first_guess <- incidence_data %>%
    mutate(date = date - first_guess_delay) %>%
    complete(date = seq.Date(minimal_date, min(date), by = "days"),
             fill = list(value = first_recorded_incidence)) %>% # left-pad with first recorded value
    complete(date = seq.Date(max(date), maximal_date, by = "days"),
             fill = list(value = last_recorded_incidence)) %>% # right-pad with last recorded value
    arrange(date) %>% 
    filter(date >=  minimal_date)
  
  original_incidence <- incidence_data %>% 
    complete(date = seq.Date(minimal_date, maximal_date, by = "days"),
             fill = list(value = 0)) %>% 
    pull(value)
  
  final_estimate <- iterate_RL(
    first_guess$value,
    original_incidence,
    delay_distribution_matrix = delay_distribution_matrix,
    max_delay = days_further_in_the_past,
    max_iterations = max_iterations,
    verbose = verbose)
  
  deconvolved_dates <- first_guess %>% pull(date)
  
  result <- tibble(date = deconvolved_dates, value = final_estimate)
  
  result <- result %>%
    filter(date <= maximal_date - first_guess_delay)
  
  return(result)
}

iterate_RL <- function(
  initial_estimate,
  original_incidence,
  delay_distribution_matrix,
  threshold_chi_squared = 1,
  max_iterations = 100,
  max_delay,
  verbose = FALSE) {
  
  current_estimate <- initial_estimate
  N <- length(current_estimate)
  N0 <- N - max_delay
  chi_squared <- Inf
  count <- 1
  
  delay_distribution_matrix <- delay_distribution_matrix[1:length(current_estimate), 1:length(current_estimate)]
  truncated_delay_distribution_matrix <- delay_distribution_matrix[(1 + max_delay):NROW(delay_distribution_matrix),, drop = F]
  
  Q_vector <- apply(truncated_delay_distribution_matrix, MARGIN = 2, sum)
  
  while(chi_squared > threshold_chi_squared & count <= max_iterations) {
    
    if (verbose) {
      cat("\t\tStep: ", count, " - Chi squared: ", chi_squared, "\n")
    }
    
    E <- as.vector(delay_distribution_matrix %*% current_estimate)
    B <- replace_na(original_incidence/E, 0)
    
    current_estimate <- current_estimate / Q_vector *  as.vector(crossprod(B, delay_distribution_matrix))
    current_estimate <- replace_na(current_estimate, 0)
    
    chi_squared <- 1/N0 * sum((E[(max_delay + 1): length(E)] - original_incidence[(max_delay + 1) : length(original_incidence)])^2/E[(max_delay + 1): length(E)], na.rm = T)
    count <- count + 1
  }
  
  return(current_estimate)
}

get_bootstrap_replicate <- function(original_time_series, block_size = 10, days_incl = 21) {
  tmp <- original_time_series
  
  # Change introduced after meeting on 19.1
  #tmp$log_value <- ifelse(tmp$value != 0, log(tmp$value), 0)
  tmp$log_value <- log(tmp$value + 1)
  
  smoothed_incidence_data <- tmp %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(log_value = 0)) %>%
    mutate(log_loess = getLOESSCases(dates = date, count_data = log_value, days_incl),
           log_diff = log_value - log_loess)
  
  log_diff_boot <- block_boot_overlap_func(smoothed_incidence_data$log_diff, block_size)
  log_smoothed_data <- smoothed_incidence_data$log_loess
  
  ts_boot <- exp(log_diff_boot + log_smoothed_data) -1
  ts_boot[ts_boot<0] <- 0
  ts_boot <- round(ts_boot)
  
  replicate <- original_time_series %>%
    complete(date = seq.Date(min(date), max(date), by = "days"), fill = list(value = 0)) %>%
    dplyr::mutate(value = ts_boot) %>%
    arrange(date)
  
  return(replicate)
}
block_boot_overlap_func <- function(ts, block_size = 10){
  
  # get the weekdays for each position of ts
  weekdays_index <- (1:length(ts)) %% 7
  weekdays_index[which(weekdays_index==0)] <- 7
  
  ts_boot <-c()
  last_day_index <- 7
  
  ###### get the ts_boot: make sure glue wrt the correct days
  while(length(ts_boot) < length(ts)){
    start_index <- sample(1:(length(ts)-block_size+1), 1)
    sampled_index <- start_index:(start_index+block_size-1)
    sampled_weekdays <- weekdays_index[sampled_index]
    
    # make sure the day related to the first sample is after the previous ts_boot
    first_day_index <- which(sampled_weekdays==last_day_index)[1] + 1
    ts_boot_index <- sampled_index[first_day_index:block_size]
    
    last_day_index <- tail(weekdays_index[ts_boot_index],1)
    
    ts_boot <- c(ts_boot, ts[ts_boot_index])
  }
  
  # take the same length as previous ts
  ts_boot <- ts_boot[1:length(ts)]
  
  return(ts_boot)
}

doAllReEstimations <- function(
  data,
  slidingWindow = 3,
  methods = c("Cori", "WallingaTeunis"),
  variationTypes = c("step", "slidingWindow"),
  all_delays,
  truncations,
  interval_ends = list(default = c("2020-04-01")),
  ...) {
  
  results_list <- list()
  
  for (source_i in unique(data$source)) {
    cat("estimating Re for data source: ", source_i, "...\n")
    for (region_i in unique(data$region)) {
      cat("  Region: ", region_i, "\n")
      
      if ("list" %in% class(interval_ends)) {
        if (!is.null(interval_ends[[region_i]])) {
          region_interval_ends <- interval_ends[[region_i]]
        } else {
          region_interval_ends <- interval_ends[["default"]]
        }
      } else if ("Date" %in% class(interval_ends)) {
        region_interval_ends <- interval_ends
      } else {
        warning(str_c(
          "no valid interval ends for region ", region_i, ". ",
          "Interval ends must be a vector of dates or a named list with names corresponding to regions",
          "(or \"default\")."
        ))
        region_interval_ends <- ""
      }
      if (length(region_interval_ends) == 0) {
        region_interval_ends <- "01-01-2020"
      }
      ## Run EpiEstim
      for (data_type_i in unique(data$data_type)) {
        subset_data <- data %>% filter(region == region_i & source == source_i & data_type == data_type_i)
        if (nrow(subset_data) == 0) {
          next
        }
        cat("    Data type: ", data_type_i, "\n")
        
        delay_i <- all_delays[[data_type_i]]
        
        for (replicate_i in unique(unique(subset_data$replicate))) {
          subset_data_rep <- subset(subset_data, subset_data$replicate == replicate_i)
          results_list <- c(results_list,
                            list(
                              doReEstimation(
                                subset_data_rep,
                                slidingWindow = slidingWindow,
                                methods = methods,
                                variationTypes = variationTypes,
                                interval_ends = region_interval_ends,
                                delays = delay_i,
                                truncations = truncations
                              )
                            )
          )
        }
      }
    }
  }
  
  return(bind_rows(results_list))
}

doReEstimation <- function(
  data_subset,
  slidingWindow = 1,
  methods,
  variationTypes,
  interval_ends = c("2020-04-01"),
  delays,
  truncations) {
  
  end_result <-  data.frame()
  
  for (method_i in methods) {
    for (variation_i in variationTypes) {
      
      if(nrow(data_subset %>% filter(local_infection == F)) > 0) {
        incidence_data_local <- data_subset %>% filter(local_infection == T) %>% pull(value)
        incidence_data_import <- data_subset %>% filter(local_infection == F) %>% pull(value)
        
        incidence_data <- data.frame(local = incidence_data_local,
                                     imported = incidence_data_import)
      } else {
        incidence_data <- data.frame(I = data_subset %>% filter(local_infection == T) %>% pull(value))
      }
      
      dates <- data_subset %>% filter(local_infection == T) %>% pull(date)
      
      offsetting <- delays[method_i]
      
      leftTrunc <- truncations$left[method_i]
      rightTrunc <- truncations$right[method_i]
      
      result <- estimateRe(
        dates = dates,
        incidenceData = incidence_data,
        windowLength =  slidingWindow,
        estimateOffsetting = offsetting,
        rightTruncation = rightTrunc,
        leftTruncation = leftTrunc,
        method = method_i,
        variationType = variation_i,
        interval_ends = interval_ends)
      if (nrow(result) > 0) {
        result$region <- unique(data_subset$region)[1]
        result$country <- unique(data_subset$country)[1]
        result$source <- unique(data_subset$source)[1]
        result$data_type <- unique(data_subset$data_type)[1]
        result$replicate <- unique(data_subset$replicate)[1]
        ## need to reorder columns in 'results' dataframe to do the same as in data
        result <- result[, c(
          "date", "region", "country", "source", "data_type", "estimate_type",
          "replicate", "value", "variable")]
        end_result <- bind_rows(end_result, result)
      }
    }
  }
  
  return(end_result)
}

estimateRe <- function(
  dates,
  incidenceData,
  estimateOffsetting = 10,
  rightTruncation = 0,
  leftTruncation = 5,
  method = "Cori",
  variationType = "slidingWindow",
  interval_ends = c("2020-03-13", "2020-03-16", "2020-03-20"),
  minimumCumul = 5,
  windowLength= 4,
  mean_si = 4.8,
  std_si  = 2.3) {
  
  offset <- 1
  cumulativeIncidence <- 0
  while (cumulativeIncidence < minimumCumul) {
    if (offset > nrow(incidenceData)) {
      return(data.frame(date = c(), variable = c(), value = c(), estimate_type = c()))
    }
    cumulativeIncidence <- cumulativeIncidence + incidenceData[offset, 1]
    offset <- offset + 1
  }
  
  ## offset needs to be at least two for EpiEstim
  offset <- max(2, offset)
  
  rightBound <- nrow(incidenceData) - (windowLength - 1)
  
  if (rightBound < offset) { ## no valid data point, return empty estimate
    return(data.frame(date = c(), variable = c(), value = c(), estimate_type = c()))
  }
  
  ## generate start and end bounds for Re estimates
  if (variationType == "step") {
    
    # index in incidenceData that corresponds to the interval_end date
    interval_end_indices <- sapply(
      interval_ends,
      function(x) {
        which(dates == as.Date(x))[1]
      }
    )
    
    #starts and end indices of the intervals (numeric vector)
    # t_start = interval_end + 1
    t_start <- c(offset, na.omit(interval_end_indices) + 1)
    t_end <- c(na.omit(interval_end_indices), nrow(incidenceData))
    
    if (offset >= nrow(incidenceData)) {
      return(data.frame(date = c(), variable = c(), value = c(), estimate_type = c()))
    }
    
    # remove intervals if the offset is greater than the
    # end of the interval
    while (offset > t_end[1]) {
      t_start <- t_start[-1]
      t_start[1] <- offset
      t_end <- t_end[-1]
    }
    
    # make sure there are no intervals beyond the length of the data
    while (t_start[length(t_start)] >= nrow(incidenceData)) {
      t_end <- t_end[-length(t_end)]
      t_start <- t_start[-length(t_start)]
    }
    
    outputDates <- dates[t_start[1]:t_end[length(t_end)]]
    
  } else if (variationType == "slidingWindow") {
    # computation intervals corresponding to every position of the
    # sliding window
    t_start <- seq(offset, rightBound)
    t_end <- t_start + windowLength - 1
    outputDates <- dates[t_end]
  } else {
    print("Unknown time variation.")
    return(data.frame(date = c(), variable = c(), value = c(), estimate_type = c()))
  }
  
  ## offset dates to account for delay between infection and recorded event (testing, hospitalization, death...)
  outputDates <- outputDates - estimateOffsetting
  
  if (method == "Cori") {
    R_instantaneous <- estimate_R(
      incidenceData,
      method = "parametric_si",
      config = make_config(
        list(
          mean_si = mean_si,
          std_si = std_si,
          t_start = t_start,
          t_end = t_end,
          mean_prior = 1)
      )
    )
  } else if (method == "WallingaTeunis") {
    R_instantaneous <- wallinga_teunis(
      incidenceData,
      method = "parametric_si",
      config = list(
        mean_si = mean_si, std_si = std_si,
        t_start = t_start,
        t_end = t_end,
        n_sim = 10)
    )
  } else {
    print("Unknown estimation method")
    return(data.frame(date = c(), variable = c(), value = c(), estimate_type = c()))
  }
  
  if (variationType == "step") {
    R_mean <- unlist(lapply(seq_along(t_start),
                            function(x) {
                              rep(R_instantaneous$R$`Mean(R)`[x], t_end[x] - t_start[x] + 1)
                            }
    ))
    R_highHPD <- unlist(lapply(seq_along(t_start),
                               function(x) {
                                 rep(R_instantaneous$R$`Quantile.0.975(R)`[x], t_end[x] - t_start[x] + 1)
                               }
    ))
    R_lowHPD <- unlist(lapply(seq_along(t_start),
                              function(x) {
                                rep(R_instantaneous$R$`Quantile.0.025(R)`[x], t_end[x] - t_start[x] + 1)
                              }
    ))
  } else {
    R_mean <- R_instantaneous$R$`Mean(R)`
    R_highHPD <- R_instantaneous$R$`Quantile.0.975(R)`
    R_lowHPD <- R_instantaneous$R$`Quantile.0.025(R)`
  }
  
  if (rightTruncation > 0) {
    if (rightTruncation >= length(outputDates)) {
      return(data.frame(date = c(), variable = c(), value = c(), estimate_type = c()))
    }
    originalLength <- length(outputDates)
    outputDates <- outputDates[-seq(originalLength, by = -1, length.out = rightTruncation)]
    R_mean <- R_mean[-seq(originalLength, by = -1, length.out = rightTruncation)]
    R_highHPD <- R_highHPD[-seq(originalLength, by = -1, length.out = rightTruncation)]
    R_lowHPD <- R_lowHPD[-seq(originalLength, by = -1, length.out = rightTruncation)]
  }
  
  if (leftTruncation > 0) {
    if (leftTruncation >= length(outputDates)) {
      return(data.frame(date = c(), variable = c(), value = c(), estimate_type = c()))
    }
    originalLength <- length(outputDates)
    outputDates <- outputDates[-seq(1, leftTruncation)]
    R_mean <- R_mean[-seq(1, leftTruncation)]
    R_highHPD <- R_highHPD[-seq(1, leftTruncation)]
    R_lowHPD <- R_lowHPD[-seq(1, leftTruncation)]
  }
  
  result <- data.frame(
    date = outputDates,
    R_mean = R_mean,
    R_highHPD = R_highHPD,
    R_lowHPD = R_lowHPD)
  
  result <- reshape2::melt(result, id.vars = "date")
  colnames(result) <- c("date", "variable", "value")
  result$estimate_type <- paste0(method, "_", variationType)
  
  return(result)
}

cleanCountryReEstimate <- function(countryEstimatesRaw, method = 'bootstrap',
                                   rename_types = T,
                                   report_sd = F,
                                   alpha=0.95){
  
  if (rename_types){
    cleanEstimate <- as_tibble(countryEstimatesRaw) %>%
      mutate(
        data_type = factor(
          data_type,
          levels = c(
            "infection_Confirmed cases",
            "infection_Confirmed cases / tests",
            "infection_Hospitalized patients",
            "infection_Deaths",
            "infection_Excess deaths"),
          labels = c(
            "Confirmed cases",
            "Confirmed cases / tests",
            "Hospitalized patients",
            "Deaths",
            "Excess deaths")))
  } else {
    cleanEstimate <- as_tibble(countryEstimatesRaw)
  }
  
  if (method == 'legacy'){
    legacy_ReEstimates <- cleanEstimate %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      dplyr::group_by(date, country, region, data_type, source, estimate_type) %>%
      dplyr::summarize(
        median_R_mean = median(R_mean),
        median_R_highHPD = median(R_highHPD),
        median_R_lowHPD = median(R_lowHPD),
        mean_R_mean = mean(R_mean),
        .groups = "keep"
      ) %>%
      dplyr::select(country, region, source, data_type, estimate_type, date,
                    median_R_mean, median_R_highHPD, median_R_lowHPD,
                    mean_R_mean) %>%
      arrange(country, region, source, data_type, estimate_type, date) %>%
      ungroup()
    ReEstimates <- legacy_ReEstimates
  } else if (method == 'bootstrap') {
    
    #low_quan <- (1-alpha)/2
    high_quan <- 1-(1-alpha)/2
    
    orig_ReEstimate <- cleanEstimate %>%
      filter(replicate == 0 ) %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      rename(median_R_mean = R_mean)
    # this is called median to be compatible with legacy code
    
    MM_ReEstimates <- cleanEstimate %>%
      filter(replicate != 0 ) %>%
      pivot_wider(names_from = "variable", values_from = "value") %>%
      dplyr::group_by(date, country, region, data_type, source, estimate_type) %>%
      dplyr::summarize(
        sd_mean = sd(R_mean), #across all bootstrap replicates
        #sd_highHPD = sd(R_highHPD), #across all bootstrap replicates
        #sd_lowHPD = sd(R_lowHPD), #across all bootstrap replicates
        .groups = "drop"
      ) %>%
      right_join(orig_ReEstimate, by = c('date', 'country', 'region',
                                         'data_type', 'source', 'estimate_type')) %>%
      dplyr::mutate(median_R_highHPD = median_R_mean + qnorm(high_quan)*sd_mean,
                    median_R_lowHPD = median_R_mean - qnorm(high_quan)*sd_mean#,
                    # R_highHPD_top = R_highHPD + qnorm(high_quan)*sd_highHPD,
                    # R_highHPD_bot = R_highHPD - qnorm(high_quan)*sd_highHPD,
                    # R_lowHPD_top = R_lowHPD + qnorm(high_quan)*sd_lowHPD,
                    # R_lowHPD_bot = R_lowHPD - qnorm(high_quan)*sd_lowHPD
      ) %>%
      mutate(median_R_highHPD = ifelse(median_R_highHPD <0, 0, median_R_highHPD),
             median_R_lowHPD = ifelse(median_R_lowHPD <0, 0, median_R_lowHPD)#,
             #R_highHPD_top = ifelse(R_highHPD_top <0, 0, R_highHPD_top),
             #R_highHPD_bot = ifelse(R_highHPD_bot <0, 0, R_highHPD_bot),
             #R_lowHPD_top = ifelse(R_lowHPD_top <0, 0, R_lowHPD_top),
             #R_lowHPD_bot = ifelse(R_lowHPD_bot <0, 0, R_lowHPD_bot)
      )
    # we add the estimate-type extension later, because simpleUnion and
    # wideHPDs still derive from this df
    
    simple_Union <- MM_ReEstimates %>%
      left_join(orig_ReEstimate, by = c('date', 'country', 'region',
                                        'data_type', 'source', 'estimate_type')) %>%
      rowwise() %>%
      mutate(median_R_mean = median_R_mean.x,
             median_R_highHPD = max(median_R_highHPD, R_highHPD.y),
             median_R_lowHPD = min(median_R_lowHPD, R_lowHPD.y))#,
    #estimate_type = paste0(estimate_type, '_simple_Union'))
    
    # wideHPDs <- MM_ReEstimates %>%
    #   mutate(median_R_highHPD = R_highHPD_top,
    #          median_R_lowHPD = R_lowHPD_bot,
    #          estimate_type = paste0(estimate_type, '_wideHPDs'))
    # 
    # MM_baggedMedian <- cleanEstimate %>%
    #   filter(replicate != 0 ) %>% 
    #   pivot_wider(names_from = "variable", values_from = "value") %>%
    #   dplyr::group_by(date, country, region, data_type, source, estimate_type) %>%
    #   dplyr::summarize(
    #     sd_mean = sd(R_mean),
    #     .groups = "drop"
    #   ) %>%
    #   left_join(legacy_ReEstimates, by = c('date', 'country', 'region', 
    #                                      'data_type', 'source', 'estimate_type')) %>%
    #   dplyr::mutate(median_R_highHPD = median_R_mean + qnorm(high_quan)*sd_mean,
    #                 median_R_lowHPD = median_R_mean - qnorm(high_quan)*sd_mean) %>%
    #   mutate(median_R_highHPD = ifelse(median_R_highHPD <0, 0, median_R_highHPD),
    #          median_R_lowHPD = ifelse(median_R_lowHPD <0, 0, median_R_lowHPD)) 
    
    # MM_baggedMean <- cleanEstimate %>%
    #   #filter(replicate != 0 ) %>%
    #   pivot_wider(names_from = "variable", values_from = "value") %>%
    #   dplyr::group_by(date, country, region, data_type, source, estimate_type) %>%
    #   dplyr::summarize(
    #     sd_mean = sd(R_mean),
    #     .groups = "drop"
    #   ) %>%
    #   left_join(legacy_ReEstimates, by = c('date', 'country', 'region',
    #                                        'data_type', 'source', 'estimate_type')) %>%
    #   dplyr::mutate(median_R_highHPD = mean_R_mean + qnorm(high_quan)*sd_mean,
    #                 median_R_lowHPD = mean_R_mean - qnorm(high_quan)*sd_mean) %>%
    #   mutate(median_R_highHPD = ifelse(median_R_highHPD <0, 0, median_R_highHPD),
    #          median_R_lowHPD = ifelse(median_R_lowHPD <0, 0, median_R_lowHPD))
    
    # bag_Union <- MM_baggedMean %>%
    #   left_join(legacy_ReEstimates, by = c('date', 'country', 'region', 
    #                                        'data_type', 'source', 'estimate_type')) %>%
    #   rowwise() %>%
    #   mutate(median_R_mean = median_R_mean.x,
    #          median_R_highHPD = max(median_R_highHPD.x, median_R_highHPD.y),
    #          median_R_lowHPD = min(median_R_lowHPD.x, median_R_lowHPD.y),
    #          estimate_type = paste0(estimate_type, '_bag_Union'))
    
    unsortedReEstimates <- bind_rows(#legacy_ReEstimates, 
      #MM_ReEstimates #%>% mutate(estimate_type = paste0(estimate_type, '_MM')),
      #MM_baggedMedian %>% mutate(estimate_type = paste0(estimate_type, '_MM_baggedMedian')),
      #MM_baggedMean #%>% mutate(estimate_type = paste0(estimate_type, '_MM_baggedMean'))
      simple_Union #, bag_Union,
      #wideHPDs
    ) 
    
    if (report_sd){
      ReEstimates <- unsortedReEstimates %>%
        dplyr::select(country, region, source, data_type, estimate_type, date,
                      median_R_mean, median_R_highHPD, median_R_lowHPD, sd_mean) %>%
        arrange(country, region, source, data_type, estimate_type, date) %>%
        ungroup()
    } else {
      ReEstimates <- unsortedReEstimates %>%
        dplyr::select(country, region, source, data_type, estimate_type, date,
                      median_R_mean, median_R_highHPD, median_R_lowHPD) %>%
        arrange(country, region, source, data_type, estimate_type, date) %>%
        ungroup()
    }
    
  }
  return(ReEstimates)
}

# # # # # # # # FUNCTION FOR Rt METHODS COMPARISON PAPER # # # # # # #
rt_huismann_function <- function( dataframe, region, weekly){
  
  colnames(dataframe) <- c("date", "predictor")
  
  # normalize data by the minimum detected value in the county
  dataframe$norm_sars <- dataframe$predictor / min(dataframe$predictor, na.rm = TRUE)
  
  # match column names in huismann method
  dataframe$region <- region

  # create config df: note this can be used to modify the serial intervals based on variant type. we do not have that, so using Huismann defaults
  config_df = expand.grid("region" = region, # for region we can replace with SW_ID  
                          'incidence_var' = c('norm_sars', 'norm_sars'),
                          'FirstGamma' = 'incubation',
                          'SecondGamma' = 'benefield' )
  
  
  # objects to store output
  deconv_ww_data <- data.frame()
  Re_ww <- data.frame()
  
  for(row_i in 1:nrow(config_df)){
    
    # huismann functions
    new_deconv_data = deconvolveIncidence(dataframe %>% filter(region == config_df[row_i, 'region']), 
                                          incidence_var = config_df[row_i, 'incidence_var'],
                                          getCountParams(as.character(config_df[row_i, 'FirstGamma'])), 
                                          getCountParams(as.character(config_df[row_i, 'SecondGamma'])),
                                          smooth_param = TRUE, n_boot = 50)
    
    new_deconv_data <- new_deconv_data %>%
      mutate(incidence_var = config_df[row_i, 'incidence_var'])
    
    ##### Get Re #####
    new_Re_ww = getReBootstrap(new_deconv_data)
    new_Re_ww <- new_Re_ww %>%
      mutate(variable = config_df[row_i, 'incidence_var'],
             region = config_df[row_i, 'region'])
    
    deconv_ww_data <- bind_rows(deconv_ww_data, new_deconv_data)
    Re_ww = bind_rows(Re_ww, new_Re_ww)
  }
  
  
  # remove unnecessary columns
  Re_ww <- Re_ww %>%
    filter(date >= min(dataframe$date))%>%
    filter(!duplicated(date))%>%
    select(-country, -source, - data_type, -variable, -region, -estimate_type) %>%
    mutate(se_rt = (median_R_highHPD - median_R_lowHPD) / 3.02
           ) %>%
    rename(mean_rt = median_R_mean,
           ll_95_rt = median_R_lowHPD,
           ul_95_rt = median_R_highHPD)

  
  if(weekly == "Yes"){
    
    Re_ww_weekly <- Re_ww %>%
      group_by(week = floor_date(date, unit = "weeks")) %>%
      summarize(
        mean_rt = mean(mean_rt, na.rm = TRUE),
        se_rt = mean(se_rt, na.rm = TRUE),
        ll_95_rt = mean(ll_95_rt, na.rm = TRUE),
        ul_95_rt = mean(ul_95_rt, na.rm = TRUE)
      ) %>%
      ungroup()
    
    return(Re_ww_weekly)
    
  } else if(weekly == "No"){
    return(Re_ww)
  }
  
}