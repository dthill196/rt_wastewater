rt_evaluation_function <- function(method_name,
                                   reference_dataframe,
                                   method_dataframe){
  
  # add week
  reference_dataframe$date <- reference_dataframe$week
  method_dataframe$date <- method_dataframe$week
  
  # edit names for reference group
  colnames(reference_dataframe) <- paste("case", colnames(reference_dataframe), sep = "_")
  dataframe <- reference_dataframe %>%
    rename(date = case_date) %>%
    left_join(method_dataframe, by = c("date"))
  
  # # rmse
  rmse_1 <- sqrt(mean((dataframe$case_mean_rt - dataframe$mean_rt)^2, na.rm = TRUE))
  
  # # # # # # # # # #
  # Metric 2 - Correlation (Pearson)
  # # # # # # # # # #
  
  cor_1 <- cor(dataframe$case_mean_rt, dataframe$mean_rt, method = "pearson", use = "na.or.complete")
  
  # # # # # # # # # #
  # Metric 3 - do timepoint for peaks coincide
  # # # # # # # # # #
  
  # function for peaks
  peaks_function <- function(dataframe,
                             reference_rt,
                             comparison_rt){
    
    # find peaks
    case_peaks <- findPeaks(reference_rt)
    fit_peaks <- findPeaks(comparison_rt)
    
    # identify rows that are peaks
    dataframe$case_peak <- ifelse(rownames(dataframe) %in% case_peaks, 1, 0)
    dataframe$fit_peak <- ifelse(rownames(dataframe) %in% fit_peaks, 1, 0)
    
    # group by week and see if they correspond
    rt_week <- dataframe %>%
      filter(case_peak >0 | fit_peak >0) %>%
      mutate(week = floor_date(date, unit = "weeks")) %>%
      group_by(week)%>%
      summarize(case_peak_sum = sum(case_peak),
                fit_peak_sum = sum(fit_peak)
      )%>%
      ungroup()
    
    rt_week$peaks_agree <- ifelse(rt_week$case_peak_sum > 0 & rt_week$fit_peak_sum >0, 1, 0)
    
    # sum and divide by the length of the df
    peaks <- sum(rt_week$peaks_agree) / nrow(rt_week)
    return(peaks)
  }
  
  peaks_fit <- peaks_function(dataframe = dataframe,
                              reference_rt = dataframe$case_mean_rt,
                              comparison_rt = dataframe$mean_rt)
  
  # # # # # # # # # #
  # Metric 4 - is the value above or below 1
  # # # # # # # # # #
  
  dataframe$case_above_1 <- ifelse(dataframe$case_mean_rt > 1, 1, 0)
  dataframe$fit_above_1 <- ifelse(dataframe$mean_rt >1, 1, 0)
  
  dataframe$above_1_agree <- ifelse(dataframe$case_above_1 == dataframe$fit_above_1, 1, 0)
  
  # percent agreement
  above_1 <- sum(dataframe$above_1_agree, na.rm = TRUE) / nrow(dataframe)
  
  # # # # # # # # # #
  # Metric 5 - Absolute difference between methods (or just mean difference)
  # # # # # # # # # #
  
  mean_abs_dif <- mean(abs(dataframe$case_mean_rt - dataframe$mean_rt), na.rm = TRUE)
  
  # # # # # # # # # #
  # Metric 6 - Sharpness of confidence intervals
  # # # # # # # # # #
  
  # ci width function (used inside sharpness function)
  ci_width_function <- function(dataframe, critical_value, mean_r, se_r){
    dataframe$ll <- mean_r - critical_value * se_r
    dataframe$ul <- mean_r + critical_value * se_r
    dataframe$ci_width <- dataframe$ul - dataframe$ll
  }
  
  # function for sharpness
  sharpness_function <- function(dataframe, mean_r, se_r){
    
    critical_value_list <- as.character(c(1.04, 1.15, 1.28, 1.44, 1.645, 1.75, 1.96, 2.05, 2.33, 2.58))
    
    datalist <- list()
    for(i in unique(critical_value_list)){
      
      # use function to calcuate widths for cis
      width_fit <- ci_width_function(
        dataframe = dataframe,
        critical_value = as.numeric(i),
        mean_r = mean_r,
        se_r = se_r 
      )
      
      # store in df as a vector
      df <- as.data.frame(width_fit)
      colnames(df) <- paste("width", i, sep = "-")
      
      # combine widths
      datalist[[i]] <- df
      
    }
    
    # combine into df
    fit_widths <- do.call(cbind, datalist)
    
    # now calculate the sharpness
    # weighted mean = sum(value * weight)/sum(weights)
    fit_widths$w_mean_sharpness <- (fit_widths$`width-1.04`*(1*(1-0.7)) +
                                      fit_widths$`width-1.15`*(1*(1-0.75))+
                                      fit_widths$`width-1.28`*(1*(1-0.8)) +
                                      fit_widths$`width-1.44`*(1*(1-0.85)) +
                                      fit_widths$`width-1.645`*(1*(1-0.9)) +
                                      fit_widths$`width-1.75`*(1*(1-0.92)) +
                                      fit_widths$`width-1.96`*(1*(1-0.95)) +
                                      fit_widths$`width-2.05`*(1*(1-0.96)) +
                                      fit_widths$`width-2.33`*(1*(1-0.98)) +
                                      fit_widths$`width-2.58`*(1*(1-0.99)) 
                                    
    )/ ((1*(1-0.7)) + 
          (1*(1-0.75)) + 
          (1*(1-0.8)) + 
          (1*(1-0.85)) + 
          (1*(1-0.9)) + 
          (1*(1-0.92)) + 
          (1*(1-0.95)) + 
          (1*(1-0.96)) + 
          (1*(1-0.98)) +
          (1*(1-0.99))
    )
    
    # weighted mean
    sharpness <- mean(fit_widths$w_mean_sharpness, na.rm = TRUE)
    
    return(sharpness)
  }
  
  sharpness_fit <- sharpness_function(dataframe = dataframe,
                                      mean_r = dataframe$mean_rt,
                                      se_r = dataframe$se_rt)
  
  
  # one function for all, return a table with all values and the method name
  metrics <- c(rmse_1, cor_1, peaks_fit, above_1, mean_abs_dif, sharpness_fit)
  metric_names <- c("RMSE", "Pearson Correlation", "Percent of peaks that coincide", "Above or below 1 agreement", "Mean absolute difference",
                    "Sharpness")
  
  eval_table <- as.data.frame(cbind( metric_names, metrics)) 
  colnames(eval_table)[2] <- method_name
  
  return(eval_table)
}
