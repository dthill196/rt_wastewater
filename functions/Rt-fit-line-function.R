# dataframe <- data_fit_line
# weekly <- "Yes"
# range <- 45
# predictor <- "lag_sars2.7avg"
rt_fit_line_function <- function(dataframe, weekly, range, predictor){

  if(weekly == "Yes"){
    
    # for loop
    range= range
    # create date range (skip first 2 weeks of data)
    date_list <- dataframe %>%
      arrange(date)%>%
      slice(8:nrow(.)) %>%
      select(date)
    
    datalist1 <- list()
    
    ww <- dataframe
    
    for(i in unique(date_list$date)){
      
      # create window for dates
      df <- ww %>%
        filter(date >= as.Date(i)-range) %>%
        filter(date <= as.Date(i))
      
      lm2 <- df %>%
        do(tidy(lm(rt ~ df[[predictor]], .)))
      
      
      df$actual_slope <- lm2$estimate[2]
      df$slope_se <- lm2$std.error[2]
      df$intercept <- lm2$estimate[1]
      df$intercept_se <- lm2$std.error[1]
      df <- df %>%
        slice(which.max(date))
      
      datalist1[[i]] <- df
    }
    
    final <-do.call(rbind, datalist1)  
    
    # calculate rt from fit and calculate
    # convert all the values from dat_fit_line copies to this scale
    final$rt_fit_line_roll_reg <- final$intercept + (final[[predictor]] * final$actual_slope)
    
    # calculate lower bound for rt
    #meanHeight - 1.96 * seHeight
    final$slope_ll <- final$actual_slope - 1.96 * final$slope_se
    final$intercept_ll <- final$intercept - 1.96 * final$intercept_se
    final$rt_fit_line_ll <- final$intercept_ll + (final[[predictor]] * final$slope_ll)
    
    # calculate upper bound for rt
    final$slope_ul <- final$actual_slope + 1.96 * final$slope_se
    final$intercept_ul <- final$intercept + 1.96 * final$intercept_se
    final$rt_fit_line_ul <- final$intercept_ul + (final[[predictor]] * final$slope_ul)
    
    # weekly means for the final data
    final_weekly <- final %>%
      group_by(week = floor_date(date, unit = "week")) %>%
      summarize(rt_fit_line_roll_reg_weekly = mean(rt_fit_line_roll_reg, na.rm = TRUE),
             rt_fit_line_ll_weekly = mean(rt_fit_line_ll, na.rm = TRUE),
             rt_fit_line_ul_weekly = mean(rt_fit_line_ul, na.rm = TRUE)
             ) %>%
      ungroup() %>%
      mutate(se_rt = (rt_fit_line_ul_weekly - rt_fit_line_ll_weekly / 3.92) ) %>%
      rename(mean_rt = rt_fit_line_roll_reg_weekly,
             ll_95_rt = rt_fit_line_ll_weekly,
             ul_95_rt = rt_fit_line_ul_weekly) %>%
      select(week, mean_rt, se_rt, ll_95_rt, ul_95_rt)
    
    return(final_weekly)
    
  } else if(weekly == "No"){
    
    # for loop
    range= 45
    # create date range (skip first week of data)
    date_list <- dataframe %>%
      arrange(date)%>%
      slice(8:nrow(.)) %>%
      select(date)
    
    datalist1 <- list()
    
    ww <- dataframe
    
    for(i in unique(date_list$date)){
      
      # create window for dates
      df <- ww %>%
        filter(date >= as.Date(i)-range) %>%
        filter(date <= as.Date(i))
      
      lm2 <- df %>%
        do(tidy(lm(rt ~ df[[predictor]], .)))
      
      
      df$actual_slope <- lm2$estimate[2]
      df$slope_se <- lm2$std.error[2]
      df$intercept <- lm2$estimate[1]
      df$intercept_se <- lm2$std.error[1]
      df <- df %>%
        slice(which.max(date))
      
      datalist1[[i]] <- df
    }
    
    final <-do.call(rbind, datalist1)  
    
    # calculate rt from fit and calculate
    # convert all the values from dat_fit_line copies to this scale
    final$rt_fit_line_roll_reg <- final$intercept + (final[[predictor]] * final$actual_slope)
    
    # calculate lower bound for rt
    #meanHeight - 1.96 * seHeight
    final$slope_ll <- final$actual_slope - 1.96 * final$slope_se
    final$intercept_ll <- final$intercept - 1.96 * final$intercept_se
    final$rt_fit_line_ll <- final$intercept_ll + (final[[predictor]] * final$slope_ll)
    
    # calculate upper bound for rt
    final$slope_ul <- final$actual_slope + 1.96 * final$slope_se
    final$intercept_ul <- final$intercept + 1.96 * final$intercept_se
    final$rt_fit_line_ul <- final$intercept_ul + (final[[predictor]] * final$slope_ul)
    
    final <- final %>%
      mutate(se_rt = (rt_fit_line_ul -rt_fit_line_ll / 3.92) ) %>%
      rename(mean_rt = rt_fit_line_roll_reg) %>%
      select(date, mean_rt, se_rt, ll_95_rt, ul_95_rt)
    
    
    return(final)
    
  }
}
