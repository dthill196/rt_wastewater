# Rt from wastewater using rolling GLM of cases ~ wastewater

rt_glm_function <- function(dataframe, predictor,case_data, range, unknown_si, mean_si, std_si, weekly){
  
  # unknown si function for case Rt
  if(unknown_si == "Yes"){
    
    # weekly Rt data for unknown si
    if(weekly == "Yes"){
      # for loop
      range= range
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
        
        # round to integer
        df$cases_2 <- round(df[[case_data]], 0)
        
        # predict in glm, poisson distribution
        model <- glm(cases_2 ~ df[[predictor]],
                     data = df,
                     family = "poisson")
        
        # calculate cases
        df$ww_pred_cases <- predict(model, newdata = df, type = "response")
        df <- df %>%
          filter(date == i)
        datalist1[[i]] <- df
      }
      
      final <-do.call(rbind, datalist1)  
      
      # prep dataframe
      glm_data <- final %>%
        select(date, ww_pred_cases) %>%
        rename(Date = date) %>%
        filter(!is.na(ww_pred_cases))
      
      # input predicted case data into the rt function with known si interval
      rt_glm_weekly <- rt_function_unknown_si(dataframe = glm_data, weekly)
      
      return(rt_glm_weekly)
      
    # daily Rt data for unknown si
    } else if(weekly == "No"){
      # for loop
      range = range
      
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
        
        # round to integer
        df$cases_2 <- round(df[[case_data]], 0)
        
        # predict in glm, poisson distribution
        model <- glm(cases_2 ~ df[[predictor]],
                     data = df,
                     family = "poisson")
        
        # calculate cases
        df$ww_pred_cases <- predict(model, newdata = df, type = "response")
        df <- df %>%
          filter(date == i)
        datalist1[[i]] <- df
      }
      
      final <-do.call(rbind, datalist1)  
      
      # prep dataframe
      glm_data <- final %>%
        select(date, ww_pred_cases) %>%
        rename(Date = date) %>%
        filter(!is.na(ww_pred_cases))
      
      # input predicted case data into the rt function with known si interval
      rt_glm <- rt_function_unknown_si(dataframe = glm_data, weekly)
      
      return(rt_glm)
    }
    
    # known si, use this one
  } else if(unknown_si == "No"){
    
    # weekly data known si
    if(weekly == "Yes"){
      # for loop
      range= range
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
        
        # round to integer
        df$cases_2 <- round(df[[case_data]], 0)
        
        # predict in glm, poisson distribution
        model <- glm(cases_2 ~ df[[predictor]],
                     data = df,
                     family = "poisson")
        
        # calculate cases
        df$ww_pred_cases <- predict(model, newdata = df, type = "response")
        df <- df %>%
          filter(date == i)
        datalist1[[i]] <- df
      }
      
      final <-do.call(rbind, datalist1)  
      
      # prep dataframe
      glm_data <- final %>%
        select(date, ww_pred_cases) %>%
        rename(Date = date) %>%
        filter(!is.na(ww_pred_cases))
      
      # input predicted case data into the rt function with known si interval
      rt_glm_weekly <- rt_function(dataframe = glm_data, mean_si, std_si, weekly)
      
      return(rt_glm_weekly)
      
    # daily data known si
    } else if(weekly == "No"){
      # for loop
      range = range
      
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
        
        # round to integer
        df$cases_2 <- round(df[[case_data]], 0)
        
        # predict in glm, poisson distribution
        model <- glm(cases_2 ~ df[[predictor]],
                     data = df,
                     family = "poisson")
        
        # calculate cases
        df$ww_pred_cases <- predict(model, newdata = df, type = "response")
        df <- df %>%
          filter(date == i)
        datalist1[[i]] <- df
      }
      
      final <-do.call(rbind, datalist1)  
      
      # prep dataframe
      glm_data <- final %>%
        select(date, ww_pred_cases) %>%
        rename(Date = date) %>%
        filter(!is.na(ww_pred_cases))
      
      # input predicted case data into the rt function with known si interval
      rt_glm <- rt_function(dataframe = glm_data, mean_si, std_si, weekly)
      
      return(rt_glm)
    }
  }
}

