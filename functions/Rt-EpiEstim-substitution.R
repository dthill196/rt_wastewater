# use rt equation -> add start date
rt_ww_sub_function <- function(dataframe, mean_si, std_si,  weekly, unknown_si){
  
  if(unknown_si == "Yes"){
    rt_ww_sub <- rt_function_unknown_si(dataframe,  weekly)
    
  } else if(unknown_si == "No"){
    rt_ww_sub <- rt_function(dataframe, mean_si, std_si, weekly)
    
  }
  
}
