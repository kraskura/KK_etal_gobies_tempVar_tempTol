
# function for extracting a specific timeframe within a dataset 

dataMutateMeans<-function(data){
  
  data$minutes.dayROUND<-as.factor(round(data$minutes.day, digits = -1))
  
  data<-data %>% 
    group_by(minutes.dayROUND, SiteID) %>% 
    mutate(mean_tempROUND = mean(TEMP)) %>% 
    dplyr:::group_by(mo, day, y, SiteID) %>%
    mutate(mean_temp = mean(TEMP),
           min_temp = min(TEMP),
           max_temp = max(TEMP),
           var_temp = var(TEMP),
           sd_temp = sd(TEMP),
           range_temp = max_temp - min_temp) %>%
    dplyr:::group_by(SiteID) %>%
    mutate(meanMEAN = round(mean(mean_temp),1),
           meanMAX = round(mean(max_temp),1), 
           meanMIN = round(mean(min_temp),1),) %>%
    as.data.frame()
  
  return(data)
  
}
  
