# Function for consistent data formatting 
transformTempData<-function(data, MeanGroupVar){
  
  # get various time frames separated out. 
  data <- transform(data,time.str = format(LocalDateTime,'%H:%M:%S'),
                    y = year(LocalDateTime),
                    mo = month(LocalDateTime),
                    day = day(LocalDateTime),
                    hour = hour(LocalDateTime),
                    minute = minute(LocalDateTime),
                    second= second(LocalDateTime))
  data$minutes.day <- 60 * data$hour + data$minute
  
  # calculate rolling average
  data <- data %>%
    dplyr::group_by(!!MeanGroupVar) %>% 
    dplyr::mutate(temp_0.5h = zoo::rollmean(TEMP, k = 3, fill = NA), # mean of 3 vals 
                  temp_1h = zoo::rollmean(TEMP, k = 6, fill = NA),# mean of 6 vals 
                  temp_2h = zoo::rollmean(TEMP, k = 12, fill = NA)) %>% # mean of 12 vals 
    dplyr::ungroup() %>% 
    as.data.frame()
  
  return(data)
}
