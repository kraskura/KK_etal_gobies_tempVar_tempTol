# function for extracting a specific timeframe within a dataset 
extractTimeFrame<-function(data, year, minMo, maxMo, minDay, maxDay){

  data.trunc<-data[c(
    c(grepl(year, data$y)  & 
        c(data$mo >= minMo & data$mo <= maxMo) &   
        c(data$day >= minDay & data$day <= maxDay))), ]
  
  return(data.trunc)
  
}
