
library(devtools)
install_github("kraskura/ggformat2")
library(here)
library(tidyverse)
library(chron)
library(lubridate)
# plotting 
library(gridExtra)
library(ggformat2) # kraskura/ggformat2
library(cowplot)

# Import functions that ensure consistent formatting ----------
source("./R/extractTimeFrame.R") # extract wanted specific time frames for plotting, or any other purpose etc. 
source("./R/transformTempData.R") # add various timeframe variables, h, min, day.min, rolling avg for TEMP ºC. 
source("./R/dataMutateMeans.R") # Add column to sumarie means, max, daily range, etc. of TEMP ºC

dataF<-read.csv("./Data/Kraskura_combined_Study_data2023.csv", header = T)
names(dataF)<-c("DateTime", "TEMP", "SiteID", "Latitude","Longitude")
dataF$DateTime<-as.POSIXct(mdy_hm(dataF$DateTime), tz = "America/Los_Angeles") # in configuartion settings this was accurately set up as PDT. not UTC timezone. - confirmed 
dataF$LocalDateTime<-dataF$DateTime # timezone: "America/Los_Angeles" or Pacific Day time "PDT"
dataF<-dataF[complete.cases(dataF),]
dataF <- transformTempData(dataF, MeanGroupVar = "SiteID")
dataF <- dataMutateMeans(dataF)
yPlot<-2019
moPlot<-7
dPlot1<-7
dPlot2<-17
data<-dataF
  
if(!dir.exists("Figures/Field temps sites")){
  dir.create("Figures/Field temps sites", recursive = T)
}

plotFieldData<-function(data, yPlot, moPlot, dPlot1, dPlot2, printplot = FALSE){

  data.plot <- extractTimeFrame(data, year = yPlot, minMo = moPlot, maxMo = moPlot, minDay = dPlot1, maxDay = dPlot2)

  if(nrow(data.plot) == 0) {
    stop ("Not data")
  }

  plotname <- paste("Figures/Field temps sites/FieldTemp_y", yPlot, "mo", moPlot,"_dStart_", dPlot1, ".png", sep ="")
  siteColors <- c("SITE2" = "#A149FA",
                  "SITE3" = "#3B44F6", "SITE1" = "black",
                  "SITE4" = "#3EC70B",
                  "TestSite" = "#00CAFF")
# "#F7EC09"

  plot.site<-ggplot(data.plot)+
      geom_line(aes(DateTime, TEMP, group = SiteID, color = SiteID), linewidth = 0.7, alpha = 0.9)+
      scale_y_continuous(limits = c(0, 40), breaks = seq(5, 40, 5))+
      scale_color_manual(values = siteColors)+
      geom_text(mapping = aes(x = as.POSIXct(data.plot$DateTime[nrow(data.plot)]) + 60*60*3,
                              y = round(mean(max_temp),2),
                              group = SiteID, label=paste(round(mean(max_temp),2), "ºC", sep ="")),
                              hjust = 0, size = 3, check_overlap = TRUE)+
      geom_text(mapping = aes(x = as.POSIXct(data.plot$DateTime[nrow(data.plot)]) + 60*60*3,
                              y = round(mean(mean_temp),2),
                              group = SiteID, label=paste(round(mean(mean_temp),2), "ºC", sep ="")),
                              hjust=0, size = 3, check_overlap = TRUE)+
      geom_text(mapping = aes(x = as.POSIXct(data.plot$DateTime[nrow(data.plot)]) + 60*60*3,
                              y = round(mean(min_temp),2),
                              group = SiteID, label=paste(round(mean(min_temp),2), "ºC", sep ="")),
                              hjust=0,  size = 3, check_overlap = TRUE)+
      scale_x_datetime(limits = c(as.POSIXct(data.plot$DateTime[1]),
                                  as.POSIXct(data.plot$DateTime[nrow(data.plot)]) + 60*60*24))
  ggformat(plot.site, title = "",
           y_title = expression(Temperature~(degree*C)),
           x_title = "",
           size_text = 10, print = printplot)
  plot.site <- plot.site + theme(axis.title.x = element_blank(), legend.position = "none")

  ggsave(filename = plotname, plot.site, width = 4, height = 2)

}

years<-c(2019, 2020, 2021, 2022)
months<-seq(1,12, 1)

for(i in 1:length(years)){
  for (j in 1:length(months)){
    print(c(i,j))
    data.plot <- extractTimeFrame(data, year = years[i], minMo = months[j], maxMo = months[j], minDay = 1, maxDay = 7)
    if(nrow(data.plot) == 0) {
      next
    }
    "here- yey!"
    plotFieldData(data = dataF, yPlot = years[i] , moPlot = months[j], dPlot1 = 1, dPlot2 = 7)
    
  }
}




 