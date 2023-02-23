
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

# Import functions that help with consistent formatting ----------
source("./R/extractTimeFrame.R") # extract wanted specific time frames for plotting, or any other purpose etc. 
source("./R/transformTempData.R") # add various timeframe variables, h, min, day.min, rolling avg for TEMP ºC. 
source("./R/dataMutateMeans.R") # Add column to sumarie means, max, daily range, etc. of TEMP ºC

if(!dir.exists("Data/Analysis source")){
  dir.create("Data/Analysis source", recursive = T)
}


# LAB treatment temps----------
labTfiles<-paste("./Data/Lab_temperatures/", list.files(path = "./Data/Lab_temperatures/"), sep = "")
data <- do.call(rbind,lapply(labTfiles, read.csv))
data <- data[data$Value < 100 & data$Value > 0, ]# filter out technical error temps

# assign each temperature probe to its appropriate treatment 
data$treatm <- ifelse(data$Label == "Probe1", '27', # ºC treatment 
                      ifelse(data$Label == "Probe2", 'V1', # consistent variation treatment 
                             ifelse(data$Label == "Probe3", '12',
                                    ifelse(data$Label == "Probe4", 'V1',
                                           ifelse(data$Label == "Probe11", '22',
                                                  ifelse(data$Label == "Probe12", '17',
                                                         ifelse(data$Label == "Probe13", 'V2', # inconsitent, spikey variation treatment 
                                                                ifelse(data$Label == "Probe14", '12',
                                                                       ifelse(data$Label == "Probe15", '22',
                                                                              ifelse(data$Label == "Probe5", '27',
                                                                                     ifelse(data$Label == "Probe6", '17',
                                                                                            ifelse(data$Label == "Probe7", 'V2',
                                                                                                   ifelse(data$Label == "Probe8", '12',
                                                                                                          ifelse(data$Label == "Probe16", 'V1',
                                                                                                                 ifelse(data$Label == "Probe17", 'V2',
                                                                                                                        ifelse(data$Label == "Probe18", '27',
                                                                                                                               ifelse(data$Label == "Probe19", '22',
                                                                                                                                      ifelse(data$Label == "Probe20", '17','AIR'))))))))))))))))))

# assign each temperature probe to its appropriate tank ID in the lab 
data$tank <- ifelse(data$Label == "Probe1", 1,
                    ifelse(data$Label == "Probe2", 2,
                           ifelse(data$Label == "Probe3", 3,
                                  ifelse(data$Label == "Probe4", 4,
                                         ifelse(data$Label == "Probe11", 5,
                                                ifelse(data$Label == "Probe12", 6,
                                                       ifelse(data$Label == "Probe13", 7,
                                                              ifelse(data$Label == "Probe14", 8,
                                                                     ifelse(data$Label == "Probe15", 9,
                                                                            ifelse(data$Label == "Probe5", 10,
                                                                                   ifelse(data$Label == "Probe6", 11,
                                                                                          ifelse(data$Label == "Probe7", 12,
                                                                                                 ifelse(data$Label == "Probe8", 13,
                                                                                                        ifelse(data$Label == "Probe16", 14,
                                                                                                               ifelse(data$Label == "Probe17", 15,
                                                                                                                      ifelse(data$Label == "Probe18", 16,
                                                                                                                             ifelse(data$Label == "Probe19", 17, 
                                                                                                                                    ifelse(data$Label == "Probe20", 18, 0))))))))))))))))))


# check data for the quality
data[which(is.na(data$Value)),] # only AIR temp, take out 
data<-data[-c(which(is.na(data$Value))),] 
data<-data[, c(2, 3, 4, 6, 7)]
names(data)<-c("Label", "TEMP", "DateTime", "treatm", "tank" )
head(data)

# data$Date = as.Date(temps.mean$ServerTime, format=c('%m-%d-%y %H:%M')) 
data$LocalDateTime<-as.POSIXct(data$DateTime, format="%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles") # server time is eastern time "America/New_York"
# class(temps.mean$LocalDateTime)

# ********************************
dataLab<-transformTempData(data, MeanGroupVar = "Label")

# these tanks did not have fish in them after indicated dates in November, exclude the temperatures 
# nrow(dataLab)
dataLab <- dataLab[-c(which((dataLab$Label == "Probe18" & dataLab$mo == 11 & dataLab$day >= 22) |
                                 (dataLab$Label == "Probe1" & dataLab$mo == 11 & dataLab$day >= 22) |
                                 (dataLab$Label == "Probe5" & dataLab$mo == 11 & dataLab$day >= 15) |
                                 (dataLab$Label == "Probe14" & dataLab$mo == 11 & dataLab$day >= 21))),]
dataLab <- dataLab[-which(dataLab$Label == "Probe16" & dataLab$mo == 11 & dataLab$day >= 20),]

# take out air temp:
dataLab<-dataLab[-c(which(dataLab$treatm == "AIR")),]
dataLab$SiteID<-"Lab"
  
# separate out static treatments & applying custom function to get some running means and daily temp stats (min, max, mean.. )
dataLab.S<-dataMutateMeans(dataLab[!c(dataLab$treatm=="V1" | dataLab$treatm=="V2"), ])
dataLab.all<-dataMutateMeans(dataLab)
# separate out variability treatments 
dataLab.V1<-dataMutateMeans(dataLab[c(dataLab$treatm=="V1"), ]) 
dataLab.V2<-dataMutateMeans(dataLab[c(dataLab$treatm=="V2"), ]) # stochastic variable treatment
dataLab.V1V2<-dataMutateMeans(dataLab[c(dataLab$treatm=="V1" | dataLab$treatm=="V2"), ]) # both 

# dataset for plotting example Oct 19 - Oct 26 (representative trends of lab treatments)
dataLab.V1.plot <- extractTimeFrame(data = dataLab.V1, year = 2019, minMo = 10, maxMo = 10, minDay = 19, maxDay = 30)  # regular variable treatment 
dataLab.V2.plot <- extractTimeFrame(data = dataLab.V2, year = 2019, minMo = 10, maxMo = 10, minDay = 19, maxDay = 30)  # stochastic variable treatment 
dataLab.S.plot <- extractTimeFrame(data = dataLab.S, year = 2019, minMo = 10, maxMo = 10, minDay = 19, maxDay = 30) 
data.plot<-extractTimeFrame(data = dataLab.all, year = 2019, minMo = 10, maxMo = 10, minDay = 19, maxDay = 25) 

# datesets of test days. 
dataLab.test.days <- extractTimeFrame(data = dataLab.all, year = 2019, minMo = 11, maxMo = 11, minDay = 11, maxDay = 31)  # variable and static treatment 
dataLab.test.daysVar <- dataLab.test.days[c(dataLab.test.days$treatm=="V1" | dataLab.test.days$treatm=="V2"), ]
# clean up tech error temps. 
dataLab.test.days<-dataLab.test.days[!c(dataLab.test.days$treatm == "27" & dataLab.test.days$TEMP<26.5), ]
write_csv(dataLab.test.days, file = "./Data/Analysis source/Lab_TestDaytemps_2019_allTreatm.csv")
write_csv(dataLab.test.daysVar, file = "./Data/Analysis source/Lab_TestDaytemps_2019_varTreatm.csv")

# Manuscript Fig 2: figure all treatment together -----
lab.SECTION<-ggplot(data.plot[!c(data.plot$treatm == "V2"), ], aes(LocalDateTime, TEMP, color = treatm, fill = treatm,
                                   group = interaction(mo, day, treatm, tank)))+
  scale_y_continuous(limits = c(3, 35), breaks = seq(5, 35, 5))+
  guides(color=guide_legend(nrow=1,byrow=TRUE))+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD"), 
                     name = "", labels = c("12ºC", "17ºC", "22ºC", "27ºC", "V1", "V2"))+
  scale_fill_manual(values=c("12" ="#00518C",  "17" = "#008A60","22" ="#DBA11C","27" ="#BE647D", "V2" = "black", "V1" = "#A3ABBD"),
                    name = "", labels = c("12ºC", "17ºC", "22ºC", "27ºC", "V1", "V2"))+
  geom_line(data = data.plot[c(data.plot$treatm == "V2"), ], aes(LocalDateTime, TEMP, color = treatm,
                                   group = interaction(mo, day, treatm, tank)), linetype = 3, linewidth = 0.6)+
  geom_line( alpha=0.6,  linewidth=1)
ggformat(lab.SECTION, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time", print=F, size_text = 12)
lab.SECTION<- lab.SECTION+theme(axis.text.x = element_text( size=11),
                                axis.title.x = element_blank(),
                                legend.position = "none")
                               # legend.position = c(0.5, 0.96),
                               # legend.background = element_rect(colour = NA, fill = NA), 
                               # legend.key.width = unit(1, "line"))

lab.SECTION


## Figure i) SUPPL : daily var all tanks, all treatm ---------------
daily.p.ALL <- ggplot(dataLab, aes(minutes.day, TEMP, group = interaction(treatm), color = treatm))+
  # geom_point(size=1, pch=".")+
  ylim(5, 35)+
  scale_color_manual(values=c("#00518C", "#008A60", "#DBA11C", "#BE647D","#3F4756", "#A3ABBD") )+
  annotate("rect", xmin=0, xmax=360, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  annotate("rect", xmin=960, xmax=1440, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  geom_line(aes(minutes.day, TEMP, group = interaction(mo, day)), linewidth=0.4, alpha=0.5)+
  facet_wrap(nrow = 7, ncol=3,
             ~factor(tank))
                     # , levels=c('1','10','16','5', '9', '17', '6', '11', '18', '3', '8', '13', '2', '4', '14', '7', '12', '15')))
ggformat(daily.p.ALL, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (min)")
daily.p.ALL <- daily.p.ALL + theme(legend.position = "top")


## Figure ii) SUPPL: static treatments ---------------
daily.p <- ggplot(dataLab.S.plot, aes(minutes.day/60, TEMP, group = interaction(treatm), color = treatm))+
  ylim(5, 35)+
  scale_color_manual(values=c("#00518C", "#008A60", "#DBA11C", "#BE647D") )+
  annotate("rect", xmin=0, xmax=360/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  annotate("rect", xmin=960/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  geom_line(aes(minutes.day/60, TEMP, group = interaction(mo, day, treatm, tank)), linewidth=0.4, alpha=0.5)
ggformat(daily.p, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 15, print = F)
daily.p <- daily.p + theme(legend.position = "none")

## Figure ii) SUPPL: V2, stochastic variable treatment; stacked trends ---------------
daily.p.V2 <- ggplot(dataLab.V2.plot, aes(minutes.day/60, TEMP, group = interaction(treatm), color = treatm))+
  ylim(5, 35)+
  scale_color_manual(values=c("#3F4756"))+
  annotate("rect", xmin=0, xmax=360/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  annotate("rect", xmin=960/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  geom_line(aes(minutes.day/60, TEMP, group = interaction(mo, day, treatm, tank)), size=0.4, alpha=0.5)
ggformat(daily.p.V2, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 12, print = F)
daily.p.V2 <- daily.p.V2 + theme(legend.position = "none")

## Figure iii) : V1, steady variable treatment; stacked trends ---------------

daily.p.V1 <- ggplot(dataLab.V1.plot, aes(minutes.day/60, TEMP, group = interaction(treatm), color = treatm))+
  ylim(5, 35)+
  scale_color_manual(values=c("#A3ABBD"))+
  annotate("rect", xmin=0, xmax=360/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  annotate("rect", xmin=960/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  geom_line(aes(minutes.day/60, TEMP, group = interaction(mo, day, treatm, tank)), size=0.4, alpha=0.5)
ggformat(daily.p.V1, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 12, print = F)
daily.p.V1 <- daily.p.V1 + theme(legend.position = "none")


## Figure iii) : V1 V2, variable treatments; stacked trends ---------------
dataLab.V1V2$treatm = factor(dataLab.V1V2$treatm, levels=c("V2", "V1"))

SMOOTHdaily.p.V1V2 <- ggplot(dataLab.V1V2, aes(minutes.day/60, TEMP, group = interaction(treatm), color = treatm))+
  ylim(0, 35)+
  scale_color_manual(values=c("#3F4756", "#A3ABBD"))+
  annotate("rect", xmin=0, xmax=360/60, ymin=-Inf, ymax=Inf, alpha=0.3, fill="grey90") +
  annotate("rect", xmin=960/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.3, fill="grey90") +
  geom_line(aes(minutes.day/60, TEMP, group = interaction(mo, day, treatm, tank)), size=0.1, alpha=0.3)+
  geom_smooth()+
  facet_wrap(.~treatm, nrow = 2)
ggformat(SMOOTHdaily.p.V1V2, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 12, print = T)
SMOOTHdaily.p.V1V2 <- SMOOTHdaily.p.V1V2 + theme(legend.position = "none")


## Figure iv) : V2, stochastic variable treatment; continuous trends ---------------
lab.SECTION.V2<-ggplot(dataLab.V2.plot, aes(LocalDateTime,  TEMP, group = interaction(mo, day, treatm, tank)))+
  ylim(5, 35)+
  geom_point( size=1, pch=".", color = "black")+
  geom_line( alpha=1, color = "black", lwd=0.5)
ggformat(lab.SECTION.V2, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time", print=F, size_text = 12)
lab.SECTION.V2<- lab.SECTION.V2+theme(legend.position = c(0.1, 0.8))

## Figure v) : V1, steady variable treatment; continuous trends ---------------
lab.SECTION.V1<-ggplot(dataLab.V1.plot, aes(LocalDateTime, TEMP, group = interaction(mo, day, treatm, tank)))+
  ylim(5, 35)+
  geom_point( size=1, pch=".", color = "black")+
  geom_line( alpha=1, color = "black", lwd=0.5)
ggformat(lab.SECTION.V1, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time", print=F, size_text = 12)
lab.SECTION.V1<- lab.SECTION.V1+theme(legend.position = c(0.1, 0.8))


## saving figures -------
plot_grid(lab.SECTION.V1, daily.p.V1,
          lab.SECTION.V2, daily.p.V2,
          ncol = 2, align = "hv",
          rel_widths = c(2, 1.2), rel_heights = c(1, 1)) %>%
  ggsave(filename = "./Figures/Figure1_Lab_variable.png", width = 6, height = 5)

plot_grid(daily.p,
          ncol = 1, align = "hv") %>%
  ggsave(filename = "./Figures/Figure1_Lab_static.png", width = 3.5, height = 4)

plot_grid(daily.p.V1,
          ncol = 1, align = "hv") %>%
  ggsave(filename = "./Figures/Figure1_Lab_V1.png", width = 3, height = 3)

plot_grid(daily.p.V2,
          ncol = 1, align = "hv") %>%
  ggsave(filename = "./Figures/Figure1_Lab_V2.png", width = 3, height = 3)


plot_grid(lab.SECTION,
          ncol = 1, align = "hv") %>%
  ggsave(filename = "./Figures/Figure3_Lab_ALL.png", width = 3, height = 3)



# Field: HOBO loggers in carp -------------

dataF<-read.csv("./Data/Kraskura_combined_Study_data2023.csv", header = T)
names(dataF)<-c("DateTime", "TEMP", "SiteID", "Latitude","Longitude")
dataF$DateTime<-as.POSIXct(mdy_hm(dataF$DateTime), tz = "America/Los_Angeles") # in configuartion settings this was accurately set up as PDT. not UTC timezone. - confirmed 
dataF$LocalDateTime<-dataF$DateTime # timezone: "America/Los_Angeles" or Pacific Day time "PDT"
dataF<-dataF[complete.cases(dataF),]

dataCarp<-transformTempData(dataF, MeanGroupVar = "SiteID")
dataCarp<-dataMutateMeans(dataCarp)

# section of data before the CT test in the field. 
# test dates: 2021-07-11, 2021-07-05, 2021-07-09, 2021-07-14, 2021-09-21, 2021-09-22. 
SiteTempsJul<- extractTimeFrame(dataCarp, year = 2021, minMo = 7, maxMo = 7, minDay = 1, maxDay = 14)
SiteTempsSep<- extractTimeFrame(dataCarp, year = 2021, minMo = 9, maxMo = 9, minDay = 5 , maxDay = 22)
SiteTempsJul<- SiteTempsJul %>% 
  filter(SiteID == "TestSite") %>% 
  write.csv(file = "./Data/Analysis source/SiteTempsJul.csv")
SiteTempsSep<- SiteTempsSep %>% 
  filter(SiteID == "TestSite") %>% 
  write.csv(file = "./Data/Analysis source/SiteTempsSep.csv")


# section of data presenting variation in 2020 January, cold mo only 
dataCarp.V1<- extractTimeFrame(dataCarp, year = 2020, minMo = 1, maxMo = 2, minDay = 1, maxDay = 31)
dataCarp.V1.1week<- extractTimeFrame(dataCarp, year = 2020, minMo = 1, maxMo = 1, minDay = 1, maxDay = 7)

# section of data presenting variation in 2019 Sept & Oct , warm mo, but site specific varition resembling the lab treatments 
dataCarp.V1V2<- extractTimeFrame(dataCarp, year = 2019, minMo = 10, maxMo = 10, minDay = 1, maxDay = 31)
dataCarp.V2.1week<- extractTimeFrame(dataCarp, year = 2019, minMo = 10, maxMo = 10, minDay = 1, maxDay = 7)

## Manuscript: Fig 2A ---------
daily.pFcont<-ggplot(dataCarp, aes(LocalDateTime, TEMP, color = as.factor(SiteID)))+
  # facet_wrap(.~SiteID)+
  # annotate("segment", x = as.POSIXct("2020-05-01"), xend = as.POSIXct("2020-10-01"), y = 32, yend = 32, color = "#BE647D", lwd=1, lty="dotted")+
  # annotate("segment", x = as.POSIXct("2021-05-01"), xend = as.POSIXct("2021-10-01"), y = 32, yend = 32, color = "#BE647D", lwd=2, lty="dotted")+
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 5))+
  annotate(geom = "rect", xmin = as.POSIXct("2019-10-01"), xmax = as.POSIXct("2019-11-01"),
           ymin = -Inf, ymax = Inf,  colour = "black", fill = "black", alpha = 0.1)+
  annotate(geom = "rect", xmin = as.POSIXct("2021-07-01"), xmax = as.POSIXct("2021-09-01"),
           ymin = -Inf, ymax = Inf,  colour = "#00CAFF", fill = "#00CAFF", alpha = 0.1)+
  annotate("segment", x = as.POSIXct("2019-08-01"), xend = as.POSIXct("2019-09-01"), y = 27, yend = 27,
           arrow = arrow(length = unit(.2,"cm")), color = "#BE647D", lwd=2, alpha=1)+
  annotate("segment", x = as.POSIXct("2019-08-01"), xend = as.POSIXct("2019-09-01"), y = 22, yend = 22,
           arrow = arrow(length = unit(.2,"cm")), color = "#DBA11C", lwd=2, alpha=1)+
  # annotate("segment", x = as.POSIXct("2021-05-01"), xend = as.POSIXct("2021-09-30"), y = 27, yend = 27, color = "#BE647D", lwd=2)+
  # annotate("segment", x = as.POSIXct("2021-05-01"), xend = as.POSIXct("2021-09-30"), y = 22, yend = 22, color = "#DBA11C", lwd=2)+
  annotate("segment", x = as.POSIXct("2019-08-01"), xend = as.POSIXct("2019-09-01"), y = 17, yend = 17,
           arrow = arrow(length = unit(.2,"cm")), color = "#008A60", lwd=2, alpha=1)+
  # annotate("segment", x = as.POSIXct("2021-02-20"), xend = as.POSIXct("2021-05-15"), y = 17, yend = 17, color = "#008A60", lwd=2)+
  annotate("segment", x = as.POSIXct("2019-08-01"), xend = as.POSIXct("2019-09-01"), y = 12, yend = 12,
           arrow = arrow(length = unit(.2,"cm")), color = "#00518C", lwd=2, alpha=1)+
  geom_point( size=0.4,  pch=".", show.legend = FALSE)+
  geom_point(dataCarp[dataCarp$SiteID == "SITE3",], mapping = aes(LocalDateTime, TEMP, color = as.factor(SiteID)),
             size=0.4,  pch=".", show.legend = FALSE)+ # to order the sites
  geom_point(dataCarp[dataCarp$SiteID == "SITE4",], mapping = aes(LocalDateTime, TEMP, color = as.factor(SiteID)),
             size=0.4,  pch=".", show.legend = FALSE)+ # to order the sites
 scale_color_manual(values = c("SITE2" =  "grey40", "SITE3" = "grey40", "SITE4" = "black", "TestSite" = "#00CAFF", "SITE1" = "grey"))
  # ylim(0, 40)
  # annotate("segment", x = as.POSIXct("2020-11-12"), xend = as.POSIXct("2021-03-20"), y = 12, yend = 12, color = "#00518C", lwd=2)
ggformat(daily.pFcont, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time", size_text = 12)

## Manuscript: Fig 2B ---------
SMOOTHdaily.p.CarpVAR<-ggplot(dataCarp.V1V2[(dataCarp.V1V2$SiteID == "SITE1"),])+
  geom_line(aes(DateTime, TEMP), linewidth = 0.5, color = "black")+
  scale_y_continuous(limits = c(3, 35), breaks = seq(5, 35, 5))+
  xlim(as.POSIXct("2019-10-01 12:30:00"), as.POSIXct("2019-11-08 12:30:00"))+
  annotate("segment", x = as.POSIXct("2019-11-01"), xend = as.POSIXct("2019-11-02"), y = 26.2, yend = 26.2,
           arrow = arrow(length = unit(.2,"cm"), ends = "first"), color = "#BE647D", lwd=2, alpha=1)+
  annotate("segment", x = as.POSIXct("2019-11-01"), xend = as.POSIXct("2019-11-02"), y = 17.1, yend = 17.1,
           arrow = arrow(length = unit(.2,"cm"), ends = "first"), color = "#008A60", lwd=2, alpha=1)+
  annotate("segment", x = as.POSIXct("2019-11-01"), xend = as.POSIXct("2019-11-02"), y = 11.1, yend = 11.1,
           arrow = arrow(length = unit(.2,"cm"), ends = "first"), color = "#00518C", lwd=2, alpha=1)+
  geom_text(data = dataCarp.V1V2[(dataCarp.V1V2$SiteID == "SITE1") & !duplicated(dataCarp.V1V2$SiteID), c("SiteID", "meanMAX") ],
            mapping = aes(x = as.POSIXct("2019-11-02 22:30:00"), y = 26.2, group = SiteID, label=paste("max: ", meanMAX, "ºC", sep ="")),
            hjust = 0, size = 4, color = "#BE647D")+
  geom_text(data =dataCarp.V1V2[(dataCarp.V1V2$SiteID == "SITE1") & !duplicated(dataCarp.V1V2$SiteID), c("SiteID", "meanMEAN") ],
           aes(x = as.POSIXct("2019-11-02 20:30:00"), y = 17.1, group = SiteID, label=paste("mean: ", meanMEAN, "ºC", sep ="")),
           hjust=0, size = 4, color = "#008A60")+
  geom_text(data =dataCarp.V1V2[(dataCarp.V1V2$SiteID == "SITE1") & !duplicated(dataCarp.V1V2$SiteID), c("SiteID", "meanMIN") ],
           aes(x = as.POSIXct("2019-11-02 20:30:00"), y = 11.1, group = SiteID, label=paste("min: ", meanMIN, "ºC", sep ="")),
           hjust=0, size = 4, color = "#00518C")
  # geom_line(lwd = 0.5, mapping = aes(x = 25, min_temp, group = SiteID), color = "grey") +
  # stat_summary(mapping = aes(x = 25, min_temp, group = SiteID), fun = "mean", geom = "point", size = 1, color = "grey")
ggformat(SMOOTHdaily.p.CarpVAR, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 12, print = T)
SMOOTHdaily.p.CarpVAR <- SMOOTHdaily.p.CarpVAR + theme(legend.position = "none", 
                                                       axis.title.x = element_blank())


## Figure i) : V1 V2, smooth, variable treatments; stacked trends ---------------
SMOOTHdaily.p.CarpV1V2 <- ggplot(dataCarp.V1V2[(dataCarp.V1V2$SiteID == "SITE1"),], aes(minutes.day/60, TEMP, group = interaction(SiteID)))+
  ylim(0, 35)+
  xlim(0, 28)+
  # scale_color_manual(values=c("#A3ABBD", "#3F4756"))+
  annotate("rect", xmin=0, xmax=420/60, ymin=-Inf, ymax=Inf, alpha=0.3, fill="grey90") +
  annotate("rect", xmin=1080/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.3, fill="grey90") +
  geom_line(aes(minutes.day/60, TEMP, group = interaction(mo, day, SiteID)), size=0.3, alpha=0.4, color = "black")+
  geom_smooth(se = TRUE, color = "black")+
  # facet_wrap(.~SiteID, nrow = 3)+
  geom_line(mapping = aes(x = 25, max_temp, group = SiteID), lwd = 0.5) +
  stat_summary(mapping = aes(x = 25, max_temp, group = SiteID), fun = "mean", geom = "point", size = 1)+
  geom_line(lwd = 0.5, mapping = aes(x = 24.5, mean_temp, group = SiteID), color = "grey50") +
  stat_summary(mapping = aes(x = 24.5, y = mean_temp, group = SiteID), fun = "mean", geom = "point", size = 1, color = "grey50")+
  geom_text(data =dataCarp.V1V2[(dataCarp.V1V2$SiteID == "SITE1") & !duplicated(dataCarp.V1V2$SiteID), c("SiteID", "meanMEAN") ],
           aes(x = 25.8, y = meanMEAN, group = SiteID, label=paste(meanMEAN)), hjust=0, size = 3, color = "black")+
  geom_text(data = dataCarp.V1V2[(dataCarp.V1V2$SiteID == "SITE1") & !duplicated(dataCarp.V1V2$SiteID), c("SiteID", "meanMAX") ],
            mapping = aes(x = 25.5, y = meanMAX, group = SiteID, label=paste(meanMAX)),
            hjust = 0, size = 3, color = "grey50")
  # geom_line(lwd = 0.5, mapping = aes(x = 25, min_temp, group = SiteID), color = "grey") +
  # stat_summary(mapping = aes(x = 25, min_temp, group = SiteID), fun = "mean", geom = "point", size = 1, color = "grey")
ggformat(SMOOTHdaily.p.CarpV1V2, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 10, print = T)
SMOOTHdaily.p.CarpV1V2 <- SMOOTHdaily.p.CarpV1V2 + theme(legend.position = "none")


## Figure ii) ---------
daily.p1.Carp.V2 <- ggplot(dataCarp.V1V2[!(dataCarp.V1V2$SiteID == "SITE2"),] , aes(minutes.day/60, TEMP))+
  annotate("rect", xmin=0, xmax=420/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  annotate("rect", xmin=1020/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  # geom_hline(yintercept = 27, color = "#BE647D", lty=1, lwd = 1, alpha = 0.7)+
  geom_hline(yintercept = 17, color = "#008A60", lty=1, lwd = 1, alpha = 0.7)+
  # geom_hline(yintercept = 12, color = "#00518C", lty=1, lwd = 1, alpha = 0.7)+
  # geom_hline(yintercept = 22, color = "#DBA11C", lty=1, lwd = 1, alpha = 0.7)+
  ylim(0, 35)+
  facet_wrap(.~SiteID, nrow = 3)+
  # geom_hline(yintercept = 32, color = "#BE647D", lty=2, lwd =1, alpha = 0.7)+
  geom_line(aes(minutes.day/60, TEMP, group = interaction(y, mo, day)), color = "grey10", size=0.2, alpha=1)
ggformat(daily.p1.Carp.V2, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 12)



## Figure iii) : variable treatment; continuous trends ---------------
Carp.SECTION.V2<-ggplot(dataCarp.V2.1week, aes(LocalDateTime,  TEMP, group = interaction(mo, day)))+
  ylim(0, 35)+
  geom_point( size=1, pch=".", color = "black")+
  # geom_hline(yintercept = 17, color = "#008A60", lty=1, lwd = 1, alpha = 0.7)+
  facet_wrap(.~SiteID, nrow = 3)+
  geom_line( alpha=1, color = "black", lwd=0.5)
ggformat(Carp.SECTION.V2, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time", print=F, size_text = 12)


## Figure iv) ---------
SMOOTHdaily.p1.Carp.V1 <- ggplot(dataCarp.V1[!(dataCarp.V1$SiteID == "SITE2"),], aes(minutes.day/60, TEMP, group = interaction(SiteID)))+
  ylim(0, 35)+
  # scale_color_manual(values=c("#A3ABBD", "#3F4756"))+
  annotate("rect", xmin=0, xmax=360/60, ymin=-Inf, ymax=Inf, alpha=0.3, fill="grey90") +
  annotate("rect", xmin=960/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.3, fill="grey90") +
  geom_line(aes(minutes.day/60, TEMP, group = interaction(mo, day, SiteID)), size=0.1, alpha=0.4, color = "black")+
  # geom_smooth(se = TRUE, color = "#926C00")+
  facet_wrap(.~SiteID, nrow = 3)
ggformat(SMOOTHdaily.p1.Carp.V1, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 12, print = T)


## Figure iv) : more stable variable treatment; continuous trends ---------------
Carp.SECTION.V1<-ggplot(dataCarp.V1.1week, aes(LocalDateTime,  TEMP, group = interaction(mo, day)))+
  ylim(0, 35)+
  geom_point( size=1, pch=".", color = "black")+
  facet_wrap(.~SiteID, nrow = 3)+
  geom_hline(yintercept = 17, color = "#008A60", lty=1, lwd = 1, alpha = 0.7)+
  geom_line( alpha=1, color = "black", lwd=0.5)
ggformat(Carp.SECTION.V1, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time", print=F, size_text = 12)
Carp.SECTION.V1<- Carp.SECTION.V1+theme(legend.position = c(0.1, 0.8))



## Figure v) ---------

daily.p1F<- ggplot(dataCarp.V1V2 , aes(minutes.day/60, TEMP))+
  ylim(0, 35)+
  facet_wrap(.~SiteID)+
  annotate("rect", xmin=0, xmax=420/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  annotate("rect", xmin=1020/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  geom_hline(yintercept = 27, color = "#BE647D", lty=1, lwd = 1)+
  geom_hline(yintercept = 17, color = "#008A60", lty=1, lwd = 1)+
  geom_hline(yintercept = 12, color = "#00518C", lty=1, lwd = 1)+
  geom_hline(yintercept = 22, color = "#DBA11C", lty=1, lwd = 1)+
  geom_hline(yintercept = 32, color = "#BE647D", lty=2, lwd =1, alpha = 0.7)+
  geom_line(aes(minutes.day/60, TEMP, group = interaction(y, mo, day)), color = "grey10", size=0.2, alpha=1)
ggformat(daily.p1F, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)")
daily.p1F <- daily.p1F + theme(legend.position = "top")



## saving figures -------
ggsave(filename = "./Figures/Figure1_Field_full.png", daily.pFcont, width = 9, height = 3)
ggsave(filename = "./Figures/Figure1_Field_site1.png", SMOOTHdaily.p.CarpVAR, width = 7, height = 3)

plot_grid(Carp.SECTION.V1, Carp.SECTION.V2,
          ncol = 2, align = "hv",
          rel_heights = c(1, 1)) %>% 
ggsave(filename = "./Figures/Figure1_Field_variable.png", width = 6, height = 5)

plot_grid(daily.p1.Carp.V1, daily.p1.Carp.V2,
          ncol = 2, align = "hv",
          rel_heights = c(1, 1)) %>% 
  ggsave(filename = "./Figures/Figure1_Field_variable.png", width = 6, height = 5)

# plot_grid(plot_grid(SMOOTHdaily.p.V1V2, NULL, nrow=2, rel_heights = c(0.7, 0.28)), SMOOTHdaily.p.CarpV1V2, 
#             align = "v", rel_widths = c(1.25, 1))


# data frames ----
dataF %>% 
  group_by(SiteID) %>% 
  summarize(minTEMP = min(TEMP), maxTEMP = max(TEMP), meanTEMP = mean(TEMP), dailyDelta = maxTEMP - minTEMP, 
            n = n(), 
            minDate = min(DateTime), 
            maxDate = max(DateTime))

# dataF %>% 
dataF$Date<-as.Date(dataF$DateTime)

minDay<-aggregate(dataF["TEMP"], by=dataF[c("Date", "SiteID")], min)
maxDay<-aggregate(dataF["TEMP"], by=dataF[c("Date", "SiteID")], max)
meanDay<-aggregate(dataF["TEMP"], by=dataF[c("Date", "SiteID")], mean)
deltaDay<-aggregate(dataF["TEMP"], by=dataF[c("Date", "SiteID")], range)
deltaDay$DELTA<-maxDay$TEMP-minDay$TEMP

minDay %>% 
  group_by(SiteID) %>% 
  summarize(mean_dailyMIN_TEMP = mean(TEMP), sd_dailyMIN_TEMP = sd(TEMP), n = n())

maxDay %>% 
  group_by(SiteID) %>% 
  summarize(mean_dailyMAX_TEMP = mean(TEMP), sd_dailyMAX_TEMP = sd(TEMP), n = n())

meanDay %>% 
  group_by(SiteID) %>% 
  summarize(mean_dailyMEAN_TEMP = mean(TEMP), sd_dailyMEAN_TEMP = sd(TEMP), n = n())

deltaDay %>% 
  group_by(SiteID) %>% 
  summarize(mean_dailyDELTA_TEMP = mean(DELTA), sd_dailyDELTA_TEMP = sd(DELTA), n = n())


  