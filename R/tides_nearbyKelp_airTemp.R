



# Tides in carp -------------
dataF<-read.csv("/Users/kristakraskura/Github_repositories/Salt-marsh-fish-physio/Data/csv/Kraskura_combined_Study_data2021.csv", header = T)
names(dataF)<-c("DateTime", "TEMP", "SiteID", "Latitude","Longitude")
dataF$DateTime<-as.POSIXct(mdy_hm(dataF$DateTime), tz = "America/Los_Angeles") # in configuartion settings this was accurately set up as PDT. not UTC timezone. - confirmed 
dataF$LocalDateTime<-dataF$DateTime # timezone: "America/Los_Angeles" or Pacific Day time "PDT"

dataCarp<-transformTempData(dataF, MeanGroupVar = "SiteID")
dataCarp<-dataCarp[complete.cases(dataCarp),]
 
  
  
  
# LTER kelp site; nearshore ocean temps:  ----------
# ************** 
# https://www.neonscience.org/tabular-time-series
# https://ourcodingclub.github.io/tutorials/time/
# https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sbc.13.28
  
  
dataKelp<-read.csv("/Users/kristakraskura/Github_repositories/LTER-data/LTER data/Temperature/Bottom_temp_all_years_20220729.csv")
dates<-dataKelp$DATE_LOCAL
time<-dataKelp$TIME_LOCAL
LocalDateTime<-paste(dates, time)
dataKelp$LocalDateTime<-as.POSIXct(strptime(LocalDateTime, "%Y-%m-%d %H:%M"), tz ="America/Los_Angeles")
dataKelp<-dataKelp[dataKelp$SITE =="CARP",] # take only Carpinteria kelp forest data 
dataKelp$SITE<-as.factor(dataKelp$SITE)
names(dataKelp)<- c("SITE", "SERIAL", "DATE_LOCAL", "TIME_LOCAL", "TEMP", "LocalDateTime")
# sanity checks
names(dataKelp)
levels(dataKelp$SITE)
length(which(dataKelp[]==-99999)) # ==0!  good all dataKelp are there, no missing values 

dataKelp<-transformTempData(dataKelp, MeanGroupVar = "SITE")

# # no year, but site, all sites together
# temp_sum_Y <- dataKelp %>%
#   dplyr:::group_by(DATE_LOCAL, YEAR) %>%
#   summarize(mean_temp = mean(TEMP_C), min_temp = min(TEMP_C), max_temp = max(TEMP_C), var_temp = var(TEMP_C),
#             sd_temp = sd(TEMP_C),
#             mean_YEAR = mean (YEAR), mean_DAY = mean (DAY), mean_MO = mean(MONTH), mean_HR = mean (HR), .groups = "keep")
# 
# temp_sum_Y$MO_DAY<-as.factor(substr(temp_sum_Y$DATE_LOCAL, start=6, stop=10))
# temp_sum_Y<-as.dataKelp.frame(temp_sum_Y)

## Figure i) ---------


dataKelp.V1<- extractTimeFrame(dataKelp, year = 2020, minMo = 1, maxMo = 2, minDay = 1, maxDay = 31)
dataKelp.V1.1week<- extractTimeFrame(dataKelp, year = 2020, minMo = 1, maxMo = 1, minDay = 1, maxDay = 7)

dataKelp.V2<- extractTimeFrame(dataKelp, year = 2019, minMo = 9, maxMo = 10, minDay = 1, maxDay = 31)
dataKelp.V2.1week<- extractTimeFrame(dataKelp, year = 2019, minMo = 10, maxMo = 10, minDay = 1, maxDay = 7)


## Figure i) ---------
daily.p1.Kelp.V2 <- ggplot(dataKelp.V2 , aes(minutes.day/60, TEMP))+
  annotate("rect", xmin=0, xmax=420/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  annotate("rect", xmin=1020/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  # geom_hline(yintercept = 27, color = "#BE647D", lty=1, lwd = 1, alpha = 0.7)+
  geom_hline(yintercept = 17, color = "#00C59E", lty=1, lwd = 1, alpha = 0.7)+
  # geom_hline(yintercept = 12, color = "#00518C", lty=1, lwd = 1, alpha = 0.7)+
  # geom_hline(yintercept = 22, color = "#DBA11C", lty=1, lwd = 1, alpha = 0.7)+
  ylim(5, 35)+
  # geom_hline(yintercept = 32, color = "#BE647D", lty=2, lwd =1, alpha = 0.7)+
  geom_line(aes(minutes.day/60, TEMP, group = interaction(y, mo, day)), color = "grey10", size=0.2, alpha=1)
ggformat(daily.p1.Kelp.V2, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 12)
# daily.p1.Kelp.V2 <- daily.p1.Kelp.V2 + theme(legend.position = "top",
#                                      axis.text.y = element_blank(),
#                                      axis.title.y = element_blank())


## Figure iv) : variable treatment; continuous trends ---------------
Kelp.SECTION.V2<-ggplot(dataKelp.V2.1week, aes(LocalDateTime,  TEMP, group = interaction(mo, day)))+
  ylim(5, 35)+
  geom_point( size=1, pch=".", color = "black")+
  geom_hline(yintercept = 17, color = "#00C59E", lty=1, lwd = 1, alpha = 0.7)+
  geom_line( alpha=1, color = "black", lwd=0.5)
ggformat(Kelp.SECTION.V2, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time", print=F, size_text = 12)
Kelp.SECTION.V2<- Kelp.SECTION.V2+theme(legend.position = c(0.1, 0.8))



## Figure i) ---------
daily.p1.Kelp.V1 <- ggplot(dataKelp.V1 , aes(minutes.day/60, TEMP))+
  annotate("rect", xmin=0, xmax=420/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  annotate("rect", xmin=1020/60, xmax=1440/60, ymin=-Inf, ymax=Inf, alpha=0.5, fill="grey90") +
  ylim(5, 35)+
  geom_hline(yintercept = 17, color = "#00C59E", lty=1, lwd = 1, alpha = 0.7)+
  # geom_hline(yintercept = 32, color = "#BE647D", lty=2, lwd =1, alpha = 0.7)+
  geom_line(aes(minutes.day/60, TEMP, group = interaction(y, mo, day)), color = "grey10", size=0.2, alpha=1)
ggformat(daily.p1.Kelp.V1, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time in the day (h)", size_text = 12)
# daily.p1.Kelp.V2 <- daily.p1.Kelp.V2 + theme(legend.position = "top",
#                                      axis.text.y = element_blank(),
#                                      axis.title.y = element_blank())

## Figure iv) : more stable variable treatment; continuous trends ---------------
Kelp.SECTION.V1<-ggplot(dataKelp.V1.1week, aes(LocalDateTime,  TEMP, group = interaction(mo, day)))+
  ylim(5, 35)+
  geom_point( size=1, pch=".", color = "black")+
  geom_hline(yintercept = 17, color = "#00C59E", lty=1, lwd = 1, alpha = 0.7)+
  geom_line( alpha=1, color = "black", lwd=0.5)
ggformat(Kelp.SECTION.V1, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time", print=F, size_text = 12)
Kelp.SECTION.V1<- Kelp.SECTION.V1+theme(legend.position = c(0.1, 0.8))


## Figure v) ---------
daily.pFcontKelp<-ggplot(dataKelp[dataKelp$y == 2019 | dataKelp$y == 2020, ], aes(LocalDateTime, TEMP))+
  geom_point( size=0.4,  pch=".", show.legend = FALSE, color = "#4D8076")+
  ylim(0, 40)+
  # facet_wrap(.~SiteID)+
  # annotate("segment", x = as.POSIXct("2020-05-01"), xend = as.POSIXct("2020-10-01"), y = 32, yend = 32, color = "#BE647D", lwd=1, lty="dotted")+
  # annotate("segment", x = as.POSIXct("2021-05-01"), xend = as.POSIXct("2021-10-01"), y = 32, yend = 32, color = "#BE647D", lwd=2, lty="dotted")+
  annotate("segment", x = as.POSIXct("2020-05-01"), xend = as.POSIXct("2020-10-01"), y = 27, yend = 27, color = "#BE647D", lwd=2)+
  annotate("segment", x = as.POSIXct("2020-05-01"), xend = as.POSIXct("2020-10-01"), y = 22, yend = 22, color = "#DBA11C", lwd=2)+
  # annotate("segment", x = as.POSIXct("2021-05-01"), xend = as.POSIXct("2021-09-30"), y = 27, yend = 27, color = "#BE647D", lwd=2)+
  # annotate("segment", x = as.POSIXct("2021-05-01"), xend = as.POSIXct("2021-09-30"), y = 22, yend = 22, color = "#DBA11C", lwd=2)+
  annotate("segment", x = as.POSIXct("2020-01-01"), xend = as.POSIXct("2020-06-01"), y = 17, yend = 17, color = "#00C59E", lwd=2)+
  # annotate("segment", x = as.POSIXct("2021-02-20"), xend = as.POSIXct("2021-05-15"), y = 17, yend = 17, color = "#00C59E", lwd=2)+
  annotate("segment", x = as.POSIXct("2019-11-01"), xend = as.POSIXct("2020-04-01"), y = 12, yend = 12, color = "#00518C", lwd=2)
# annotate("segment", x = as.POSIXct("2020-11-12"), xend = as.POSIXct("2021-03-20"), y = 12, yend = 12, color = "#00518C", lwd=2)
ggformat(daily.pFcontKelp, title = "", y_title = expression(Temperature~(degree*C)), x_title = "Time", size_text = 12)



# ------
daily.kelp<-ggplot(dataKelp[dataKelp$mo == 9 & c(dataKelp$y == 2019 | dataKelp$y == 2020 | dataKelp$y == 2021),], aes(LocalDateTime, TEMP))+
  ylim(0, 35)+
  geom_point( size=1, pch=".", color = "black")+
  geom_line( alpha=1, color = "black", lwd=0.5)+
  facet_wrap(.~y, scales = "free")
ggformat(daily.kelp,
         title = "",
         y_title = expression(Temperature~(degree*C)),
         x_title = "Time", print=T, size_text = 12)





