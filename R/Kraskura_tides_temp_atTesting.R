
library(chron)
library(lubridate)
library(birk)

library(ggformat2) # from kraskura/ggformat github package. 
library(ggsci)
library(data.table)
library(zoo)
library(ggrepel)

library(tidyverse)
library(lattice)

library(lme4)
library(car)
library(emmeans)
library(lmerTest) 
library(factoextra) # pca visualization
library(FactoMineR) # pca 
library(here)


tide<-read.csv(file = "./Data/Tide_data_Jul_Sep2021.csv")
temp.all<-read.csv("./Data/Field_temps_CombinedData2021.csv")
temp.lab<-read.csv("./Data/Lab_TestDaytemps_2019.csv")
ct.tests<-read.csv(file = "./Data/CTtests_Lab2019_Field2021.csv")
ct.tests<-ct.tests[c(ct.tests$treatment == "V1" | ct.tests$treatment == "V2" | is.na(ct.tests$treatment)),]

# combine field and lab temps 
head(temp.all)
head(temp.lab)
temp.lab <- temp.lab[, c("LocalDateTime", "TEMP", "SiteID", "treatm", "tank")]
temp.all$DateTime<-mdy_hm(temp.all$DateTime_GMT, tz =  "UTC") 
temp.lab$DateTime<-ymd_hms(temp.lab$LocalDateTime, tz =  "UTC") - hours(8)

head(temp.all)
head(temp.lab)

colnames(temp.lab)<-colnames(temp.all)
temp.all <- rbind(temp.all, temp.lab)
temp.all$DateTime <- ymd_hms(temp.all$DateTime)

# make tide data organised
tide$DateTime<-mdy_hm(paste(tide$Date, tide$Time, sep = " "), tz =  "UTC")
tide$test.round<-"JUL"
tide[substr(tide$Date, start = 1, stop=1) == "9","test.round"]<-"SEPT"
tide$test.round<-as.factor(tide$test.round)

any(is.na(temp.all$DateTime))
testDays <- mdy(ct.tests$DateTime)

ct.tests$DateTime<-mdy_hms(paste(ct.tests$Date, ct.tests$TimeDay, sep = " "), tz =  "UTC")
# take out only one line per CT test 
ct.tests <- ct.tests[!duplicated(ct.tests$TestID),]

temp.tests <- temp.all[grepl(paste(testDays, collapse = "|"), as.Date(temp.all$DateTime)), ]
tide.tests <- tide[grepl(paste(testDays, collapse = "|"), as.Date(tide$DateTime)), ]

temp.tests2 <- temp.all[c(c(c(as.Date(temp.all$DateTime) >= as.Date("2021-07-01") & as.Date(temp.all$DateTime) <= as.Date("2021-07-20")) | 
                        c(as.Date(temp.all$DateTime) >= as.Date("2021-09-10") & as.Date(temp.all$DateTime) <= as.Date("2021-09-30"))) &
                          temp.all$SiteID == "TestSite") | 
                        c(c(as.Date(temp.all$DateTime) >= as.Date("2019-11-10") & as.Date(temp.all$DateTime) <= as.Date("2019-11-30")) &
                          temp.all$SiteID == "Lab"),]

tide.tests2 <- tide[c(c(as.Date(tide$DateTime) >= as.Date("2021-07-01") & as.Date(tide$DateTime) <= as.Date("2021-07-20")) | 
                  c(as.Date(tide$DateTime) >= as.Date("2021-09-10") & as.Date(tide$DateTime) <= as.Date("2021-09-30"))),]

# sum_daily_temp<-temp.tests2 %>%
#   dplyr:::group_by(as.Date(DateTime)) %>%
#   summarize(min_T = min(Temp), max_T = max(Temp), mean_T = mean(Temp)) %>%
#   as.data.frame()

ct.tests$Tide.m<--999
ct.tests$Env.Temp.Test<-NA
ct.tests$DateTimeTest<-ymd_hms("2022-01-12 00:00:00")
for(i in 1:nrow(ct.tests)){
  ct.tests$Tide.m[i]<-tide[which.closest(tide$DateTime, ct.tests$DateTime[i]), "Prediction_m"]
}

for(i in 1:nrow(ct.tests)){
  if(!ct.tests$Date[i] == "2021-07-05" & ct.tests$Field_Lab[i] == "FIELD"){ # dont have env temp recordings for this date
    ct.tests$Env.Temp.Test[i]<-temp.tests2[which.closest(temp.tests2$DateTime, ct.tests$DateTime[i]), "Temp"]
    ct.tests$DateTimeTest[i]<-temp.tests2[which.closest(temp.tests2$DateTime, ct.tests$DateTime[i]), "DateTime"]
  }
  
  if(!ct.tests$Date[i] == "2019-11-14" & ct.tests$Field_Lab[i] == "LAB"){ # dont have env temp recordings for this date
    ct.tests$Env.Temp.Test[i]<-temp.tests2[which.closest(temp.tests2$DateTime, ct.tests$DateTime[i]), "Temp"]
    ct.tests$DateTimeTest[i]<-temp.tests2[which.closest(temp.tests2$DateTime, ct.tests$DateTime[i]), "DateTime"]# print(i)
  }
  print(c(ct.tests$Env.Temp[i],ct.tests$Temp_test_start[i]))
}

# as.data.frame(ct.tests)
# write.csv(file = "./Field_CTtests_meansLongWenv.csv", ct.tests, row.names = FALSE)

temp.p<-ggplot(data = temp.tests, aes(DateTime, Temp, group = Longitude, color = factor(SiteID)))+
  geom_line()+
  scale_x_datetime()+
  scale_y_continuous(limits = c(0,40))+
  scale_color_aaas()+
  # scale_linetype_manual(values = c(3, 2, 1))+
  facet_wrap(.~as.Date(DateTime), scales = "free_x", nrow=4)
ggformat(temp.p, title = "", y_title = "Temperature (ºC)", x_title = "Date Time", size_text = 10)

tides.p<-ggplot(data = tide.tests, aes(DateTime, Prediction_m))+
  geom_line()+
  scale_x_datetime()+
  facet_wrap(.~as.Date(DateTime), scales = "free_x")
ggformat(tides.p, title = "", y_title = "Tide (m)", x_title = "Date Time")

# # combo plot
nightRanges <- data.frame(
  from=ymd_hms(c('2021-07-05 20:00:00', '2021-07-09 20:00:00','2021-07-11 20:00:00', '2021-07-14 20:00:00',
                '2021-09-21 20:00:00', '2021-09-22 20:00:00',
                '2021-07-05 00:00:01', '2021-07-09 00:00:01', '2021-07-11 00:00:01', '2021-07-14 00:00:01',
                '2021-09-21 00:00:01', '2021-09-22 00:00:01')),
  to=ymd_hms(c('2021-07-05 23:59:59', '2021-07-09 23:59:59',  '2021-07-11 23:59:59', '2021-07-14 23:59:59',
               '2021-09-21 23:59:59', '2021-09-22 23:59:59',
               '2021-07-05 06:00:00', '2021-07-09 06:00:00','2021-07-11 06:00:00', '2021-07-14 06:00:00',
               '2021-09-21 06:00:00', '2021-09-22 06:00:00')),
  DateTime = as.Date(c('2021-07-05 20:00:00', '2021-07-09 20:00:00','2021-07-11 20:00:00', '2021-07-14 20:00:00',
                       '2021-09-21 20:00:00', '2021-09-22 20:00:00',
                       '2021-07-05 20:00:00', '2021-07-09 20:00:00','2021-07-11 20:00:00', '2021-07-14 20:00:00',
                       '2021-09-21 20:00:00', '2021-09-22 20:00:00')))

# temp.tests

# ct.tests$DateTime
# tide2$DateTime
# temp.tests$DateTime

t.p<-ggplot()+
  geom_rect(data = nightRanges, aes(xmin = from, xmax = to, ymin = -Inf, ymax = Inf), alpha = 0.1)+
  # geom_line(data = tide.tests, aes(DateTime, Prediction_m + 20, group = NULL, linetype = NULL),
            # color = "grey30", alpha=0.4, linewidth=1)+
  geom_line(data = temp.tests, aes(DateTime, Temp, group = Longitude, color = factor(SiteID)), linewidth = 0.5)+
  scale_y_continuous(limits = c(14,40), sec.axis = sec_axis(~ (. - 20), name = "Tide (m)")) +
  # scale_linetype_manual(values = c(3, 2, 1))+
  facet_wrap(.~as.Date(DateTime), scales = "free_x") +
             # labeller = as_labeller(c(`2021-07-05` = "July 5",
             #                                `2021-07-09`= "July 9",
             #                                `2021-07-11`= "July 11",
             #                                `2021-07-14`= "July 14",
             #                                `2021-09-21`= "September 21",
             #                                `2021-09-22` = "September 22"))
# )+
  scale_x_datetime(date_breaks = "4 hours", date_labels = "%H:%M")+
  geom_point(data = ct.tests[ct.tests$Test=="CTmin" , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "blue", color = "black", pch=23, size=3.5)+
  geom_point(data = ct.tests[ct.tests$Test=="CTmax" , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "red3", color = "red3", pch=21, size=2)
ggformat(t.p, title = "", y_title = "Temperature (ºC)", x_title = "Time of the day (24h)", size_text = 10)
t.p<-t.p+theme(axis.text.x = element_text( size=8))


# ggsave(t.p, filename = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Carpinteria marsh fish /DataAnalysis/PLOTS/Fig_Conceptual_fieldPart.png", width = 8, height = 4)

###  tides ------
t.p.tide2.Jul<-ggplot(data = tide2[c(as.Date(tide2$DateTime) >= as.Date("2021-07-03") & as.Date(tide2$DateTime) <= as.Date("2021-07-12")),])+
  geom_line(aes(x=DateTime, y=Prediction_m), color = "grey", alpha=0.4, lwd=1)+
  # geom_line(data = temp.tests2, aes(DateTime, Temp), lwd=0.5)+
  # scale_x_datetime(date_breaks = "1 day", date_labels = "%m-%d")+
  geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmin") & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-01") & as.Date(ct.tests$DateTime) <= as.Date("2021-07-20")), "DateTimeTest"],
             color = "blue", alpha = 0.5)+
  geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmax") & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-01") & as.Date(ct.tests$DateTime) <= as.Date("2021-07-20")), "DateTimeTest"],
             color = "red", alpha = 0.5)
ggformat(t.p.tide2.Jul, title = "", y_title = "Tide (m)", x_title = "Date")
t.p.tide2.Jul<-t.p.tide2.Jul+theme(axis.text.x = element_text( size=8))
# t.p.tide2.Jul

t.p.tide2.Sep<-ggplot(data = tide2[c(as.Date(tide2$DateTime) >= as.Date("2021-09-19") & as.Date(tide2$DateTime) <= as.Date("2021-09-24")),])+
  geom_line(aes(x=DateTime, y=Prediction_m), color = "grey", alpha=0.4, lwd=1)+
  # geom_line(data = temp.tests2, aes(DateTime, Temp), lwd=0.5)+
  # scale_x_datetime(date_breaks = "1 day", date_labels = "%m-%d")+
  geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmin") & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-31")), "DateTime"], color = "blue", alpha = 0.5)+
  geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmax") & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-31")), "DateTime"], color = "red", alpha = 0.5)
ggformat(t.p.tide2.Sep, title = "", y_title = "Tide (m)", x_title = "Date")
t.p.tide2.Sep<-t.p.tide2.Sep + theme(axis.text.x = element_text( size=8),
                                                 axis.text.y = element_blank(), 
                                                 axis.title.y = element_blank())



### temperatures  ------

t.p.temp.tests2.Jul<-ggplot(data = temp.tests2[c(as.Date(temp.tests2$DateTime) >= as.Date("2021-07-03") & as.Date(temp.tests2$DateTime) <= as.Date("2021-07-14")),])+
  geom_line(aes(DateTime, Temp), linewidth = 0.8)+
  scale_y_continuous(limits = c(12, 35), breaks = seq(12, 35, 3))+
  scale_x_datetime(date_breaks = "1 day", date_labels = "%b-%d")+
  geom_point(data = ct.tests[ct.tests$Test=="CTmin" & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-01") & as.Date(ct.tests$DateTime) <= as.Date("2021-07-20")) , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "dodgerblue", color = "dodgerblue", pch=21, size=2)+
  geom_point(data = ct.tests[ct.tests$Test=="CTmax" & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-01") & as.Date(ct.tests$DateTime) <= as.Date("2021-07-20")) , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "red3", color = "red3", pch=21, size=2)
  # geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmin") & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-01") & as.Date(ct.tests$DateTime) <= as.Date("2021-07-20")), "DateTime"], color = "blue", alpha = 0.5)+
  # geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmax") & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-01") & as.Date(ct.tests$DateTime) <= as.Date("2021-07-20")), "DateTime"], color = "red", alpha = 0.5)
ggformat(t.p.temp.tests2.Jul, title = "", y_title = "Temperature (ºC)", x_title = "Date", print = F, size_text = 12)
t.p.temp.tests2.Jul<-t.p.temp.tests2.Jul + theme(axis.text.x = element_text(size=8),
                                                 axis.text.y = element_blank(), 
                                                 axis.title.y = element_blank())
# t.p.temp.tests2.Jul


t.p.temp.tests2.Sep<-ggplot(data = temp.tests2[c(as.Date(temp.tests2$DateTime) >= as.Date("2021-09-20") & as.Date(temp.tests2$DateTime) <= as.Date("2021-09-22")),])+
  geom_line(aes(DateTime, Temp), linewidth = 0.8)+
  scale_y_continuous(limits = c(12, 35), breaks = seq(12, 35, 3))+
  scale_x_datetime(date_breaks = "1 day", date_labels = "%b-%d")+
  geom_point(data = ct.tests[ct.tests$Test=="CTmin" & c(as.Date(ct.tests$DateTime) >= as.Date("2021-09-20") & as.Date(ct.tests$DateTime) <= as.Date("2021-09-22")) , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "dodgerblue", color = "dodgerblue", pch=21, size=2)+
  geom_point(data = ct.tests[ct.tests$Test=="CTmax" & c(as.Date(ct.tests$DateTime) >= as.Date("2021-09-20") & as.Date(ct.tests$DateTime) <= as.Date("2021-09-22")) , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "red3", color = "red3", pch=21, size=2)
  # geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmin") & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-01") & as.Date(ct.tests$DateTime) <= as.Date("2021-07-20")), "DateTime"], color = "blue", alpha = 0.5)+
  # geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmax") & c(as.Date(ct.tests$DateTime) >= as.Date("2021-07-01") & as.Date(ct.tests$DateTime) <= as.Date("2021-07-20")), "DateTime"], color = "red", alpha = 0.5)
ggformat(t.p.temp.tests2.Sep, title = "", y_title = "Temperature (ºC)", x_title = "Date", print = F, size_text = 12)
t.p.temp.tests2.Sep<-t.p.temp.tests2.Sep + theme(axis.text.x = element_text( size=8),
                                                 axis.text.y = element_blank(), 
                                                 axis.title.y = element_blank(), 
                                                 axis.title.x = element_blank())
# t.p.temp.tests2.Sep


t.p.temp.lab.tests<-ggplot(data = temp.tests2[c(as.Date(temp.tests2$DateTime) >= as.Date("2019-11-9") & as.Date(temp.tests2$DateTime) <= as.Date("2019-11-13")),])+
  geom_line(aes(DateTime, Temp, group  = Longitude), lwd=0.5, alpha=0.2)+
  scale_y_continuous(limits = c(12, 35), breaks = seq(12, 35, 3))+
  stat_smooth(aes(DateTime, Temp), color = "black", linewidth = 0.8)+
  scale_x_datetime(date_breaks = "1 day", date_labels = "%b-%d")+
  geom_point(data = ct.tests[ct.tests$Test=="CTmin" & c(as.Date(ct.tests$DateTime) >= as.Date("2019-11-9") & as.Date(ct.tests$DateTime) <= as.Date("2019-11-13")) , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "dodgerblue", color = "dodgerblue", pch=21, size=2)+
  geom_point(data = ct.tests[ct.tests$Test=="CTmax" & c(as.Date(ct.tests$DateTime) >= as.Date("2019-11-9") & as.Date(ct.tests$DateTime) <= as.Date("2019-11-13")) , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "red3", color = "red3", pch=21, size=2)
  # geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmin") & c(as.Date(ct.tests$DateTime) >= as.Date("2019-11-9") & as.Date(ct.tests$DateTime) <= as.Date("2019-11-13")), "DateTime"], color = "blue", alpha = 0.5)+
  # geom_vline(xintercept = ct.tests[c(ct.tests$Test=="CTmax") & c(as.Date(ct.tests$DateTime) >= as.Date("2019-11-9") & as.Date(ct.tests$DateTime) <= as.Date("2019-11-13")), "DateTime"], color = "red", alpha = 0.5)
ggformat(t.p.temp.lab.tests, title = "", y_title = "Temperature (ºC)", x_title = "Date", print = F, size_text = 12)
t.p.temp.lab.tests <- t.p.temp.lab.tests + theme(axis.text.x = element_text( size=8),
                                                 axis.text.y = element_text( size=8),
                                                 axis.title.x = element_blank())
# t.p.temp.lab.tests


t.p.temp.lab.tests2<-ggplot(data = temp.tests2[c(as.Date(temp.tests2$DateTime) >= as.Date("2019-11-24") & as.Date(temp.tests2$DateTime) <= as.Date("2019-11-25")),])+
  geom_line(aes(DateTime, Temp, group  = Longitude), lwd = 0.5, alpha = 0.2)+
  scale_y_continuous(limits = c(12, 35), breaks = seq(12, 35, 3))+
  stat_smooth(aes(DateTime, Temp), color = "black", linewidth = 0.8)+
  scale_x_datetime(date_breaks = "1 day", date_labels = "%b-%d")+
  geom_point(data = ct.tests[ct.tests$Test=="CTmin" & c(as.Date(ct.tests$DateTime) >= as.Date("2019-11-24") & as.Date(ct.tests$DateTime) <= as.Date("2019-11-25")) , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "dodgerblue", color = "dodgerblue", pch=21, size=2)+
  geom_point(data = ct.tests[ct.tests$Test=="CTmax" & c(as.Date(ct.tests$DateTime) >= as.Date("2019-11-24x") & as.Date(ct.tests$DateTime) <= as.Date("2019-11-25")) , ], mapping = aes(x=DateTime, y=Env.Temp.Test), fill = "red3", color = "red3", pch=21, size=2)
ggformat(t.p.temp.lab.tests2, title = "", y_title = "Temperature (ºC)", x_title = "Date", print = F, size_text = 12)
t.p.temp.lab.tests2<-t.p.temp.lab.tests2 + theme(axis.text.x = element_text( size=8),
                                                 axis.text.y = element_blank(), 
                                                 axis.title.y = element_blank(),
                                                 axis.title.x = element_blank())
# t.p.temp.lab.tests2

cowplot::plot_grid(t.p.temp.lab.tests, t.p.temp.lab.tests2, t.p.temp.tests2.Jul, t.p.temp.tests2.Sep,
                  nrow = 1,
                  align = "h", 
                  rel_widths = c(0.3, 0.2, 0.8, 0.3)) %>% 
  ggsave(filename = "./Figures/Figure2_Lab_Field_tests_temps.png", width = 10, height = 2)

# Lab temp matching time frames and testes






