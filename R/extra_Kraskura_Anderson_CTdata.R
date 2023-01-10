

library(ggformat2) # from kraskura/ggformat github package. 
library(ggsci)
library(chron)
library(lubridate)
library(data.table)
library(zoo)
library(ggrepel)
library(lattice)
library(lme4)
library(car)
library(emmeans)
library(merTools)
library(lmerTest) 
library(tidyverse)

library(here)

# function used to calculate ∆BIC scores and re-order based on the lowest score. 
BICdelta<-function(BICtable){
  BIC.t <- BICtable [order(BICtable$BIC), ]
  BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
  return( BIC.t)
}

# set colors 
cols.klab<- c("#00518C", "#008A60", "#DBA11C", "#BE647D", "black", "#A3ABBD")  # 27, 22, 17, 12, V1, V2
# cols.tide<-c("#3D459F", "#BE647D", "#001041")


# DATA PREP -----
## 1. CTmin and CTmax -----------
data<-read.csv("./Data/CTtests_Lab2019_Field2021.csv")
data<-data %>%
  mutate_if(is.integer, as.numeric) %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(Date = as.Date(Date, "%m/%d/%y"),
         TimeDay = as.times(TimeDay),
         TimeDay2 = as.numeric(paste(substr(TimeDay, start = 1, stop=2), substr(TimeDay,start = 4, stop = 5), sep="")),
         DateTime = ymd_hms(paste(Date, TimeDay, sep = " "), tz =  "America/Los_Angeles")) %>% 
  mutate(Date = ymd(Date), 
         K_cond = mass_mg*100/TL_mm^3)
data$temp<-as.numeric(as.character(data$temp_treatment)) # NA are variable treatments 
data$tank<-factor(data$tank) # NA are variable treatments 

# import field temps to get daily min and max temps.
fieldTempsJul<- read.csv("./Data/SiteTempsJul.csv")
fieldTempsJul<- fieldTempsJul %>% 
  dplyr::select(max_temp, min_temp, mean_temp, range_temp, y, day, DateTime) %>% 
  distinct(day, .keep_all = TRUE) %>% 
  mutate(max.Env.Temp = rollmeanr(max_temp, 3, na.pad = T), 
         min.Env.Temp = rollmeanr(min_temp, 3, na.pad = T), 
         mean.Env.Temp = rollmeanr(mean_temp, 3, na.pad = T), 
         delta.T = max.Env.Temp - min.Env.Temp,
         Date = ymd(as.character(substr(DateTime, start = 1, stop=10)))) %>% 
  filter(day == 9 | day == 11 | day == 14) %>% 
  dplyr::select(max.Env.Temp, min.Env.Temp, mean.Env.Temp, delta.T, day, Date)
  
fieldTempsSep<- read.csv("./Data/SiteTempsSep.csv")
fieldTempsSep<- fieldTempsSep %>% 
  dplyr::select(max_temp, min_temp, mean_temp, range_temp, y, day, DateTime) %>% 
  distinct(day, .keep_all = TRUE) %>% 
  mutate(max.Env.Temp = rollmeanr(max_temp, 3, na.pad = T), 
         min.Env.Temp = rollmeanr(min_temp, 3, na.pad = T), 
         mean.Env.Temp = rollmeanr(mean_temp, 3, na.pad = T), 
         delta.T = max.Env.Temp - min.Env.Temp,
         Date = ymd(as.character(substr(DateTime, start = 1, stop=10)))) %>% 
  filter(day == 21 | day == 22) %>% 
  dplyr::select(max.Env.Temp, min.Env.Temp, mean.Env.Temp, delta.T, day, Date)
   
fieldTemps<-rbind(fieldTempsSep, fieldTempsJul)
for ( i in 1: nrow(data)){
  if(is.na(data$max.Env.Temp[i]) & !data$Date[i] == "2021-07-05"){
    data$max.Env.Temp[i]<-fieldTemps[which(fieldTemps$Date == data$Date[i]), "max.Env.Temp"]
    data$min.Env.Temp[i]<-fieldTemps[which(fieldTemps$Date == data$Date[i]), "min.Env.Temp"]
    data$mean.Env.Temp[i]<-fieldTemps[which(fieldTemps$Date == data$Date[i]), "mean.Env.Temp"]
    data$delta.T[i]<-fieldTemps[which(fieldTemps$Date == data$Date[i]), "delta.T"]
  }
}     


data.sum<-data %>%
  dplyr::group_by(TestID, TestID2, Test, temp) %>%
  summarise(mean_CT = mean(temp_tolerance, na.rm = TRUE), sd_CT = sd(temp_tolerance, na.rm = TRUE), min_CT = min(temp_tolerance, na.rm = TRUE), max_CT = max(temp_tolerance, na.rm = TRUE),
            mean_mass = mean(mass_mg, na.rm = TRUE), sd_mass = sd(mass_mg, na.rm = TRUE), min_mass = min(mass_mg, na.rm = TRUE), max_mass = max(mass_mg, na.rm = TRUE),
            mean_TL = mean(TL_mm, na.rm = TRUE), sd_TL = sd(TL_mm, na.rm = TRUE), min_TL = min(TL_mm, na.rm = TRUE), max_TL = max(TL_mm, na.rm = TRUE),
            n=length(!is.na(temp_tolerance))) %>% 
  as.data.frame()



## 2. Weights (lab 2019 only): -----------
data2<-read.csv("./Data/Lab_Weights_GobiesFAll2019.csv")
data2$tank<-as.factor(as.character(data2$tank))
data2<-data2[complete.cases(data2),]
data2$date<-as.POSIXct(strptime(dates(as.character(data2$date)), format = "%m/%d/%y", tz = "America/Los_Angeles"))# assign time points of weighing the fish 
data2$days_bw_measures <- c(data2$date - (as.POSIXct(strptime("10/18/19", format = "%m/%d/%y", tz = "America/Los_Angeles")))) / (60*60*24) 

data2 <- data2 %>% 
  mutate(K_cond = mass_mg*100/TL_mm^3,
         GR_mg_d = NA, 
         timepoint = ifelse(days_bw_measures == 0, 1, # timepoint 1 = all fish measured on oct 18 ("10/18/19")
                            ifelse(round(days_bw_measures) == 28, 2, 3))) %>% # timepoint 2 = all fish measured on Nov 15 ("11/15/19") (28.04167 days in between)
  mutate(timepoint = factor(timepoint)) # the original is in seconds, adjust from there) %>% # calculate condition factor K= 100weight / (L^3)
data2$temp_treatment[data2$temp_treatment == "var1"]<-"V1"
data2$temp_treatment[data2$temp_treatment == "var2"]<-"V2"

data2.sum.all<-data2 %>%
  dplyr::group_by(treatment, temp_treatment, tank, timepoint) %>% 
  summarise(mean_mg=mean(mass_mg), sd_mg=sd(mass_mg), min_mg=min(mass_mg), max_mg=max(mass_mg),
            mean_mm=mean(TL_mm), sd_mm=sd(TL_mm), min_mm=min(TL_mm), max_mm=max(TL_mm), n=length(!is.na(mass_mg)))
data2.sum<- data2.sum.all %>% 
  dplyr::filter(timepoint!="3") %>%  # take out the individuals that were measured three times (all post Nov 15), CT tests because those were not fed, unfair comparison
  as.data.frame() 

# long to wide format based no the timepoints 
# don't use the ones after CT tests because those were not fed, unfair comparison
timepoint1<-data2.sum[data2.sum$timepoint==1,]
timepoint2<-data2.sum[data2.sum$timepoint==2,]
data2.sum.w<-merge(timepoint1, timepoint2, by = "tank" , all.x = TRUE)
data2.sum.w$days_bw_measures.y<-28 # days between timepoints 1 and 2
data2.sum.w<-data2.sum.w[complete.cases(data2.sum.w),] # two tanks that did not have any fish left due to mortality, tank 10 (HOT), and tank 14 (V1)
data2.sum.w$GR_mg_d<-(as.numeric(data2.sum.w$mean_mg.y)-as.numeric(data2.sum.w$mean_mg.x)) / as.numeric(data2.sum.w$days_bw_measures.y)
data2.sum.w$GR_mm_d<-(as.numeric(data2.sum.w$mean_mm.y)-as.numeric(data2.sum.w$mean_mm.x)) / as.numeric(data2.sum.w$days_bw_measures.y)

# don't include tanks that have < 10 fish at the time of weighing
# take only static treatments 
growth.mod.data<-data2.sum.w %>% 
  dplyr::filter(temp_treatment.x=="12" | temp_treatment.x=="17" | temp_treatment.x=="22" | temp_treatment.x=="27") %>% 
  mutate_if(is.integer, as.numeric) %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(temp = as.numeric(as.character(temp_treatment.x))) 

growth.mod.mean<-data2.sum.w %>% 
  mutate_if(is.character, as.factor) %>%
  dplyr::group_by(temp_treatment.x) %>%
  # mutate_if(is.integer, as.numeric) %>% 
  # mutate(temp = as.numeric(as.character(temp_treatment.x))) %>% 
  summarise(mean_GR=mean(GR_mg_d), sd_GR=sd(GR_mg_d), min_GR=min(GR_mg_d), max_GR=max(GR_mg_d),
            mean_GRmm=mean(GR_mm_d), sd_GRmm=sd(GR_mm_d), min_GRmm=min(GR_mm_d), max_GRmm=max(GR_mm_d), n=length(!is.na(GR_mg_d))) %>% 
  mutate(temp = as.numeric(as.character(temp_treatment.x))) %>% # NAs on V1 and V2
  as.data.frame()


# STATS --------
## i) NOT used growth mg  ----------------------
# mod.poly.growth <- lm(GR_mg_d ~ poly(temp,2), data = growth.mod.data)
# plot(mod.poly.growth)
predict(mod.poly.growth)
summary(mod.poly.growth)

pred.datamg<-as.data.frame(expand.grid(temp =seq(12, 27, 0.2)))
pred.datamg$pred.GR<-predict(mod.poly.growth, newdata = pred.datamg)
pred.datamg$pred.GR.se<-predict(mod.poly.growth, newdata = pred.datamg, se.fit = TRUE)[[2]]

## ii) NOT used CTminmax: field only CT tests ----------------------
# are CT results different between times of measurements, days, etc.?
data.minF<-data[c(data$Field_Lab == "FIELD" & data$Test=="CTmin"),] # all Field CTmax
data.maxF<-data[c(data$Field_Lab == "FIELD" & data$Test=="CTmax"),] # all Field CTmin
modF1_min<-lm(temp_tolerance ~  TestID, data = data.minF)
modF1_max<-lm(temp_tolerance ~  TestID, data = data.maxF)
# plot(modF1_min)
# plot(modF1_max)

# >> yes very different 
car::Anova(modF1_min)
car::Anova(modF1_max)


## iii) NOT used: all tests (Variables and statics; categorical treatment predictor)----------
data.max<-data[c(data$Field_Lab == "LAB" & data$Test=="CTmax"),] # all LAB CTmax
# mod1max.all<-lmer(temp_tolerance ~ treatment + (1|tank), REML = FALSE, data = data.max) # leads to singular fit
mod1max.all<-lm(temp_tolerance ~ treatment , data = data.max) # leads to singular fit
plot(mod1max.all)
# qqmath(mod1max.all)
# summary(mod1max.all)

data.min<-data[c(data$Field_Lab == "LAB" & data$Test=="CTmin" ),] # all LAB CTmin
# mod1min.all<-lmer(temp_tolerance ~  treatment  + (1|tank), REML = FALSE, data = data.min) # leads to singular fit
mod1min.all<-lm(temp_tolerance ~  treatment , data = data.min) # leads to singular fit
plot(mod1min.all) 
# qqmath(mod1min.all)
# summary(mod1min.all)

# Type II anovas 
car::Anova(mod1max.all)
car::Anova(mod1min.all)
# post hocs 
emmeans(mod1max.all, "treatment")
emmeans(mod1min.all, "treatment")
emmeans(mod1max.all, pairwise ~ treatment)
emmeans(mod1min.all, pairwise ~ treatment)


## iv). FIELD + LAB: Variable treatments only  -----------

data.VARmin<-data[c(data$Field_Lab == "FIELD" | data$treatment == "V1" | data$treatment == "V2") & data$Test == "CTmin",  ]
data.VARmax<-data[c(data$Field_Lab == "FIELD" | data$treatment == "V1" | data$treatment == "V2") & data$Test == "CTmax",  ]

# CTmin -- All variability tests together, different explanatory variables:
mod.envT.CTmin.D1<-lmer(temp_tolerance ~ Temp_test_start + (1|TestID), data =data.VARmin, REML = FALSE)
mod.envT.CTmin.D2<-lmer(temp_tolerance ~ delta.T + (1|TestID), data =data.VARmin, REML = FALSE) 
mod.envT.CTmin.D3<-lmer(temp_tolerance ~ TimeDay2 + (1|TestID),data =data.VARmin, REML = FALSE)

mod.envT.CTmin.D4<-lmer(temp_tolerance ~ Temp_test_start + delta.T + (1|TestID), data =data.VARmin, REML = FALSE)
mod.envT.CTmin.D5<-lmer(temp_tolerance ~ Temp_test_start + TimeDay2 + (1|TestID), data =data.VARmin, REML = FALSE)
mod.envT.CTmin.D6<-lmer(temp_tolerance ~ delta.T + TimeDay2 + (1|TestID), data =data.VARmin, REML = FALSE) 
mod.envT.CTmin.D7<-lmer(temp_tolerance ~ Temp_test_start + delta.T + TimeDay2 + (1|TestID), data =data.VARmin, REML = FALSE) 

mod.envT.CTmin.D8<-lmer(temp_tolerance ~ min.Env.Temp + (1|TestID), data = data.VARmin, REML = FALSE)
mod.envT.CTmin.D9<-lmer(temp_tolerance ~ max.Env.Temp + (1|TestID), data = data.VARmin, REML = FALSE)
mod.envT.CTmin.D10<-lmer(temp_tolerance ~ mean.Env.Temp + (1|TestID), data = data.VARmin, REML = FALSE)
mod.envT.CTmin.D11<-lmer(temp_tolerance ~ mean.Env.Temp + delta.T + (1|TestID), data = data.VARmin, REML = FALSE)
mod.envT.CTmin.D12<-lmer(temp_tolerance ~ mean.Env.Temp + Temp_test_start + (1|TestID), data = data.VARmin, REML = FALSE)

BICdelta(BIC(mod.envT.CTmin.D1, mod.envT.CTmin.D2, mod.envT.CTmin.D3,
    mod.envT.CTmin.D4, mod.envT.CTmin.D5, mod.envT.CTmin.D6, mod.envT.CTmin.D7,
    mod.envT.CTmin.D8, mod.envT.CTmin.D9, mod.envT.CTmin.D10, 
    mod.envT.CTmin.D11, mod.envT.CTmin.D12)) # 


# CTmax -- All variability tests together, different explanatory variables:
mod.envT.CTmax.D1<-lmer(temp_tolerance ~ Temp_test_start + (1|TestID), data =data.VARmax, REML = FALSE)
mod.envT.CTmax.D2<-lmer(temp_tolerance ~ delta.T + (1|TestID), data =data.VARmax, REML = FALSE) 
mod.envT.CTmax.D3<-lmer(temp_tolerance ~ TimeDay2 + (1|TestID),data =data.VARmax, REML = FALSE)

mod.envT.CTmax.D4<-lmer(temp_tolerance ~ Temp_test_start + delta.T + (1|TestID), data =data.VARmax, REML = FALSE)
mod.envT.CTmax.D5<-lmer(temp_tolerance ~ Temp_test_start + TimeDay2 + (1|TestID), data =data.VARmax, REML = FALSE)
mod.envT.CTmax.D6<-lmer(temp_tolerance ~ delta.T + TimeDay2 + (1|TestID), data =data.VARmax, REML = FALSE) 
mod.envT.CTmax.D7<-lmer(temp_tolerance ~ Temp_test_start + delta.T + TimeDay2 + (1|TestID), data =data.VARmax, REML = FALSE) 

mod.envT.CTmax.D8<-lmer(temp_tolerance ~ min.Env.Temp + (1|TestID), data = data.VARmax, REML = FALSE)
mod.envT.CTmax.D9<-lmer(temp_tolerance ~ max.Env.Temp + (1|TestID), data = data.VARmax, REML = FALSE)
mod.envT.CTmax.D10<-lmer(temp_tolerance ~ mean.Env.Temp + (1|TestID), data = data.VARmax, REML = FALSE)
mod.envT.CTmax.D11<-lmer(temp_tolerance ~ mean.Env.Temp + delta.T + (1|TestID), data = data.VARmax, REML = FALSE)
mod.envT.CTmax.D12<-lmer(temp_tolerance ~ mean.Env.Temp + Temp_test_start + (1|TestID), data = data.VARmax, REML = FALSE)


BICdelta(BIC(mod.envT.CTmax.D1, mod.envT.CTmax.D2, mod.envT.CTmax.D3,
    mod.envT.CTmax.D4, mod.envT.CTmax.D5, mod.envT.CTmax.D6, mod.envT.CTmax.D7,
    mod.envT.CTmax.D8, mod.envT.CTmax.D9, mod.envT.CTmax.D10, 
    mod.envT.CTmax.D11, mod.envT.CTmax.D12)) # mod.envT.CTmax

# Stats 
summary(mod.envT.CTmax.D10)
summary(mod.envT.CTmin.D10)

car:::Anova(mod.envT.CTmin.D10, "II")
car:::Anova(mod.envT.CTmax.D10, "II")

pred.dataVAR<-as.data.frame(expand.grid( mean.Env.Temp = seq(12, 30, 2)))
pred.dataVAR$predCTmin<-predict(mod.envT.CTmin.D10, newdata = pred.dataVAR, re=NA)
pred.dataVAR$predCTmax<-predict(mod.envT.CTmax.D10, newdata = pred.dataVAR, re=NA)

## v). [suppl 3] size (mm) and CT -----
# field only
modF.min.scaling<-lmer(temp_tolerance ~  TL_mm + (1|TestID), data = data.minF)
modF.min.scaling0<-lmer(temp_tolerance ~ 1 + (1|TestID), data = data.minF)
modF.max.scaling<-lmer(temp_tolerance ~  TL_mm + (1|TestID), data = data.maxF)
modF.max.scaling0<-lmer(temp_tolerance ~  1 + (1|TestID), data = data.maxF)

modF.min.scaling<-lmer(temp_tolerance ~  TL_mm + (1|TestID), data = data.minF)
modF.min.scaling0<-lmer(temp_tolerance ~ 1 + (1|TestID), data = data.minF)
modF.max.scaling<-lmer(temp_tolerance ~  TL_mm + (1|TestID), data = data.maxF)
modF.max.scaling0<-lmer(temp_tolerance ~  1 + (1|TestID), data = data.maxF)

BIC(modF.min.scaling, modF.min.scaling0)
BIC(modF.max.scaling, modF.max.scaling0)
car::Anova(modF.max.scaling)
car::Anova(modF.min.scaling)


# Figures -------------
## [extra 3]: LAB - Growth mg  ----------
plot_growth <- ggplot()+
  geom_line(data = pred.datamg, mapping = aes(x = temp, y = pred.GR), color = "grey30", lwd = 0.9)+
  geom_line(data = pred.datamg, mapping = aes(x = temp, y = pred.GR+pred.GR.se), color = "grey30", lwd = 0.2, lty=2)+
  geom_line(data = pred.datamg, mapping = aes(x = temp, y = pred.GR-pred.GR.se), color = "grey30", lwd = 0.2, lty=2)+
  geom_line(data = growth.mod.data, aes(x = temp, GR_mg_d, color = temp_treatment.x), alpha=0.3)+
  geom_point(data = growth.mod.data, aes(x = temp, GR_mg_d, color = temp_treatment.x, fill = temp_treatment.x),
             alpha= 0.3, pch = 22, size=2)+
  geom_point(data = growth.mod.mean, aes(x = temp, mean_GR, color = temp_treatment.x, fill = temp_treatment.x),
             alpha= 1, pch = 22, size=3)+
  geom_errorbar(data = growth.mod.mean, mapping = aes(x = temp, ymin = mean_GR-sd_GR, ymax = mean_GR+sd_GR, color = temp_treatment.x),
                linewidth = 0.2, width = 0.1)+
  # geom_point(data2.sum.w[c(!(is.na(data2.sum.w$GR_mg_d)) & c(data2.sum.w$temp_treatment.x=="V1")),], 
  #            mapping = aes(x = 20, y = GR_mg_d), color = "#A3ABBD", fill = "#A3ABBD", pch = 22, size=2, alpha = 0.3)+
  # geom_point(data2.sum.w[c(!(is.na(data2.sum.w$GR_mg_d)) & c(data2.sum.w$temp_treatment.x=="V2")),], 
  #            mapping = aes(x = 11.8, y = GR_mg_d), color = "black", fill = "black",  pch = 22, size=2, alpha = 0.3)+
  geom_hline(yintercept = 0, linetype="solid", color="red", alpha=0.1)+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD") )+
  scale_fill_manual(values=c("12" ="#00518C",  "17" = "#008A60","22" ="#DBA11C","27" ="#BE647D", "V2" = "black", "V1" = "#A3ABBD") )+
  theme_classic()+
  # ylim(-15, 40)+
  xlim(10,30)
ggformat(plot_growth, x_title = "Temperature treatment", y_title = expression(Growth~rate~(mg~d^-1)), print=FALSE)
plot_growth <- plot_growth + theme (legend.position = "none")
# ggsave(plot_growth, filename = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Carpinteria marsh fish /DataAnalysis/PLOTS/Lab_growth_rate.png", width = 4, height = 4)

## Size all individuals --------
plot_K<- ggplot(data2, aes(temp_treatment, K_cond, colour=timepoint))+
  geom_point(alpha=1, position=position_dodge(width=0.5))+
  geom_boxplot(alpha=0.3, position=position_dodge(width=0.5))+
  scale_color_manual(values=c("grey40", "dodgerblue", "coral"))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = "dashed")+
  ylim(0, 5)
ggformat(plot_K, x_title = "Temperature treatment", y_title = expression(Condition~Factor~(K:~100*'%'*mg/TL^3)), print=FALSE)

plot_mg<- ggplot(data2, aes(temp_treatment, mass_mg, color = timepoint))+
  geom_point(alpha=1, position=position_dodge(width=0.5))+
  geom_boxplot(alpha=0.3, position=position_dodge(width=0.5))+
  scale_color_manual(values=c("grey40", "dodgerblue", "coral"))+
  geom_vline(xintercept = 4.5, linetype="dashed")+
  theme_classic()
ggformat(plot_mg, x_title = "Temperature treatment", y_title = "Mass (mg)", print=FALSE)
plot_mg<-plot_mg+theme(legend.position = "none") 

plot_mm<- ggplot(data2, aes(temp_treatment, TL_mm, colour=timepoint))+
  geom_point(alpha=1, position=position_dodge(width=0.5))+
  scale_color_manual(values=c("grey40", "dodgerblue", "coral"))+
  geom_boxplot(alpha=0.3, position=position_dodge(width=0.5))+
  geom_vline(xintercept = 4.5, linetype="dashed")+
  theme_classic()
ggformat(plot_mm, x_title = "Temperature treatment", y_title = "Total length (mm)", print=FALSE)
plot_mm<-plot_mm+theme(legend.position = "none") 

## 2. LAB: Figures CTmax: --------
plot_CTmax <- ggplot(data.max, aes(temp_treatment, temp_tolerance, fill = temp_treatment, color = temp_treatment, group = temp_treatment ))+
  geom_boxplot(alpha=0.5, width = 1)+
  geom_point(alpha=0.8, pch=21,  size=3)+
  geom_vline(xintercept = 4.5, linetype="dashed")+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD") )+
  scale_fill_manual(values=c("12" ="#00518C",  "17" = "#008A60","22" ="#DBA11C","27" ="#BE647D", "V2" = "black", "V1" = "#A3ABBD") )+
  theme_classic()+
  ylim(30, 42)
ggformat(plot_CTmax, x_title = "Temperature treatment", y_title = expression(CT[max]~(degree*C)), print=FALSE)
plot_CTmax<- plot_CTmax+theme(legend.position = "top")

plot_CTmin<- ggplot(data.min, aes(temp_treatment, temp_tolerance, fill = temp_treatment, color = temp_treatment, group = temp_treatment ))+
  geom_boxplot(alpha=0.5, width = 1)+
  geom_point(alpha=0.8, pch=21,  size=3)+
  geom_vline(xintercept = 4.5, linetype="dashed")+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD") )+
  scale_fill_manual(values=c("12" ="#00518C",  "17" = "#008A60","22" ="#DBA11C","27" ="#BE647D", "V2" = "black", "V1" = "#A3ABBD") )+
  theme_classic()+
  ylim(0, 12)
ggformat(plot_CTmin, x_title = "Temperature treatment", y_title = expression(CT[min]~(degree*C)), print=FALSE)
plot_CTmin<- plot_CTmin+theme(legend.position = "top")

# conditions factor --------
# 
plot_KCTtests<- ggplot(data, aes(TestID, K_cond, color = Date.1))+
  geom_hline(yintercept = c(0.5, 1), linetype = "dashed")+
  geom_point(alpha=1, position=position_dodge(width=0.5))+
  geom_boxplot(alpha=0.3, position=position_dodge(width=0.5))+
  geom_vline(xintercept = 4.5, linetype="dashed")+
  theme_classic()
ggformat(plot_KCTtests, x_title = "Temperature treatment", y_title = expression(Condition~Factor~(K:~100*'%'*mg/TL^3)), print=TRUE)
plot_KCTtests <- plot_KCTtests +theme(axis.text.x = element_text(angle = 45, size=6))

# Field correl, delta ----
# 
CTplotFIELD.delta<-ggplot(data=data, aes(y=temp_tolerance, x=delta.T, shape = Field_Lab, color = temp_treatment, fill = temp_treatment, label = TestID))+
  geom_point(size=2)+
  scale_shape_manual(values = c(1, 21))+
  # geom_errorbar( data=dat.ctF.l[dat.ctF.l$test == "CTmin",], mapping = aes(ymin = mean_CT - sd_CT, ymax = mean_CT + sd_CT),  width=0.2, size=0.5, alpha=0.8)+
  theme_classic()+
  # geom_text(check_overlap = T)+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD", "NA" = "darkgreen") )+
  scale_fill_manual(values=c("12" ="#00518C",  "17" = "#008A60","22" ="#DBA11C","27" ="#BE647D", "V2" = "black", "V1" = "#A3ABBD", "NA" = "darkgreen") )
ggformat(CTplotFIELD.delta, x_title = expression(Delta~Daily~Temperature~(degree*C)), y_title = expression(Temperature~tolerance~(degree*C)), print = TRUE)
# CTplotFIELD.delta<-CTplotFIELD.delta+theme(legend.position = c(0.2, 0.5), 
#                                              legend.title = element_text(), 
#                                              legend.background = element_rect(fill = "white", color = "black"))+
#   guides(fill=guide_legend(title=expression(Delta~daily~degree*C)))
#  *************************************

CTplotFIELD3<-ggplot(data=data, aes(y=temp_tolerance, x=min.Env.Temp, shape = Field_Lab,color = temp_treatment, fill = temp_treatment))+
  geom_point(size=2)+
  scale_shape_manual(values = c(1, 21))+
  theme_classic()+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD", "NA" = "darkgreen") )+
  scale_fill_manual(values=c("12" ="#00518C",  "17" = "#008A60","22" ="#DBA11C","27" ="#BE647D", "V2" = "black", "V1" = "#A3ABBD", "NA" = "darkgreen") )
ggformat(CTplotFIELD3, x_title = expression(Min~Daily~Temperature~(degree*C)), y_title = expression(Temperature~tolerance~(degree*C)), print = FALSE)



# misc
# 
# delta and variabitliy field 
PLOT1<-ggplot()+
  # geom_line(data=pred.data, mapping = aes(y=predCTmax, x=delta.T), lty=2, color = "grey", lwd = 0.5)+
  geom_abline(slope = fixef(mod.envT.CTmin.D1)[2], intercept = fixef(mod.envT.CTmin.D1)[1], lty=1, color = "grey30", lwd = 0.7)+
  geom_point(data=dat.V[dat.V$Test == "CTmax",], mapping = aes(x=delta.T, y=CT_degC, fill =mean.T, color = mean.T), size=2, pch=21, alpha = 0.8, show.legend = FALSE)+
  geom_errorbar( data=dat.sum.V[dat.sum.V$Test == "CTmax",], mapping = aes(x=delta.T, ymin = mean_CT - se_CT, ymax = mean_CT + se_CT), color = "grey30", width=0.1, size=0.4, alpha=0.8, show.legend = FALSE)+
  geom_point(data=dat.sum.V[dat.sum.V$Test == "CTmax",], mapping = aes(x=delta.T, y=mean_CT, fill =mean.T, label = n), color = "black", size=3, pch=21, show.legend = FALSE)+
  geom_point(data=dat.V[dat.V$Test == "CTmin",], mapping = aes(x=delta.T, y=CT_degC, fill =mean.T, color = mean.T), size=2, pch=23, alpha = 0.8, show.legend = F)+
  geom_errorbar( data=dat.sum.V[dat.sum.V$Test == "CTmin",], mapping = aes(x=delta.T, ymin = mean_CT - se_CT, ymax = mean_CT + se_CT), color = "grey30", width=0.1, size=0.4, alpha=0.8, show.legend = FALSE)+
  geom_point(data=dat.sum.V[dat.sum.V$Test == "CTmin",], mapping = aes(x=delta.T, y=mean_CT, fill =mean.T, label = n), color = "black", size=3, pch=23, show.legend = FALSE)+
  scale_color_gradient2(low = "#5302a3", high = "#f0f921", mid = "#f48849", midpoint = 22 , guide = "none")+
  scale_fill_gradient2(low = "#5302a3", high = "#f0f921", mid = "#f48849", midpoint = 22 , guide = "none")+
 
  geom_point(dat.V[ c(dat.V$treatm=="V1") & c(dat.V$Test=="CTmax"),], mapping = aes(x=delta.T, y=CT_degC,), color = "black", fill = "black", alpha=1, pch = 21, size=2, show.legend = FALSE)+
  geom_point(dat.V[ c(dat.V$treatm=="V2") & c(dat.V$Test=="CTmax"),], mapping = aes(x=delta.T, y=CT_degC,), color = "#A3ABBD", fill = "#A3ABBD", alpha=1, pch = 21, size=2, show.legend = FALSE)+
  geom_point(dat.V[ c(dat.V$treatm=="V1") & c(dat.V$Test=="CTmin"),], mapping = aes(x=delta.T, y=CT_degC,), color = "black", fill = "black", alpha=1, pch = 23, size=2, show.legend = FALSE)+
  geom_point(dat.V[ c(dat.V$treatm=="V2") & c(dat.V$Test=="CTmin"),], mapping = aes(x=delta.T, y=CT_degC,), color = "#A3ABBD", fill = "#A3ABBD", alpha=1, pch = 23, size=2, show.legend = FALSE)+
  geom_point(data=dat.sum.V[c(dat.sum.V$treatm=="V1" & dat.sum.V$Test == "CTmax"),], mapping = aes(x=delta.T, y=mean_CT), color = "black", fill = "black", size=3, pch=21, show.legend = FALSE)+
  geom_point(data=dat.sum.V[c(dat.sum.V$treatm=="V2" & dat.sum.V$Test == "CTmax"),], mapping = aes(x=delta.T, y=mean_CT), color = "black", fill = "#A3ABBD", size=3, pch=21, show.legend = FALSE)+
  geom_point(data=dat.sum.V[c(dat.sum.V$treatm=="V1" & dat.sum.V$Test == "CTmin"),], mapping = aes(x=delta.T, y=mean_CT), color = "black", fill = "black", size=3, pch=23, show.legend = FALSE)+
  geom_point(data=dat.sum.V[c(dat.sum.V$treatm=="V2" & dat.sum.V$Test == "CTmin"),], mapping = aes(x=delta.T, y=mean_CT), color = "black", fill = "#A3ABBD", size=3, pch=23, show.legend = FALSE)+
  theme_classic()
ggformat(PLOT1, x_title = expression(Delta~Daily~Temperature~(degree*C)), y_title = expression(Temperature~tolerance~(degree*C)), print = FALSE)
PLOT1<-PLOT1 +theme(legend.position = c(0.2, 0.5), 
                legend.title = element_text(), 
                legend.background = element_rect(fill = "white", size = 0.1, color = "black"))+
  guides(fill=guide_legend(title=expression(Mean~daily~degree*C)))
PLOT1


PLOT2<-ggplot()+
  # geom_line(data=pred.data, mapping = aes(y=predCTmax, x=Start.T), lty=2, color = "grey", lwd = 0.5)+
  # geom_point(data=dat[dat$Test == "CTmax",], mapping = aes(x=Start.T, y=CT_degC, fill =delta.T, color = delta.T), size=2, pch=21, alpha = 0.8, show.legend = FALSE)+
  geom_errorbar( data=data.sum[dat.sum$Test == "CTmax",], mapping = aes(x=Start.T, ymin = mean_CT - se_CT, ymax = mean_CT + se_CT), color = "grey30", width=0.1, size=0.4, alpha=0.8, show.legend = FALSE)+
  geom_point(data=dat.sum[dat.sum$Test == "CTmax",], mapping = aes(x=Start.T, y=mean_CT, fill =Locat, label = n), color = "black", size=3, pch=21, show.legend = FALSE)+
  # geom_point(data=dat[dat$Test == "CTmin",], mapping = aes(x=Start.T, y=CT_degC, fill =delta.T, color = delta.T), size=2, pch=23, alpha = 0.8, show.legend = FALSE)+
  geom_errorbar( data=dat.sum[dat.sum$Test == "CTmin",], mapping = aes(x=Start.T, ymin = mean_CT - se_CT, ymax = mean_CT + se_CT), color = "grey30", width=0.1, size=0.4, alpha=0.8, show.legend = FALSE)+
  geom_point(data=dat.sum[dat.sum$Test == "CTmin",], mapping = aes(x=Start.T, y=mean_CT, fill =Locat, label = n), color = "black", size=3, pch=23, show.legend = FALSE)+
  scale_fill_manual(values = c("#EA573D", "white"))+
  
  # geom_point(dat[ c(dat$treatm=="V1") & c(dat$Test=="CTmax"),], mapping = aes(x=Start.T, y=CT_degC,), color = "black", fill = "black", alpha=1, pch = 21, size=2, show.legend = FALSE)+
  # geom_point(dat[ c(dat$treatm=="V2") & c(dat$Test=="CTmax"),], mapping = aes(x=Start.T, y=CT_degC,), color = "#A3ABBD", fill = "#A3ABBD", alpha=1, pch = 21, size=2, show.legend = FALSE)+
  # geom_point(dat[ c(dat$treatm=="V1") & c(dat$Test=="CTmin"),], mapping = aes(x=Start.T, y=CT_degC,), color = "black", fill = "black", alpha=1, pch = 23, size=2, show.legend = FALSE)+
  # geom_point(dat[ c(dat$treatm=="V2") & c(dat$Test=="CTmin"),], mapping = aes(x=Start.T, y=CT_degC,), color = "#A3ABBD", fill = "#A3ABBD", alpha=1, pch = 23, size=2, show.legend = FALSE)+
  geom_point(data=dat.sum[c(dat.sum$treatm=="V1" & dat.sum$Test == "CTmax"),], mapping = aes(x=Start.T, y=mean_CT), color = "black", fill = "black", size=3, pch=21, show.legend = FALSE, stroke=1)+
  geom_point(data=dat.sum[c(dat.sum$treatm=="V2" & dat.sum$Test == "CTmax"),], mapping = aes(x=Start.T, y=mean_CT), color = "black", fill = "#A3ABBD", size=3, pch=21, show.legend = FALSE, stroke=1)+
  geom_point(data=dat.sum[c(dat.sum$treatm=="V1" & dat.sum$Test == "CTmin"),], mapping = aes(x=Start.T, y=mean_CT), color = "black", fill = "black", size=3, pch=23, show.legend = FALSE, stroke=1)+
  geom_point(data=dat.sum[c(dat.sum$treatm=="V2" & dat.sum$Test == "CTmin"),], mapping = aes(x=Start.T, y=mean_CT), color = "black", fill = "#A3ABBD", size=3, pch=23, show.legend = FALSE, stroke=1)+
  theme_classic()
ggformat(PLOT2, x_title = expression(Test~Start~Temperature~(degree*C)), y_title = expression(Temperature~tolerance~(degree*C)), print = FALSE)
PLOT2<-PLOT2 +theme(legend.position = c(0.14, 0.5), 
                    legend.title = element_text(), 
                    legend.background = element_rect(fill = "white", size = 0.1, color = "black"))+
  guides(fill=guide_legend(title=expression(Delta~Field~degree*C)))
PLOT2


cowplot:::plot_grid(PLOT1, PLOT2,
                    labels = "AUTO", 
                    nrow =1,
                    ncol=2,
                    align = "hv", 
                    label_x = c(0.3, 0.22),
                    label_y = c(0.9)) %>% 
  ggsave( filename = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Carpinteria marsh fish /DataAnalysis/PLOTS/Figure3_CT_VAR.png", width = 8, height = 4)
# ggsave(plot_CTsizeF, filename = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Carpinteria marsh fish /DataAnalysis/PLOTS/Field_CTmaxCT_size.png", width = 7, height = 3.5)



# saving plots ############# 
setwd("/Users/kristakraskura/Desktop/BOX/UCSB/Research/Carpinteria - marsh fish /DataAnalysis/PLOTS/")

png("Plot1_mass_weight_3tp_Dec4.png", width=6, height = 8 , res=300, units = "in")
grid.arrange(plot_mg, plot_mm, nrow=2, ncol=1, heights=c(1,0.9))
dev.off()

png("Plot2_temp_toler_Dec4.png", width= 8, height = 4, res=300, units = "in")
grid.arrange(plot_CTmax, plot_CT, nrow=1, ncol=2, widths = c(1,0.95))
dev.off()

png("plot3_cond_fact_Dec4.png", width=5, height = 3.5, res=300, units = "in")
grid.arrange(plot_K)
dev.off()

png("Plot4_CTall_Dec4.png", width= 5, height = 3.5, res=300, units = "in")
print(CTplot)
dev.off()

png("Plot5_CTsize_Dec4.png", width=7 , height =4 , res=300, units = "in")
print(plot_CTsize)
dev.off()



# merge
CTplot<-
  ggplot()+
  geom_line(data = pred.data, mapping = aes(y = ctmax_pred, x = as.numeric(temp)),
            alpha=1, colour="grey30", lty=1)+
  geom_line(data = pred.data, mapping = aes(y = CT_pred, x = as.numeric(temp)),
            alpha=1, colour="grey30", lty=1)+
  # geom_line(data = pred.data, mapping = aes(y = CT_predERmax, x = as.numeric(temp)),
  #           alpha=1, colour="grey30", lwd=0.2, lty="dashed")+
  # geom_line(data = pred.data, mapping = aes(y = CT_predERmin, x = as.numeric(temp)),
  #           alpha=1, colour="grey30", lwd=0.2, lty="dashed")+
  # geom_line(data = pred.data, mapping = aes(y = ctmax_predERmax, x = as.numeric(temp)),
  #           alpha=1, colour="grey30", lwd=0.2, lty="dashed")+
  # geom_line(data = pred.data, mapping = aes(y = ctmax_predERmin, x = as.numeric(temp)),
  #           alpha=1, colour="grey30", lwd=0.2, lty="dashed" )+
  geom_point(data=data.static.min, aes(temp, temp_tolerance,  color = treatm, fill = treatm), pch=23, size=2, alpha=1, show.legend = FALSE)+
  geom_point(data=data.static.max, aes(temp, temp_tolerance,  color = treatm, fill = treatm), pch=21, size=2, alpha=1, show.legend = FALSE)+
  geom_point(data[c(!(is.na(data$temp_tolerance)) & c(data$temp=="V1") & c(data$Test=="CTmax")),], mapping = aes(x = inters.V1ctmax.x, y = temp_tolerance), color = "black", fill = "black", alpha=1, pch = 21, size=2, show.legend = FALSE)+
  geom_point(data[c(!(is.na(data$temp_tolerance)) & c(data$temp=="V2") & c(data$Test=="CTmax")),], mapping = aes(x = inters.V2ctmax.x, y = temp_tolerance), color = "#A3ABBD", fill = "#A3ABBD", alpha=1, pch = 23, size=2, show.legend = FALSE)+
  geom_point(data[c(!(is.na(data$temp_tolerance)) & c(data$temp=="V1") & c(data$Test=="CTmin")),], mapping = aes(x = inters.V1CT.x, y = temp_tolerance), color = "black", fill = "black", alpha=1, pch = 21, size=2, show.legend = FALSE)+
  geom_point(data[c(!(is.na(data$temp_tolerance)) & c(data$temp=="V2") & c(data$Test=="CTmin")),], mapping = aes(x = inters.V2CT.x, y = temp_tolerance), color = "#A3ABBD", fill = "#A3ABBD", alpha=1, pch = 23, size=2, show.legend = FALSE)+
  geom_point(data=datsum, aes(as.numeric(temp), mean_CT, group = treatment,  fill = treatment), color = "black", size=3, alpha=1, pch=23)+
  geom_point(data=datsum, aes(as.numeric(temp), mean_CTmax, group = treatment,  fill = treatment), color = "black", size=3, alpha=1, pch=21)+
  geom_point(mapping=aes(x = inters.V1CT.x, y=inters.V1CT.y), pch=23,  size=4, alpha=1, color = "black", fill = "#A3ABBD")+
  geom_point(mapping=aes(x = inters.V1ctmax.x, y=inters.V1ctmax.y), pch=21,  size=4, alpha=1, color = "black", fill = "#A3ABBD")+
  geom_point(mapping=aes(x = inters.V2CT.x, y=inters.V2CT.y), pch=23,  size=4, alpha=1, color = "black", fill = "black")+
  geom_point(mapping=aes(x = inters.V2ctmax.x, y=inters.V2ctmax.y), pch=21,  size=4, alpha=1, color = "black", fill = "black")+# geom_point(data[c(!(is.na(data$temp_tolerance)) & c(data$temp=="V1")),],
  ylim(0, 45)+
  xlim(10, 30)+
  scale_fill_manual(values=c("COLD" = "#00518C", "HOT" = "#BE647D",  "MED" = "#DBA11C", "LOW" ="#008A60",  "V1" = "#A3ABBD", "V2" =  "black") )+
  scale_color_manual(values=c("COLD" = "#00518C", "HOT" = "#BE647D",  "MED" = "#DBA11C", "LOW" ="#008A60",  "V1" = "#A3ABBD", "V2" =  "black") )+
  theme_classic()
  # geom_line(aes(temp, mean_CTmax, group = timepoint))
  # annotate("text", x=1, y=44, colour="darkred", size=4, label="CTmax", hjust=0)+
  # annotate("text", x=1, y=41, colour="dodgerblue", size=4, label="CTmin", hjust=0)
ggformat(CTplot, x_title = "Temperature treatment", y_title = expression(Temperature~tolerance~(degree*C)), print = TRUE)
CTplot <- CTplot + theme(legend.position = "none")


# where on the regresion line is ct min?
# y = int+x(m)
# -xm = int - y
# xm = y - int
# x = (y - int) /m 

# inters.V1CT.x<-(dat.CTmin$mean_CT[5] - fixef(mod1min)[1])/ fixef(mod1min)[2]
# inters.V1CT.y<-(dat.CTmin$mean_CT[5])
# inters.V1ctmax.x<-(dat.CTmax$mean_CTmax[5] - fixef(mod1)[1])/ fixef(mod1)[2]
# inters.V1ctmax.y<-(dat.CTmax$mean_CTmax[5])
# 
# inters.V2CT.x<-(dat.CTmin$mean_CT[6] - fixef(mod1min)[1])/ fixef(mod1min)[2]
# inters.V2CT.y<-(dat.CTmin$mean_CT[6])
# inters.V2ctmax.x<-(dat.CTmax$mean_CTmax[6] - fixef(mod1)[1])/ fixef(mod1)[2]
# inters.V2ctmax.y<-(dat.CTmax$mean_CTmax[6])
