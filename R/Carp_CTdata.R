
library(ggformat2) # from kraskura/ggformat github package. 
library(chron)
library(lubridate)
library(lattice) # for qqmath
library(zoo)
library(ggrepel)
library(lme4) # linear models
library(emmeans) # post hocs
library(merTools)# for predicting CI intervals
# library(lmerTest) 
library(tidyverse)

library(here)

# function used to calculate ∆BIC scores and re-order based on the lowest score. 
BICdelta<-function(BICtable){
  BIC.t <- BICtable [order(BICtable$BIC), ]
  BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
  return( BIC.t)
}

# set colors 
cols.klab<- c("#00518C", "#008A60", "#DBA11C", "#BE647D", "#A3ABBD","black", "#00CAFF")  # 27, 22, 17, 12, V1, V2, "field"
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

# lab tests
data.min<-data[c(data$Field_Lab == "LAB" & data$Test=="CTmin" ),] # all LAB CTmin
data.max<-data[c(data$Field_Lab == "LAB" & data$Test=="CTmax"),] # all LAB CTmax
data.static.max<-data[c(data$Field_Lab == "LAB" & data$Test=="CTmax" & data$Stat_Var == "static"),] # used for CT min and max stats 
data.static.min<-data[c(data$Field_Lab == "LAB" & data$Test=="CTmin" & data$Stat_Var == "static"),] # used for CT min and max stats 

# variables 
data.VARmin<-data[data$Stat_Var == "variable" & data$Test == "CTmin",  ]
data.VARmax<-data[data$Stat_Var == "variable" & data$Test == "CTmax",  ]

#field tests:
data.minF<-data[c(data$Field_Lab == "FIELD" & data$Test=="CTmin"),] # all Field CTmax
data.maxF<-data[c(data$Field_Lab == "FIELD" & data$Test=="CTmax"),] # all Field CTmin
summary(data.minF)
summary(data.maxF)

# import field temps to get daily min and max temps.
fieldTempsJul<- read.csv("./Data/Analysis source/SiteTempsJul.csv")
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
  
fieldTempsSep<- read.csv("./Data/Analysis source/SiteTempsSep.csv")
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
            mean_mm=mean(TL_mm), sd_mm=sd(TL_mm), min_mm=min(TL_mm), max_mm=max(TL_mm), n=length(!is.na(mass_mg))) %>% 
  as.data.frame()

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


# STATS:----- 

boxplot(temp_tolerance ~ TestID2, data = data[data$Test == "CTmin",])
boxplot(temp_tolerance ~ TestID2, data = data[data$Test == "CTmax",])

## 1. [suppl 2] LAB: growth ----- 
# mass in mg <<!  this is less reliable because fish are so small and the water residual error can overpower the mass. 
# total length in mm 
mod.poly.growthmm <- lm(GR_mm_d ~ poly(temp,2), data = growth.mod.data)
predict(mod.poly.growthmm)
summary(mod.poly.growthmm)

pred.datamm<-as.data.frame(expand.grid(temp =seq(12, 27, 0.2)))
pred.datamm$pred.GRmm<-predict(mod.poly.growthmm, newdata = pred.datamm)
pred.datamm$pred.GRmm.se<-predict(mod.poly.growthmm, newdata = pred.datamm, se.fit = TRUE)[[2]]

## 2. [suppl 1a] LAB & STATIC ONLY  (continuous temp predictor) --------
# H: CTmin and CTmax increase with increasing static acclimation temperatures.
# ***********************************************
# mixed model using tank as a random effect
# data<-data[data$Phys.cond =="NOTFED", ] # only fish that were not fed. 
# dataFED<-data[data$Phys.cond =="FED", ] 

# fish that were FED are included in this run. 
mod1<-lmer(temp_tolerance ~  temp + (1|tank), REML = FALSE, data = data.static.max)
mod1.b<-lmer(temp_tolerance ~  temp + TL_mm + (1|tank), REML = FALSE, data = data.static.max)

BIC(mod1, mod1.b) # with mass better 
qqmath(mod1.b)
summary(mod1.b)
# plot(resid(mod1.b))

mod1min<-lmer(temp_tolerance ~  temp + (1|tank), REML = FALSE, data = data.static.min) # singular fit for the random effects
mod1min.b<-lmer(temp_tolerance ~  temp + TL_mm + (1|tank), REML = FALSE, data = data.static.min) # singular fit for the random effects
# plot(mod1min)
qqmath(mod1min.b)
summary(mod1min.b)
# plot(resid(mod1min.b))
# singular fit is effectively the same as simple linear regression (lm)

# type II anovas:  REPORTED
car::Anova(mod1.b, type = "II") #CTmax
car::Anova(mod1min.b, "II") # CTmin

# for predicting, plotting the data
pred.dataCT<-as.data.frame(expand.grid(temp =seq(12, 27, 0.2), TL_mm = 36))
pred.dataCT$ctmax_pred<-predict(mod1.b,newdata = pred.dataCT, re=NA)
pred.dataCT$ctmin_pred<-predict(mod1min.b,newdata = pred.dataCT, re=NA)

predCTmaxStatic<-predictInterval(mod1.b, data.static.max, which = "fixed",include.resid.var = TRUE, level = 0.95, n.sims = 1000, type="linear.prediction")
predCTminStatic<-predictInterval(mod1min.b, data.static.min, which = "fixed",include.resid.var = TRUE, level = 0.95, n.sims = 1000, type="linear.prediction")
predCTmaxStatic$temp<-data.static.max$temp
predCTminStatic$temp<-data.static.min$temp

## 3. [suppl 1b] ALL TESTS: ANOVA between all tests -----------
# H: CTmin and CTmax differ between treatments (lab) and tests (field)
# ***********************************************
# testID2 - treatments, or days tested. 
mod.CTmin<-lm(temp_tolerance ~ treatment, data = data[data$Test== "CTmin", ])
mod.CTmin.b<-lm(temp_tolerance ~ TL_mm + treatment, data = data[data$Test== "CTmin", ])
BICdelta(BIC(mod.CTmin, mod.CTmin.b))
car::Anova(mod.CTmin.b, type = "II")
# plot(mod.CTmin.b)
  
mod.CTmax<-lm(temp_tolerance ~ treatment, data = data[data$Test== "CTmax", ])
mod.CTmax.b<-lm(temp_tolerance ~ TL_mm + treatment, data = data[data$Test== "CTmax", ])
BICdelta(BIC(mod.CTmax, mod.CTmax.b))
car::Anova(mod.CTmax.b, type = "II")
# plot(mod.CTmax.b)

contrast(emmeans(mod.CTmax.b, ~ treatment), "pairwise")
contrast(emmeans(mod.CTmin.b, ~ treatment), "pairwise")

## 4. [suppl 1b] FIELD vs LAB: variable treatments -----------
# H: CTmin and CTmax is not different between lab and field variable treatments
# ***********************************************
mod.var.CTmin<-lmer(temp_tolerance ~ Field_Lab + (1|treatment), data = data.VARmin, REML = FALSE)
mod.var.CTmin.b<-lmer(temp_tolerance ~ Field_Lab + TL_mm + (1|treatment), data = data.VARmin, REML = FALSE)
BICdelta(BIC(mod.var.CTmin, mod.var.CTmin.b))
car::Anova(mod.var.CTmin.b, type = "II")

mod.var.CTmax<-lmer(temp_tolerance ~ Field_Lab + (1|treatment), data = data.VARmax, REML = FALSE)
mod.var.CTmax.b<-lmer(temp_tolerance ~ Field_Lab + TL_mm + (1|treatment), data = data.VARmax, REML = FALSE)
BICdelta(BIC(mod.var.CTmax, mod.var.CTmax.b))
car::Anova(mod.var.CTmax.b, type = "II")

## 5. [main 2]  ALL TESTS: continuous , static and variable.  max, mean, min, delta, start temp, start daytime -----
# H: [exploratory] what recent thermal history variable explains CTmin and CT max best?
# CT min tests 
model.CTmin.0<-lmer(temp_tolerance ~ 1  + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE)
model.CTmin.1<-lmer(temp_tolerance ~ max.Env.Temp + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.1.b<-lmer(temp_tolerance ~ max.Env.Temp + TL_mm + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.2<-lmer(temp_tolerance ~ min.Env.Temp + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.2.b<-lmer(temp_tolerance ~ min.Env.Temp + TL_mm + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.3<-lmer(temp_tolerance ~ mean.Env.Temp + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.3.b<-lmer(temp_tolerance ~ mean.Env.Temp +TL_mm +  (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.4<-lmer(temp_tolerance ~ Temp_test_start + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.4.b<-lmer(temp_tolerance ~ Temp_test_start + TL_mm + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.5<-lmer(temp_tolerance ~ delta.T + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.5.b<-lmer(temp_tolerance ~ delta.T + TL_mm + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.6<-lmer(temp_tolerance ~ TimeDay2 + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.6.b<-lmer(temp_tolerance ~ TimeDay2 + TL_mm + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 

model.CTmin.7<-lmer(temp_tolerance ~ max.Env.Temp + delta.T + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.7.c<-lmer(temp_tolerance ~ max.Env.Temp + delta.T + Stat_Var + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.7.b<-lmer(temp_tolerance ~ max.Env.Temp + delta.T + TL_mm + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) # best 
model.CTmin.7.d<-lmer(temp_tolerance ~ max.Env.Temp + delta.T + TL_mm + Stat_Var + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) # best 

model.CTmin.8<-lmer(temp_tolerance ~ max.Env.Temp + Temp_test_start + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.9<-lmer(temp_tolerance ~ delta.T + Temp_test_start + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.10<-lmer(temp_tolerance ~ max.Env.Temp + delta.T + TimeDay2 + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE) 
model.CTmin.11<-lmer(temp_tolerance ~ max.Env.Temp + Temp_test_start + delta.T + (1|treatment), data = data[data$Test== "CTmin", ], REML = FALSE)

### [suppl TABL 2] ----------
# very convincing that size matters, add as a covariate also below 
BICdelta(BIC(model.CTmin.0, model.CTmin.1, model.CTmin.2, model.CTmin.3, model.CTmin.4, model.CTmin.5, model.CTmin.6,
             model.CTmin.1.b, model.CTmin.2.b, model.CTmin.3.b, model.CTmin.4.b, model.CTmin.5.b, model.CTmin.6.b,
             model.CTmin.7,model.CTmin.7.d, model.CTmin.7.b, model.CTmin.7.c,
             model.CTmin.8, model.CTmin.9, model.CTmin.10, model.CTmin.11))


### [suppl TABL 2] ----------
car::Anova(model.CTmin.1.b, type = "II") #CTmax

# CT max tests 
model.CTmax.0<-lmer(temp_tolerance ~ 1  + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE)
model.CTmax.1<-lmer(temp_tolerance ~ max.Env.Temp + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.1.b<-lmer(temp_tolerance ~ max.Env.Temp + TL_mm + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.2<-lmer(temp_tolerance ~ min.Env.Temp + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.2.b<-lmer(temp_tolerance ~ min.Env.Temp + TL_mm + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 

model.CTmax.3<-lmer(temp_tolerance ~ mean.Env.Temp + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.3.b<-lmer(temp_tolerance ~ mean.Env.Temp + TL_mm + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) # <<< BEST
model.CTmax.3.c<-lmer(temp_tolerance ~ mean.Env.Temp + TL_mm + delta.T + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 

model.CTmax.4<-lmer(temp_tolerance ~ Temp_test_start + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.4.b<-lmer(temp_tolerance ~ Temp_test_start + TL_mm + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.5<-lmer(temp_tolerance ~ delta.T + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.5.b<-lmer(temp_tolerance ~ delta.T + TL_mm + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.6<-lmer(temp_tolerance ~ TimeDay2 + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.6.b<-lmer(temp_tolerance ~ TimeDay2 + TL_mm + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 

model.CTmax.7<-lmer(temp_tolerance ~ max.Env.Temp + delta.T + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.7.b<-lmer(temp_tolerance ~ max.Env.Temp + delta.T + TL_mm + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 

model.CTmax.8<-lmer(temp_tolerance ~ max.Env.Temp + Temp_test_start + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.9<-lmer(temp_tolerance ~ delta.T + Temp_test_start + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.10<-lmer(temp_tolerance ~ max.Env.Temp + delta.T + TimeDay2 + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE) 
model.CTmax.11<-lmer(temp_tolerance ~ max.Env.Temp + Temp_test_start + delta.T + (1|treatment), data = data[data$Test== "CTmax", ], REML = FALSE)

### [suppl TABL 2] ----------
# very convincing that size matters, add as a covariate also below 
BICdelta(BIC(model.CTmax.0, model.CTmax.1, model.CTmax.2, model.CTmax.3, model.CTmax.4, model.CTmax.5, model.CTmax.6,
             model.CTmax.1.b, model.CTmax.2.b, model.CTmax.3.b,model.CTmax.3.c, model.CTmax.4.b, model.CTmax.5.b, model.CTmax.6.b,
             model.CTmax.7,model.CTmax.7.b, 
             model.CTmax.8, model.CTmax.9, model.CTmax.10, model.CTmax.11))
### [suppl TABL 2] ----------
car::Anova(model.CTmax.3.b, type = "II") #CTmax


# model residuals, performance, etc. 
plot(model.CTmax.3.b) # good 
qqmath(model.CTmax.3.b)
hist(resid(model.CTmax.3.b), breaks = 50) # normal, good; few outliers

plot(model.CTmin.1.b) # good 
qqmath(model.CTmin.1.b)
hist(resid(model.CTmin.1.b), breaks = 50) # normal, good; few outliers

### variables for plotting -----
# mean env temp
CTmax.int<-round(fixef(model.CTmax.3.b)[1], 2) 
CTmax.slopeMean<-round(fixef(model.CTmax.3.b)[2], 2)
CTmax.slopeTL<-round(fixef(model.CTmax.3.b)[3], 2)
CTmax.n<-unlist(summary(model.CTmax.3.b)[[3]][2])[1]

# max env temp & delta 
CTmin.int<-round(fixef(model.CTmin.1.b)[1], 2)
CTmin.slopeMax<-round(fixef(model.CTmin.1.b)[2], 2)
CTmin.slopeTL<-round(fixef(model.CTmin.1.b)[3], 2)
CTmin.n<-unlist(summary(model.CTmin.1.b)[[3]][2])[1]

data.ctmin.predict<-as.data.frame(expand.grid(max.Env.Temp = c(12, 34), delta.T = c(7.844), TL_mm = 30.82))
data.ctmin.predict$pred.ctmin<- predict(model.CTmin.7.b, newdata =data.ctmin.predict, re.form = NA)

data.ctmax.predict<-as.data.frame(expand.grid(mean.Env.Temp = c(12, 34), delta.T = c(7.844), TL_mm = 30.82))
data.ctmax.predict$pred.ctmax<- predict(model.CTmax.3.b, newdata =data.ctmax.predict, re.form = NA)

# Data Summaries, stats: --------
###  get one representative value for a treatment 
# a dataframe with one mean per temp treatment 

dat<-data2  %>% 
  filter(timepoint == "1" | timepoint == "2") %>% 
  group_by(treatment, temp_treatment, timepoint) %>% 
  summarise(mean_mg=mean(mass_mg), sd_mg=sd(mass_mg), min_mg=min(mass_mg), max_mg=max(mass_mg),
            mean_mm=mean(TL_mm), sd_mm=sd(TL_mm), min_mm=min(TL_mm), max_mm=max(TL_mm), n=length(mass_mg))

dat.CTmin<-data %>% 
  filter(Test=="CTmin") %>% 
  # filter(Field_Lab=="FIELD") %>% 
  # filter(Stat_Var=="variable") %>%
  group_by(treatment, temp, Phys.cond, Field_Lab, TestID) %>% 
  summarise(mean_CT=mean(temp_tolerance),
            sd_CT=sd(temp_tolerance),
            min_CT=min(temp_tolerance),
            max_CT=max(temp_tolerance),
            n=length(temp_tolerance), 
            dailyMin = mean(min.Env.Temp), 
            dailyMax = mean(max.Env.Temp), 
            dailyMean = mean(mean.Env.Temp), 
            dailyRange = mean(mean.Env.Temp), 
            startTemp = mean(Temp_test_start)) %>% 
  as.data.frame()

dat.CTmax<-data %>% 
  filter(Test=="CTmax") %>% 
  # filter(Field_Lab=="FIELD") %>% 
  # filter(Stat_Var=="variable") %>%
  group_by(treatment, temp, Phys.cond,  Field_Lab, TestID) %>% 
  summarise(mean_CT=mean(temp_tolerance),
            sd_CT=sd(temp_tolerance),
            min_CT=min(temp_tolerance),
            max_CT=max(temp_tolerance),
            n=length(temp_tolerance),
            dailyMin = mean(min.Env.Temp), 
            dailyMax = mean(max.Env.Temp), 
            dailyMean = mean(mean.Env.Temp), 
            dailyRange = mean(mean.Env.Temp), 
            startTemp = mean(Temp_test_start)) %>% 
  as.data.frame()
  
data %>% 
  filter(Test=="CTmax") %>% 
  filter(Field_Lab=="LAB") %>% 
  filter(Stat_Var=="variable") %>% 
  summarise(mean = mean(temp_tolerance),
            min = min(temp_tolerance),
            max = max(temp_tolerance),
            CV = sd(temp_tolerance)/mean,
            n = n())

data %>% 
  filter(Test=="CTmin") %>% 
  filter(Field_Lab=="LAB") %>% 
  filter(Stat_Var=="variable") %>% 
  summarise(mean = mean(temp_tolerance),
            min = min(temp_tolerance),
            max = max(temp_tolerance),
            CV = sd(temp_tolerance)/mean,
            n = n())

dat.ct<-rbind(dat.CTmin, dat.CTmax)

data.ctmax<-data[data$Test == "CTmax", ]
data.ctmin<-data[data$Test == "CTmin", ]

summary(data.ctmax)
summary(data.ctmin)

length(levels(factor(data.ctmax$TestID)))
length(levels(factor(data.ctmin$TestID)))

n.tests.max.F<-length(levels(factor(data.maxF$TestID))) # field only 
n.tests.min.F<-length(levels(factor(data.minF$TestID))) # field only 
n.tests.max.L<-length(levels(factor(data.max$TestID))) # lab only 
n.tests.min.L<-length(levels(factor(data.min$TestID))) # lab only 

n.indiv.max.F<-nrow(data.maxF) # field only 
n.indiv.min.F<-nrow(data.minF) # field only 
n.indiv.max.L<-nrow(data.max) # lab only 
n.indiv.min.L<-nrow(data.min) # lab only 

# number of individuals timepoint 1:
sum(data2.sum.all[data2.sum.all$timepoint==1, "n" ])
# n field:
n.indiv.max.F + n.indiv.min.F

# CTmin and CTmax by static treatment
dat.CTmin.lab<-data.min %>% 
  filter(Stat_Var=="static") %>%  # comment out to have lab variable treatments here
  group_by(treatment) %>% 
  summarise(mean_CT=mean(temp_tolerance),
            sd_CT=sd(temp_tolerance),
            min_CT=min(temp_tolerance),
            max_CT=max(temp_tolerance),
            n=length(temp_tolerance), 
            diff_treatment = mean(temp - mean_CT), 
            CV = sd_CT / mean_CT) %>% 
  arrange(desc(-mean_CT)) %>% 
  mutate(diff_CT=mean_CT-lag(mean_CT))

dat.CTmax.lab<-data.max %>% 
  filter(Stat_Var=="static") %>%  # comment out to have lab variable treatments here
  group_by(treatment) %>% 
  summarise(mean_CT=mean(temp_tolerance),
            sd_CT=sd(temp_tolerance),
            min_CT=min(temp_tolerance),
            max_CT=max(temp_tolerance),
            n=length(temp_tolerance),
            diff_treatment = mean(mean_CT- temp),
            CV = sd_CT / mean_CT) %>% 
  arrange(desc(-mean_CT)) %>% 
  mutate(diff_CT=mean_CT-lag(mean_CT))


# FIGURES:  ----------
## 3. [main 1] Figures main, all CTs with correlations of env. temp-----
## [main 1a] start temp *************************************
CTplotFIELD<-ggplot(data=data, aes(y=temp_tolerance, x=Temp_test_start, shape = Stat_Var, alpha = Stat_Var, size = Stat_Var, fill = temp_treatment, color = temp_treatment))+
  geom_point()+
  scale_shape_manual(values = c(5, 21), name = "")+
  scale_size_manual(values = c(2, 3))+
  scale_alpha_manual(values = c(1, 0.5))+
  theme_classic()+
  guides(size = "none", alpha = "none")+
  annotate(geom = "text", y = 46.5, x = 17, hjust = 0, color = "black", size = 3.5,
         label = bquote(Field ~ n[tests] == ~ .(n.tests.max.F) ~ "(" * n[Indiv] == ~ .(n.indiv.max.F) * ")"))+
  annotate(geom = "text", y = 43.5, x = 17, hjust = 0, color = "black", size = 3.5,
         label = bquote(Lab ~ n[tests] == ~ .(n.tests.max.L) ~ "(" * n[Indiv] == ~ .(n.indiv.max.L) * ")"))+
  annotate(geom = "text", y = -4, x = 17, hjust = 0, color = "black", size = 3.5,
         label = bquote(Field ~ n[tests] == ~ .(n.tests.min.F) ~ "(" * n[Indiv] == ~ .(n.indiv.min.F) * ")"))+
  annotate(geom = "text", y = -7, x = 17, hjust = 0, color = "black", size = 3.5,
         label = bquote(Lab ~ n[tests] == ~ .(n.tests.min.L) ~ "(" * n[Indiv] == ~ .(n.indiv.min.L) * ")"))+
  scale_fill_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD", "FIELD" = "#00CAFF"),
                    name = "", labels = c("12ºC", "17ºC", "22ºC", "27ºC", "Field", "V2", "V1"))+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD", "FIELD" = "#00CAFF"),
                     name = "", labels = c("12ºC", "17ºC", "22ºC", "27ºC", "Field", "V2", "V1"))+
  scale_y_continuous(limits = c(-8, 47), breaks = c(seq(0, 45, 5)))
ggformat(CTplotFIELD, x_title = expression(Test~start~temperature~(degree*C)), y_title = expression(Temperature~tolerance~(degree*C)), print = F, size_text = 12)
CTplotFIELD <- CTplotFIELD + theme(legend.position = "none", 
                                     legend.title = element_text(), 
                                     axis.title.y = element_blank())
CTplotFIELD

## [main 1b] mean temp *************************************
CTplotFIELD2<-ggplot(data=data)+
  geom_line(data = data.ctmax.predict,
            mapping = aes(y = pred.ctmax, x = mean.Env.Temp, alpha = NULL, size = NULL, shape = NULL, fill = NULL, color =NULL),
            lty = "dashed", color = "black", linewidth=0.5)+
  geom_point(aes(y=temp_tolerance, x=mean.Env.Temp, shape = Stat_Var, alpha = Stat_Var, size = Stat_Var, fill = temp_treatment, color = temp_treatment))+
  scale_shape_manual(values = c(5, 21), name = "")+
  scale_size_manual(values = c(2, 3))+
  scale_alpha_manual(values = c(1, 0.5))+
  theme_classic()+
  guides(size = "none", alpha = "none")+
  theme_classic()+
  scale_fill_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D","V1" = "#A3ABBD", "V2" = "black",  "FIELD" = "#00CAFF"),
                     name = "", labels = c("12ºC", "17ºC", "22ºC", "27ºC", "Field", "V1", "V2"))+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V1" = "#A3ABBD","V2" = "black",  "FIELD" = "#00CAFF"),
                     name = "", labels = c("12ºC", "17ºC", "22ºC", "27ºC","Field", "V1", "V2"))+
  scale_y_continuous(limits = c(-8, 47), breaks = c(seq(0, 45, 5)))+
  annotate(geom = "text", y = 44, x = 15, hjust = 0, color = "black", size = 3.5,
         label = bquote( CT[max] == ~ .(CTmax.int) ~ "+" ~ .(CTmax.slopeMean) * T[mean] ~ + .(CTmax.slopeTL) * "TL"))
ggformat(CTplotFIELD2, x_title = expression(Mean~daily~T~(degree*C)), y_title = expression(Temperature~tolerance~(degree*C)), print = F, size_text = 12)
CTplotFIELD2<-CTplotFIELD2+theme(legend.position = c(0.88, 0.5),
                               legend.key.size = unit(0.2, "cm"),
                               legend.margin = margin(-0.8,0,0,0, unit="cm"))
CTplotFIELD2

# correlations with means, min, max 
## [main 1c] max temp *************************************
CTplotFIELD1<-ggplot(data=data, aes(y=temp_tolerance, x=max.Env.Temp,  shape = Stat_Var, alpha = Stat_Var, size = Stat_Var, fill = temp_treatment, color = temp_treatment))+
  geom_line(data = data.ctmin.predict, mapping = aes(y = pred.ctmin, x = max.Env.Temp, size = NULL, alpha = NULL, shape = NULL, fill = NULL, color =NULL), lty = "dashed", color = "black", size=0.5)+
  geom_point()+
  scale_shape_manual(values = c(5, 21), name = "Location")+
  scale_size_manual(values = c(2, 3))+
  scale_alpha_manual(values = c(1, 0.5))+
  guides(size = "none", alpha = "none")+
  # geom_errorbar( data=dat.ctF.l[dat.ctF.l$test == "CTmin",], mapping = aes(ymin = mean_CT - sd_CT, ymax = mean_CT + sd_CT),  width=0.2, size=0.5, alpha=0.8)+
  theme_classic()+
  scale_fill_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD", "FIELD" = "#00CAFF"), name = "Treatment ºC" )+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD", "FIELD" = "#00CAFF"), name = "Treatment ºC" )+
  scale_y_continuous(limits = c(-8, 47), breaks = c(seq(0, 45, 5)))+
  annotate(geom = "text", y = -5, x = 11.5, hjust = 0, color = "black", size = 3.5,
         label = bquote( CT[min] == ~ .(CTmin.int) ~ "+" ~ .(CTmin.slopeMax) * T[max] ~ .(CTmin.slopeTL) * "TL"))
 ggformat(CTplotFIELD1, x_title = expression(Max~daily~T~(degree*C)), y_title = expression(Temperature~tolerance~(degree*C)), print = TRUE, size_text = 12)
CTplotFIELD1 <- CTplotFIELD1 + theme(legend.position = "none", 
                                     legend.title = element_text(), 
                                     axis.title.y = element_blank())

cowplot:::plot_grid(CTplotFIELD2, CTplotFIELD1,
                    labels = "AUTO", 
                    nrow =1,
                    ncol=2,
                    align = "hv", 
                    label_size = 16,
                    label_x = c(0.18, 0.18),
                    label_y = c(0.9, 0.9),
                    rel_widths = c(1,1)) %>% 
ggsave(filename = "Figures/Figure3_EnvCorrelations.png", width = 8, height = 4)
# ggsave(plot_CTsizeF, filename = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Carpinteria marsh fish /DataAnalysis/PLOTS/Field_CTmaxCT_size.png", width = 7, height = 3.5)
 
## 1. [not in ms]: LAB - Growth, mm------
plot_growth2 <- ggplot()+
  geom_line(data = pred.datamm, mapping = aes(x = temp, y = pred.GRmm), color = "grey30", lwd = 0.9)+
  geom_line(data = pred.datamm, mapping = aes(x = temp, y = pred.GRmm+pred.GRmm.se), color = "grey30", lwd = 0.2, lty=2)+
  geom_line(data = pred.datamm, mapping = aes(x = temp, y = pred.GRmm-pred.GRmm.se), color = "grey30", lwd = 0.2, lty=2)+
  geom_line(data = growth.mod.data, aes(x = temp, GR_mm_d, color = temp_treatment.x), alpha=0.3)+
  geom_point(data = growth.mod.data, aes(x = temp, GR_mm_d, color = temp_treatment.x, fill = temp_treatment.x), alpha= 0.3, pch = 22, size=2)+
  # geom_text(data = growth.mod.data, aes(x = temp, GR_mm_d, label = paste("n=",n.y, sep="")), nudge_x = 0.6, size = 3)+
  geom_text(data = NULL, aes(x = c(12, 12, 12, 17, 17, 17, 22, 22, 22, 27, 27), y = c(0.02, 0.0125, 0.005, 0.02, 0.0125, 0.005, 0.02, 0.0125, 0.005,0.02, 0.0125)),
            label = c(14, 23, 21, 22, 24, 18,22, 25, 14, 11, 6), size = 3)+
  geom_text(data = NULL, aes(x = 12, y = 0.03), label = "n (tank)", size = 3)+
  geom_point(data = growth.mod.mean, aes(x = temp, mean_GRmm, color = temp_treatment.x, fill = temp_treatment.x), alpha= 1, pch = 22, size=3)+
  geom_errorbar(data = growth.mod.mean, mapping = aes(x = temp, ymin = mean_GRmm-sd_GRmm, ymax = mean_GRmm+sd_GRmm, color = temp_treatment.x), size=0.2, width = 0.1)+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V" = "#A3ABBD") )+
  scale_fill_manual(values=c("12" ="#00518C",  "17" = "#008A60","22" ="#DBA11C","27" ="#BE647D", "V2" = "black", "V1" = "#A3ABBD") )+
  theme_classic()+
  xlim(10,30)+
  ylim(0, 0.2)
ggformat(plot_growth2, x_title = "Temperature treatment", y_title = expression(Growth~rate~(mm~d^-1)), print=F, size_text = 12)
plot_growth2 <- plot_growth2 + theme (legend.position = "none")

plot_growth2B <- ggplot()+
  geom_point(data2.sum.w[c(!(is.na(data2.sum.w$GR_mm_d)) & c(data2.sum.w$temp_treatment.x=="V1")),],
             mapping = aes(x = temp_treatment.x, y = GR_mm_d), color = "#A3ABBD", fill = "#A3ABBD", pch = 22, size=2, alpha = 0.3)+
  geom_point(data2.sum.w[c(!(is.na(data2.sum.w$GR_mm_d)) & c(data2.sum.w$temp_treatment.x=="V2")),],
             mapping = aes(x = temp_treatment.x, y = GR_mm_d), color = "black", fill = "black",  pch = 22, size=2, alpha = 0.3)+
  geom_point(data = growth.mod.mean[growth.mod.mean$temp_treatment.x == "V1" | growth.mod.mean$temp_treatment.x == "V2",],
             aes(x = temp_treatment.x, mean_GRmm, color = temp_treatment.x, fill = temp_treatment.x),
             alpha= 1, pch = 22, size=3)+
  geom_errorbar(data = growth.mod.mean[growth.mod.mean$temp_treatment.x == "V1" | growth.mod.mean$temp_treatment.x == "V2",],
                mapping = aes(x = temp_treatment.x, ymin = mean_GRmm-sd_GRmm, ymax = mean_GRmm+sd_GRmm, color = temp_treatment.x),
                size=0.2, width = 0.1)+
  # geom_text(data = data2.sum.w, aes(x = temp_treatment.x, GR_mm_d, label = paste("n=",n.y, sep="")), nudge_x = -0.6, size = 3)+
  geom_text(data = NULL, aes(x = c(1,1, 2,2,2), y = c(0.02, 0.0125, 0.02, 0.0125, 0.005)),
            label = c(20, 16, 15, 17, 17), size = 3)+
  scale_color_manual(values=c("12" ="#00518C", "17" = "#008A60","22" ="#DBA11C", "27" = "#BE647D", "V2" = "black", "V1" = "#A3ABBD") )+
  scale_fill_manual(values=c("12" ="#00518C",  "17" = "#008A60","22" ="#DBA11C","27" ="#BE647D", "V2" = "black", "V1" = "#A3ABBD") )+
  theme_classic()+
  # ylim(-15, 40)+
  ylim(0, 0.2)
ggformat(plot_growth2B, x_title = "Temperature treatment", y_title = expression(Growth~rate~(mm~d^-1)), print=F, size_text = 12)
plot_growth2B <- plot_growth2B + theme (legend.position = "none", 
                                        axis.title.y = element_blank(),
                                        axis.text.y = element_blank(),
                                        axis.title.x = element_blank())

plot.size1 <- cowplot::plot_grid(plot_growth2, plot_growth2B,
                                 labels = c("AUTO"),
                                 rel_widths = c(1, 0.2),
                                 align = "h")
ggsave(plot.size1, filename = "./Figures/FigureS2_growth.png", width = 6, height = 4.5)


## 2. [suppl 1b]: Figures CTmax CT between-test variation: --------
data$TestID2<-factor(data$TestID2, levels = c("LAB-12", "LAB-17", "LAB-22", "LAB-27", "LAB-V1", "LAB-V2", 
                     "FIELD-44382", "FIELD-44386", "FIELD-44388", "FIELD-44391",
                     "FIELD-44460", "FIELD-44461"))
data.sum$TestID2<-factor(data.sum$TestID2, levels = c("LAB-12", "LAB-17", "LAB-22", "LAB-27", "LAB-V1", "LAB-V2", 
                     "FIELD-44382", "FIELD-44386", "FIELD-44388", "FIELD-44391",
                     "FIELD-44460", "FIELD-44461"))

plot_CTmaxmin <- ggplot(data[!c(data$Field_Lab == "LAB" & !c(data$temp_treatment=="V1" | data$temp_treatment=="V2")),],
                        aes(x = TestID2, y = temp_tolerance, fill = mean.Env.Temp, color = mean.Env.Temp))+
  geom_text_repel(data = data.sum[!c(data.sum$TestID2 == "LAB-27" | data.sum$TestID2 == "LAB-22" | data.sum$TestID2 == "LAB-17" | data.sum$TestID2 == "LAB-12") & data.sum$Test == "CTmin",],
                  aes(x = TestID2, y = mean_CT, label = n, fill = NULL), min.segment.length = 0.1, nudge_x = 0.2, direction = "y", hjust = "left", size = 2, color = "black")+
  geom_text_repel(data = data.sum[!c(data.sum$TestID2 == "LAB-27" | data.sum$TestID2 == "LAB-22" | data.sum$TestID2 == "LAB-17" | data.sum$TestID2 == "LAB-12") & data.sum$Test == "CTmax",],
                  aes(x = TestID2, y = mean_CT, label = n, fill = NULL), min.segment.length =  0.1, nudge_x = 0.2, direction = "y", hjust = "left", size = 2, color = "black")+
  geom_point(alpha=0.7, size=1.5, pch=21)+
  geom_point(data.sum[!c(data.sum$TestID2 == "LAB-27" | data.sum$TestID2 == "LAB-22" | data.sum$TestID2 == "LAB-17" | data.sum$TestID2 == "LAB-12"),],
             mapping = aes(y = mean_CT, x = TestID2, fill = NULL), alpha=1, size=3, color = "black", pch=21)+
  geom_vline(xintercept = c(2.5), linetype="dashed")+
  scale_fill_viridis_c(option = "D", name = "Mean T (ºC)")+
  scale_color_viridis_c(option = "D", guide = NULL)+
  scale_x_discrete(labels = c("LAB-V1" = "V1", "LAB-V2" = "V2", "FIELD-44382" = "Jul-5", "FIELD-44386" = "Jul-9", "FIELD-44388" = "Jul-11", "FIELD-44391" = "Jul-14",
                     "FIELD-44460" = "Sep-21", "FIELD-44461" = "Sep-22"))+
  theme_classic()+
  ylim(0, 42)
ggformat(plot_CTmaxmin, x_title = "Temperature treatment", y_title = expression(Temperature~tolerance~(degree*C)), size_text = 12,  print=FALSE)
plot_CTmaxmin <- plot_CTmaxmin + theme(legend.position = c(0.87, 0.5),
                                     axis.title.y = element_blank())

###  ---------
plot_CTmaxminLAB <- ggplot(data[c(data$Field_Lab == "LAB" & !c(data$temp_treatment=="V1" | data$temp_treatment=="V2")),],
                        aes(x = temp, y = temp_tolerance, fill = TestID2, color = TestID2))+
  geom_line(data = predCTmaxStatic, aes(x = temp, y = upr, fill = NULL, color = NULL), color = "grey30", lwd = 0.2, lty=2 )+
  geom_line(data = predCTmaxStatic, aes(x = temp, y = fit, fill = NULL, color = NULL), color = "grey30", lwd = 0.2, lty=1 )+
  geom_line(data = predCTmaxStatic, aes(x = temp, y = lwr, fill = NULL, color = NULL), color = "grey30", lwd = 0.2, lty=2)+
  geom_line(data = predCTminStatic, aes(x = temp, y = upr, fill = NULL, color = NULL), color = "grey30", lwd = 0.2, lty=2 )+
  geom_line(data = predCTminStatic, aes(x = temp, y = fit, fill = NULL, color = NULL), color = "grey30", lwd = 0.2, lty=1 )+
  geom_line(data = predCTminStatic, aes(x = temp, y = lwr, fill = NULL, color = NULL), color = "grey30", lwd = 0.2, lty=2)+# intercept = fixef(mod1)[1], slope = fixef(mod1)[2])+
  geom_text_repel(data = data.sum[c(data.sum$TestID2 == "LAB-27" | data.sum$TestID2 == "LAB-22" | data.sum$TestID2 == "LAB-17" | data.sum$TestID2 == "LAB-12") & data.sum$Test == "CTmin",],
                  aes(x = temp, y = mean_CT, label = n, fill = NULL), min.segment.length = 0.1, nudge_x = 1, direction = "y", hjust = "left", size = 2, color = "black")+
  geom_text_repel(data = data.sum[c(data.sum$TestID2 == "LAB-27" | data.sum$TestID2 == "LAB-22" | data.sum$TestID2 == "LAB-17" | data.sum$TestID2 == "LAB-12") & data.sum$Test == "CTmax",],
                  aes(x = temp, y = mean_CT, label = n, fill = NULL), min.segment.length =  0.1, nudge_x = 1, direction = "y", hjust = "left", size = 2, color = "black")+
  geom_point(alpha=0.7, size=1.5, pch=21)+
  geom_point(data.sum[c(data.sum$TestID2 == "LAB-27" | data.sum$TestID2 == "LAB-22" | data.sum$TestID2 == "LAB-17" | data.sum$TestID2 == "LAB-12"),],
             mapping = aes(y = mean_CT, x = temp, fill = NULL), alpha=1, size=3, color = "black", pch=21)+
  scale_color_manual(values=c("LAB-12" ="#00518C", "LAB-17" = "#008A60","LAB-22" ="#DBA11C", "LAB-27" = "#BE647D") )+
  scale_fill_manual(values=c("LAB-12" ="#00518C",  "LAB-17" = "#008A60","LAB-22" ="#DBA11C","LAB-27" ="#BE647D") )+
  theme_classic()+
  ylim(0, 42)+
  xlim(11, 28)
ggformat(plot_CTmaxminLAB, x_title = "Temperature treatment", y_title = expression(Temperature~tolerance~(degree*C)), size_text = 12,  print=FALSE)
plot_CTmaxminLAB <- plot_CTmaxminLAB + theme(legend.position = "none")

cowplot::plot_grid(plot_CTmaxminLAB, plot_CTmaxmin,
                   labels = c('A', 'B'),
                   rel_widths = c(0.5, 1),
                   align = "hv") %>% 
ggsave( filename = "./Figures/FigureS2_CTresults.png", width = 8, height = 4)



## 3. [suppl 3] scaling plot ------
plot_CTsize_mg<-ggplot(data=data, aes(y=temp_tolerance , x=mass_mg,
                                   group = interaction(Test, temp),
                                   fill = Field_Lab, color = Field_Lab)) +
  geom_point( colour="black", pch=21, size=2)+
  facet_wrap(.~Test, scale="free")+
  theme_classic()
ggformat(plot_CTsize_mg, x_title = "Body mass (mg)", y_title = expression(Temperature~tolerance~(degree*C)), print = TRUE)

plot_CTsize_mm<-ggplot(data=data, aes(y=temp_tolerance , x=TL_mm,
                                   group = interaction(Test, temp),
                                   fill = Field_Lab, color = Field_Lab)) +
  geom_point( colour="black", pch=21, size=2)+
  facet_wrap(.~Test, scale="free")+
  theme_classic()
ggformat(plot_CTsize_mm, x_title = "Body Lentgth (mm)", y_title = expression(Temperature~tolerance~(degree*C)), print = TRUE, size_text = 12)
ggsave( filename = "./Figures/Figure3_scaling_TLmm.png", width = 8, height = 4)





