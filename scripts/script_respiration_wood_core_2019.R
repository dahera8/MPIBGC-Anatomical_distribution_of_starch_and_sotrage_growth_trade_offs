#### estimation of the stem respiration patterns of the incubated wood cores in 2019

library(tidyverse)
library(ggplot2)
library(dplyr)
library(dplyrUtil, quietly = TRUE)
library(rstatix)

setwd("/Users/_dherrera/MPIBGC-Anatomical_distribution_of_starch_and_sotrage_growth_trade_offs/")


###### processing data from May 2019 #####
data_incubation_chambers=read.csv2("data/resp_from_incubation_flasks_May2019.csv", header = T, sep= ";")

data_incubation_chambers=data_incubation_chambers[seq(1, nrow(data_incubation_chambers), 2), ]
data_incubation_chambers=data_incubation_chambers%>%separate(FlaskID, c("FlaskID", "Flask_vol"))
data_incubation_chambers$Flask_vol=as.numeric(data_incubation_chambers$Flask_vol)


### transform ppm to mg/m3

presure= 1013 ## pressure at the moment of measurement hpa
Temp= 25 ## temperature at the moment of measurement c
vol_mol=22.4 ## volume of 1 mol at 1 atm at O Celcious
M=44 ## molecular weight of CO2
data_incubation_chambers$CO2_mg_m3=data_incubation_chambers$ppmCO2*(M/vol_mol)*(273/(273+Temp))*(presure/1013)

## caltulating the amount of CO2 in mg, volume of the flask is in cc so m3 from the CO2_mg_m3 has to be transformer to cc
data_incubation_chambers$CO2_mg=data_incubation_chambers$CO2_mg_m3*1e-6*(data_incubation_chambers$Flask_vol*2)

data_incubation_chambers$C_mg=data_incubation_chambers$CO2_mg*0.27
data_incubation_chambers$CO2_mg_ml=data_incubation_chambers$CO2_mg_m3*1e-6

## calculating respiration rate in mg per hour 
## the incubation lasted for 36 hours
data_incubation_chambers$CO2_mg_hr=data_incubation_chambers$CO2_mg/36

ggplot(data_incubation_chambers, aes(x=ID_sp, y=CO2_mg_hr))+geom_boxplot()
ggplot(data_incubation_chambers, aes(x=ID_sp, y=ppmCO2))+geom_boxplot()

## standardize per volume then we will obtain mg/hr*m3

probe_volume_m3=(pi*(0.5/2)^2)*6/10000

data_incubation_chambers$CO2_mg_hr_m3=data_incubation_chambers$CO2_mg_hr/probe_volume_m3

ggplot(data_incubation_chambers, aes(x=ID_sp, y=CO2_mg_hr_m3))+geom_boxplot()

#################################################################################################################################
#### processing data for respiration from the incubation in July 2018

data_CO2_July2018=read.csv2("data/13C_results_jul2018 .csv", header = T, sep = ";")
data_CO2_July2018=data_CO2_July2018[1:37,]
data_CO2_July2018$CO2_mg_m3=(44/22.4)*(273/(273+Temp))*(presure/1013)*data_CO2_July2018$ppm
data_CO2_July2018$CO2_mg=data_CO2_July2018$CO2_mg_m3*1e-6*160 ## 160 cc is from the sum of the two flasks 115 and 45 cc. 
data_CO2_July2018$CO2_mg_hr=data_CO2_July2018$CO2_mg/12 ## 12 hours of incubation. 

data_CO2_July2018$month="July"
data_CO2_July2018$Date="218-07-01"
data_CO2_July2018$TreeID=data_CO2_July2018$ID2
data_CO2_July2018$ID_sp=data_CO2_July2018$ID3
data_CO2_July2018$ppmCO2=data_CO2_July2018$ppm

ggplot(data_CO2_July2018, aes(x=ID3, y=ppm))+geom_boxplot()

ggplot(data_CO2_July2018, aes(x=ID3, y=CO2_mg_hr))+geom_boxplot()+
  geom_boxplot(data=data_incubation_chambers, aes(x=ID_sp, y=CO2_mg_hr), col="red")


### merging Jul 2018 and may 2019
for_merge_CO2_Jul2018=data_CO2_July2018%>%select(ID_sp, ppmCO2, CO2_mg_hr, month, Date, TreeID)
for_merge_CO2_May2019=data_incubation_chambers%>%select(ID_sp, ppmCO2, CO2_mg_hr, month, Date, TreeID)

for_merge_CO2_May2019$ID_sp=as.character(for_merge_CO2_May2019$ID_sp)
for_merge_CO2_Jul2018$ID_sp=as.character(for_merge_CO2_Jul2018$ID_sp)

for_merge_CO2_May2019$ID_sp[which(for_merge_CO2_May2019$ID_sp=="D.mic")]="Dacryodes microcarpa"
for_merge_CO2_May2019$ID_sp[which(for_merge_CO2_May2019$ID_sp=="S.gui")]="Sacoglotis guianensis"
for_merge_CO2_May2019$ID_sp[which(for_merge_CO2_May2019$ID_sp=="O.leu")]="Ocotea leucoxylon"

for_merge_CO2_Jul2018$ID_sp[which(for_merge_CO2_Jul2018$ID_sp=="Ocotea leocoxylon ")]="Ocotea leucoxylon"
for_merge_CO2_Jul2018$ID_sp[which(for_merge_CO2_Jul2018$ID_sp=="Sacoglotis guianensis ")]="Sacoglotis guianensis"


CO2_incubations=rbind(for_merge_CO2_Jul2018, for_merge_CO2_May2019)
CO2_incubations$month=factor(CO2_incubations$month, levels = c("May", "July"))
CO2_incubations$season=ifelse(CO2_incubations$month=="July", "Dry season", "Wet season")
CO2_incubations$season=factor(CO2_incubations$season, levels = c("Wet season", "Dry season"))
CO2_incubations$ID_sp2[CO2_incubations$ID_sp=="Ocotea leucoxylon"]="O. leucoxylon\n(evergreen/parenchyma-storing sepecies)"
CO2_incubations$ID_sp2[CO2_incubations$ID_sp=="Dacryodes microcarpa"]="D. microcarpa\n(semi-deciduous/fiber-storing sepecies)"
CO2_incubations$ID_sp2[CO2_incubations$ID_sp=="Sacoglotis guianensis"]="S. guianensis\n(semi-deciduous/parenchyma-storing sepecies)"

CO2_incubations$storage_strategy[CO2_incubations$ID_sp2=="O. leucoxylon\n(evergreen/parenchyma-storing sepecies)"]="Parenchyma_storage"
CO2_incubations$storage_strategy[CO2_incubations$ID_sp2=="D. microcarpa\n(semi-deciduous/fiber-storing sepecies)"]="Fiber_storage"
CO2_incubations$storage_strategy[CO2_incubations$ID_sp2=="S. guianensis\n(semi-deciduous/parenchyma-storing sepecies)"]="Parenchyma_storage"

### visualizing the comparisons #####
png(file="outputs/CO2_eflux_from_core_incubation.png",
    width = 1200,
    height = 800)
ggplot(CO2_incubations, aes(x=season, y=CO2_mg_hr, col=storage_strategy))+geom_boxplot(size=2)+
   ylab(expression(increment~ core~ CO[2]~eflux (mg/hr/cm^3)))+
  xlab("")+
  facet_wrap(~ID_sp2)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        strip.text.x=element_text(size = 17, margin = margin(10), vjust = 1, face = "italic"))
dev.off()

#### using anova for testing differences withing species and between sepcies
AOV_CO2_eflux=aov(log(CO2_mg_hr)~ID_sp*month, data=CO2_incubations)
plot(AOV_CO2_eflux)
summary(AOV_CO2_eflux)
TukeyHSD(AOV_CO2_eflux)


### using wilcoxon test for comparins respiration rates within species
str(CO2_incubations)

CO2_pair_comparisons=CO2_incubations%>%group_by(ID_sp)%>%
  mapGroups(function(df){
    unique_ID=as.integer(names(which(table(df$TreeID)>1)))
    df=df%>%filter(TreeID%in%unique_ID)
    Wc_test=wilcox_test(CO2_mg_hr~month, data=df, paired = T)
    return(Wc_test)
  })%>%ungroup()


