#### estimate the growth rates for the species with dendromiters in Tanguro for the year 2019

### libraries
library(tidyverse)
library(dplyr)
library(dplyrUtil, quietly = TRUE)
library(ggplot2)

###reading the data of parenchyma and starch out of the ImageJ files
path="/Users/_dherrera/Documents/balzan_project/Manuscript_submissions/NSC_seasonality_and_growth/Journal_of_ecology/revision2/scripts_and_data/"
setwd(path)

source("scripts/R_function_stat_smooth.R")
source("scripts/functions.R")

# fun <- function(x, ...){
#   c(m=mean(x, ...),sd=sd(x, ...), l=length(x))
# }

#### estimating growth rates

monthly_growth=read.csv2("data/dendrometro_data.csv", sep = ",")

monthly_growth2=monthly_growth%>%gather(date, dbh, dbh_Aug.18, dbh_Sep.18, dbh_Oct.18,
                                        dbh_Nov.18, dbh_Dec.18, dbh_Jan.19, dbh_Feb.19, dbh_Mar.19, 
                                        dbh_Apr.19,dbh_May.19, dbh_Jun.19, dbh_Jul.19, dbh_Aug.19, dbh_Sep.19, dbh_Oct.19, dbh_Nov.19,
                                        dbh_Dec.19, dbh_Jan.20, dbh_Feb.20, dbh_Mar.20, dbh_May.20)
monthly_growth2=monthly_growth2[,-c(4, 8:27)]
monthly_growth2=monthly_growth2%>%separate(date, c("measurement", "date"), sep="_")

lev=unique(monthly_growth2$date)
monthly_growth2$date2=factor(monthly_growth2$date, levels=lev ,
                             labels=lev)
monthly_growth2$date3=as.numeric(monthly_growth2$date2)

monthly_growth2$dbh=as.numeric(monthly_growth2$dbh)
str(monthly_growth2)



ggplot(monthly_growth2, aes(x=date3, y=dbh, col=Scientific_name))+geom_point()


monthly_growth2$rate_growth=NA
monthly_growth2$rate_growth_summ_three_months_prior=NA
monthly_growth2$rate_growth_mean_three_months_prior=NA
monthly_growth2$rate_growth_summ_three_months_in=NA
monthly_growth2$rate_growth_mean_three_months_in=NA
monthly_growth2$rate_growth_summ_three_months_after=NA
monthly_growth2=monthly_growth2 %>% group_by(num_placa)%>% 
  mapGroups(function(groups_ID){
    for(i in 2:nrow(groups_ID)){
      groups_ID$rate_growth[1]=0
      groups_ID$rate_growth[i]=groups_ID$dbh[i]-groups_ID$dbh[i-1]
    }
    return(groups_ID)
  })%>% ungroup()

### excluding outliers and possible mistakes (we consider that growth or shrinkage of more than 1 cm per month could have been a reading mistake and then we eliminated this data)
monthly_growth2$rate_growth[which(monthly_growth2$rate_growth>=1 | monthly_growth2$rate_growth<=-1)]=NA

ggplot(monthly_growth2, aes(x=date2, y=rate_growth, col=Scientific_name))+geom_point()+facet_grid(~Scientific_name)


# estimating mean averages every three months to reduce the effect of water dynamics in the growth data 
monthly_growth2=monthly_growth2 %>% group_by(num_placa)%>% 
  mapGroups(function(groups_ID){
    uplimit=nrow(groups_ID)-2
for(i in 3:uplimit){
  groups_ID$rate_growth_summ_three_months_prior[i]=groups_ID$rate_growth[i]+groups_ID$rate_growth[i-1]+groups_ID$rate_growth[i-2]
  groups_ID$rate_growth_mean_three_months_prior[i]=mean(c(groups_ID$rate_growth[i],groups_ID$rate_growth[i-1], groups_ID$rate_growth[i-2]))
  groups_ID$rate_growth_summ_three_months_in[i]=groups_ID$rate_growth[i]+groups_ID$rate_growth[i-1]+groups_ID$rate_growth[i+1]
  groups_ID$rate_growth_mean_three_months_in[i]=mean(c(groups_ID$rate_growth[i],groups_ID$rate_growth[i-1],groups_ID$rate_growth[i+1]))
  groups_ID$rate_growth_summ_three_months_after[i]=groups_ID$rate_growth[i]+groups_ID$rate_growth[i+1]+groups_ID$rate_growth[i+2]
}
    return(groups_ID)
  })%>% ungroup()


#### visualzing the data ######
ggplot(monthly_growth2, aes(x=date3, y=rate_growth_summ_three_months_prior, col=Scientific_name))+geom_point()+facet_grid(~Scientific_name)
ggplot(monthly_growth2, aes(x=date2, y=rate_growth_summ_three_months_prior, col=Scientific_name))+geom_point()+facet_grid(~Scientific_name)
ggplot(monthly_growth2, aes(x=date3, y=rate_growth_mean_three_months_prior, col=Scientific_name))+geom_point()+facet_grid(~Scientific_name)
ggplot(monthly_growth2, aes(x=date3, y=rate_growth_summ_three_months_in, col=Scientific_name))+geom_point()+facet_grid(~Scientific_name)
ggplot(monthly_growth2, aes(x=date3, y=rate_growth_mean_three_months_in, col=Scientific_name))+geom_point()+facet_grid(~Scientific_name)
ggplot(monthly_growth2, aes(x=date3, y=rate_growth_summ_three_months_after, col=Scientific_name))+geom_point()+facet_grid(~Scientific_name)


### calculating mean averages per species #####
monthly_growth_species_means=monthly_growth2 %>%group_by(Scientific_name, date2)%>% 
  summarise(mean_growth_rate=mean(rate_growth, na.rm=T),
            sd_growth_rate=sd(rate_growth, na.rm=T),
            mean_growth_three_months_in=mean(rate_growth_mean_three_months_in, na.rm=T))

level=unique(monthly_growth_species_means$date2)

monthly_growth_species_means$date3=as.numeric(factor(monthly_growth_species_means$date2, levels=level, labels=level))
ggplot(monthly_growth_species_means, aes(x=date2, y=mean_growth_rate))+geom_point(aes(col=Scientific_name))+geom_line(aes(x=date3, col=Scientific_name))

ggplot(monthly_growth_species_means, aes(x=date2, y=mean_growth_three_months_in))+geom_point(aes(col=Scientific_name))+geom_line(aes(x=date3, col=Scientific_name))+
  facet_grid(~Scientific_name)



monthly_growth2$ID_sp=monthly_growth2$comon_name
comon_name_key=unique(monthly_growth2$comon_name)
ID_sp_key=c("Sac", "ab", "Dac")

library(plyr)
monthly_growth2$ID_sp=mapvalues(monthly_growth2$ID_sp,comon_name_key,ID_sp_key)

monthly_growth2$ID_sp2=monthly_growth2$comon_name
comon_name_key=unique(monthly_growth2$comon_name)
ID_sp_key2=c("S. guianensis\n(semi-deciduous/parenchyma-storing-species)",
             "O. leucoxylon\n(evergreen/parenchyma-storing-species)",
             "D. microcarpa\n(semi-deciduous/fiber-storing-species)")

monthly_growth2$ID_sp2=mapvalues(monthly_growth2$ID_sp2,comon_name_key,ID_sp_key2)

detach("package:plyr", unload=TRUE)

monthly_growth2$ID=paste(monthly_growth2$ID_sp, monthly_growth2$num_placa, sep="")

monthly_growth2=monthly_growth2%>%separate(date, c("month2", "year"))
monthly_growth2$month=paste(monthly_growth2$month2, monthly_growth2$year, sep = "")

monthly_growth2$month_numeric=match(monthly_growth2$month2, month.abb)
monthly_growth2$year=as.numeric(monthly_growth2$year)+2000
monthly_growth2$date=paste(monthly_growth2$year, monthly_growth2$month_numeric, "01", sep = "-")
monthly_growth2$date2=as.POSIXct(strptime(as.character(monthly_growth2$date), format = "%Y-%m-%d"))


monthly_means_growt=monthly_growth2%>%filter(date3%in%c(7:19))%>%group_by(Scientific_name, month_numeric)
monthly_means_growt$Scientific_name=as.character(monthly_means_growt$Scientific_name)
monthly_means_growt$Scientific_name[which(monthly_means_growt$Scientific_name=="Sacoglotis guianensis ")]="Sacoglottis guianensis"
monthly_means_growt$Scientific_name[which(monthly_means_growt$Scientific_name=="Ocotea leocoxylon ")]="Ocotea leucoxylon"


#dev.off()
str(monthly_means_growt)
ID_sp2=unique(monthly_means_growt$ID_sp2)
png("outputs/growth_seasonality.png", width = 1400, height = 800)
ggplot(monthly_means_growt)+
  # geom_rect(aes(xmin = 6, xmax = 9, ymin = -Inf, ymax = Inf), fill = "grey86",
  #           alpha = 0.01) +
  # geom_rect(aes(xmin = 1, xmax = 4, ymin = -Inf, ymax = Inf), fill = "gray25",
  #           alpha = 0.01) +
  geom_boxplot(aes(x=as.factor(month_numeric), y=rate_growth_mean_three_months_in, col=Scientific_name), alpha=0.2)+
  geom_point(aes(x=month_numeric, y=rate_growth_mean_three_months_in, col=Scientific_name), alpha=0.2)+
  stat_smooth(aes(x=month_numeric, y=rate_growth_mean_three_months_in, col=Scientific_name), se=T)+
  scale_x_discrete("month", breaks = seq(1,12, 2), 
                   labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  ylab("Mean monthly growth rate (cm)")+
  ylim(-0.1,0.25)+
  #scale_y_continuous("Growth rate (cm/month)", limits = c(-0.3, 0.3))+
  scale_color_manual(values=c("orangered4", "aquamarine4","aquamarine4"))+
  facet_grid(~ID_sp2)+
   theme_bw()+
   theme(axis.text.x=element_text(size = rel(3)),
         axis.text.y=element_text(size = rel(3)),
         axis.title.x = element_text(size = rel(3)),
         axis.title.y = element_text(size = rel(3)),
         strip.text = element_text(size = rel(2)),
         legend.position="none",
         plot.title = element_text(size = rel(3), face = "italic"),
         strip.text.x=element_text(size = 20, margin = margin(10), vjust = 1, face = "italic"))
dev.off()

### Comparison between three months accumulated growth withing each species 

str(monthly_means_growt)
monthly_means_growt$month_numeric=factor(monthly_means_growt$month_numeric)
monthly_means_growt$month_numeric2=as.factor(monthly_means_growt$month)

Wilcoxon_monthly_growth=monthly_means_growt%>%group_by(Scientific_name)%>%
  mapGroups(function(df){
    Wc_test=wilcox_test(rate_growth_mean_three_months_in~month_numeric2, data=df, paired = T, detailed = T)
    return(Wc_test)
  })%>%ungroup()

Wilcoxon_monthly_growth%>%filter(p.adj<0.05)

save(monthly_growth2, file="data/growth_data_tanguro_2018_2020.RData")
