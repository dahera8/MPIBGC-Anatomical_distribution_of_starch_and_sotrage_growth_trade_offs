
###packages needed 
library(tidyverse)
library(dplyr)
library(dplyrUtil, quietly = TRUE)
library(ggplot2)
library(broom)
library(ggpubr)
library(rstatix)
library(reshape2)
###reading the data of parenchyma and starch out of the ImageJ files

path="/Users/_dherrera/Documents/balzan_project/Manuscript_submissions/NSC_seasonality_and_growth/Journal_of_ecology/revision2/scripts_and_data/"
setwd(path)

source("scripts/R_function_stat_smooth.R")
source("scripts/functions.R")
load("data/growth_data_tanguro_2018_2020.RData") ##growth data calculated calculated in the script "Growth_data"

###### Calucaltion for the evergreen/parenchyma-storing-species O. leucoxylon #####


load("data/Ab_woodStarch.Rda")

############## Calculating mean statistics per individual, month and sampling depth ####################

## mean statistics by tree 

Ab_woodStarch_means_by_tree=Ab_woodStarch%>%group_by(ID, depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))


### mean statistics by species 

Ab_woodStarch_means_by_specie=Ab_woodStarch%>%group_by(depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

Ab_woodStarch_means_by_tree$month=factor(Ab_woodStarch_means_by_tree$month, 
                                         levels = c("May19", "Aug19", "Nov19", "Feb20"))
Ab_woodStarch_means_by_specie$month=factor(Ab_woodStarch_means_by_specie$month, 
                                           levels = c("May19", "Aug19", "Nov19", "Feb20"))

Ab_woodStarch_means_by_specie$month_season=Ab_woodStarch_means_by_specie$month

library(plyr)

months=unique(Ab_woodStarch_means_by_specie$month)
months_season=c("May 19 (transition)", "Aug 19 (Dry season)", "Nov 19 (transition)", "Feb 20 (wet season)")

Ab_woodStarch_means_by_specie$month_season=mapvalues(Ab_woodStarch_means_by_specie$month_season,
                                                     months,
                                                     months_season)
detach("package:plyr", unload=TRUE)

Ab_woodStarch_means_by_specie$month_season=factor(Ab_woodStarch_means_by_specie$month_season,
                                                  labels = c("May 19 (transition)", 
                                                             "Aug 19 (Dry season)", 
                                                             "Nov 19 (transition)", 
                                                             "Feb 20 (wet season)"))

png("outputs/radial_starch_content_O-leu_scale_seasons_withlengend.png", width = 800, height = 600)
ggplot()+
  geom_line(data=Ab_woodStarch_means_by_specie , aes(x=depth_mm, y=mean_starch_area, col=month_season), size=1)+
  geom_ribbon(data=Ab_woodStarch_means_by_specie, aes(x=depth_mm, 
                                                      ymin=mean_starch_area-sd_starch_area,
                                                      ymax=mean_starch_area+sd_starch_area,
                                                      fill= month_season),
              #fill="gray80",
              alpha=0.1)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("starch %")+
  labs(color="Seasons", fill="Seasons")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(2.5)),
        axis.title.y = element_text(size = rel(2.5)),
        strip.text = element_text(size = rel(3)),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(3)),
        legend.title = element_text(size = rel(3)),
        #legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


## statistical comparisons of the mean starch area at different depth ranges (ANOVA) and wilcoxon

Ab_pairwise_comparison=Ab_woodStarch_means_by_tree
Ab_pairwise_comparison=Ab_pairwise_comparison[-c(1:24), c(1,2,3,5)]
Ab_pairwise_comparison$time=as.factor(Ab_pairwise_comparison$month)
Ab_pairwise_comparison=Ab_pairwise_comparison[order(Ab_pairwise_comparison$month),]
rownames(Ab_pairwise_comparison)=NULL

Ab_pairwise_comparison$ID=as.factor(Ab_pairwise_comparison$ID)
AB_Anva2=Ab_pairwise_comparison%>%group_by(depth_mm)%>%nest()%>%
  mutate(ANOVA=map(data, ~ aov(mean_starch_area~month, data=.x)),
         tukey_comparison=map(ANOVA, TukeyHSD),
         tidied=map(ANOVA, tidy),
         tidied2=map(tukey_comparison, tidy), 
  )


names(AB_Anva2$tidied)=AB_Anva2$depth_mm
names(AB_Anva2$tidied2)=AB_Anva2$depth_mm
AB_Anva2$tidied. # this shows the p_values of the comparisons between months at each depth 
AB_Anva2$tidied2 # this is the posthoch Tikey comparison to see what months differed in starch concentration at each month
estimates_ANOVA_by_depth=AB_Anva2%>%unnest(tidied)

#wilcoxon comparison for mean starch area (starch concentrations) at different depth ranges in the wood cores
Ab_pairwise_comparison2=filter(Ab_pairwise_comparison, ID!="ab34288") ## We had to remove this tree because we lost one month of measurements for this individuals and Wilcoxon  test requires identical sample size for each group. 
Ab_pairwise_comparison2$ID=factor(Ab_pairwise_comparison2$ID)
Ab_pairwise_comparison2$depth_mm2=as.factor(Ab_pairwise_comparison2$depth_mm)
table(Ab_pairwise_comparison2$ID)

str(Ab_pairwise_comparison2)
AB_Anva3=Ab_pairwise_comparison2%>%group_by(depth_mm2)%>%
  wilcox_test(mean_starch_area~month, data=., paired = T)
  
### Figure to show the mean starch are at each depth for each sampling date 

ggplot(Ab_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_boxplot(data=Ab_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_point(aes(col=ID))+
  geom_point(data=Ab_woodStarch_means_by_specie, aes(x=month, y=mean_starch_area), col="black")+
  facet_wrap(~depth_mm, scales = "free")



### estimate the maximum and relative seasonal range in starch concentrations during the year of measurements in each wood depth

Ab_woodStarch_radial_activity=Ab_woodStarch_means_by_tree%>%group_by(ID, depth_mm)%>%
  summarise(activ_mean=max(mean_starch_area)-min(mean_starch_area),
            rel_activ_mean=(max(mean_starch_area)-min(mean_starch_area))/max(mean_starch_area)*100,
            activ_median=max(median_starch_are)-min(median_starch_are),
            activ_IQR=max(Int_quart_range)-min(Int_quart_range),
            max_starch=max(mean_starch_area),
            max_month=.data$month[max(which(.data$mean_starch_area==max(.data$mean_starch_area)))],
            min_starch=min(mean_starch_area),
            min_month=.data$month[min(which(.data$mean_starch_area==min(.data$mean_starch_area)))]
  )

Ab_woodStarch_radial_activity=na.omit(Ab_woodStarch_radial_activity)
#Ab_woodStarch_radial_activity=Ab_woodStarch_radial_activity%>%filter(ID!="ab26820")

Ab_woodStarch_radial_Activity_by_species=Ab_woodStarch_radial_activity%>%group_by(depth_mm)%>%
  summarise(mean_active_mean=mean(activ_mean),
            rel_mean_active_mean=mean(rel_activ_mean, na.rm=T),
            SD_rel_mean_active_mean=sd(rel_activ_mean, na.rm=T),
            mean_active_median=mean(activ_median),
            mean_active_IQR=mean(activ_IQR),
            SD_active_mean=sd(activ_mean),
            SD_active_median=sd(activ_median),
            Sd_active_IQR=sd(activ_IQR), 
            max=max(max_starch),
            mean_max=mean(max_starch),
            max_month_mode=names(which.max(table(max_month))),
            min=min(min_starch),
            mean_min=mean(min_starch),
            min_month_mode=names(which.max(table(min_month)))
  )


## ANOVA to compare changes in the starch change at different wood depths 

# check outliers 
Ab_woodStarch_radial_activity%>%group_by(depth_mm)%>%identify_outliers(activ_mean)

#normality assumption 
Ab_woodStarch_radial_activity%>%group_by(depth_mm)%>%shapiro_test(activ_mean)
ggqqplot(Ab_woodStarch_radial_activity, "activ_mean", facet.by = "depth_mm")

str(Ab_woodStarch_radial_activity)
seasonal_amplitud_anova=aov(activ_mean~as.factor(depth_mm), data=Ab_woodStarch_radial_activity)
summary(seasonal_amplitud_anova)
#plot(seasonal_amplitud_anova)
TukeyHSD(seasonal_amplitud_anova)

### visulization
png("outputs/metabolic_activity_O.leu.png", width = 800, height=600)
ggplot(data=Ab_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=mean_active_mean))+
  geom_line(data=Ab_woodStarch_radial_activity, aes(x=depth_mm, y=activ_mean, factor=ID), col="gray95")+
  geom_line( col="turquoise3")+
  geom_ribbon( aes(ymin=mean_active_mean-SD_active_mean, ymax=mean_active_mean+SD_active_mean,
                   fill= "var"),
               fill="turquoise3",
               alpha=0.3)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("seasonal amplitud starch %")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

png("outputs/relative_sorage_change_O.leu.png", width = 800, height = 600)
ggplot(data=Ab_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=rel_mean_active_mean))+
  geom_line(data=Ab_woodStarch_radial_activity, aes(x=depth_mm, y=rel_activ_mean, factor=ID), col="gray95")+
  geom_line( col="turquoise3")+
  geom_ribbon( aes(ymin=rel_mean_active_mean-SD_rel_mean_active_mean, ymax=rel_mean_active_mean+SD_rel_mean_active_mean,
                   fill= "var"),
               fill="turquoise3",
               alpha=0.3)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("relative seasonal amplitud")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

## mean content of starch in the entire core 
## wood density was measured in previous studies for the species in the area: Ocotea leucoxylon has a mean wood density value of 0.45 g/cm3

Ab_woodStarch_means_by_tree$segment_volume_cm3=pi*(0.5/2)^2*0.5 ### this in in cm3

Ab_woodStarch_means_by_tree$grams_wood=Ab_woodStarch_means_by_tree$segment_volume_cm3*0.45
Ab_woodStarch_means_by_tree$grams_starch=Ab_woodStarch_means_by_tree$grams_wood*Ab_woodStarch_means_by_tree$mean_starch_area/100
Ab_woodStarch_means_by_tree$sd_grams_starch=Ab_woodStarch_means_by_tree$grams_wood*Ab_woodStarch_means_by_tree$sd_starch_area/100

### estimating the mean starch wood core content and standard deviation in the entire wood core
Ab_woodStarch_means_entire_tree=Ab_woodStarch_means_by_tree%>%group_by(month, ID)%>%
  summarise(mean_starch_percentage=mean(mean_starch_area),
            median_starch_percentage=median(median_starch_are),
            mean_IQR=mean(Int_quart_range),
            starch_woodcore_content_by_mean=sum(mean_starch_area),
            starch_woodcore_content_by_median=sum(median_starch_are),
            starch_woodcore_content_grams=sum(grams_starch),
            sd_starch_woodcore_content_gr=sd(grams_starch),
            sd2_starch_woodcore_content_gr=mean(sd_grams_starch))

Ab_woodStarch_means_entire_tree$starch_woodcore_content_mg=Ab_woodStarch_means_entire_tree$starch_woodcore_content_grams*1000
Ab_woodStarch_means_entire_tree$sd_starch_woodcore_content_mg=Ab_woodStarch_means_entire_tree$sd2_starch_woodcore_content_gr*1000

#### estimating the man and standard deviation per species 
Ab_woodStarch_means_entire_tree_by_specie=Ab_woodStarch_means_entire_tree%>%group_by(month)%>%
  summarise(mean_starch_percentage=mean(mean_starch_percentage),
            sd_mean_starch_percentage=sd(mean_starch_percentage),
            median_starch_percentage=mean(median_starch_percentage),
            sd_median_starch_percentage=sd(median_starch_percentage),
            mean_IQR=mean(mean_IQR), 
            mean_starch_wood_content_by_mean=mean(starch_woodcore_content_by_mean),
            sd_starch_wood_content_by_mean=sd(starch_woodcore_content_by_mean),
            mean_starch_wood_content_by_median=mean(starch_woodcore_content_by_median),
            sd_starch_wood_content_by_medial=sd(starch_woodcore_content_by_median),
            mean_starch_wood_content_grams=mean(starch_woodcore_content_grams),
            sd_mean_starch_wood_content_grams=sd(starch_woodcore_content_grams),
            mean_starch_wood_content_mg=mean(starch_woodcore_content_mg),
            sd_mean_starch_wood_content_mg=sd(starch_woodcore_content_mg),
            sd2_mean_starch_wood_content_mg=sum(sd_starch_woodcore_content_mg)
  )

### note that the standard deviation calculated from the individuals sd is bigger than the standard deviation of starch content in each individual with the uncertainty propagated. 

png("outputs/Ab_starch_wood_content_comp_month.png", width = 800, height = 600)
ggplot(Ab_woodStarch_means_entire_tree, aes(x=month, y=starch_woodcore_content_mg))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_point(data=Ab_woodStarch_means_entire_tree_by_specie, aes(x=month, y=mean_starch_wood_content_mg),  col="black", pch=17, size=3)+
  xlab("")+
  ylab("Starch woodcore content (mg)")+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


### pairwise comparison with paired samples for the content of sugar in the entire wood core.

Ab_pairwise_comparison_entire_tree=Ab_woodStarch_means_entire_tree
Ab_pairwise_comparison_entire_tree$month2=as.factor(as.numeric(Ab_pairwise_comparison_entire_tree$month))
Ab_pairwise_comparison_entire_tree=Ab_pairwise_comparison_entire_tree[-c(which(Ab_pairwise_comparison_entire_tree$ID=="ab26820")),] 
names(Ab_pairwise_comparison_entire_tree)
Ab_pairwise_comparison_entire_tree=Ab_pairwise_comparison_entire_tree[ ,c(2,13,11)]
Ab_pairwise_comparison_entire_tree$ID=as.factor(Ab_pairwise_comparison_entire_tree$ID)
str(Ab_pairwise_comparison_entire_tree)


wilcox_test(starch_woodcore_content_mg~month2, data=Ab_pairwise_comparison_entire_tree, paired = T)


#### Comparisons between starch content and annual growth.  

monthly_growth2_summarized=monthly_growth2%>%filter(date3%in%c(8:20))%>%group_by(ID)%>%
  summarise(annual_growth_2019=sum(rate_growth),
            mean_monthly_growth_2019=mean(rate_growth),
            max_monthly_growth_2019=max(rate_growth),
            min_monthly_growth_2019=min(rate_growth)
  )


Ab_woodStarch_entire_tree=merge(Ab_woodStarch_means_entire_tree, monthly_growth2, by=c("ID", "month"))
Ab_woodStarch_entire_tree=merge(Ab_woodStarch_entire_tree, monthly_growth2_summarized, by=c("ID"))

str(Ab_woodStarch_entire_tree)
Ab_lm_annual_growth_starch=Ab_woodStarch_entire_tree%>%group_by(month)%>%nest()%>%
  mutate(lm_Ab_starch_annual_growth=map(data, ~lm(annual_growth_2019~starch_woodcore_content_mg, data = .x)),
         tidied_lm=map(lm_Ab_starch_annual_growth, tidy))

statistics_ab_lm=Ab_lm_annual_growth_starch%>%unnest(tidied_lm)
#plot(Ab_lm_annual_growth_starch$lm_Ab_starch_annual_growth[[4]])

Ab_woodStarch_entire_tree$significant_rel=as.character(Ab_woodStarch_entire_tree$month)
key_months_1=unique(as.character(Ab_woodStarch_entire_tree$month))
ID_key_signif2=c("Not_significant",
                 "Not_significant",
                 "Not_significant",
                 "Not_significant")

library(plyr)

Ab_woodStarch_entire_tree$significant_rel=mapvalues(Ab_woodStarch_entire_tree$significant_rel,
                                                    key_months_1,
                                                    ID_key_signif2)
detach("package:plyr", unload=TRUE)

png("outputs/annual_growth_vs_starch_content_O_leu.png",
    width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree, aes(x=starch_woodcore_content_mg, y=annual_growth_2019))+
  geom_point(color="turquoise2")+
  geom_smooth(aes(linetype=significant_rel), method = "lm", se=F, color="turquoise2", linetype = "dashed")+
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=8, xpos = 1,
                             ypos = 1.1, xpos2 = 1,
                             ypos2= 1)+
  xlab("starch content (mg)")+ ylab("annual growth 2019 (cm)")+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  facet_wrap(~month)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


# Calculation of the changes in starch content between months. 

Ab_woodStarch_entire_tree=Ab_woodStarch_entire_tree[order(Ab_woodStarch_entire_tree$ID, Ab_woodStarch_entire_tree$date3),]
Ab_woodStarch_entire_tree$starch_change=NA
Ab_woodStarch_entire_tree$starch_change_percent=NA
Ab_woodStarch_entire_tree$rel_starch_change=NA

Ab_woodStarch_entire_tree2=Ab_woodStarch_entire_tree%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 2:nrow(x)){
      x$starch_change[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1]
      x$starch_change_percent[i]=x$mean_starch_percentage[i]-x$mean_starch_percentage[i-1]
      x$rel_starch_change[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1])/x$starch_woodcore_content_mg[i]
    }
    return(x)
  }
  )%>%ungroup()


png("outputs/Ab_starch_mass_change_comp_month.png", width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (mg)")+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


png("outputs/Ab_starch_mass_change_rel_comp_month.png", width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=rel_starch_change))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  ylim(-4, 2)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

png("outputs/Ab_starch_mass_change_percent_month.png", width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change_percent))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  ylim(-4, 2)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


## pairwise comparisons with Wilcox test to compare starch changes within the species 

Ab_woodStarch_entire_tree2_wilcox_pair_comp=filter(Ab_woodStarch_entire_tree2, month!="May19")
Ab_woodStarch_entire_tree2_wilcox_pair_comp$month=factor(Ab_woodStarch_entire_tree2_wilcox_pair_comp$month,
                                                         levels = c("Aug19", "Nov19", "Feb20"))
Ab_woodStarch_entire_tree2_wilcox_pair_comp$month2=as.factor(as.numeric(Ab_woodStarch_entire_tree2_wilcox_pair_comp$month))


Ab_woodStarch_entire_tree2_wilcox_pair_comp=Ab_woodStarch_entire_tree2_wilcox_pair_comp%>%
  filter(ID!="ab26820") ### this tree is excluded because it was not measured in February. 
names(Ab_woodStarch_entire_tree2_wilcox_pair_comp)
Ab_woodStarch_entire_tree2_wilcox_pair_comp=Ab_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "starch_change", "rel_starch_change")]
Ab_woodStarch_entire_tree2_wilcox_pair_comp$ID=as.factor(Ab_woodStarch_entire_tree2_wilcox_pair_comp$ID)
str(Ab_woodStarch_entire_tree2_wilcox_pair_comp)

wilcox_test(starch_change~month2, data=Ab_woodStarch_entire_tree2_wilcox_pair_comp, paired = T)
wilcox_test(rel_starch_change~month2, data=Ab_woodStarch_entire_tree2_wilcox_pair_comp, paired = T)


### linear regression models to account for the relationship growth (with time lag) vs starch changes


Ab_woodStarch_entire_tree2=filter(Ab_woodStarch_entire_tree2,month!="May19")
Ab_woodStarch_entire_tree2$month_changes=as.character(Ab_woodStarch_entire_tree2$month)

str(Ab_woodStarch_entire_tree2)
Ab_lm_change_starch_vs_growth=Ab_woodStarch_entire_tree2%>%group_by(month)%>%nest()%>%
  mutate(lm_Ab_schange_starch_vs_growth=map(data, ~lm(rate_growth_summ_three_months_after~starch_change, data = .x)),
         tidied_lm=map(lm_Ab_schange_starch_vs_growth, tidy))

statistics_ab_lm_change_starhc_vs_growth=Ab_lm_change_starch_vs_growth%>%unnest(tidied_lm)
#plot(Ab_lm_annual_growth_starch$lm_Ab_starch_annual_growth[[2]])


month_chages_key=unique(as.character(Ab_woodStarch_entire_tree2$month))
ID_key=c("starch (May19-Aug19)\ngrowth (Aug19-Nov19)", 
         "starch (Aug19-Nov19)\ngrowth (Nov19-Feb20)", 
         "starch (Nov19-Feb20)\n growth (Feb20-May20)")

Ab_woodStarch_entire_tree2$significant_rel=as.character(Ab_woodStarch_entire_tree2$month)
ID_key_signif=c("Not_significant",
                "significant",
                "Not_significant")

library(plyr)
Ab_woodStarch_entire_tree2$month_changes=mapvalues(Ab_woodStarch_entire_tree2$month_changes,
                                                   month_chages_key,
                                                   ID_key)
Ab_woodStarch_entire_tree2$significant_rel=mapvalues(Ab_woodStarch_entire_tree2$significant_rel,
                                                   month_chages_key,
                                                   ID_key_signif)
detach("package:plyr", unload=TRUE)

Ab_woodStarch_entire_tree2$month_changes=factor(Ab_woodStarch_entire_tree2$month_changes,
                                                levels = c("starch (May19-Aug19)\ngrowth (Aug19-Nov19)", 
                                                           "starch (Aug19-Nov19)\ngrowth (Nov19-Feb20)", 
                                                           "starch (Nov19-Feb20)\n growth (Feb20-May20)"))

png("outputs/starch_change_and_growth_O_leu.png",
    width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree2, aes(x=starch_change, y=rate_growth_summ_three_months_after))+
  geom_point(color="turquoise2", size=4, alpha=0.5)+
  geom_smooth(aes(linetype=significant_rel),method = "lm", se=F, color="turquoise2")+
  scale_linetype_manual(values=c("dashed", "solid"))+ 
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=7, xpos = -15,
                             ypos = 0.7, xpos2 = -15,
                             ypos2= 0.65)+
  #geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  xlab("starch change (mg)")+ ylab("accumulated growth (cm)")+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  facet_wrap(~month_changes)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),
        plot.title = element_text(size = rel(2.5), face = "italic"),
        strip.text.x=element_text(size = 20, margin = margin(10), vjust = 1))
dev.off()



##### Calculations for the semi-deciduous/parenchyma-storing species Sacoglotis guianensis####

load("data/Sac_woodStarch.Rda")
############## Calculating mean statistics per individual, month and sampling depth ####################

## mean statistics by tree 

Sac_woodStarch_means_by_tree=Sac_woodStarch%>%group_by(ID, depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

### mean statistics by species 

Sac_woodStarch_means_by_specie=Sac_woodStarch%>%group_by(depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

Sac_woodStarch_means_by_tree$month=factor(Sac_woodStarch_means_by_tree$month, 
                                          levels = c("May19", "Aug19", "Nov19", "Feb20"))
Sac_woodStarch_means_by_specie$month=factor(Sac_woodStarch_means_by_specie$month, 
                                            levels = c("May19", "Aug19", "Nov19", "Feb20"))

png("outputs/radial_starch_content_S-gui.png", width = 900, height = 600)
ggplot()+
  geom_line(data=filter(Sac_woodStarch_means_by_specie,month!="Jan18" && month !="Jul18") , aes(x=depth_mm, y=mean_starch_area, col=month), size=1)+
  #facet_wrap(~month)+
  #geom_line(data=filter(Sac_woodStarch_means_by_tree, month!="Jan18" && month !="Jul18"), aes(x=depth_mm, y=mean_starch_area, col=ID), col="gray95")+
  geom_ribbon(data=filter(Sac_woodStarch_means_by_specie,month!="Jan18" && month !="Jul18"), aes(x=depth_mm, 
                                                                                                 ymin=abs(mean_starch_area-sd_starch_area),
                                                                                                 ymax=mean_starch_area+sd_starch_area,
                                                                                                 fill= month),
              #fill="gray80",
              alpha=0.1)+
  
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("starch %")+
  #xlim(0,110)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

## statistical comparisons at different depth ranges (ANOVA)

Sac_pairwise_comparison=Sac_woodStarch_means_by_tree
Sac_pairwise_comparison=Sac_pairwise_comparison[, c(1,2,3,5)]
Sac_pairwise_comparison$time=as.factor(Sac_pairwise_comparison$month)
Sac_pairwise_comparison=Sac_pairwise_comparison[order(Sac_pairwise_comparison$month),]
rownames(Sac_pairwise_comparison)=NULL

Sac_pairwise_comparison$ID=as.factor(Sac_pairwise_comparison$ID)
Sac_Anva2=Sac_pairwise_comparison%>%group_by(depth_mm)%>%nest()%>%
  mutate(ANOVA=map(data, ~ aov(mean_starch_area~month, data=.x)),
         tukey_comparison=map(ANOVA, TukeyHSD),
         tidied=map(ANOVA, tidy),
         tidied2=map(tukey_comparison, tidy), 
  )

names(Sac_Anva2$tidied)=Sac_Anva2$depth_mm
names(Sac_Anva2$tidied2)=Sac_Anva2$depth_mm
Sac_Anva2$tidied
Sac_Anva2$tidied2
estimates_ANOVA_by_depth=Sac_Anva2%>%unnest(tidied)

#wilcoxon comparison for mean starch area (starch concentrations) at different depth ranges in the wood cores####

#sac_pairwise_comparison2=filter(sac_pairwise_comparison, ID!="sac34288") ## We had to remove this tree because we lost one month of measurements for this individuals and Wilcoxon  test requires identical sample size for each group. 
# Sac_pairwise_comparison2=filter(Sac_woodStarch_means_by_tree, ID!="Sac1054")
# Sac_pairwise_comparison2=filter(Sac_woodStarch_means_by_tree, depth_mm2)
# 
# Sac_pairwise_comparison2$ID=factor(Sac_pairwise_comparison2$ID)
# Sac_pairwise_comparison2$depth_mm2=as.factor(Sac_pairwise_comparison2$depth_mm)
# table(Sac_pairwise_comparison2$ID)
# with(Sac_pairwise_comparison2, table(ID, depth_mm2))
# 
# str(Sac_pairwise_comparison2)
# sac_Anva3=Sac_pairwise_comparison2%>%group_by(depth_mm2)%>%
#   wilcox_test(mean_starch_area~month, data=., paired = T)

#####

ggplot(Sac_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_boxplot(data=Sac_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_point(aes(col=ID))+
  geom_point(data=Sac_woodStarch_means_by_specie, aes(x=month, y=mean_starch_area), col="black")+
  facet_wrap(~depth_mm)


### estimate the maximum and relative seasonal range in starch mass during the year of measurements in each wood depth

Sac_woodStarch_radial_activity=Sac_woodStarch_means_by_tree%>%group_by(ID, depth_mm)%>%
  summarise(activ_mean=max(mean_starch_area)-min(mean_starch_area),
            rel_activ_mean=(max(mean_starch_area)-min(mean_starch_area))/max(mean_starch_area)*100,
            activ_median=max(median_starch_are)-min(median_starch_are),
            activ_IQR=max(Int_quart_range)-min(Int_quart_range),
            max_starch=max(mean_starch_area),
            max_month=.data$month[max(which(.data$mean_starch_area==max(.data$mean_starch_area)))],
            min_starch=min(mean_starch_area),
            min_month=.data$month[min(which(.data$mean_starch_area==min(.data$mean_starch_area)))]
  )

Sac_woodStarch_radial_activity=na.omit(Sac_woodStarch_radial_activity)

Sac_woodStarch_radial_Activity_by_species=Sac_woodStarch_radial_activity%>%group_by(depth_mm)%>%
  summarise(mean_active_mean=mean(activ_mean, na.rm=T),
            rel_mean_active_mean=mean(rel_activ_mean, na.rm=T),
            SD_rel_mean_active_mean=sd(rel_activ_mean, na.rm=T),
            mean_active_median=mean(activ_median),
            mean_active_IQR=mean(activ_IQR),
            SD_active_mean=sd(activ_mean),
            SD_active_median=sd(activ_median),
            Sd_active_IQR=sd(activ_IQR), 
            max=max(max_starch),
            mean_max=mean(max_starch),
            max_month_mode=names(which.max(table(max_month))),
            min=min(min_starch),
            mean_min=mean(min_starch),
            min_month_mode=names(which.max(table(min_month)))
  )

## ANOVA to compare changes in the maximum starch seasonal amplitud at different wood depths 

# check outliers 
Sac_woodStarch_radial_activity%>%group_by(depth_mm)%>%identify_outliers(activ_mean)

#normality assumption 
Sac_woodStarch_radial_activity%>%filter(depth_mm<115)%>%group_by(depth_mm)%>%shapiro_test(activ_mean)
ggqqplot(Sac_woodStarch_radial_activity, "activ_mean", facet.by = "depth_mm")

str(Sac_woodStarch_radial_activity)
seasonal_amplitud_anova=aov(activ_mean~as.factor(depth_mm), data=Sac_woodStarch_radial_activity)
summary(seasonal_amplitud_anova)
#plot(seasonal_amplitud_anova)
TukeyHSD(seasonal_amplitud_anova)

png("outputs/metabolic_activity_Sac.png", width = 800, height=600)
ggplot(data=Sac_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=mean_active_mean))+
  geom_line(data=Sac_woodStarch_radial_activity, aes(x=depth_mm, y=activ_mean, factor=ID), col="gray95")+
  geom_line( col="turquoise3")+
  geom_ribbon( aes(ymin=mean_active_mean-SD_active_mean, ymax=mean_active_mean+SD_active_mean,
                   fill= "var"),
               fill="turquoise3",
               alpha=0.3)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("seasonal amplitud starch %")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

png("outputs/relative_sorage_change_Sac.png", width = 900, height = 600)
ggplot(data=Sac_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=rel_mean_active_mean))+
  geom_line(data=Sac_woodStarch_radial_activity, aes(x=depth_mm, y=rel_activ_mean, factor=ID), col="gray95")+
  geom_line( col="turquoise3")+
  geom_ribbon( aes(ymin=rel_mean_active_mean-SD_rel_mean_active_mean, ymax=rel_mean_active_mean+SD_rel_mean_active_mean,
                   fill= "var"),
               fill="turquoise3",
               alpha=0.3)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("relative seasonal amplitud")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

## mean content of starch in the entire core 
## wood density 0.81 g/cm3 was taken from tprevious records of wood densitty 

Sac_woodStarch_means_by_tree$segment_volume_cm3=pi*(0.5/2)^2*0.5 ### this in in cm3

Sac_woodStarch_means_by_tree$grams_wood=Sac_woodStarch_means_by_tree$segment_volume_cm3*0.8
Sac_woodStarch_means_by_tree$grams_starch=Sac_woodStarch_means_by_tree$grams_wood*Sac_woodStarch_means_by_tree$mean_starch_area/100
Sac_woodStarch_means_by_tree$sd_grams_starch=Sac_woodStarch_means_by_tree$grams_wood*Sac_woodStarch_means_by_tree$sd_starch_area/100

Sac_woodStarch_means_entire_tree=Sac_woodStarch_means_by_tree%>%group_by(month, ID)%>%
  summarise(mean_starch_percentage=mean(mean_starch_area),
            median_starch_percentage=mean(median_starch_are),
            mean_IQR=mean(Int_quart_range),
            starch_woodcore_content_by_mean=sum(mean_starch_area),
            starch_woodcore_content_by_median=sum(median_starch_are),
            starch_woodcore_content_grams=sum(grams_starch),
            #starch_woodcore_content_mg=sum(mg_starch),
            sd_starch_woodcore_content_gr=sd(grams_starch),
            sd2_starch_woodcore_content_gr=mean(sd_grams_starch))


Sac_woodStarch_means_entire_tree$starch_woodcore_content_mg=Sac_woodStarch_means_entire_tree$starch_woodcore_content_grams*1000
Sac_woodStarch_means_entire_tree$sd_starch_woodcore_content_mg=Sac_woodStarch_means_entire_tree$sd2_starch_woodcore_content_gr*1000

Sac_woodStarch_means_entire_tree_by_specie=Sac_woodStarch_means_entire_tree%>%group_by(month)%>%
  summarise(mean_starch_percentage=mean(mean_starch_percentage),
            sd_mean_starch_percentage=sd(mean_starch_percentage),
            median_starch_percentage=mean(median_starch_percentage),
            sd_median_starch_percentage=sd(median_starch_percentage),
            mean_IQR=mean(mean_IQR), 
            mean_starch_wood_content_by_mean=mean(starch_woodcore_content_by_mean),
            sd_starch_wood_content_by_mean=sd(starch_woodcore_content_by_mean),
            mean_starch_wood_content_by_median=mean(starch_woodcore_content_by_median),
            sd_starch_wood_content_by_medial=sd(starch_woodcore_content_by_median),
            mean_starch_wood_content_grams=mean(starch_woodcore_content_grams),
            sd_mean_starch_wood_content_grams=sd(starch_woodcore_content_grams),
            mean_starch_wood_content_mg=mean(starch_woodcore_content_mg),
            sd_mean_starch_wood_content_mg=sd(starch_woodcore_content_mg),
            sd2_mean_starch_wood_content_mg=sum(sd_starch_woodcore_content_mg)
  )
Sac_woodStarch_means_entire_tree$starch_woodcore_content_mg=Sac_woodStarch_means_entire_tree$starch_woodcore_content_grams*1000

png("outputs/Sac_starch_wood_content_comp_month.png", width = 800, height = 600)
ggplot(Sac_woodStarch_means_entire_tree, aes(x=month, y=starch_woodcore_content_mg))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_point(data=Sac_woodStarch_means_entire_tree_by_specie, aes(x=month, y=mean_starch_wood_content_mg),  col="black", pch=17, size=3)+
  xlab("")+
  ylab("Starch woodcore content (mg)")+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


### pairwised comparison with paired samples for the content of sugar in the entire tree. The differences between months are lost mainly because the influence of the lack 
## of seasonality deeper in the stemwood or the contraty seasonality, and then all the averaging and the integration. 

Sac_pairwise_comparison_entire_tree=Sac_woodStarch_means_entire_tree
Sac_pairwise_comparison_entire_tree$month2=as.factor(as.numeric(Sac_pairwise_comparison_entire_tree$month))
#Sac_pairwise_comparison_entire_tree=Sac_pairwise_comparison_entire_tree[-c(which(Sac_pairwise_comparison_entire_tree$ID=="Sac26820")),]
Sac_pairwise_comparison_entire_tree=Sac_pairwise_comparison_entire_tree[-c(which(Sac_pairwise_comparison_entire_tree$ID=="Sac1054")),]
names(Sac_pairwise_comparison_entire_tree)
Sac_pairwise_comparison_entire_tree=Sac_pairwise_comparison_entire_tree[ ,c(2,11,13)]
Sac_pairwise_comparison_entire_tree$ID=as.factor(Sac_pairwise_comparison_entire_tree$ID)
str(Sac_pairwise_comparison_entire_tree)

Sac_pwc_starch_mg <- Sac_pairwise_comparison_entire_tree %>%
  pairwise_t_test(starch_woodcore_content_mg ~ month2, paired = T,
                  p.adjust.method = "BH"
  )
Sac_pwc_starch_mg


wilcox_test(starch_woodcore_content_mg~month2, data=Sac_pairwise_comparison_entire_tree, paired = T)


#### add growth data to Sac_woodStarch_content_woodcore and compare with total starch mass in the entire woodcore. 

monthly_growth2_summarized=monthly_growth2%>%filter(date3%in%c(8:20))%>%group_by(ID)%>%
  summarise(annual_growth_2019=sum(rate_growth),
            mean_monthly_growth_2019=mean(rate_growth),
            max_monthly_growth_2019=max(rate_growth),
            min_monthly_growth_2019=min(rate_growth)
  )


Sac_woodStarch_entire_tree=merge(Sac_woodStarch_means_entire_tree, monthly_growth2, by=c("ID", "month"))
Sac_woodStarch_entire_tree=merge(Sac_woodStarch_entire_tree, monthly_growth2_summarized, by=c("ID"))

str(Sac_woodStarch_entire_tree)
Sac_lm_annual_growth_starch=Sac_woodStarch_entire_tree%>%group_by(month)%>%nest()%>%
  mutate(lm_Sac_starch_annual_growth=map(data, ~lm(annual_growth_2019~starch_woodcore_content_mg, data = .x)),
         tidied_lm=map(lm_Sac_starch_annual_growth, tidy))

statistics_Sac_lm=Sac_lm_annual_growth_starch%>%unnest(tidied_lm)
#plot(Sac_lm_annual_growth_starch$lm_Sac_starch_annual_growth[[4]])


Sac_woodStarch_entire_tree$significant_rel=as.character(Sac_woodStarch_entire_tree$month)
key_months_1=unique(as.character(Sac_woodStarch_entire_tree$month))
ID_key_signif2=c("Not_significant",
                 "Not_significant",
                 "significant",
                 "significant")

library(plyr)

Sac_woodStarch_entire_tree$significant_rel=mapvalues(Sac_woodStarch_entire_tree$significant_rel,
                                                    key_months_1,
                                                    ID_key_signif2)
detach("package:plyr", unload=TRUE)

png("outputs/annual_growth_vs_starch_content_Sac.png",
    width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree, aes(x=starch_woodcore_content_mg, y=annual_growth_2019))+
  geom_point(color="turquoise2")+
  geom_smooth(aes(linetype=significant_rel), method = "lm", se=F, color="turquoise2")+
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=8, xpos = 1,
                             ypos = 2.4, xpos2 = 1,
                             ypos2= 2.1)+
  xlab("starch content (mg)")+ ylab("annual growth 2019 (cm)")+
  ylim(0, 2.5)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  facet_wrap(~month)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


# Calculating the changes in the total starch mass in the entire wood core between seasons. 

Sac_woodStarch_entire_tree=Sac_woodStarch_entire_tree[order(Sac_woodStarch_entire_tree$ID, Sac_woodStarch_entire_tree$date3),]
Sac_woodStarch_entire_tree$starch_change=NA
Sac_woodStarch_entire_tree$starch_change_percent=NA
Sac_woodStarch_entire_tree$rel_starch_change=NA

Sac_woodStarch_entire_tree2=Sac_woodStarch_entire_tree%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 2:nrow(x)){
      x$starch_change[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1]
      x$starch_change_percent[i]=x$mean_starch_percentage[i]-x$mean_starch_percentage[i-1]
      x$rel_starch_change[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1])/x$starch_woodcore_content_mg[i]
    }
    return(x)
  }
  )%>%ungroup()


## pairwise comparison between changes in the starch content between seasons 

png("outputs/Sac_starch_mass_change_comp_month.png", width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (mg)")+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


png("outputs/Sac_starch_mass_change_rel_comp_month.png", width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=rel_starch_change))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  #ylim(-4, 2)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

png("outputs/Sac_starch_mass_change_percentage_month.png", width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change_percent))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  #ylim(-4, 2)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


### wilcoxon comparisson with paired samples to compare the changes in total starch mass in the wood core between different months. 
Sac_woodStarch_entire_tree2_wilcox_pair_comp=filter(Sac_woodStarch_entire_tree2, month!="May19")
Sac_woodStarch_entire_tree2_wilcox_pair_comp$month=factor(Sac_woodStarch_entire_tree2_wilcox_pair_comp$month,
                                                          levels = c("Aug19", "Nov19", "Feb20"))
Sac_woodStarch_entire_tree2_wilcox_pair_comp$month2=as.factor(as.numeric(Sac_woodStarch_entire_tree2_wilcox_pair_comp$month))
names(Sac_woodStarch_entire_tree2_wilcox_pair_comp)
Sac_woodStarch_entire_tree3_wilcox_pair_comp=Sac_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "starch_change")]
Sac_woodStarch_entire_tree4_wilcox_pair_comp=Sac_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "rel_starch_change")]
Sac_woodStarch_entire_tree3_wilcox_pair_comp$ID=as.factor(Sac_woodStarch_entire_tree3_wilcox_pair_comp$ID)
Sac_woodStarch_entire_tree4_wilcox_pair_comp$ID=as.factor(Sac_woodStarch_entire_tree4_wilcox_pair_comp$ID)
str(Sac_woodStarch_entire_tree3_wilcox_pair_comp)
str(Sac_woodStarch_entire_tree4_wilcox_pair_comp)

wilcox_test(starch_change~month2, data=Sac_woodStarch_entire_tree3_wilcox_pair_comp, paired = T)

wilcox_test(rel_starch_change~month2, data=Sac_woodStarch_entire_tree4_wilcox_pair_comp, paired = T)


### linear models for the relationship between accumulated three months of growth (with a time lag) vs starch changes

Sac_woodStarch_entire_tree2=filter(Sac_woodStarch_entire_tree2,month!="May19")
Sac_woodStarch_entire_tree2$month_changes=as.character(Sac_woodStarch_entire_tree2$month)


str(Sac_woodStarch_entire_tree2)
Sac_lm_change_starch_vs_growth=Sac_woodStarch_entire_tree2%>%group_by(month)%>%nest()%>%
  mutate(lm_Sac_schange_starch_vs_growth=map(data, ~lm(rate_growth_summ_three_months_after~starch_change, data = .x)),
         tidied_lm=map(lm_Sac_schange_starch_vs_growth, tidy))

statistics_Sac_lm_change_starhc_vs_growth=Sac_lm_change_starch_vs_growth%>%unnest(tidied_lm)
#plot(Sac_lm_annual_growth_starch$lm_Sac_starch_annual_growth[[2]])


month_chages_key=unique(as.character(Sac_woodStarch_entire_tree2$month))
ID_key=c("starch (May19-Aug19)\ngrowth (Aug19-Nov19)", 
         "starch (Aug19-Nov19)\ngrowth (Nov19-Feb20)", 
         "starch (Nov19-Feb20)\n growth (Feb20-May20)")


library(plyr)
Sac_woodStarch_entire_tree2$month_changes=mapvalues(Sac_woodStarch_entire_tree2$month_changes,
                                                   month_chages_key,
                                                   ID_key)
detach("package:plyr", unload=TRUE)

Sac_woodStarch_entire_tree2$month_changes=factor(Sac_woodStarch_entire_tree2$month_changes,
                                                levels = c("starch (May19-Aug19)\ngrowth (Aug19-Nov19)", 
                                                           "starch (Aug19-Nov19)\ngrowth (Nov19-Feb20)", 
                                                           "starch (Nov19-Feb20)\n growth (Feb20-May20)"))

png("outputs/three_month_growth_after_vs_starch_change_S_gui.png",
    width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree2, aes(x=starch_change, y=rate_growth_summ_three_months_after))+
  geom_point(color="turquoise2", size=4, alpha=0.5)+
  #geom_text(aes(label=ID))+
  geom_smooth(linetype=c("dashed"), method = "lm", se=F, color="turquoise")+
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=7, xpos = -75,
                             ypos = 0.66, xpos2 = -75,
                             ypos2= 0.63)+
  xlab("starch change (mg)")+ylab("accumulated growth (cm)")+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  ylim(c(0.0, 0.7))+
  geom_vline(xintercept = 0, linetype="dashed")+
  facet_wrap(~month_changes)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3)),
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),
        legend.position="none",
        plot.title = element_text(size = rel(2.5), face = "italic"),
        strip.text.x = element_text(size= 20, margin = margin(10), vjust = 1),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()



############################################################################################################################################################


##### Calculations for the semi-deciduous/fiber storing species Dacryodes microcarpa####

load("data/Dac_woodStarch.Rda")

############## Calculating mean statistics per individual, month and sampling depth ####################

## mean statistics by tree 

Dac_woodStarch_means_by_tree=Dac_woodStarch%>%group_by(ID, depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

### mean statistics by species 

Dac_woodStarch_means_by_specie=Dac_woodStarch%>%group_by(depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

Dac_woodStarch_means_by_tree$month=factor(Dac_woodStarch_means_by_tree$month, 
                                          levels = c("May19", "Aug19", "Nov19", "Feb20"))
Dac_woodStarch_means_by_specie$month=factor(Dac_woodStarch_means_by_specie$month, 
                                            levels = c("May19", "Aug19", "Nov19", "Feb20"))

png("outputs/radial_starch_content_Dmic.png", width = 900, height = 600)
ggplot()+
  geom_line(data=filter(Dac_woodStarch_means_by_specie,month!="Jan18" && month !="Jul18") , aes(x=depth_mm, y=mean_starch_area, col=month), size=1)+
  #facet_wrap(~month)+
  #geom_line(data=filter(Dac_woodStarch_means_by_tree, month!="Jan18" && month !="Jul18"), aes(x=depth_mm, y=mean_starch_area, col=ID), col="gray95")+
  geom_ribbon(data=filter(Dac_woodStarch_means_by_specie,month!="Jan18" && month !="Jul18"), aes(x=depth_mm, 
                                                                                                 ymin=mean_starch_area-sd_starch_area,
                                                                                                 ymax=mean_starch_area+sd_starch_area,
                                                                                                 fill= month),
              #fill="gray80",
              alpha=0.1)+
  
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  xlab("depth (mm)")+ylab("starch %")+
  #xlim(0,110)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

## statistical comparisons at different depth ranges (ANOVA)

Dac_pairwise_comparison=Dac_woodStarch_means_by_tree
Dac_pairwise_comparison=Dac_pairwise_comparison[, c(1,2,3,5)]
Dac_pairwise_comparison$time=as.factor(Dac_pairwise_comparison$month)
Dac_pairwise_comparison=Dac_pairwise_comparison[order(Dac_pairwise_comparison$month),]
rownames(Dac_pairwise_comparison)=NULL

Dac_pairwise_comparison$ID=as.factor(Dac_pairwise_comparison$ID)
Dac_Anva2=Dac_pairwise_comparison%>%group_by(depth_mm)%>%nest()%>%
  mutate(ANOVA=map(data, ~ aov(mean_starch_area~month, data=.x)),
         tukey_comparison=map(ANOVA, TukeyHSD),
         tidied=map(ANOVA, tidy),
         tidied2=map(tukey_comparison, tidy), 
  )

names(Dac_Anva2$tidied)=Dac_Anva2$depth_mm
names(Dac_Anva2$tidied2)=Dac_Anva2$depth_mm
Dac_Anva2$tidied
Dac_Anva2$tidied2
estimates_ANOVA_by_depth=Dac_Anva2%>%unnest(tidied)

table(Dac_pairwise_comparison$ID)

#wilcoxon comparison for mean starch area (starch concentrations) at different depth ranges in the wood cores
Dac_pairwise_comparison2=filter(Dac_pairwise_comparison, ID!="Dac32001") ## We had to remove this tree because we lost one month of measurements for this individuals and Wilcoxon  test requires identical sample size for each group. 
Dac_pairwise_comparison2$ID=factor(Dac_pairwise_comparison2$ID)
Dac_pairwise_comparison2$depth_mm2=as.factor(Dac_pairwise_comparison2$depth_mm)
table(Dac_pairwise_comparison2$ID)

str(Dac_pairwise_comparison2)
Dac_Anva3=Dac_pairwise_comparison2%>%group_by(depth_mm2)%>%
  wilcox_test(mean_starch_area~month, data=., paired = T)


ggplot(Dac_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_boxplot(data=Dac_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_point(aes(col=ID))+
  geom_point(data=Dac_woodStarch_means_by_specie, aes(x=month, y=mean_starch_area), col="black")+
  facet_wrap(~depth_mm)


### estimate the maximum and relative seasonal range in starch mass during the year of measurements in each wood depth

Dac_woodStarch_radial_activity=Dac_woodStarch_means_by_tree%>%group_by(ID, depth_mm)%>%
  summarise(activ_mean=max(mean_starch_area)-min(mean_starch_area),
            rel_activ_mean=(max(mean_starch_area)-min(mean_starch_area))/max(mean_starch_area)*100,
            activ_median=max(median_starch_are)-min(median_starch_are),
            activ_IQR=max(Int_quart_range)-min(Int_quart_range),
            max_starch=max(mean_starch_area),
            max_month=.data$month[max(which(.data$mean_starch_area==max(.data$mean_starch_area)))],
            min_starch=min(mean_starch_area),
            min_month=.data$month[min(which(.data$mean_starch_area==min(.data$mean_starch_area)))]
  )

Dac_woodStarch_radial_activity=na.omit(Dac_woodStarch_radial_activity)

Dac_woodStarch_radial_Activity_by_species=Dac_woodStarch_radial_activity%>%group_by(depth_mm)%>%
  summarise(mean_active_mean=mean(activ_mean, na.rm=T),
            rel_mean_active_mean=mean(rel_activ_mean, na.rm=T),
            SD_rel_mean_active_mean=sd(rel_activ_mean, na.rm=T),
            mean_active_median=mean(activ_median),
            mean_active_IQR=mean(activ_IQR),
            SD_active_mean=sd(activ_mean),
            SD_active_median=sd(activ_median),
            Sd_active_IQR=sd(activ_IQR), 
            max=max(max_starch),
            mean_max=mean(max_starch),
            max_month_mode=names(which.max(table(max_month))),
            min=min(min_starch),
            mean_min=mean(min_starch),
            min_month_mode=names(which.max(table(min_month)))
  )


## ANOVA to compare changes in the starch change at different wood depths 

# check outliers 
Dac_woodStarch_radial_activity%>%group_by(depth_mm)%>%identify_outliers(activ_mean)

#normality assumption 
Dac_woodStarch_radial_activity%>%filter(depth_mm<65)%>%group_by(depth_mm)%>%shapiro_test(activ_mean)
ggqqplot(Dac_woodStarch_radial_activity, "activ_mean", facet.by = "depth_mm")

str(Dac_woodStarch_radial_activity)
seasonal_amplitud_anova=aov(activ_mean~as.factor(depth_mm), data=Dac_woodStarch_radial_activity)
summary(seasonal_amplitud_anova)
#plot(seasonal_amplitud_anova)
TukeyHSD(seasonal_amplitud_anova)


png("outputs/metabolic_activity_Dac.png", width = 800, height=600)
ggplot(data=Dac_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=mean_active_mean))+
  geom_line(data=Dac_woodStarch_radial_activity, aes(x=depth_mm, y=activ_mean, factor=ID), col="gray95")+
  geom_line( col="orangered2")+
  geom_ribbon( aes(ymin=mean_active_mean-SD_active_mean, ymax=mean_active_mean+SD_active_mean,
                   fill= "var"),
               fill="orangered2",
               alpha=0.3)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  xlab("depth (mm)")+ylab("seasonal amplitud starch %")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

png("outputs/relative_sorage_change_Dac.png", width = 900, height = 600)
ggplot(data=Dac_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=rel_mean_active_mean))+
  geom_line(data=Dac_woodStarch_radial_activity, aes(x=depth_mm, y=rel_activ_mean, factor=ID), col="gray95")+
  geom_line( col="orangered2")+
  geom_ribbon( aes(ymin=rel_mean_active_mean-SD_rel_mean_active_mean, ymax=rel_mean_active_mean+SD_rel_mean_active_mean,
                   fill= "var"),
               fill="orangered2",
               alpha=0.3)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  xlab("depth (mm)")+ylab("relative seasonal amplitud (%)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

## mean content of starch in the entire core 
## wood density 0.51 g/cm3 was taken from previous measurements of trees from this species in the study area

Dac_woodStarch_means_by_tree$segment_volume_cm3=pi*(0.5/2)^2*0.5 ### this in in cm3

Dac_woodStarch_means_by_tree$grams_wood=Dac_woodStarch_means_by_tree$segment_volume_cm3*0.51
Dac_woodStarch_means_by_tree$grams_starch=Dac_woodStarch_means_by_tree$grams_wood*Dac_woodStarch_means_by_tree$mean_starch_area/100
Dac_woodStarch_means_by_tree$sd_grams_starch=Dac_woodStarch_means_by_tree$grams_wood*Dac_woodStarch_means_by_tree$sd_starch_area/100


Dac_woodStarch_means_entire_tree=Dac_woodStarch_means_by_tree%>%group_by(month, ID)%>%
  summarise(mean_starch_percentage=mean(mean_starch_area),
            median_starch_percentage=median(median_starch_are),
            mean_IQR=mean(Int_quart_range),
            starch_woodcore_content_by_mean=sum(mean_starch_area),
            starch_woodcore_content_by_median=sum(median_starch_are),
            starch_woodcore_content_grams=sum(grams_starch),
            sd_starch_woodcore_content_gr=sd(grams_starch),
            sd2_starch_woodcore_content_gr=mean(sd_grams_starch))

Dac_woodStarch_means_entire_tree$starch_woodcore_content_mg=Dac_woodStarch_means_entire_tree$starch_woodcore_content_grams*1000
Dac_woodStarch_means_entire_tree$sd_starch_woodcore_content_mg=Dac_woodStarch_means_entire_tree$sd2_starch_woodcore_content_gr*1000



Dac_woodStarch_means_entire_tree_by_specie=Dac_woodStarch_means_entire_tree%>%group_by(month)%>%
  summarise(mean_starch_percentage=mean(mean_starch_percentage),
            sd_mean_starch_percentage=sd(mean_starch_percentage),
            median_starch_percentage=mean(median_starch_percentage),
            sd_median_starch_percentage=sd(median_starch_percentage),
            mean_IQR=mean(mean_IQR), 
            mean_starch_wood_content_by_mean=mean(starch_woodcore_content_by_mean),
            sd_starch_wood_content_by_mean=sd(starch_woodcore_content_by_mean),
            mean_starch_wood_content_by_median=mean(starch_woodcore_content_by_median),
            sd_starch_wood_content_by_medial=sd(starch_woodcore_content_by_median),
            mean_starch_wood_content_grams=mean(starch_woodcore_content_grams),
            sd_mean_starch_wood_content_grams=sd(starch_woodcore_content_grams),
            mean_starch_wood_content_mg=mean(starch_woodcore_content_mg),
            sd_mean_starch_wood_content_mg=sd(starch_woodcore_content_mg),
            sd2_mean_starch_wood_content_mg=sum(sd_starch_woodcore_content_mg)
  )


png("outputs/Dac_starch_wood_content_comp_month.png", width = 800, height = 600)
ggplot(Dac_woodStarch_means_entire_tree, aes(x=month, y=starch_woodcore_content_mg))+geom_boxplot(fill="orangered2",  alpha=0.3)+
  geom_point(aes(col=ID))+
  geom_point(data=Dac_woodStarch_means_entire_tree_by_specie, aes(x=month, y=mean_starch_wood_content_mg),  col="black", pch=17, size=3)+
  xlab("")+
  ylab("Starch woodcore content (mg)")+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


### pairwise comparison with paired samples for the content of sugar in the entire tree. 

Dac_pairwise_comparison_entire_tree=Dac_woodStarch_means_entire_tree
Dac_pairwise_comparison_entire_tree$month2=as.factor(as.numeric(Dac_pairwise_comparison_entire_tree$month))
Dac_pairwise_comparison_entire_tree=Dac_pairwise_comparison_entire_tree[-c(which(Dac_pairwise_comparison_entire_tree$ID=="Dac32001")),]
names(Dac_pairwise_comparison_entire_tree)
Dac_pairwise_comparison_entire_tree=Dac_pairwise_comparison_entire_tree[ ,c(2,11,13)]
Dac_pairwise_comparison_entire_tree$ID=as.factor(Dac_pairwise_comparison_entire_tree$ID)
str(Dac_pairwise_comparison_entire_tree)


Dac_pwc_starch_mg <- Dac_pairwise_comparison_entire_tree %>%
  pairwise_t_test(starch_woodcore_content_mg ~ month2, paired = T,
                  p.adjust.method = "BH"
  )
Dac_pwc_starch_mg


wilcox_test(starch_woodcore_content_mg~month2, data=Dac_pairwise_comparison_entire_tree, paired = T)


#### add growth data to Dac_woodStarch_content_woodcore and compare with total starch, mean starch, starch in the first two cm... 

monthly_growth2_summarized=monthly_growth2%>%filter(date3%in%c(8:20))%>%group_by(ID)%>%
  summarise(annual_growth_2019=sum(rate_growth),
            mean_monthly_growth_2019=mean(rate_growth),
            max_monthly_growth_2019=max(rate_growth),
            min_monthly_growth_2019=min(rate_growth)
  )


Dac_woodStarch_entire_tree=merge(Dac_woodStarch_means_entire_tree, monthly_growth2, by=c("ID", "month"))
Dac_woodStarch_entire_tree=merge(Dac_woodStarch_entire_tree, monthly_growth2_summarized, by=c("ID"))

str(Dac_woodStarch_entire_tree)
Dac_lm_annual_growth_starch=Dac_woodStarch_entire_tree%>%group_by(month)%>%nest()%>%
  mutate(lm_Dac_starch_annual_growth=map(data, ~lm(annual_growth_2019~starch_woodcore_content_mg, data = .x)),
         tidied_lm=map(lm_Dac_starch_annual_growth, tidy))

statistics_Dac_lm=Dac_lm_annual_growth_starch%>%unnest(tidied_lm)
#plot(Dac_lm_annual_growth_starch$lm_Dac_starch_annual_growth[[4]])

Dac_woodStarch_entire_tree$significant_rel=as.character(Dac_woodStarch_entire_tree$month)
key_months_1=unique(as.character(Dac_woodStarch_entire_tree$month))
ID_key_signif2=c("Not_significant",
                 "Not_significant",
                 "Not_significant",
                 "Not_significant")

library(plyr)

Dac_woodStarch_entire_tree$significant_rel=mapvalues(Dac_woodStarch_entire_tree$significant_rel,
                                                     key_months_1,
                                                     ID_key_signif2)
detach("package:plyr", unload=TRUE)

png("outputs/annual_growth_vs_starch_content_Dac.png",
    width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree, aes(x=starch_woodcore_content_mg, y=annual_growth_2019))+
  geom_point(color="orangered2")+
  geom_smooth(method = "lm", se=F, color="orangered2", linetype="dashed")+
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=8, xpos = 1,
                             ypos = 2.4, xpos2 = 1,
                             ypos2= 2.1)+
  xlab("starch content (mg)")+ ylab("annual growth 2019 (cm)")+
  ylim(0, 2.5)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  facet_wrap(~month)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


####calculating the changes in total starch mass in the entire wood core between seasons. 

Dac_woodStarch_entire_tree=Dac_woodStarch_entire_tree[order(Dac_woodStarch_entire_tree$ID, Dac_woodStarch_entire_tree$date3),]
Dac_woodStarch_entire_tree$starch_change=NA
Dac_woodStarch_entire_tree$starch_change_percent=NA
Dac_woodStarch_entire_tree$rel_starch_change=NA

Dac_woodStarch_entire_tree2=Dac_woodStarch_entire_tree%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 2:nrow(x)){
      x$starch_change[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1]
      x$starch_change_percent[i]=x$mean_starch_percentage[i]-x$mean_starch_percentage[i-1]
      x$rel_starch_change[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1])/x$starch_woodcore_content_mg[i]
    }
    return(x)
  }
  )%>%ungroup()



png("outputs/Dac_starch_mass_change_comp_month.png", width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change))+geom_boxplot(fill="orangered2",  alpha=0.3)+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (mg)")+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


png("outputs/Dac_starch_mass_change_rel_comp_month.png", width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=rel_starch_change))+geom_boxplot(fill="orangered2",  alpha=0.3)+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  #ylim(-4, 2)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

png("outputs/Dac_starch_mass_change_percentage_month.png", width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change_percent))+geom_boxplot(fill="orangered2",  alpha=0.3)+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  #ylim(-4, 2)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

## pairwise comparison using wilxocon singed-rank test of the changes in total starch content between seasons 

Dac_woodStarch_entire_tree2_wilcox_pair_comp=filter(Dac_woodStarch_entire_tree2, month!="May19")
Dac_woodStarch_entire_tree2_wilcox_pair_comp$month=factor(Dac_woodStarch_entire_tree2_wilcox_pair_comp$month,
                                                          levels = c("Aug19", "Nov19", "Feb20"))
Dac_woodStarch_entire_tree2_wilcox_pair_comp$month2=as.factor(as.numeric(Dac_woodStarch_entire_tree2_wilcox_pair_comp$month))
names(Dac_woodStarch_entire_tree2_wilcox_pair_comp)
Dac_woodStarch_entire_tree3_wilcox_pair_comp=Dac_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "starch_change")]
Dac_woodStarch_entire_tree4_wilcox_pair_comp=Dac_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "rel_starch_change")]
Dac_woodStarch_entire_tree3_wilcox_pair_comp$ID=as.factor(Dac_woodStarch_entire_tree3_wilcox_pair_comp$ID)
Dac_woodStarch_entire_tree4_wilcox_pair_comp$ID=as.factor(Dac_woodStarch_entire_tree4_wilcox_pair_comp$ID)
str(Dac_woodStarch_entire_tree3_wilcox_pair_comp)
str(Dac_woodStarch_entire_tree4_wilcox_pair_comp)

wilcox_test(starch_change~month2, data=Dac_woodStarch_entire_tree3_wilcox_pair_comp, paired = T)

wilcox_test(rel_starch_change~month2, data=Dac_woodStarch_entire_tree4_wilcox_pair_comp, paired = T)


### linear models for the relationship growth vs starch changes

Dac_woodStarch_entire_tree2=filter(Dac_woodStarch_entire_tree2,month!="May19")
Dac_woodStarch_entire_tree2$month_changes=as.character(Dac_woodStarch_entire_tree2$month)

str(Dac_woodStarch_entire_tree2)
Dac_lm_change_starch_vs_growth=Dac_woodStarch_entire_tree2%>%group_by(month)%>%nest()%>%
  mutate(lm_Dac_schange_starch_vs_growth=map(data, ~lm(rate_growth_summ_three_months_prior~starch_change, data = .x)),
         tidied_lm=map(lm_Dac_schange_starch_vs_growth, tidy))

statistics_Dac_lm_change_starhc_vs_growth=Dac_lm_change_starch_vs_growth%>%unnest(tidied_lm)
#plot(Dac_lm_annual_growth_starch$lm_Dac_starch_annual_growth[[2]])

month_chages_key=unique(as.character(Dac_woodStarch_entire_tree2$month))
ID_key=c("starch (May19-Aug19)\ngrowth (May19-Aug19)",
         "starch (Aug19-Nov19)\ngrowth (Aug19-Nov20)", 
         "starch (Nov19-Feb20)\n growth (Nov19-Feb20)"
         )

Dac_woodStarch_entire_tree2$significant_rel=as.character(Dac_woodStarch_entire_tree2$month)
ID_key_dac_signif=c("Not_significant",
                    "Not_significant",
                    "significant")

library(plyr)
Dac_woodStarch_entire_tree2$month_changes=mapvalues(Dac_woodStarch_entire_tree2$month_changes,
                                                    month_chages_key,
                                                    ID_key)

Dac_woodStarch_entire_tree2$significant_rel=mapvalues(Dac_woodStarch_entire_tree2$significant_rel,
                                                    month_chages_key,
                                                    ID_key_dac_signif)
detach("package:plyr", unload=TRUE)

Dac_woodStarch_entire_tree2$month_changes=factor(Dac_woodStarch_entire_tree2$month_changes,
                                                 levels = c("starch (May19-Aug19)\ngrowth (May19-Aug19)",
                                                            "starch (Aug19-Nov19)\ngrowth (Aug19-Nov20)", 
                                                            "starch (Nov19-Feb20)\n growth (Nov19-Feb20)"
                                                            ))

png("outputs/three_month_growth_after_vs_starch_change_Dac.png",
    width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree2, aes(x=starch_change, y=rate_growth_summ_three_months_prior))+
  geom_point(color="orangered2", size=4, alpha=0.5)+
  geom_smooth(aes(linetype=significant_rel), method = "lm", se=F, color="orangered")+
  scale_linetype_manual(values=c("dashed", "solid")) +
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=7, xpos = -25,
                             ypos = 0.6, xpos2 = -25,
                             ypos2= 0.57)+
  ylim(-0.1,0.60)+
  xlab("starch change (mg)")+ylab("accumulated growth (cm)")+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  geom_vline(xintercept = 0, linetype="dashed")+
  facet_wrap(~month_changes)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3)),
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),
        legend.position="none",
        plot.title = element_text(size = rel(2.5), face = "italic"),
        strip.text.x = element_text(size= 20, margin = margin(10), vjust = 1),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

#########################################################################################################################################################

### merging the three species together for comparisons between species 


woodstarch_all_species=rbind(Ab_woodStarch_entire_tree,
                             Sac_woodStarch_entire_tree,
                             Dac_woodStarch_entire_tree
                            )
woodstarch_all_species$date=woodstarch_all_species$month

date_vector=c("2019-05-01",
              "2019-08-01",
              "2019-11-01",
              "2020-02-01")
month_vector=unique(woodstarch_all_species$month)

library(plyr)
woodstarch_all_species$date=mapvalues(woodstarch_all_species$date,
                                      month_vector,
                                      date_vector)

woodstarch_all_species$date=as.POSIXct(strptime(as.character(woodstarch_all_species$date),
                                                format = "%Y-%m-%d"))

woodstarch_all_species$species=gsub("[[:digit:]]","", woodstarch_all_species$ID)

species_tag=unique(woodstarch_all_species$species)
species_name=c("O.leucoxylon",
               "S. guianensis",
               "D. microcarpa"
               )

woodstarch_all_species$species=mapvalues(woodstarch_all_species$species,
                                         species_tag,
                                         species_name)

woodstarch_all_species$species_trait=woodstarch_all_species$species
species_tag=unique(woodstarch_all_species$species)
species_name2=c("O. leucoxylon\n(evergreen/parenchyma-storing-species)",
               "S. guianensis\n(semi-deciduous/parenchyma-storing-species)",
               "D. microcarpa\n(semi-deciduous/fiber-storing-species)"
)

woodstarch_all_species$species_trait=mapvalues(woodstarch_all_species$species_trait,
                                         species_tag,
                                         species_name2)


detach("package:plyr", unload=TRUE)

woodstarch_all_species$storage_strategy=ifelse(woodstarch_all_species$species%in%c("D. microcarpa", "T. burserifolia", "T. glaziovii", "T. guianensis"), "Fiber storage", "Parenchyma storage")
woodstarch_all_species$storage_lipids=ifelse(woodstarch_all_species$species%in%c("D. microcarpa", "T. burserifolia", "T. glaziovii", "V. vismiifolia"), "Lipid storage", "No lipid storage")
woodstarch_all_species$leaf_habit=ifelse(woodstarch_all_species$species%in%c("D. microcarpa", "V. vismiifolia", "S. guianensis", "T. guianensis"), "semi-deciduous", "evergreen")


load("data/climatic_Data_2018_2020.RData")


## calculating mean values per species 
woodstarch_all_species_mean_species=woodstarch_all_species%>%group_by(month, species)%>%
  summarise(mean_starch_percentage2=mean(mean_starch_percentage),
            sd_mean_starch=sd(mean_starch_percentage),
            sterr_mean_starch =var(mean_starch_percentage)/sqrt(length(mean_starch_percentage)),
            median_starch_percentage=mean(median_starch_percentage),
            mean_IQR=mean(mean_IQR),
            starch_woodcore_content_by_mean2=mean(starch_woodcore_content_by_mean),
            sd_starch_woodcore_content_by_mean=sd(starch_woodcore_content_by_mean),
            starch_woodcore_content_by_median=mean(starch_woodcore_content_by_median),
            mean_starch_content_mg=mean(starch_woodcore_content_mg),
            sd_starch_content_mg=sd(starch_woodcore_content_mg),
            date=first(.data$date),
            storage_strategy=first(.data$storage_strategy,),
            storage_lipids=first(.data$storage_lipids),
           
  )

### adding significance differences from the comparisons withing the species done above for each species confidence level is 90%
woodstarch_all_species_mean_species$cdl=c("a", "a","a",
                                          "b", "a", "a",
                                          "a", "a", "ab",
                                          "ab", "a", "b")
 
str(woodstarch_all_species_mean_species)
names(woodstarch_all_species)


sp_trait_key=unique(woodstarch_all_species$species_trait)
woodstarch_all_species_mean_species$species_trait=woodstarch_all_species_mean_species$species
species_key=c("O.leucoxylon", "S. guianensis", "D. microcarpa")

library(plyr)
woodstarch_all_species_mean_species$species_trait=mapvalues(woodstarch_all_species_mean_species$species_trait,
                                               species_key,
                                               sp_trait_key)  
detach("package:plyr", unload=TRUE)

png("outputs/total_starch_seasonality_D_S_O_separated_mg_boxplot.png",width = 2000, height = 800)
ggplot(woodstarch_all_species)+
  geom_area(data=climatic_data_months_2018_2020%>%filter(Year=="2019"| Year=="2020"),
            aes(x=date, y=precip_month/3), fill="grey80")+
  geom_boxplot(aes(x=date, y=starch_woodcore_content_mg, fill=storage_strategy, factor=as.factor(date)), )+
  geom_point(aes(x=date, y=starch_woodcore_content_mg))+
  geom_text(data=woodstarch_all_species_mean_species, aes(label=cdl, x=date, y= mean_starch_content_mg+(sd_starch_content_mg+20)), size=8)+
  #geom_errorbar(aes(x=date, ymin=mean_starch_content_mg-sd_starch_content_mg,
  #                 ymax=mean_starch_content_mg+sd_starch_content_mg, factor=species),
  #             position=position_dodge(), width=3000000)+
  scale_x_time(name="",
               breaks = c(unique(c(climatic_data_months_2018_2020%>%filter(Year=="2019"| Year=="2020"))$date))[c(3,4,9,10,15,16, 21, 22)],
               labels = c("Feb19", "Feb20", "May19", "May20", "Aug19", "Aug20", "Nov19", "Nov20"),
               limits = as.POSIXct(strptime(c("2019-01-01", "2020-06-01"), format = "%Y-%m-%d")))+
  #ylab("Starch content (g)")+
  scale_y_continuous(name="Starch mass (mg)",
                     sec.axis=sec_axis(~.*3, name="Rainfall"), limits = c(0,150))+
  #scale_fill_manual("Species")+ #values=c("Brown","darkolivegreen3", "aquamarine4"))+
  #ylim(0,250)+
  facet_wrap(~species_trait, scales = "free")+
  theme_classic()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2.4), face = "italic"),
        legend.position = "none",
        legend.title = element_text(size=rel(3)),
        legend.text = element_text(size = rel(3))
  )
dev.off()

#comparissons between storage strategies grouped at different leaf habit levels 
#checking normality assumption 
ggqqplot(woodstarch_all_species, "starch_woodcore_content_mg", facet.by = c("leaf_habit","storage_strategy"))
lm_allsp=aov(starch_woodcore_content_mg~storage_strategy*leaf_habit, data=
               woodstarch_all_species)
summary(lm_allsp)
TukeyHSD(lm_allsp)

ggplot(woodstarch_all_species)+geom_boxplot(aes(x=storage_strategy, y=starch_woodcore_content_mg))+
  facet_wrap(~leaf_habit)


##### comparison between the changes in starch content between all the months examined. 
woodstarch_all_species2=woodstarch_all_species
woodstarch_all_species2=woodstarch_all_species2%>%filter(ID!="ab26820" & ID!="Sac1054" & ID!="Dac32001")

woodstarch_all_species2$starch_change=NA
woodstarch_all_species2$rel_starch_change=NA
woodstarch_all_species2$changed_months=as.character("NA")

woodstarch_all_species2=woodstarch_all_species2%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 2:nrow(x)){
      x$starch_change[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1]
      x$rel_starch_change[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1])/x$starch_woodcore_content_mg[i]
      x$changed_months[i]=paste(x$month[i-1],"-",x$month[i], sep = "")
       }
    return(x)
  }
  )%>%ungroup()


woodstarch_all_species2$starch_change2=NA
woodstarch_all_species2$rel_starch_change2=NA
woodstarch_all_species2$changed_months2=as.character("NA")

woodstarch_all_species2=woodstarch_all_species2%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 3:nrow(x)){
      x$starch_change2[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-2]
      x$rel_starch_change2[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-2])/x$starch_woodcore_content_mg[i]
      x$changed_months2[i]=paste(x$month[i-2],"-", x$month[i], sep = "") 
       }
    return(x)
  }
  )%>%ungroup()


woodstarch_all_species2$starch_change3=NA
woodstarch_all_species2$rel_starch_change3=NA
woodstarch_all_species2$changed_months3=as.character("NA")

woodstarch_all_species2=woodstarch_all_species2%>%group_by(ID)%>%
  mapGroups(function(x){
      x$starch_change3[4]=x$starch_woodcore_content_mg[4]-x$starch_woodcore_content_mg[4-3]
      x$rel_starch_change3[4]=(x$starch_woodcore_content_mg[4]-x$starch_woodcore_content_mg[4-3])/x$starch_woodcore_content_mg[4]
      x$changed_months3[4]=paste(x$month[4-3],"-", x$month[4], sep = "")     
      return(x)
  }
  )%>%ungroup()

#### melting the columns together

names(woodstarch_all_species2)
id_vars=c("ID", "month", "starch_woodcore_content_mg", "Scientific_name", "species_trait", "storage_strategy", "storage_lipids", "leaf_habit")

woodstarch_all_species2_melt=woodstarch_all_species2%>%melt(id.vars = id_vars, measure.vars =c("starch_change","starch_change2", "starch_change3"),
                                                              value.name = "starch_change" 
                                                              )
  woodstarch_all_species2_melt2=woodstarch_all_species2%>%melt(id.vars = id_vars, measure.vars =c("rel_starch_change","rel_starch_change2", "rel_starch_change3"), 
                                                                    value.name = "rel_starch_change")
  woodstarch_all_species2_melt3=woodstarch_all_species2%>%melt(id.vars = id_vars, measure.vars =c("changed_months","changed_months2", "changed_months3"), 
                                                                    value.name = "changed_months")

  woodstarch_all_species2_melt=cbind(woodstarch_all_species2_melt,rel_starch_change=woodstarch_all_species2_melt2$rel_starch_change, changed_months=woodstarch_all_species2_melt3$changed_months)
  
  woodstarch_all_species2_melt=woodstarch_all_species2_melt[!is.na(woodstarch_all_species2_melt$starch_change),]
  
  unique(woodstarch_all_species2_melt$changed_months)
  woodstarch_all_species2_melt$changed_months=factor( woodstarch_all_species2_melt$changed_months, 
                                                      levels = c("May19-Aug19",
                                                                 "Aug19-Nov19",
                                                                 "Nov19-Feb20",
                                                                 "May19-Nov19",
                                                                 "Aug19-Feb20",
                                                                 "May19-Feb20"))
  
  #### comparison using paired samples (wilcoxon signed rank test) to compare total starch change in the entire wood core between all considered time periods and species######
  str(woodstarch_all_species2_melt)
  
  table(woodstarch_all_species2_melt$ID)
  
  woodstarch_all_species2_melt$wilcox_factor=paste(woodstarch_all_species2_melt$Scientific_name, woodstarch_all_species2_melt$changed_months, sep= "_")
  woodstarch_all_species2_melt$wilcox_factor=as.factor(woodstarch_all_species2_melt$wilcox_factor)

  wilcox_1=wilcox_test(rel_starch_change~wilcox_factor, data=woodstarch_all_species2_melt)

  ### calculating the non-parametric confidence interva for the mean of the total starch mass in the entire wood core at each time period considered.
  ### this is for visually assessing differences from zero and differences between groups. 
  library(boot) 
  
  # function to obtain the mean
  Bmean <- function(data, indices) {
    d <- data[indices] # allows boot to select sample 
    return(mean(d))
  } 
  
  woodstarch_all_species2_melt_boot_ci=woodstarch_all_species2_melt%>%filter(rel_starch_change>-2)%>%group_by(changed_months, Scientific_name)%>%nest()%>%
    mutate(boot_results=map(data, ~ boot(data=.x$rel_starch_change, statistic=Bmean, R=1000)),
           boostrap_var=map(boot_results, ~var(.x$t)),
           vector_bootstrap=unlist(boostrap_var),
           ci_boot=map(boot_results, boot.ci, var.t0=vector_bootstrap))

   names(woodstarch_all_species2_melt_boot_ci$ci_boot)=c(paste(woodstarch_all_species2_melt_boot_ci$Scientific_name,woodstarch_all_species2_melt_boot_ci$changed_months))
  
  mean_groupe=woodstarch_all_species2_melt%>%group_by(changed_months, Scientific_name)%>%
    summarise(mean_groups=mean(rel_starch_change))
  
  boot_results_unnest=woodstarch_all_species2_melt_boot_ci$boot_results
  boot_results_unnest_df=lapply(seq_along(boot_results_unnest),
                                function(i){
                                  unlist(boot_results_unnest[[i]][1])
                                })
  boot_results_unnest_df=as.data.frame(do.call(rbind, boot_results_unnest_df))

  ci_boot_unnest=woodstarch_all_species2_melt_boot_ci$ci_boot
  ci_boot_unnest_df= lapply(seq_along(ci_boot_unnest), function(i){
    unlist(ci_boot_unnest[[i]][7])
  })
  ci_boot_unnest_df=do.call("rbind", ci_boot_unnest_df)
  ci_boot_unnest_df=as.data.frame(ci_boot_unnest_df)
  ci_boot_unnest_df=cbind(ci_boot_unnest_df, sp=woodstarch_all_species2_melt_boot_ci$Scientific_name, month_changes=woodstarch_all_species2_melt_boot_ci$changed_months, 
                          mean=boot_results_unnest_df)

  sp_names=c("D. microcarpa\n(semi-deciduous/fiber-storing-species)",
             "S. guianensis\n(semi-deciduous/parenchyma-storing-species)",
             "O. leucoxylon\n(evergreen/parenchyma-storing-species)"
            )
  names(sp_names)=unique(ci_boot_unnest_df$sp)
  
  ci_boot_unnest_df$storage_strategy=ifelse(ci_boot_unnest_df$sp=="Dacryodes microcarpa", "Fiber_storage", "parenchyma_storage")
  
png("outputs/CI_rel_starch_changes_between_months_D_S_O_separated_boxplot.png",width = 2600, height = 800)
ggplot(ci_boot_unnest_df, aes(x=month_changes, y=t0, col=storage_strategy))+
  geom_point(size=6)+
  geom_errorbar(aes(ymin=bca4, ymax=bca5), size=3)+
  #ylim(-3, 2)+
  xlab("")+ylab("Mean relative starch change (%)")+
  geom_hline(yintercept = 0, linetype="dashed", color="red")+
  facet_wrap(~sp, labeller = labeller(sp=sp_names))+
  theme_classic()+
  theme(axis.text.x=element_text(size = rel(2.4)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3), face = "italic"),
        legend.position = "none",
        legend.title = element_text(size=rel(3)),
        legend.text = element_text(size = rel(3))
  )
dev.off()
  
 
png("outputs/rel_starch_changes_between_months_D_S_O_separated_boxplot.png",width = 1500, height = 2000)
ggplot(woodstarch_all_species2_melt)+
  geom_boxplot(aes(x=changed_months, y=rel_starch_change, fill=storage_strategy, alpha=0.3))+
  geom_hline(yintercept = 0, linetype="dashed", color="red")+
  xlab("")+ylab("Relative starch change (%)")+
  ylim(-4,1)+
  facet_wrap(~species_trait, scales = "free", ncol=1)+
  theme_classic()+
  theme(axis.text.x=element_text(size = rel(3.5)),
        axis.text.y=element_text(size = rel(6)),
        axis.title.x = element_text(size = rel(6)),
        axis.title.y = element_text(size = rel(6)),
        strip.text = element_text(size = rel(5), face = "italic"),
        legend.position = "none",
        legend.title = element_text(size=rel(3)),
        legend.text = element_text(size = rel(3))
  )
dev.off()
