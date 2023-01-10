library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(cowplot)
library(ggrepel)

#12/27/22

##Figure S1 is in 'apprate_to_soilconc.R'
##Figure S2
lab_bb<-read.csv('data_in/lab_bb.csv')#should be 1158
field_bb<-read.csv('data_in/field_bb.csv')  #should be 30611

#read in final versions of the estimated dermal dose 
dd_pyra<-read.csv('data_out/pyraclostrobin_final_DD.csv')#should be 1158
mod<-read.csv('data_in/As_Modifier.csv')
param<-read.csv('data_in/parameters.csv')

#let's extract the pyraclostrobin data from the field dataset
pyra_f<-field_bb[field_bb$chemical == 'Pyraclostrobin', ]  
pyra_f<-na.omit(pyra_f)
pyra_f$unit<-'ug/g' 

#let's get the trifloxistrobin data 
strob_l<-lab_bb[lab_bb$chemical == 'trifloxystrobin', ]  
strob_l<-na.omit(strob_l)
strob_l$unit<-'ug/g' 
#get apprate in ug/cm2
strob_l$app_rate_ug_cm2<-strob_l$app_rate_g_cm2*1000000



## We'll make a boxplot that will show the differences in lab, field, and estimated tiss. concn with
## an application rate of 0.302 ug/cm2 and using the OG bodyweights
effects<-dd_pyra
effects<-effects[effects$Study %in% "Cusaac_2017c" | effects$Study %in% "Cusaac_2017d",]
#I don't need the modifier here; we're going to only extract soil-only applications
#let's calculate a new dermal dose based on the 0.302 ug/cm2 application rate used in the trifloxistrobin
pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.91*24 #from comptox profile
dt<-param[1,2]
move_rate<-param[3,2]
bioavail<-param[5,2]
dsa<-param[9,2]*effects$M_Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_deg<-log(2)/hl*8/move_rate #with adjusted time of 8h to match lab dataset
soil_concs<-((0.302*10)/16000)*1000 #1cm mixing depth
dermal_dose<-(soil_concs^soil_concs_deg * kp_pyra * (dsa/dt) * derm_frac * bioavail)/effects$M_Body_Weight_g

#we've taken the estimated DDs based on the body weights, an application rate of 0.302, and 8h durations. 
effects$dermaldose<-dermal_dose

dd_est<-cbind(effects$dermaldose,"Estimate")
dd_l<-cbind(strob_l$tissue_conc_ugg,"Lab")
dd_f<-cbind(pyra_f$tissue, "Field")

final_df<-as.data.frame(rbind(dd_est,dd_l,dd_f))
names(final_df)<-c("DD","Type")

final_df$aq<-as.numeric(final_df$DD)* (0.05*(10^3.99))



ratio_all<-ggplot(final_df, aes(x=as.factor(Type), y=(as.numeric(DD)), fill=Type, color=Type)) +
  geom_boxplot()+
  #facet_wrap(.~Type, scales = "free")+
  labs(y = expression(paste('Log (Tissue Concentration) [ ug ', g^-1, ' ]')), 
       x = NULL)+ # adds in labels for x and y axis
  
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(margin = margin(t = 10, r = 0, b = , l = 0), size=14,face="bold"),
        axis.title.y=element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=14,face="bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=1, size=12))+
  theme(legend.position = "none")
ratio_all


mean(strob_l$body_weight_g)

