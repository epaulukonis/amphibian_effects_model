library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(cowplot)

##Script for comparison of field and lab-based tissue concentrations to estimated tissue-based concentrations
lab_bb<-read.csv('data_in/lab_bb.csv')#should be 1158
field_bb<-read.csv('data_in/field_bb.csv')  #should be 30611

#pyraclostrobin has a relatively high soil sorption factor, i.e., the bioavailable fraction of pyra
#rapidly decreases after soil contact -  however, water content in soil prolongs the period of potential soil exposure
#cusaac et al. 2017 performed a second experiment where they exposed toads 0, 60, and 120 minutes after 
#spraying, and also included with/without simulated
#they found that mortality was eliminated when not exposed immediately; 2.5and 5.0 ug/cm2 application rates

#we don't obviously have application rates or body weights for these field data
#'3-5g' in weight

#we could use this as an example of field-level tissue concentrations
#note that further experiments are needed to better capture low application rates or more appropriately, 
#survival rates and dermal doses after real-life field conditions?



#let's look at see what chemicals we've got to compare 
sort(unique(field_bb$chemical)) 
sort(unique(lab_bb$chemical))



#field
#get pyraclsotrobin data from field data
pyra_f<-field_bb[field_bb$chemical == 'Pyraclostrobin', ]  
pyra_f<-na.omit(pyra_f)
pyra_f$unit<-'ug/g' 

ggplot(pyra_f, aes(y = Species)) +
  geom_density_ridges(
    aes(x=tissue,fill=paste(Species, Source)),
    scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Species") +
  xlab("Tissue Concentration (ug/g)")+
  theme_ridges()



#lab
#get trifloxystrobin pesticide from lab data
strob_l<-lab_bb[lab_bb$chemical == 'trifloxystrobin', ]  
strob_l<-na.omit(strob_l)
strob_l$unit<-'ug/g' 

ggplot(strob_l, aes(y = species)) +
  geom_density_ridges(
    aes(x=tissue_conc_ugg,fill=paste(species, source)),
    scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Species") +
  xlab("Tissue Concentration (ug/g)")+
  theme_ridges()


names(effectsp)
pyra_est<-effectsp[,c(1,5:6,8,11,21)]

names(strob_l)
names(pyra_est)

unique(pyra_est$Application_Rate)
unique(strob_l$app_rate_g_cm2)  #need to convert to ug/cm2
strob_l$Application_Rate<-strob_l$app_rate_g_cm2*1000000

mean(strob_l$tissue_conc_ugg)
pyra_est_lowdose<-pyra_est %>% filter(Application_Rate < 1)

ggplot(pyra_est_lowdose, aes(x = dermaldose, y = Application_Rate))+
  geom_point(colour=Application_Rate)

#compare dermal doses calculated at similar app rates and the 8 hour trifloxystrobin next week


#### toxicity data
# pyraclostrobin
effectsp<-read.csv('data_in/Headline_updated.csv')
param<-read.csv('data_in/parameters.csv')
effectsp$half_life<-(effectsp$Duration_h*log(2))/log(1/effectsp$Survival) #need to modify survival by duration, using an exponential growth curve
effectsp$adj_sur_96<-round(1/(2^(96/effectsp$half_life)),3)

pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.91*24 #from comptox profile
dt<-param[1,2]
move_rate<-param[3,2]
bioavail<-param[5,2]
dsa<-param[9,2]*effectsp$M_Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_deg<-log(2)/hl*96/move_rate #with adjusted time of 96h to match adjusted survival
soil_concs<-((effectsp$Application_Rate*10)/16000)*1000 #1cm mixing depth; or change depending on soil?

#this calculates dermal dose 
dermal_dose<-(soil_concs^soil_concs_deg * kp_pyra * (dsa/dt) * derm_frac * bioavail)/effectsp$M_Body_Weight_g
head(dermal_dose) #take a look
effectsp$dermaldose<-dermal_dose
effectsp$Mortality<-1-effectsp$adj_sur_96


# glyphosate
effectsg<-read.csv('data_in/Glyphosate_updated.csv')
param<-read.csv('data_in/parameters.csv')#need to modify survival by duration, using an exponential growth curve
effectsg$adj_sur_96<-effectsg$Survival

pmolweight<-169.1 #g/mol #comptox
plogKow<- -3.12 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.46*24 #from comptox profile
dt<-param[1,2]
move_rate<-param[3,2]
bioavail<-param[5,2]
dsa<-param[9,2]*effectsg$Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_deg<-log(2)/hl*96/move_rate #with adjusted time of 96h to match adjusted survival
soil_concs<-((effectsg$Application_Rate*10)/16000)*1000 #1cm mixing depth; or change depending on soil?

#this calculates dermal dose 
dermal_dose<-(soil_concs^soil_concs_deg * kp_pyra * (dsa/dt) * derm_frac * bioavail)/effectsg$Body_Weight_g
head(dermal_dose) #take a look
effectsg$dermaldose<-dermal_dose
effectsg$Mortality<-1-effectsg$adj_sur_96



field<-ggplot(pyra_f, aes(y = Species)) +
  geom_density_ridges(aes(x=tissue,fill=paste(Species)), stat = "binline", bins = 50, scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.01,.25)) +
  coord_cartesian(clip = "off") +
  ylab("Species") +
  xlab("Tissue Concentration (ug/g)")+
  guides(fill=guide_legend(title=""))+
  theme_ridges()
field


lab<-ggplot(lab_bb, aes(y = species)) +
  geom_density_ridges(aes(x=tissue_conc_ugg,fill=paste(species)),stat = "binline", bins = 50,scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Species") +
  xlab("Tissue Concentration (ug/g)")+
  guides(fill=guide_legend(title=""))+
  # theme(axis.title.y = element_text( size=12, face="bold"), 
  #       axis.title.x = element_text( size=12, face="bold"),
  #       legend.title = element_text(size=14, face="bold"),
  #       legend.text = element_text( size=12, face="bold"),
  #       )+
  theme_ridges()
lab

est<-ggplot(effectsg, aes(y = Species)) +
  geom_density_ridges(aes(x=dermaldose,fill=paste(Species)),stat = "binline", bins = 50,scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Species") +
  xlab("Tissue Concentration (ug/g)")+
  guides(fill=guide_legend(title=""))+
  # theme(axis.title.y = element_text( size=12, face="bold"), 
  #       axis.title.x = element_text( size=12, face="bold"),
  #       legend.title = element_text(size=14, face="bold"),
  #       legend.text = element_text( size=12, face="bold"),
  #       )+
  theme_ridges()
est



plot_grid(field,calc)
#plot_grid(field,calc, labels = c('Field Body Burdens', 'Calculated Body Burdens'), label_size = 12)
