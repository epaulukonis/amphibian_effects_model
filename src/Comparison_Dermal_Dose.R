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
#note that further experiments need to better capture low application rates or more appropriately, 
#survival rates and dermal doses after real-life field conditions?

#we could also consider additional soil mixing rates?

#let's look at see what chemicals we've got to compare 
sort(unique(field_bb$chemical)) 
sort(unique(lab_bb$chemical))

#field
sort(unique(field_bb$chemical)) ##glyphosate not on here
pyra_f<-field_bb[field_bb$chemical == 'Pyraclostrobin', ]  
pyra_f<-na.omit(pyra_f)
pyra_f$unit<-'ug/g' ##these are old and I forgot to update them!

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

#toxicity data
effects<-read.csv('data_in/Headline_updated.csv')
param<-read.csv('data_in/parameters.csv')

pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
dt<-param[1,2]
move_rate<-param[3,2]
hl<-4.91*24 #from comptox profile
bioavail<-param[5,2]
dsa<-param[9,2]*effects$M_Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_deg<-log(2)/hl*96/move_rate #with adjusted time of 96h to match adjusted survival
soil_concs<-((effects$Application_Rate*10)/16000)*1000 #1cm mixing depth; or change depending on soil?

#this calculates dermal dose 
dermal_dose<-(soil_concs^soil_concs_deg * kp_pyra * (dsa/dt) * derm_frac * bioavail)/effects$M_Body_Weight_g
head(dermal_dose) #take a look
effects$dermaldose<-dermal_dose

soils<-effects[effects$Study == 'Cusaac_2017c'|effects$Study == 'Cusaac_2017d', ] 
spray<-effects[c(1:30,39:47),]
#spray<-effects[effects$Study != 'Cusaac_2017c'|effects$Study != 'Cusaac_2017d', ] 

median(soils$dermaldose)
median(spray$dermaldose)
mean(pyra_f$tissue)
min(effects$dermaldose)


unique(soils$Species)
unique(spray$Species)
unique(pyra_f$Species)

unique(pyra$Source)



field<-ggplot(pyra_f, aes(y = Species)) +
  geom_density_ridges(aes(x=tissue,fill=paste(Species)), stat = "binline", bins = 100, scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.01,.25)) +
  coord_cartesian(clip = "off") +
  ylab("Species") +
  xlab("Tissue Concentration (ug/g)")+
  guides(fill=guide_legend(title=""))+
  theme_ridges()
field

calc<-ggplot(effects, aes(y = Species)) +
  geom_density_ridges(aes(x=dermaldose,fill=paste(Species)),stat = "binline", bins = 100,scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Species") +
  xlab("Tissue Concentration (ug/g)")+
  guides(fill=guide_legend(title=""))+
  theme_ridges()
calc



plot_grid(field,calc, labels = c('Field Body Burdens', 'Calculated Body Burdens'), label_size = 12)

