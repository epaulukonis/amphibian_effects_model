library(ggplot2)

#let's test and see how close we get to real soil concentrations using the weir method
#we'll assume an even mixing depth of 1cm, and also test out some other depths 
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')
lab <- read.csv("data_in/lab_bb.csv")
head(lab)
names(lab)

soil<-na.omit(lab[,c(1,2,4,8)])
unique(lab$chemical)

##all pesticides
#let's also alter the app rate units to match our dataset
soil$app_rate_ug_cm2<-soil$app_rate_g_cm2*1000000
#convert from g/cm2 to mg/m2
soil$conc_mgm2<-(soil$app_rate_ug_cm2*0.001*10000)

#1 cm
soil$conc_ugg_1cm<-(soil$conc_mgm2/16000)*1000
#10 cm
soil$conc_ugg_10cm<-(soil$conc_mgm2/160000)*1000
#2.54 cm (1 inch)
soil$conc_ugg_1inch<-(soil$conc_mgm2/40640)*1000
#5.08 cm (2 inch)
soil$conc_ugg_2inch<-(soil$conc_mgm2/81280)*1000
#1mm
soil$conc_ugg_1mm<-(soil$conc_mgm2/1600)*1000
#5mm
soil$conc_ugg_5mm<-(soil$conc_mgm2/8000)*1000
#7mm
soil$conc_ugg_7mm<-(soil$conc_mgm2/11200)*1000
#9mm
soil$conc_ugg_9mm<-(soil$conc_mgm2/14400)*1000
#8mm
soil$conc_ugg_8mm<-(soil$conc_mgm2/12800)*1000

names(soil)
soil_e<-soil[,c(4:15)]
soil_t<-tidyr::gather(soil_e, "depth","conc",1,4:12)
max(soil_e$soil_conc_ugg)

#scatter plot
p <- ggplot(soil_t, aes(app_rate_ug_cm2, conc,color=depth))+
  geom_point() +
  ylim(0,200)+
  xlim(0,100)+
 # geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  ylab("Calculated Soil Sonc (ug/g)") +
  xlab("Application rate (ug/cm2)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key=element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12),
        axis.title.x = element_text(face='bold'))
p

#keep 1cm as mixing depth


##individual pesticides
prop <- dplyr::filter(soil, chemical == 'atrazine')

#let's also alter the app rate units to match our dataset
prop$app_rate_ug_cm2<-prop$app_rate_g_cm2*1000000
#convert from g/cm2 to mg/m2
prop$conc_mgm2<-(prop$app_rate_ug_cm2*0.001*10000)

#1 cm
prop$conc_ugg_1cm<-(prop$conc_mgm2/16000)*1000
#10 cm
prop$conc_ugg_10cm<-(prop$conc_mgm2/160000)*1000
#2.54 cm (1 inch)
prop$conc_ugg_1inch<-(prop$conc_mgm2/40640)*1000
#5.08 cm (2 inch)
prop$conc_ugg_2inch<-(prop$conc_mgm2/81280)*1000
#1mm
prop$conc_ugg_1mm<-(prop$conc_mgm2/1600)*1000
#5mm
prop$conc_ugg_5mm<-(prop$conc_mgm2/8000)*1000
#7mm
prop$conc_ugg_7mm<-(prop$conc_mgm2/11200)*1000
#9mm
prop$conc_ugg_9mm<-(prop$conc_mgm2/14400)*1000
#8mm
prop$conc_ugg_8mm<-(prop$conc_mgm2/12800)*1000

prop_e<-prop[,c(4:15)]
prop_t<-tidyr::gather(prop_e, "depth","conc",1,4:12)
max(prop_e$soil_conc_ugg)

#scatter plot
p <- ggplot(prop_t, aes(app_rate_ug_cm2, conc,color=depth))+
  geom_point() +
  ylim(0,50)+
  xlim(22,25)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  ylab("Calculated Soil Conc (ug/g)") +
  xlab("Application rate (ug/cm2)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key=element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12),
        axis.title.x = element_text(face='bold'))
p


##boxplot with error between estimates
#do 1cm, 1mm, 1 inch for all chemicals

soil$app_rate_ug_cm2<-soil$app_rate_g_cm2*1000000
#convert from g/cm2 to mg/m2
prop$conc_mgm2<-(prop$app_rate_ug_cm2*0.001*10000)

#1 cm
soil$conc_ugg_1cm<-(soil$conc_mgm2/16000)*1000
#1mm
soil$conc_ugg_1mm<-(soil$conc_mgm2/1600)*1000
#1 inch
soil$conc_ugg_1inch<-(soil$conc_mgm2/40640)*1000

names(soil)
soil<-soil[,c(3:9)]
soil_edit<-tidyr::gather(soil, "depth","conc",5:7)

soil_edit$error<-soil_edit$soil_conc_ugg-soil_edit$conc
names(soil_edit)
unique(soil_edit$depth)

t_1cm <- dplyr::filter(soil_edit, depth == 'conc_ugg_1cm')



p <- ggplot(t_1cm, aes(x=chemical, y=error)) + 
  geom_boxplot()
p






