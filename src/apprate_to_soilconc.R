library(ggplot2)

#let's test and see how close we get to real soil concentrations using the weir method
#we'll assume an even mixing depth of 1cm, and also test out some other depths 
lab <- read.csv("data_in/lab_bb.csv")
head(lab)
names(lab)

soil<-na.omit(lab[,c(1,2,9)])
unique(lab$chemical)

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
soil<-soil[,c(3:4,6:14)]
soil_t<-tidyr::gather(soil, "depth","conc",1,3:11)
max(soil$soil_conc_ugg)

#scatter plot
p <- ggplot(soil_t, aes(app_rate_ug_cm2, conc,color=depth))+
  geom_point() +
  ylim(0,200)+
  xlim(0,100)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  ylab("Calculated Soil Sonc (ug/g)") +
  xlab("Application rate (ug/cm2)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key=element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12),
        axis.title.x = element_text(face='bold'))
p

#keep 1cm as mixing depth
