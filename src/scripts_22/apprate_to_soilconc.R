library(ggplot2)
library(cowplot)

#let's test and see how close we get to real soil concentrations using the weir method
#we'll assume an even mixing depth of 1cm, and also test out some other depths 
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')
lab <- read.csv("data_in/lab_bb.csv")
head(lab)
names(lab)

soil<-na.omit(lab[,c(1,2,4,8)])
#unique(soil$chemical)


###scatter plot lm all pesticides----
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
  ylab("Calculated Soil Conc (ug/g)") +
  xlab("Application rate (ug/cm2)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key=element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12),
        axis.title.x = element_text(face='bold'))
p

#keep 1cm as mixing depth


###scatter plot lm individual pesticides----
unique(soil$chemical)
prop <- dplyr::filter(soil, chemical == 'trifloxystrobin')

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

p <- ggplot(prop_t, aes(app_rate_ug_cm2, conc, color=depth, shape= depth))+
  geom_point() +
  ylim(0,0.45)+
  xlim()+
  scale_shape_manual(values=c(1:15))+
  scale_color_manual(values=colorBlindGrey8)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  ylab("Calculated Soil Conc (ug/g)") +
  xlab("Application rate (ug/cm2)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key=element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12),
        axis.title.x = element_text(face='bold'))
p


avg_conc<-prop_t %>% group_by(depth) %>% summarise(avg = mean(conc))
#0.302 application rate results in a soil concentration of 0.260 for trifloxystrobin
# the closest estimated mixing depth was 7mm, with a concentration of 0.270


###boxplot with error between estimates----
#do 1cm, 1mm, 5mm, 1 inch for all chemicals
soil$app_rate_ug_cm2<-soil$app_rate_g_cm2*1000000
#convert from g/cm2 to mg/m2
soil$conc_mgm2<-(soil$app_rate_ug_cm2*0.001*10000)

#1 cm
soil$conc_ugg_1cm<-(soil$conc_mgm2/16000)*1000
#1mm
soil$conc_ugg_1mm<-(soil$conc_mgm2/1600)*1000
#1 inch
soil$conc_ugg_1inch<-(soil$conc_mgm2/40640)*1000
#5mm
soil$conc_ugg_5mm<-(soil$conc_mgm2/8000)*1000

names(soil)
soil<-soil[,c(3:10)]
soil_edit<-tidyr::gather(soil, "depth","conc",5:8)

unique(soil_edit$depth)

soil_edit$Ratio<-soil_edit$conc/soil_edit$soil_conc_ugg #original/calculated
names(soil_edit)
unique(soil_edit$depth)

CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

soil_edit$Chemical<-sapply(soil$chemical, CapStr)


t_1cm <- dplyr::filter(soil_edit, depth == 'conc_ugg_1cm')
#t_1cm$look<-log2(t_1cm$Ratio)
p1 <- ggplot(t_1cm, aes(x=Chemical, y=log2(Ratio))) + 
  geom_boxplot(aes(fill=Chemical))+
  # scale_y_continuous(breaks=seq(0,100, 25))+
  ggtitle("1 cm Depth") +
  ylab("Log2(Ratio)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16),
        legend.position = 'none')
p1



t_5mm <- dplyr::filter(soil_edit, depth == 'conc_ugg_5mm')
p2 <- ggplot(t_5mm, aes(x=Chemical, y=log2(Ratio))) + 
  geom_boxplot(aes(fill=Chemical))+
  scale_y_continuous(breaks=seq(-100,175, 50))+
  ggtitle("5 mm Depth") +
  ylab("Log2(Ratio)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16),
        legend.position = 'none')
p2


t_1mm <- dplyr::filter(soil_edit, depth == 'conc_ugg_1mm')
p3 <- ggplot(t_1mm, aes(x=Chemical, y=log2(Ratio))) + 
  geom_boxplot(aes(fill=Chemical))+
  scale_y_continuous(breaks=seq(-500,50,100))+
  ggtitle("1 mm Depth") +
  ylab("Log2(Ratio)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16),
        legend.position = 'none')
p3


t_1inch <- dplyr::filter(soil_edit, depth == 'conc_ugg_1inch')
p4 <- ggplot(t_1inch, aes(x=Chemical, y=log2(Ratio))) + 
  geom_boxplot(aes(fill=Chemical))+
  # scale_y_continuous(breaks=seq(-500,50,100))+
  ggtitle("1 inch Depth") +
  ylab("Log2(Ratio)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16),
        legend.position = 'none')
p4


plot_grid(p1,p4,p2,p3)

#error < 0 mean that calculated soil concentrations at the mixing depth were higher than the actual soil concentration
#indicating that that mixing depth over-predicts the concentration

#error > 0 mean that calculated soil concentrations at the mixing depth were higher than the actual soil concentration
#indicating that that mixing depth under-predicts the concentration

