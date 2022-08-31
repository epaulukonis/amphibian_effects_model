#create plots for Paulukonis et al. 2021
library(Rcpp)
library(RcppGSL)
library(RcppEigen)
library(forcats)
library(plotly)
library(shiny)
library(scales)
library(ToxicR)
library(nlme)
library(ggplot2)
library(ggridges)
library(gridExtra)
library(grid)
library(readr)
library(bayestestR)
library(dplyr)
library(tidyverse)
library(cowplot)

set.seed(6379)


setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')
#load('myEnvironment_analysis.RData')


##this script contains plot and data analysis for the BMDS MCMC outputs
#all plots here are part of the manuscript
#some of the script will be similar to the comparison plots found in 'Comparison_BBMD_BMDS.R'



############################Plots---- 
#### Figure 1

#Acute mortality by body weight (g) for terrestrial stage anurans for 
#a) pyraclostrobin and b) glyphosate exposures. For pyraclostrobin, soil-only exposures are denoted in red. 


mortp<-ggplot(effectsp, aes(x=M_Body_Weight_g, y=Mortality)) +
  geom_point(size=2) 
mortp


#### Figure 2
mod<-read.csv('data_in/As_Modifier.csv')

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
effectsp <-effectsp[order(effectsp$Study),]


fam<-as.data.frame(cbind(effectsp$Family, effectsp$Study, effectsp$Application_Rate))
names(fam)<-c('Family', 'Study', "Application")
mod_As<- merge(fam,mod, by  = "Family") 
mod_As<-mod_As[order(mod_As$Study),]

#adjust with modifier
bw_effect<-effectsp$M_Body_Weight_g
exp<-bw_effect^mod_As$Exponent 
as<-mod_As$Modifier*exp
derm_d<-effectsp[c(1:39),23] #only do the direct exposures
app_d<-effectsp[c(1:39),9] #only do the direct exposures
as<-as[1:39] #only do the direct exposures
derm_original<-derm_d*((as*as.numeric(app_d))/2) #modify the calculated dermal dose by the SA and application rate product
derm_original<-c(derm_original,effectsp[c(40:47),23]) 
effectsp$dermaldose<-derm_original

#note; this is with no distinction between the zero doses and dropped experiments; 47 records


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

bw_effect<-effectsg$Body_Weight_g
exp<-bw_effect^mod_As$Exponent 
as_d<-mod_As$Modifier*exp
derm_d<-effectsg$dermaldose 
app_d<-as.numeric(effectsg$Application_Rate)
derm_original<-derm_d*((as_d*app_d)/2)

effectsg$dermaldose<-derm_original


#run pyraclostrobin in Data_Simulation.R
effectsp$label<-paste(effectsp$Study, effectsp$Species, sep="")
effectsp <-effectsp[order(effectsp$label),]
unique(effectsp$label) # get unique combos of study and species

# Figure 1: application rate and 96hr mortality figure?
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#762a83", "#000000","#7fbc41")

pp<-ggplot(effectsp, aes(x=dermaldose, y=Mortality, colour= Study_Type,shape=label)) +
  geom_point(size=2) + 
  scale_shape_manual(name = "Study and Species", 
    labels = c("Belden et al. 2010, Anaxyrus cognatus","Bruhl et al. 2013, Rana temporaria","Cusaac et al. 2015, Anaxyrus woodhousii", "Cusaac et al. 2015, Acris blanchardi",
               "Cusaac et al. 2016a, Acris blanchardi", "Cusaac et al. 2016b, Acris blanchardi",
               "Cusaac et al. 2017a, Anaxyrus cognatus", "Cusaac et al. 2017a, Anaxyrus woodhousii", "Cusaac et al. 2017b, Anaxyrus cognatus",
               "Cusaac et al. 2017c, Anaxyrus cognatus", "Cusaac et al. 2017d, Anaxyrus cognatus"),
    values=c(3:4,7:12,15,17,19))+
  # scale_color_manual( name = "Study and Species",
  #                     labels = c("Belden et al. 2010, Anaxyrus cognatus","Bruhl et al. 2013, Rana temporaria","Cusaac et al. 2015, Anaxyrus woodhousii", "Cusaac et al. 2015, Acris blanchardi",
  #                                "Cusaac et al. 2016a, Acris blanchardi", "Cusaac et al. 2016b, Acris blanchardi",
  #                                "Cusaac et al. 2017a, Anaxyrus cognatus", "Cusaac et al. 2017a, Anaxyrus woodhousii", "Cusaac et al. 2017b, Anaxyrus cognatus",
  #                                "Cusaac et al. 2017c, Anaxyrus cognatus", "Cusaac et al. 2017d, Anaxyrus cognatus"),
  #   values=colorBlindGrey8)+
  # # geom_smooth(method=lm)+
  ggtitle("Adjusted 96hr Mortality vs. Estimated Dermal Dose, Pyraclostrobin") +
  ylab("Mortality") +
  xlab("Dose (ug/g)")+
  #scale_x_continuous(limits = c(0, 50))+
  #scale_color_manual(values='#2ca25f')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16))
        #legend.position = 'none') #first run without this to get legend
pp

#next, run glyphosate in Data_Simulation.R
effectsg$label<-paste(effectsg$Study, effectsg$Species, sep="")
effectsg <-effectsg[order(effectsg$label),]
unique(effectsg$label) # get unique combos of study and species

# Figure 1: application rate and 96hr mortality figure?
pg<-ggplot(effectsg, aes(x=dermaldose, y=Mortality, color=label, shape =  label)) +
  geom_point(size=2) + 
  scale_shape_manual(name = "Study and Species", 
                     labels = c("Bernal et al. 2009, Centrolene prosblepon","Bernal et al. 2009, Engystomops pustulosus", 
                                "Bernal et al. 2009, Pristimantis taeniatus","Bernal et al. 2009, Rhinella granulosa",
                                "Bernal et al. 2009, Rhinella marina","Bernal et al. 2009, Rhinella roqueana",
                                "Bernal et al. 2009, Scinax ruber","Meza-Joya et al. 2013, Eleutherodactylus johnstonei"),
    values=c(3:4,7:12))+
  scale_color_manual(name = "Study and Species", 
                     labels = c("Bernal et al. 2009, Centrolene prosblepon","Bernal et al. 2009, Engystomops pustulosus", 
                                "Bernal et al. 2009, Pristimantis taeniatus","Bernal et al. 2009, Rhinella granulosa",
                                "Bernal et al. 2009, Rhinella marina","Bernal et al. 2009, Rhinella roqueana",
                                "Bernal et al. 2009, Scinax ruber","Meza-Joya et al. 2013, Eleutherodactylus johnstonei"),
                     values=colorBlindGrey8)+
  # geom_smooth(method=lm)+
  ggtitle("Adjusted 96hr Mortality vs. Calculated Dermal Dose, Glyphosate") +
  ylab("Mortality") +
  xlab("Tissue Concentration (ug/g)")+
  # scale_color_manual(values='#2ca25f')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16))
        #legend.position = 'none') #first run without this to get legend

pg

grid.arrange(pp,pg, nrow = 2, ncol = 1)



#### Figure 3 done in MCMC_Pyra.R and MCMC_Glypho.R
#### Figure 4 BMDs with curve
ggplot(bmds_est, aes(y = Model)) +
  geom_density_ridges(aes(x=BMDSEstimates, fill = order), 
                      stat = "binline", bins = 100, scale = 0.9,
                      alpha = .6, ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.5,1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Fitted BMD with Upper (BMDU) and Lower (BMDL) Estimates",
    y = "Model",
    title = "Benchmark Dose (BMD) Estimates by Model" ) +
  theme_ridges()

#just for bmds
ggplot(bmds, aes(x=BMDSEstimates, group=Model)) + 
  geom_density(aes(fill=Model), alpha=0.2)+
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#56B4E9", "#009E73"))+
  xlim(-2.5, 1)+
  xlab('BMDS')+
  ylab('Density')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12))



mcmc_ma$exp_n1$
plot(mcmc_ma$exp_n1)
plot(mcmc_q$exp_n1)
plot(mcmc_ma$exp_n1$Individual_Model_8)
derm_seq <-as.data.frame(seq(0,1.4, length=1000))  #sequence from 0 to 1.4 by x, 1000 length
colnames(derm_seq)[1]<-'dermdose'


save.image(file='myEnvironment_analysis.RData') 



#Figure X - plot of posterior probability weights across 9 models 
# ggplot(postw, aes(x = PosteriorProbs, y = Model, group = Model, fill = Model)) +
#   geom_density_ridges(stat = "binline", bins = 100, scale = 1.5) +
#   scale_y_discrete(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   coord_cartesian(clip = "off") +
#   ylab("Model") +
#   xlab("Posterior Probability by Model for 1000 Simulated Studies")+
#   theme_ridges()


