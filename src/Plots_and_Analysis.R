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



setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')
#load('myEnvironment_analysis.RData')


##this script contains plot and data analysis for the BMDS MCMC outputs
#all plots here are part of the manuscript
#some of the script will be similar to the comparison plots found in 'Comparison_BBMD_BMDS.R'



#### Pyraclostrobin Analysis ----
sims<-read.csv('data_out/BMDS_headline_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T) #split by simulation
mod_list<-sort(c('hill','log-logistic','logistic','log-probit','weibull','qlinear','probit','multistage','gamma')) #order the models by name
mcmc_ma = lapply(by_s, function(y) ma_dichotomous_fit(y[,2],y[,4],y[,3], fit_type = "mcmc")) #apply ma function over list


##first nested level of function output
mcmc_ma$exp_n1$ma_bmd
#ma_bmd -  cdf of the BMD? first used to calculate the CL of the bmd 
mcmc_ma$exp_n1$bmd
#median bmd, bmdl, bmdu
mcmc_ma$exp_n1$posterior_probs
#prob weights for all models

##second nested level of function output
mcmc_ma$exp_n1$Individual_Model_1$mcmc_result #30,000 rows of param samples 
post<-lapply(mcmc_ma, function (x) x['posterior_probs']) #pull out posterior probabilities 
postw<-as.data.frame(unlist(post))
postw<-tibble::rownames_to_column(postw, "Simulation")
postw$Simulation = substr(postw$Simulation,1,nchar(postw$Simulation)-1) #remove extraneous characters
colnames(postw)[2]<-'PosteriorProbs'
mod_list_f<-rep(mod_list, times=1000)
postw$Model<-mod_list_f #add column for model
#write.csv(postw, 'data_out/PosteriorWeights.csv')

#look at median of PW across all models 
tab_med<-postw %>%
  group_by(Model) %>% 
  summarize(med=median(PosteriorProbs))
print(tab_med)
write.csv(tab_med,'data_out/BMDS_posteriorweight_median_headline.csv')

#pull model average BMD, bmdl, bmdu
bmdsorder<-c('bmds','bmdl','bmdu')
bmds_order<-rep(bmdsorder, times=1000)
bmds<-lapply(mcmc_ma, function (x) x['bmd']) #pull out bmds and bmdls
bmds<-as.data.frame(unlist(bmds))
bmds<-tibble::rownames_to_column(bmds, "Simulation")
bmds$Simulation = substr(bmds$Simulation,1,nchar(bmds$Simulation)-1)
colnames(bmds)[2]<-'BMDSEstimates'
bmds$order<-bmds_order
bmds$Model<-'ModelAverage'

#get bmds average
tab_bmds<-bmds%>%
  group_by(order) %>% 
  summarize(med=median(BMDSEstimates), SD=sd(BMDSEstimates))
tab_bmds


#run mcmc over top 3 models
mcmc_ll<- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "probit",fit_type = "mcmc"))
mcmc_p <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-probit",fit_type = "mcmc"))
mcmc_q <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-logistic",fit_type = "mcmc"))

#try plotting with built-in plot function
# out_mcmc_result<-mcmc_ll$exp_n1$mcmc_result
# out_mcmc_fittedmodel<-mcmc_ll$exp_n1$fitted_model
# out_mcmc_prior <-mcmc_ll$exp_n1$prior
# out_bmd_samples<-as.data.frame(mcmc_ll$exp_n1$mcmc_result$BMD_samples)
# plot(mcmc_ll$exp_n1)
# df<-mcmc_ll$exp_n1
# plot(test$V1~test$V2) #doesn't seem to be original plot; goes above 0.5 (plotly plot doesn't)
plot(mcmc_ma$exp_n500)
plot(mcmc_ma$exp_n1$Individual_Model_1)

ggplot_fit<-plot(mcmc_ma$exp_n500)
ggplot_fit + scale_x_continuous(trans = "pseudo_log") + theme_classic()

MAdensity_plot(mcmc_ma$exp_n1)+ggtitle("Density Plot MCMC")

#organize and put outputs in dataframe
ll<-lapply(mcmc_ll, function (x) x['bmd'])
ll<-as.data.frame(unlist(ll))
ll<-tibble::rownames_to_column(ll, "Simulation")
ll$Simulation = substr(ll$Simulation,1,nchar(ll$Simulation)-1)
colnames(ll)[2]<-'BMDSEstimates'
ll$order<-bmds_order
ll$Model<-'Probit'

p<-lapply(mcmc_p, function (x) x['bmd'])
p<-as.data.frame(unlist(p))
p<-tibble::rownames_to_column(p, "Simulation")
p$Simulation = substr(p$Simulation,1,nchar(p$Simulation)-1)
colnames(p)[2]<-'BMDSEstimates'
p$order<-bmds_order
p$Model<-'Log-Probit'

q<-lapply(mcmc_q, function (x) x['bmd'])
q<-as.data.frame(unlist(q))
q<-tibble::rownames_to_column(q, "Simulation")
q$Simulation = substr(q$Simulation,1,nchar(q$Simulation)-1)
colnames(q)[2]<-'BMDSEstimates'
q$order<-bmds_order
q$Model<-'Log-Logistic'

bmds_est<-rbind(bmds,ll,p,q)
write.csv(bmds_est,'data_out/BMDS_headline_top3.csv')


##### Glyphosate Analysis ----
sims<-read.csv('data_out/BMDS_glyphosate_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T) #split by simulation
mod_list<-sort(c('hill','log-logistic','logistic','log-probit','weibull','qlinear','probit','multistage','gamma')) #order the models by name
mcmc_ma_g = lapply(by_s, function(y) ma_dichotomous_fit(y[,2],y[,4],y[,3], fit_type = "mcmc")) #apply ma function over list




##first nested level of function output
mcmc_ma$exp_n1$ma_bmd
#ma_bmd -  cdf of the BMD? first used to calculate the CL of the bmd 
mcmc_ma$exp_n1$bmd
#median bmd, bmdl, bmdu
mcmc_ma$exp_n1$posterior_probs
#prob weights for all models

##second nested level of function output
mcmc_ma$exp_n1$Individual_Model_1$mcmc_result #30,000 rows of param samples 


post<-lapply(mcmc_ma, function (x) x['posterior_probs']) #pull out posterior probabilities 
postw<-as.data.frame(unlist(post))
postw<-tibble::rownames_to_column(postw, "Simulation")
postw$Simulation = substr(postw$Simulation,1,nchar(postw$Simulation)-1) #remove extraneous characters
colnames(postw)[2]<-'PosteriorProbs'
mod_list_f<-rep(mod_list, times=1000)
postw$Model<-mod_list_f #add column for model
#write.csv(postw, 'data_out/PosteriorWeights.csv')

#look at median of PW across all models 
tab_med<-postw %>%
  group_by(Model) %>% 
  summarize(med=median(PosteriorProbs))
print(tab_med)
write.csv(tab_med,'data_out/BMDS_posteriorweight_median_glyphosate.csv')

#pull model average BMD, bmdl, bmdu
bmdsorder<-c('bmds','bmdl','bmdu')
bmds_order<-rep(bmdsorder, times=1000)
bmds<-lapply(mcmc_ma_g, function (x) x['bmd']) #pull out bmds and bmdls
bmds<-as.data.frame(unlist(bmds))
bmds<-tibble::rownames_to_column(bmds, "Simulation")
bmds$Simulation = substr(bmds$Simulation,1,nchar(bmds$Simulation)-1)
colnames(bmds)[2]<-'BMDSEstimates'
bmds$order<-bmds_order
bmds$Model<-'ModelAverage'


#get bmds average
tab_bmds<-bmds%>%
  group_by(order) %>% 
  summarize(med=median(BMDSEstimates), SD=sd(BMDSEstimates))
tab_bmds


#run mcmc over top 3 models
mcmc_ll<- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "gamma",fit_type = "mcmc"))
mcmc_p <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-probit",fit_type = "mcmc"))
mcmc_q <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "weibull",fit_type = "mcmc"))


#try plotting with built-in plot function
# out_mcmc_result<-mcmc_ll$exp_n1$mcmc_result
# out_mcmc_fittedmodel<-mcmc_ll$exp_n1$fitted_model
# out_mcmc_prior <-mcmc_ll$exp_n1$prior
# out_bmd_samples<-as.data.frame(mcmc_ll$exp_n1$mcmc_result$BMD_samples)
# plot(mcmc_ll$exp_n1)
# df<-mcmc_ll$exp_n1
# plot(test$V1~test$V2) #doesn't seem to be original plot; goes above 0.5 (plotly plot doesn't)

ggplot_fit<-plot(mcmc_ma_g$exp_n1)
ggplot_fit + scale_x_continuous(trans = "pseudo_log") + theme_classic()

MAdensity_plot(mcmc_ma$exp_n1)+ggtitle("Density Plot MCMC")


#organize and put outputs in dataframe
ll<-lapply(mcmc_ll, function (x) x['bmd'])
ll<-as.data.frame(unlist(ll))
ll<-tibble::rownames_to_column(ll, "Simulation")
ll$Simulation = substr(ll$Simulation,1,nchar(ll$Simulation)-1)
colnames(ll)[2]<-'BMDSEstimates'
ll$order<-bmds_order
ll$Model<-'Gamma'

p<-lapply(mcmc_p, function (x) x['bmd'])
p<-as.data.frame(unlist(p))
p<-tibble::rownames_to_column(p, "Simulation")
p$Simulation = substr(p$Simulation,1,nchar(p$Simulation)-1)
colnames(p)[2]<-'BMDSEstimates'
p$order<-bmds_order
p$Model<-'Log-Probit'

q<-lapply(mcmc_q, function (x) x['bmd'])
q<-as.data.frame(unlist(q))
q<-tibble::rownames_to_column(q, "Simulation")
q$Simulation = substr(q$Simulation,1,nchar(q$Simulation)-1)
colnames(q)[2]<-'BMDSEstimates'
q$order<-bmds_order
q$Model<-'Weibull'

bmds_est<-rbind(bmds,ll,p,q)
write.csv(bmds_est,'data_out/BMDS_headline_top3.csv')


############################Plots---- 
#toxicity data

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

bw_effect<-effects$Body_Weight_g
exp<-bw_effect^mod_As$Exponent 
as_d<-mod_As$Modifier*exp
derm_d<-effects$dermaldose 
app_d<-as.numeric(effects$Application_Rate)
derm_original<-derm_d*((as_d*app_d)/2)

effectsg$dermaldose<-derm_original


#combine both datasets
#names(effectsg)[11]<-"Body_Weight_g"
#effects<-rbind(effectsg[,c(1:18,20:22)], effectsg[,1:21])

# Figure 1: application rate and 96hr mortality figure?
# Figure 2: Distribution of posterior probabilities across 9 models to show why we picked them
# Figure 3: KDEs of BMDs across top 3 models
# Figure 4: DR curve for top 3 models and 95% credibility interval + median, plotted with original dataset

#run pyraclostrobin in Data_Simulation.R
effectsp<-effects
effectsp$label<-paste(effectsp$Study, effectsp$Species, sep="")
effectsp <-effectsp[order(effectsp$label),]
unique(effectsp$label) # get unique combos of study and species

# Figure 1: application rate and 96hr mortality figure?
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#762a83", "#000000","#7fbc41")

pp<-ggplot(effectsp, aes(x=dermaldose, y=Mortality, colour=label,shape=label)) +
  geom_point(size=2) + 
  scale_shape_manual(name = "Study and Species", 
    labels = c("Belden et al. 2010, Anaxyrus cognatus","Bruhl et al. 2013, Rana temporaria","Cusaac et al. 2015, Anaxyrus woodhousii", "Cusaac et al. 2015, Acris blanchardi",
               "Cusaac et al. 2016a, Acris blanchardi", "Cusaac et al. 2016b, Acris blanchardi",
               "Cusaac et al. 2017a, Anaxyrus cognatus", "Cusaac et al. 2017a, Anaxyrus woodhousii", "Cusaac et al. 2017b, Anaxyrus cognatus",
               "Cusaac et al. 2017c, Anaxyrus cognatus", "Cusaac et al. 2017d, Anaxyrus cognatus"),
    values=c(3:4,7:12,15,17,19))+
  scale_color_manual( name = "Study and Species",
                      labels = c("Belden et al. 2010, Anaxyrus cognatus","Bruhl et al. 2013, Rana temporaria","Cusaac et al. 2015, Anaxyrus woodhousii", "Cusaac et al. 2015, Acris blanchardi",
                                 "Cusaac et al. 2016a, Acris blanchardi", "Cusaac et al. 2016b, Acris blanchardi",
                                 "Cusaac et al. 2017a, Anaxyrus cognatus", "Cusaac et al. 2017a, Anaxyrus woodhousii", "Cusaac et al. 2017b, Anaxyrus cognatus",
                                 "Cusaac et al. 2017c, Anaxyrus cognatus", "Cusaac et al. 2017d, Anaxyrus cognatus"),
    values=colorBlindGrey8)+
  # geom_smooth(method=lm)+
  ggtitle("Adjusted 96hr Mortality vs. Calculated Dermal Dose, Pyaclostrobin") +
  ylab("Mortality") +
  xlab("Tissue Concentration (ug/g)")+
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
effectsg<-effects 
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


#Figure 2 - plot of posterior probability weights across 9 models 
ggplot(postw, aes(x = PosteriorProbs, y = Model, group = Model, fill = Model)) +
  geom_density_ridges(stat = "binline", bins = 100, scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Model") +
  xlab("Posterior Probability by Model for 1000 Simulated Studies")+
  theme_ridges()


#Figure 3 - KDE of benchmark doses 
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


#Figure 4 - we need to pull out the individual DR plots for the top models?
mcmc_ma$exp_n1$
plot(mcmc_ma$exp_n1)
plot(mcmc_q$exp_n1)
plot(mcmc_ma$exp_n1$Individual_Model_8)
derm_seq <-as.data.frame(seq(0,1.4, length=1000))  #sequence from 0 to 1.4 by x, 1000 length
colnames(derm_seq)[1]<-'dermdose'


save.image(file='myEnvironment_analysis.RData') 


#Andy extra analysis----
#do laplace and mle estimation as well, for Andy
D1 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "logistic",fit_type = "laplace"))
D2 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "logistic",fit_type = "mle"))
D3 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "logistic",fit_type = "mcmc"))

E1 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-probit",fit_type = "laplace"))
E2 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-probit",fit_type = "mle"))
E3 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-probit",fit_type = "mcmc"))

H1 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "weibull",fit_type = "laplace"))
H2 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "weibull",fit_type = "mle"))
H3 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "weibull",fit_type = "mcmc"))

I1 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-logistic",fit_type = "laplace"))
I2 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-logistic",fit_type = "mle"))
I3 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-logistic",fit_type = "mcmc"))

J1 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "qlinear",fit_type = "laplace"))
J2 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "qlinear",fit_type = "mle"))
J3 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "qlinear",fit_type = "mcmc"))

K1 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "probit",fit_type = "laplace"))
K2 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "probit",fit_type = "mle"))
K3 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "probit",fit_type = "mcmc"))

L1 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "multistage",fit_type = "laplace"))
L2 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "multistage",fit_type = "mle"))
L3 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "multistage",fit_type = "mcmc"))


M1 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "gamma",fit_type = "laplace"))
M2 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "gamma",fit_type = "mle"))
M3 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "gamma",fit_type = "mcmc"))


N1 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "hill",fit_type = "laplace"))
#Note I generally wouldn't run the hill MLE model!
N2 = lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "hill",fit_type = "mcmc"))


save.image(file='myEnvironment_andy_analysis.RData') 
load('myEnvironment_andy_analysis.RData')
output<-list(D1,D2,D3,E1,E2,E3,H1,H2,H3,I1,I2,I3,J1,J2,J3,K1,K2,K3,L1,L2,L3,M1,M2,M3,N1,N2)

out<-sapply(output,function(x)sapply(x, function(y) y['bmd']))
out<-bmdsorder<-c('bmds','bmdl','bmdu')
bmds_order<-rep(bmdsorder, times=26000)

model.names<-mod_list<-c('logistic','log-probit','weibull','log-logistic','qlinear','probit','multistage','gamma','hill')
model_order<-as.data.frame(rep(model.names,each=3))
model_order<-as.data.frame(model_order[1:26,])

model.types<-c('laplace','mle','mcmc')
model_types<-as.data.frame(rep(model.types,times=9))
model_types<-as.data.frame(model_types[1:26,])

models<-cbind(model_order, model_types)
models$fin_names<-paste(models$`model_order[1:26, ]`,"-",models$`model_types[1:26, ]`)
models<-models[,3]

test<-as.data.frame(out)
colnames(test)<-models
test$exp<-row.names(test)
testy<-tidyr::gather(test, "model", "n", 1:26)
testf<-testy %>% unnest(n) %>% group_by(model)
testf$dose_measure<-bmds_order

write.csv(testf, 'data_out/bmds_single_models.csv')
