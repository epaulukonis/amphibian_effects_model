### April 2022 Analysis

###Evaluating MCMC MA runs, and pulling out 'top' model using BMD


###5/25/22 USE LAPLACE CAUSE THE ALPHA IS WEIRD

library(Rcpp)
library(RcppGSL)
library(RcppEigen)
library(forcats)
library(plotly)
library(shiny)
library(scales)
library(ToxicR) ## 
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
library(drc)

set.seed(6379)
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')

#simulated datasets
sims<-read.csv('data_out/BMDS_headline_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T) #split by simulation

##og dataset
effects<-read.csv('data_in/Headline_updated.csv') #run code in Data_Sim to clean up the effects and remove rows we don't want
mod<-read.csv('data_in/As_Modifier.csv')
param<-read.csv('data_in/parameters.csv')

####OG MA analysis----
# 
# mod_list<-sort(c('hill','log-logistic','logistic','log-probit','weibull','qlinear','probit','multistage','gamma')) #order the models by name
# mcmc_ma = lapply(by_s, function(y) ma_dichotomous_fit(y[,2],y[,4],y[,3], fit_type = "mcmc")) #apply ma function over list
# 
# post<-lapply(mcmc_ma, function (x) x['posterior_probs']) #pull out posterior probabilities 
# postw<-as.data.frame(unlist(post))
# postw<-tibble::rownames_to_column(postw, "Simulation")
# postw$Simulation = substr(postw$Simulation,1,nchar(postw$Simulation)-1) #remove extraneous characters
# colnames(postw)[2]<-'PosteriorProbs'
# mod_list_f<-rep(mod_list, times=1000)
# postw$Model<-mod_list_f #add column for model
# 
# #look at median of PW across all models 
# tab_med<-postw %>%
#   group_by(Model) %>% 
#   summarize(med=median(PosteriorProbs))
# print(tab_med)
#write.csv(tab_med,'data_out/BMDS_posteriorweight_median_headline.csv')

##the log-logistic has the highest posterior probability weight


##here we model the ll model (top)
#model the log-logistic single model using the simulated sets
ll_fit <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type="log-logistic",fit_type="laplace"))


# #pull model average BMD, bmdl, bmdu
bmdsorder<-c('bmds','bmdl','bmdu')
bmds_order<-rep(bmdsorder, times=1000)
bmds<-lapply(ll_fit, function (x) x['bmd']) #pull out bmds and bmdls
bmds<-as.data.frame(unlist(bmds))
bmds<-tibble::rownames_to_column(bmds, "Simulation")
bmds$Simulation = substr(bmds$Simulation,1,nchar(bmds$Simulation)-1)
colnames(bmds)[2]<-'BMDSEstimates'
bmds$order<-bmds_order
bmds$Model<-'log_logistic'


#get parameters 
# para<-lapply(ll_fit, function (x) x['parameters'])
# para<-as.data.frame(unlist(para))
# para<-tibble::rownames_to_column(para, "Exp")
# colnames(para)[2]<-'Value'
# para_list_f<-rep(c("p1","p2","p3"), times=1000)
# para$Parameters<-para_list_f


## fit og pyraclostrobin data ----
#first modify so that the survival is adjusted for 96 hours across all studies
effects$half_life<-(effects$Duration_h*log(2))/log(1/effects$Survival) #need to modify survival by duration, using an exponential growth curve
effects$adj_sur_96<-round(1/(2^(96/effects$half_life)),3)
effects_sub<-effects[(effects$Study =='Cusaac_2015' & effects$Species == 'Anaxyrus woodhousii'),]
effects<-effects[!(effects$Study == 'Cusaac_2017a' | effects$Study == 'Cusaac_2015'), ] 
effects<-rbind(effects,effects_sub)
#effects<-effects[order(effects$Study),]
#group by study, species, and application rate
effects<-effects %>%
  group_by(Study) %>% #study
  group_by(Species) %>% #species
  group_by(Application_Rate)%>% #application rate
  mutate(SE = SD/sqrt(length(effects))) #mutate to add standard error to account for variability between studies
effects<-as.data.frame(effects)
#effects<-na.omit(effects)
controls<-effects[effects$Application_Rate == 0, ]  #take out all controls
control<-effects[1,] # take out single row which will contain final control
#function for weighted mean
weight_mean<-function(data,vector,weight){
  round(with(data, sum(vector*weight)/sum(weight)),2)
}
control$Survival<-weight_mean(controls,controls$Survival, controls$N_Exp) #survival
control$N_Exp<-sum(controls$N_Exp)#use sum of the N used in dose for final N_exp
effects<-effects[!effects$Application_Rate == 0, ] #remove application rate of 0
effects<-rbind(effects,control)
effects <-effects[order(effects$Study),]


#then, format modifier DF
fam<-as.data.frame(cbind(effects$Family, effects$Study, effects$Application_Rate))
names(fam)<-c('Family', 'Study', "Application")
mod_As<- merge(fam,mod, by  = "Family") 
mod_As<-mod_As[order(mod_As$Study),]

#this is for the in-house dermal doses, we'll calculate it this way:
pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.91*24 #from comptox profile
dt<-param[1,2]
move_rate<-param[3,2]
bioavail<-param[5,2]
dsa<-param[9,2]*effects$M_Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_deg<-log(2)/hl*96/move_rate #with adjusted time of 96h to match adjusted survival
soil_concs<-((effects$Application_Rate*10)/16000)*1000 #1cm mixing depth; or change depending on soil?
dermal_dose<-(soil_concs^soil_concs_deg * kp_pyra * (dsa/dt) * derm_frac * bioavail)/effects$M_Body_Weight_g
#head(dermal_dose) #take a look
effects$dermaldose<-dermal_dose
effects$Mortality<-1-effects$adj_sur_96
bw_effect<-effects$M_Body_Weight_g
exp<-bw_effect^mod_As$Exponent 
as<-mod_As$Modifier*exp
derm_d<-effects[c(1:23),23] #only do the direct exposures
app_d<-effects[c(1:23),9] #only do the direct exposures
as<-as[1:23] #only do the direct exposures
derm_original<-derm_d*((as*as.numeric(app_d))/2) #modify the calculated dermal dose by the SA and application rate product
derm_original<-c(derm_original,effects[c(24:29),23])
effects$dermaldose<-derm_original

effects<-effects[with(effects, order(Study, Species, Application_Rate)), ]
#effects<-effects[order(effects$dermaldose),]
effects$death<-round(effects$Mortality*effects$N_Exp,0)
effects$Exp<-ifelse(effects$Method=="Overspray",1,0)
#check order
print(effects$death)
by_s[[748]]$Effect



og_data<-matrix(0,nrow=nrow(effects), ncol=4)
colnames(og_data) <- c("Dose","N","Incidence", "Exp")
og_data[,1] <- effects$dermaldose
og_data[,2] <- effects$N_Exp
og_data[,3] <- effects$death
og_data[,4] <- effects$Exp


#### time to get the upper and lower bounds out of the 1000 curves ----
#post-process all data for credible intervals
d<-og_data[,1]
nd<-length(by_s[[1]]$Dose)
nsims <- 1000
mortality_df <- matrix(ncol = nsims, nrow = nd)

for(i in 1:length(by_s)){
  parms <- ll_fit[[i]]$parameters
  g <- 1/(1+exp(-parms[1])); 
  a <- parms[2];
  b <- parms[3]; 
  mortality_df[,i] <- g + (1-g)*(1/(1+exp(-a-b*log(d)))) #or by d or by_s[[i]]$Dose
} 
View(mortality_df)


percentiles_df <- matrix(ncol = 2, nrow = nd)
for(i in 1:nd){
  mort_percentiles <- quantile(mortality_df[i,], probs = c(0.025,0.975), names=T)
  percentiles_df[i,1] <- mort_percentiles[[1]]
  # percentiles_df[i,2] <- mort_percentiles[[2]]
  # percentiles_df[i,3] <- mort_percentiles[[3]]
  # percentiles_df[i,4] <- mort_percentiles[[4]]
  # percentiles_df[i,5] <- mort_percentiles[[5]]
  # percentiles_df[i,6] <- mort_percentiles[[6]]
  percentiles_df[i,2] <- mort_percentiles[[2]]
}
dim(percentiles_df)

df_percentiles <- data.frame(x=d, val= as.vector(percentiles_df), 
                             variable=rep(paste0("category", 1:2), each=29))


#### plot curve with original data as well as upper and lower bounds (2.5% and 97.5% percentiles of mortality), with OG points ----
ll_laplace<-single_dichotomous_fit(og_data[,1],og_data[,3],og_data[,2],model_type="log-logistic",fit_type="laplace")


parms <- ll_laplace$parameters
g <- 1/(1+exp(-parms[1])); 
a <- (parms[2]);
b <- parms[3]; 
d <- effects$dermaldose
rval_og <- g + (1-g)*(1/(1+exp(-a-b*log(d))))
df<-as.data.frame(cbind(rval_og,d))
names(df)<-c("effect","dose")

lerror<-df_percentiles[df_percentiles$variable == "category1",2]
uerror<-df_percentiles[df_percentiles$variable == "category2",2]

sort(df$effect)


og_data<-as.data.frame(og_data)

my_breaks_m<-c(0,10,20,30,40,50,55)
main<-ggplot() + 
 #geom_density(data=bmds_df, aes(x=bmds))+
  geom_line(data = df, aes(x=dose, y=effect))+
  geom_ribbon(aes(x = df$dose, ymin =lerror, ymax = uerror), alpha = .2) +
  geom_point(data=og_data, aes(x=Dose,y=(Incidence/N), colour=Exp, fill=Exp))+
  ggtitle("Log-logistic Pyraclostrobin Curve") +
  ylab("Mortality") +
  xlab("Dose (ug/g)")+
  scale_x_continuous(limits=c(0,55),breaks=my_breaks_m,expand = c(0, 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"), 
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
      axis.text.y = element_text(size=12, face='bold'),
      axis.title.x = element_text(size=14, face='bold'),
      axis.title.y = element_text(size=14, face='bold'),
      plot.title = element_text(face = 'bold', size = 16), legend.position = 'none') #first run without this to get legend

main


#range_webice<-c(2.165,172.912)

# main+geom_point(aes(x=6,y=0.50),colour="red")+
#   geom_pointrange(aes(x=6, y=0.50,xmin=2.165, xmax=172.912))

ld50<-read.csv('data_in/ld50_pyra.csv')
#ld50<-ld50[!ld50$LD50>100,]

cons<-sort(ld50$LD50)[1:5]
ld50c<-ld50[ld50$LD50 %in% cons,]
frog<-6.288521517
frog<-ld50[ld50$LD50 %in% frog,]


# main +
#   geom_boxplot(data=ld50c, aes(x=LD50, y=Mortality), width=0.05) +
#   geom_vline(xintercept=frog$LD50, linetype='dotted', col = 'red')+
#   # geom_point(data=frog, aes(x=LD50,y=Mortality),colour="black", shape=18, size=4) +
#   geom_point(aes(x=6,y=0.50),colour="red") 
#   
  #geom_pointrange(aes(x=6, y=0.50,xmin=2.165, xmax=172.912))

#bmds_est<-bmds[bmds$order =='bmds',]


# nsims <- 1000
# bmds_df <- matrix(nrow = nsims, ncol = 2)
# for(i in 1:length(ll_fit)){
#   parms <- ll_fit[[i]]$parameters
#   g <- 1/(1+exp(-parms[1])); 
#   a <- parms[2];
#   b <- parms[3]; 
#   d <- bmds_est[i,2]
#  bmds_df[i,2] <- g + (1-g)*(1/(1+exp(-a-b*log(d)))) #or by d or by_s[[i]]$Dose
#  bmds_df[i,1]<-d
# } 
# 
# bmds_df<-as.data.frame(bmds_df)
# names(bmds_df)<-c("bmds","mortality")

my_breaks<-c(0,1,2,3,4,5,10,20,30,40,50, 55)


# parms <- ll_laplace$parameters
# g <- 1/(1+exp(-parms[1])); 
# a <- (parms[2]);
# b <- parms[3]; 
# d <- bmds$BMDSEstimates
# rval_og <- g + (1-g)*(1/(1+exp(-a-b*log(d))))
# df_bmd<-as.data.frame(cbind(rval_og,d))
# names(df_bmd)<-c("bmdseffect","bmdsdose")


bmds_hist_py<-
  ggplot(data=bmds,aes(x=BMDSEstimates, y=..count..,group=order,fill=order))+
  geom_histogram(bins=300)+
  ylab("Density BMDs") +
  xlab("Dose (ug/g)")+
  scale_x_continuous(limits=c(0,55),breaks=my_breaks, expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16), legend.position = "none")

bmds_hist_py


final_glypho<-main +
  geom_boxplot(data=ld50c, aes(x=LD50, y=Mortality), width=0.05) +
  geom_vline(xintercept=frog[1,1], linetype='dotted', col = 'red')+
  geom_vline(xintercept=frog[2,1], linetype='dotted', col = 'red')+
  # geom_point(data=frog, aes(x=LD50,y=Mortality),colour="black", shape=18, size=4) +
  geom_point(aes(x=0.017,y=0.50),colour="red") 


grid.arrange(final_glypho, bmds_hist_py,nrow=2,heights=c(4.5,1.5))



