##script for simulation and analysis of dermal dose and effects
library(dplyr)
library(tidyr)

load('myEnvironment_simulation.RData')
effects<-read.csv('data_in/Headline_updated.csv')
param<-read.csv('data_in/parameters.csv') #go back and take weighted mean; sum alive / sum assessed 
nsims<-1000
set.seed(6379) 

#group by study, species, and application rate
effects<-effects %>%
  group_by(Study) %>% #study
  group_by(Species) %>% #species
  group_by(Application_Rate)%>% #application rate
  mutate(SE = SD/sqrt(length(effects))) #mutate to add standard error to account for variability between studies
effects<-as.data.frame(effects)


str(effects)
#N_exp is int
#adj_sur_96 is num


controls<-effects[effects$Application_Rate == 0, ]  #take out all controls
control<-effects[5,] # take out single row which will contain final control
#function for weighted mean
weight_mean<-function(data,vector,weight){
  round(with(data, sum(vector*weight)/sum(weight)),2)
}
control$Survival<-weight_mean(controls,controls$Survival, controls$N_Exp) #survival
control$N_Exp<-sum(controls$N_Exp)#use sum of the N used in dose for final N_exp
effects<-effects[effects$Application_Rate != 0, ]
effects<-rbind(effects,control)
effects <-effects[order(effects$Study),]


##let's bootstrap some data based on our original dataset
#Simulate body weights (BWs)- we need 1000 simulations of the X different BWs
by_s<-split(effects, list(effects$Study,effects$Species,effects$Application_Rate), drop=T) #split by these groups, use SE
bw_sim<-replicate(nsims,
                  {normalize <- lapply(by_s, function(y) rnorm(nrow(y), mean=y[[11]], sd=y[[19]]))
                  bw_sim<-as.data.frame(unlist(normalize))
                  })
bw_sim<-do.call(cbind.data.frame, bw_sim)
colnames(bw_sim)[1:nsims]<-names<-paste0("BW",1:nsims,"")


#simulate survival using a binomial distribution with modified duration
effects$half_life<-(effects$Duration_h*log(2))/log(1/effects$Survival)#need to modify survival by duration, using an exponential growth curve
effects$adj_sur_96<-round(1/(2^(96/effects$half_life)),3)

survival_sim <- matrix(data=NA,nrow=nrow(effects),ncol=nsims)
colnames(survival_sim)[1:nsims]<-paste0("Sur",1:nsims,"")
for(i in 1:nrow(effects)){
  survival_sim[i,]<-rbinom(nsims, effects[i,17], effects[i,21]) #trials is number of observations for each experiment, prob is adjusted prob of survival for 96 hours
  survival_sim[i,]<-round(survival_sim[i,]/effects[i,17],3) #divide the number surviving by total trials to find the proportion that actually survived 
}
survival_sim<-as.data.frame(survival_sim)

#simulate the exposure parameters from Purucker et al. 2021
exposure_sims <- matrix(data=NA,nrow=nsims,ncol=6)
colnames(exposure_sims) <- c("dt_mean","movement_rate_mean","bioavail_mean",
                             "dermal_sa_slope_mean","dermal_sa_exponent_mean",
                             "dermal_fraction_mean")
exposure_sims[,1]<-rnorm(nsims, mean=param[1,2], sd=param[2,2])
exposure_sims[,2]<-rpois(nsims, param[3,2])                     #recall this is a poisson distribution, so no SD
exposure_sims[,3]<-rnorm(nsims, mean=param[5,2], sd=param[6,2])
exposure_sims[,4]<-rnorm(nsims, mean=param[7,2], sd=param[8,2])
exposure_sims[,5]<-rnorm(nsims, mean=param[9,2], sd=param[10,2])
exposure_sims[,6]<-rnorm(nsims, mean=param[11,2], sd=param[12,2])
exposure_sims<-as.data.frame(exposure_sims)


#calculate dermal dosage 
pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.91*24 #from comptox profile

#calculate soil concentrations; application rate will remain the same, the exponent value changes with exposure parameter simulations
soil_concs_degs<-as.data.frame(log(2)/hl*96/exposure_sims$movement_rate_mean) #use 96 as duration, as we adjusted all survivals for a 96h timeframe
soil_concs<-((effects$Application_Rate*10)/16000)*1000 #1cm #converts it into units of ugg #this comes from Weir
my_list<-list()
for (i in 1:nsims){
  dsa<-exposure_sims[i,5]*bw_sim[,i]^exposure_sims[i,5]
  dermal_dose<-(soil_concs^soil_concs_degs[i,] * kp_pyra * (dsa/exposure_sims[i,1]) * exposure_sims[i,6] * exposure_sims[i,3])/bw_sim[,i]
  my_list[[i]]<-dermal_dose
}
my_list[3] #check the list outputs to be sure
derm<-as.data.frame(t(do.call(rbind.data.frame, my_list)))
colnames(derm)<-paste0("dermdose",1:nsims,"")
row.names(derm)<-NULL
#outputs are in ugg

#format the dermal dose estimates and survival estimates to use in the BBMD app
#should be in order with effects, so that each experiment N matches with each study
mort_n<-as.data.frame(sapply(survival_sim, function(x) effects$N_Exp - round(x*effects$N_Exp,0))) #number of individuals dying per exp
#we are treating each batch like its own separate study; essentially, each set of doses has a corresponding set of survival data
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
exp_n<-as.data.frame(rep.col(effects$N_Exp,nsims))
colnames(exp_n)<-paste0("exp_n",1:nsims,"")
#need ID, dose, N, and mortality
#for 100 sets of 39 rows
exp_n<-gather(exp_n,"set","N",1:nsims)
mort_n<-gather(mort_n,"set","Effect", 1:nsims)
derm_dose<-gather(derm,"set","Dose",1:nsims)

BMDS_Run3<-cbind(exp_n,mort_n,derm_dose) #name the output whatever run you'd like
BMDS_Run3<-BMDS_Run3[,c(1,6,2,4)]

#BBMD requires each ID group to be in descending order
BMDS_Run3<-BMDS_Run3 %>%
  group_by(set) %>% #set of sims
  arrange(Dose, .by_group=T) #this will rearrange the output in terms of order of sims, but will still work fine
#we need to split up the 1000 set simulation into 10 sets of 100
bbmd<-split(BMDS_Run3, (as.numeric(rownames(BMDS_Run3))-1) %/% 3900)
lapply(bbmd,dim) #check dimensions of the output; should be 3900*4 for each set of 10 (0-9)

#quick test if needed
# test<-bbmd
# testy<-test[[3]] #check that it's formatted with 100 sets of unique sims
# unique(testy$set) #should be 100

colnames<-c("Dataset Index", "Dose", "N","Effect")
bbmd<-lapply(bbmd, setNames, colnames)
##for BBMD format
# for(i in names(bbmd)){
#   write.csv(bbmd[[i]], paste0(i,"bbmdrun.csv"), row.names=FALSE)
# }

save.image(file='myEnvironment_simulation.RData') 
write.csv(BMDS_Run3,'data_out/BMDS_Run3.csv') #save whole file for bmds analysis

