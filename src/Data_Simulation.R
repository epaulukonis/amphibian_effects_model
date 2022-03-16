##script for simulation and analysis of dermal dose and effects
library(dplyr)
library(tidyr)
library(stats)
library(multcomp)
library(DescTools)
library(reshape2)
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')



# we run through the analysis and simulation with Pyraclostrobin and Glyphosate
# evaluating statistical significance using cochran-armitage test (Neurhauser et al. 1999)
# modify the exposure duration using a decay function to 96 hours


mod<-read.csv('data_in/As_Modifier.csv')


#### Pyraclostrobin ----
#note to self; specific column numbers will change depending on dataset
effects<-read.csv('data_in/Headline_updated.csv')
param<-read.csv('data_in/parameters.csv') #go back and take weighted mean; sum alive / sum assessed 
nsims<-1000
#set.seed(6379) 

#first modify so that the survival is adjusted for 96 hours across all studies
effects$half_life<-(effects$Duration_h*log(2))/log(1/effects$Survival) #need to modify survival by duration, using an exponential growth curve
effects$adj_sur_96<-round(1/(2^(96/effects$half_life)),3)
print(names(effects))
df_sig<-effects[,c(1,5,9,11,18)] #make sure this output matches
df_sig<-na.omit(df_sig)
df_sig$mort<-round(df_sig$N_Exp - (df_sig$Survival*df_sig$N_Exp),0)
df_sig$alive<-df_sig$N_Exp-df_sig$mort
head(df_sig)
#split data into list by study/species for cochran armitage trend test
by_s_stat<-split(df_sig, list(df_sig$Study, df_sig$Species), drop=T)
by_s_stat<-lapply(by_s_stat, function (x) x[c('Application_Rate','mort','alive')]) 

#change format to matrix to run CA trend test to test for trends in increasing significance
ch0=list()
for (i in seq_along(by_s_stat)) {
a=t(by_s_stat[[i]])
colnames(a)<-a[1,]
a<-a[c(3,2),] #need to have no mort first row, mort second
row.names(a)<-c(0,1)
a<-as.matrix(a)
ch0[[i]]<-a
}

#run CA test over matrix list; format needs to match format in ?Cochran-Armitage 'dose' ex
CA_increasing_trend<-lapply(ch0,function(x) CochranArmitageTest(x,'increasing')) 
CA_increasing_trend
#10,5,1 fail to reject the null hypothesis for pyraclostrobin
#non fail to reject the null hypothesis for glyphosate 

#for pryaclostrobin
#1 is cusaac 2015 acris blanchardi
#5 is cusaac 2017a anaxyrus cognatus
#10 is cusaac 2017a anaxyrus woodhousii
#low sample sizes for two (5,10), other two not much response from exposure, also lower application rates
#remove 4 studies for a total of 25 points including weighted control

effects_sub<-effects[(effects$Study =='Cusaac_2015' & effects$Species == 'Acris blanchardi'),]
effects<-effects[!(effects$Study == 'Cusaac_2017a' | effects$Study == 'Cusaac_2015'), ] 
effects<-rbind(effects,effects_sub)
effects<-effects[order(effects$Study),]

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
min(effects$M_Body_Weight_g)
min(effects$SD)
effects <-effects[order(effects$Study),]


###Body weights
#Simulate body weights (BWs)- we need 1000 simulations of the X different BWs
by_s<-split(effects, list(effects$Study,effects$Species,effects$Application_Rate), drop=T) #split by study, species, application rate use SE
by_s<-by_s[order(names(by_s))]
print(names(by_s[[1]]))

bw_sim<-replicate(nsims,
                  {normalize <- lapply(by_s, function(y) rnorm(nrow(y), mean=y[[12]], sd=y[[13]]))
                  bw_sim<-as.data.frame(unlist(normalize))
                  })

bw_sim<-do.call(cbind.data.frame, bw_sim)
colnames(bw_sim)[1:nsims]<-names<-paste0("BW",1:nsims,"")




###survival
#simulate survival using a binomial distribution with modified duration
survival_sim <- matrix(data=NA,nrow=nrow(effects),ncol=nsims)
colnames(survival_sim)[1:nsims]<-paste0("Sur",1:nsims,"")
names(effects)
for(i in 1:nrow(effects)){
  survival_sim[i,]<-rbinom(nsims, effects[i,17], effects[i,20]) #trials is number of observations for each experiment, prob is adjusted prob of survival for 96 hours
  survival_sim[i,]<-round(survival_sim[i,]/effects[i,17],3) #divide the number surviving by total trials to find the proportion that actually survived 
}
survival_sim<-as.data.frame(survival_sim)

###exposure parameters
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


###pyraclostrobin
#calculate dermal dosage 
pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.91*24 #from comptox profile

#calculate soil concentrations; application rate will remain the same, the exponent value changes with exposure parameter simulations
soil_concs_degs<-as.data.frame(log(2)/hl*96/exposure_sims$movement_rate_mean) #use 96 as duration, as we adjusted all survivals for a 96h timeframe
soil_concs<-((effects$Application_Rate*10)/16000)*1000 #1cm depth #converts it into final units of ug/g (tissue conc/bw)#this comes from Weir
my_list<-list()
for (i in 1:nsims){
  dsa<-exposure_sims[i,5]*bw_sim[,i]^exposure_sims[i,5]
  dermal_dose<-(soil_concs^soil_concs_degs[i,] * kp_pyra * (dsa/exposure_sims[i,1]) * exposure_sims[i,6] * exposure_sims[i,3])/bw_sim[,i]
  my_list[[i]]<-dermal_dose
}



my_list[9] #check the list outputs to be sure
derm<-as.data.frame(t(do.call(rbind.data.frame, my_list)))
colnames(derm)<-paste0("dermdose",1:nsims,"")
row.names(derm)<-NULL


#this calculates just the soil exposure, NOT overspray; we need to adjust it using the method I describe in the manuscript
#I.e, use allometric relationship, x that by the application rate, and that is the factor that should be used to modify the tissue concentration

#For each output; multiply the simulated bw by the As value; 
#then adjust the dermal dose by that much. 100% of As, 50% As, and 2/3 As.

fam<-as.data.frame(cbind(effects$Family, effects$Study))
names(fam)<-c('Family', 'Study')
mod_As<- merge(fam,mod, by  = "Family") 
mod_As<-mod_As[order(mod_As$Study),]


bw_sim_as<-bw_sim
exp<-bw_sim_as^mod_As$Exponent 
as<-mod_As$As_Modifyer*exp
colnames(as)[1:nsims]<-names<-paste0("As",1:nsims,"")

##100% of surface area 
derm_adj_r<-derm*as

#Check ratio:
ratio_100<-derm_adj_r/derm
ratio_100$Family<-effects$Family #add family
BCF_100_m<- ratio_100 %>% 
  group_by(Family) %>%
  summarise(across(dermdose1:dermdose1000,mean, na.rm = TRUE))
print(rowMeans(BCF_100_m[2:1001]))


##50% of surface area
derm_adj_r<-derm*(as/2)

#Check ratio:
ratio_50<-derm_adj_r/derm
ratio_50$Family<-effects$Family #add family
BCF_50_m<- ratio_50 %>% 
  group_by(Family) %>%
  summarise(across(dermdose1:dermdose1000, mean, na.rm = TRUE))
print(rowMeans(BCF_50_m[2:1001]))

##2/3 of surface area
derm_adj_r<-derm*(as*0.666)

#Check ratio:
ratio_66<-derm_adj_r/derm
ratio_66$Family<-effects$Family #add family
BCF_66_m<- ratio_66 %>% 
  group_by(Family) %>%
  summarise(across(dermdose1:dermdose1000, mean, na.rm = TRUE))
print(rowMeans(BCF_66_m[2:1001]))

####

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

BMDS_headline_fin<-cbind(exp_n,mort_n,derm_dose) #name the output whatever run you'd like
BMDS_headline_fin<-BMDS_headline_fin[,c(1,6,2,4)]

#order by dose for clarity
BMDS_headline_fin<-BMDS_headline_fin %>%
  group_by(set) %>% #set of sims
  arrange(Dose, .by_group=T) #this will rearrange the output in terms of order of sims, but will still work fine


#save.image(file='myEnvironment_simulation.RData') 
write.csv(BMDS_headline_fin,'data_out/BMDS_headline_fin_132022.csv') #save whole file for bmds analysis

#calculate means for paper here
#mortality
mort_sim<-1-survival_sim
mort_n_forpaper<-gather(mort_sim,"set","Effect", 1:nsims)
mean(mort_n_forpaper$Effect)

#dose
mean(BMDS_headline_fin$Dose)




#### Glyphosate ----
#note to self; specific column numbers will change depending on dataset
effects<-read.csv('data_in/Glyphosate_updated.csv')
param<-read.csv('data_in/parameters.csv') #go back and take weighted mean; sum alive / sum assessed 
nsims<-1000
#set.seed(6379) 

#the duration for all of the glyphosate studies is already set to 96 hours; not need to calculate it using exponential growth curve
effects$adj_sur_96<-effects$Survival
print(names(effects))
df_sig<-effects[,c(1,5,9,11,18)] #make sure this output matches
df_sig<-na.omit(df_sig)
df_sig$mort<-round(df_sig$N_Exp - (df_sig$Survival*df_sig$N_Exp),0)
df_sig$alive<-df_sig$N_Exp-df_sig$mort
head(df_sig)
#split data into list by study/species for cochran armitage trend test
by_s_stat<-split(df_sig, list(df_sig$Study, df_sig$Species), drop=T)
by_s_stat<-lapply(by_s_stat, function (x) x[c('Application_Rate','mort','alive')]) 

#change format to matrix to run CA trend test to test for trends in increasing significance
ch0=list()
for (i in seq_along(by_s_stat)) {
  a=t(by_s_stat[[i]])
  colnames(a)<-a[1,]
  a<-a[c(3,2),] #need to have no mort first row, mort second
  row.names(a)<-c(0,1)
  a<-as.matrix(a)
  ch0[[i]]<-a
}

#run CA test over matrix list; format needs to match format in ?Cochran-Armitage 'dose' ex
CA_increasing_trend<-lapply(ch0,function(x) CochranArmitageTest(x,'increasing')) 
CA_increasing_trend
#none fail to reject the null hypothesis for glyphosate 


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
effects<-na.omit(effects)
#effects <-effects[order(effects$Study),]
min(effects$Body_Weight_g)
max(effects$Body_Weight_g)
min(effects$SD)

min(effects$Survival)
max(effects$Survival)

unique(effects$Species)

min(effects$Application_Rate)
max(effects$Application_Rate)


#format data for sampling from lognormal distribution
effects$bw_m<-log(effects$Body_Weight_g^2 / sqrt(effects$SD^2 + effects$Body_Weight_g^2))
effects$bw_sd<-sqrt(log(1 + (effects$SD^2 / effects$Body_Weight_g^2)))


##let's bootstrap some data based on our original dataset
#Simulate body weights (BWs)- we need 1000 simulations of the X different BWs
by_s<-split(effects, list(effects$Study, effects$Species,effects$Application_Rate), drop=T) #split by study, species, application rate
print(names(by_s[[19]]))

bw_sim<-replicate(nsims,
                  {normalize <- lapply(by_s, function(y) rnorm(nrow(y), mean=y[[12]], sd=y[[13]]))
                  bw_sim<-as.data.frame(unlist(normalize))
                  })
bw_sim<-do.call(cbind.data.frame, bw_sim)
colnames(bw_sim)[1:nsims]<-names<-paste0("BW",1:nsims,"")

#simulate survival using a binomial distribution with modified duration
print(names(effects))
survival_sim <- matrix(data=NA,nrow=nrow(effects),ncol=nsims)
colnames(survival_sim)[1:nsims]<-paste0("Sur",1:nsims,"")
for(i in 1:nrow(effects)){
  survival_sim[i,]<-rbinom(nsims, effects[i,17], effects[i,19]) #trials is number of observations for each experiment, prob is adjusted prob of survival for 96 hours
  survival_sim[i,]<-round(survival_sim[i,]/effects[i,17],3) #divide the number surviving by total in exp to find the proportion that actually survived 
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


#glyphosate
#calculate dermal dosage 
pmolweight<-169.1 #g/mol #comptox
plogKow<- -3.12 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.46*24 #from comptox profile


#calculate soil concentrations; application rate will remain the same, the exponent value changes with exposure parameter simulations
soil_concs_degs<-as.data.frame(log(2)/hl*96/exposure_sims$movement_rate_mean) #use 96 as duration, as we adjusted all survivals for a 96h timeframe
soil_concs<-((effects$Application_Rate*10)/16000)*1000 #1cm depth #converts it into final units of ug/g (tissue conc/bw)#this comes from Weir
my_list<-list()
for (i in 1:nsims){
  dsa<-exposure_sims[i,5]*bw_sim[,i]^exposure_sims[i,5]
  dermal_dose<-(soil_concs^soil_concs_degs[i,] * kp_pyra * (dsa/exposure_sims[i,1]) * exposure_sims[i,6] * exposure_sims[i,3])/bw_sim[,i]
  my_list[[i]]<-dermal_dose
}
my_list[9] #check the list outputs to be sure
derm<-as.data.frame(t(do.call(rbind.data.frame, my_list)))
colnames(derm)<-paste0("dermdose",1:nsims,"")
row.names(derm)<-NULL
#outputs are in ugg


fam<-as.data.frame(cbind(effects$Family, effects$Study))
names(fam)<-c('Family', 'Study')
mod_As<- merge(fam,mod, by  = "Family") 
mod_As<-mod_As[order(mod_As$Study),]


bw_sim_as<-bw_sim
exp<-bw_sim_as^mod_As$Exponent 
as<-mod_As$As_Modifyer*exp
colnames(as)[1:nsims]<-names<-paste0("As",1:nsims,"")

##100% of surface area 
derm_adj_r<-derm*as

#Check ratio:
ratio_100<-derm_adj_r/derm
ratio_100$Family<-effects$Family #add family
BCF_100_m<- ratio_100 %>% 
  group_by(Family) %>%
  summarise(across(dermdose1:dermdose1000,mean, na.rm = TRUE))
print(rowMeans(BCF_100_m[2:1001]))


##50% of surface area
derm_adj_r<-derm*(as/2)

#Check ratio:
ratio_50<-derm_adj_r/derm
ratio_50$Family<-effects$Family #add family
BCF_50_m<- ratio_50 %>% 
  group_by(Family) %>%
  summarise(across(dermdose1:dermdose1000, mean, na.rm = TRUE))
print(rowMeans(BCF_50_m[2:1001]))

##2/3 of surface area
derm_adj_r<-derm*(as*0.666)

#Check ratio:
ratio_66<-derm_adj_r/derm
ratio_66$Family<-effects$Family #add family
BCF_66_m<- ratio_66 %>% 
  group_by(Family) %>%
  summarise(across(dermdose1:dermdose1000, mean, na.rm = TRUE))
print(rowMeans(BCF_66_m[2:1001]))


####

#format the dermal dose estimates and survival estimates to use in the BBMD app
#should be in order with effects, so that each experiment N matches with each study
mort_n<-as.data.frame(sapply(survival_sim, function(x) effects$N_Exp - round(x*effects$N_Exp,0))) #number of individuals dying per exp
#we are treating each batch like its own separate study; essentially, each set of doses has a corresponding set of survival data
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
exp_n<-as.data.frame(rep.col(effects$N_Exp, nsims))
colnames(exp_n)<-paste0("exp_n",1:nsims,"")
#need ID, dose, N, and mortality
#for 100 sets of 39 rows
exp_n<-gather(exp_n,"set","N",1:nsims)
mort_n<-gather(mort_n,"set","Effect", 1:nsims)
derm_dose<-gather(derm,"set","Dose",1:nsims)

BMDS_glyphosate_fin<-cbind(exp_n,mort_n,derm_dose) #name the output whatever run you'd like
BMDS_glyphosate_fin<-BMDS_glyphosate_fin[,c(1,6,2,4)]

#order by dose for clarity
BMDS_glyphosate_fin<-BMDS_glyphosate_fin %>%
  group_by(set) %>% #set of sims
  arrange(Dose, .by_group=T) #this will rearrange the output in terms of order of sims, but will still work fine


#save.image(file='myEnvironment_simulation.RData') 
write.csv(BMDS_glyphosate_fin,'data_out/BMDS_glyphosate_fin_132022.csv') #save whole file for bmds analysis

#calculate means for paper here
#mortality
mort_sim<-1-survival_sim
mort_n_forpaper<-gather(mort_sim,"set","Effect", 1:nsims)
mean(mort_n_forpaper$Effect)

#dose
mean(BMDS_glyphosate_fin$Dose)



