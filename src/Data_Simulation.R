##script for simulation and analysis of dermal dose and effects
library(dplyr)
library(tidyr)
library(stats)
library(multcomp)
library(DescTools)
library(reshape2)

#load('myEnvironment_simulation.RData')
#note to self; specific column numbers will change depending on dataset
effects<-read.csv('data_in/Glyphosate_updated.csv')
param<-read.csv('data_in/parameters.csv') #go back and take weighted mean; sum alive / sum assessed 
nsims<-1000
set.seed(6379) 

#we run through the analysis and simulation with Pyraclostrobin and Glyphosate, plus hopefully one other pesticide

##evaluating statistical significance using cochran-armitage test (Neurhauser et al. 1999)
#try fitting with 95 hour mortality un-rounded

#let's put together the columns we need to do the CA test
#first modify so that the survival is adjusted for 96 hours across all studies
# effects$half_life<-(effects$Duration_h*log(2))/log(1/effects$Survival)#need to modify survival by duration, using an exponential growth curve
# effects$adj_sur_96<-round(1/(2^(96/effects$half_life)),3)
names(effects)
df_sig<-effects[,c(1,5,8,10,17)] #make sure this output matches
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
CA_increasing_trend<-lapply(ch0,function(x) CochranArmitageTest(x,'increasing')) #run CA test over matrix list; format needs to match format in ?Cochran-Armitage 'dose' ex
CA_increasing_trend
#10,9,5,1 fail to reject the null hypothesis for pyraclostrobin
#non fail to reject the null hypothesis for glyphosate 

#for pryaclostrobin
#1 is cusaac 2015 acris blanchardi
#5 is cusaac 2017a anaxyrus cognatus
#9 is cusaac 2015 anaxyrus woodhousii
#10 is cusaac 2017a anaxyrus woodhousii
#low sample sizes for two (5,10), other two not much response from exposure, also lower application rates
#remove 4 studies for a total of 25 points including weighted control
#effects<-effects[!(effects$Study =='Cusaac_2015'| effects$Study == 'Cusaac_2017a'), ] #pyraclostrobin


#group by study, species, and application rate
effects<-effects %>%
  group_by(Study) %>% #study
  group_by(Species) %>% #species
  group_by(Application_Rate)%>% #application rate
  mutate(SE = SD/sqrt(length(effects))) #mutate to add standard error to account for variability between studies
effects<-as.data.frame(effects)
effects<-na.omit(effects)
controls<-effects[effects$Application_Rate == 0, ]  #take out all controls
control<-effects[1,] # take out single row which will contain final control
#function for weighted mean
weight_mean<-function(data,vector,weight){
  round(with(data, sum(vector*weight)/sum(weight)),2)
}
control$Survival<-weight_mean(controls,controls$Survival, controls$N_Exp) #survival
control$N_Exp<-sum(controls$N_Exp)#use sum of the N used in dose for final N_exp
effects<-effects[effects$Application_Rate != 0, ]
effects<-rbind(effects,control)
#effects <-effects[order(effects$Study),]


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
survival_sim <- matrix(data=NA,nrow=nrow(effects),ncol=nsims)
colnames(survival_sim)[1:nsims]<-paste0("Sur",1:nsims,"")
for(i in 1:nrow(effects)){
  survival_sim[i,]<-rbinom(nsims, effects[i,17], effects[i,10]) #trials is number of observations for each experiment, prob is adjusted prob of survival for 96 hours
  survival_sim[i,]<-round(survival_sim[i,]/effects[i,17],3) #divide the number surviving by total trials to find the proportion that actually survived 
}
survival_sim<-as.data.frame(survival_sim)
#let's walk through this with Tom; 

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


#pyraclostrobin
#calculate dermal dosage 
# pmolweight<-387.8 #g/mol #comptox
# plogKow<- 4.44 #comptox
# kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
# hl<-4.91*24 #from comptox profile

#glyphosate
#calculate dermal dosage 
pmolweight<-169.1 #g/mol #comptox
plogKow<- -3.12 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.46*24 #from comptox profile

#calculate soil concentrations; application rate will remain the same, the exponent value changes with exposure parameter simulations
soil_concs_degs<-as.data.frame(log(2)/hl*96/exposure_sims$movement_rate_mean) #use 96 as duration, as we adjusted all survivals for a 96h timeframe
soil_concs<-((effects$Application_Rate*10)/16000)*1000 #1cm depth #converts it into units of ugg #this comes from Weir
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

BMDS_glyphosate_fin<-cbind(exp_n,mort_n,derm_dose) #name the output whatever run you'd like
BMDS_glyphosate_fin<-BMDS_glyphosate_fin[,c(1,6,2,4)]

#order by dose for clarity
BMDS_glyphosate_fin<-BMDS_glyphosate_fin %>%
  group_by(set) %>% #set of sims
  arrange(Dose, .by_group=T) #this will rearrange the output in terms of order of sims, but will still work fine


###only needed for BBMD, BMDS can do it in one go
#we need to split up the 1000 set simulation into 10 sets of 100
# bbmd<-split(BMDS_Run3, (as.numeric(rownames(BMDS_Run3))-1) %/% 3900)
# lapply(bbmd,dim) #check dimensions of the output; should be 3900*4 for each set of 10 (0-9)

#quick test if needed
# test<-bbmd
# testy<-test[[3]] #check that it's formatted with 100 sets of unique sims
# unique(testy$set) #should be 100

# colnames<-c("Dataset Index", "Dose", "N","Effect")
# bbmd<-lapply(bbmd, setNames, colnames)
##for BBMD format
# for(i in names(bbmd)){
#   write.csv(bbmd[[i]], paste0(i,"bbmdrun.csv"), row.names=FALSE)
# }

save.image(file='myEnvironment_simulation.RData') 
write.csv(BMDS_glyphosate_fin,'data_out/BMDS_glyphosate_fin.csv') #save whole file for bmds analysis

#scraps ###---- 
#tests trend in proportions? in this example, there does seem to be a trend in proportions? hmm
smokers  <- c( 83, 90, 129, 70 )
patients <- c( 86, 93, 136, 82 )
prop.test(smokers, patients)
prop.trend.test(smokers, patients)
prop.trend.test(smokers, patients, c(0,0,0,1))
?prop.trend.test
#survival and study? 
# we want to test for statistical significance between individual studies
#do fairly simple lm between models? #ancova?
dose <- as.data.frame(matrix(c(10,9,10,7, 0,1,0,3), byrow=TRUE, nrow=2, dimnames=list(resp=0:1, dose=0:3)))
dose
CochranArmitageTest(dose, "increasing")
CochranArmitageTest(dose)
CochranArmitageTest(dose, "decreasing")
?CochranArmitageTest
# not exactly the same as in package coin:
# independence_test(tumor ~ dose, data = lungtumor, teststat = "quad")
lungtumor <- data.frame(dose = rep(c(0, 1, 2), c(40, 50, 48)),
                        tumor = c(rep(c(0, 1), c(38, 2)),
                                  rep(c(0, 1), c(43, 7)),
                                  rep(c(0, 1), c(33, 15))))


tab <- table(lungtumor$dose, lungtumor$tumor)
CochranArmitageTest(tab, 'increasing')
str(tab)

out_1<-t(by_s_stat[[1]])
out_1
colnames(out_1)<-out[1,]
out<-out[2:3,]
row.names(out)<-c(0,1)
out<-as.matrix(out)
out
dose

CochranArmitageTest(out,'increasing')