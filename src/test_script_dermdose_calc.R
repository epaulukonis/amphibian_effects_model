##script for simulation and analysis of dermal dose and effects 
library(ggplot2)
library(dplyr)
library(tibble)
library(fitdistrplus)


effects<-read.csv('data_in/Headline_updated.csv')
param<-read.csv('data_in/parameters.csv')


#Body weight ----
#desired output is 100 simulations of distributions of bw for 5 different studies
#recall that we need to convert SD to SE
#we'll assume a normal distribution here
nsims<-100
#calculate SE
effects_SE<-effects %>%
  group_by(Study) %>%
  mutate(SE = SD/sqrt(length(effects))) 
levels(effects_SE$Study)
#we need 100 simulations of the 45 different bws, but we'll simulate through each study
#because of this, our BW simulation is slightly more unique, because mean BWs come from each study
by_s<-split(effects_SE, effects_SE$Study) 
bw_sim<-replicate(nsims,
{normalize <- lapply(by_s, function(y) rnorm(nrow(y), mean=y[[11]], sd=y[[19]]))
  bw_sim<-as.data.frame(unlist(normalize))
  bw_sim<-tibble::rownames_to_column(bw_sim, "row_names")
                })
bw_sim<-do.call(cbind.data.frame, bw_sim)
bw_sim<-bw_sim[,c(1,seq(2, ncol(bw_sim), by=2))] 
colnames(bw_sim)[2:100]<-names<-paste0("BW",1:100,"")
colnames(bw_sim)[1]<-'Study'

#Survival ---- 
#need to modify survival by duration, using an exponential growth curve
#then use a binomial distribution to estimate probability of survival
effects$half_life<-(effects$Duration_h*log(2))/log(1/effects$Survival)
effects$adj_sur_96<-1/(2^(96/effects$half_life))

survival_sim <- matrix(data=NA,nrow=47,ncol=100)
colnames(survival_sim)[1:100]<-paste0("Sur",1:100,"")
for(i in 1:nrow(effects)){
  survival_sim[i,]<-rbinom(nsims, effects[i,17], effects[i,20]) #trials is number of observations for each experiment, prob is adjusted prob of survival for 96 hours
  
}
survival_sim<-as.data.frame(round((survival_sim/30), 3))


#Exposure Parameters----
exposure_sims <- matrix(data=NA,nrow=100,ncol=6)
colnames(exposure_sims) <- c("dt_mean","movement_rate_mean","bioavail_mean",
                                  "dermal_sa_slope_mean","dermal_sa_exponent_mean",
                                  "dermal_fraction_mean")
exposure_sims[,1]<-rnorm(nsims, mean=param[1,2], sd=param[2,2])
exposure_sims[,2]<-rnorm(nsims, mean=param[3,2], sd=param[4,2])
exposure_sims[,3]<-rnorm(nsims, mean=param[5,2], sd=param[6,2])
exposure_sims[,4]<-rnorm(nsims, mean=param[7,2], sd=param[8,2])
exposure_sims[,5]<-rnorm(nsims, mean=param[9,2], sd=param[10,2])
exposure_sims[,6]<-rnorm(nsims, mean=param[11,2], sd=param[12,2])
exposure_sims<-as.data.frame(exposure_sims)



#Calculation of dermal dose and comparison to survival simulations----
#calculate one time values for exposure - these do not change for Headline
pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.91*24 #from comptox profile

#calculate soil concentrations; application rate will remain the same, the exponent value changes with exposure parameter simulations
soil_concs_deg<-as.data.frame(log(2)/hl*96/exposure_sims$movement_rate_mean) #use 96 as duration, as we adjusted all survivals for a 96h timeframe
soil_concs<-((effects$Application_Rate*10)/16000)*1000 #1cm 

my_list<-list()
for (j in 2:ncol(bw_sim)){
  for (i in 1:nrow(exposure_sims)){
    for (t in 1:nrow(soil_concs_deg)){
dsa<-exposure_sims[i,1]*bw_sim[,j]^exposure_sims[i,5]
dermal_dose<-(soil_concs^soil_concs_deg[t,] * kp_pyra * (dsa/exposure_sims[i,1]) * exposure_sims[i,6] * exposure_sims[i,3])/bw_sim[,j]
my_list[[j]]<-dermal_dose
}}}
my_list[57] #check the list outputs, first one will be empty
my_list<-my_list[-1]
names(my_list)<-paste0("dermdose",1:100,"")
test<-melt(my_list)
#try using gather to put the final columns together

#we should have 100 distinct simulations of 47 dermal dose estimates to compare to our 100 survival estimates

##Process description notes ----
#our goal: to fill in the curve for these parameters by simulating the information we don't have
#we have several steps here to do that

#things we could simulate, even if we have no individual body data
#body weights and survival and exposure parameters - means
#we'll simulate using the mean and sd

#fit a curve to the numbers we have so far
#dermal dose vs. survival -  that's the end goal here, where we're plotting those
#100 simulations to start
#normal distribution for all corresponding distributions, minus binomial for survival (0-1)
#first choice to vary is exposure variables 
#second choice vary means of bw (from some n) so by group mean (meaning 7 curves)
#assign sampling distribution - by group
#third choice is to vary survival as an output (binomial distribution)

#replacing body weights
#replacing survival
#new column for dose
#think about putting it in a data frame
#100 snap shots of response curve
#SO we're simulating all 3 of these aspects 100 times, for each study
#and then we're gonig to be modifying the dose-reponse curve using the parameters 

#ok SO we have to modify the exposure variables; we have to use sampling distributions 
#because those are for individuals; SO modify the sds for both BWs and exposure 
#because we need to account for group variability, as the SD and exposure variables are for individual variability
#simulated survival 
#simulated body weights
#and simulated exposure
#I need to ask him how to do this

#we DO need to adjust for duration
#align survival metrics to be at 96 hours
#add in adjustment for ones above/below 
#so we have to fit an exponential decay model to those ones in order to sync everything to a 96hr mortality
#build e-lambda t 

#official imputation data method
#go ahead and use the means of adults/juveniles of other species to fill in for cusaac 2015/bruhl
#sd use the coeff of variant approach, and approach - NAH
#JK; so jsut use the means/sd of similar species at terrestrial stage

#dose response thing 

#at this point we're not doing a Bayesian analysis, we're just doing a monte carlo analysis
#Khan -  workshop for SRA; we still need to determine how to fit those curves


#first; we'll impute the data, and record how/why we did that
#we'll start with body weights - ask Tom about normal distribution; what did he mean by 'first round' is uniform?
#then we'll try survival, but we'll have to adjust the duration rate following an exponential curve
#then we'll do exposure parameters 

#to ask tom: 
#1) you mentioned running a uniform distribution for the'first' round; can you clarify? 
#if we have mean and sd for these, why would we do a uniform (assuming normal/binomial)?
#2)for the adjusting of the SDs; I know you mentioned it, but I didn't catch all of it. We need to account for individual variability
#we obviously don't have individual BWs to calculate population standard deviation, can you elaborate on a good approach?


#answer: 1) you do NOT need to do a uniform; do rnorm or rbinom (not quite sure how to do this yet)
#answer: 2) use standard error to get group variability, because we don't have individual variability

#recall that we have to develop these simulations for each study
#for BW/survival
#100 sims x 5 studies x the number of individual means per batch?
#for exposure parameters
#100 sims x 6
#ok, so I need to modify the script to run 100 x 45

#what we'd then do is I believe is plot estimated dermal dose with new body weights 

#Tom wants to know what the package does
#for httk; what sort of data exists

#some sort of math function -  
#x axis o-some exposure
# y is 0-1
#at 0 our survival is 1
#dose response curve packages?
#robust estimate of dose response estimate that we collected
#could do a couple different generations 

#100 sets of simulations
#each time the dose response curve
#also simulating a monte carlo curve through
#top 50 parameters; slope parameter the first time was between 0 and -5
#probability 
#amphibians 100 times, draw from binomial
#10 amphibians, probability is 0.8 

#rbinom(100, # trials are n(row), prob is the survival for that row)
#probability 
#ok, so now we're saying that many actually survived 28/30

#comparison 
#creating 1000 dose reponse curves
#and then generating those points

#covential fitting tehcniques
#least  squares or AIC 
#Initital Calculation of Dermal Dose and Comparison to Survival----
##for now let's assume 1cm is our best bet for even mixing, assuming a 1cm soil depth
#we're treating headline and headline amp the same, per results from Cusaac 2015 and Cusaac 2017
#plot mortality as a function of body burden
#bb =  cSoil*e^(-kt)*kp*(SA/dt)Saf*DAF*h/BW
#SA= 10*bw^0.667
#Kp
#dt *
#Saf *
#BW
#H (bioavail) *
#DAF *
#cSoil is a conversion from application rate to estimate soil concentration

##calculate parameters, add in parameters from ABC
#for Kp: I need mol-weight and logKow info of Pyraclostrobin and Metconazole

#we'll be modeling added toxicity via the Monte Carlo route, and maybe the ABC route
#for now, I will treat them the same, using the molecular info for pyra as the primary source of data

#pyra
pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox

#metconazole
# mmolweight<-319.8 #g/mol #comptox
# mlogKow<- 3.77 #comptox

kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
# kp_met =  10^(-2.72+(0.71*mlogKow)-(0.0061*mmolweight))

dt<-param[1,2]
move_rate<-param[3,2]
hl<-4.91*24 #from comptox profile
bioavail<-param[5,2]
dsa<-param[9,2]*effects$M_Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_degradation<-log(2)/hl*effects$Duration_h/move_rate
soil_concs<-((effects$Application_Rate*10)/16000)*1000 #1cm 

#this calculates dermal dose 
dermal_dose<-(soil_concs^soil_concs_degradation * kp_pyra * (dsa/dt) * derm_frac * bioavail)/effects$M_Body_Weight_g
head(dermal_dose) #take a look
effects$dermaldose<-dermal_dose

fm1<-lm(Survival ~ dermaldose, data = effects)
summary(fm1)

plot(Survival ~ dermaldose, data=effects, type="p",
     xlab="Dermal Dose (ug/g)", ylab="Survival (%)")  

plot(Application_Rate ~ dermaldose, data=effects, type="p",
     xlab="Derm Dose", ylab="App Rate")  

ggplot(effects, aes(x=Application_Rate, y=Survival, group=Study)) +
  geom_point(aes(color=Study))

#Scraps ####---- 
#scrap code
# norm<-function(x,y){
#   rnorm(nsims, mean=x,sd=y)
# }


# create function
# sampletest <- function(df){
#   df %>%
#     group_by(Study) %>%
#     do({
#       norm = rnorm(nsims, df$M_Body_Weight_g, df$SE)
#       data.frame(norm)
#     }) 
# }
# normalize <- apply(effects, 1, function(x) rnorm(nsims, mean=x[11], sd=x[19]))

#bw is 11
#se is 19

# by_s[[5]][[11]]
# 
# 
# ?do.call
# normalize <- apply(dataFrameApply, 1, function(x) rnorm(x[1], mean=x[2], sd=x[3]))
# out<-do.call(rbind.data.frame, lapply(by_s,function(x){
#   rnorm(nrow(x[[]]), mean=x[[11]], sd=x[[19]])
# }))
# 
# 
# 
# nrow(by_s[[4]])
# 
# 
# normalize <- lapply(by_s, function(y) rnorm(nrow(y), mean=y[[11]], sd=y[[19]]))
# df <- data.frame(matrix(unlist(normalize), nrow=length(normalize), ncol=1))
# testy<-as.data.frame(unlist(normalize))
# testy<-tibble::rownames_to_column(testy, "row_names")
# 
# test[[1]]
# test[[2]]
# normalize <- lapply(by_s, function(y) rnorm(nrow(y), mean=y[[11]], sd=y[[19]]))
# names(normalize)
# 
# out<-do.call(rbind.data.frame, lapply(by_s,function(x){
#   rnorm(nsims, effects_SE$M_Body_Weight_g, effects_SE$SE)
#   
# }))
# 
# normalize <- apply(dataFrameApply, 1, function(x) rnorm(x[1], mean=x[2], sd=x[3]))
# require(plyr)
# func <- function(df){
#   return(data.frame(rnorm(nrow(.), effects_SE$M_Body_Weight_g, effects_SE$SE)))
# }
# test<-ddply(effects_SE, .(Study), func)
# 
# 
# #ddply(data.frame, variable(s), function, optional arguments)
# 
# #100*45 = 4500 of bw
# #100*45 = 4500 of survival
# #100*6 = 600 parameters
# 
# BW_sim<-ddply(effects_SE, .(Study), func)
# colnames(BW_sim)[2]<-'simulated_BW'
##test that it works
# testy<-effects_SE[1:4,]
# set.seed(1)
# norm = rnorm(nsims, testy$M_Body_Weight_g, testy$SE)
# hist(norm)
# test_s <- fitdist(norm, "norm")
# 
# sim_pa<-function(x,output){
#   for(t in 1:nrow(x)){
#     out<-return(rnorm(nsims, mean=x[t,2], sd=x[t+1,2]))
#   }
# }
# out<-sim_pa(param)
# out
# 
# 
# 
# 
# # R function
# f = function(x, output) {
#   # x is the row of type Character
#   # access element in first column
#   name = x[1]
#   # access element in second column
#   income = x[3]
#   #your code to process x
#   cat(name, income, "\n")
# }
# 
# 
# #apply(X, MARGIN, FUN, …)
# apply(celebrities, 1, f) 

# func <- function(df){
#   return(data.frame(rnorm(nrow(.), effects_SE$M_Body_Weight_g, effects_SE$SE)))
# }
