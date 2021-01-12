##script for simulation and analysis of dermal dose and effects
library(ggplot2)
library(dplyr)
library(tidyr)
setwd('C:\\Users\\eliza\\Dropbox\\amphib_modeling_manuscript\\manuscript_pop_modeling\\amphibian_effects_model')
load('myEnvironment.RData')

effects<-read.csv('data_in/Headline_updated.csv')
param<-read.csv('data_in/parameters.csv')
nsims<-1000

#Body weight ----
#desired output is 100 simulations of distributions of bw for every combo of study*species*dose
#recall that we need to convert SD to SE
#we'll assume a normal distribution here
#calculate SE
#split by study*species(dose)
effects<-effects %>%
  group_by(Study) %>% #study
   group_by(Species) %>% #species
  group_by(Application_Rate)%>% #application rate
  mutate(SE = SD/sqrt(length(effects)))
effects<-as.data.frame(effects)

#it's probably best if I take out all rows with zero except for one at this point
#39 total rows
#while we are treating simulations of BW within each group of study/species/application, 
#we will be treating each simulated batch of dose-response as its own set of individual records
#because we have accounted for variability in the simulated BWs

zero_r<-effects[1,] #pull out single 0 dose; will function as our threshold
effects<-effects[effects$Application_Rate != 0, ]  
effects<-rbind(effects,zero_r)
effects <-effects[order(effects$Study),]

#we need 100 simulations of the X different bws, but we'll simulate through each study*species*dose
#because of this, our BW simulation is slightly more complicated - this will account for variation between 
#astudies as well as species
by_s<-split(effects, list(effects$Study,effects$Species,effects$Application_Rate), drop=T)
bw_sim<-replicate(nsims,
                  {normalize <- lapply(by_s, function(y) rnorm(nrow(y), mean=y[[11]], sd=y[[19]]))
                  bw_sim<-as.data.frame(unlist(normalize))
                  #bw_sim<-tibble::rownames_to_column(bw_sim, "row_names")
                  })
bw_sim<-do.call(cbind.data.frame, bw_sim)
#bw_sim<-bw_sim[,c(1,seq(2, ncol(bw_sim), by=2))] 
colnames(bw_sim)[1:nsims]<-names<-paste0("BW",1:nsims,"")
#colnames(bw_sim)[1]<-'Study'

#Survival ---- 
#need to modify survival by duration, using an exponential growth curve
#then use a binomial distribution to estimate probability of survival
effects$half_life<-(effects$Duration_h*log(2))/log(1/effects$Survival)
effects$adj_sur_96<-round(1/(2^(96/effects$half_life)),3)
#effectsx<-effects # for other dermal calc

survival_sim <- matrix(data=NA,nrow=nrow(effects),ncol=nsims)
colnames(survival_sim)[1:nsims]<-paste0("Sur",1:nsims,"")
for(i in 1:nrow(effects)){
  survival_sim[i,]<-rbinom(nsims, effects[i,17], effects[i,21]) #trials is number of observations for each experiment, prob is adjusted prob of survival for 96 hours
survival_sim[i,]<-round(survival_sim[i,]/effects[i,17],3) #divide the number surviving by total trials to find the proportion that actually survived 
}
survival_sim<-as.data.frame(survival_sim)


#Exposure Parameters----
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
 

#Calculation of dermal dose----
#calculate one time values for exposure - these do not change for Headline or Headline AMP
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
my_list[458] #check the list outputs
derm<-as.data.frame(t(do.call(rbind.data.frame, my_list)))
colnames(derm)<-paste0("dermdose",1:nsims,"")
row.names(derm)<-NULL
#outputs are in ugg


# format the dermal dose estimates and survival estimates to use in the BBMD app
# should be in order with effects, so that each experiment N matches with each study
mort_n<-as.data.frame(sapply(survival_sim, function(x) effects$N_Exp - round(x*effects$N_Exp,0))) #number of individuals dying per exp
#we are treating each batch like its own seperate study; essentially, each set of doses has a corresponding set of survival data
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

BBMD_Run1<-cbind(exp_n,mort_n,derm_dose)
BBMD_Run1<-BBMD_Run1[,c(1,6,2,4)]

#BBMD requires each ID group to be in descending order
BBMD_Run1<-BBMD_Run1 %>%
  group_by(set) %>% #set of sims
  arrange(Dose, .by_group=T) #this will rearrange the output in terms of order of sims, but will still work fine
#we need to split up the 1000 set simulation into 10 sets of 100
bbmd<-split(BBMD_Run1, (as.numeric(rownames(BBMD_Run1))-1) %/% 3900)
lapply(bbmd,dim)
test<-bbmd
testy<-test[[3]] #check that it's formatted with 100 sets of unique sims
unique(testy$set) #should be 100
colnames<-c("Dataset Index", "Dose", "N","Effect")
bbmd<-lapply(bbmd, setNames, colnames)
for(i in names(bbmd)){
  write.csv(bbmd[[i]], paste0(i,"bbmdrun.csv"), row.names=FALSE)
}


save.image(file='myEnvironment.RData')


#Survial prior estimates and improving DR model - WIP ----
slope<-runif(nsims, min=-5, max=0) #changing it to -5
intercept<-runif(nsims, min=0.5,max=1.5)

sl<-slope[1]
int<-intercept[1]
test<-int+sl*derm[,1] #prior model, survival; using simple lm equation for now


survival_simprior <- matrix(data=NA,nrow=47,ncol=nsims)
colnames(survival_simprior)[1:nsims]<-paste0("Sur2",1:nsims,"")
#sur_sum_abs_diff <- matrix(data=NA, nrow=100, ncol=1)

for (j in 1:ncol(derm)){
  # for (i in 1:nrow(sur_sum_abs_diff)){
  survival_simprior[,j]<-intercept[j] + slope[j]*derm[,j]
  sur_sum_abs_diff<-sum(abs(survival_simprior[,j] - survival_sim[,j])) 
}
#fix this; but not at the moment
# abs(0.2-0.8)
# sum_abs_differences <- sum(abs(survival_simprior[,1] - survival_sim[,1])) 
# sur_sum_abs_diff <- matrix(data=NA, nrow=nsims, ncol=1)

#median of that 
#example:
#find those that are greater are 900, fit slope and intercept to those
#simulate again
#only accept those that have a sum of abso diff less tha 900
#second time, you keep going until you have 100 that beat your old score
#sort by score and take the top 25%
#higher scores mean greater differences between observed and predicted
#lower scores mean lower differences between observed and predicted



#Caluclate original dermal dose for comparison, and output new dataset for BBMD test----
pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
dt<-param[1,2]
move_rate<-param[3,2]
hl<-4.91*24 #from comptox profile
bioavail<-param[5,2]
dsan<-param[9,2]*effectsx$M_Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_deg<-log(2)/hl*96/move_rate
soil_concsn<-((effectsx$Application_Rate*10)/16000)*1000 #1cm 

#this calculates dermal dose 
dermal_dose<-(soil_concsn^soil_concs_deg * kp_pyra * (dsan/dt) * derm_frac * bioavail)/effectsx$M_Body_Weight_g
head(dermal_dose) #take a look
effectsx$dermaldose<-dermal_dose #ugg
fm1<-lm(adj_sur_96 ~ dermaldose, data = effectsx)
summary(fm1)
plot(adj_sur_96 ~ dermaldose, data=effectsx)

p<-ggplot(effectsx, aes(x=dermaldose, y=adj_sur_96, group=Study)) +
  geom_point(aes(color=Study))+
  ylab("Adjusted 96hr Survival") +
  xlab("Dermal Dose (ug/g)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key=element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12),
        )
p


#this calculates application rate v. study
p2<-ggplot(effectsx, aes(x=Application_Rate, y=adj_sur_96, group=Study)) +
  geom_point(aes(color=Study, shape = Study), size=3.5)+
  scale_shape_manual(values=c(0:4,15:17))+
  ggtitle("Mean 96hr Survival as a function of Application Rate by Study") +
  ylab("96hr Survival Rate") +
  xlab("Application Rate (ug a.i./cm2)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key=element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size = 13)
  )
p2



#prep data for BBMD
#effects$N_s<-round(effects$adj_sur_96*effects$N_Exp,0)
#write.csv(effects,'doseresponse_data.csv')
#go ahead and eliminate duplicate dermal dose rows
#length(unique(dermal_dose))
#length(unique(effects$adj_sur_96))
#length(dermal_dose)
#issue for BBMD platform; althought we have a unique set of dermal dose and survival values, some dermal doses
#are the same due to body weight and application rate (exp: soil; )
# dr<-read.csv('doseresponse_data.csv')
# length(unique(dr$dermaldose)) #pull out only unique values for now, just to look
# out<-group_by(dr, dermaldose) %>% slice(1)
# write.csv(out,'doseresponse_data_fix.csv')

#Process description notes ----
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
#100 dose response curves
#slope and intercept are two initial parameters
#generate uniform priors 
#slope 0 to -10 (dose - x axis)
#Intercept 0.5 to 1.5 (survival -  y axis)
#I believe I can just generate the intercept and slope, and then run them through DRC to generate priors
#what I think I'm doing is simply generating a set of slope/intercept values 
#and creating a new set of prior distributions where I estimate survival based on those values
#we then compare that estimated survival to the survival outputs that I generated earlier...?

#at heart, what we're doing is estimating new survival outputs 
#given a new slope/intercept in the dose response curve, using the dermal dose data that we generate
#we're then comparing that to our original survival outputs
#dumb question; why not just use the original model? is that because we don't have enough data points?
#I think it's because we don't have enough data points
#yi =  survival
#xi =  dose
#b0 =  y intercept
#b1 =  slope

#yi = b0 + -b1*x1# linear model equation
#we could go ahead and generate these 100 slopes and intercepts and then take a look at the sum of absolute deviations

#we would then use this equation:
#dermal_sum_abs_differences[generation,i] <- sum(abs(dermal_dose_amphib[,i] - amphib_concs))
#subtracting the estimating survival from the actual survival estimates that we generated
#Sum of absolute differences attempts to help use identofy whether the function we have used n

#surival = intercept + -slope x dose linear model with negative slope
#47 pairs of those observations
#predict simple linear model 
#sum mean deviation or sum of squared deviations
#run a bunch of different models at the same time!
#model selection aspect

#for each set of those we are generating slopes and intercepts
#sum of absolute deviations
#47 doses - for that model you would fit that
#model prediction
#dose response predicts 
#47 deviations 

#when dose = 0, we would expect survival to be at least 0.5
#once I output the lms; what's the next step?
#I've found a package that allows me to estimate dose response curves 
#100 dose-response cloud
#go ahead
#2parameters in regression model
#priors on slope and intercept
#slope is around 1
#slope flat to negative
#squared deviations or as in abc model, use sum of absolute 
#intercept that's close to 1
#slope that's negative 
#intercept uniform from 0.5-1.5
#slope uniform 0--10
#draw priors from these 100
#ok so we're generating 100 simulations of the uniform 
#using the drc package
#then we compare these to our predicted observations 
#then we'll take the inputs 

#we should have 100 distinct simulations of 47 dermal dose estimates to compare to our 100 survival estimates
#let's loop through and calculate some lms for this set, and look at the drc package
# 
# dermd<-tidyr::gather(derm,"simd","vald",1:100)
# surs<-tidyr::gather(survival_sim,"sims","vals",1:100)
# fin<-cbind(dermd,surs)
# by_sim<-split(fin, fin$simd)
# out<-purrr::map(by_sim, ~lm(vals ~ vald, data = .x) )
# summary(out[[67]]) # take a look at results
#wait, was this even necessary?


#further down the line we will use a covariate approach
#assumption of indepedence when that isn't necessarily true

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
dsa<-param[9,2]*effectsx$M_Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_deg<-log(2)/hl*96/move_rate #with adjusted time of 96h to match adjusted survival
soil_concs<-((effectsx$Application_Rate*10)/16000)*1000 #1cm 

#this calculates dermal dose 
dermal_dose<-(soil_concs^soil_concs_deg * kp_pyra * (dsa/dt) * derm_frac * bioavail)/effectsx$M_Body_Weight_g
head(dermal_dose) #take a look
effectsx$dermaldose<-dermal_dose


fm1<-lm(adj_sur_96 ~ dermaldose, data = effects)
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

# lmp <- function (modelobject) {
#   f <- summary(modelobject)$fstatistic
#   p <- pf(f[1],f[2],f[3],lower.tail=F)
#   attributes(p) <- NULL
#   return(p)
# }
# summ_sims<-sapply(out,lmp) # pulls out p values for comparison 
# min(summ_sims)
# max(summ_sims)
# mean(summ_sims)
#compare to original data
# fm1<-lm(adj_sur_96 ~ dermaldose, data = effects)
# summary(fm1)
#   
# mean(effects$M_Body_Weight_g) #this seems pretty similar
# mean(bw_sim$BW7) #pull out random BW
# 
# mean(effects$adj_sur_96) # this isn't as far off
# mean(surs$vals) #recall that this is 4700 variables
# 
# mean(effects$dermaldose) #this does not seem so similar; so it seems to be an issue with the calculation
# mean(dermd$vald) #note that does give mean of 4700 vs 47 values, but trend appears to be there for 



#let's try comparing a few outputs
# dsa<-exposure_sims[1,5]*bw_sim[,2]^exposure_sims[1,5]
# dermal_dosen<-(soil_concs^soil_concs_degs[1,] * kp_pyra * (dsa/exposure_sims[1,1]) * exposure_sims[1,6] * exposure_sims[1,3])/bw_sim[,2]
# dsa
# dsan
# #the issue was the dsa

#effects <-effects[order(effects$Study),] #order by study to match output of split if not already ordered

##
#update.packages(ask=FALSE)
#library(fitdistrplus)
#library(drc)
#library(purrr)

#calculate dsa 
# dsa<-list()
# for (j in 2:ncol(bw_sim)){
#   for (i in 1:nrow(exposure_sims)){
#     dsa[[j]]<-exposure_sims[i,5]*bw_sim[,j]^exposure_sims[i,5]
#   }}
# dsa[[3]]
# 
# my_list<-list()
# lapply(2:ncol(bw_sim), function(j) 
#   lapply(1:nrow(exposure_sims), function(i) 
#     lapply(1:nrow(soil_concs_degs), function(t) 
#       my_list[[j]]<-(soil_concs^soil_concs_degs[t,] * kp_pyra * ((exposure_sims[i,5]*bw_sim[,j]^exposure_sims[i,5])/exposure_sims[i,1]) * exposure_sims[i,6] * exposure_sims[i,3])/bw_sim[,j]
#     )))
# my_list<-list()
# for (j in 2:ncol(bw_sim)){
#   for (i in 1:nrow(exposure_sims)){
#     for (t in 1:nrow(soil_concs_degs)){
#       dsa<-exposure_sims[i,5]*bw_sim[,j]^exposure_sims[i,5]
#       dermal_dose<-(soil_concs^soil_concs_degs[t,] * kp_pyra * (dsa/exposure_sims[i,1]) * exposure_sims[i,6] * exposure_sims[i,3])/bw_sim[,j]
#       my_list[[j]]<-dermal_dose
#     }}}

