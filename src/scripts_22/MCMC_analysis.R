### April 2022 Analysis

###Evaluating MCMC MA runs, and pulling out 'top' model using BMD

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
#get dataset 
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')

sims<-read.csv('data_out/BMDS_headline_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T) #split by simulation



###note, only uncomment if need to revisit MA
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


#pull out parameters for all 1000 models
#then, we'll manually do a log-logistic DR curve on them with the parameters
#we'll look for the BMDS and across the 1000, we'll get the high and low estimates of the curve (using quantile)
#we'll plot the DR curves with BMDs corresponding to the median, 2.5, and 97.5% 
#the quantiles are used as the upper/lower bounds on the curve



#model the log-logistic single model using the simulated sets
ll_fit <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type="log-logistic",fit_type="mcmc"))


#pull model average BMD, bmdl, bmdu
bmdsorder<-c('bmds','bmdl','bmdu')
bmds_order<-rep(bmdsorder, times=1000)
bmds<-lapply(ll_fit, function (x) x['bmd']) #pull out bmds and bmdls
bmds<-as.data.frame(unlist(bmds))
bmds<-tibble::rownames_to_column(bmds, "Simulation")
bmds$Simulation = substr(bmds$Simulation,1,nchar(bmds$Simulation)-1)
colnames(bmds)[2]<-'BMDSEstimates'
bmds$order<-bmds_order
bmds$Model<-'log_logistic'


#note; quantile does not return the specific row values, so I used a 'closest' function to pick the datasets
quantile_sims<-quantile(bmds_bmd$BMDSEstimates, probs = c(0.025,0.5,0.975), names=T)
print(quantile_sims)

five_per<-bmds_bmd[which.min(abs(1.842985 -bmds_bmd$BMDSEstimates)),]
fifty_per<-bmds_bmd[which.min(abs(2.520026-bmds_bmd$BMDSEstimates)),]
ninetyfive_per<-bmds_bmd[which.min(abs(2.995469 -bmds_bmd$BMDSEstimates)),]
perc<-rbind(five_per,fifty_per,ninetyfive_per)
print(perc) #these are the data rows to pull to use the data to plot the curves


#get parameters 
para<-lapply(ll_fit, function (x) x['parameters'])
para<-as.data.frame(unlist(para))
para<-tibble::rownames_to_column(para, "Exp")
colnames(para)[2]<-'Value'
para_list_f<-rep(c("p1","p2","p3"), times=1000)
para$Parameters<-para_list_f

##just test with a basic set of parameters for now

##this would not run when I input the parameters themselves; llogistic or LL.3 don't like my fixed parameters 
alpha=para[1,2]
beta=para[2,2]
gamma=para[3,2]

g=gamma
a=alpha
b=beta

dr<-by_s[[1]]
colnames(dr)[c(2,4)]<-c('dose','response')
dr$effect<-dr$response/dr$N

##none of this would work for me (except the basic version where we fit the data)
test_ll <- drm(response/N~dose, data=dr, fct = LL.3u(), type='binomial')
plot(test_ll, type='all')
out<-LL.3u(upper=1, fixed = c(para[1,2], para[2,2], para[3,2]), names = c("a", "b", "g"))
test_ll <- drm(response/N~dose, data=dr, fct = LL.3u(upper=1, fixed=c(alpha,beta,gamma)), type='binomial', )
plot(test_ll)
test_ll <- drm(response/N~dose, set, weights = N, data = dr, 
               fct = LL.3(fixed = c(alpha, beta, gamma), names = c("a", "b", "g")), type="binomial",
              )



#here's the equation for the actual log logistic model from the BMDS User's Manual
log_logistic<-function(x){
(g + ((1-g)/1+exp(-a-b*log(x))) )

}
p(dose) = log_logistic(dose) #this is the predicted probability distribution of the dose based on the fitted model



#here's an example of just running the basic logistic regression using GLM
model_list<-list()
for(model in 1:3){
  model=1
  dr = by_s[[model]]
  dr$response<-dr$Effect/dr$N
    
  log_dr<-glm(response ~ Dose, family=binomial(link='logit'), data=dr )
  model_list[[model]]<-log_dr
  
}


log_dr<-glm(response ~ Dose, family=binomial(link='logit'), data=dr )
dr2<-by_s[[4]]
dr2$response<-dr2$Effect/dr2$N

pred_dr<-predict(log_dr, newdata=dr2, se.fit=TRUE, type="link")
dr$p <- plogis(pred_dr$fit)
dr$lower <- plogis(pred_dr$fit - 1.96*pred_dr$se.fit)
dr$upper <- plogis(pred_dr$fit + 1.96*pred_dr$se.fit)
plot(response ~ Dose, data=dr2, type="l", ylim=c(0,1),xlab="Dose", ylab="Response")
points(response ~ Dose, dr2)
lines(response ~ Dose, dr, lty=2)
lines(response ~ Dose, dr, lty=2)

####scrap notes:

##just plot the fitted models (using some function for logistic regression), and stack as 95/5% CI
##we have the parameters, just need to put them in a function

#check range assumptions in this implementation of the drm function
#maybe log transformed?

#fitting and not predict function 

#(gamma+(1-gamma))/(1+exp{-[alpha+betaln(x)]}) #actual ll equation
#notes:
#see BMDS guidance

#not sure why gammma is wonky; beta and alpha look ok? (gamma between 0 and 1, but isn't)

##might need some guidance, but let me figure it out with the log-logistic model
##pulling top 2.5, 97.5 and 50% quantile BMDS values as upper and lower
##challenge: haven't plotted from scratch before, probably not as complicated as I'm making it

##other challenge: the guidance online is confusing

