### April 2022 Analysis

###Evalauting MCMC MA runs, and pulling out 'top' model using LD50 or BMD


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
library(drc)

set.seed(6379)
#get dataset 
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')

sims<-read.csv('data_out/BMDS_headline_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T) #split by simulation
mod_list<-sort(c('hill','log-logistic','logistic','log-probit','weibull','qlinear','probit','multistage','gamma')) #order the models by name
mcmc_ma = lapply(by_s, function(y) ma_dichotomous_fit(y[,2],y[,4],y[,3], fit_type = "mcmc")) #apply ma function over list

post<-lapply(mcmc_ma, function (x) x['posterior_probs']) #pull out posterior probabilities 
postw<-as.data.frame(unlist(post))
postw<-tibble::rownames_to_column(postw, "Simulation")
postw$Simulation = substr(postw$Simulation,1,nchar(postw$Simulation)-1) #remove extraneous characters
colnames(postw)[2]<-'PosteriorProbs'
mod_list_f<-rep(mod_list, times=1000)
postw$Model<-mod_list_f #add column for model

#look at median of PW across all models 
tab_med<-postw %>%
  group_by(Model) %>% 
  summarize(med=median(PosteriorProbs))
print(tab_med)
#write.csv(tab_med,'data_out/BMDS_posteriorweight_median_headline.csv')


#pull out parameters for all 1000 models
#then, we'll manually do a log-logistic DR curve on them with the parameters
#we'll look for the LD50 and across the 1000, we'll get the high and low estimates of the curve
#we'll plot the median, 2.5, and 2.75 CI curves



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


bmds_bmd<-bmds[bmds$order == 'bmds',]
quantile_sims<-as.data.frame(quantile(bmds_bmd$BMDSEstimates, probs = c(0.05,0.5,0.95), names=F))


para<-lapply(ll_fit, function (x) x['parameters'])
para<-as.data.frame(unlist(para))
para<-tibble::rownames_to_column(para, "Exp")
colnames(para)[2]<-'Value'
para_list_f<-rep(c("p1","p2","p3"), times=1000)
para$Parameters<-para_list_f


dr<-by_s[[1]]
colnames(dr)[c(2,4)]<-c('dose','response')

test_ll <- drm(response/N~dose, data=dr,, fct = LL.3u(), type='binomial')
plot(test_ll, type='all')
LL.3u(upper=1, fixed = c(para[1,2], para[2,2], para[3,2]), names = c("b", "d", "e"))


#issue: is the drm function dichotomous ?

ll_fit$exp_n1$options

#user defined
alpha=para[1,2]
beta=para[2,2]
gamma=para[3,2]
g=gamma


log(g/(1-g))

test_ll <- drm(response/N~dose, set, weights = N, data = dr, fct = LL.3(fixed = c(alpha, beta, gamma), names = c("b", "d", "e")), type="binomial")
plot(test_ll)

#check range assumptions in this implementation of the drm function
#maybe log transformed?

#fitting and not predict function 


#(gamma+(1-gamma))/(1+exp{-[alpha+betaln(x)]}) #actual ll equation
#notes:
#see BMDS guidance

#not sure why gammma is wonky; beta and alpha look ok? (gamma between 0 and 1, but isn't)


