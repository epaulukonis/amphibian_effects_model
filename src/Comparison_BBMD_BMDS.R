library(scales)
library(ToxicR)
library(nlme)
library(ggridges)
library(dplyr)
getwd()
set.seed(25395) #match seed to BBMD for comparison

#BBMDS using ToxicR - upload csv of runs ----
#current iteration is the non-zero mortality at zero dose
sims<-read.csv('data_in/BMDS_test_run.csv')
by_s<-split(sims, list(sims$set), drop=T) #split by simulation
mod_list<-sort(c('hill','log-logistic','logistic','log-probit','weibull','qlinear','probit','multistage','gamma')) #order the models by name
mcmc_ma = lapply(by_s, function(y) ma_dichotomous_fit(y[,2],y[,4],y[,3], fit_type = "mcmc")) #apply ma function over list
post<-lapply(mcmc_ma, function (x) x['posterior_probs']) #pull out posterior probabilities 
postw<-as.data.frame(unlist(post))
postw<-tibble::rownames_to_column(postw, "Simulation")
postw$Simulation = substr(postw$Simulation,1,nchar(postw$Simulation)-1)
colnames(postw)[2]<-'PosteriorProbs'
mod_list_f<-rep(mod_list, times=1000)
postw$Model<-mod_list_f #add column for model


mydir = "C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/BBMD Outputs/batch1000/batch_set_2/csv"
modsbmr = list.files(path=mydir, pattern="*bmrs.csv", full.names=TRUE)
modsbmr<- plyr::ldply(modsbmr, read_csv)
modsbmr<-modsbmr[!(modsbmr$posterior_weight==1 ),]
colnames(modsbmr)[4]<-"Model"