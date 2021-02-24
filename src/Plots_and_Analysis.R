#create plots for Paulukonis et al. 2021
library(ggplot2)
library(ggridges)
library(readr)
library(bayestestR)
library(dplyr)
library(tidyverse)
#setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Paulukonis_Documents/manuscript_pop_modeling/amphibian_effects_model/amphibian_effects_model')
load('myEnvironment_analysis.RData')

##this script contains plot and data analysis for the BMDS MCMC outputs
#all plots here are part of the manuscript
#some of the script will be similar to the comparison plots found in 'Comparison_BBMD_BMDS.R'

####Analysis----\
sims<-read.csv('data_out/BMDS_run3.csv') #read in the compiled simulation data from 'Data_Simulation.R'
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
write.csv(postw, 'data_out/PosteriorWeights.csv')

#pull out 3 models to look at bmds
tab_med<-postw %>%
  group_by(Model) %>% 
  summarize(med=median(PosteriorProbs))
#probit
#log-logistic
#qlinear
write.csv(tab_med,'data_out/BMDS_posteriorweight_median.csv')

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


#run mcmc over top 3 models
mcmc_ll<- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-logistic",fit_type = "mcmc"))
mcmc_p <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "probit",fit_type = "mcmc"))
mcmc_q <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "qlinear",fit_type = "mcmc"))

#organize and put outputs in dataframe
ll<-lapply(mcmc_ll, function (x) x['bmd'])
ll<-as.data.frame(unlist(ll))
ll<-tibble::rownames_to_column(ll, "Simulation")
ll$Simulation = substr(ll$Simulation,1,nchar(ll$Simulation)-1)
colnames(ll)[2]<-'BMDSEstimates'
ll$order<-bmds_order
ll$Model<-'Log-Logistic'

p<-lapply(mcmc_p, function (x) x['bmd'])
p<-as.data.frame(unlist(p))
p<-tibble::rownames_to_column(p, "Simulation")
p$Simulation = substr(p$Simulation,1,nchar(p$Simulation)-1)
colnames(p)[2]<-'BMDSEstimates'
p$order<-bmds_order
p$Model<-'Probit'

q<-lapply(mcmc_p, function (x) x['bmd'])
q<-as.data.frame(unlist(q))
q<-tibble::rownames_to_column(q, "Simulation")
q$Simulation = substr(q$Simulation,1,nchar(q$Simulation)-1)
colnames(q)[2]<-'BMDSEstimates'
q$order<-bmds_order
q$Model<-'QuantalLinear'

bmds_est<-rbind(bmds,ll,p,q)
write.csv(bmds_est,'data_out/BMDS_bmdsbmdlbmdu_top3.csv')

save.image(file='myEnvironment_analysis.RData') 


####Plots----
# Figure 1: application rate and 96hr mortality figure?
# Figure 2: Distribution of posterior probabilities across 9 models to show why we picked them
# Figure 3: KDEs of BMDs across top 3 models
# Figure 4: DR curve for top 3 models and 95% credibility interval + median, plotted with original dataset



#Figure 2 - plot of posterior probability - BMDS
mydir = "C:/Users/epauluko/Dropbox/amphib_modeling_manuscript/manuscript_pop_modeling_OLD/amphibian_effects_model/data_in/csv"
modsbmr = list.files(path=mydir, pattern="*bmrs.csv", full.names=TRUE)
modsbmr<- plyr::ldply(modsbmr, read_csv)
modsbmr<-modsbmr[!(modsbmr$posterior_weight==1 ),]
colnames(modsbmr)[4]<-"Model"

ggplot(modsbmr, aes(x = posterior_weight, y = Model, group = Model, fill = Model)) +
  geom_density_ridges(stat = "binline", bins = 100, scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Model") +
  xlab("Posterior Probability by Model for 1000 Simulated Studies")+
  theme_ridges()

#choose 3 models that are above threshold of 0.1 or greater
modsp = list.files(path=mydir, pattern="*models.csv", full.names=TRUE) #filter csv with posterior values
modsp <- plyr::ldply(modsp, read_csv)
modsbmr %>%
  group_by(Model) %>% 
  summarize(med=median(posterior_weight)) #returns median value of posterior weights across 1000 sets
mods3<-filter(modsp, grepl('LogProbit|QuantalLinear|DichotomousHill', model)) #pull out top 3 models


#Figure 3 - KDE of benchmark doses 
#all models 

#top 3 
ggplot(bmds_est, aes(x=BMDSEstimates, group=Model)) + 
  geom_density(aes(fill=Model), alpha=0.2)+
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#56B4E9", "#009E73"))+
  xlim(-2.5, 1)+
  xlab('BMDS')+
  ylab('Density')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12))



#Figure 4 - plot confidence interval and median of 3 models, overlay original dataset on top
#simulate some doses
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
