library(Rcpp)
library(RcppGSL)
library(RcppEigen)
library(ggplot2)
library(ggridges)
library(forcats)
library(plotly)
library(shiny)
library(scales)
library(ToxicR)
library(nlme)
library(ggridges)
library(dplyr)

##this script contains code to compare the BMDS ToxicR package for MCMC dose response analysis and the BBMD dose response analysis app

getwd()
set.seed(25395) #match seed to BBMD for comparison - this for non-zero mort
set.seed(31443) #match seed to BBMD for comparison - this for zero mort

##note; make sure that the data input matches either batch 1 (non-zero) or batch 2 (zero) mortality. 
#current iteration is batch 2, the zero mortality at zero dose

###BMDS ----
sims<-read.csv('data_in/BMDS_test_run_2.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T) #split by simulation
mod_list<-sort(c('hill','log-logistic','logistic','log-probit','weibull','qlinear','probit','multistage','gamma')) #order the models by name
mcmc_ma = lapply(by_s, function(y) ma_dichotomous_fit(y[,2],y[,4],y[,3], fit_type = "mcmc")) #apply ma function over list
post<-lapply(mcmc_ma, function (x) x['posterior_probs']) #pull out posterior probabilities 
postw<-as.data.frame(unlist(post))
postw<-tibble::rownames_to_column(postw, "Simulation")
postw$Simulation = substr(postw$Simulation,1,nchar(postw$Simulation)-1) #remove extraneous characters

mcmc_ma[1]
#plot BMDS posterior probs from the model average
colnames(postw)[2]<-'PosteriorProbs'
mod_list_f<-rep(mod_list, times=1000)
postw$Model<-mod_list_f #add column for model
#plot of all models
ggplot(postw, aes(x = PosteriorProbs, y = Model, group = Model, fill = Model)) +
  geom_density_ridges(stat = "binline", bins = 100, scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Model") +
  xlab("Posterior Distribution Weight of 1000 Simulated Studies")+
  theme_ridges()

#look at median posterior prob to determine top 3 BMDS models
tab_med<-postw %>%
  group_by(Model) %>% 
  summarize(med=median(PosteriorProbs))
write.csv(tab_med,'data_out/BMDS_median_zero.csv')
#for here, using median, it's log-probit, probit, and weibull

#pull out bmds and bmdu for model average
bmdsorder<-c('bmds','bmdl','bmdu')
bmds_order<-rep(bmdsorder, times=1000)
bmds<-lapply(mcmc_ma, function (x) x['bmd'])
bmds<-as.data.frame(unlist(bmds))
bmds<-tibble::rownames_to_column(bmds, "Simulation")
bmds$Simulation = substr(bmds$Simulation,1,nchar(bmds$Simulation)-1)
colnames(bmds)[2]<-'BMDSEstimates'
bmds$order<-bmds_order
bmds$Model<-'ModelAverage'

#run mcmc over top 3 models
mcmc_log<- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "log-probit",fit_type = "mcmc"))
mcmc_p <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "probit",fit_type = "mcmc"))
mcmc_w <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type = "weibull",fit_type = "mcmc"))

#organize and put outputs in dataframe
log<-lapply(mcmc_log, function (x) x['bmd'])
log<-as.data.frame(unlist(log))
log<-tibble::rownames_to_column(log, "Simulation")
log$Simulation = substr(log$Simulation,1,nchar(log$Simulation)-1)
colnames(log)[2]<-'BMDSEstimates'
log$order<-bmds_order
log$Model<-'Log-Probit'

p<-lapply(mcmc_p, function (x) x['bmd'])
p<-as.data.frame(unlist(p))
p<-tibble::rownames_to_column(p, "Simulation")
p$Simulation = substr(p$Simulation,1,nchar(p$Simulation)-1)
colnames(p)[2]<-'BMDSEstimates'
p$order<-bmds_order
p$Model<-'Probit'

w<-lapply(mcmc_p, function (x) x['bmd'])
w<-as.data.frame(unlist(w))
w<-tibble::rownames_to_column(w, "Simulation")
w$Simulation = substr(w$Simulation,1,nchar(w$Simulation)-1)
colnames(w)[2]<-'BMDSEstimates'
w$order<-bmds_order
w$Model<-'Weibull'

bmds_est<-rbind(bmds,log,p,w)
bmds_est_f<-bmds_est %>% filter(order == 'bmds') 

#BMDS KDE plot
ggplot(bmds_est_f, aes(x=BMDSEstimates, group=Model)) + 
  geom_density(aes(fill=Model), alpha=0.2)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9","#D55E00", "#009E73"))+
  xlim(-2.5, 1)+
  xlab('BMDS')+
  ylab('Density')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12))

#rm(list = ls()) #clear global if you'd like

###BBMD----
mydir = "C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/BBMD Outputs/batch1000/batch_set_2/csv"
modsbmr = list.files(path=mydir, pattern="*bmrs.csv", full.names=TRUE) #filter csv with bmds values
modsbmr<- plyr::ldply(modsbmr, read.csv)
modwse<-modsbmr[!(modsbmr$posterior_weight==1 ),]#remove any posterior weights of 1
colnames(modwse)[4]<-"Model"

modsp = list.files(path=mydir, pattern="*models.csv", full.names=TRUE) #filter csv with posterior values
modsp <- plyr::ldply(modsp, read.csv)

#plot BBMD posterior probs 
ggplot(modwse, aes(x = posterior_weight, y = Model, group = Model, fill = Model)) +
  geom_density_ridges(stat = "binline", bins = 100, scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Model") +
  xlab("Posterior Distribution Weight of 100 Simulated Studies")+
  theme_ridges()

#look at median posterior prob to determine top 3 BBMD models
tab_med<-modsbmr %>%
  group_by(model) %>% 
  summarize(med=median(posterior_weight)) #returns median value of posterior weights across 1000 sets
write.csv(tab_med,'data_out/BBMD_median_zero.csv')

bbmd_est<-filter(modsbmr, grepl('Model average|LogProbit', model)) #pull out top 3 models

#BBMD KDE plot
ggplot(bbmd_est, aes(x=bmd, group=model)) + 
  geom_density(aes(fill=model), alpha=0.3)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  xlim(0, 1.5)+
  xlab('BBMD')+
  ylab('Density')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12))
