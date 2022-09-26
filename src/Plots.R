#create plots for Paulukonis et al. 2021
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

set.seed(6379)


setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')
#load('myEnvironment_analysis.RData')


##this script contains plot and data analysis for the BMDS MCMC outputs
#all plots here are part of the manuscript
#some of the script will be similar to the comparison plots found in 'Comparison_BBMD_BMDS.R'
glypho<-read.csv("data_out/glyphosate_final_DD.csv")
pyra<-read.csv("data_out/pyraclostrobin_final_DD.csv")

############################ Main Plots---- 
##Figure 1 is workflow


### Figure 2 - plot of posterior probability weights across 9 models 
sims<-read.csv('data_out/BMDS_headline_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T)
mod_list<-sort(c('hill','log-logistic','logistic','log-probit','weibull','qlinear','probit','multistage','gamma')) #order the models by name
mcmc_ma = lapply(by_s, function(y) ma_dichotomous_fit(y[,2],y[,4],y[,3], fit_type = "mcmc")) #apply ma function over list

post<-lapply(mcmc_ma, function (x) x['posterior_probs']) #pull out posterior probabilities
postw<-as.data.frame(unlist(post))
postw<-tibble::rownames_to_column(postw, "Simulation")
postw$Simulation = substr(postw$Simulation,1,nchar(postw$Simulation)-1) #remove extraneous characters
colnames(postw)[2]<-'PosteriorProbs'
mod_list_f<-rep(mod_list, times=1000)
postw$Model<-mod_list_f #add column for model


ggplot(postw, aes(x = PosteriorProbs, y = Model, group = Model, fill = Model)) +
  geom_density_ridges(stat = "binline", bins = 75, scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Model") +
  xlab("Posterior Probability Weight by Model - 1000 Simulations, Pyraclostrobin")+
  theme(legend.position="none",  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16))
 


sims<-read.csv('data_out/BMDS_glyphosate_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T)
mod_list<-sort(c('hill','log-logistic','logistic','log-probit','weibull','qlinear','probit','multistage','gamma')) #order the models by name
mcmc_ma = lapply(by_s, function(y) ma_dichotomous_fit(y[,2],y[,4],y[,3], fit_type = "mcmc")) #apply ma function over list

post<-lapply(mcmc_ma, function (x) x['posterior_probs']) #pull out posterior probabilities
postw<-as.data.frame(unlist(post))
postw<-tibble::rownames_to_column(postw, "Simulation")
postw$Simulation = substr(postw$Simulation,1,nchar(postw$Simulation)-1) #remove extraneous characters
colnames(postw)[2]<-'PosteriorProbs'
mod_list_f<-rep(mod_list, times=1000)
postw$Model<-mod_list_f #add column for model


ggplot(postw, aes(x = PosteriorProbs, y = Model, group = Model, fill = Model)) +
  geom_density_ridges(stat = "binline", bins = 75, scale = 2) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Model") +
  xlab("Posterior Probability Weight by Model - 1000 Simulations, Glyphosate")+
  theme(legend.position="none",  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16))




#### Figure 3 done in MCMC_Pyra.R and MCMC_Glypho.R




####### Supplementary ----

#Figure S1 is in apprate_to_soilconc.R

#Figure S2
lab_bb<-read.csv('data_in/lab_bb.csv')#should be 1158
field_bb<-read.csv('data_in/field_bb.csv')  #should be 30611

#Figure S3


#Figure S4? BMDS or effect by body weight
ggplot(bmds_est, aes(y = Model)) +
  geom_density_ridges(aes(x=BMDSEstimates, fill = order), 
                      stat = "binline", bins = 100, scale = 0.9,
                      alpha = .6, ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.5,1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Fitted BMD with Upper (BMDU) and Lower (BMDL) Estimates",
    y = "Model",
    title = "Benchmark Dose (BMD) Estimates by Model" ) +
  theme_ridges()

#just for bmds
ggplot(bmds, aes(x=BMDSEstimates, group=Model)) + 
  geom_density(aes(fill=Model), alpha=0.2)+
  scale_fill_manual(values=c("#E69F00", "#D55E00", "#56B4E9", "#009E73"))+
  xlim(-2.5, 1)+
  xlab('BMDS')+
  ylab('Density')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 12))



mcmc_ma$exp_n1$
  plot(mcmc_ma$exp_n1)
plot(mcmc_q$exp_n1)
plot(mcmc_ma$exp_n1$Individual_Model_8)
derm_seq <-as.data.frame(seq(0,1.4, length=1000))  #sequence from 0 to 1.4 by x, 1000 length
colnames(derm_seq)[1]<-'dermdose'



