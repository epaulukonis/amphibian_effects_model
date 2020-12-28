#create plots for Paulukonis et al. 2021
library(ggplot2)
library(ggridges)
#plot a) histogram of the 8 posterior weights using ggridges 
#plot b) histogram of the 8 BMD with CI (BMDU/DL) for 
#plot c) DR curve using 1000 generated doses with parameters for top 4 models, plotted with original data


##plot a)
modws<-read.csv('C:/Users/eliza/Dropbox/amphib_modeling_manuscript/manuscript_pop_modeling/bbmd_analysis/bbmd1_modelweights.csv')
modwse<-modws[!(modws$posterior_weight==1 ),]
colnames(modwse)[4]<-"Model"

ggplot(modwse, aes(x = posterior_weight, y = Model, group = Model, fill = Model)) +
  geom_density_ridges(stat = "binline", bins = 100, scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Model") +
  xlab("Posterior Distribution Weight of 100 Simulated Studies")+
  theme_ridges()


##plot b)
modbmds<-read.csv('C:/Users/eliza/Dropbox/amphib_modeling_manuscript/manuscript_pop_modeling/bbmd_analysis/bbmd1_modelparams.csv')
