#create plots for Paulukonis et al. 2021
library(ggplot2)
library(ggridges)
library(readr)
library(bayestestR)
library(dplyr)
getwd()
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Paulukonis_Documents/manuscript_pop_modeling/amphibian_effects_model/amphibian_effects_model')

##this script contains plot and data analysis for the BMDS MCMC outputs
#all plots here are part of the manuscript
#some of the script will be similar to the comparison plots found in 'Comparison_BBMD_BMDS.R'

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


