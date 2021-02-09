#create plots for Paulukonis et al. 2021
library(ggplot2)
library(ggridges)
library(readr)
library(bayestestR)
library(dplyr)
getwd()
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Paulukonis_Documents/manuscript_pop_modeling/amphibian_effects_model/amphibian_effects_model')

#plot of posterior probability
mydir = "C:/Users/eliza/Dropbox/amphib_modeling_manuscript/manuscript_pop_modeling/amphibian_effects_model/data_in/csv"
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


#plot of BMDs with upper and lower limit
test<-tidyr::gather(modsbmr, "mark","n",7:9)
test$mark <- factor(test$mark, levels = c("bmdl", "bmd", "bmdu"))
ggplot(test, aes(y = Model)) +
  geom_density_ridges(aes(x=n, fill = mark), 
                      stat = "binline", bins = 100, scale = 0.9,
                      alpha = .6, ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(limits = c(-1,5.0), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Fitted BMD with Upper (BMDU) and Lower (BMDL) Estimates",
    y = "Model",
    title = "Benchmark Dose (BMD) Estimates by Model" ) +
  theme_ridges()

#determine top model by median
#choose 3 models that are above threshold of 0.1 or greater
modsp = list.files(path=mydir, pattern="*models.csv", full.names=TRUE) #filter csv with posterior values
modsp <- plyr::ldply(modsp, read_csv)
modsbmr %>%
  group_by(Model) %>% 
  summarize(med=median(posterior_weight)) #returns median value of posterior weights across 1000 sets
mods3<-filter(modsp, grepl('LogProbit|QuantalLinear|DichotomousHill', model)) #pull out top 3 models

#max(derm_dose$Dose) #sequence from 0 to 1.4 by x, 1000 length
derm_seq <-as.data.frame(seq(0,1.4, length=1000))
colnames(derm_seq)[1]<-'dermdose'


