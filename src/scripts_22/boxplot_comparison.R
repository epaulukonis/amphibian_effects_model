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




#### Glyphosate ----

#script for comparing LD50 DDs as boxpltos across Web-ICE, SSD estimates, and our OG data 1000 sims
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')

#simulated datasets
sims<-read.csv('data_out/BMDS_glyphosate_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T) #split by simulation


#Web-ICE raw data 
ld50<-read.csv('data_in/ld50_glypho.csv')


gam_fit <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type="gamma",fit_type="laplace"))

#### Get the LD50 estiamtes for all simulated datasets;
#first, estimate the mortalities associated with each of the datasets and their corresponding doses
#what I chose to do was the pull the value closest to 0.50 mortality for each dataset, and the dose associated with that

nd<-39
nsims <- 1000
mortality_df <- as.data.frame(matrix(ncol = nsims, nrow = nd))
mortality_df_ld50<-as.data.frame(matrix(ncol = nsims, nrow = 1))
dose_df<-as.data.frame(matrix(ncol = nsims, nrow = nd))
dose_df_ld50<-as.data.frame(matrix(ncol = nsims, nrow = 1))
dose_df_quant<-as.data.frame(matrix(ncol = nsims, nrow = 1))

for(i in 1:length(by_s)){
  parms <- gam_fit[[i]]$parameters
  g <-  1/(1+exp(-parms[1])); 
  a <- parms[2];
  b <- parms[3];
  d<-by_s[[i]]$Dose
  mortality_df[,i] <- g + (1-g)*pgamma(b*d,a,1)
  dose_df[,i]<-d
  mortality_df_ld50[,i]<-which(abs(mortality_df[,i] - 0.50) == min(abs(mortality_df[,i]  - 0.50)))
  dose_df_ld50[,i]<-dose_df[mortality_df_ld50[,i],]

  dose_df_quant[,i]<-quantile(dose_df[,i],probs = c(0.50), names=T)
  
} 

#from the original DF: our LD50 was 0.017
toxicR<-as.data.frame(t(dose_df_ld50))
toxicR$Name<-"ToxicR"
names(toxicR)[1]<-"Dose"


#now let's get the original LD50s from web-ICE
names(ld50)<-c("Dose","Name")
ld50$Name<-"Web-ICE"


#to get Alex's data, you'll need to run the x_amphibian.rmd to get the full set of 1000 simulated lines and pull out the
#doses associated with the HC5
hc5dat<-bootdat %>% 
  group_by(variable) %>%
  arrange(abs(value - 0.05)) %>%
  slice(1)
  
#HC50 (50 percent of species impacted)
hc50dat<-bootdat %>% 
  group_by(variable) %>%
  arrange(abs(value - 0.50)) %>%
  slice(1)

hc50dat<-hc50dat[,1:2]
names(hc50dat)<-c("Dose","Name")
hc50dat$Name<-"SSD"
hc50dat$Dose<-hc50dat$Dose* (0.05*(10^-3.4)) # convet to LC50* BCF


ld50_df<-rbind(toxicR,ld50, hc50dat)
ggplot(data=ld50_df,aes(x=Name, y=log(Dose),group=Name, fill=Name))+
  geom_boxplot()



#use HC5s
hc5dat<-hc5dat[,1:2]
names(hc5dat)<-c("Dose","Name")
hc5dat$Name<-"SSD"
hc5dat$Dose<-hc5dat$Dose* (0.05*(10^-3.4)) # convet to LC50* BCF

ld50_df<-rbind(toxicR,ld50, hc5dat)
ld50_df$Name<-factor(ld50_df$Name, levels = c("SSD","Web-ICE","ToxicR"))

ggplot(data=ld50_df,aes(x=Name, y=log(Dose),group=Name, fill=Name))+
  geom_boxplot()

#OG prediction SSD (using Alex's model estimate)
ld50_alex<-as.data.frame(Gly.Pred.HC.Data$HC.Value)
names(ld50_alex)<-"Dose"

ld50_alex$Name<-"SSD"
ld50_alex$Dose<-ld50_alex$Dose* (0.05*(10^-3.4)) # convet to LC50* BCF

ld50_df<-rbind(toxicR,ld50, ld50_alex)
ld50_df$Name<-factor(ld50_df$Name, levels = c("SSD","Web-ICE","ToxicR"))
ggplot(data=ld50_df,aes(x=Name, y=log(Dose),group=Name, fill=Name))+
  geom_boxplot()+
  labs(y = expression(paste('Log Glyphosate LD50 [ ug ', g^-1, ' ]')), 
       x = 'Method',
       title = paste0("Comparison of LD50 estimates across 3 Methods")) +
  theme_minimal()



#### Pyraclostrobin ----
#script for comparing LD50 DDs as boxpltos across Web-ICE, SSD estimates, and our OG data 1000 sims
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')

#simulated datasets
sims<-read.csv('data_out/BMDS_headline_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T) #split by simulation

#Web-ICE raw data 
ld50<-read.csv('data_in/ld50_pyra.csv')

ll_fit <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type="log-logistic",fit_type="laplace"))

nd<-29
nsims <- 1000
mortality_df <- as.data.frame(matrix(ncol = nsims, nrow = nd))
mortality_df_ld50<-as.data.frame(matrix(ncol = nsims, nrow = 1))
dose_df<-as.data.frame(matrix(ncol = nsims, nrow = nd))
dose_df_ld50<-as.data.frame(matrix(ncol = nsims, nrow = 1))
dose_df_quant<-as.data.frame(matrix(ncol = nsims, nrow = 1))


for(i in 1:length(by_s)){
  parms <- ll_fit[[i]]$parameters
  g <- 1/(1+exp(-parms[1])); 
  a <- parms[2];
  b <- parms[3]; 
  d<-by_s[[i]]$Dose
  mortality_df[,i] <- g + (1-g)*(1/(1+exp(-a-b*log(d)))) #or by d or by_s[[i]]$Dose
  dose_df[,i]<-d
  mortality_df_ld50[,i]<-which(abs(mortality_df[,i] - 0.50) == min(abs(mortality_df[,i]  - 0.50)))
  dose_df_ld50[,i]<-dose_df[mortality_df_ld50[,i],]
  dose_df_quant[,i]<-quantile(dose_df[,i],probs = c(0.50), names=T)
  
} 

#from the original DF: our LD50 was 0.017
toxicR<-as.data.frame(t(dose_df_ld50))
toxicR$Name<-"ToxicR"
names(toxicR)[1]<-"Dose"


#now let's get the original LD50s from web-ICE
names(ld50)<-c("Dose","Name")
ld50$Name<-"Web-ICE"


#to get Alex's data, you'll need to run the x_amphibian.rmd to get the full set of 1000 simulated lines and pull out the
#doses associated with the HC5
hc5dat<-bootdat %>% 
  group_by(variable) %>%
  arrange(abs(value - 0.05)) %>%
  slice(1)

#HC50 (50 percent of species impacted)
hc50dat<-bootdat %>% 
  group_by(variable) %>%
  arrange(abs(value - 0.50)) %>%
  slice(1)

hc50dat<-hc50dat[,1:2]
names(hc50dat)<-c("Dose","Name")
hc50dat$Name<-"SSD"
hc50dat$Dose<-hc50dat$Dose*  (0.05*(10^3.99)) # convet to LC50* BCF


ld50_df<-rbind(toxicR,ld50, hc50dat)
ggplot(data=ld50_df,aes(x=Name, y=log(Dose),group=Name, fill=Name))+
  geom_boxplot()

#HC5s
hc5dat<-hc5dat[,1:2]
names(hc5dat)<-c("Dose","Name")
hc5dat$Name<-"SSD"
hc5dat$Dose<-hc5dat$Dose*  (0.05*(10^3.99))# convet to LC50* BCF

ld50_df<-rbind(toxicR,ld50, hc5dat)
ld50_df$Name<-factor(ld50_df$Name, levels = c("SSD","Web-ICE","ToxicR"))

ggplot(data=ld50_df,aes(x=Name, y=log(Dose),group=Name, fill=Name))+
  geom_boxplot()

#OG prediction SSD (using Alex's model estimate)
ld50_alex<-as.data.frame(pyra.pred.hc.data$HC.Value)
names(ld50_alex)<-"Dose"

ld50_alex$Name<-"SSD"
ld50_alex$Dose<-ld50_alex$Dose* (0.05*(10^3.99)) # convet to LC50* BCF

ld50_df<-rbind(toxicR,ld50, ld50_alex)
ld50_df$Name<-factor(ld50_df$Name, levels = c("SSD","Web-ICE","ToxicR"))

pyraplot<-ggplot(data=ld50_df,aes(x=Name, y=log(Dose),group=Name, fill=Name))+
  geom_boxplot()+
  labs(y = expression(paste('Log Pyraclostrobin LD50 [ ug ', g^-1, ' ]')), 
       x = 'Method',
       title = paste0("Comparison of LD50 estimates across 3 Methods")) +
  theme_minimal()

