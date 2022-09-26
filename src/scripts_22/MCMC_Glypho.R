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
setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')

#simulated datasets
sims<-read.csv('data_out/BMDS_glyphosate_fin_32122.csv') #read in the compiled simulation data from 'Data_Simulation.R'
by_s<-split(sims, list(sims$set), drop=T) #split by simulation



##og dataset
effects<-effects<-read.csv('data_in/Glyphosate_updated.csv') #run code in Data_Sim to clean up the effects and remove rows we don't want
mod<-read.csv('data_in/As_Modifier.csv')
param<-read.csv('data_in/parameters.csv')

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
#  summarize(med=median(PosteriorProbs))
#  print(tab_med)
#write.csv(tab_med,'data_out/BMDS_posteriorweight_median_headline.csv')
 
 #gamma 

##here we model the ll model (top)
#model the log-logistic single model using the simulated sets
gam_fit <- lapply(by_s, function(y) single_dichotomous_fit(y[,2],y[,4],y[,3],model_type="gamma",fit_type="laplace"))


# #pull model average BMD, bmdl, bmdu
bmdsorder<-c('bmds','bmdl','bmdu')
bmds_order<-rep(bmdsorder, times=1000)
bmds<-lapply(gam_fit, function (x) x['bmd']) #pull out bmds and bmdls
bmds<-as.data.frame(unlist(bmds))
bmds<-tibble::rownames_to_column(bmds, "Simulation")
bmds$Simulation = substr(bmds$Simulation,1,nchar(bmds$Simulation)-1)
colnames(bmds)[2]<-'BMDSEstimates'
bmds$order<-bmds_order
bmds$Model<-'gamma'


## fit og glyphosate data ----
#the duration for all of the glyphosate studies is already set to 96 hours; not need to calculate it using exponential growth curve
effects$adj_sur_96<-effects$Survival
print(names(effects))
df_sig<-effects[,c(1,5,9,11,18)] #make sure this output matches
df_sig<-na.omit(df_sig)
df_sig$mort<-round(df_sig$N_Exp - (df_sig$Survival*df_sig$N_Exp),0)
df_sig$alive<-df_sig$N_Exp-df_sig$mort
head(df_sig)
#split data into list by study/species for cochran armitage trend test
by_s_stat<-split(df_sig, list(df_sig$Study, df_sig$Species), drop=T)
by_s_stat<-lapply(by_s_stat, function (x) x[c('Application_Rate','mort','alive')]) 

#change format to matrix to run CA trend test to test for trends in increasing significance
ch0=list()
for (i in seq_along(by_s_stat)) {
  a=t(by_s_stat[[i]])
  colnames(a)<-a[1,]
  a<-a[c(3,2),] #need to have no mort first row, mort second
  row.names(a)<-c(0,1)
  a<-as.matrix(a)
  ch0[[i]]<-a
}

#run CA test over matrix list; format needs to match format in ?Cochran-Armitage 'dose' ex
CA_increasing_trend<-lapply(ch0,function(x) CochranArmitageTest(x,'increasing')) 
CA_increasing_trend
#none fail to reject the null hypothesis for glyphosate 


#group by study, species, and application rate
effects<-effects %>%
  group_by(Study) %>% #study
  group_by(Species) %>% #species
  group_by(Application_Rate)%>% #application rate
  mutate(SE = SD/sqrt(length(effects))) #mutate to add standard error to account for variability between studies
effects<-as.data.frame(effects)
#effects<-na.omit(effects)
controls<-effects[effects$Application_Rate == 0, ]  #take out all controls
control<-effects[1,] # take out single row which will contain final control
#function for weighted mean
weight_mean<-function(data,vector,weight){
  round(with(data, sum(vector*weight)/sum(weight)),2)
}
control$Survival<-weight_mean(controls,controls$Survival, controls$N_Exp) #survival
control$N_Exp<-sum(controls$N_Exp)#use sum of the N used in dose for final N_exp
effects<-effects[!effects$Application_Rate == 0, ] #remove application rate of 0
effects<-rbind(effects,control)
effects$adj_sur_96<-effects$Survival
effects<-na.omit(effects)
effects <-effects[order(effects$Study),]


fam<-as.data.frame(cbind(effects$Family, effects$Study, effects$Application_Rate))
names(fam)<-c('Family', 'Study', "Application")
mod_As<- merge(fam,mod, by  = "Family") 
mod_As<-mod_As[order(mod_As$Study),]


pmolweight<-169.1 #g/mol #comptox
plogKow<- -3.12 #comptox
kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
hl<-4.46*24 #from comptox profile
dt<-param[1,2]
move_rate<-param[3,2]
bioavail<-param[5,2]
dsa<-param[9,2]*effects$Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_deg<-log(2)/hl*96/move_rate #with adjusted time of 96h to match adjusted survival
soil_concs<-((effects$Application_Rate*10)/16000)*1000 #1cm mixing depth; or change depending on soil?
dermal_dose<-(soil_concs^soil_concs_deg * kp_pyra * (dsa/dt) * derm_frac * bioavail)/effects$Body_Weight_g
head(dermal_dose) #take a look
effects$dermaldose<-dermal_dose
effects$Mortality<-1-effects$adj_sur_96

bw_effect<-effects$Body_Weight_g
exp<-bw_effect^mod_As$Exponent 
as_d<-mod_As$Modifier*exp


# 50% SA
derm_d<-effects$dermaldose 
app_d<-as.numeric(effects$Application_Rate)
derm_original<-derm_d*((as_d*app_d)/2) #modify the calculated dermal dose by the SA and application rate product
TCR_org_50<-as.data.frame(derm_original/effects$dermaldose) #this calculates the ratio between the direct (overspray) and indirect (soil)
names(TCR_org_50)<-"TCR"
TCR_org_50$Family<-effects$Family #add family
TCR_org_50<-na.omit(TCR_org_50)
TCR_org_50_r<- TCR_org_50 %>% 
  group_by(Family) %>%
  summarise(TCR_r = mean(TCR))
print(TCR_org_50_r$TCR_r)
print(TCR_org_50_r$Family)

effects$dermaldose<-derm_original

effects<-effects[with(effects, order(Study, Species, Application_Rate)), ]
effects$death<-round(effects$Mortality*effects$N_Exp,0)
#effects<-effects[order(effects$dermaldose),]
effects$Exp<-ifelse(effects$Method=="Overspray",1,0)
print(effects$death)
by_s[[4]]$Effect


write.csv(effects,'data_out/glyphosate_final_DD.csv')



og_data<-matrix(0,nrow=nrow(effects), ncol=5)
colnames(og_data) <- c("Dose","N","Incidence", "Exp","BW")
og_data[,1] <- effects$dermaldose
og_data[,2] <- effects$N_Exp
og_data[,3] <- effects$death
og_data[,4] <- effects$Exp
og_data[,5] <- effects$Body_Weight_g

round(1-(og_data[,3]/og_data[,2]),3)
print(sims[1:39,4])
print(effects$death)

max(sims[,2])


#### time to get the upper and lower bounds out of the 1000 curves ----
#post-process all data for credible intervals
d<-og_data[,1]
nd<-length(by_s[[1]]$Dose)
nsims <- 1000
mortality_df <- matrix(ncol = nsims, nrow = nd)


# para<-lapply(gam_fit, function (x) x['parameters'])
# para<-as.data.frame(unlist(para))
# para<-tibble::rownames_to_column(para, "Exp")
# colnames(para)[2]<-'Value'
# para_list_f<-rep(c("p1","p2","p3"), times=1000)
# para$Parameters<-para_list_f
# 
# b<-para[para$Parameters=='p3',]
# max(b$Value)
# mean(b$Value)


for(i in 1:length(by_s)){
  parms <- gam_fit[[i]]$parameters
  g <-  1/(1+exp(-parms[1])); 
  a <- parms[2];
  b <- parms[3]; 
  mortality_df[,i] <- g + (1-g)*pgamma(b*d,a,1)
} #

#14 23
# 25


percentiles_df <- matrix(ncol = 3, nrow = nd)
for(i in 1:nd){
  mort_percentiles <- quantile(mortality_df[i,], probs = c(0.025,0.50, 0.975), names=T)
  percentiles_df[i,1] <- mort_percentiles[[1]]
  # percentiles_df[i,2] <- mort_percentiles[[2]]
  # percentiles_df[i,3] <- mort_percentiles[[3]]
  percentiles_df[i,2] <- mort_percentiles[[2]]
  # percentiles_df[i,5] <- mort_percentiles[[5]]
  # percentiles_df[i,6] <- mort_percentiles[[6]]
  percentiles_df[i,3] <- mort_percentiles[[3]]
}
dim(percentiles_df)

df_percentiles <- data.frame(x=d, val= as.vector(percentiles_df), 
                             variable=rep(paste0("category", 1:3), each=39))

#### plot curve with original data as well as upper and lower bounds (2.5% and 97.5% percentiles of mortality), with OG points ----

gam_laplace<-single_dichotomous_fit(og_data[,1],og_data[,3],og_data[,2],model_type="gamma",fit_type="laplace")

parms <- gam_laplace$parameters
g <-  1/(1+exp(-parms[1])); 
a <- parms[2];
b <- parms[3]; 
rval_og <- g + (1-g)*pgamma(b*d,a,1)
df<-as.data.frame(cbind(rval_og,d))
names(df)<-c("effect","dose")



lerror<-df_percentiles[df_percentiles$variable == "category1",2]
median<-df_percentiles[df_percentiles$variable == "category2",2]
uerror<-df_percentiles[df_percentiles$variable == "category3",2]

med_data<-as.data.frame(cbind(median,df$dose))

og_data<-as.data.frame(og_data)
og_data[,5]<-as.character(og_data[,5])


my_breaks_m<-c(-10,-5,-3,-2.7)

main<-ggplot() + 
  #geom_density(data=bmds_df, aes(x=bmds))+
  geom_line(data = df, aes(x=log(dose), y=effect))+
  geom_ribbon(aes(x = log(df$dose), ymin =lerror, ymax = uerror), alpha = .2) +
  geom_point(data=og_data, aes(x=log(Dose),y=(Incidence/N), colour=Exp, fill=Exp))+
  ggtitle("Gamma Glyphosate Curve") +
  ylab("Mortality") +
  xlab("log(Dose (ug/g))")+
  scale_x_continuous(limits=c(-10,-2.7),breaks=my_breaks_m,expand = c(0, 0)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16), legend.position = 'none') #first run without this to get legend

main

ld50<-read.csv('data_in/ld50_glypho.csv')
#ld50<-ld50[!ld50$LD50>100,]
cons<-sort(ld50$LD50)[1:5]
ld50c<-ld50[ld50$LD50 %in% cons,]
frog<-c(0.002980385, 0.000480939)

min(cons)
max(cons)
mean(cons)

frog<-ld50[ld50$LD50 %in% frog,]

my_breaks<-c(-10,-5,-3,-2.7)

bmds_hist_gly<-
  ggplot(data=bmds,aes(x=log(BMDSEstimates), y=..count..,group=order, fill=order))+
  geom_histogram(bins=100)+
  ylab("Density BMDs") +
  xlab("log(Dose (ug/g))")+
  scale_x_continuous(limits=c(-10,-2.7), breaks=round(my_breaks,2), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        axis.title.x = element_text(size=14, face='bold'),
        axis.title.y = element_text(size=14, face='bold'),
        plot.title = element_text(face = 'bold', size = 16), legend.position="none")

bmds_hist_gly




final_glypho<-main +
  geom_boxplot(data=ld50c, aes(x=log(LD50), y=Mortality), width=0.05) +
  geom_vline(xintercept=log(frog[1,1]), linetype='dotted', col = 'red')+
  geom_vline(xintercept=log(frog[2,1]), linetype='dotted', col = 'red')+
  # geom_point(data=frog, aes(x=LD50,y=Mortality),colour="black", shape=18, size=4) +
  geom_point(aes(x=log(0.017),y=0.50),colour="red") 


grid.arrange(final_glypho, bmds_hist_gly,nrow=2,heights=c(4.5,1.5))




# ggplot(data = df, aes(x=dose, y=effect)) + 
#   geom_line()+
#   geom_line(data=med_data,aes(x=V2,y=median,colour='blue'))+
#   geom_ribbon(aes(ymin =lerror, ymax = uerror), alpha = .2) +
#   geom_point(data=og_data, aes(x=Dose,y=(Incidence/N), colour=BW, fill=BW, shape=BW))+
#   scale_shape_manual(values=c(1:7))+
#   
#   ggtitle("Gamma model fit, with 2.5 and 97.5 percentiles as upper and lower bound: Glyphosate") +
#   ylab("Mortality") +
#   xlab("Estimated Dermal Dose (ug/g)")+
#   #scale_x_continuous(expand = c(0, 0)) + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"), 
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size= 12, face='bold'),
#         axis.text.y = element_text(size=12, face='bold'),
#         axis.title.x = element_text(size=14, face='bold'),
#         axis.title.y = element_text(size=14, face='bold'),
#         plot.title = element_text(face = 'bold', size = 16)) #first run without this to get legend

#look at max of points (larger x-axis for simulations could explain )