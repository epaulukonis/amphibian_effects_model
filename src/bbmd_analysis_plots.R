#create plots for Paulukonis et al. 2021
library(ggplot2)
library(ggridges)
library(readr)
library(bayestestR)
library(dplyr)

#plot a) histogram of the 8 posterior weights using ggridges 
#plot b) histogram of the 8 BMD with CI (BMDU/DL) for 
#plot c) DR curve using 1000 generated doses with parameters for top 3 models, plotted with original data


##plot a)
setwd("C:/Users/eliza/Dropbox/amphib_modeling_manuscript/manuscript_pop_modeling/amphibian_effects_model/data_in")
mydir = "csv"
modsbmr = list.files(path=mydir, pattern="*bmrs.csv", full.names=TRUE)
modsbmr<- plyr::ldply(modsbmr, read_csv)
modsbmr<-modsbmr[!(modsbmr$posterior_weight==1 ),]
#colnames(modsbmr)[4]<-"Model"

ggplot(modsbmr, aes(x = posterior_weight, y = Model, group = Model, fill = Model)) +
  geom_density_ridges(stat = "binline", bins = 100, scale = 1.5) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  ylab("Model") +
  xlab("Posterior Distribution Weight of 1000 Simulated Studies")+
  theme_ridges()


## plot b)
# ggplot(modsbmr, aes(x = bmd, y = Model, group = Model, fill = Model)) +
#   geom_density_ridges(stat = "binline", bins = 100, scale = 1.5) +
#   scale_y_discrete(expand = c(0, 0)) +
#   scale_x_continuous(limits = c(0,1.25), breaks=c(0, 0.25, 0.50, 0.75, 1.0), expand = c(0, 0)) +
#   coord_cartesian(clip = "off") +
#   ylab("Model") +
#   xlab("Estimated Benchmark Dose of 1000 Simulated Studies")+
#   theme_ridges()
# 
# 
# ggplot(modsbmr, aes(x = bmdl, y = Model, group = Model, fill = Model)) +
#   geom_density_ridges(stat = "binline", bins = 100, scale = 1.5) +
#   scale_y_discrete(expand = c(0, 0)) +
#   scale_x_continuous(limits = c(0,1.25), breaks=c(0, 0.25, 0.50, 0.75, 1.0), expand = c(0, 0)) +
#   coord_cartesian(clip = "off") +
#   ylab("Model") +
#   xlab("Estimated Lower Benchmark Dose of 1000 Simulated Studies")+
#   theme_ridges()

colnames(modsbmr)
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


##plot c)
modsp = list.files(path=mydir, pattern="*models.csv", full.names=TRUE)
modsp <- plyr::ldply(modsp, read_csv)
#determine top model by median
modsbmr %>%
  group_by(model) %>% 
  summarize(med=median(posterior_weight))
#choose 3 models that are above threshold of 0.1 or greater

mods3<-filter(modsp, grepl('LogProbit|QuantalLinear|DichotomousHill', model))

#we'll loop through the 3 different models and populate the survival outputs for each of those
#we will sequence a range of values to represent dermal dose, and build 1000 curves
#get 89% credible interval
#median

max(derm_dose$Dose) #sequence from 0 to 1.4 by x, 1000 length
derm_seq <-as.data.frame(seq(0,1.4,length.out = 1000))
colnames(derm_seq)[1]<-'dermdose'

####QuantalLinear
QL<-mods3 %>%
  filter(model == 'QuantalLinear') %>%
  select(-c( model, seed, ppp, sd, p025, p50,p975,n_eff,rhat)) %>%
  tidyr:: spread(parameter, mean)
QL_s<-list()
for (i in 1:nrow(QL)){
  en<-exp(1)
  QL_s[i]<- QL[i,2] + (1-QL[i,2])*(1-en^(-QL[i,3]*derm_seq))
}
QL_s_f<-as.data.frame(t(do.call(rbind.data.frame, QL_s)))
colnames(QL_s_f)<-paste0("mortality",1:ncol(QL_s_f),"")
row.names(QL_s_f)<-NULL
#median curve over rows
QL_med_vals <- as.data.frame(apply(QL_s_f, 1, median, na.rm = T))
#CI curves over rows
QL_CI<-apply(as.matrix(QL_s_f), 1, ci, ci=0.95, method = "HDI")
CI_high<-matrix(data=NA, nrow=1000, ncol=1)
CI_low<-matrix(data=NA, nrow=1000, ncol=1)
for (i in 1:length(QL_CI)){
  CI_high[[i]]<-QL_CI[[i]][["CI_high"]]
  CI_low[[i]]<-QL_CI[[i]][["CI_low"]]
}
QL_CI_high<-as.data.frame(CI_high)
QL_CI_low<-as.data.frame(CI_low)
QL<-cbind(QL_med_vals,QL_CI_high,QL_CI_low)
colnames(QL)<-c('Median','CI_2.5','CI_97.5')
QL_s<-tidyr::gather(QL,"Curve","mortality",1:ncol(QL))
QL_s<-cbind(QL_s,derm_seq)

ggplot(QL_s, aes(dermdose, mortality, group = Curve)) +
  geom_line(aes(color=Curve))+
  ylab("Mortality") +
  xlab("Dermal Dose Estimates")+
  coord_cartesian(clip = "off") +
  theme_ridges()



#####LogProbit
LP<-mods3 %>%
  filter(model == 'LogProbit') %>%
  select(-c( model, seed, ppp, sd, p025, p50,p975,n_eff,rhat)) %>%
  tidyr:: spread(parameter, mean)

f(dose) = a + (1-a) * x * (c+b*(log(dose)))
# x is the cumulative distribution function of the standard normal distribution, denoted x 
LP<-mods3 %>%
  filter(model == 'LogProbit') 


#DichotomousHill
DH<-mods3 %>%
  filter(model == 'DichotomousHill') %>%
  select(-c( model, seed, ppp, sd, p025, p50,p975,n_eff,rhat)) %>%
  tidyr:: spread(parameter, mean)
DH_s<-list()
for (i in 1:nrow(DH)){
  en<-exp(1)
  DH_s[i]<- DH[i,2] * DH[i,5] + ((DH[i,2]-(DH[i,2] * DH[i,5]))/(1+en^((-DH[i,4])*(-DH[i,3])*log(derm_seq))))
}

t1<-a*g +((a-(a*g))/(1+en^(-c*-b*log(derm_seq))))



a<-DH[1,2]
b<-DH[1,3]
c<-DH[1,4]
g<-DH[1,5]


t1<-(DH[1,2] * DH[1,5]) + ((DH[1,2]-(DH[1,2] * DH[1,5]))/(1+en^((-DH[1,4])*(-DH[1,3])*log(derm_seq))))

t1<-DH[1,2] * DH[1,5] + ((DH[1,2]-(DH[1,2] * DH[1,5]))/1+en^((-DH[1,4])*(-DH[1,3])*log(derm_seq)))


DH_s_f<-as.data.frame(t(do.call(rbind.data.frame, DH_s)))
colnames(DH_s_f)<-paste0("mortality",1:ncol(DH_s_f),"")
row.names(DH_s_f)<-NULL
# QL_s_f$dose<-derm_seq
# QL_s_f<-QL_s_f[,c(1001,1:1000)]
DH_s_ft<-tidyr::gather(DH_s_f,"set","mortality",1:ncol(DH_s_f))
DH_s_ft<-cbind(DH_s_ft,derm_seq)
max(DH_s_ft$mortality)
ggplot(DH_s_ft, aes(dermdose, mortality, colour = mortality)) + geom_line()



ci_func <- function(x) sapply(x, ci)
mtcars  %>% filter(mpg > 20) %>%  myFunc()
sapply(data,ci,ci=0.95)


posterior <- distribution_normal(1000)

# Compute HDI and ETI
ci_hdi <- ci(posterior, method = "HDI")
ci_eti <- ci(posterior, method = "ETI")
#because these methods give the same for normal distribution, we can choose one or another. also pull out median


####Scraps ----
# #Median Set
# dose_m<-mean(dermal_dose)
# dose_sd<-sd(dermal_dose)
# sim_d<-rnorm(1000,dose_m, dose_sd)
# en<-exp(1)
# QL_S_med<- med$a + (1-med$a)*1-en^(-med$b*sim_d)
# #take a look
# plot(QL_S_med~sim_d)
# 
# #lower CI
# dose_m<-mean(dermal_dose)
# dose_sd<-sd(dermal_dose)
# sim_d<-rnorm(1000,dose_m, dose_sd)
# en<-exp(1)
# QL_S_med<- med$a + (1-med$a)*1-en^(-med$b*sim_d)
# #take a look
# plot(QL_S_med~sim_d)
# 
# #upper CI
# dose_m<-mean(dermal_dose)
# dose_sd<-sd(dermal_dose)
# sim_d<-rnorm(1000,dose_m, dose_sd)
# en<-exp(1)
# QL_S_med<- med$a + (1-med$a)*1-en^(-med$b*sim_d)
# #take a look
# plot(QL_S_med~sim_d)

# data<-QL[,c(2:3)]
# CI<-sapply(QL_s_f,ci, ci=0.95)
# med<-as.data.frame(t(sapply(data, function(x) median(x))))
# CI_25<-unlist(CI[2,])
# CI_95<-unlist(CI[3,]) 