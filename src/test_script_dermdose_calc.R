##test script for bb calculation using the 

pest<-read.csv('data_in/Headline_test.csv')
param<-read.csv('data_in/parameters.csv')

##for now let's assume 1cm is our best bet for even mixing, assuming a 1cm soil depth
##let's grab the Headline as is and go with the BW we have
pyra<-pest[1:10,]

#plot mortality as a function of body burden
#bb =  cSoil*e^(-kt)*kp*(SA/dt)Saf*DAF*h/BW
#SA= 10*bw^0.667
#Kp
#dt *
#Saf *
#BW
#H (bioavail) *
#DAF *
#cSoil is a conversion from application rate to estimate soil concentration

##calculate parameters, add in parameters from ABC
#for Kp: I need mol-weight and logKow info of Pyraclostrobin and Metconazole

#treat seperately - we'll be modeling added toxicity via the ABC route

#pyra
pmolweight<-387.8 #g/mol #comptox
plogKow<- 4.44 #comptox
#metconazole
mmolweight<-319.8 #g/mol #comptox
mlogKow<- 3.77 #comptox

kp_pyra =  10^(-2.72+(0.71*plogKow)-(0.0061*pmolweight))
kp_met =  10^(-2.72+(0.71*mlogKow)-(0.0061*mmolweight))

dt<-param[1,2]
move_rate<-param[3,2]
hl<-4.91*24 #from comptox profile
bioavail<-param[5,2]
dsa<-param[9,2]*pyra$Body_Weight_g^param[9,2] 
derm_frac<-param[11,2]
soil_concs_degradation<-log(2)/hl*pyra$Duration_h/move_rate
soil_concs<-((pyra$Application_Rate*10)/16000)*1000 #1cm 

dermal_dose<-(soil_concs^soil_concs_degradation * kp * (dsa/dt) * derm_frac * bioavail)/pyra$Body_Weight_g
dermal_dose 
pyra$dermaldose<-dermal_dose

fm1<-lm(Survival ~ dermaldose, data = pyra)
summary(fm1)

plot(Survival ~ dermaldose, data=pyra, type="p",
     xlab="Dermal Dose (ug/g)", ylab="Survival (%)")  

#move x axis to 0
#maybe we know the survival rate; we pick the application rate;
#bootstrap simulations
#each time we have a draw 

#1 cm even 
#sample size - if he has individual body weights that'd be great
#if not; we could always use the other body weights as a proxy for simulation 