#scraps ###---- 
#tests trend in proportions? in this example, there does seem to be a trend in proportions? hmm
smokers  <- c( 83, 90, 129, 70 )
patients <- c( 86, 93, 136, 82 )
prop.test(smokers, patients)
prop.trend.test(smokers, patients)
prop.trend.test(smokers, patients, c(0,0,0,1))
?prop.trend.test
#survival and study? 
# we want to test for statistical significance between individual studies
#do fairly simple lm between models? #ancova?
dose <- as.data.frame(matrix(c(10,9,10,7, 0,1,0,3), byrow=TRUE, nrow=2, dimnames=list(resp=0:1, dose=0:3)))
dose
CochranArmitageTest(dose, "increasing")
CochranArmitageTest(dose)
CochranArmitageTest(dose, "decreasing")
?CochranArmitageTest
# not exactly the same as in package coin:
# independence_test(tumor ~ dose, data = lungtumor, teststat = "quad")
lungtumor <- data.frame(dose = rep(c(0, 1, 2), c(40, 50, 48)),
                        tumor = c(rep(c(0, 1), c(38, 2)),
                                  rep(c(0, 1), c(43, 7)),
                                  rep(c(0, 1), c(33, 15))))


tab <- table(lungtumor$dose, lungtumor$tumor)
CochranArmitageTest(tab, 'increasing')
str(tab)

out_1<-t(by_s_stat[[1]])
out_1
colnames(out_1)<-out[1,]
out<-out[2:3,]
row.names(out)<-c(0,1)
out<-as.matrix(out)
out
dose

CochranArmitageTest(out,'increasing')


#format the dermal dose estimates and survival estimates to use in the BBMD app
#should be in order with effects, so that each experiment N matches with each study
mort_n<-as.data.frame(sapply(survival_sim, function(x) effects$N_Exp - round(x*effects$N_Exp,0))) #number of individuals dying per exp
#we are treating each batch like its own separate study; essentially, each set of doses has a corresponding set of survival data
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
exp_n<-as.data.frame(rep.col(effects$N_Exp,nsims))
colnames(exp_n)<-paste0("exp_n",1:nsims,"")
#need ID, dose, N, and mortality
#for 100 sets of 39 rows
exp_n<-gather(exp_n,"set","N",1:nsims)
mort_n<-gather(mort_n,"set","Effect", 1:nsims)
derm_dose<-gather(derm,"set","Dose",1:nsims)

BMDS_headline_fin<-cbind(exp_n,mort_n,derm_dose) #name the output whatever run you'd like
BMDS_headline_fin<-BMDS_headline_fin[,c(1,6,2,4)]

#order by dose for clarity
BMDS_headline_fin<-BMDS_headline_fin %>%
  group_by(set) %>% #set of sims
  arrange(Dose, .by_group=T) #this will rearrange the output in terms of order of sims, but will still work fine








#let's test that I can make this work



###only needed for BBMD, BMDS can do it in one go
#we need to split up the 1000 set simulation into 10 sets of 100
# bbmd<-split(BMDS_Run3, (as.numeric(rownames(BMDS_Run3))-1) %/% 3900)
# lapply(bbmd,dim) #check dimensions of the output; should be 3900*4 for each set of 10 (0-9)

#quick test if needed
# test<-bbmd
# testy<-test[[3]] #check that it's formatted with 100 sets of unique sims
# unique(testy$set) #should be 100

# colnames<-c("Dataset Index", "Dose", "N","Effect")
# bbmd<-lapply(bbmd, setNames, colnames)
##for BBMD format
# for(i in names(bbmd)){
#   write.csv(bbmd[[i]], paste0(i,"bbmdrun.csv"), row.names=FALSE)
# }







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


