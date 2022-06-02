##just test with a basic set of parameters for now

##this would not run when I input the parameters themselves; llogistic or LL.3 don't like my fixed parameters 
alpha=para[1,2]
beta=para[2,2]
gamma=para[3,2]

g=gamma
a=alpha
b=beta

dr<-by_s[[1]]
colnames(dr)[c(2,4)]<-c('dose','response')
dr$effect<-dr$response/dr$N

##none of this would work for me (except the basic version where we fit the data)
test_ll <- drm(response/N~dose, data=dr, fct = LL.3u(), type='binomial')
plot(test_ll, type='all')
out<-LL.3u(upper=1, fixed = c(para[1,2], para[2,2], para[3,2]), names = c("a", "b", "g"))
test_ll <- drm(response/N~dose, data=dr, fct = LL.3u(upper=1, fixed=c(alpha,beta,gamma)), type='binomial', )
plot(test_ll)
test_ll <- drm(response/N~dose, set, weights = N, data = dr, 
               fct = LL.3(fixed = c(alpha, beta, gamma), names = c("a", "b", "g")), type="binomial",
)






#here's the equation for the actual log logistic model from the BMDS User's Manual
log_logistic<-function(x){
  (g + ((1-g)/1+exp(-a-b*log(x))) )
  
}
p(dose) = log_logistic(dose) #this is the predicted probability distribution of the dose based on the fitted model



#here's an example of just running the basic logistic regression using GLM
model_list<-list()
for(model in 1:3){
  model=1
  dr = by_s[[model]]
  dr$response<-dr$Effect/dr$N
  
  log_dr<-glm(response ~ Dose, family=binomial(link='logit'), data=dr )
  model_list[[model]]<-log_dr
  
}


log_dr<-glm(response ~ Dose, family=binomial(link='logit'), data=dr )
dr2<-by_s[[4]]
dr2$response<-dr2$Effect/dr2$N

pred_dr<-predict(log_dr, newdata=dr2, se.fit=TRUE, type="link")
dr$p <- plogis(pred_dr$fit)
dr$lower <- plogis(pred_dr$fit - 1.96*pred_dr$se.fit)
dr$upper <- plogis(pred_dr$fit + 1.96*pred_dr$se.fit)
plot(response ~ Dose, data=dr2, type="l", ylim=c(0,1),xlab="Dose", ylab="Response")
points(response ~ Dose, dr2)
lines(response ~ Dose, dr, lty=2)
lines(response ~ Dose, dr, lty=2)

####scrap notes:

##just plot the fitted models (using some function for logistic regression), and stack as 95/5% CI
##we have the parameters, just need to put them in a function

#check range assumptions in this implementation of the drm function
#maybe log transformed?

#fitting and not predict function 

#(gamma+(1-gamma))/(1+exp{-[alpha+betaln(x)]}) #actual ll equation
#notes:
#see BMDS guidance

#not sure why gammma is wonky; beta and alpha look ok? (gamma between 0 and 1, but isn't)

##might need some guidance, but let me figure it out with the log-logistic model
##pulling top 2.5, 97.5 and 50% quantile BMDS values as upper and lower
##challenge: haven't plotted from scratch before, probably not as complicated as I'm making it

##other challenge: the guidance online is confusing


## example using fit 505 from the list
# parms <- ll_fit[[505]]$parameters
# data<-by_s[[505]]
# 
# g <- 1/(1+exp(-parms[1])); 
# a <- parms[2];
# b <- parms[3]; 
# d <- data$Dose
# rval_505 <- g + (1-g)*(1/(1+exp(-a-b*log(d))))
# 
# df<-as.data.frame(cbind(rval_505,d))
# names(df)<-c("effect","dose")
# ggplot(data = df, aes(x=dose, y=effect)) + geom_line()
# plot(ll_fit[[505]])



##here we extract the quantiles of the BMD; this did not ending up being that great
#note; quantile does not return the specific row values, so I used a 'closest' function to pick the datasets
# quantile_sims<-quantile(bmds$BMDSEstimates, probs = c(0.025,0.5,0.975), names=T)
# print(quantile_sims)
# 
# five_per<-bmds_bmd[which.min(abs(1.842985 -bmds_bmd$BMDSEstimates)),]
# fifty_per<-bmds_bmd[which.min(abs(2.520026-bmds_bmd$BMDSEstimates)),]
# ninetyfive_per<-bmds_bmd[which.min(abs(2.995469 -bmds_bmd$BMDSEstimates)),]
# perc<-rbind(five_per,fifty_per,ninetyfive_per)
# print(perc) #these are the data rows to pull to use the data to plot the curves

geom_line(aes(y=effect,size=1))+
  geom_ribbon(aes(y = effect, ymin =lerror, ymax = uerror), alpha = .2) +
  geom_point(data=og_data, aes(x=Dose,y=(Incidence/N)))



ggplot(data = df_percentiles, aes(x=x, y=val)) + geom_line(aes(colour=variable))






ll_laplace<- single_dichotomous_fit(og_data[,"Dose"],og_data[,"Incidence"],
                                    og_data[,"N"],model_type="log-logistic", fit_type="laplace")

og_ma<- ma_dichotomous_fit(og_data[,"Dose"],og_data[,"Incidence"],
                           og_data[,"N"],fit_type="mcmc")
plot(og_ma)



testy<-ll_fit$exp_n1$data
probs <- (0.5+testy[,2,drop=T])/(1.0 + testy[,3,drop=T])
se <- sqrt(probs*(1-probs)/testy[,3,drop=T])
uerror <- apply(cbind(probs*0+1,probs+se),1,min)
lerror <- apply(cbind(probs*0,probs-se),1,max)
doses =testy[,1,drop=T]
dose = c(doses,doses)
Response = c(uerror,lerror)



og_data<-as.data.frame(og_data)
ggplot(data = df, aes(x=dose, y=effect)) + geom_line() + geom_point(data=og_data, aes(x=Dose,y=(Incidence/N)))


ggplot(data = df, aes(x=dose)) + 
  geom_line(aes(y = effect), size = 1) +
  geom_ribbon(aes(y = effect, ymin =lerror, ymax = uerror), alpha = .2) +
  geom_point(data=og_data, aes(x=Dose,y=(Incidence/N)))

