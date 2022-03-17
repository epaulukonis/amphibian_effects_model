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