


setwd('C:/Users/epauluko/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/GitHub/amphibian_effects_model')
lor<-within(read.csv('data_in/lor_soil.csv'),{
  Category<-as.factor(Category)
})





hsb2 <- within(read.csv("https://stats.idre.ucla.edu/stat/data/hsb2.csv"), {
  race <- as.factor(race)
  schtyp <- as.factor(schtyp)
  prog <- as.factor(prog)
})

attach(lor)


# null hypothesis: there will be no difference between categories 

summary(aov(Survival ~ Category + Application_Rate))
