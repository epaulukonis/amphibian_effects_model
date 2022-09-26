library(janitor)
library(tidyr)
setwd("C:\\Users\\EPAULUKO\\OneDrive - Environmental Protection Agency (EPA)\\Profile\\Documents\\Paulukonis_Documents\\manuscript_doseresponse_modeling")

df<-read.csv("TableS1.csv")

df<-df%>%clean_names()


unique(df$species_tested)
unique(df$study)
unique(df$chemical)
