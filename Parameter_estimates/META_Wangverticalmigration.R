###Vertical Migration Validation RScript for updated vertical migration parameters used in GLOWORM-META
#Olivia Ingle - olivia.ingle@liverpool.ac.uk (Last updated: 07/02/2025)
data <- read.csv("META_verticalmigration_data.csv")

##C. oncophora Parameters
Coop_logm1 <- log(data$Cooperia.m1)
mod.coop1 = lm(Coop_logm1 ~ data$Cooperia.Rainfall) 
print(summary(mod.coop1)) 
###Intercept = -3.35002, 0.05559

##O. ostertagia Parameters
Ost_logm1 <- log(data$Ostertagia.m1)
mod.ost1 = lm(Ost_logm1 ~ data$Ostertagia.Rainfall) 
print(summary(mod.ost1)) 
###Intercept = -3.29608, 0.04594 