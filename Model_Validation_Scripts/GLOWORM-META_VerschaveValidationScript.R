# Author: Olivia Ingle - olivia.ingle@liverpool.ac.uk /// Hannah Vineer - hannah.vineer@liverpool.ac.uk
# Created: 08/08/2024

##GLOWORM-META Validation Script for Ostertagia ostertagi and Cooperia oncophora for the Verschave thesis data for the Belgium herds. 

##############Source Files/Packages
# Load in packages and required source files -----------------------------------
packages = c("deSolve", "geosphere", "chron", "ncdf4", "imputeTS", "forecast", "gplots", "tidyr", "rstudioapi") # Package names


# Check for, and install missing package
installed_packages = packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE)) # Load packages

setwd(dirname(getActiveDocumentContext()$path)) ###saves outputs to the folder where your CSV is located
getwd()

# Source model functions -------------------------------------------------------
source('Model_functions/livestockfuns.R') # liveweight, faeces production and dry matter intake
source('Model_functions/weatherfuns.R') # EOBS data extraction function
source('Model_functions/parainterp.R') # Interpolation function for parasitological input
source('Model_functions/gloworm_meta.R') # GLOWORM-META function
source('Model_functions/initialvalues.r') # sets initial conditions
source('Model_functions/ginparms.r') # species-specific parameters

# Read in validation data ####
Verschave <- read.csv("ValidationData/VerschaveData.csv")
Versch_Herd1 = subset(x = Verschave, subset = Verschave$HerdID=="Herd1")
Versch_Herd2 = subset(x = Verschave, subset = Verschave$HerdID=="Herd2")
Versch_Herd3 = subset(x = Verschave, subset = Verschave$HerdID=="Herd3")
Versch_Herd4 = subset(x = Verschave, subset = Verschave$HerdID=="Herd4")
Versch_Herd5 = subset(x = Verschave, subset = Verschave$HerdID=="Herd5")
Versch_Herd6 = subset(x = Verschave, subset = Verschave$HerdID=="Herd6")
Versch_Herd7 = subset(x = Verschave, subset = Verschave$HerdID=="Herd7")

#Average herbage biomass over the simulation (kilograms dry matter per hectare) **This is estimated as no value given***
kgDM = 2000

# Host species:
# Choose from:
# 1 = sheep
# 2 = cow
host = 2




##############Herd Runs

#Herd 1 Set-up #####
#Seperate Ostertagia and Cooperia data
Versch_Herd1Ost = subset(x = Versch_Herd1, subset = Versch_Herd1$Nematode=="Ostertagia")
Versch_Herd1Coop = subset(x = Versch_Herd1, subset = Versch_Herd1$Nematode=="Cooperia")
# Enter start and end dates:
ss = '2012-04-20'
ee = '2012-11-09'

# Enter site location:
lat = 51.2860833
lon = 3.2677981

#Enter temperature and rainfall data, using the eobs() function defined above: 
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
precip = eobspoint(data = 'rr_ens_mean_0.1deg_reg_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)

# Host age at beginning of simulation (in days) (estimated 20 months of age **unsure on exact age of cows**)
host_age = 600

# Add events 
# Treatment added on the 2012-09-07. 
##treatment was moxidectin pour-on: 100% effective against all stages and has a residual activity of 35 days
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(141, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat
# Estimate weights and use this to estimate dry matter intake and faeces production
# Generate a vector of ages to match the length of the simulation
host.age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1)
# Use this to make a vector of estimated host weight using the liveweight function
lw = if(host == 1) {lw_sheep(age = host.age)} else {lw_cattle(age = host.age)}

# Use the weight to estimate dry matter intake
DMI = dmi(lwt = lw, age = host.age, host_species = host)

# Create a vector of DM values for each day of the simulation
kgDMha = rep(kgDM, length(DMI))

# The amount of faeces produced is not important for the simulation but is used in posthoc estimates of faecal egg counts from the simulation output.
f = faeces(lwt = lw, host_species = host)
###No exact stocking rate given --- used Herd 3 averages instead
stocking_rate = c(rep(5.2, length(seq.Date(as.Date(ss), as.Date('2012-08-03'), 'day'))), rep(3, length(temp)-length(seq.Date(as.Date(ss), as.Date('2012-08-03'), 'day'))))
movement_dates = c(ss, '2012-08-01', ee)  #must include the start and end dates
grazing_plots = c(1,0,0) # must include the start and end locations


#Herd 1 Ostertagia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.5, Preadult_in_host = 13877.04) ###calculated using the eggs to worms calculation. Used the second timepoint egg count. 


Herd1Oster = gloworm_meta(start = ss, end = ee, lat = lat, 
                     temp = temp, precip = precip, 
                     statevars = initial_values, host = 2, nematode = 3,
                     stocking_rate = stocking_rate, graze = grazing_plots,
                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                     lwt = lw, faeces = f, eventdat = eventdat)
Herd1Oster_days <- 1:nrow(Herd1Oster) + 20##shifts FEC to start date for graphing
plot(Herd1Oster_days, Herd1Oster[,'FEC'], ylim = c(0, 500), xlim = c(1, length(Versch_Herd1Ost$Day) + 243),  main = 'Verschave Herd 1', xaxt = 'n', xlab = 'Month', ylab = expression(paste('FEC ', italic('O. ostertagia'))), type ='l', col = 'grey')
axis(side = 1, at = c(1, 30, 61, 91, 122, 153, 183, 213, 243), labels = c('Apr','May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd1Ost$day+20), ui = Versch_Herd1Ost$mean+((Versch_Herd1Ost$upperCI-Versch_Herd1Ost$mean)), li = Versch_Herd1Ost$mean-((Versch_Herd1Ost$mean-Versch_Herd1Ost$lowerCI)), add = TRUE, y = Versch_Herd1Ost$mean, pch = 16)


validation_output1ostertagia = data.frame(day = Versch_Herd1Ost$day, 
                                simulated = Herd1Oster[,'FEC'][Versch_Herd1Ost$day],
                                observed.mean = Versch_Herd1Ost$mean,
                                observed.uci = Versch_Herd1Ost$mean-((Versch_Herd1Ost$mean-Versch_Herd1Ost$lowerCI)),
                                observed.lci = Versch_Herd1Ost$mean+((Versch_Herd1Ost$upperCI-Versch_Herd1Ost$mean)),
                                species = rep('Ostertagia', length(Versch_Herd1Ost$day)),
                                study = rep('Verschave', length(Versch_Herd1Ost$day)), 
                                group = rep("Herd 1", length(Versch_Herd1Ost$day)))
dev.off()

#Herd 1 Cooperia####

initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.4,  Preadult_in_hostA = 294.5377) ###calculated using the eggs to worms calculation. Used the second timepoint egg count. 

Herd1Coop = gloworm_meta(start = ss, end = ee, lat = lat, 
                          temp = temp, precip = precip, 
                          statevars = initial_values, host = 2, nematode = 4,
                          stocking_rate = stocking_rate, graze = grazing_plots,
                          movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                          lwt = lw, faeces = f, eventdat = eventdat)
Herd1Coopdays <- 1:nrow(Herd1Coop) + 20
plot(Herd1Coopdays, Herd1Coop[,'FEC'], ylim = c(0, 800), xlim = c(1, length(Versch_Herd1Coop$Day) + 243), xaxt = 'n', main = 'Verschave Herd 1', xlab = 'Month', ylab = expression(paste('FEC ', italic('C. oncophora'))), type ='l', col = 'grey')
axis(side = 1, at = c(1, 30, 61, 91, 122, 153, 183, 213, 243), labels = c('Apr','May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd1Coop$day+20), ui = Versch_Herd1Coop$mean+((Versch_Herd1Coop$upperCI-Versch_Herd1Coop$mean)), li = Versch_Herd1Coop$mean-((Versch_Herd1Coop$mean-Versch_Herd1Coop$lowerCI)), add = TRUE, y = Versch_Herd1Coop$mean, pch = 16)

validation_output1cooperia = data.frame(day = Versch_Herd1Coop$day, 
                                          simulated = Herd1Coop[,'FEC'][Versch_Herd1Coop$day],
                                          observed.mean = Versch_Herd1Coop$mean,
                                          observed.uci = Versch_Herd1Coop$mean-((Versch_Herd1Coop$mean-Versch_Herd1Coop$lowerCI)),
                                          observed.lci = Versch_Herd1Coop$mean+((Versch_Herd1Coop$upperCI-Versch_Herd1Coop$mean)),
                                          species = rep('Cooperia', length(Versch_Herd1Coop$day)),
                                          study = rep('Verschave', length(Versch_Herd1Coop$day)), 
                                          group = rep("Herd 1", length(Versch_Herd1Coop$day)))

dev.off()



#Herd 2 Set-up ####
Versch_Herd2Ost = subset(x = Versch_Herd2, subset = Versch_Herd2$Nematode=="Ostertagia")
Versch_Herd2Coop = subset(x = Versch_Herd2, subset = Versch_Herd2$Nematode=="Cooperia")

# Enter start and end dates:
ss = "2012-05-15"
ee = "2012-09-21"

# Enter site location:
lat = 51.2846915
lon = 4.6641447

#Enter temperature and rainfall data, using the eobs() function defined above: 
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
precip = eobspoint(data = 'rr_ens_mean_0.1deg_reg_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)

# Host age at beginning of simulation (in days) (estimated 6 months of age - Verschave thesis)
host_age = 180

# Add events 
# Treatment added on the 2012-09-07. 
#treatment was Doramectin injectable formulation 100% effective against all stages and has a residual activity of 1-2wks against Co
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(119, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

# Estimate weights and use this to estimate dry matter intake and faeces production
# Generate a vector of ages to match the length of the simulation
host.age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1)
# Use this to make a vector of estimated host weight using the liveweight function
lw = if(host == 1) {lw_sheep(age = host.age)} else {lw_cattle(age = host.age)}

# Use the weight to estimate dry matter intake
DMI = dmi(lwt = lw, age = host.age, host_species = host)

# Create a vector of DM values for each day of the simulation
kgDMha = rep(kgDM, length(DMI))

# The amount of faeces produced is not important for the simulation but is used in posthoc estimates of faecal egg counts from the simulation output.
f = faeces(lwt = lw, host_species = host)

###No exact stocking rate given --- used Herd 3 averages instead
stocking_rate = c(rep(10, length(seq.Date(as.Date(ss), as.Date('2012-09-06'), 'day'))), rep(10, length(temp)-length(seq.Date(as.Date(ss), as.Date('2012-09-06'), 'day'))))
movement_dates = c(ss, '2012-09-07', '2012-09-08', ee)  #must include the start and end dates
grazing_plots = c(1,0,0,0) # must include the start and end locations

#Herd 2 Ostertagia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.1, Preadult_in_host = 3205.178) ###calculated using the eggs to worms calculation. Used the second timepoint egg count. 
Herd2Oster = gloworm_meta(start = ss, end = ee, lat = lat, 
                          temp = temp, precip = precip, 
                          statevars = initial_values, host = 2, nematode = 3,
                          stocking_rate = stocking_rate, graze = grazing_plots,
                          movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                          lwt = lw, faeces = f, eventdat = eventdat)
Herd2Osterdays <- 1:nrow(Herd2Oster) + 15
plot(Herd2Osterdays , Herd2Oster[,'FEC'], ylim = c(0, 500), xlim = c(-30, length(Versch_Herd2Ost$day) + 213), xaxt = 'n', main = 'Verschave Herd 2', xlab = 'Month', ylab = expression(paste('FEC ', italic('O. ostertagia'))), type ='l', col = 'grey')
axis(side = 1, at = c(-30, 1, 31, 61, 91, 122, 153, 183, 213), labels = c('Apr', 'May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd2Ost$day+15), ui = Versch_Herd2Ost$mean+((Versch_Herd2Ost$upperCI-Versch_Herd2Ost$mean)), li = Versch_Herd2Ost$mean-((Versch_Herd2Ost$mean-Versch_Herd2Ost$lowerCI)), add = TRUE, y = Versch_Herd2Ost$mean, pch = 16)

validation_output2ostertagia = data.frame(day = Versch_Herd2Ost$day, 
                                          simulated = Herd2Oster[,'FEC'][Versch_Herd2Ost$day],
                                          observed.mean = Versch_Herd2Ost$mean,
                                          observed.uci = Versch_Herd2Ost$mean-((Versch_Herd2Ost$mean-Versch_Herd2Ost$lowerCI)),
                                          observed.lci = Versch_Herd2Ost$mean+((Versch_Herd2Ost$upperCI-Versch_Herd2Ost$mean)),
                                          species = rep('Ostertagia', length(Versch_Herd2Ost$day)),
                                          study = rep('Verschave', length(Versch_Herd2Ost$day)), 
                                          group = rep("Herd 2", length(Versch_Herd2Ost$day)))
dev.off()

#Herd 2 Cooperia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.2,  Preadult_in_host = 1790.729) ###calculated using the eggs to worms calculation. Used the second timepoint egg count. 

Herd2Coop = gloworm_meta(start = ss, end = ee, lat = lat, 
                         temp = temp, precip = precip, 
                         statevars = initial_values, host = 2, nematode = 4,
                         stocking_rate = stocking_rate, graze = grazing_plots,
                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                         lwt = lw, faeces = f, eventdat = eventdat)

Herd2Coopdays <- 1:nrow(Herd2Coop) + 45
plot(Herd2Coopdays, Herd2Coop[,'FEC'], ylim = c(0, 800), xaxt = 'n', xlim = c(1, length(Versch_Herd2Coop$day) + 243), main = 'Verschave Herd 2',  xlab = 'Month', ylab = expression(paste(italic('C. oncophora'), 'FEC')),  type ='l', col = 'grey')
axis(side = 1, at = c(1, 30, 61, 91, 122, 153, 183, 213, 243), labels = c('Apr','May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd2Coop$day+45), ui = Versch_Herd2Coop$mean+((Versch_Herd2Coop$upperCI-Versch_Herd2Coop$mean)), li = Versch_Herd2Coop$mean-((Versch_Herd2Coop$mean-Versch_Herd2Coop$lowerCI)), add = TRUE, y = Versch_Herd2Coop$mean, pch = 16)

validation_output2cooperia = data.frame(day = Versch_Herd2Coop$day, 
                                        simulated = Herd2Coop[,'FEC'][Versch_Herd2Coop$day],
                                        observed.mean = Versch_Herd2Coop$mean,
                                        observed.uci = Versch_Herd2Coop$mean-((Versch_Herd2Coop$mean-Versch_Herd2Coop$lowerCI)),
                                        observed.lci = Versch_Herd2Coop$mean+((Versch_Herd2Coop$upperCI-Versch_Herd2Coop$mean)),
                                        species = rep('Cooperia', length(Versch_Herd2Coop$day)),
                                        study = rep('Verschave', length(Versch_Herd2Coop$day)), 
                                        group = rep("Herd 2", length(Versch_Herd2Coop$day)))

dev.off()


#Herd 3 Set-up####
#Seperate Ostertagia and Cooperia data
Versch_Herd3Ost = subset(x = Versch_Herd3, subset = Versch_Herd3$Nematode=="Ostertagia")
Versch_Herd3Coop = subset(x = Versch_Herd3, subset = Versch_Herd3$Nematode=="Cooperia")

# Enter start and end dates:
ss = "2013-06-14"
ee = "2013-10-15"

# Enter site location:
lat = 51.1169567
lon = 3.5585687

# Enter temperature and rainfall data, using the eobs() function defined above:
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
precip = eobspoint(data = 'rr_ens_mean_0.1deg_reg_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)

# Host age at beginning of simulation (in days) (estimated 21 months of age - Verschave thesis)
host_age = 630


# Add events 
# NO TREATMENTS
eventdat = NULL
eventdat

# Estimate weights and use this to estimate dry matter intake and faeces production
# Generate a vector of ages to match the length of the simulation
host.age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1)
# Use this to make a vector of estimated host weight using the liveweight function
lw = if(host == 1) {lw_sheep(age = host.age)} else {lw_cattle(age = host.age)}

# Use the weight to estimate dry matter intake
DMI = dmi(lwt = lw, age = host.age, host_species = host)

# Create a vector of DM values for each day of the simulation
kgDMha = rep(kgDM, length(DMI))

# The amount of faeces produced is not important for the simulation but is used in posthoc estimates of faecal egg counts from the simulation output.
f = faeces(lwt = lw, host_species = host)

###Stocking rate  
stocking_rate = c(rep(5.2, length(seq.Date(as.Date(ss), as.Date('2013-08-20'), 'day'))), rep(2.5, length(temp)-length(seq.Date(as.Date(ss), as.Date('2013-08-20'), 'day'))))
movement_dates = c(ss, '2013-08-03', '2013-09-20', '2013-09-21', ee)  #must include the start and end dates
grazing_plots = c(1,1,0,0,0) # must include the start and end locations

#Herd 3 Ostertagia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.3,  Preadult_in_host = 1433.585) 

Herd3Oster = gloworm_meta(start = ss, end = ee, lat = lat, 
                      temp = temp, precip = precip, 
                      statevars = initial_values, host = 2, nematode = 3,
                      stocking_rate = stocking_rate, graze = grazing_plots,
                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                      lwt = lw, faeces = f, eventdat = eventdat)

Herd3Osterday <- 1:nrow(Herd3Oster) + 14
plot(Herd3Oster[,'FEC'], ylim = c(0, 500), xlim = c(-60, length(Versch_Herd3Ost$day) + 183), xaxt = 'n', main = 'Verschave Herd 3', xlab = 'Month', ylab = expression(paste('FEC ', italic('O. ostertagia'))), type ='l', col = 'grey')
axis(side = 1, at = c(-60, -30,1,31, 61, 91, 122, 153, 183), labels = c('Apr', 'May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd3Ost$day), ui = Versch_Herd3Ost$mean+((Versch_Herd3Ost$upperCI-Versch_Herd3Ost$mean)), li = Versch_Herd3Ost$mean-((Versch_Herd3Ost$mean-Versch_Herd3Ost$lowerCI)), add = TRUE, y = Versch_Herd3Ost$mean, pch = 16)



validation_output3ostertagia = data.frame(day = Versch_Herd3Ost$day, 
                                          simulated = Herd3Oster[,'FEC'][Versch_Herd3Ost$day],
                                          observed.mean = Versch_Herd3Ost$mean,
                                          observed.uci = Versch_Herd3Ost$mean-((Versch_Herd3Ost$mean-Versch_Herd3Ost$lowerCI)),
                                          observed.lci = Versch_Herd3Ost$mean+((Versch_Herd3Ost$upperCI-Versch_Herd3Ost$mean)),
                                          species = rep('Ostertagia', length(Versch_Herd3Ost$day)),
                                          study = rep('Verschave', length(Versch_Herd3Ost$day)), 
                                          group = rep("Herd 3", length(Versch_Herd3Ost$day)))
dev.off()

#Herd 3 Cooperia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.5,  Adult_in_host = 10) ###calculated using the eggs to worms calculation. Used the second timepoint egg count. 
Herd3Coop = gloworm_meta(start = ss, end = ee, lat = lat, 
                         temp = temp, precip = precip, 
                         statevars = initial_values, host = 2, nematode = 4,
                         stocking_rate = stocking_rate, graze = grazing_plots,
                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                         lwt = lw, faeces = f, eventdat = eventdat)

Herd3Coopdays <- 1:nrow(Herd3Coop) +74
plot(Herd3Coopdays, Herd3Coop[,'FEC'], ylim = c(0, 800), xlim = c(1, length(Versch_Herd3Coop$day) + 243), main = 'Verschave Herd 3', xaxt = 'n', xlab = 'Month', ylab = expression(paste(italic('C. oncophora'), 'FEC')),  type ='l', col = 'grey')
axis(side = 1, at = c(1, 30, 61, 91, 122, 153, 183, 213, 243), labels = c('Apr','May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd3Coop$day+74), ui = Versch_Herd3Coop$mean+((Versch_Herd3Coop$upperCI-Versch_Herd3Coop$mean)), li = Versch_Herd3Coop$mean-((Versch_Herd3Coop$mean-Versch_Herd3Coop$lowerCI)), add = TRUE, y = Versch_Herd3Coop$mean, pch = 16)

validation_output3cooperia = data.frame(day = Versch_Herd3Coop$day, 
                                        simulated = Herd3Coop[,'FEC'][Versch_Herd3Coop$day],
                                        observed.mean = Versch_Herd3Coop$mean,
                                        observed.uci = Versch_Herd3Coop$mean-((Versch_Herd3Coop$mean-Versch_Herd3Coop$lowerCI)),
                                        observed.lci = Versch_Herd3Coop$mean+((Versch_Herd3Coop$upperCI-Versch_Herd3Coop$mean)),
                                        species = rep('Cooperia', length(Versch_Herd3Coop$day)),
                                        study = rep('Verschave', length(Versch_Herd3Coop$day)), 
                                        group = rep("Herd 3", length(Versch_Herd3Coop$day)))



#Herd 4 Set-up####
#Seperate Ostertagia and Cooperia data
Versch_Herd4Ost = subset(x = Versch_Herd4, subset = Versch_Herd4$Nematode=="Ostertagia")
Versch_Herd4Coop = subset(x = Versch_Herd4, subset = Versch_Herd4$Nematode=="Cooperia")

# Enter start and end dates:
ss = '2013-04-20'
ee = '2013-11-26'

# Enter site location:
lat = 50.85
lon = 3.636457

# Enter temperature and rainfall data, using the eobs() function defined above:
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
precip = eobspoint(data = 'rr_ens_mean_0.1deg_reg_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)
# Add events 
# NO TREATMENTS
eventdat = NULL
eventdat

# Estimate weights and use this to estimate dry matter intake and faeces production
# Generate a vector of ages to match the length of the simulation
host.age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1)
# Use this to make a vector of estimated host weight using the liveweight function
lw = if(host == 1) {lw_sheep(age = host.age)} else {lw_cattle(age = host.age)}

# Use the weight to estimate dry matter intake
DMI = dmi(lwt = lw, age = host.age, host_species = host)

# Create a vector of DM values for each day of the simulation
kgDMha = rep(kgDM, length(DMI))

# The amount of faeces produced is not important for the simulation but is used in posthoc estimates of faecal egg counts from the simulation output.
f = faeces(lwt = lw, host_species = host)

#stocking rates
stocking_rate = c(rep(13, length(seq.Date(as.Date(ss), as.Date('2013-08-20'), 'day'))), rep(5, length(temp)-length(seq.Date(as.Date(ss), as.Date('2013-08-20'), 'day'))))
movement_dates = c(ss, '2013-08-20', '2013-11-24', ee)  #must include the start and end dates
grazing_plots = c(1,0,0,0) # must include the start and end locations

#Herd 4 Ostertagia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.3,  Preadult_in_host = 213.1807) ###calculated using the eggs to worms calculation. Used the second timepoint egg count. 

Herd4Oster = gloworm_meta(start = ss, end = ee, lat = lat, 
                      temp = temp, precip = precip, 
                      statevars = initial_values, host = 2, nematode = 3,
                      stocking_rate = stocking_rate, graze = grazing_plots,
                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                      lwt = lw, faeces = f, eventdat = eventdat)

Herd4Osterdays <- 1:nrow(Herd4Oster) + 15
plot(Herd4Osterdays , Herd4Oster[,'FEC'], ylim = c(0, 500), xlim = c(-30, length(Versch_Herd4Ost$day) + 213), xaxt = 'n', main = 'Verschave Herd 4', xlab = 'Month', ylab = expression(paste('FEC ', italic('O. ostertagia'))), type ='l', col = 'grey')
axis(side = 1, at = c(-30, 1, 31, 61, 91, 122, 153, 183, 213), labels = c('Apr', 'May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd4Ost$day+15), ui = Versch_Herd4Ost$mean+((Versch_Herd4Ost$upperCI-Versch_Herd4Ost$mean)), li = Versch_Herd4Ost$mean-((Versch_Herd4Ost$mean-Versch_Herd4Ost$lowerCI)), add = TRUE, y = Versch_Herd4Ost$mean, pch = 16)


validation_output4ostertagia = data.frame(day = Versch_Herd4Ost$day, 
                                          simulated = Herd4Oster[,'FEC'][Versch_Herd4Ost$day],
                                          observed.mean = Versch_Herd4Ost$mean,
                                          observed.uci = Versch_Herd4Ost$mean-((Versch_Herd4Ost$mean-Versch_Herd4Ost$lowerCI)),
                                          observed.lci = Versch_Herd4Ost$mean+((Versch_Herd4Ost$upperCI-Versch_Herd4Ost$mean)),
                                          species = rep('Ostertagia', length(Versch_Herd4Ost$day)),
                                          study = rep('Verschave', length(Versch_Herd4Ost$day)), 
                                          group = rep("Herd 4", length(Versch_Herd4Ost$day)))
dev.off()

#Herd 4 Cooperia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.6,  Preadult_in_host = 10.37104) ###calculated using the eggs to worms calculation. Used the second timepoint egg count. 

Herd4Coop = gloworm_meta(start = ss, end = ee, lat = lat, 
                         temp = temp, precip = precip, 
                         statevars = initial_values, host = 2, nematode = 4,
                         stocking_rate = stocking_rate, graze = grazing_plots,
                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                         lwt = lw, faeces = f, eventdat = eventdat)

Herd4Coopdays <- 1:nrow(Herd4Coop) + 45
plot(Herd4Coopdays, Herd4Coop[,'FEC'], ylim = c(0, 800), xaxt = 'n', xlim = c(1, length(Versch_Herd4Coop$day) + 243), main = 'Verschave Herd 4',  xlab = 'Month', ylab = expression(paste(italic('C. oncophora'), 'FEC')),  type ='l', col = 'grey')
axis(side = 1, at = c(1, 30, 61, 91, 122, 153, 183, 213, 243), labels = c('Apr','May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd4Coop$day+45), ui = Versch_Herd4Coop$mean+((Versch_Herd4Coop$upperCI-Versch_Herd4Coop$mean)), li = Versch_Herd4Coop$mean-((Versch_Herd4Coop$mean-Versch_Herd4Coop$lowerCI)), add = TRUE, y = Versch_Herd4Coop$mean, pch = 16)

validation_output4cooperia = data.frame(day = Versch_Herd4Coop$day, 
                                        simulated = Herd4Coop[,'FEC'][Versch_Herd4Coop$day],
                                        observed.mean = Versch_Herd4Coop$mean,
                                        observed.uci = Versch_Herd4Coop$mean-((Versch_Herd4Coop$mean-Versch_Herd4Coop$lowerCI)),
                                        observed.lci = Versch_Herd4Coop$mean+((Versch_Herd4Coop$upperCI-Versch_Herd4Coop$mean)),
                                        species = rep('Cooperia', length(Versch_Herd4Coop$day)),
                                        study = rep('Verschave', length(Versch_Herd4Coop$day)), 
                                        group = rep("Herd 4", length(Versch_Herd4Coop$day)))

dev.off()


#Herd 5 Set-up####
#Seperate Ostertagia and Cooperia data
Versch_Herd5Ost = subset(x = Versch_Herd5, subset = Versch_Herd5$Nematode=="Ostertagia")
Versch_Herd5Coop = subset(x = Versch_Herd5, subset = Versch_Herd5$Nematode=="Cooperia")

# Enter start and end dates:
ss = "2013-06-07"
ee = "2013-10-23"

# Enter site location:
lat = 51.161453
lon = 3.6421456

# Enter temperature and rainfall data, using the eobs() function defined above:
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
precip = eobspoint(data = 'rr_ens_mean_0.1deg_reg_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)

# Host age at beginning of simulation (in days) (estimated 15 months of age - Verschave thesis)
host_age = 330

# Add events 
# NO TREATMENTS
eventdat = NULL
eventdat

# Estimate weights and use this to estimate dry matter intake and faeces production
# Generate a vector of ages to match the length of the simulation
host.age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1)
# Use this to make a vector of estimated host weight using the liveweight function
lw = if(host == 1) {lw_sheep(age = host.age)} else {lw_cattle(age = host.age)}

# Use the weight to estimate dry matter intake
DMI = dmi(lwt = lw, age = host.age, host_species = host)

# Create a vector of DM values for each day of the simulation
kgDMha = rep(kgDM, length(DMI))

# The amount of faeces produced is not important for the simulation but is used in posthoc estimates of faecal egg counts from the simulation output.
f = faeces(lwt = lw, host_species = host)

###Stocking rates from Verschave thesis, 2015. -- starting stocking rate high -- assuming when in barn stocking rate is lower. 
stocking_rate = c(rep(62.4, length(seq.Date(as.Date(ss), as.Date('2013-10-05'), 'day'))), rep(20, length(temp)-length(seq.Date(as.Date(ss), as.Date('2013-10-05'), 'day'))))
movement_dates = c(ss, '2013-10-05','2013-10-14', ee)  #must include the start and end dates
grazing_plots = c(1,0,0,0) # must include the start and end locations

#Herd 5 Ostertagia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.2,  Preadult_in_host = 310.081) ###calculated using the eggs to worms calculation. Used the third timepoint egg count. 

Herd5Oster = gloworm_meta(start = ss, end = ee, lat = lat, 
                          temp = temp, precip = precip, 
                          statevars = initial_values, host = 2, nematode = 3,
                          stocking_rate = stocking_rate, graze = grazing_plots,
                          movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                          lwt = lw, faeces = f, eventdat = eventdat)

Herd5Osterdays <- 1:nrow(Herd5Oster) + 7
plot(Herd5Osterdays, Herd5Oster[,'FEC'], ylim = c(0, 500), xlim = c(-60, length(Versch_Herd5Ost$day) + 183),  xaxt = 'n', main = 'Verschave Herd 5', xlab = 'Month', ylab = expression(paste('FEC ', italic('O. ostertagia'))), type ='l', col = 'grey')
axis(side = 1, at = c(-60, -30,1,31, 61, 91, 122, 153, 183), labels = c('Apr', 'May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd5Ost$day+7), ui = Versch_Herd5Ost$mean+((Versch_Herd5Ost$upperCI-Versch_Herd5Ost$mean)), li = Versch_Herd5Ost$mean-((Versch_Herd5Ost$mean-Versch_Herd5Ost$lowerCI)), add = TRUE, y = Versch_Herd5Ost$mean, pch = 16)


validation_output5ostertagia = data.frame(day = Versch_Herd5Ost$day, 
                                          simulated = Herd5Oster[,'FEC'][Versch_Herd5Ost$day],
                                          observed.mean = Versch_Herd5Ost$mean,
                                          observed.uci = Versch_Herd5Ost$mean-((Versch_Herd5Ost$mean-Versch_Herd5Ost$lowerCI)),
                                          observed.lci = Versch_Herd5Ost$mean+((Versch_Herd5Ost$upperCI-Versch_Herd5Ost$mean)),
                                          species = rep('Ostertagia', length(Versch_Herd5Ost$day)),
                                          study = rep('Verschave', length(Versch_Herd5Ost$day)), 
                                          group = rep("Herd 5", length(Versch_Herd5Ost$day)))
dev.off()

#Herd 5 Cooperia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.4,  Preadult_in_host = 10.37104)

Herd5Coop = gloworm_meta(start = ss, end = ee, lat = lat, 
                         temp = temp, precip = precip, 
                         statevars = initial_values, host = 2, nematode = 4,
                         stocking_rate = stocking_rate, graze = grazing_plots,
                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                         lwt = lw, faeces = f, eventdat = eventdat)

Herd5Coopdays <- 1:nrow(Herd5Coop) + 45
plot(Herd5Coopdays, Herd5Coop[,'FEC'], ylim = c(0, 800), xaxt = 'n', xlim = c(1, length(Versch_Herd5Coop$day) + 243), main = 'Verschave Herd 5',  xlab = 'Month', ylab = expression(paste(italic('C. oncophora'), 'FEC')),  type ='l', col = 'grey')
axis(side = 1, at = c(1, 30, 61, 91, 122, 153, 183, 213, 243), labels = c('Apr','May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd5Coop$day+45), ui = Versch_Herd5Coop$mean+((Versch_Herd5Coop$upperCI-Versch_Herd5Coop$mean)), li = Versch_Herd5Coop$mean-((Versch_Herd5Coop$mean-Versch_Herd5Coop$lowerCI)), add = TRUE, y = Versch_Herd5Coop$mean, pch = 16)

validation_output5cooperia = data.frame(day = Versch_Herd5Coop$day, 
                                        simulated = Herd5Coop[,'FEC'][Versch_Herd5Coop$day],
                                        observed.mean = Versch_Herd5Coop$mean,
                                        observed.uci = Versch_Herd5Coop$mean-((Versch_Herd5Coop$mean-Versch_Herd5Coop$lowerCI)),
                                        observed.lci = Versch_Herd5Coop$mean+((Versch_Herd5Coop$upperCI-Versch_Herd5Coop$mean)),
                                        species = rep('Cooperia', length(Versch_Herd5Coop$day)),
                                        study = rep('Verschave', length(Versch_Herd5Coop$day)), 
                                        group = rep("Herd 5", length(Versch_Herd5Coop$day)))
dev.off()



#Herd 6 Set-up####
#Seperate Ostertagia and Cooperia data
Versch_Herd6Ost = subset(x = Versch_Herd6, subset = Versch_Herd6$Nematode=="Ostertagia")
Versch_Herd6Coop = subset(x = Versch_Herd6, subset = Versch_Herd6$Nematode=="Cooperia")

# Enter start and end dates:
ss = '2013-06-20'
ee = '2013-11-30'

# Enter site location:
lat = 51.1563506
lon = 4.0054557

# Enter temperature and rainfall data, using the eobs() function defined above:
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
precip = eobspoint(data = 'rr_ens_mean_0.1deg_reg_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)

# Host age at beginning of simulation (in days) (estimated 15 months of age - Verschave thesis)
host_age = 450

# Add events 
# NO TREATMENTS
eventdat = NULL
eventdat

# Estimate weights and use this to estimate dry matter intake and faeces production
# Generate a vector of ages to match the length of the simulation
host.age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1)
# Use this to make a vector of estimated host weight using the liveweight function
lw = if(host == 1) {lw_sheep(age = host.age)} else {lw_cattle(age = host.age)}

# Use the weight to estimate dry matter intake
DMI = dmi(lwt = lw, age = host.age, host_species = host)

# Create a vector of DM values for each day of the simulation
kgDMha = rep(kgDM, length(DMI))

# The amount of faeces produced is not important for the simulation but is used in posthoc estimates of faecal egg counts from the simulation output.
f = faeces(lwt = lw, host_species = host)

#stocking rates
stocking_rate = c(rep(3.7, length(seq.Date(as.Date(ss), as.Date('2013-08-20'), 'day'))), rep(3.7, length(temp)-length(seq.Date(as.Date(ss), as.Date('2013-08-20'), 'day'))))
movement_dates = c(ss, '2013-11-19', ee)  #must include the start and end dates
grazing_plots = c(1,0,0) # must include the start and end locations

#Herd 6 Ostertagia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.3,  Preadult_in_host = 787.9234) ###calculated using the eggs to worms calculation. Used the third timepoint egg count. 

Herd6Oster = gloworm_meta(start = ss, end = ee, lat = lat, 
                          temp = temp, precip = precip, 
                          statevars = initial_values, host = 2, nematode = 3,
                          stocking_rate = stocking_rate, graze = grazing_plots,
                          movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                          lwt = lw, faeces = f, eventdat = eventdat)

Herd6Osterdays <- 1:nrow(Herd6Oster) + 7
plot(Herd6Osterdays, Herd6Oster[,'FEC'], ylim = c(0, 500), xlim = c(-60, length(Versch_Herd6Ost$day) + 183),  xaxt = 'n', main = 'Verschave Herd 6', xlab = 'Month', ylab = expression(paste('FEC ', italic('O. ostertagia'))), type ='l', col = 'grey')
axis(side = 1, at = c(-60, -30,1,31, 61, 91, 122, 153, 183), labels = c('Apr', 'May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd6Ost$day+7), ui = Versch_Herd6Ost$mean+((Versch_Herd6Ost$upperCI-Versch_Herd6Ost$mean)), li = Versch_Herd6Ost$mean-((Versch_Herd6Ost$mean-Versch_Herd6Ost$lowerCI)), add = TRUE, y = Versch_Herd6Ost$mean, pch = 16)


validation_output6ostertagia = data.frame(day = Versch_Herd6Ost$day, 
                                          simulated = Herd6Oster[,'FEC'][Versch_Herd6Ost$day],
                                          observed.mean = Versch_Herd6Ost$mean,
                                          observed.uci = Versch_Herd6Ost$mean-((Versch_Herd6Ost$mean-Versch_Herd6Ost$lowerCI)),
                                          observed.lci = Versch_Herd6Ost$mean+((Versch_Herd6Ost$upperCI-Versch_Herd6Ost$mean)),
                                          species = rep('Ostertagia', length(Versch_Herd6Ost$day)),
                                          study = rep('Verschave', length(Versch_Herd6Ost$day)), 
                                          group = rep("Herd 6", length(Versch_Herd6Ost$day)))

dev.off()
#Herd 6 Cooperia####
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.5,  Preadult_in_host = 24.89051) ###calculated using the eggs to worms calculation. Used the third timepoint egg count. 

Herd6Coop = gloworm_meta(start = ss, end = ee, lat = lat, 
                         temp = temp, precip = precip, 
                         statevars = initial_values, host = 2, nematode = 4,
                         stocking_rate = stocking_rate, graze = grazing_plots,
                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                         lwt = lw, faeces = f, eventdat = eventdat)

Herd6Coopdays <- 1:nrow(Herd6Coop) + 45
plot(Herd6Coopdays, Herd6Coop[,'FEC'], ylim = c(0, 800), xaxt = 'n', xlim = c(1, length(Versch_Herd6Coop$day) + 243), main = 'Verschave Herd 6',  xlab = 'Month', ylab = expression(paste(italic('C. oncophora'), 'FEC')),  type ='l', col = 'grey')
axis(side = 1, at = c(1, 30, 61, 91, 122, 153, 183, 213, 243), labels = c('Apr','May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd6Coop$day+45), ui = Versch_Herd6Coop$mean+((Versch_Herd6Coop$upperCI-Versch_Herd6Coop$mean)), li = Versch_Herd6Coop$mean-((Versch_Herd6Coop$mean-Versch_Herd6Coop$lowerCI)), add = TRUE, y = Versch_Herd6Coop$mean, pch = 16)

validation_output6cooperia = data.frame(day = Versch_Herd6Coop$day, 
                                        simulated = Herd6Coop[,'FEC'][Versch_Herd6Coop$day],
                                        observed.mean = Versch_Herd6Coop$mean,
                                        observed.uci = Versch_Herd6Coop$mean-((Versch_Herd6Coop$mean-Versch_Herd6Coop$lowerCI)),
                                        observed.lci = Versch_Herd6Coop$mean+((Versch_Herd6Coop$upperCI-Versch_Herd6Coop$mean)),
                                        species = rep('Cooperia', length(Versch_Herd6Coop$day)),
                                        study = rep('Verschave', length(Versch_Herd6Coop$day)), 
                                        group = rep("Herd 6", length(Versch_Herd6Coop$day)))


#Herd 7 Set-up####
#Seperate Ostertagia and Cooperia data
Versch_Herd7Ost = subset(x = Versch_Herd7, subset = Versch_Herd7$Nematode=="Ostertagia")
Versch_Herd7Coop = subset(x = Versch_Herd7, subset = Versch_Herd7$Nematode=="Cooperia")

# Enter start and end dates:
ss = "2013-06-13"
ee = "2013-09-10"

# Enter site location:
lat = 51.2054297
lon = 3.5456105

# Enter temperature and rainfall data, using the eobs() function defined above:
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
precip = eobspoint(data = 'rr_ens_mean_0.1deg_reg_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)

# Host age at beginning of simulation (in days) (estimated 15 months of age - Verschave thesis)
host_age = 300


# Add events 
# Moxidectin pour on added 19/08/2013
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(68, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

# Estimate weights and use this to estimate dry matter intake and faeces production
# Generate a vector of ages to match the length of the simulation
host.age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1)
# Use this to make a vector of estimated host weight using the liveweight function
lw = if(host == 1) {lw_sheep(age = host.age)} else {lw_cattle(age = host.age)}

# Use the weight to estimate dry matter intake
DMI = dmi(lwt = lw, age = host.age, host_species = host)

# Create a vector of DM values for each day of the simulation
kgDMha = rep(kgDM, length(DMI))

# The amount of faeces produced is not important for the simulation but is used in posthoc estimates of faecal egg counts from the simulation output.
f = faeces(lwt = lw, host_species = host)

#stocking rate
stocking_rate = c(rep(25.5, length(seq.Date(as.Date(ss), as.Date('2013-09-05'), 'day'))), rep(25.5, length(temp)-length(seq.Date(as.Date(ss), as.Date('2013-09-05'), 'day'))))
movement_dates = c(ss, '2013-09-05', ee)  #must include the start and end dates
grazing_plots = c(1,0,0) # must include the start and end locations

#Herd 7 Ostertagia####
initial_values = init.vals(immunity_host = 0.3,  Preadult_in_host = 775.2025) ###calculated using the eggs to worms calculation. Used the third timepoint egg count. 

Herd7Oster = gloworm_meta(start = ss, end = ee, lat = lat, 
                          temp = temp, precip = precip, 
                          statevars = initial_values, host = 2, nematode = 3,
                          stocking_rate = stocking_rate, graze = grazing_plots,
                          movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                          lwt = lw, faeces = f, eventdat = eventdat)

Herd7Osterdays <- 1:nrow(Herd7Oster) + 13
plot(Herd7Osterdays, Herd7Oster[,'FEC'], ylim = c(0, 500), xlim = c(-60, length(Versch_Herd7Ost$day) + 183),  xaxt = 'n', main = 'Verschave Herd 7', xlab = 'Month', ylab = expression(paste('FEC ', italic('O. ostertagia'))), type ='l', col = 'grey')
axis(side = 1, at = c(-60, -30,1,31, 61, 91, 122, 153, 183), labels = c('Apr', 'May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd7Ost$day+7), ui = Versch_Herd7Ost$mean+((Versch_Herd7Ost$upperCI-Versch_Herd7Ost$mean)), li = Versch_Herd7Ost$mean-((Versch_Herd7Ost$mean-Versch_Herd7Ost$lowerCI)), add = TRUE, y = Versch_Herd7Ost$mean, pch = 16)


validation_output7ostertagia = data.frame(day = Versch_Herd7Ost$day, 
                                          simulated = Herd7Oster[,'FEC'][Versch_Herd7Ost$day],
                                          observed.mean = Versch_Herd7Ost$mean,
                                          observed.uci = Versch_Herd7Ost$mean-((Versch_Herd7Ost$mean-Versch_Herd7Ost$lowerCI)),
                                          observed.lci = Versch_Herd7Ost$mean+((Versch_Herd7Ost$upperCI-Versch_Herd7Ost$mean)),
                                          species = rep('Ostertagia', length(Versch_Herd7Ost$day)),
                                          study = rep('Verschave', length(Versch_Herd7Ost$day)), 
                                          group = rep("Herd 7", length(Versch_Herd7Ost$day)))
dev.off()
#Herd 7 Cooperia####

initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.3,  Preadult_in_host = 95.64178) ###calculated using the eggs to worms calculation. Used the third timepoint egg count. 

Herd7Coop = gloworm_meta(start = ss, end = ee, lat = lat, 
                         temp = temp, precip = precip, 
                         statevars = initial_values, host = 2, nematode = 4,
                         stocking_rate = stocking_rate, graze = grazing_plots,
                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                         lwt = lw, faeces = f, eventdat = eventdat)

Herd7Coopdays <- 1:nrow(Herd7Coop) + 45
plot(Herd7Coopdays, Herd7Coop[,'FEC'], ylim = c(0, 800), xaxt = 'n', xlim = c(1, length(Versch_Herd7Coop$day) + 243), main = 'Verschave Herd 7',  xlab = 'Month', ylab = expression(paste(italic('C. oncophora'), 'FEC')),  type ='l', col = 'grey')
axis(side = 1, at = c(1, 30, 61, 91, 122, 153, 183, 213, 243), labels = c('Apr','May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
plotCI(x = (Versch_Herd7Coop$day+45), ui = Versch_Herd7Coop$mean+((Versch_Herd7Coop$upperCI-Versch_Herd7Coop$mean)), li = Versch_Herd7Coop$mean-((Versch_Herd7Coop$mean-Versch_Herd7Coop$lowerCI)), add = TRUE, y = Versch_Herd7Coop$mean, pch = 16)


dev.off()

validation_output7cooperia = data.frame(day = Versch_Herd7Coop$day, 
                                        simulated = Herd7Coop[,'FEC'][Versch_Herd7Coop$day],
                                        observed.mean = Versch_Herd7Coop$mean,
                                        observed.uci = Versch_Herd7Coop$mean-((Versch_Herd7Coop$mean-Versch_Herd7Coop$lowerCI)),
                                        observed.lci = Versch_Herd7Coop$mean+((Versch_Herd7Coop$upperCI-Versch_Herd7Coop$mean)),
                                        species = rep('Cooperia', length(Versch_Herd7Coop$day)),
                                        study = rep('Verschave', length(Versch_Herd7Coop$day)), 
                                        group = rep("Herd 7", length(Versch_Herd7Coop$day)))



##OstertagiaValidation Graphs####

ostertagia.validation <- rbind(validation_output1ostertagia, validation_output2ostertagia, validation_output3ostertagia, validation_output4ostertagia, validation_output5ostertagia , validation_output6ostertagia, validation_output7ostertagia )
osterlm <- lm(ostertagia.validation$observed.mean ~ 0 + ostertagia.validation$simulated)

plot(ostertagia.validation$simulated, ostertagia.validation$observed.mean,
     ylab = 'observed FEC', main = expression(paste(italic('O. ostertagia'), 'Validation ')), xlab = 'simulated FEC', pch = 16)
abline(a = 0, b = 1)  # y = x line
abline(a = 0, b = osterlm$coefficients[1], lty = 2)  # fitted model line

write.table(ostertagia.validation, file = '/Users/oliviaingle/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/BRACE R/Validations/FinalValidationScripts/Verschave_OstertagiaValidation.csv', sep = ',', row.names=FALSE, col.names=TRUE)

par(mfrow=c(3,3))
tiff('Ost_Verschave_longitudinal.tiff', width = 10, height = 10, units = 'in', res = 300)

dev.off()
tiff('Ost_Verschave_sims_obvs.tiff', width = 5, height = 5, units = 'in', res = 300)

dev.off()

##Cooperia Validations Graphs####
cooperia.validation <- rbind(validation_output1cooperia, validation_output2cooperia, validation_output3cooperia, validation_output4cooperia, validation_output5cooperia, validation_output6cooperia, validation_output7cooperia)
cooplm <- lm(cooperia.validation$observed.mean ~ 0 + cooperia.validation$simulated)

plot(cooperia.validation$simulated, cooperia.validation$observed.mean,
     ylab = 'observed FEC', main = expression(paste(italic('C. oncophora'), 'Validation ')), xlab = 'simulated FEC', pch = 16)
abline(a = 0, b = 1)  # y = x line
abline(a = 0, b = cooplm$coefficients[1], lty = 2)  # fitted model line

write.table(cooperia.validation, file = '/Users/oliviaingle/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/R/BRACE R/Validations/FinalValidationScripts/Verschave_CooperiaValidation.csv', sep = ',', row.names=FALSE, col.names=TRUE)

par(mfrow=c(3,3))
tiff('Coop_Verschave_longitudinal.tiff', width = 10, height = 10, units = 'in', res = 300)


dev.off()
tiff('Coop_Verschave_sims_obvs.tiff', width = 5, height = 5, units = 'in', res = 300)

dev.off()



