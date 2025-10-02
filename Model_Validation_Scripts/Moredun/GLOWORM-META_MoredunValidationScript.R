##Author: Olivia Ingle - olivia.ingle@liverpool.ac.uk // Hannah Vineer - hannah.vineer@liverpool.ac.uk
# Created 13/08/2024

##GLOWORM-META Vallidation Script for Ostertagia ostertagi and Cooperia oncophora for the Verschave thesis data for the Belgium herds. 
#Moredun data uses the same lon/lats, eventdat, host age and ss/ee dates. 
#Requires different stocking rates and kgdm. 


##############Source Files/Packages
# Load in packages and required source files -----------------------------------
packages = c("deSolve", "geosphere", "chron", "ncdf4", "imputeTS", "forecast", "gplots", "tidyr", "rstudioapi") # Package names


# Check for, and install missing package
installed_packages = packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE)) # Load packages

setwd(dirname(getActiveDocumentContext()$path)) ###saves outputs to the folder where your Rscript is located
getwd()

# Source model functions -------------------------------------------------------
source('Model_functions/livestockfuns.R') # liveweight, faeces production and dry matter intake
source('Model_functions/weatherfuns.R') # EOBS data extraction function
source('Model_functions/parainterp.R') # Interpolation function for parasitological input
source('Model_functions/gloworm_meta.R') # GLOWORM-META function
source('Model_functions/initialvalues.r') # sets initial conditions
source('Model_functions/ginparms.r') # species-specific parameters
library(dplyr)
library(tidyr)
# Read in validation data and set up basic code####
Moredun <- read.csv("ValidationData/MoredunData.csv")
Moredun_Ostertagia = subset(x = Moredun, subset = Moredun$Nematode=="Ostertagia")
Moredun_Cooperia = subset(x = Moredun, subset = Moredun$Nematode=="Cooperia")

##seperate Ostertagia data into Oster200, Oster500 for graphing later. 
Moredun_Ostertagia200 <- Moredun_Ostertagia %>%
  filter(CowID == "Oster200")

Moredun_Ostertagia500 <- Moredun_Ostertagia %>%
  filter(CowID == "Oster500")

##seperate Cooperia data into Coop300, Coop400 for graphing later. 
Moredun_Cooperia300 <- Moredun_Cooperia %>%
  filter(CowID == "Coop300")

Moredun_Cooperia400 <- Moredun_Cooperia %>%
  filter(CowID == "Coop400")


# Host species:
# Choose from:
# 1 = sheep
# 2 = cow
host = 2

# Enter start and end dates:
ss = '2023-05-24'
ee = '2023-09-30'

# Enter site location: ###Moredun Research Institute
lat = 55.8573407
lon = -3.1981503

# Enter temperature and rainfall data, using the eobs() function defined above:
pre = read.csv("/Validation_data/Moredun/moredunprecipmaysept.csv")
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
precip = pre$mm


# Host age at beginning of simulation (in days) (estimated 1 year old **unsure on exact age of cows**)
host_age = 365

# Add events 
# No events recorded for these cows. 
eventdat = NULL

eventdat


#Ostertagia####

#Stocking Rate
##Ostertagia field = 2263.1 m2
##=0.2263 hectares
##Stocking rate is cows/hectare (2/0.2263 = 8.84 )
# Host management dates and rates

stocking_rate = c(rep(9, length(seq.Date(as.Date(ss), as.Date('2023-08-25'), 'day'))), rep(9, length(temp)-length(seq.Date(as.Date(ss), as.Date('2023-08-25'), 'day'))))
movement_dates = c(ss, '2023-08-25', ee)  #must include the start and end dates
grazing_plots = c(1,1,1)


# Enter starting parasite population size and level of acquired immunity...
# Because we now have the init.vals() function above, we only need to specify state variables that should be non-zero. 
###Ostertagia immunity does not develop in the same way as Cooperia immunity, so lower in this case. See Ploeger et al, 1995. 
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.3, Adult_in_hostA = 20000)

#Average herbage biomass over the simulation (kilograms dry matter per hectare)
kgDM = 3852

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


ostertagia = gloworm_meta(start = ss, end = ee, lat = lat, 
                         temp = temp, precip = precip, 
                         statevars = initial_values, host = 2, nematode = 3,
                         stocking_rate = stocking_rate, graze = grazing_plots,
                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                         lwt = lw, faeces = f, eventdat = eventdat)

tiff('Oster_Moredun_longitudinal.tiff', width = 8, height = 8, units = 'in', res = 300)

plot(ostertagia[,'FEC'], ylim = c(0, 500), xlim = c(-60, length(Moredun_Ostertagia$day) + 182), main = 'Moredun Ostertagia', xaxt = 'n', xlab = 'Month', ylab = expression(paste('FEC ', italic('O. ostertagia'))), type ='l', col = 'grey')
axis(side = 1, at = c(-60, -30,1,31, 61, 91, 122, 153, 183), labels = c('Apr', 'May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
points(Moredun_Ostertagia$day, Moredun_Ostertagia$mean, pch = 16)

dev.off()

Moredun_Ostertagia$simulated = ostertagia[Moredun_Ostertagia$day, 'FEC']

##seperate Ostertagia data into Oster200, Oster500 for graphing later. 
Moredun_Ostertagia200 <- Moredun_Ostertagia %>%
  filter(CowID == "Oster200")

Moredun_Ostertagia500 <- Moredun_Ostertagia %>%
  filter(CowID == "Oster500")

##save as validations

ostertagia.validation200 = data.frame(day = Moredun_Ostertagia200$day, 
                                simulated = Moredun_Ostertagia200$simulated, 
                                observed = Moredun_Ostertagia200$mean,
                                observed.uci = 'NULL',
                                observed.lci = 'NULL',
                                species = rep('Ostertagia', length(Moredun_Ostertagia200$day)),
                                study = rep('Moredun', length(Moredun_Ostertagia200$day)), 
                                group = rep("Oster200", length(Moredun_Ostertagia200$day)))

ostertagia.validation500 = data.frame(day = Moredun_Ostertagia500$day, 
                                   simulated = Moredun_Ostertagia500$simulated, 
                                   observed = Moredun_Ostertagia500$mean,
                                   observed.uci = 'NULL',
                                   observed.lci = 'NULL',
                                   species = rep('Ostertagia', length(Moredun_Ostertagia500$day)),
                                   study = rep('Moredun', length(Moredun_Ostertagia500$day)), 
                                   group = rep("Oster500", length(Moredun_Ostertagia500$day)))

validation_outputoster <- rbind(ostertagia.validation200, ostertagia.validation500)

lm_validation1 = lm(validation_outputoster$observed ~ 0 + validation_outputoster$simulated)
summary(lm_validation1)

tiff('Oster_Moredun_obvs.tiff', width = 8, height = 8, units = 'in', res = 300)

plot(validation_outputoster$simulated, validation_outputoster$observed, xlim = c(0, 300), ylim = c(0, 300),
     ylab = 'observed FEC', main = 'MoredunOstertagia', xlab = 'simulated FEC', pch = 16)
abline(a = 0, b = 1)  # y = x line
abline(a = 0, b = lm_validation1$coefficients[1], lty = 2)  # fitted model line

dev.off()

#Cooperia####
##Cooperia field = 2231.6 m2
##=0.551440369 acres / 0.22316 hectares
##Stocking rate is cows/hectare (2/0.22316 = 9 )
# Host management dates and rates
stocking_rate = c(rep(2, length(seq.Date(as.Date(ss), as.Date('2023-08-25'), 'day'))), rep(2, length(temp)-length(seq.Date(as.Date(ss), as.Date('2023-08-25'), 'day'))))
movement_dates = c(ss, '2023-08-25', ee)  #must include the start and end dates
grazing_plots = c(1,1,1) # must include the start and end locations

initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.3, Preadult_in_hostA = 20000)

#Average herbage biomass over the simulation (kilograms dry matter per hectare)
kgDM = 3054

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

cooperia = gloworm_meta(start = ss, end = ee, lat = lat, 
                          temp = temp, precip = precip, 
                          statevars = initial_values, host = 2, nematode = 4,
                          stocking_rate = stocking_rate, graze = grazing_plots,
                          movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                          lwt = lw, faeces = f, eventdat = eventdat)

Moredun_Cooperia$simulated = cooperia[Moredun_Cooperia$day, 'FEC']

##seperate Ostertagia data into Oster200, Oster500 for graphing later. 
Moredun_Cooperia300 <- Moredun_Cooperia %>%
  filter(CowID == "Coop300")

Moredun_Cooperia400 <- Moredun_Cooperia %>%
  filter(CowID == "Coop400")

tiff('Coop_Moredun_longitudinal.tiff', width = 8, height = 8, units = 'in', res = 300)
#Plot out validation
plot(cooperia[,'FEC'], ylim = c(0, 500), xlim = c(-60, length(Moredun_Cooperia$day) + 182), main = 'Moredun Cooperia', xaxt = 'n', xlab = 'Month', ylab = expression(paste('FEC ', italic('C. oncophora'))), type ='l', col = 'grey')
axis(side = 1, at = c(-60, -30,1,31, 61, 91, 122, 153, 183), labels = c('Apr', 'May','Jun', 'Jul', 'Aug', 'Sept','Oct', 'Nov', 'Dec'))
points(Moredun_Cooperia300$day, Moredun_Cooperia300$mean, pch = 16, col = 'black')
points(Moredun_Cooperia400$day, Moredun_Cooperia400$mean, pch = 16, col = 'black')


dev.off()


##save as validations

cooperia.validation300 = data.frame(day = Moredun_Cooperia300 $day, 
                                      simulated = Moredun_Cooperia300 $simulated, 
                                      observed = Moredun_Cooperia300 $mean,
                                      observed.uci = 'NULL',
                                      observed.lci = 'NULL',
                                      species = rep('Cooperia', length(Moredun_Cooperia300 $day)),
                                      study = rep('Moredun', length(Moredun_Cooperia300 $day)), 
                                      group = rep("Coop300", length(Moredun_Cooperia300 $day)))

cooperia.validation400 = data.frame(day = Moredun_Cooperia400 $day, 
                                      simulated = Moredun_Cooperia400 $simulated, 
                                      observed = Moredun_Cooperia400 $mean,
                                      observed.uci = 'NULL',
                                      observed.lci = 'NULL',
                                      species = rep('Cooperia', length(Moredun_Cooperia400$day)),
                                      study = rep('Moredun', length(Moredun_Cooperia400$day)), 
                                      group = rep("coop400", length(Moredun_Cooperia400$day)))

validation_outputcoop <- rbind(cooperia.validation300, cooperia.validation400)

lm_validation2 = lm(validation_outputcoop$observed ~ 0 + validation_outputcoop$simulated)
summary(lm_validation2)

tiff('Coop_Moredun_obvs.tiff', width = 8, height = 8, units = 'in', res = 300)

plot(validation_outputcoop$simulated, validation_outputcoop$observed, xlim = c(0, 300), ylim = c(0, 300),
     ylab = 'observed FEC', main = 'Moredun Cooperia', xlab = 'simulated FEC', pch = 16)
abline(a = 0, b = 1)  # y = x line
abline(a = 0, b = lm_validation2$coefficients[1], lty = 2)  # fitted model line

dev.off()
