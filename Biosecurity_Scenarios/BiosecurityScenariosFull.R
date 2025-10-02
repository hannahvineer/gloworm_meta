###Paper 1 Biosecurity Scenarios for Scotland, UK lons/lats 
##Created by: Olivia Ingle olivia.ingle@liverpool.ac.uk - 23/01/2025
##Last Modified: Olivia Ingle 

###Peniciuk Scotland and Cahors, France lon/lats. 
##Ostertagia and Cooperia simulated
#Spring and Autumn start dates simulated


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
source('Model_Functions/livestockfuns.R') # liveweight, faeces production and dry matter intake
source('Model_Functions/weatherfuns.R') # EOBS data extraction function
source('Model_Functions/parainterp.R') # Interpolation function for parasitological input
source('Model_Functions/gloworm_meta.R') # GLOWORM-META function
source('Model_Functions/initialvalues.r') # sets initial conditions
source('Model_Functions/ginparms.r') # species-specific parameters
# Sping SCOTLAND constant model inputs####
# Enter start and end dates:
ss = '2020-05-01'
ee = '2021-04-30'

# Enter site location: ###Peniciuk
lat = 55.8573407
lon = -3.1981503

# Enter temperature and rainfall data, using the eobs() function defined above:
precip = eobspoint('rr_ens_mean_0.1deg_reg_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)

#Cows on hard pasture (0) and then put onto pasture 1 month later
kgDM = 2000

# Host species:
# Choose from:
# 1 = sheep
# 2 = cow
host = 2


#set initial values = adults in host = 10,000. immunity is low as calves 6 months old. 
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.01, Preadult_in_host = 2000, Adult_in_hostA = 8000)

# Host age at beginning of simulation (in days) (estimated 6 months old)
host_age = 180

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

# No Quarantine, No Treatment, animals straight to Pasture 1 ####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-05-08', ee)  #must include the start and end dates
grazing_plots = c(1,1,1)

##Ostertagia

SPRINGNoQ_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 3,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

plot(SPRINGNoQ_OSTER[,'L3p_A'], ylim = c(0, 2000000), xlim = range(month_positions), 
     main = 'No Quarantine', xaxt = 'n', xlab = 'Month', 
     ylab = expression(paste(italic('O. ostertagia'), ' L3p_A')), 
     type = 'l', col = 'grey')
axis(side = 1, at = month_positions, labels = month_labels)
text(230, 1980000, "No Effective Treatment", col = "black")

##Cooperia 

SPRINGNoQ_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                               temp = temp, precip = precip, 
                               statevars = initial_values, host = 2, nematode = 4,
                               stocking_rate = stocking_rate, graze = grazing_plots,
                               movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                               lwt = lw, faeces = f, eventdat = eventdat)

plot(SPRINGNoQ_COOP[,'L3p_A'], ylim = c(0, 2000000), xlim = range(month_positions), 
     main = 'No Quarantine', xaxt = 'n', xlab = 'Month', 
     ylab = expression(paste(italic('O. ostertagia'), ' L3p_A')), 
     type = 'l', col = 'grey')
axis(side = 1, at = month_positions, labels = month_labels)
text(230, 1980000, "No Effective Treatment", col = "black")


# 7 Days Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-05-08', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

SPRING0week_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

SPRING0week_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 99% Effective Treatment####

eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING99week_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99week_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                temp = temp, precip = precip, 
                                statevars = initial_values, host = 2, nematode = 4,
                                stocking_rate = stocking_rate, graze = grazing_plots,
                                movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                lwt = lw, faeces = f, eventdat = eventdat)


# 7 Days Quarantine - 80% Effective Treatment####

eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat


## Ostertagia

SPRING80week_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80week_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING60week_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60week_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)


# 14 Days Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-05-14', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

SPRING0twoweek_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 3,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

SPRING0twoweek_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                temp = temp, precip = precip, 
                                statevars = initial_values, host = 2, nematode = 4,
                                stocking_rate = stocking_rate, graze = grazing_plots,
                                movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING99twoweek_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99twoweek_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING80twoweek_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80twoweek_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING60twoweek_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60twoweek_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-06-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

## Ostertagia

SPRING0month_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING0month_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING99month_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99month_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING80month_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80month_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING60month_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60month_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-07-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

## Ostertagia

SPRING0twomonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING0twomonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING99twomonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99twomonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)


# 2 Months Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING80twomonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80twomonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING60twomonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60twomonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 3 Months Quarantine - No Effective Treatment####

eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-07-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-07-01'), 'day'))))
movement_dates = c(ss, '2020-08-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

SPRING0threemonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 3,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

SPRING0threemonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                temp = temp, precip = precip, 
                                statevars = initial_values, host = 2, nematode = 4,
                                stocking_rate = stocking_rate, graze = grazing_plots,
                                movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                lwt = lw, faeces = f, eventdat = eventdat)

# 3 Months Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

SPRING99threemonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99threemonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)








# 3 Months Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

SPRING80threemonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80threemonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)
# 3 Months Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

SPRING60threemonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60threemonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)
# Winter date constant inputs####
# Enter start and end dates:
ss = '2020-10-01'
ee = '2021-09-30'

#set initial values = adults in host = 10,000. immunity is low as calves 6 months old. 
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.01, Preadult_in_host = 4000, Adult_in_hostA = 6000)

# No Quarantine, No Treatment, incoming animals straight to pasture ####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-10-08', ee)  #must include the start and end dates
grazing_plots = c(1,1,1)

##Ostertagia

WINTERNoQ_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                               temp = temp, precip = precip, 
                               statevars = initial_values, host = 2, nematode = 3,
                               stocking_rate = stocking_rate, graze = grazing_plots,
                               movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                               lwt = lw, faeces = f, eventdat = eventdat)

##Cooperia 

WINTERNoQ_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                              temp = temp, precip = precip, 
                              statevars = initial_values, host = 2, nematode = 4,
                              stocking_rate = stocking_rate, graze = grazing_plots,
                              movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                              lwt = lw, faeces = f, eventdat = eventdat)




# 7 Days Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-10-08', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

WINTER0week_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 3,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

WINTER0week_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                temp = temp, precip = precip, 
                                statevars = initial_values, host = 2, nematode = 4,
                                stocking_rate = stocking_rate, graze = grazing_plots,
                                movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER99week_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99week_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat


## Ostertagia

WINTER80week_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80week_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER60week_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60week_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-10-14', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

WINTER0twoweek_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 3,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

WINTER0twoweek_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 4,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER99twoweek_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99twoweek_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)






# 14 Days Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER80twoweek_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80twoweek_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER60twoweek_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60twoweek_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-11-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

## Ostertagia

WINTER0month_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER0month_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER99month_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99month_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER80month_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80month_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER60month_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60month_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-12-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

## Ostertagia

WINTER0twomonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 3,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER0twomonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 4,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER99twomonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99twomonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)


# 2 Months Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER80twomonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80twomonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER60twomonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60twomonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)


















# 3 Months Quarantine - No Effective Treatment####

eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2021-01-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2021-01-01'), 'day'))))
movement_dates = c(ss, '2021-01-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

WINTER0threemonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 3,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)

WINTER0threemonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 4,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

# 3 Months Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

WINTER99threemonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99threemonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)








# 3 Months Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

WINTER80threemonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80threemonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)
# 3 Months Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

WINTER60threemonths_OSTER = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60threemonths_COOP = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)




















###SPRING EUROPEAN SCENARIOS####
# Enter start and end dates:
ss = '2020-05-01'
ee = '2021-04-30'

# Enter site location: ###Cahors
lat = 44.4
lon = 1.6

# Enter temperature and rainfall data, using the eobs() function defined above:
precip = eobspoint('rr_ens_mean_0.1deg_reg_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)
temp = eobspoint(data = 'tg_ens_mean_0.1deg_reg_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)

#Cows on hard pasture (0) and then put onto pasture 1 month later
kgDM = 2000

# Host species:
# Choose from:
# 1 = sheep
# 2 = cow
host = 2


#set initial values = adults in host = 10,000. immunity is low as calves 6 months old. 
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.01, Preadult_in_host = 2000, Adult_in_hostA = 8000)

# Host age at beginning of simulation (in days) (estimated 6 months old)
host_age = 180

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
# No Quarantine, No Treatment, animals straight to Pasture 1 ####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-05-08', ee)  #must include the start and end dates
grazing_plots = c(1,1,1)

##Ostertagia

SPRINGNoQ_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                               temp = temp, precip = precip, 
                               statevars = initial_values, host = 2, nematode = 3,
                               stocking_rate = stocking_rate, graze = grazing_plots,
                               movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                               lwt = lw, faeces = f, eventdat = eventdat)

##Cooperia 

SPRINGNoQ_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                              temp = temp, precip = precip, 
                              statevars = initial_values, host = 2, nematode = 4,
                              stocking_rate = stocking_rate, graze = grazing_plots,
                              movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                              lwt = lw, faeces = f, eventdat = eventdat)


# 7 Days Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-05-08', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

SPRING0week_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 3,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

SPRING0week_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                temp = temp, precip = precip, 
                                statevars = initial_values, host = 2, nematode = 4,
                                stocking_rate = stocking_rate, graze = grazing_plots,
                                movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 99% Effective Treatment####

eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING99week_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99week_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)


# 7 Days Quarantine - 80% Effective Treatment####

eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat


## Ostertagia

SPRING80week_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80week_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING60week_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60week_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)


# 14 Days Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-05-14', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

SPRING0twoweek_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 3,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

SPRING0twoweek_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 4,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING99twoweek_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99twoweek_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING80twoweek_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80twoweek_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING60twoweek_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60twoweek_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-06-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

## Ostertagia

SPRING0month_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING0month_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING99month_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99month_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING80month_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80month_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING60month_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60month_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-06-01'), 'day'))))
movement_dates = c(ss, '2020-07-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

## Ostertagia

SPRING0twomonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 3,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING0twomonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 4,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING99twomonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99twomonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)


# 2 Months Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING80twomonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80twomonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

SPRING60twomonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60twomonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)

# 3 Months Quarantine - No Effective Treatment####

eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-07-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-07-01'), 'day'))))
movement_dates = c(ss, '2020-08-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

SPRING0threemonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 3,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)

SPRING0threemonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 4,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

# 3 Months Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

SPRING99threemonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING99threemonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)

# 3 Months Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

SPRING80threemonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING80threemonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)
# 3 Months Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

SPRING60threemonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

SPRING60threemonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)



##WINTER EUROPEAN SCENARIOS####
ss = '2020-10-01'
ee = '2021-09-30'

#set initial values = adults in host = 10,000. immunity is low as calves 6 months old. 
initial_values = init.vals()
initial_values = init.vals(immunity_host = 0.01, Preadult_in_host = 4000, Adult_in_hostA = 6000)

# No Quarantine, No Treatment, incoming animals straight to pasture ####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-10-08', ee)  #must include the start and end dates
grazing_plots = c(1,1,1)

##Ostertagia

WINTERNoQ_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                               temp = temp, precip = precip, 
                               statevars = initial_values, host = 2, nematode = 3,
                               stocking_rate = stocking_rate, graze = grazing_plots,
                               movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                               lwt = lw, faeces = f, eventdat = eventdat)

##Cooperia 

WINTERNoQ_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                              temp = temp, precip = precip, 
                              statevars = initial_values, host = 2, nematode = 4,
                              stocking_rate = stocking_rate, graze = grazing_plots,
                              movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                              lwt = lw, faeces = f, eventdat = eventdat)




# 7 Days Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-10-08', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

WINTER0week_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 3,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

WINTER0week_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                temp = temp, precip = precip, 
                                statevars = initial_values, host = 2, nematode = 4,
                                stocking_rate = stocking_rate, graze = grazing_plots,
                                movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER99week_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99week_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat


## Ostertagia

WINTER80week_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80week_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 7 Days Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER60week_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60week_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-10-14', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

WINTER0twoweek_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 3,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

WINTER0twoweek_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 4,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER99twoweek_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99twoweek_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER80twoweek_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80twoweek_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 14 Days Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER60twoweek_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 3,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60twoweek_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                    temp = temp, precip = precip, 
                                    statevars = initial_values, host = 2, nematode = 4,
                                    stocking_rate = stocking_rate, graze = grazing_plots,
                                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                    lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-11-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

## Ostertagia

WINTER0month_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 3,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER0month_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = 2, nematode = 4,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER99month_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99month_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER80month_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80month_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 1 Month Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER60month_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                   temp = temp, precip = precip, 
                                   statevars = initial_values, host = 2, nematode = 3,
                                   stocking_rate = stocking_rate, graze = grazing_plots,
                                   movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                   lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60month_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                  temp = temp, precip = precip, 
                                  statevars = initial_values, host = 2, nematode = 4,
                                  stocking_rate = stocking_rate, graze = grazing_plots,
                                  movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                  lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - No Effective Treatment####
eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2020-10-01'), 'day'))))
movement_dates = c(ss, '2020-12-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

## Ostertagia

WINTER0twomonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 3,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER0twomonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                     temp = temp, precip = precip, 
                                     statevars = initial_values, host = 2, nematode = 4,
                                     stocking_rate = stocking_rate, graze = grazing_plots,
                                     movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                     lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER99twomonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99twomonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)


# 2 Months Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER80twomonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80twomonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)

# 2 Months Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat

## Ostertagia

WINTER60twomonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 3,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60twomonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                      temp = temp, precip = precip, 
                                      statevars = initial_values, host = 2, nematode = 4,
                                      stocking_rate = stocking_rate, graze = grazing_plots,
                                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                      lwt = lw, faeces = f, eventdat = eventdat)

# 3 Months Quarantine - No Effective Treatment####

eventdat = NULL
eventdat

stocking_rate = c(rep(1, length(seq.Date(as.Date(ss), as.Date('2021-01-01'), 'day'))), rep(1, length(temp)-length(seq.Date(as.Date(ss), as.Date('2021-01-01'), 'day'))))
movement_dates = c(ss, '2021-01-01', ee)  #must include the start and end dates
grazing_plots = c(1,2,2)

WINTER0threemonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 3,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)

WINTER0threemonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                       temp = temp, precip = precip, 
                                       statevars = initial_values, host = 2, nematode = 4,
                                       stocking_rate = stocking_rate, graze = grazing_plots,
                                       movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                       lwt = lw, faeces = f, eventdat = eventdat)

# 3 Months Quarantine - 99% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.01, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

WINTER99threemonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER99threemonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)

# 3 Months Quarantine - 80% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.2, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

WINTER80threemonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER80threemonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)
# 3 Months Quarantine - 60% Effective Treatment####
eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(1, 3), 
  value = rep(0.4, 3), 
  method = rep('mult', 3))
eventdat
## Ostertagia

WINTER60threemonths_OSTER_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                         temp = temp, precip = precip, 
                                         statevars = initial_values, host = 2, nematode = 3,
                                         stocking_rate = stocking_rate, graze = grazing_plots,
                                         movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                         lwt = lw, faeces = f, eventdat = eventdat)

## Cooperia

WINTER60threemonths_COOP_EU = gloworm_meta(start = ss, end = ee, lat = lat, 
                                        temp = temp, precip = precip, 
                                        statevars = initial_values, host = 2, nematode = 4,
                                        stocking_rate = stocking_rate, graze = grazing_plots,
                                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                        lwt = lw, faeces = f, eventdat = eventdat)




















##Graphs####
#These can be changed to which sequence is needed. For these we selected the 7 day vs 30 day plots. 
dev.off()

# Set up the plotting area for a 2x2 grid (4 plots total)
par(mfrow = c(2,2), 
    mar = c(0,0,0,0),
    oma = c(4,4,1,1))

#Ostertagia Spring

month_positions <- c(1, 31, 62, 93, 124, 154, 185, 215, 246, 277, 305, 336, 366)
month_labels <- c('May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')

# Set up the plot area with 14 Days data for EU
plot(SPRING99week_OSTER_EU[,'L3p_A'], ylim = c(0, 1500000), xlim = range(month_positions), 
     main = NULL, xaxt = 'n', xlab = "", 
     ylab = expression(paste(italic('O. ostertagia'), ' L3p_A')), 
     type = 'l', lty = 1, col = "#FFB1A8", lwd = 2) # Solid dark orange for 14 Days 99% Effective

# Add quarantine end lines for EU
abline(v = 8, col = "black", lwd = 2)
text(x = 10,  y = max(c(600000, 680000, 98000)), labels = "7 Days Quarantine End", cex = 1.2, srt = 90, pos = 4, col = "black")
abline(v = 31, col = "black", lwd = 2)
text(x = 33,  y = max(c(600000, 680000, 98000)), labels = "30 Days Quarantine End", cex = 1.2, srt = 90, pos = 4, col = "black")
text(x = 33,  y = 1450000, labels = expression(paste(italic('O. ostertagia'), ' Spring')), cex = 1.4, pos = 4, col = "black")


# Add lines for EU treatments
lines(SPRINGNoQ_OSTER_EU[,'L3p_A'], lty = 5, col = "#9F4A3E", lwd = 3) # Dot-dash for No Treatments
lines(SPRING99week_OSTER_EU[,'L3p_A'], lty = 1, col = "#FFB1A8", lwd = 3) # Solid light orange for 7 Days 99% Effective
lines(SPRING60week_OSTER_EU[,'L3p_A'], lty = 2, col = "#FFB1A8", lwd = 3) # Dashed for 60%
lines(SPRING80week_OSTER_EU[,'L3p_A'], lty = 3, col = "#FFB1A8", lwd = 3) # Dotted for 80%

lines(SPRING99month_OSTER_EU[,'L3p_A'], lty = 1, col = "#C65E57", lwd = 3) # Solid dark orange for 7 Days 99% Effective
lines(SPRING60month_OSTER_EU[,'L3p_A'], lty = 2, col = "#C65E57", lwd = 3) # Dashed for 60%
lines(SPRING80month_OSTER_EU[,'L3p_A'], lty = 3, col = "#C65E57", lwd = 3) # Dotted for 80%

# Add lines for UK treatments
lines(SPRINGNoQ_OSTER[,'L3p_A'], lty = 5, col = "#006C75", lwd = 3) # Dot-dash for No Treatments
lines(SPRING99week_OSTER[,'L3p_A'], lty = 1, col = "#A9ECEE", lwd = 3) # Solid light cyan for 7 Days 99% Effective
lines(SPRING60week_OSTER[,'L3p_A'], lty = 2, col = "#A9ECEE", lwd = 3) # Dashed for 60%
lines(SPRING80week_OSTER[,'L3p_A'], lty = 3, col = "#A9ECEE", lwd = 3) # Dotted for 80%

lines(SPRING99month_OSTER[,'L3p_A'], lty = 1, col = "#008A91", lwd = 3) # Solid dark cyan for 7 Days 99% Effective
lines(SPRING60month_OSTER[,'L3p_A'], lty = 2, col = "#008A91", lwd = 3) # Dashed for 60%
lines(SPRING80month_OSTER[,'L3p_A'], lty = 3, col = "#008A91", lwd = 3) # Dotted for 80%

#Ostertagia Winter

month_positions <- c(1, 31, 62, 93, 124, 154, 185, 215, 246, 277, 305, 336, 366)
month_labels <- c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct')

# Set up the plot area with 14 Days data for EU
plot(WINTER99week_OSTER_EU[,'L3p_A'], ylim = c(0, 1500000), xlim = range(month_positions), 
     main = NULL, yaxt = 'n', xaxt = 'n', xlab = "", ylab = "",
     type = 'l', lty = 1, col = "#FFB1A8", lwd = 2) # Solid dark orange for 14 Days 99% Effective

# Add quarantine end lines for EU
abline(v = 8, col = "black", lwd = 2)
text(x = 10,  y = max(c(600000, 680000, 98000)), labels = "7 Days Quarantine End", cex = 1.2, srt = 90, pos = 4, col = "black")
abline(v = 31, col = "black", lwd = 2)
text(x = 33,  y = max(c(600000, 680000, 98000)), labels = "30 Days Quarantine End", cex = 1.2, srt = 90, pos = 4, col = "black")
text(x = 33,  y = 1450000, labels = expression(paste(italic('O. ostertagia'), ' Autumn')), cex = 1.4, pos = 4, col = "black")


# Add lines for EU treatments
lines(WINTERNoQ_OSTER_EU[,'L3p_A'], lty = 5, col = "#9F4A3E", lwd = 3) # Dot-dash for No Treatments
lines(WINTER99week_OSTER_EU[,'L3p_A'], lty = 1, col = "#FFB1A8", lwd = 3) # Solid light orange for 7 Days 99% Effective
lines(WINTER60week_OSTER_EU[,'L3p_A'], lty = 2, col = "#FFB1A8", lwd = 3) # Dashed for 60%
lines(WINTER80week_OSTER_EU[,'L3p_A'], lty = 3, col = "#FFB1A8", lwd = 3) # Dotted for 80%

lines(WINTER99month_OSTER_EU[,'L3p_A'], lty = 1, col = "#C65E57", lwd = 3) # Solid dark orange for 7 Days 99% Effective
lines(WINTER60month_OSTER_EU[,'L3p_A'], lty = 2, col = "#C65E57", lwd = 3) # Dashed for 60%
lines(WINTER80month_OSTER_EU[,'L3p_A'], lty = 3, col = "#C65E57", lwd = 3) # Dotted for 80%

# Add lines for UK treatments
lines(WINTERNoQ_OSTER[,'L3p_A'], lty = 5, col = "#006C75", lwd = 3) # Dot-dash for No Treatments
lines(WINTER99week_OSTER[,'L3p_A'], lty = 1, col = "#A9ECEE", lwd = 3) # Solid light cyan for 7 Days 99% Effective
lines(WINTER60week_OSTER[,'L3p_A'], lty = 2, col = "#A9ECEE", lwd = 3) # Dashed for 60%
lines(WINTER80week_OSTER[,'L3p_A'], lty = 3, col = "#A9ECEE", lwd = 3) # Dotted for 80%

lines(WINTER99month_OSTER[,'L3p_A'], lty = 1, col = "#008A91", lwd = 3) # Solid dark cyan for 7 Days 99% Effective
lines(WINTER60month_OSTER[,'L3p_A'], lty = 2, col = "#008A91", lwd = 3) # Dashed for 60%
lines(WINTER80month_OSTER[,'L3p_A'], lty = 3, col = "#008A91", lwd = 3) # Dotted for 80%
legend("topright", 
       legend = c("No Biosecurity Measures (UK)", "No Biosecurity Measures (EU)", 
                  "60% Effective", "80% Effective", "99% Effective", 
                  "7 Days (UK)", "30 Days (UK)", "7 Days (EU)", "30 Days (EU)"), 
       lty = c(5, 5, 2, 3, 1, NA, NA, NA, NA),  
       col = c("#006C75", "#9F4A3E", "black", "black", "black", 
               "#FFB1A8", "#C65E57", "#A9ECEE", "#008A91"),  
       lwd = c(2, 2, 2, 2, 2, NA, NA, NA, NA),  
       pch = c(NA, NA, NA, NA, NA, 15, 15, 15, 15),
       y.intersp = 0.6,             # Reduce vertical spacing
       x.intersp = 0.5,             # Reduce space between symbols and text
       adj = 0,                     # Align text more to the left
       pt.cex = c(1, 1, 1, 1, 1, 2, 2, 2, 2),  
       cex = 1,  
       bty = "o",                   # Keep the box, remove unnecessary padding
       box.lwd = 0.5,                # Make box thinner for a tighter fit
       box.col = "black",            # Ensure the box is visible but not oversized
       inset = c(-0.3, 0))


#COOPERIA Spring

month_positions <- c(1, 31, 62, 93, 124, 154, 185, 215, 246, 277, 305, 336, 366)
month_labels <- c('May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May')

# Set up the plot area with 14 Days data for EU
plot(SPRING99week_COOP_EU[,'L3p_A'], ylim = c(0, 1500000), xlim = range(month_positions), 
     main = NULL, xaxt = 'n', xlab = 'Month', 
     ylab = expression(paste(italic('C. oncophora'), ' L3p_A')), 
     type = 'l', lty = 1, col = "#FFB1A8", lwd = 2) # Solid dark orange for 14 Days 99% Effective

# Add axis labels for UK
axis(side = 1, at = month_positions, labels = month_labels)

# Add quarantine end lines for EU
abline(v = 8, col = "black", lwd = 2)
text(x = 10,  y = max(c(600000, 680000, 98000)), labels = "7 Days Quarantine End", cex = 1.2, srt = 90, pos = 4, col = "black")
abline(v = 31, col = "black", lwd = 2)
text(x = 33,  y = max(c(600000, 680000, 98000)), labels = "30 Days Quarantine End", cex = 1.2, srt = 90, pos = 4, col = "black")
text(x = 33,  y = 1450000, labels = expression(paste(italic('C. oncophora'), ' Spring')), cex = 1.4, pos = 4, col = "black")


# Add lines for EU treatments
lines(SPRINGNoQ_COOP_EU[,'L3p_A'], lty = 5, col = "#9F4A3E", lwd = 3) # Dot-dash for No Treatments
lines(SPRING99week_COOP_EU[,'L3p_A'], lty = 1, col = "#FFB1A8", lwd = 3) # Solid light orange for 7 Days 99% Effective
lines(SPRING60week_COOP_EU[,'L3p_A'], lty = 2, col = "#FFB1A8", lwd = 3) # Dashed for 60%
lines(SPRING80week_COOP_EU[,'L3p_A'], lty = 3, col = "#FFB1A8", lwd = 3) # Dotted for 80%

lines(SPRING99month_COOP_EU[,'L3p_A'], lty = 1, col = "#C65E57", lwd = 3) # Solid dark orange for 7 Days 99% Effective
lines(SPRING60month_COOP_EU[,'L3p_A'], lty = 2, col = "#C65E57", lwd = 3) # Dashed for 60%
lines(SPRING80month_COOP_EU[,'L3p_A'], lty = 3, col = "#C65E57", lwd = 3) # Dotted for 80%

# Add lines for UK treatments
lines(SPRINGNoQ_COOP[,'L3p_A'], lty = 5, col = "#006C75", lwd = 3) # Dot-dash for No Treatments
lines(SPRING99week_COOP[,'L3p_A'], lty = 1, col = "#A9ECEE", lwd = 3) # Solid light cyan for 7 Days 99% Effective
lines(SPRING60week_COOP[,'L3p_A'], lty = 2, col = "#A9ECEE", lwd = 3) # Dashed for 60%
lines(SPRING80week_COOP[,'L3p_A'], lty = 3, col = "#A9ECEE", lwd = 3) # Dotted for 80%

lines(SPRING99month_COOP[,'L3p_A'], lty = 1, col = "#008A91", lwd = 3) # Solid dark cyan for 7 Days 99% Effective
lines(SPRING60month_COOP[,'L3p_A'], lty = 2, col = "#008A91", lwd = 3) # Dashed for 60%
lines(SPRING80month_COOP[,'L3p_A'], lty = 3, col = "#008A91", lwd = 3) # Dotted for 80%


#COOPERIA WINTER

month_positions <- c(1, 31, 62, 93, 124, 154, 185, 215, 246, 277, 305, 336, 366)
month_labels <- c('Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct')

# Set up the plot area with 14 Days data for EU
plot(WINTER99week_COOP_EU[,'L3p_A'], ylim = c(0, 1500000), xlim = range(month_positions), 
     main = NULL, xaxt = 'n', xlab = 'Month', yaxt = 'n',
     ylab = "", 
     type = 'l', lty = 1, col = "#FFB1A8", lwd = 2) # Solid dark orange for 14 Days 99% Effective

# Add axis labels for UK
axis(side = 1, at = month_positions, labels = month_labels)

# Add quarantine end lines for EU
abline(v = 8, col = "black", lwd = 2)
text(x = 10,  y = max(c(600000, 680000, 98000)), labels = "7 Days Quarantine End", cex = 1.2, srt = 90, pos = 4, col = "black")
abline(v = 31, col = "black", lwd = 2)
text(x = 33,  y = max(c(600000, 680000, 98000)), labels = "30 Days Quarantine End", cex = 1.2, srt = 90, pos = 4, col = "black")
text(x = 33,  y = 1450000, labels = expression(paste(italic('C. oncophora'), ' Autumn')), cex = 1.4, pos = 4, col = "black")


# Add lines for EU treatments
lines(WINTERNoQ_COOP_EU[,'L3p_A'], lty = 5, col = "#9F4A3E", lwd = 3) # Dot-dash for No Treatments
lines(WINTER99week_COOP_EU[,'L3p_A'], lty = 1, col = "#FFB1A8", lwd = 3) # Solid light orange for 7 Days 99% Effective
lines(WINTER60week_COOP_EU[,'L3p_A'], lty = 2, col = "#FFB1A8", lwd = 3) # Dashed for 60%
lines(WINTER80week_COOP_EU[,'L3p_A'], lty = 3, col = "#FFB1A8", lwd = 3) # Dotted for 80%

lines(WINTER99month_COOP_EU[,'L3p_A'], lty = 1, col = "#C65E57", lwd = 3) # Solid dark orange for 7 Days 99% Effective
lines(WINTER60month_COOP_EU[,'L3p_A'], lty = 2, col = "#C65E57", lwd = 3) # Dashed for 60%
lines(WINTER80month_COOP_EU[,'L3p_A'], lty = 3, col = "#C65E57", lwd = 3) # Dotted for 80%

# Add lines for UK treatments
lines(WINTERNoQ_COOP[,'L3p_A'], lty = 5, col = "#006C75", lwd = 3) # Dot-dash for No Treatments
lines(WINTER99week_COOP[,'L3p_A'], lty = 1, col = "#A9ECEE", lwd = 3) # Solid light cyan for 7 Days 99% Effective
lines(WINTER60week_COOP[,'L3p_A'], lty = 2, col = "#A9ECEE", lwd = 3) # Dashed for 60%
lines(WINTER80week_COOP[,'L3p_A'], lty = 3, col = "#A9ECEE", lwd = 3) # Dotted for 80%

lines(WINTER99month_COOP[,'L3p_A'], lty = 1, col = "#008A91", lwd = 3) # Solid dark cyan for 7 Days 99% Effective
lines(WINTER60month_COOP[,'L3p_A'], lty = 2, col = "#008A91", lwd = 3) # Dashed for 60%
lines(WINTER80month_COOP[,'L3p_A'], lty = 3, col = "#008A91", lwd = 3) # Dotted for 80%


















