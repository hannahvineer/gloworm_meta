# Unit testing GLOWORM-META

# Author: Hannah Vineer - hannah.vineer@liverpool.ac.uk
# Created: 03/08/2024
# Last modified: 03/08/2024
# Last modified by: Hannah Vineer

# Load in packages and required source files -----------------------------------
packages = c("testthat", "deSolve", "geosphere", "chron", "ncdf4", "imputeTS", "forecast", "gplots", "tidyr") # Package names

# Check for, and install missing package
installed_packages = packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE)) # Load packages

# Source model functions -------------------------------------------------------
source('Model_functions/livestockfuns.R') # liveweight, faeces production and dry matter intake
source('Model_functions/weatherfuns.R') # EOBS data extraction function
source('Model_functions/parainterp.R') # Interpolation function for parasitological input
source('Model_functions/gloworm_meta.R') # GLOWORM-META function
source('Model_functions/initialvalues.r') # sets initial conditions
source('Model_functions/ginparms.R') # species-specific parameters

# Start and end dates, and lat/lon for AFBI Hillsborough
ss = "2019-04-19"
ee = "2020-04-19"
lon = -6.078
lat = 54.45

# Enter temperature and rainfall data, using the eobs function sourced above:
temp = read.csv(file = 'Validation_data/McFarland/published_data/Raw temperature data EOBS.csv')[,2]
precip = read.csv(file = 'Validation_data/McFarland/published_data/Raw precipitation data EOBS.csv')[,2]

#Average herbage biomass over the simulation (kilograms dry matter per hectare)
kgDM = 3000

# Host species:
# Choose from:
# 1 = sheep
# 2 = cow
host = 2

# Host age at beginning of simulation (in days)
host_age = 365 # 10-12 month calves. Use upper age limit as these are high performance meat breeds

# Events
#1) One doramectin treatment 78% effective on 28th May, 40 days after heifer turn out. 

eventdat = data.frame(
  var = c('Pa', 'P', 'A'), 
  time = rep(40, 3), 
  value = rep(0.22, 3), 
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

stocking_rate = c(rep(40, length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))), rep(20, length(temp)-length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))))
movement_dates = c(ss, ee)
grazing_plots_reps = data.frame(start = c(1:6), end = c(1:6))
# Ostertagia simulations =======================================================

# Nematode:
# Choose from:
# 1 = H. contortus
# 2 = T. circumcincta
# 3 = O. ostertagi
# 4 = C. oncophora
nematode = 3

# Enter starting parasite population size and level of acquired immunity...
# Need to estimate adult population based on egg counts
# assuming immunity = 0.5 as they are 10-12 months of age...
# FEC ~100 epg (Fig 1. McFarland et al., 2022, states decrease in initial FECs, so assuming downward trajectory since April = ~100 epg)
ginparms = gin(nematode = nematode, temp = temp, precip = precip, photoperiod = daylength(lat = lat, doy = seq(as.Date(ss), as.Date(ee), "days")))
fecundity_init = exp(ginparms$max.lambda-(ginparms$max.lambda-ginparms$min.lambda)*0.25)
total_eggs_per_day_per_head = 100*f[1]
est_adults_from_FEC = total_eggs_per_day_per_head/fecundity_init

# Because we now have the init.vals() function above, we only need to specify state variables that should be non-zero. In this case we have newborn lambs and clean pasture so use the defaults for all state variables.
initial_values = init.vals(immunity_host = 0.25, q = 0.22, Adult_in_hostA = est_adults_from_FEC)


# Run the model once to give the expected output
expected_ost = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate, graze = grazing_plots_reps[1,],
                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)


# Loop through the patches - and check if output is identical.
# Differences would suggest that errors were introduced when code was duplicated

for(i in 2:length(grazing_plots_reps[,1])) {
test_ost = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate, graze = grazing_plots_reps[i,],
                    movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)
print(i)
stopifnot(identical(test_ost[,"FEC"], expected_ost[,"FEC"]))
}

# If the loop above stops, you should be able to identify which patch threw 
# the error based on the iteration number.