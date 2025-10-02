# Biosecurity simulations - UK
# Uses two sampling matrices generated using the the latin hypercube samples (see "LHS_scenario..." file)

# Load in packages and required source files -----------------------------------
packages = c("dplyr", "rnaturalearth", "sf", "ggplot2", "sensitivity", "deSolve", "geosphere", "chron", "ncdf4", "imputeTS", "forecast", "gplots", "tidyr") # Package names

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


# Load data and generate 'sa' object -------------------------------------------

n = 1000
X1 = data.frame(matrix(runif(5 * n), nrow = n))
colnames(X1) = c("quar_t", "quar_loc", "purch_doy", "fecs", "lon_lat")
X2 = data.frame(matrix(runif(5 * n), nrow = n))
colnames(X2) = c("quar_t", "quar_loc", "purch_doy", "fecs", "lon_lat")

# Generate the "sa" class object using the Sobol function in sensitivity
# Use this after running the simulations
sa = soboljansen(model=NULL, X1 = X1, X2 = X2, nboot=100)
summary(sa)
saveRDS(sa, file = "Biosecurity_simulations/sa.RDS")
scenarios = as.data.frame(sa$X)

# Quarantine time
quar_t = c(1, 84) # days = 0-12 weeks
# Quarantine location
quar_loc = c(0,1)
# Day of year of purchase
purch_doy = c(1, 365)
# Initial FECs
fecs = c(0, 500)
# Lat/lon
UK = ne_countries(country = "United Kingdom", returnclass = "sf", scale = 'medium')
plot(UK)
locations = st_sample(x = UK, size = 5000)
plot(st_geometry(UK))
plot(locations, add = T)
points_df = st_coordinates(locations) %>% as.data.frame()
# Check which locations have available EOBS data
points_df$EOBS_available = FALSE
for (i in 1:5000){
  print(i)
  # Start and end dates, lat and lon
  ss = "2022-01-01"
  ee = "2022-01-10"
  lon = points_df$X[i]
  lat = points_df$Y[i]
  temp = eobspoint(data = 'Biosecurity_simulations/tg_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
  if(!anyNA(temp)) {points_df$EOBS_available[i] = TRUE}
}
points_df = points_df[which(points_df$EOBS_available==TRUE),]
length(which(points_df$EOBS_available==TRUE))
write.csv(points_df, file = "Biosecurity_simulations/random_locations_UK_5000.csv")

# Rescale the sa parameters for use in the simulations
rescale_scenarios <- function(x, newMax, newMin){(x - min(x))/(max(x)-min(x)) * (newMax - newMin) + newMin}
scenarios[,1] = round(rescale_scenarios(x = scenarios[,1], newMin = quar_t[1], newMax = quar_t[2]))
scenarios[,2] = round(scenarios[,2])
scenarios[,3] = round(rescale_scenarios(x = scenarios[,3], newMin = purch_doy[1], newMax = purch_doy[2]))
scenarios[,4] = round(rescale_scenarios(x = scenarios[,4], newMin = fecs[1], newMax = fecs[2]))
scenarios[,5] = round(rescale_scenarios(x = scenarios[,5], newMin = 1, newMax = length(points_df[,1])))
colnames(scenarios) = c("quar_t", "quar_loc", "purch_doy", "fecs", "lon_lat")

# Get the locations and select based on the indices in the scenarios object
locs = read.csv("Biosecurity_simulations/random_locations_UK_5000.csv")
latitudes = locs$Y[scenarios$lon_lat]
longitudes = locs$X[scenarios$lon_lat]

# Source model functions -------------------------------------------------------
source('Model_functions/livestockfuns.R') # liveweight, faeces production and dry matter intake
source('Model_functions/weatherfuns.R') # EOBS data extraction function
source('Model_functions/parainterp.R') # Interpolation function for parasitological input
source('Model_functions/gloworm_meta.R') # GLOWORM-META function
source('Model_functions/initialvalues.r') # sets initial conditions
source('Model_functions/ginparms.R') # species-specific parameters

# Common environmental input ---------------------------------------------------
sim_duration = 182 # days - running the model for 6 months

#Average herbage biomass over the simulation (kilograms dry matter per hectare)
kgDM = 2000

# Host species:
# Choose from:
# 1 = sheep
# 2 = cow
host = 2

# Host age at beginning of simulation (in days)
host_age = 180

# # Events
# #1) One doramectin treatment 78% effective on 28th May, 40 days after heifer turn out. 
# 
# eventdat = data.frame(
#   var = c('Pa', 'P', 'A'), 
#   time = rep(40, 3), 
#   value = rep(0.22, 3), 
#   method = rep('mult', 3))
# 
# eventdat
eventdat = NULL

# Ostertagia -------------------------------------------------------------------

# Nematode:
# Choose from:
# 1 = H. contortus
# 2 = T. circumcincta
# 3 = O. ostertagi
# 4 = C. oncophora
nematode = 3

ost_out = data.frame(L3_patchA_6m = rep(NA, length(scenarios[,1])), 
                     cum_L3_patchB = rep(NA, length(scenarios[,1])),
                     A_turnout = rep(NA, length(scenarios[,1])),
                     A_end = rep(NA, length(scenarios[,1])))

for (i in 1:length(scenarios[,1])){
  print(i)

    # Start and end dates, lat and lon
    ss = as.Date(scenarios$purch_doy[i], origin = "2022-01-01")
    ee = as.Date(scenarios$purch_doy[i]+sim_duration, origin = "2022-01-01") # run for 6 months after purchase date
    lon = longitudes[i]
    lat = latitudes[i]
    
    # Enter temperature and rainfall data, using the eobs function sourced above:
    temp = eobspoint(data = 'Biosecurity_simulations/tg_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
    precip = eobspoint(data = 'Biosecurity_simulations/rr_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)
    # # some of the lats/lons fall slightly outside of the EOBS grids. 
    # #"next" allows us to skip over those scenarios and keep the look going
    # if(anyNA(temp)) next 
    # if(anyNA(precip)) next
    # 
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
    
    
    stocking_rate = c(rep(2, length(seq.Date(as.Date(ss), as.Date(ee), 'day'))))
    movement_dates = c(ss, as.Date(scenarios$purch_doy[i]+scenarios$quar_t[i], origin = "2022-01-01"), ee)
    grazing_plots = c(scenarios$quar_loc[i], 2, 2)
    
    ginparms = gin(nematode = nematode, temp = temp, precip = precip, photoperiod = daylength(lat = lat, doy = seq(as.Date(ss), as.Date(ee), "days")))
    fecundity_init = exp(ginparms$max.lambda-(ginparms$max.lambda-ginparms$min.lambda)*0.25)
    total_eggs_per_day_per_head = scenarios$fecs[i]*f[1]
    est_adults_from_FEC = total_eggs_per_day_per_head/fecundity_init
    
    # Because we now have the init.vals() function above, we only need to specify state variables that should be non-zero. In this case we have newborn lambs and clean pasture so use the defaults for all state variables.
    initial_values = init.vals(immunity_host = 0.25, q = 0.001, Adult_in_hostA = est_adults_from_FEC)
    
    ostA = gloworm_meta(start = ss, end = ee, lat = lat, 
                        temp = temp, precip = precip, 
                        statevars = initial_values, host = host, nematode = nematode,
                        stocking_rate = stocking_rate, graze = grazing_plots,
                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                        lwt = lw, faeces = f, eventdat = eventdat)
    
    # Quantitites of interest are:
    # residual contamination of the quarantine area 6 months after purchase
    # additional exposure of main herd to L3 contributed by purchased cattle
    # infection levels of purchased cattle at turnout
    # infection levels of purchased cattle 6 months later
    
    ost_out$L3_patchA_6m[i] = ostA[sim_duration,"L3p_A"]
    ost_out$cum_L3_patchB[i] = sum(ostA[,"L3p_B"], na.rm = T)
    ost_out$A_turnout[i] = ostA[scenarios$quar_t[i],"A"]
    ost_out$A_end[i] = ostA[sim_duration,"A"]
}

write.csv(x = cbind(ost_out, scenarios), file = "Biosecurity_simulations/Ostertagia_QoIs.csv")

# GLobal sensitivity analysis to find the most influential factors
# Dubplicate the sa object for each QoI
sa1_ost = sa
sa2_ost = sa
sa3_ost = sa
sa4_ost = sa

# Function to plot SI how I want it
plotSI = function(datS, datT, x_labs, pch_cex, x_axis = TRUE, y_axis = TRUE, legend = TRUE) {
  plot((1:nrow(datT))+0.1, datT$original.1, type = "p", ylim = c(0,1), xlim = c(0.6, 5.4), xaxt = "n", yaxt = "n", xlab = "", pch = 17, cex = pch_cex)
  segments((1:nrow(datT))+0.1, datT$`min. c.i..1`, (1:nrow(datT))+0.1, datT$`max. c.i..1`)
  points((1:nrow(datS))-0.1, datS$original.1, pch = 16, cex = pch_cex)
  segments((1:nrow(datS))-0.1, datS$`min. c.i..1`, (1:nrow(datS))-0.1, datS$`max. c.i..1`)
  if(x_axis == TRUE) axis(1, at = 1:nrow(datS), labels = x_labs)
  if(y_axis == TRUE) axis(2, at = seq(0,1,0.1))
  if(y_axis == TRUE) mtext("Sobol' Index", side = 2, line = 2)
  if(legend == TRUE) legend("topright", legend = c("Main effect", "Total effect"), pch = c(16,17))
}
# Model output is Z-transformed to standardise before computing the Sobol indices

# L3 contamination of quanratine ground after 6 months
tell(sa1_ost, scale(ost_out[,1])) # tell the sa object the model output so that it can compute the sensitivity indices
print(sa1_ost)

plotSI(as.data.frame(sa1_ost$S), as.data.frame(sa1_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Residual contamination", line = 0.5)

saveRDS(sa1_ost, file = "Biosecurity_simulations/sa1_ost.RDS")

# Cumulative L3 on patch B
tell(sa2_ost, scale(ost_out[,2])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa2_ost$S), as.data.frame(sa2_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Additional exposure", line = 0.5)

saveRDS(sa2_ost, file = "Biosecurity_simulations/sa2_ost.RDS")

# Adult burden at time of turnout with main herd
tell(sa3_ost, scale(ost_out[,3])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa3_ost$S), as.data.frame(sa3_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Burden at turnout", line = 0.5)

saveRDS(sa3_ost, file = "Biosecurity_simulations/sa3_ost.RDS")

# Adult burden after 6 months
tell(sa4_ost, scale(ost_out[,4])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa4_ost$S), as.data.frame(sa4_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Residual burden", line = 0.5)

saveRDS(sa4_ost, file = "Biosecurity_simulations/sa4_ost.RDS")

tiff('Biosecurity_simulations/Ostertagia_SA.tiff', width = 7, height = 7, units = 'in', res = 300)

par(mfrow = c(2,2), 
    mar = c(0,0,0,0),
    oma = c(3,4,1,1))

plotSI(as.data.frame(sa1_ost$S), as.data.frame(sa1_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2, x_axis = FALSE, legend = FALSE)
text(0.5, 0.95, "Residual contamination", pos = 4, cex = 1.2)

plotSI(as.data.frame(sa2_ost$S), as.data.frame(sa2_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2, x_axis = FALSE, y_axis = FALSE)
text(0.5, 0.95, "Additional exposure", pos = 4, cex = 1.2)

plotSI(as.data.frame(sa3_ost$S), as.data.frame(sa3_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2, legend = FALSE)
text(0.5, 0.95, "Burden at turnout", pos = 4, cex = 1.2)

plotSI(as.data.frame(sa4_ost$S), as.data.frame(sa4_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2, y_axis = FALSE, legend = FALSE)
text(0.5, 0.95, "Residual burden", pos = 4, cex = 1.2)

dev.off()

# Cooperia -------------------------------------------------------------------

# Nematode:
# Choose from:
# 1 = H. contortus
# 2 = T. circumcincta
# 3 = O. ostertagi
# 4 = C. oncophora
nematode = 4

coop_out = data.frame(L3_patchA_6m = rep(NA, length(scenarios[,1])), 
                     cum_L3_patchB = rep(NA, length(scenarios[,1])),
                     A_turnout = rep(NA, length(scenarios[,1])),
                     A_end = rep(NA, length(scenarios[,1])))

for (i in 1:length(scenarios[,1])){
  print(i)
  
  # Start and end dates, lat and lon
  ss = as.Date(scenarios$purch_doy[i], origin = "2022-01-01")
  ee = as.Date(scenarios$purch_doy[i]+sim_duration, origin = "2022-01-01") # run for 6 months after purchase date
  lon = longitudes[i]
  lat = latitudes[i]
  
  # Enter temperature and rainfall data, using the eobs function sourced above:
  temp = eobspoint(data = 'Biosecurity_simulations/tg_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
  precip = eobspoint(data = 'Biosecurity_simulations/rr_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)
  # # some of the lats/lons fall slightly outside of the EOBS grids. 
  # #"next" allows us to skip over those scenarios and keep the look going
  # if(anyNA(temp)) next 
  # if(anyNA(precip)) next
  # 
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
  
  
  stocking_rate = c(rep(2, length(seq.Date(as.Date(ss), as.Date(ee), 'day'))))
  movement_dates = c(ss, as.Date(scenarios$purch_doy[i]+scenarios$quar_t[i], origin = "2022-01-01"), ee)
  grazing_plots = c(scenarios$quar_loc[i], 2, 2)
  
  ginparms = gin(nematode = nematode, temp = temp, precip = precip, photoperiod = daylength(lat = lat, doy = seq(as.Date(ss), as.Date(ee), "days")))
  fecundity_init = exp(ginparms$max.lambda-(ginparms$max.lambda-ginparms$min.lambda)*0.5)
  total_eggs_per_day_per_head = scenarios$fecs[i]*f[1]
  est_adults_from_FEC = total_eggs_per_day_per_head/fecundity_init
  
  # Because we now have the init.vals() function above, we only need to specify state variables that should be non-zero. In this case we have newborn lambs and clean pasture so use the defaults for all state variables.
  initial_values = init.vals(immunity_host = 0.5, q = 0.001, Adult_in_hostA = est_adults_from_FEC)
  
  coopA = gloworm_meta(start = ss, end = ee, lat = lat, 
                      temp = temp, precip = precip, 
                      statevars = initial_values, host = host, nematode = nematode,
                      stocking_rate = stocking_rate, graze = grazing_plots,
                      movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                      lwt = lw, faeces = f, eventdat = eventdat)
  
  # Quantitites of interest are:
  # residual contamination of the quarantine area 6 months after purchase
  # additional exposure of main herd to L3 contributed by purchased cattle
  # infection levels of purchased cattle at turnout
  # infection levels of purchased cattle 6 months later
  
  coop_out$L3_patchA_6m[i] = coopA[sim_duration,"L3p_A"]
  coop_out$cum_L3_patchB[i] = sum(coopA[,"L3p_B"], na.rm = T)
  coop_out$A_turnout[i] = coopA[scenarios$quar_t[i],"A"]
  coop_out$A_end[i] = coopA[sim_duration,"A"]
}

write.csv(x = cbind(coop_out, scenarios), file = "Biosecurity_simulations/Cooperia_QoIs.csv")

# GLobal sensitivity analysis to find the most influential factors
# Dubplicate the sa object for each QoI
sa1_coop = sa
sa2_coop = sa
sa3_coop = sa
sa4_coop = sa

# Model output is Z-transformed to standardise before computing the Sobol indices

# L3 contamination of quanratine ground after 6 months
tell(sa1_coop, scale(coop_out[,1])) # tell the sa object the model output so that it can compute the sensitivity indices
print(sa1_coop)

plotSI(as.data.frame(sa1_coop$S), as.data.frame(sa1_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Residual contamination", line = 0.5)

saveRDS(sa1_coop, file = "Biosecurity_simulations/sa1_coop.RDS")

# Cumulative L3 on patch B
tell(sa2_coop, scale(coop_out[,2])) # tell the sa object the model output so that it can compute the sensitivity indices
print(sa2_coop)

plotSI(as.data.frame(sa2_coop$S), as.data.frame(sa2_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Additional exposure", line = 0.5)

saveRDS(sa2_coop, file = "Biosecurity_simulations/sa2_coop.RDS")

# Adult burden at time of turnout with main herd
tell(sa3_coop, scale(coop_out[,3])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa3_coop$S), as.data.frame(sa3_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Burden at turnout", line = 0.5)

saveRDS(sa3_coop, file = "Biosecurity_simulations/sa3_coop.RDS")

# Adult burden after 6 months
tell(sa4_coop, scale(coop_out[,4])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa4_coop$S), as.data.frame(sa4_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Residual burden", line = 0.5)

saveRDS(sa4_coop, file = "Biosecurity_simulations/sa4_coop.RDS")

tiff('Biosecurity_simulations/Cooperia_SA.tiff', width = 7, height = 7, units = 'in', res = 300)

par(mfrow = c(2,2), 
    mar = c(0,0,0,0),
    oma = c(3,4,1,1))

plotSI(as.data.frame(sa1_coop$S), as.data.frame(sa1_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2, x_axis = FALSE, legend = FALSE)
text(0.5, 0.95, "Residual contamination", pos = 4, cex = 1.2)

plotSI(as.data.frame(sa2_coop$S), as.data.frame(sa2_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2, x_axis = FALSE, y_axis = FALSE)
text(0.5, 0.95, "Additional exposure", pos = 4, cex = 1.2)

plotSI(as.data.frame(sa3_coop$S), as.data.frame(sa3_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2, legend = FALSE)
text(0.5, 0.95, "Burden at turnout", pos = 4, cex = 1.2)

plotSI(as.data.frame(sa4_coop$S), as.data.frame(sa4_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2, y_axis = FALSE, legend = FALSE)
text(0.5, 0.95, "Residual burden", pos = 4, cex = 1.2)

dev.off()

# Useful links:
# https://stats.stackexchange.com/questions/43504/interpreting-results-from-sobol-sensitivity-analysis-in-r 
# https://uc-ebook.org/docs/html/3_sensitivity_analysis_the_basics.html
# https://towardsdatascience.com/sobol-indices-to-measure-feature-importance-54cedc3281bc

# Visualisations after SA ------------------------------------------------------

dat = read.csv("Biosecurity_simulations/Ostertagia_QoIs.csv")
head(dat)

# Append the coordinates
locs = read.csv("Biosecurity_simulations/random_locations_UK_5000.csv")
dat$latitudes = locs$Y[scenarios$lon_lat]
dat$longitudes = locs$X[scenarios$lon_lat]

tiff('Biosecurity_simulations/Ostertagia_exposure_vs_params.tiff', width = 7, height = 7, units = 'in', res = 300)
par(mfrow = c(2,2), mar = c(4,4,1,1))
plot(dat$cum_L3_patchB~dat$quar_t, xlab = "Quarantine time (days)", ylab = "Additional exposure")
plot(dat$cum_L3_patchB~dat$purch_doy, xlab = "Day of purchase (day 1 = 1st January)", ylab = "Additional exposure")
plot(dat$cum_L3_patchB~dat$fecs, xlab = "FECs on day of purchase", ylab = "Additional exposure")
plot(dat$cum_L3_patchB~dat$latitudes, xlab = "Purchasing farm latitude", ylab = "Additional exposure")
text(x = 59.5, y = 1e08, "FIC", col = "darkgrey")
arrows(59.5, 0.7e08, 59.5, 1e03, col = "darkgrey", length = 0.1)
dev.off()

tiff('Biosecurity_simulations/Ostertagia_residualL3_vs_params.tiff', width = 7, height = 7, units = 'in', res = 300)
par(mfrow = c(2,2), mar = c(4,4,1,1))
plot(dat$L3_patchA_6m~dat$quar_t, xlab = "Quarantine time (days)", ylab = "Residual contamination")
plot(dat$L3_patchA_6m~dat$purch_doy, xlab = "Day of purchase (day 1 = 1st January)", ylab = "Residual contamination")
plot(dat$L3_patchA_6m~dat$fecs, xlab = "FECs on day of purchase", ylab = "Residual contamination")
plot(dat$L3_patchA_6m~dat$latitudes, xlab = "Purchasing farm latitude", ylab = "Residual contamination")
text(x = 59.5, y = 1e06, "FIC", col = "darkgrey")
arrows(59.5, 0.7e06, 59.5, 1e05, col = "darkgrey", length = 0.1)
dev.off()
