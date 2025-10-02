# Biosecurity simulations - Europe
# Uses two sampling matrices generated using the the latin hypercube samples (see "LHS_scenario..." file)

# Load in packages and required source files -----------------------------------
packages = c("dplyr", "rnaturalearth", "sf", "ggplot2", "sensitivity", "deSolve", "geosphere", "chron", "ncdf4", "imputeTS", "forecast", "gplots", "tidyr"

             
) # Package names

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

n = 3000
X1 = data.frame(matrix(runif(5 * n), nrow = n))
colnames(X1) = c("quar_t", "quar_loc", "purch_doy", "fecs", "lon_lat")
X2 = data.frame(matrix(runif(5 * n), nrow = n))
colnames(X2) = c("quar_t", "quar_loc", "purch_doy", "fecs", "lon_lat")

# Generate the "sa" class object using the Sobol function in sensitivity
# Use this after running the simulations
sa = soboljansen(model=NULL, X1 = X1, X2 = X2, nboot=100)
summary(sa)
saveRDS(sa, file = "20241115sa_Europen3000s15000.RDS")
scenarios = as.data.frame(sa$X)

# Quarantine time
quar_t = c(1, 90) # days = 0-12 weeks
# Quarantine location
quar_loc = c(0,1)
# Day of year of purchase
purch_doy = c(1, 365)
# Initial FECs
fecs = c(0, 900)
# Lat/lon
Eur = ne_countries(continent = "europe", returnclass = "sf", scale = 'medium')
plot(Eur[1])
Eur = sf::st_crop(Eur, xmin = -20, xmax = 45,
                  ymin = 30, ymax = 73)
plot(Eur[1])
locations = st_sample(x = Eur, size = 15000)

plot(st_geometry(Eur))
plot(locations, add = T)
points_df = st_coordinates(locations) %>% as.data.frame()

write.csv(points_df, file = "20241115originallocationsn3000s15000.csv")
points_df <- read.csv("20241115originallocationsn3000s15000.csv")

# Initialize tracking columns in points_df for validation stages
points_df$EOBS_available <- FALSE
points_df$DevSuccess_available <- FALSE
points_df$Precipitation_available <- FALSE

# Data frames to capture points_df state at each step for tracking
valid_locations_temp_precip <- data.frame()  # After temperature and precipitation check
valid_locations_dev_success <- data.frame()  # After dev.success check for both species
valid_locations_hmig <- data.frame()         # After h.mig check for both species

# Loop through each location to check for valid data
for (i in 1:nrow(points_df)) {
  ss = "2022-01-01"
  ee = "2022-04-01"
  lon = points_df$X[i]
  lat = points_df$Y[i]
  
  # Print current location index and coordinates for tracking
  cat("Processing location", i, "with coordinates: (", lon, ", ", lat, ")\n")
  
  # Retrieve temperature and precipitation data
  temp <- eobspoint(data = 'tg_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
  precip <- eobspoint(data = 'rr_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)
  
  # Step 1: Check that temperature and precipitation data are not all NAs
  if (all(is.na(temp)) || all(is.na(precip))) {
    cat("Data removed due to all NAs in temperature or precipitation at location", i, "\n")
    next
  }
  points_df$EOBS_available[i] <- TRUE
  
  # Add this to valid locations for temperature and precipitation
  valid_locations_temp_precip <- rbind(valid_locations_temp_precip, points_df[i,])
  
  # Step 2: Calculate dev.1 for Ostertagia and Cooperia and check for sufficient valid values
  dev_1_ost <- pmin(1, pmax(0, -0.07258 + 0.00976 * temp))
  dev_1_coop <- pmin(1, pmax(0, -0.08547 + 0.01021 * temp))
  
  # Ensure at least 2 valid values exist in both dev_1_ost and dev_1_coop
  if (sum(!is.na(dev_1_ost)) < 2 || sum(!is.na(dev_1_coop)) < 2) {
    cat("Data removed due to insufficient valid points in dev.1 for both species at location", i, "\n")
    next
  }
  
  # Step 2b: Calculate dev.success for both species using the moving average only if enough points exist
  devsuccess_ost <- tryCatch({
    as.numeric(na.ma(ma(dev_1_ost, order = 7, centre = FALSE), k = 7))
  }, error = function(e) NA)
  
  devsuccess_coop <- tryCatch({
    as.numeric(na.ma(ma(dev_1_coop, order = 7, centre = FALSE), k = 7))
  }, error = function(e) NA)
  
  # Check that devsuccess is valid (not all NAs or zeros) for both species
  if (all(is.na(devsuccess_ost)) || all(devsuccess_ost == 0) || 
      all(is.na(devsuccess_coop)) || all(devsuccess_coop == 0)) {
    cat("Data removed due to dev.success condition not met for both species at location", i, "\n")
    next
  }
  points_df$DevSuccess_available[i] <- TRUE
  
  # Add this to valid locations after dev.success check
  valid_locations_dev_success <- rbind(valid_locations_dev_success, points_df[i,])
  
  # Step 3: Calculate h.mig for both species and check that h.mig is non-NA and non-zero
  h_mig_ost <- pmax(0, pmin(1, exp(-3.31415 + 0.04637 * precip)))
  h_mig_coop <- pmax(0, pmin(1, exp(-3.29608 + 0.04594 * precip)))
  
  # Ensure h.mig for both species is valid (not all NA or zero)
  if (all(is.na(h_mig_ost)) || all(h_mig_ost == 0) || 
      all(is.na(h_mig_coop)) || all(h_mig_coop == 0)) {
    cat("Data removed due to h.mig being NA or zero for both species at location", i, "\n")
    next
  }
  points_df$Precipitation_available[i] <- TRUE
  
  # Add this to valid locations after h.mig check
  valid_locations_hmig <- rbind(valid_locations_hmig, points_df[i,])
  
  # If all checks are passed, data at this location is valid for both species
  cat("Location", i, "passed all checks.\n")
}

# Save CSVs at each stage
write.csv(valid_locations_temp_precip, file = "14updatedvalid_locations_temp_precipn3000s15000.csv", row.names = FALSE)
write.csv(valid_locations_dev_success, file = "14updatedvalid_locations_dev_succesn3000s15000.csv", row.names = FALSE)
write.csv(valid_locations_hmig, file = "14updatedvalid_locations_hmign3000s15000.csv", row.names = FALSE)


locations <- read.csv("14updatedvalid_locations_hmign3000s15000.csv")


# Optionally, filter out locations that pass all checks and update `points_df`
updatedvalid_locations_final <- points_df[points_df$Precipitation_available == TRUE &
                                     points_df$DevSuccess_available == TRUE, ]

cat("Number of valid locations remaining after all checks: ", nrow(updatedvalid_locations_final), "\n")


write.csv(updatedvalid_locations_final, file = "14updatedvalid_locations_finaln3000s15000.csv")


###creating graphs showing locations####

originallocations <- read.csv("20241115originallocationsn3000s15000.csv")
checkedeobs <- read.csv("14updatedvalid_locations_temp_precipn3000s15000.csv")
stage1 <- read.csv("14updatedvalid_locations_dev_succesn3000s15000.csv")
stage2 <- read.csv("14updatedvalid_locations_hmign3000s15000.csv")
#stage3 <- read.csv("20241112_valid_locations_finaln4000s20000.csv")
cat("Number of valid locations remaining after all checks: ", nrow(checkedeobs), "\n")
cat("Number of valid locations remaining after all checks: ", nrow(stage1), "\n")

# Get the European map
Eur = ne_countries(continent = "europe", returnclass = "sf", scale = 'medium')

Eur = sf::st_crop(Eur, xmin = -20, xmax = 45, ymin = 30, ymax = 73)

originallocations <- st_as_sf(originallocations, coords = c("X", "Y"), crs = 4326)

##Initial Selection Sampling of Locations
dev.off()
tiff('20241114_originallocationsn2000s10000.tiff', width = 7, height = 7, units = 'in', res = 300)

ggplot() +
  geom_sf(data = Eur, fill = "lightgray", color = "white") +  # Plot Europe
  geom_sf(data = originallocations, color = "black", size = 2, shape = 1) +  
  labs(title = "Original Locations in Europe") +
  theme_minimal() +
  theme(legend.position = "none")

dev.off()

##After EOBS has checked for data presence

checkedeobs <- st_as_sf(checkedeobs, coords = c("X", "Y"), crs = 4326)
dev.off()
tiff('20241115_eobscheckn3000s15000.tiff', width = 7, height = 7, units = 'in', res = 300)

ggplot() +
  geom_sf(data = Eur, fill = "lightgray", color = "white") +  # Plot Europe
  geom_sf(data = checkedeobs, color = "black", size = 2, shape = 1) +  
  labs(title = "Checked EOBs") +
  theme_minimal() +
  theme(legend.position = "none")

##Stage 1 -- checking weather data is not all NA's

stage1 <- st_as_sf(stage1, coords = c("X", "Y"), crs = 4326)

dev.off()
tiff('20241115_tempcheckn3000s15000.tiff', width = 7, height = 7, units = 'in', res = 300)


ggplot() +
  geom_sf(data = Eur, fill = "lightgray", color = "white") +  # Plot Europe
  geom_sf(data = stage1, color = "black", size = 2, shape = 1) +  
  labs(title = "Weather data availability in Europe") +
  theme_minimal() +
  theme(legend.position = "none")

sum(is.na(stage1$X))  # Should be zero after Stage 1
sum(is.na(stage2$Y))  # Should be zero after Stage 1

##Stage 2 -- checking temperature reaches requirements for both Cooperia and Ostertagia 

stage2 <- st_as_sf(stage2, coords = c("X", "Y"), crs = 4326)

dev.off()
tiff('20241115_precipcheckn3000s15000.tiff', width = 7, height = 7, units = 'in', res = 300)


ggplot() +
  geom_sf(data = Eur, fill = "lightgray", color = "white") +  # Plot Europe
  geom_sf(data = stage2, color = "black", size = 2, shape = 1) +  
  labs(title = "Weather data availability in Europe") +
  theme_minimal() +
  theme(legend.position = "none")

##Stage 3 -- checking valid locations requirements for both Cooperia and Ostertagia 

#stage3 <- st_as_sf(stage3, coords = c("X", "Y"), crs = 4326)

#dev.off()
#tiff('20241113_finalchecksn4000s20000.tiff', width = 7, height = 7, units = 'in', res = 300)


#ggplot() +
#  geom_sf(data = Eur, fill = "lightgray", color = "white") +  # Plot Europe
#  geom_sf(data = stage3, color = "black", size = 2, shape = 1) +  
#  labs(title = "Weather data availability in Europe") +
#  theme_minimal() +
#  theme(legend.position = "none")


#######

# Rescale scenario parameters using only valid points
rescale_scenarios <- function(x, newMax, newMin) {
  (x - min(x)) / (max(x) - min(x)) * (newMax - newMin) + newMin
}

# Rescale scenario parameters
scenarios[,1] <- round(rescale_scenarios(x = scenarios[,1], newMin = quar_t[1], newMax = quar_t[2]))
scenarios[,2] <- round(scenarios[,2])
scenarios[,3] <- round(rescale_scenarios(x = scenarios[,3], newMin = purch_doy[1], newMax = purch_doy[2]))
scenarios[,4] <- round(rescale_scenarios(x = scenarios[,4], newMin = fecs[1], newMax = fecs[2]))

# Set lon_lat to refer to valid indices in valid_locations_final
num_points <- nrow(updatedvalid_locations_final)
scenarios[,5] <- pmin(num_points, pmax(1, round(rescale_scenarios(x = scenarios[,5], newMin = 1, newMax = num_points))))

# Rename columns for clarity
colnames(scenarios) <- c("quar_t", "quar_loc", "purch_doy", "fecs", "lon_lat")

# Match locations based on updated lon_lat in valid_locations_final
locs <- updatedvalid_locations_final # Use filtered points directly
latitudes <- locs$Y[scenarios$lon_lat]
longitudes <- locs$X[scenarios$lon_lat]


# Filter valid locations
#valid_locations_final <- points_df[points_df$Precipitation_available == TRUE &
#                                     points_df$DevSuccess_available == TRUE, ]
#cat("Number of valid locations remaining after all checks: ", nrow(valid_locations_final), "\n")
#write.csv(valid_locations_final, file = "20241108_valid_locations_final.csv")

# Rescale scenario parameters using only valid points
#rescale_scenarios <- function(x, newMax, newMin) {
#  (x - min(x)) / (max(x) - min(x)) * (newMax - newMin) + newMin
#}

# Rescale scenario parameters
#scenarios[,1] <- round(rescale_scenarios(x = scenarios[,1], newMin = quar_t[1], newMax = quar_t[2]))
#scenarios[,2] <- round(scenarios[,2])  # Quarantine location (already integer)
#scenarios[,3] <- round(rescale_scenarios(x = scenarios[,3], newMin = purch_doy[1], newMax = purch_doy[2]))
#scenarios[,4] <- round(rescale_scenarios(x = scenarios[,4], newMin = fecs[1], newMax = fecs[2]))

# Rescale lon_lat to match valid locations and get indices
#num_valid_points <- nrow(valid_locations_final)
#scenarios[, 5] <- pmin(num_valid_points, pmax(1, round(rescale_scenarios(scenarios[, 5], 1, num_valid_points))))

# Add valid latitudes and longitudes to scenarios
#scenarios$latitude <- valid_locations_final$Y[scenarios$lon_lat]
#scenarios$longitude <- valid_locations_final$X[scenarios$lon_lat]

# OPTIONAL: Save to CSV (for debugging)
#write.csv(scenarios, "filtered_scenarios.csv", row.names = FALSE)

# OPTIONAL: View the first few rows of the filtered scenarios
#head(scenarios)

head(scenarios)

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

for (i in 1:length(scenarios[, 1])) {
  print(i)
  
  # Wrap the entire scenario processing logic in tryCatch
  tryCatch({
    # Start and end dates, lat and lon
    ss = as.Date(scenarios$purch_doy[i], origin = "2022-01-01")
    ee = as.Date(scenarios$purch_doy[i] + sim_duration, origin = "2022-01-01") # run for 6 months after purchase date
    lon = longitudes[i]
    lat = latitudes[i]
    
    
    temp = eobspoint(data = 'tg_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
    precip = eobspoint(data = 'rr_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)
    
    # Estimate weights and use this to estimate dry matter intake and faeces production
    host.age = seq(host_age, host_age + as.numeric(as.Date(ee) - as.Date(ss)), 1)
    lw = if (host == 1) { lw_sheep(age = host.age) } else { lw_cattle(age = host.age) }
    
    DMI = dmi(lwt = lw, age = host.age, host_species = host)
    
    # Create a vector of kgDM values for each day of the simulation
    kgDMha = rep(kgDM, length(DMI))
    f = faeces(lwt = lw, host_species = host)
    
    stocking_rate = c(rep(2, length(seq.Date(as.Date(ss), as.Date(ee), 'day'))))
    movement_dates = c(ss, as.Date(scenarios$purch_doy[i] + scenarios$quar_t[i], origin = "2022-01-01"), ee)
    grazing_plots = c(scenarios$quar_loc[i], 2, 2)
    
    
    ginparms = gin(nematode = nematode, temp = temp, precip = precip, photoperiod = daylength(lat = lat, doy = seq(as.Date(ss), as.Date(ee), "days")))
    fecundity_init = exp(ginparms$max.lambda-(ginparms$max.lambda-ginparms$min.lambda)*0.25)
    total_eggs_per_day_per_head = scenarios$fecs[i]*f[1]
    est_adults_from_FEC = total_eggs_per_day_per_head/fecundity_init
    
    # Initial values for the model
    initial_values = init.vals(immunity_host = 0.25, q = 0.001, Adult_in_hostA = est_adults_from_FEC)
    
    # Use gloworm_meta to compute ostA
    ostA = gloworm_meta(start = ss, end = ee, lat = lat,
                        temp = temp, precip = precip,
                        statevars = initial_values, host = host, nematode = nematode,
                        stocking_rate = stocking_rate, graze = grazing_plots,
                        movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                        lwt = lw, faeces = f, eventdat = eventdat)
    
    # Proceed with quantities of interest
    ost_out$L3_patchA_6m[i] = ostA[sim_duration, "L3p_A"]
    ost_out$cum_L3_patchB[i] = sum(ostA[, "L3p_B"], na.rm = TRUE)
    ost_out$A_turnout[i] = ostA[scenarios$quar_t[i], "A"]
    ost_out$A_end[i] = ostA[sim_duration, "A"]
    
  }, error = function(e) {
    message(paste("Error in scenario", i, ":", e$message))  # Log the error message
    return(NULL)  # Return NULL in case of error
  })
}

write.csv(x = cbind(ost_out, scenarios), file = "20241115Ostertagia_QoIs_Europen3000s15000.csv")


# Global sensitivity analysis to find the most influential factors
# Duplicate the sa object for each QoI
sa1_ost = sa
sa2_ost = sa
sa3_ost = sa
sa4_ost = sa

plotSI = function(datS, datT, x_labs, pch_cex, x_axis = TRUE, y_axis = TRUE, legend = TRUE) {
  plot((1:nrow(datT)) + 0.1, datT$original, type = "p", ylim = c(0, 1), xlim = c(0.6, 5.4), xaxt = "n", yaxt = "n", xlab = "", pch = 17, cex = pch_cex)
  segments((1:nrow(datT)) + 0.1, datT$`min. c.i.`, (1:nrow(datT)) + 0.1, datT$`max. c.i.`)
  points((1:nrow(datS)) - 0.1, datS$original, pch = 16, cex = pch_cex)
  segments((1:nrow(datS)) - 0.1, datS$`min. c.i.`, (1:nrow(datS)) - 0.1, datS$`max. c.i.`)
  if (x_axis == TRUE) axis(1, at = 1:nrow(datS), labels = x_labs)
  if (y_axis == TRUE) axis(2, at = seq(0, 1, 0.1))
  if (y_axis == TRUE) mtext("Sobol' Index", side = 2, line = 2)
  if (legend == TRUE) legend("topright", legend = c("Main effect", "Total effect"), pch = c(16, 17))
}

any(is.na(ost_out[,1])) || any(is.infinite(ost_out[,1]))
which(is.na(ost_out[,1]) | is.infinite(ost_out[,1]))
for (col in 1:ncol(ost_out)) {
  Interpolate NA and Inf values using linear interpolation
  ost_out[, col] <- approx(1:nrow(ost_out), ost_out[, col], xout = 1:nrow(ost_out), method = "linear", rule = 2)$y
  }

for (col in 1:ncol(ost_out)) {
  # Extract the column
  column_data <- ost_out[, col]
  
  # Identify non-missing and finite values
  valid_idx <- which(!is.na(column_data) & is.finite(column_data))
  
  # Check if there are at least two valid points
  if (length(valid_idx) > 1) {
    # Perform linear interpolation
    column_data <- approx(x = valid_idx,
                          y = column_data[valid_idx],
                          xout = 1:nrow(ost_out),
                          method = "linear",
                          rule = 2)$y
  }
  
  # Replace the column with interpolated data, keeping original where interpolation isn't possible
  ost_out[, col] <- column_data
}


tell(sa1_ost, scale(ost_out[,1])) 
print(sa1_ost)

plotSI(as.data.frame(sa1_ost$S), as.data.frame(sa1_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Residual contamination", line = 0.5)

saveRDS(sa1_ost, file = "20241114_sa1_ost_Europen3000s15000.RDS")

# Cumulative L3 on patch B
tell(sa2_ost, scale(ost_out[,2])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa2_ost$S), as.data.frame(sa2_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Additional exposure", line = 0.5)

saveRDS(sa2_ost, file = "20241115_sa2_ost_Europen3000s15000.RDS")

# Adult burden at time of turnout with main herd
tell(sa3_ost, scale(ost_out[,3])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa3_ost$S), as.data.frame(sa3_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Burden at turnout", line = 0.5)

saveRDS(sa3_ost, file = "20241154_sa3_ost_Europen3000s15000.RDS")

# Adult burden after 6 months
tell(sa4_ost, scale(ost_out[,4])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa4_ost$S), as.data.frame(sa4_ost$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Residual burden", line = 0.5)

saveRDS(sa4_ost, file = "20241115_sa4_ost_Europen3000s15000.RDS")


dev.off()
tiff('20241115_Ostertagia_SA_Europen3000s15000.tiff', width = 7, height = 7, units = 'in', res = 300)

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
  temp = eobspoint(data = 'tg_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
  precip = eobspoint(data = 'rr_ens_mean_0.25deg_reg_2011-2023_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)
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

write.csv(x = cbind(coop_out, scenarios), file = "20241115_Cooperia_QoIs_Europen3000s15000.csv")

coop_out <- read.csv("20241115_Cooperia_QoIs_Europen3000s15000.csv")

#for (col in 1:ncol(coop_out)) {
#  Interpolate NA and Inf values using linear interpolation
#coop_out[, col] <- approx(1:nrow(coop_out), coop_out[, col], xout = 1:nrow(coop_out), method = "linear", rule = 2)$y
#}

any(is.na(coop_out[,1])) || any(is.infinite(coop_out[,1]))
which(is.na(coop_out[,1]) | is.infinite(coop_out[,1]))

for (col in 1:ncol(coop_out)) {
  # Extract the column
  column_data <- coop_out[, col]
  
  # Identify non-missing and finite values
  valid_idx <- which(!is.na(column_data) & is.finite(column_data))
  
  # Check if there are at least two valid points
  if (length(valid_idx) > 1) {
    # Perform linear interpolation
    column_data <- approx(x = valid_idx,
                          y = column_data[valid_idx],
                          xout = 1:nrow(coop_out),
                          method = "linear",
                          rule = 2)$y
  }
  
  # Replace the column with interpolated data, keeping original where interpolation isn't possible
 coop_out[, col] <- column_data
}
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

saveRDS(sa1_coop, file = "20241115sa1_coop_Europen3000s15000.RDS")

# Cumulative L3 on patch B
tell(sa2_coop, scale(coop_out[,2])) # tell the sa object the model output so that it can compute the sensitivity indices
print(sa2_coop)

plotSI(as.data.frame(sa2_coop$S), as.data.frame(sa2_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Additional exposure", line = 0.5)

saveRDS(sa2_coop, file = "20241115sa2_coop_Europen3000s15000.RDS")

# Adult burden at time of turnout with main herd
tell(sa3_coop, scale(coop_out[,3])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa3_coop$S), as.data.frame(sa3_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Burden at turnout", line = 0.5)

saveRDS(sa3_coop, file = "20241115sa3_coop_Europen3000s15000.RDS")

# Adult burden after 6 months
tell(sa4_coop, scale(coop_out[,4])) # tell the sa object the model output so that it can compute the sensitivity indices

plotSI(as.data.frame(sa4_coop$S), as.data.frame(sa4_coop$T), c("Q(t)", "Q(loc)", "Date", "FEC", "Location"), pch_cex = 1.2)
title("Residual burden", line = 0.5)

saveRDS(sa4_coop, file = "20241115sa4_coop_Europen3000s15000.RDS")
dev.off()


tiff('20241115Cooperia_SA_Europen3000s15000.tiff', width = 7, height = 7, units = 'in', res = 300)

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

library(ggplot2)
library(gridExtra)

plot1 <-ggplot(coop_out, aes(x = quar_t, y = cum_L3_patchB)) +
  geom_point(alpha = 0.7) +
  labs(title = "Additional Exposure vs Duration of Quarantine",
       x = "Duration of Quarantine (days)",
       y = "Additional Exposure (cum_L3_patchB)") +
  scale_x_continuous(breaks = seq(0, 90, by = 10)) +  # Set x-axis breaks for quar_t (0-90)
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),   # Adjust title size
    axis.title = element_text(size = 9),                  # Adjust axis titles
    axis.text = element_text(size = 8),                   # Adjust axis labels
    plot.margin = margin(5, 5, 5, 5)                      # Adjust plot margins (top, right, bottom, left)
  )

plot2 <- ggplot(coop_out, aes(x = purch_doy, y = cum_L3_patchB)) +
  geom_point(alpha = 0.7) +
  labs(title = "Additional Exposure vs Purchase Date",
       x = "Purchase Day (day 1 = 1st January)",
       y = "Additional Exposure (cum_L3_patchB)") + # Set x-axis breaks for quar_t (0-90)
  theme_minimal()  + 
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),   # Adjust title size
    axis.title = element_text(size = 9),                  # Adjust axis titles
    axis.text = element_text(size = 8),                   # Adjust axis labels
    plot.margin = margin(5, 5, 5, 5)                      # Adjust plot margins (top, right, bottom, left)
  )

plot3 <- ggplot(coop_out, aes(x = fecs, y = cum_L3_patchB)) +
  geom_point(alpha = 0.7) +
  labs(title = "Additional Exposure vs Initial FECs",
       x = "FECs on day of purchase",
       y = "Additional Exposure (cum_L3_patchB)") + # Set x-axis breaks for quar_t (0-90)
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),   # Adjust title size
    axis.title = element_text(size = 9),                  # Adjust axis titles
    axis.text = element_text(size = 8),                   # Adjust axis labels
    plot.margin = margin(5, 5, 5, 5)                      # Adjust plot margins (top, right, bottom, left)
  )


# Merge the two data frames based on the common identifier 'lon_lat' and 'X.1'
merged_data <- merge(coop_out, locations, by.x = "lon_lat", by.y = "X.1", all.x = TRUE)

# Add the latitude column 'Y' from 'locations' as a new column 'lats' in merged_data
merged_data$lats <- merged_data$Y

# Check the merged data to confirm that 'lats' has been added correctly
head(merged_data)

# Plot residual contamination against latitude
plot4 <- ggplot(merged_data, aes(x = lats, y = cum_L3_patchB)) + 
  geom_point(alpha = 0.7) +  # Scatter plot with transparency
  labs(title = "Additional Exposure vs Farm Location",
       x = "Farm Location Latitude",
       y = "Additional Exposure (cum_L3_patchB)") +
  theme_minimal()  + 
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),   # Adjust title size
    axis.title = element_text(size = 9),                  # Adjust axis titles
    axis.text = element_text(size = 8),                   # Adjust axis labels
    plot.margin = margin(5, 5, 5, 5)                      # Adjust plot margins (top, right, bottom, left)
  )
dev.off()

tiff('Cooperia_AdditionalExposure.tiff', width = 7, height = 7, units = 'in', res = 300)
grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
dev.off()
