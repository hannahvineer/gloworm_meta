# Author: Olivia Ingle - olivia.ingle@liverpool.ac.uk
# Created: 04/07/2024
# Last modified: 02/08/2024
# Last modified by: Hannah Vineer

# Load in packages and required source files -----------------------------------
packages = c("deSolve", "geosphere", "chron", "ncdf4", "imputeTS", "forecast", "gplots", "tidyr") # Package names

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

# McFarland Validation ---------------------------------------------------------

# 4 herds: A-D
# Parasitism prior to the start of the trial unknown.
# 10-12 months of age at turnout (start of trial).
# Clean pasture at turnout.
# Turnout 19th April 2019.
# Doramectin treatment 28th May 2019, with approx. 78% FECR across the 4 groups.
# Housed 15th October 2019.
# Each group rotiationally grazed 6 unshared fields.
# 19th September 2019, groups A, B and D reduced to 5 individuals, group C reduced to 6 individuals (from 11 per herd).
# Stocking rates variable around 40 cattle/ha, reducing to around 20/ha in September (supplementary material, McFarland et al. 2022)
# 3000 kg DM ha-1
# FECs (from supp material) were approx 83% Ostertagia.

# All herds ====================================================================

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

# Read in validation data ####
McF_dat = read.csv("Validation_data/McFarland/published_data/Raw FEC data Mendeley Data.csv")

# Replace column headers with the DAY OF TRIAL to assist with plotting later on
FEC_dates = c(3:length(colnames(McF_dat)))
for (i in 3:length(colnames(McF_dat))) {
  FEC_dates[i-2] = strsplit(colnames(McF_dat)[i], split = "G.")[[1]][2]
}
FEC_dates = as.Date(FEC_dates, "%d.%m.%Y")
FEC_days = which(seq.Date(from = as.Date(ss), to = as.Date(ee), by = "days") %in% FEC_dates)
colnames(McF_dat)[3:12] = as.character(FEC_days)
long = McF_dat %>% pivot_longer(cols = '40':'180', names_to = 'day', values_to = 'FEC') %>% as.data.frame()
long[,'day'] = as.numeric(long[,'day'])
herdA = subset(x = long, subset = long$Group_ID=="A")
herdB = subset(x = long, subset = long$Group_ID=="B")
herdC = subset(x = long, subset = long$Group_ID=="C")
herdD = subset(x = long, subset = long$Group_ID=="D")

# Herd-specific data ===========================================================

# Host management dates and rates
stocking_rate_A = c(rep(40, length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))), rep(20, length(temp)-length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))))
movement_dates_A = c(ss, '2019-04-28', '2019-05-06', '2019-05-21', '2019-05-30', '2019-06-03', '2019-06-06', '2019-06-09', '2019-06-13', '2019-06-17', '2019-06-21', '2019-06-24', '2019-06-27', '2019-07-02', '2019-07-08', '2019-07-14', '2019-07-23', '2019-07-28', '2019-08-01', '2019-08-07', '2019-08-10', '2019-08-12', '2019-08-14', '2019-08-17', '2019-08-19', '2019-08-22', '2019-08-25', '2019-08-30', '2019-09-02', '2019-09-06', '2019-09-08', '2019-09-12', '2019-09-15', '2019-09-17', '2019-09-19', '2019-09-25', '2019-09-29', '2019-10-03', '2019-10-06', '2019-10-10', '2019-10-11', ee)
grazing_plots_A = c(1, 2, 5, 1, 2, 3, 4, 5, 6, 1, 3, 2, 4, 5, 6, 1, 3, 2, 4, 5, 6, 1, 4, 3, 2, 4, 5, 6, 1, 3, 2, 4, 5, 6, 1, 3, 2, 6, 5, 4, 0, 0)

stocking_rate_B = c(rep(40, length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))), rep(20, length(temp)-length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))))
movement_dates_B = c(ss, '2019-04-23', '2019-05-07', '2019-05-26', '2019-05-31', '2019-06-05', '2019-06-09', '2019-06-14', '2019-06-19', '2019-06-24', '2019-06-30', '2019-07-04', '2019-07-09', '2019-07-15', '2019-07-20', '2019-07-28', '2019-08-05', '2019-08-07', '2019-08-11', '2019-08-14', '2019-08-16', '2019-08-20', '2019-08-23', '2019-08-26', '2019-08-31', '2019-09-02', '2019-09-06', '2019-09-10', '2019-09-13', '2019-09-15', '2019-09-18', '2019-09-23', '2019-09-28', '2019-10-03', '2019-10-07', '2019-10-10', '2019-10-14', '2019-10-15',ee)
grazing_plots_B = c(1, 2, 3, 1, 2, 4, 5, 6, 2, 3, 1, 4, 5, 6, 2, 3, 1, 4, 5, 6, 2, 3, 1, 4, 6, 2, 5, 3, 1, 4, 6, 2, 3, 5, 1, 4, 6, 0, 0)

stocking_rate_C = c(rep(40, length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))), rep(20, length(temp)-length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))))
movement_dates_C = c(ss, '2019-04-28', '2019-05-06', '2019-05-15', '2019-05-26', '2019-05-31', '2019-06-03', '2019-06-06', '2019-06-09', '2019-06-12', '2019-06-17', '2019-06-19', '2019-06-23', '2019-06-25', '2019-06-28', '2019-07-03', '2019-07-09', '2019-07-18', '2019-07-22', '2019-07-26', '2019-08-01', '2019-08-05', '2019-08-08', '2019-08-12', '2019-08-15', '2019-08-16', '2019-08-19', '2019-08-24', '2019-08-27', '2019-09-02', '2019-09-06', '2019-09-09', '2019-09-12', '2019-09-15', '2019-09-17',  '2019-09-22', '2019-09-28', '2019-10-03', '2019-10-06', '2019-10-07', '2019-10-09', '2019-10-11', '2019-10-13', '2019-10-14', ee)
grazing_plots_C = c(1, 2, 3, 1, 2, 4, 5, 6, 3, 1, 2, 5, 4, 3, 6, 2, 1, 4, 3, 5, 2, 6, 1, 4, 3, 5, 2, 6, 1, 4, 3, 5, 2, 6, 4, 1, 5, 3, 2, 6, 4, 1, 5, 0, 0)

stocking_rate_D = c(rep(40, length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))), rep(20, length(temp)-length(seq.Date(as.Date(ss), as.Date('2019-09-19'), 'day'))))
movement_dates_D = c(ss, '2019-04-28', '2019-05-06', '2019-05-20', '2019-06-01', '2019-06-03', '2019-06-05', '2019-06-10', '2019-06-13', '2019-06-17', '2019-06-21', '2019-06-24',      '2019-06-28', '2019-07-02', '2019-07-08', '2019-07-15', '2019-07-23', '2019-07-28', '2019-08-02', '2019-08-07', '2019-08-10', '2019-08-12', '2019-08-15', '2019-08-17', '2019-08-19', '2019-08-23', '2019-08-29', '2019-09-03', '2019-09-06', '2019-09-09', '2019-09-12', '2019-09-15', '2019-09-17', '2019-09-19', '2019-09-25', '2019-09-29', '2019-10-02', '2019-10-06', '2019-10-10', '2019-10-11', ee)
grazing_plots_D = c(1, 2, 5, 1, 2, 3, 4, 5, 6, 1, 3, 2, 4, 5, 6, 1, 3, 2, 4, 5, 6, 1, 4, 3, 2, 4, 6, 1, 3, 2, 4, 5, 6, 1, 3, 2, 6, 5, 4, 0, 0)


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

# Herd A ####
ostA = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate_A, graze = grazing_plots_A,
                    movements = movement_dates_A, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)

# plot.ts(ostA[,1:10])
# plot.ts(ostA[,11:20])
# plot.ts(ostA[,21:30])
# plot.ts(ostA[,31:40])
# plot.ts(ostA[,41:50])
# plot.ts(ostA[,51:56])

# par(mfrow = c(1,1))
# boxplot(herdA$FEC*0.83 ~ herdA$day, at = herdA$day[1:10], boxwex = 5, xlim = c(1,365), xaxt = 'n', ylab = 'FEC (epg)') 
# axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
# mtext("Days after turnout", 1, 3)
# lines(ostA[,'FEC'])
# abline(v = 40, lty = 1, col = 'lightgrey')
# text(x = 35, y = 200, col = 'grey', labels = "Doramectin Tx", srt = 90)
# abline(v = 180, lty = 1, col = 'lightgrey')
# text(x = 175, y = 200, col = 'grey', labels = "Housed", srt = 90)

herdA$simulated = ostA[herdA$day, 'FEC']

validation_output = data.frame(day = herdA$day, 
                               simulated = herdA$simulated, 
                               observed = herdA$FEC*0.83,
                               species = rep('Ostertagia', length(herdA$day)),
                               study = rep('HerdA', length(herdA$day)))

# Herd B ####

ostB = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate_B, graze = grazing_plots_B,
                    movements = movement_dates_B, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)

# plot.ts(ostB[,1:10])
# plot.ts(ostB[,11:20])
# plot.ts(ostB[,21:30])
# plot.ts(ostB[,31:40])
# plot.ts(ostB[,41:50])
# plot.ts(ostB[,51:56])
# 
# par(mfrow = c(1,1))
# boxplot(herdB$FEC*0.83 ~ herdB$day, at = herdB$day[1:10], boxwex = 5, xlim = c(1,365), xaxt = 'n', ylab = 'FEC (epg)') 
# axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
# mtext("Days after turnout", 1, 3)
# lines(ostB[,'FEC'])
# abline(v = 40, lty = 1, col = 'lightgrey')
# text(x = 35, y = 200, col = 'grey', labels = "Doramectin Tx", srt = 90)
# abline(v = 180, lty = 1, col = 'lightgrey')
# text(x = 175, y = 200, col = 'grey', labels = "Housed", srt = 90)

herdB$simulated = ostB[herdB$day, 'FEC']

# Save output for statistical analysis
validation_output = rbind(validation_output, data.frame(day = herdB$day, 
                               simulated = herdB$simulated, 
                               observed = herdB$FEC*0.83,
                               species = rep('Ostertagia', length(herdB$day)),
                               study = rep('HerdB', length(herdB$day))))

# Herd C ####

ostC = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate_C, graze = grazing_plots_C,
                    movements = movement_dates_C, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)

# plot.ts(ostC[,1:10])
# plot.ts(ostC[,11:20])
# plot.ts(ostC[,21:30])
# plot.ts(ostC[,31:40])
# plot.ts(ostC[,41:50])
# plot.ts(ostC[,51:56])
# 
# par(mfrow = c(1,1))
# boxplot(herdC$FEC*0.83 ~ herdC$day, at = herdC$day[1:10], boxwex = 5, xlim = c(1,365), xaxt = 'n', ylab = 'FEC (epg)') 
# axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
# mtext("Days after turnout", 1, 3)
# lines(ostC[,'FEC'])
# abline(v = 40, lty = 1, col = 'lightgrey')
# text(x = 35, y = 100, col = 'grey', labels = "Doramectin Tx", srt = 90)
# abline(v = 180, lty = 1, col = 'lightgrey')
# text(x = 175, y = 100, col = 'grey', labels = "Housed", srt = 90)

herdC$simulated = ostC[herdC$day, 'FEC']

# Save output for statistical analysis
validation_output = rbind(validation_output, data.frame(day = herdC$day, 
                                                        simulated = herdC$simulated, 
                                                        observed = herdC$FEC*0.83,
                                                        species = rep('Ostertagia', length(herdC$day)),
                                                        study = rep('HerdC', length(herdC$day))))

# Herd D ####

ostD = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate_D, graze = grazing_plots_D,
                    movements = movement_dates_D, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)

# plot.ts(ostD[,1:10])
# plot.ts(ostD[,11:20])
# plot.ts(ostD[,21:30])
# plot.ts(ostD[,31:40])
# plot.ts(ostD[,41:50])
# plot.ts(ostD[,51:56])
# 
# par(mfrow = c(1,1))
# boxplot(herdD$FEC*0.83 ~ herdD$day, at = herdD$day[1:10], boxwex = 5, xlim = c(1,365), xaxt = 'n', ylab = 'FEC (epg)') 
# axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
# mtext("Days after turnout", 1, 3)
# lines(ostD[,'FEC'])
# abline(v = 40, lty = 1, col = 'lightgrey')
# text(x = 35, y = 200, col = 'grey', labels = "Doramectin Tx", srt = 90)
# abline(v = 180, lty = 1, col = 'lightgrey')
# text(x = 175, y = 200, col = 'grey', labels = "Housed", srt = 90)

herdD$simulated = ostD[herdD$day, 'FEC']

# Save output for statistical analysis
validation_output = rbind(validation_output, data.frame(day = herdD$day, 
                                                        simulated = herdD$simulated, 
                                                        observed = herdD$FEC*0.83,
                                                        species = rep('Ostertagia', length(herdD$day)),
                                                        study = rep('HerdD', length(herdD$day))))
write.csv(validation_output, file = "Ost_McFarland_validation_output.csv")

# Stats ####
lm_validation_McF = lm(validation_output$observed.mean ~ 0 + validation_output$simulated)
summary(lm_validation_McF)
plot(lm_validation_McF) # some deviance from normality in the residuals. This would be because the FEC data are likely to be negative binomially distributed. 

plot(validation_output$observed.mean ~ validation_output$simulated, ylab = 'observed FEC', xlab = 'simulated FEC', pch = 21)
abline(a = 0, b = 1)
abline(a = 0, b = lm_validation_McF$coefficients[1], lty = 2)

capture.output(summary(lm_validation_McF), file = 'Ost_McFarland_stats.txt')
cor_validation_McF = cor.test(validation_output$simulated, validation_output$observed.mean, method = "spearman")
capture.output(cor_validation_McF, file = 'Ost_McFarland_stats.txt', append = T)

# Final plots ####
dev.off()
tiff('Ost_McFarland_longitudinal.tiff', width = 8, height = 8, units = 'in', res = 300)
par(mfrow = c(2,2))
boxplot(herdA$FEC*0.83 ~ herdA$day, at = herdA$day[1:10], boxwex = 10, xlim = c(1,365), ylim = c(0, 500), xaxt = 'n', ylab = 'FEC (epg)', xlab = "Days after turnout") 
axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
lines(ostA[,'FEC'])
abline(v = 40, lty = 1, col = 'lightgrey')
text(x = 30, y = 400, col = 'grey', labels = "Doramectin Tx", srt = 90)
abline(v = 180, lty = 1, col = 'lightgrey')
text(x = 170, y = 400, col = 'grey', labels = "Housed", srt = 90)
text(x = 300, y = 450, labels = "Herd A")

boxplot(herdB$FEC*0.83 ~ herdB$day, at = herdB$day[1:10], boxwex = 10, xlim = c(1,365), ylim = c(0, 500), xaxt = 'n', ylab = 'FEC (epg)', xlab = "Days after turnout") 
axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
lines(ostB[,'FEC'])
abline(v = 40, lty = 1, col = 'lightgrey')
text(x = 30, y = 400, col = 'grey', labels = "Doramectin Tx", srt = 90)
abline(v = 180, lty = 1, col = 'lightgrey')
text(x = 170, y = 400, col = 'grey', labels = "Housed", srt = 90)
text(x = 300, y = 450, labels = "Herd B")

boxplot(herdC$FEC*0.83 ~ herdC$day, at = herdC$day[1:10], boxwex = 10, xlim = c(1,365), ylim = c(0, 500), xaxt = 'n', ylab = 'FEC (epg)', xlab = "Days after turnout") 
axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
lines(ostC[,'FEC'])
abline(v = 40, lty = 1, col = 'lightgrey')
text(x = 30, y = 400, col = 'grey', labels = "Doramectin Tx", srt = 90)
abline(v = 180, lty = 1, col = 'lightgrey')
text(x = 170, y = 400, col = 'grey', labels = "Housed", srt = 90)
text(x = 300, y = 450, labels = "Herd C")

boxplot(herdD$FEC*0.83 ~ herdD$day, at = herdD$day[1:10], boxwex = 10, xlim = c(1,365), ylim = c(0, 500), xaxt = 'n', ylab = 'FEC (epg)', xlab = "Days after turnout") 
axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
lines(ostD[,'FEC'])
abline(v = 40, lty = 1, col = 'lightgrey')
text(x = 30, y = 400, col = 'grey', labels = "Doramectin Tx", srt = 90)
abline(v = 180, lty = 1, col = 'lightgrey')
text(x = 170, y = 400, col = 'grey', labels = "Housed", srt = 90)
text(x = 300, y = 450, labels = "Herd D")

dev.off()

tiff('Ost_McFarland_sim_obs.tiff', width = 5, height = 5, units = 'in', res = 300)

plot(validation_output$observed.mean ~ validation_output$simulated, ylab = 'observed FEC', xlab = 'simulated FEC', pch = 21)
abline(a = 0, b = 1)
abline(a = 0, b = lm_validation_McF$coefficients[1], lty = 2)

dev.off()

# Cooperia simulations =========================================================

# Nematode:
# Choose from:
# 1 = H. contortus
# 2 = T. circumcincta
# 3 = O. ostertagi
# 4 = C. oncophora
nematode = 4

# Enter starting parasite population size and level of acquired immunity...
# Need to estimate adult population based on egg counts
# assuming immunity = 0.5 as they are 10-12 months of age...
# FEC ~100 epg (Fig 1. McFarland et al., 2022, states decrease in initial FECs, so assuming downward trajectory since April = ~100 epg)
ginparms = gin(nematode = nematode, temp = temp, precip = precip, photoperiod = daylength(lat = lat, doy = seq(as.Date(ss), as.Date(ee), "days")))
fecundity_init = exp(ginparms$max.lambda-(ginparms$max.lambda-ginparms$min.lambda)*0.8)
total_eggs_per_day_per_head = 100*f[1]*0.17
est_adults_from_FEC = total_eggs_per_day_per_head/fecundity_init

# Because we now have the init.vals() function above, we only need to specify state variables that should be non-zero. In this case we have newborn lambs and clean pasture so use the defaults for all state variables.
initial_values = init.vals(immunity_host = 0.8, q = 0.22, Adult_in_hostA = est_adults_from_FEC)

# Herd A ####
coopA = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate_A, graze = grazing_plots_A,
                    movements = movement_dates_A, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)

# plot.ts(coopA[,1:10])
# plot.ts(coopA[,11:20])
# plot.ts(coopA[,21:30])
# plot.ts(coopA[,31:40])
# plot.ts(coopA[,41:50])
# plot.ts(coopA[,51:56])
# 
# par(mfrow = c(1,1))
# boxplot(herdA$FEC*0.17 ~ herdA$day, at = herdA$day[1:10], boxwex = 5, xlim = c(1,365), xaxt = 'n', ylab = 'FEC (epg)') 
# axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
# mtext("Days after turnout", 1, 3)
# lines(coopA[,'FEC'])
# abline(v = 40, lty = 1, col = 'lightgrey')
# text(x = 35, y = 200, col = 'grey', labels = "Doramectin Tx", srt = 90)
# abline(v = 180, lty = 1, col = 'lightgrey')
# text(x = 175, y = 200, col = 'grey', labels = "Housed", srt = 90)

herdA$simulated = coopA[herdA$day, 'FEC']

validation_output = data.frame(day = herdA$day, 
                               simulated = herdA$simulated, 
                               observed = herdA$FEC*0.17,
                               species = rep('Cooperia', length(herdA$day)),
                               study = rep('HerdA', length(herdA$day)))

# Herd B ####

coopB = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate_B, graze = grazing_plots_B,
                    movements = movement_dates_B, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)

# plot.ts(coopB[,1:10])
# plot.ts(coopB[,11:20])
# plot.ts(coopB[,21:30])
# plot.ts(coopB[,31:40])
# plot.ts(coopB[,41:50])
# plot.ts(coopB[,51:56])
# 
# par(mfrow = c(1,1))
# boxplot(herdB$FEC*0.17 ~ herdB$day, at = herdB$day[1:10], boxwex = 5, xlim = c(1,365), xaxt = 'n', ylab = 'FEC (epg)') 
# axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
# mtext("Days after turnout", 1, 3)
# lines(coopB[,'FEC'])
# abline(v = 40, lty = 1, col = 'lightgrey')
# text(x = 35, y = 200, col = 'grey', labels = "Doramectin Tx", srt = 90)
# abline(v = 180, lty = 1, col = 'lightgrey')
# text(x = 175, y = 200, col = 'grey', labels = "Housed", srt = 90)

herdB$simulated = coopB[herdB$day, 'FEC']

# Save output for statistical analysis
validation_output = rbind(validation_output, data.frame(day = herdB$day, 
                                                        simulated = herdB$simulated, 
                                                        observed = herdB$FEC*0.17,
                                                        species = rep('Cooperia', length(herdB$day)),
                                                        study = rep('HerdB', length(herdB$day))))

# Herd C ####

coopC = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate_C, graze = grazing_plots_C,
                    movements = movement_dates_C, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)

# plot.ts(coopC[,1:10])
# plot.ts(coopC[,11:20])
# plot.ts(coopC[,21:30])
# plot.ts(coopC[,31:40])
# plot.ts(coopC[,41:50])
# plot.ts(coopC[,51:56])
# 
# par(mfrow = c(1,1))
# boxplot(herdC$FEC*0.17 ~ herdC$day, at = herdC$day[1:10], boxwex = 5, xlim = c(1,365), xaxt = 'n', ylab = 'FEC (epg)') 
# axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
# mtext("Days after turnout", 1, 3)
# lines(coopC[,'FEC'])
# abline(v = 40, lty = 1, col = 'lightgrey')
# text(x = 35, y = 100, col = 'grey', labels = "Doramectin Tx", srt = 90)
# abline(v = 180, lty = 1, col = 'lightgrey')
# text(x = 175, y = 100, col = 'grey', labels = "Housed", srt = 90)

herdC$simulated = coopC[herdC$day, 'FEC']

# Save output for statistical analysis
validation_output = rbind(validation_output, data.frame(day = herdC$day, 
                                                        simulated = herdC$simulated, 
                                                        observed = herdC$FEC*0.17,
                                                        species = rep('Cooperia', length(herdC$day)),
                                                        study = rep('HerdC', length(herdC$day))))

# Herd D ####

coopD = gloworm_meta(start = ss, end = ee, lat = lat, 
                    temp = temp, precip = precip, 
                    statevars = initial_values, host = host, nematode = nematode,
                    stocking_rate = stocking_rate_D, graze = grazing_plots_D,
                    movements = movement_dates_D, kgDMha = kgDMha, DMI = DMI,
                    lwt = lw, faeces = f, eventdat = eventdat)

# plot.ts(coopD[,1:10])
# plot.ts(coopD[,11:20])
# plot.ts(coopD[,21:30])
# plot.ts(coopD[,31:40])
# plot.ts(coopD[,41:50])
# plot.ts(coopD[,51:56])
# 
# par(mfrow = c(1,1))
# boxplot(herdD$FEC*0.17 ~ herdD$day, at = herdD$day[1:10], boxwex = 5, xlim = c(1,365), xaxt = 'n', ylab = 'FEC (epg)') 
# axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
# mtext("Days after turnout", 1, 3)
# lines(coopD[,'FEC'])
# abline(v = 40, lty = 1, col = 'lightgrey')
# text(x = 35, y = 200, col = 'grey', labels = "Doramectin Tx", srt = 90)
# abline(v = 180, lty = 1, col = 'lightgrey')
# text(x = 175, y = 200, col = 'grey', labels = "Housed", srt = 90)

herdD$simulated = coopD[herdD$day, 'FEC']

# Save output for statistical analysis
validation_output = rbind(validation_output, data.frame(day = herdD$day, 
                                                        simulated = herdD$simulated, 
                                                        observed = herdD$FEC*0.17,
                                                        species = rep('Cooperia', length(herdD$day)),
                                                        study = rep('HerdD', length(herdD$day))))
write.csv(validation_output, file = "Coop_McFarland_validation_output.csv")

# Stats ####
lm_validation_McF = lm(validation_output$observed.mean ~ 0 + validation_output$simulated)
summary(lm_validation_McF)
plot(lm_validation_McF) # some deviance from normality in the residuals. This would be because the FEC data are likely to be negative binomially distributed. 

plot(validation_output$observed.mean ~ validation_output$simulated, ylab = 'observed FEC', xlab = 'simulated FEC', pch = 21)
abline(a = 0, b = 1)
abline(a = 0, b = lm_validation_McF$coefficients[1], lty = 2)

capture.output(summary(lm_validation_McF), file = 'Coop_McFarland_stats.txt')
cor_validatoin_McF = cor.test(validation_output$simulated, validation_output$observed.mean, method = "spearman")
capture.output(cor_validatoin_McF, file = 'Coop_McFarland_stats.txt', append = T)

# Final plots ####
dev.off()
tiff('Coop_McFarland_longitudinal.tiff', width = 8, height = 8, units = 'in', res = 300)
par(mfrow = c(2,2))
boxplot(herdA$FEC*0.17 ~ herdA$day, at = herdA$day[1:10], boxwex = 10, xlim = c(1,365), ylim = c(0, 100), xaxt = 'n', ylab = 'FEC (epg)', xlab = "Days after turnout") 
axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
lines(coopA[,'FEC'])
abline(v = 40, lty = 1, col = 'lightgrey')
text(x = 30, y = 80, col = 'grey', labels = "Doramectin Tx", srt = 90)
abline(v = 180, lty = 1, col = 'lightgrey')
text(x = 170, y = 90, col = 'grey', labels = "Housed", srt = 90)
text(x = 300, y = 90, labels = "Herd A")

boxplot(herdB$FEC*0.17 ~ herdB$day, at = herdB$day[1:10], boxwex = 10, xlim = c(1,365), ylim = c(0, 100), xaxt = 'n', ylab = 'FEC (epg)', xlab = "Days after turnout") 
axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
lines(coopB[,'FEC'])
abline(v = 40, lty = 1, col = 'lightgrey')
text(x = 30, y = 80, col = 'grey', labels = "Doramectin Tx", srt = 90)
abline(v = 180, lty = 1, col = 'lightgrey')
text(x = 170, y = 90, col = 'grey', labels = "Housed", srt = 90)
text(x = 300, y = 90, labels = "Herd B")

boxplot(herdC$FEC*0.17 ~ herdC$day, at = herdC$day[1:10], boxwex = 10, xlim = c(1,365), ylim = c(0, 100), xaxt = 'n', ylab = 'FEC (epg)', xlab = "Days after turnout") 
axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
lines(coopC[,'FEC'])
abline(v = 40, lty = 1, col = 'lightgrey')
text(x = 30, y = 80, col = 'grey', labels = "Doramectin Tx", srt = 90)
abline(v = 180, lty = 1, col = 'lightgrey')
text(x = 170, y = 90, col = 'grey', labels = "Housed", srt = 90)
text(x = 300, y = 90, labels = "Herd C")

boxplot(herdD$FEC*0.17 ~ herdD$day, at = herdD$day[1:10], boxwex = 10, xlim = c(1,365), ylim = c(0, 100), xaxt = 'n', ylab = 'FEC (epg)', xlab = "Days after turnout") 
axis(1, at = c(1, 50, 100, 150, 200, 250, 300, 350), labels = c(1, 50, 100, 150, 200, 250, 300, 350))
lines(coopD[,'FEC'])
abline(v = 40, lty = 1, col = 'lightgrey')
text(x = 30, y = 80, col = 'grey', labels = "Doramectin Tx", srt = 90)
abline(v = 180, lty = 1, col = 'lightgrey')
text(x = 170, y = 90, col = 'grey', labels = "Housed", srt = 90)
text(x = 300, y = 90, labels = "Herd D")

dev.off()

tiff('Coop_McFarland_sim_obs.tiff', width = 5, height = 5, units = 'in', res = 300)

plot(validation_output$observed.mean ~ validation_output$simulated, ylab = 'observed FEC', xlab = 'simulated FEC', pch = 21)
abline(a = 0, b = 1)
abline(a = 0, b = lm_validation_McF$coefficients[1], lty = 2)

dev.off()
