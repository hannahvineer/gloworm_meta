# Vertical migration parameter definitions for Ostertagia and Cooperia
# Current parameter definition is for sheep parasites
# New data available from Moredun

# read in Moredun data
dat = read.csv("Parameter_estimation/Moredun_vmig_data.csv")
head(dat)

# Subset by species
dat$Proportion = dat$Proportion*5 # recovery rates may be as low as 20% (Wang et al., 2021 https://doi.org/10.1186/s13071-021-05101-w)
coop = dat[which(dat$Species=="Coop"),]
ost = dat[which(dat$Species=="Ost"),]
min(ost$Proportion)
max(ost$Proportion)
min(coop$Proportion)
max(coop$Proportion)

# Quick visualisation
plot(coop$Temp, coop$Proportion, col = coop$Day)
plot(coop$Temp, coop$Count, col = coop$Day)

plot(ost$Temp, ost$Proportion, col = ost$Day)
plot(coop$Temp, coop$Count, col = coop$Day)

# Are any of the temperatures significantly different?

source("Parameter_estimation/pairwise_ks_test.R")
# The above function was taken directly from the source code of the DataScienceR package as the package wasn't available for this version of R
pairwise_ks_test(coop$Count, as.character(coop$Temp), n_min = 4)
# 20 and 25 degrees are significantly higher than 15 degrees - Cooperia
pairwise_ks_test(ost$Count, as.character(ost$Temp), n_min = 4)
# 25 degrees is significantly higher than 10 degrees - Ostertagia
# Overall, no clear pattern to the data, and data from all groups could plausibly have come from the same distribution, except 20-25 degrees which are significantly higher than some other groups. 

# Vertical migration parameter for sheep species is monotonic, peaking at ~20 degrees
# Logical based on our understanding of metabolics and mortality rates being lowest at around this temperature
x <- seq(1, 40)
y = exp(-5.48240 + 0.45392*x - 0.01252*(x^2))
plot(x,y, type = "l")

# Proportions on herbage are much higher than observed for Cooperia and Ostertagia
points(coop$Temp, coop$Proportion, col = "red")
points(ost$Temp, ost$Proportion, col = "blue")

# Scale this parameter to fit Ostertagia and Cooperia
# But first check whether Coop and Ost data also came from the same distribution
ks.test(ost$Count, coop$Count)
# Significant, therefore evidence to support Coop and Ost having different migration rates

# Minimise sum of squared error to fit the sheep vertical migration model to Coop and Ost
# Code adapted from: https://jootse84.github.io/notes/optim-in-R

sumSqMin = function(prop, temp, par = 1) {
  y = exp(-5.48240 + 0.45392*x - 0.01252*(x^2))*par
  sum((prop - y[temp])^2)
}

coop.scaling.factor <- optim(par = c(0.5), fn = sumSqMin, temp = coop$Temp, prop = coop$Proportion, method = "Brent", lower = 0, upper = 1)

ost.scaling.factor <- optim(par = c(0.5), fn = sumSqMin, temp = ost$Temp, prop = ost$Proportion, method = "Brent", lower = 0, upper = 1)

# Plot the new estimates using the scaling factors
x <- seq(1, 40)
y = exp(-5.48240 + 0.45392*x - 0.01252*(x^2))
plot(x,y, type = "l")
points(coop$Temp, coop$Proportion, col = "red")
points(ost$Temp, ost$Proportion, col = "blue")
lines(x, y*coop.scaling.factor$par, col = "red")
lines(x, y*ost.scaling.factor$par, col = "blue")

# Scaling factors for v.mig:
coop.scaling.factor$par
ost.scaling.factor$par

# Vertical migration rates for initial validation

coop.v.mig = exp(-5.48240 + 0.45392*x - 0.01252*(x^2))*0.459
ost.v.mig = exp(-5.48240 + 0.45392*x - 0.01252*(x^2))*0.268
