# GLOWORM-META model

### A model to simulate the weather-dependent population dynamics of ruminant gastrointestinal nematodes, incorporating pasture metapopulations
## Currently parameterised for cattle GINs only

## Model phylogeny

GLOWORM-FL[^1] simulated the weather-dependent population dynamics of the free-living stages of ruminant gastrointestinal nematodes (GINs) and was parameterised and validated for *Haemonchus contortus, Teladorsagia circumcincta* and *Ostertagia ostertagi* in a temperate climate (UK). The model incorporated the behaviour of nematodes in pasture, resulting in an improved fit to field data than previously published models.  
Model parameters were further adapted for *O. ostertagi* and *Cooperia oncophora* and validated in Calgary, Canada [^2].

GLOWORM-PARA[^3] simulated the population dynamics of the parasitic stages of ruminant GINs, and was parameterised and validated for *O. ostertagi* and *C. oncophora* in first-season grazing calves in Belgium. 

GLOWORM-META[^4] combines GLOWORM-FL and GLOWORM-PARA and extended the model to include a pasture metapopulation. The model was further parameterised and validated for *O. ostertagi* and *C. oncophora* in first and second-season grazers in Ireland, Belgium and Scotland, managed in both set-stocking and rotational grazing systems. 

## Running the model in R

**Install and load necessary packages:**

    packages = c("deSolve", "geosphere", "chron", "ncdf4", "imputeTS", "forecast", "gplots", "tidyr", "rstudioapi", "openmeteo") # Package names
  
    installed_packages = packages %in% rownames(installed.packages())
    if (any(installed_packages == FALSE)) {
      install.packages(packages[!installed_packages])
    }
    
    invisible(lapply(packages, library, character.only = TRUE)) # Load packages

**Source model functions:**

    source('Model_Functions/livestockfuns.R') # liveweight, faeces production and dry matter intake
    source('Model_Functions/weatherfuns.R') # EOBS data extraction function
    source('Model_Functions/parainterp.R') # Interpolation function for parasitological input
    source('Model_Functions/gloworm_meta.R') # GLOWORM-META function
    source('Model_Functions/initialvalues.r') # sets initial conditions
    source('Model_Functions/ginparms.r') # species-specific parameters

**Enter start and end dates, and location (lat/lon in decimal degrees):**

    ss = '2020-05-01'
    ee = '2020-12-31'
    lat = 55.8573407
    lon = -3.1981503

**Import weather data:**
Use this to import daily mean temperature and rainfall data for the simulations. If you're using E-OBS gridded data, you may use the *eobspoint* function sourced above. Otherwise, source demo data from OpenMeteo:

    wdat = weather_history(location = c(lat, lon), start = ss, end = ee, daily = list("temperature_2m_max", "temperature_2m_min", "precipitation_sum"))
    precip = wdat$daily_precipitation_sum
    temp = apply(X = wdat[,2:3], MARGIN = 1, FUN = mean)

**Specify host and GIN species for simulation:**

    host = 2 # 1 = ovine, 2 = bovine
    nematode = 3 # 1 = *H. contortus*, 2 = *T. circumcincta*, 3 = *O. ostertagi*, 4 = *C. oncphora*

**Specify farm management data:**
    
    kgDM = 2000

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

Generate vectors of daily data for stocking rates (hosts per hectare), movement dates (must include the start and end dates of the simulations), and grazing plots.
Plot 0 = zero-grazing (housed), plots 1-6 = pastures 1-6.
In this example, stocking rate is 15 hosts per hectare.
Hosts begin on pasture 0 (housed), move to pasture 1 on day 7, and to pasture 2 on day 120, where they remain there until the end of the simulation.
    
    stocking_rate = c(rep(15, length(seq.Date(as.Date(ss), as.Date(ee), 'day'))))
    movement_dates = c(ss, as.character(as.Date(ss)+7), as.character(as.Date(ss)+120), ee)  #must include the start and end dates - no movements in this example
    grazing_plots = c(0,1,2,2)

**Generate initial conditions for the state variables:**

Use the init.vals function sourced above to specify numbers of individuals of each life cycle stage on each pasture and in the host.
Options and defaults are described in "Model_Functions/initialvalues.r". 
    
    initial_values = init.vals(immunity_host = 0.01, Preadult_in_host = 2000, Adult_in_hostA = 8000)
  

**Add events:**

Treatments are simulated as events. You must generate a data frame with a column defining the state variable affected, a column indicating the time step of the event, a column indicating the value, and a column indicating the method of the event. 

    # No Quarantine, No Treatment, animals straight to Pasture 1 ####
    eventdat = NULL
    eventdat

In the example below taken from the gloworm-ar documentation, we are applying a treatment to the lambs every 28 days, which removes the pre-adult (Pa), hypobiotic (P) and adult (A) nematodes in the host, according to their genotype and the percentage reductions expected for each genotype which we specified above. The method "mult" multiplies each state variable by the specified value at the specified time. This code is not run here, and is included for information only.

    '''eventdat_SS = data.frame(
      var = c(rep('Pa_B_SS', 1), rep('P_B_SS', 1), rep('A_B_SS', 1)), 
      time = rep(141, 3), 
      value = rep(1-reduction_SS, 1*3), 
      method = rep('mult', 1*3))
    eventdat_RS = data.frame(
      var = c(rep('Pa_B_RS', 1), rep('P_B_RS', 1), rep('A_B_RS', 1)), 
      time = rep(141, 3), 
      value = rep(1-reduction_RS, 1*3), 
      method = rep('mult', 1*3))
    eventdat_RR = data.frame(
      var = c(rep('Pa_B_RR', 1), rep('P_B_RR', 1), rep('A_B_RR', 1)), 
      time = rep(141, 3), 
      value = rep(1-reduction_RR, 1*3), 
      method = rep('mult', 1*3))
    eventdat= rbind(eventdat_RR, eventdat_RS, eventdat_SS) # bind the three subpopn events
    eventdat = eventdat[order(eventdat$time),] # prevents warnings after sims
    eventdat'''

**Run the simulation:**

    sim = gloworm_meta(start = ss, end = ee, lat = lat, 
                                 temp = temp, precip = precip, 
                                 statevars = initial_values, host = host, nematode = nematode,
                                 stocking_rate = stocking_rate, graze = grazing_plots,
                                 movements = movement_dates, kgDMha = kgDMha, DMI = DMI,
                                 lwt = lw, faeces = f, eventdat = eventdat)
            
[^1]: Rose et al., 2015. https://doi.org/10.1016/j.ecolmodel.2014.11.033
[^2]: Wang et al., 2022. https://doi.org/10.1016/j.vetpar.2022.109777
[^3]: Rose Vineer et al., 2020. https://doi.org/10.1016/j.ijpara.2019.11.005
[^4]: Ingle et al., 2025 (pending submission - will be available in preprint early 2025)
[^5]: Leathwick et al., 1995. https://doi.org/10.1016/0020-7519(95)00059-3
