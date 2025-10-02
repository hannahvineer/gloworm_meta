#' GLOWORM full life cycle model with pasture metapopulations
#'
#' This model combines the model of the free-living stages on pasture (Rose et al., 2015) and parasitic stages in the host (Rose Vineer et al., 2020) to simulate the complete life cycle of several GI nematode species infecting either cattle or sheep.
#' The free-living module is replicated to allow for 6 metapopulations of free-living stages, representing 6 grazing plots/paddocks/fields.
#' A single group of hosts can move between the 6 pastures to simulate rotational grazing.
#' To use this function to simulate ONLY the free-living stages (representing GLOWORM-FL, Rose et al., 2015), set "graze" to 0 throughout and only enter initial conditions for the free-living stages. To use this function to simulate ONLY the parasitic stages (representing GLOWORM-PARA, Rose Vineer et al., 2020), set "graze" to 0 throughout and only set initial conditions for the parasitic state variables.
#'
#' Rose et al., 2015. GLOWORM-FL: A simulation model of the effects of climate and climate change on the free-living stages of gastro-intestinal nematode parasites of ruminants. Ecological Modelling, 297, 232-245. https://doi.org/10.1016/j.ecolmodel.2014.11.033
#'
#' Rose Vineer et al., 2020. GLOWORM-PARA: a flexible framework to simulate the population dynamics of the parasitic phase of gastrointestinal nematodes infecting grazing livestock, International Journal for Parasitology. 50, 133-144. https://doi.org/10.1016/j.ijpara.2019.11.005
#'
#' @param start start date for simulations in the format "YYYY-MM-DD"
#' @param end end date for simulations in the format "YYYY-MM-DD"
#' @param lat latitude of the simulated site in decimal degrees. This is used to estimate daylength using the geosphere::daylength function
#' @param temp daily time series of mean temperature. See ?gloworm::eobspoint for a possible method to obtain these data
#' @param precip daily time series of total rainfall. See ?gloworm::eobspoint for a possible method to obtain these data
#' @param statevars data frame defining the initial conditions for state variables. Must be produced using gloworm::init.vals(). See ?gloworm::init.vals for details.
#' @param host 1 = sheep, 2 = cattle
#' @param nematode 1 = Haemonchus contortus, 2 = Teladorsagia circumcinta, 3 = Ostertagia ostertagi, 4 = Cooperia oncophora
#' @param stocking_rate Vector of the daily values for stocking rate, individuals per hectare (not LSU)
#' @param graze Vector of values indicating where animals are grazing. Values 1-6 represent grazing plots A-F. Values of 0 can be used to represent housed animals, animals removed from the simulation, or to isolate the free-living and parasitic parts of the life cycle (i.e. prevent pasture contamination and prevent ingestion of L3).
#' @param movements Vector of movement dates in the format "YYYY-MM-DD". Must include the start and end dates of the simulations, which must match ss and ee
#' @param kgDMha Vector of daily values of kg dry matter per hectare of grazing
#' @param DMI Vector of daily dry matter intake values, in kg per day. You may use ?gloworm::dmi to produce this vector if empirical data are not available
#' @param lwt Vector of daily liveweight values, in kg. You may use ?gloworm::lw_cattle or ?gloworm::lw_sheep to produce this vector if empirical data are not available
#' @param faeces Vector of daily faeces production per host, in grams per day. You may use ?gloworm::faeces to produce this vector if empirical data are not available
#' @param eventdat Event dataframe. See ?deSolve::lsoda for details. Note that available state variables for the GLOWORM model are: E, L, L3f, L3p (with the suffixes "_A", "_B", "_C", "_D", "_E" and "_F" for the different grazing plots), P, Pa, A and r (eggs, L1/L2 in faeces, L3 in faeces, L3 on pasture, Preadult parasitic stages, Arrested preadult parasitic stages, Adult, and host immunity/resistance).
#' @return Data frame with daily simulated values for each state variable.
#' @examples
#' # Not run
#' ss = "2010-01-01"
#' ee = "2010-12-31"
#' lat = 56
#' lon = -4.15
#' temp = gloworm::eobspoint(data = 'Weather_data/tg_ens_mean_0.1deg_reg_1995-2010_v29.0e.nc', var = 'tg', lat = lat, lon = lon, start = ss, end = ee)
#' precip = gloworm::eobspoint(data = 'Weather_data/rr_ens_mean_0.1deg_reg_1995-2010_v29.0e.nc', var = 'rr', lat = lat, lon = lon, start = ss, end = ee)
#' stocking_rate = rep(50, length(temp))
#' grazing_plots = c(1,2,3) # these sheep move from grazing plot A, to B, to C, represented by 1,2,3 respectively.
#' movement_dates = c(ss, "2010-06-01", ee)
#' host = 1 # sheep
#' host_age = -90 # simulation starts 90 days before lambing
#' kgDMha = rep(1200, length(temp))
#' lw = if(host == 1) {
#'    lw_sheep(age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1))
#'      } else {
#'        lw_cattle(age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1))
#'        }
#' f = faeces(lwt = lw, host_species = host)
#' DMI = dmi(lwt = lw, age = seq(host_age, host_age + as.numeric(as.Date(ee)-as.Date(ss)), 1), host_species = host)
#'
#' # Add events: periparturient rise from day of lambing
#' # Simple representation of 500 epg, 1kg faecal output per day
#' eventdat = data.frame(
#'   var = rep('E_A', 85),
#'   time = seq(90, 90+84, 1),
#'   value = parainterp(para = c(500*1000*stocking_rate[1], 0), dates = c(ss, as.character(as.Date(ss)+84)), method = 'linear'),
#'   method = rep('add', 85))
#'
#' # Run the model simulation
#' gloworm_meta(start = ss, end = ee, lat = lat, temp = temp, precip = precip,
#'        statevars = init.vals(), host = host, nematode = 1,
#'        stocking_rate = stocking_rate, graze = grazing_plots, movements = movement_dates,
#'        kgDMha = kgDMha, DMI = DMI, lwt = lw, faeces = f, eventdat = eventdat)
#' @importFrom deSolve lsoda
#' @importFrom geosphere daylength
#' @importFrom imputeTS na.ma
#' @importFrom forecast ma
#' @export

gloworm_meta = function(start, end, lat, temp, precip, statevars,
                   host, nematode,
                   stocking_rate, graze, movements, kgDMha,
                   DMI, lwt, faeces,
                   eventdat = NULL) {

  # This code checks the lengths of environmental variables are correct and gets the species-specific parameters
  date.range = seq(as.Date(start), as.Date(end), "days")
  global.t = seq(1, length(date.range))
  photoperiod = daylength(lat = lat, doy = date.range)

  # Check if weather data match length of simulation
  if (length(date.range) != length(temp))
    stop("Error: length of temperature data and dates do not match. Ensure climatic data correspond to the start and end dates (one entry per day)")
  if (length(date.range) != length(precip))
    stop("Error: length of rainfall data and dates do not match. Ensure climatic data correspond to the start and end dates (one entry per day)")
  if (length(date.range) != length(stocking_rate))
    stop("Error: length of stocking_rate variable and dates do not match. Ensure rates (animals per hectare) correspond to the start and end dates (one entry per day)")
  if (length(date.range) != length(kgDMha))
    stop("Error: length of kgDMha (herbage biomass) variable and dates do not match. Ensure biomass corresponds to the start and end dates (one entry per day)")

  # The dry matter intake is converted to the proportion of available dry matter consumed each day (based on the user input value) and then to the daily instantaneous intake rate.
  DMI[is.nan(DMI)] = 0
  pDMI = (DMI/kgDMha)
  rDMI = -log(1-pDMI)

  # Get GI nematode parameters (constant and forced variables)
  ginparms = gin(nematode = nematode, temp = temp, precip = precip, photoperiod = photoperiod)

  # Interpolation functions for external forcing
  dry_matter_intake = approxfun(rDMI, method = "linear")
  # Create interpolation functions for the grazing locations
  g.interp = approxfun(y = graze, x = c(which(as.character(date.range) %in% as.character(movements))), method = 'constant')
  g = g.interp(seq(1, length(date.range)))
  grazeA = approxfun(ifelse(g==1, 1, 0), method = "constant")
  grazeB = approxfun(ifelse(g==2, 1, 0), method = "constant")
  grazeC = approxfun(ifelse(g==3, 1, 0), method = "constant")
  grazeD = approxfun(ifelse(g==4, 1, 0), method = "constant")
  grazeE = approxfun(ifelse(g==5, 1, 0), method = "constant")
  grazeF = approxfun(ifelse(g==6, 1, 0), method = "constant")
  density = approxfun(stocking_rate, method = 'constant')
  eggC <- approxfun(x = ginparms$egg.correction, method = 'constant', rule = 2)
  dev1rate <- approxfun(x = ginparms$dev.1, method = "linear", rule = 2)
  mu1rate <- approxfun(x = ginparms$mu.1, method = "linear", rule = 2)
  mu2rate <- approxfun(x = ginparms$mu.2, method = "linear", rule = 2)
  mu3rate <- approxfun(x = ginparms$mu.3, method = "linear", rule = 2)
  mu4rate <- approxfun(x = ginparms$mu.4, method = "linear", rule = 2)
  mu5rate <- approxfun(x = ginparms$mu.5, method = "linear", rule = 2)
  hmigrate <- approxfun(x = ginparms$h.mig, method = "linear", rule = 2)
  vmigrate <- approxfun(x = ginparms$v.mig, method = "linear", rule = 2)
  hrate = approxfun(ginparms$hypobiosis, method = "linear")
  rrate = approxfun(ginparms$resume, method = 'linear')

  parms = c(rho = ginparms$rho.response,
            sigma = ginparms$sigma.decay,
            dev2 = ginparms$dev.2,
            maxmu6 = ginparms$max.mu6,
            minmu6 = ginparms$min.mu6,
            mu7 = ginparms$mu7,
            maxmu8 = ginparms$max.mu8,
            minmu8 = ginparms$min.mu8,
            minlambda = ginparms$min.lambda,
            maxlambda = ginparms$max.lambda)

# Make sure you apply the egg correction (dessiccation) to eggs input as events too!
for (i in 1:length(eventdat$value[which(eventdat$var=='E_A')])) {
  eventdat$value[which(eventdat$var=='E_A')][i] = ginparms$egg.correction[which(eventdat$var=='E_A')][i]*eventdat$value[which(eventdat$var=='E_A')][i]
}
for (i in 1:length(eventdat$value[which(eventdat$var=='E_B')])) {
  eventdat$value[which(eventdat$var=='E_B')][i] = ginparms$egg.correction[which(eventdat$var=='E_B')][i]*eventdat$value[which(eventdat$var=='E_B')][i]
}
for (i in 1:length(eventdat$value[which(eventdat$var=='E_C')])) {
  eventdat$value[which(eventdat$var=='E_C')][i] = ginparms$egg.correction[which(eventdat$var=='E_C')][i]*eventdat$value[which(eventdat$var=='E_C')][i]
}
for (i in 1:length(eventdat$value[which(eventdat$var=='E_D')])) {
  eventdat$value[which(eventdat$var=='E_D')][i] = ginparms$egg.correction[which(eventdat$var=='E_D')][i]*eventdat$value[which(eventdat$var=='E_D')][i]
}
for (i in 1:length(eventdat$value[which(eventdat$var=='E_E')])) {
  eventdat$value[which(eventdat$var=='E_E')][i] = ginparms$egg.correction[which(eventdat$var=='E_E')][i]*eventdat$value[which(eventdat$var=='E_E')][i]
}
for (i in 1:length(eventdat$value[which(eventdat$var=='E_F')])) {
  eventdat$value[which(eventdat$var=='E_F')][i] = ginparms$egg.correction[which(eventdat$var=='E_F')][i]*eventdat$value[which(eventdat$var=='E_F')][i]
}

  # Get initial values

  y = c(E_A = statevars$Eggs_plotA * ginparms$egg.correction[1],
        L_A = statevars$L1L2_plotA,
        L3f_A = statevars$L3faeces_plotA,
        L3p_A = statevars$L3herbage_plotA * (1/ginparms$v.mig[1]),
        E_B = statevars$Eggs_plotB * ginparms$egg.correction[1],
        L_B = statevars$L1L2_plotB,
        L3f_B = statevars$L3faeces_plotB,
        L3p_B = statevars$L3herbage_plotB * (1/ginparms$v.mig[1]),
        E_C = statevars$Eggs_plotC * ginparms$egg.correction[1],
        L_C = statevars$L1L2_plotC,
        L3f_C = statevars$L3faeces_plotC,
        L3p_C = statevars$L3herbage_plotC * (1/ginparms$v.mig[1]),
        E_D = statevars$Eggs_plotD * ginparms$egg.correction[1],
        L_D = statevars$L1L2_plotD,
        L3f_D = statevars$L3faeces_plotD,
        L3p_D = statevars$L3herbage_plotD * (1/ginparms$v.mig[1]),
        E_E = statevars$Eggs_plotE * ginparms$egg.correction[1],
        L_E = statevars$L1L2_plotE,
        L3f_E = statevars$L3faeces_plotE,
        L3p_E = statevars$L3herbage_plotE * (1/ginparms$v.mig[1]),
        E_F = statevars$Eggs_plotF * ginparms$egg.correction[1],
        L_F = statevars$L1L2_plotF,
        L3f_F = statevars$L3faeces_plotF,
        L3p_F = statevars$L3herbage_plotF * (1/ginparms$v.mig[1]),
        P = statevars$Preadult_in_hostA,
        Pa = statevars$Arrested_in_hostA,
        A = statevars$Adult_in_hostA,
        r = statevars$immunity_hostA)

  gloworm_meta_mod = function (t, y, parms) {

    with(as.list(c(y, parms)), {

      # specify the parameter to take from the external forcing
      # climate dependant rates
      dev1 = dev1rate(t)
      mu1 = mu1rate(t)
      mu2 = mu2rate(t)
      mu3 = mu3rate(t)
      mu4 = mu4rate(t)
      mu5 = mu5rate(t)
      m1 = hmigrate(t)
      m2 = vmigrate(t)
      correction = eggC(t)
      # Other seasonal forcing
      h = hrate(t)    # hypobiosis
      h2 = rrate(t) # resumed development
      intake = dry_matter_intake(t)  # mean dry matter intake (kgDM)
      stock.rate = density(t)
      gA = grazeA(t)
      gB = grazeB(t)
      gC = grazeC(t)
      gD = grazeD(t)
      gE = grazeE(t)
      gF = grazeF(t)

      # Estimate immune-dependent parameters
      mu6 = minmu6+(maxmu6-minmu6)*r
      mu8 = minmu8+(maxmu8-minmu8)*r
      lambda = maxlambda+(minlambda-maxlambda)*r
      decay = sigma

      # Calculate the derivatives
      # Freeliving (per unit area)

      #Plot A
      dE_A = -(dev1 * 2 + mu1) * E_A + lambda*A*correction*stock.rate*gA
      dL_A = -(dev1 * 2 + mu2) * L_A + (dev1 * 2) * E_A
      dL3f_A = -(mu3 + m1) * L3f_A + (dev1 * 2) * L_A
      dL3p_A = -mu4 * (L3p_A * (1 - m2)) - mu5 * (L3p_A * m2) + m1 * L3f_A - (L3p_A * m2)*intake*stock.rate*gA

      #Plot B
      dE_B = -(dev1 * 2 + mu1) * E_B + lambda*A*correction*stock.rate*gB
      dL_B = -(dev1 * 2 + mu2) * L_B + (dev1 * 2) * E_B
      dL3f_B = -(mu3 + m1) * L3f_B + (dev1 * 2) * L_B
      dL3p_B = -mu4 * (L3p_B * (1 - m2)) - mu5 * (L3p_B * m2) + m1 * L3f_B - (L3p_B * m2)*intake*stock.rate*gB

      #Plot C
      dE_C = -(dev1 * 2 + mu1) * E_C + lambda*A*correction*stock.rate*gC
      dL_C = -(dev1 * 2 + mu2) * L_C + (dev1 * 2) * E_C
      dL3f_C = -(mu3 + m1) * L3f_C + (dev1 * 2) * L_C
      dL3p_C = -mu4 * (L3p_C * (1 - m2)) - mu5 * (L3p_C * m2) + m1 * L3f_C - (L3p_C * m2)*intake*stock.rate*gC

      #Plot D
      dE_D = -(dev1 * 2 + mu1) * E_D + lambda*A*correction*stock.rate*gD
      dL_D = -(dev1 * 2 + mu2) * L_D + (dev1 * 2) * E_D
      dL3f_D = -(mu3 + m1) * L3f_D + (dev1 * 2) * L_D
      dL3p_D = -mu4 * (L3p_D * (1 - m2)) - mu5 * (L3p_D * m2) + m1 * L3f_D - (L3p_D * m2)*intake*stock.rate*gD

      #Plot E
      dE_E = -(dev1 * 2 + mu1) * E_E + lambda*A*correction*stock.rate*gE
      dL_E = -(dev1 * 2 + mu2) * L_E + (dev1 * 2) * E_E
      dL3f_E = -(mu3 + m1) * L3f_E + (dev1 * 2) * L_E
      dL3p_E = -mu4 * (L3p_E * (1 - m2)) - mu5 * (L3p_E * m2) + m1 * L3f_E - (L3p_E * m2)*intake*stock.rate*gE

      #Plot F
      dE_F = -(dev1 * 2 + mu1) * E_F + lambda*A*correction*stock.rate*gF
      dL_F = -(dev1 * 2 + mu2) * L_F + (dev1 * 2) * E_F
      dL3f_F = -(mu3 + m1) * L3f_F + (dev1 * 2) * L_F
      dL3p_F = -mu4 * (L3p_F * (1 - m2)) - mu5 * (L3p_F * m2) + m1 * L3f_F - (L3p_F * m2)*intake*stock.rate*gF

      # Parasitic (per host)

      dP = - dev2*P - mu6*P + L3p_A*m2*intake*gA + L3p_B*m2*intake*gB + L3p_C*m2*intake*gC + L3p_D*m2*intake*gD + L3p_E*m2*intake*gE + L3p_F*m2*intake*gF
      dPa = dev2*h*P  - (h2+mu7)*Pa
      dA = dev2*(1-h)*P  + h2*Pa - mu8*A

      # Immunity
      dr = -r*decay + rho*(L3p_A*m2*intake*gA + L3p_B*m2*intake*gB + L3p_C*m2*intake*gC + L3p_D*m2*intake*gD + L3p_E*m2*intake*gE + L3p_F*m2*intake*gF)*(1-r)

      return(list(c(dE_A = dE_A, dL_A = dL_A, dL3f_A = dL3f_A, dL3p_A = dL3p_A,
                    dE_B = dE_B, dL_B = dL_B, dL3f_B = dL3f_B, dL3p_B = dL3p_B,
                    dE_C = dE_C, dL_C = dL_C, dL3f_C = dL3f_C, dL3p_C = dL3p_C,
                    dE_D = dE_D, dL_D = dL_D, dL3f_D = dL3f_D, dL3p_D = dL3p_D,
                    dE_E = dE_E, dL_E = dL_E, dL3f_E = dL3f_E, dL3p_E = dL3p_E,
                    dE_F = dE_F, dL_F = dL_F, dL3f_F = dL3f_F, dL3p_F = dL3p_F,
                    dP = dP, dPa =dPa, dA = dA,
                    dr = dr)))
    })
  }


  sol = lsoda(y = y, times = global.t, func = gloworm_meta_mod, parms = parms, events = list(data = eventdat))

  # Post hoc calculations
  fecundity = exp(ginparms$max.lambda-(ginparms$max.lambda-ginparms$min.lambda)*sol[,"r"])
  total.eggs = sol[,"A"]*fecundity
  FEC = total.eggs/faeces
  preadult.mortality = ginparms$min.mu6+(ginparms$max.mu6-ginparms$min.mu6)*sol[,"r"]
  adult.mortality = ginparms$min.mu8+(ginparms$max.mu8-ginparms$min.mu8)*sol[,"r"]
  L3s_A = sol[, "L3p_A"] * (1 - (vmigrate(global.t)))
  L3h_A = sol[, "L3p_A"] * vmigrate(global.t)
  L3_per_kgDM_A = L3h_A/kgDMha
  L3s_B = sol[, "L3p_B"] * (1 - (vmigrate(global.t)))
  L3h_B = sol[, "L3p_B"] * vmigrate(global.t)
  L3_per_kgDM_B = L3h_B/kgDMha
  L3s_C = sol[, "L3p_C"] * (1 - (vmigrate(global.t)))
  L3h_C = sol[, "L3p_C"] * vmigrate(global.t)
  L3_per_kgDM_C = L3h_C/kgDMha
  L3s_D = sol[, "L3p_D"] * (1 - (vmigrate(global.t)))
  L3h_D = sol[, "L3p_D"] * vmigrate(global.t)
  L3_per_kgDM_D = L3h_D/kgDMha
  L3s_E = sol[, "L3p_E"] * (1 - (vmigrate(global.t)))
  L3h_E = sol[, "L3p_E"] * vmigrate(global.t)
  L3_per_kgDM_E = L3h_E/kgDMha
  L3s_F = sol[, "L3p_F"] * (1 - (vmigrate(global.t)))
  L3h_F = sol[, "L3p_F"] * vmigrate(global.t)
  L3_per_kgDM_F = L3h_F/kgDMha
  L3i = ((sol[,'L3p_A']+sol[,'L3p_B']+sol[,'L3p_C']+sol[,'L3p_D']+sol[,'L3p_E']+sol[,'L3p_F']) * ginparms$v.mig)*pDMI*ginparms$grazing.correction

  # Bind all data together
  sol = cbind(sol, L3s_A, L3h_A, L3_per_kgDM_A, L3s_B, L3h_B, L3_per_kgDM_B, L3s_C, L3h_C, L3_per_kgDM_C, L3s_D, L3h_D, L3_per_kgDM_D, L3s_E, L3h_E, L3_per_kgDM_E, L3s_F, L3h_F, L3_per_kgDM_F, L3i, preadult.mortality, adult.mortality, fecundity, total.eggs, FEC, lwt, DMI, faeces)

  return(sol)
}
