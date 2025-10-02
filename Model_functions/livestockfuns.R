#' Age-based liveweight estimates for sheep
#'
#' Estimates the liveweight of sheep based on on logistic growth with a birth weight of 5kg, mature weight of 65kg and daily growth rate of 0.015 which produced 38kg lambs at 7 months of age. These defaults can be adjusted to suit larger/smaller and quicker/slower growing breeds.
#' @param birth birth weight in kg. Default = 5kg
#' @param mature mature weight in kg. Default = 65kg
#' @param r daily growth rate. Default 0.015
#' @param age age in days
#' @return Vector of daily estimated liveweight in kg
#' @examples 
#' lw_sheep(age = c(1:100)) # estimate liveweights from 1-100 days of age
#' @export

lw_sheep = function(birth = 5, mature = 65, r = 0.015, age) {
  
  wt.dyn = function(t, wt.init, wt.par) {
    with(as.list(c(wt.init, wt.par)), {
      
      dW = r*W*(1-(W/K))
      
      return(list(c(dW = dW)))
    })
  }
  wt.init = c(W = birth)
  wt.par = c(K = mature, r = 0.015)
  wt.sol = lsoda(y = wt.init, times = seq(1, 365*6, 1), func = wt.dyn, parms = wt.par)
  colnames(wt.sol) = c('age in days', 'weight in kg')
  
  weight = rep(NA, length(age))
  for(i in 1:length(age)){
    if(age[i]>0) weight[i]=wt.sol[age[i],2] else(weight[i] = 0)
  }
  return(weight)
  
}

#' Age-based liveweight estimates for cattle
#'
#' Estimates the liveweight of based on the fixed effects model defined by Cue et al., (2012) (coefficients in Table 3) for growing heifers up to 24 months of age. After 24 months, liveweight is fixed at the 24 month weight. Browns swiss is the default breed, as it is of intermediate between Holstein and Ayreshire.
#' 
#' Cue et al., 2012. Growth modeling of dairy heifers in Quebec based on random regression. Canadian Journal of Animal Science, 92, 33-47. https://doi.org/10.4141/cjas2011-083
#' @param age age in days
#' @param breed name of breed. Options are "brownswiss", "holstein" and "ayreshire"
#' @return Vector of daily estimated liveweight in kg
#' @examples 
#' lw_cattle(age = c(1:100)) # estimate liveweights from 1-100 days of age, with default model for Brown Swiss breed
#' @export


lw_cattle = function(age, breed = "brownswiss") {
  
  if(breed == "ayreshire") {wt = 35.147 + 23.0215 * (age/30) -0.14053 * (age/30)^2} else {
    if(breed == "brownswiss") {wt = 35.5273 + 24.1564 * (age/30) -0.08871 * (age/30)^2} else {
      if(breed == "holstein") {wt = 36.5146 + 27.3524 * (age/30) -0.1203 * (age/30)^2} else {
        print("Error: check spelling of breed. Should be one of ayreshire, browsnswiss or holstein")}
    } 
  }
  
  if(breed == "ayreshire") {maxwt = 35.147 + 23.0215 * 24 -0.14053 * 24^2} else {
    if(breed == "brownswiss") {maxwt = 35.5273 + 24.1564 * 24 -0.08871 * 24^2} else {
      if(breed == "holstein") {maxwt = 36.5146 + 27.3524 * 24 -0.1203 * 24^2} else {
        print("Error: check spelling of breed. Should be one of ayreshire, browsnswiss or holstein")}
    } 
  }
  for(i in 1:length(age)){
    if(age[i]<1) (wt[i] = 0)
    if(age[i]>720) (wt[i] = maxwt)
  }
  return(wt)
}


#' Liveweight-based dry matter intake estimates for sheep and cattle
#'
#' Estimates the daily dry matter intake for sheep and cattle based on a daily vector of liveweights. 
#' To estimate liveweights see ?lw_sheep or ?lw_cattle. 
#' Cattle intake is based on MAFF (1975) (also see Ingvartsen (1994) if other options are required for DMI estimation) and AHDB estimates of DMI as a percentage of body weight. 
#' For cattle less than 200kg, DMI is 3 percent of liveweight (AHDB). For cattle between 201 and 300kg, DMI is 2.75 percent of liveweight. And for cattle over 301kg, DMI is 2.5 percent liveweight (AHDB) corrected for milk production (MAFF, 1975). 
#' Milk argument indicates milk yield per day (kg). Default milk = 0. Sheep intake is based on Greer et al. (2009) but modified to increase lamb intake rate from zero at birth to weight-based estimate at 51 days of age, 
#' thus accounting for the gradual increase in dry matter intake of suckling lambs.
#' 
#' Greer et al., 2009. Development and field evaluation of a decision support model for anthelmintic treatments as part of a targeted selective treatment (TST) regime in lambs. Veterinary Parasitology, 164, 12-20. https://doi.org/10.1016/j.vetpar.2009.04.017
#' 
#' Ingvartsen 1994. Models of voluntary food intake in cattle. Livestock Production Science, 39, 19-38. https://doi.org/10.1016/0301-6226(94)90149-X 
#' 
#' MAFF 1975. Energy allowances and feeding systems for ruminants. Technical Bulletin 33. https://wellcomecollection.org/works/ey8cqebf/items?canvas=4 
#' @param lwt liveweight in kg
#' @param age age in days
#' @param host_species 1 (sheep) or 2 (cattle)
#' @param milk milk yield per day in kg. Default milk = 0
#' @return Vector of daily estimated dry matter intake in kg
#' @examples 
#' dmi(lwt = lw_cattle(c(1:100)), age = c(1:100), host_species = 2) # uses the lw_cattle() function to estimate lwt
#' @export

dmi = function(lwt, age, host_species, milk = 0) {
  
  if(host_species==1) {
    DMI = 0.0545*(lwt^0.86)
    if(any(age<51)) {
      interp_DMI = approxfun(y = c(0, DMI[which(age==50)[1]]), x = c(1, 50), method = 'linear')
      DMI[which(age>0&age<=50)] = interp_DMI(age[which(age>0&age<=50)])
    }
  }
  
  if(host_species==2) {
    DMI = rep(NA, length(lwt))
    for(i in 1:length(lwt)) {
      if(lwt[i]<=200) {DMI[i]=(0.03*(lwt[i]))} else {
        if(lwt[i]<=300) {DMI[i]=(0.0275*(lwt[i]))} else {
          DMI[i] = 0.025*(lwt[i])+0.1*milk}
      
      }
    }
  }
  
  return(DMI)
  
}



#' Estimated faeces production for cattle and sheep
#'
#' Estimates the amount of faeces produced per day, in grams, for sheep or cattle. Sheep estimate is based on Singleton et al. (2011). Cattle extimate is based on Nennich et al., (2005). (equation 15 for heifers, and 18 for calves). Optimised for dairy cattle. Corrected for the proportion of urine in cattle manure (Verschave et al., 2014, based on Nennich et al., 2005 and Masse et al., 2014).
#' 
#' Masse et al., 2014. Effect of corn dried distiller grains with solubles (DDGS) in dairy cow diets on manure bioenergy production potential. Animals, 4, 82-92. https://doi.org/10.3390/ani4010082 
#' 
#' Nennich et al., 2005. Prediction of Manure and Nutrient Excretion from Dairy Cattle. Journal of Dairy Science, 88, 3721-3733. https://doi.org/10.3168/jds.S0022-0302(05)73058-7
#' 
#' Singleton et al., 2011. A mechanistic model of developing immunity to Teladorsagia circumcincta infection in lambs. Parasitology. 138, 322-332. http://doi.org/10.1017/S0031182010001289
#' 
#' Verschave et al., 2014. The parasitic phase of Ostertagia ostertagi: quantification of the main life history traits through systematic review and meta-analysis. International Journal for Parasitology, 44, 1091-1104. https://doi.org/10.1016/j.ijpara.2014.08.006
#' @param lwt liveweight in kg
#' @param host_species 1 (sheep) or 2 (cattle)
#' @return Vector of daily estimated faeces produced in grams
#' @examples 
#' faeces(lwt = lw_sheep(age = c(1:100)), host_species = 1) # uses lw_sheep() to estimate lwt
#' @export

faeces = function(lwt, host_species) {
  if(host_species==1) { f = lwt*20}
  if(host_species==2) {
    f = ifelse(lwt<=250, (lwt*0.0811)*(1-0.34), (lwt*0.0181+17.8)*(1-0.34))
    f = f*1000
  }
  return(f)
}