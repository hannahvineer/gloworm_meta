#' Get GI nematode life history parameters for the GLOWORM models
#'
#' This function defines the life history parameters needed as input for the gloworm models, dependent on nematode species and environmental conditions.
#' This function can be used alone to extract the parameters, but it is also called by the GLOWORM model functions, and therefore you do not have to use this function prior to running the models.
#'
#' @param nematode 1 = Haemonchus contortus, 2 = Teladorsagia circumcinta, 3 = Ostertagia ostertagi, 4 = Cooperia oncophora
#' @param temp daily time series of mean temperature. See ?gloworm::eobspoint for a possible method to obtain these data
#' @param precip daily time series of total rainfall. See ?gloworm::eobspoint for a possible method to obtain these data
#' @param photoperiod daily time series of daylight hours. See ?geosphere::daylength for a possible method to obtain these data
#' @return Data frame with species-specific life history parameters, to be used as input for the GLOWORM model functions
#' @examples
#' gin(nematode = 2, temp = rep(10, 10), precip = rep(5, 10), photoperiod = rep(12, 10))
#' @importFrom imputeTS na.ma
#' @importFrom forecast ma
#' @export

gin = function(nematode, temp, precip, photoperiod) {


  if(nematode==1) {
    # H. contortus [COMING SOON]

    (sapply(ls(),function(x)get(x),simplify=F,USE.NAMES=T))

  } else {


    if(nematode==2) {
      # T. circumcincta [COMING SOON]

      (sapply(ls(),function(x)get(x),simplify=F,USE.NAMES=T))


    } else {


      if(nematode==3) {
        # O ostertagia

        dev.1 = pmin(1, pmax(0, -0.07258 + 0.00976*temp))
        mu.1 = pmin(1, exp(-4.38278 -0.10640*temp + 0.00540*(temp^2)))
        mu.2 = pmin(1, exp(-4.38278 -0.10640*temp + 0.00540*(temp^2)))
        mu.4 = pmin(1, exp(-5.743  - 0.1755 *temp + 0.0068 *(temp^2) -0.00003*(temp^3))) # pmin(1, exp(-6.388  - 0.2681 *temp + 0.01633 *(temp^2) -0.00016*(temp^3)))
        mu.3 = pmin(1, exp(-4.864  - 0.206 *temp + 0.0048 *(temp^2) +0.00008*(temp^3))) # pmin(1, mu.4*10)
        mu.5 = mu.3
        h.mig = ifelse(precip>=5, pmax(0, pmin(1, exp(-3.29608 + 0.04594*precip))),0)
        v.mig = (pmax(0, exp(-5.48240 + 0.45392*temp - 0.01252*(temp^2))*0.268))

        egg.correction = rep(1, length(temp))

        dev.2 = 0.04077336
        max.mu6 = -log(0.23)/17 # PPP of 17 days, low establishment rate of 0.05% in Verschave
        min.mu6 = 0.054 # -log(0.27)/17 # PPP 17 days, mean establishment rate of 0.27
        mu7 = 0.002
        min.mu8 = 0.028 # mean adult mortality rate verschave
        max.mu8 = 0.08 # high adult mortality rate verschave
        min.h = 0.02
        max.h = 0.3 # min max range in Verschave
        min.lambda = 4.58
        max.lambda = 4.96
        rho.response = 5.981e-5
        sigma.decay = 0.002
        pAf = 0.5

        devsuccess = as.numeric(na.ma(ma(dev.1, order = 7,centre = FALSE), k = 7))
        hypobiosis = max.h - ((max.h-min.h)/(max(devsuccess, na.rm = TRUE)))*(devsuccess) # (max(devsuccess, na.rm = TRUE)-min(devsuccess, na.rm = TRUE)))*(devsuccess-min(devsuccess, na.rm = TRUE)) #max.h - ((max.h-min.h)/(max(photoperiod)-min(photoperiod)))*(photoperiod-min(photoperiod))
        resume = (1/(max(devsuccess, na.rm = TRUE)-min(devsuccess, na.rm = TRUE)))*(devsuccess-min(devsuccess, na.rm = TRUE))
        eggC <- approxfun(x = egg.correction, method = 'constant', rule = 2)
        dev1rate <- approxfun(x = dev.1, method = "linear", rule = 2)
        mu1rate <- approxfun(x = mu.1, method = "linear", rule = 2)
        mu2rate <- approxfun(x = mu.2, method = "linear", rule = 2)
        mu3rate <- approxfun(x = mu.3, method = "linear", rule = 2)
        mu4rate <- approxfun(x = mu.4, method = "linear", rule = 2)
        mu5rate <- approxfun(x = mu.5, method = "linear", rule = 2)
        hmigrate <- approxfun(x = h.mig, method = "linear", rule = 2)
        vmigrate <- approxfun(x = v.mig, method = "linear", rule = 2)
        hrate = approxfun(hypobiosis, method = "linear")
        rrate = approxfun(resume, method = 'linear')

        print('Ostertagia parameters sourced')

        (sapply(ls(),function(x)get(x),simplify=F,USE.NAMES=T))


      } else {

        
        if(nematode==4) {
          #C oncophora
          
          dev.1 = pmax(0, -0.0401 + 0.00821 * temp)
          dev.1 = pmin(1, dev.1)
          mu.1 = pmin(1, exp(-4.38278 - 0.1064 * temp + 0.0054 * (temp^2)))
          mu.2 = pmin(1, exp(-4.38278 - 0.1064 * temp + 0.0054 * (temp^2)))
          mu.4 = pmin(1, exp(-5.822 - 0.1749 * temp + 0.0070 * (temp^2) + 0.00002 * (temp^3)))
          mu.3 = pmin(1, exp(-4.726 - 0.204 * temp + 0.0044 * (temp^2) + 0.00009 * (temp^3)))
          mu.5 = mu.3
          h.mig = ifelse(precip>=6, pmax(0, pmin(1, exp(-3.35002 + 0.05559*precip))),0)
          v.mig = (pmax(0, exp(-5.48240 + 0.45392*temp - 0.01252*(temp^2))*0.459)) ##vmig updated with recovery rate. 
          
          egg.correction = rep(1, length(temp))
          
          dev.2 = 0.04077336
          max.mu6 = 0.12 # high pre-adult mortality in Verschave thesis https://biblio.ugent.be/publication/6986324
          min.mu6 = 0.044 # mean pa mortality in Verschave thesis
          mu7 = 0.002
          min.mu8 = 0.04 # mean Verschave thesis
          max.mu8 = 0.12 # high Verschave thesis
          min.h = 0
          max.h = 0.06
          min.lambda = 6.44
          max.lambda = 7.30 
          rho.response = 1.316e-4
          sigma.decay = 0.002
          pAf = 0.5
          
          devsuccess = as.numeric(na.ma(ma(dev.1, order = 7,centre = FALSE), k = 7))
          hypobiosis = max.h - ((max.h-min.h)/(max(devsuccess, na.rm = TRUE)))*(devsuccess) # (max(devsuccess, na.rm = TRUE)-min(devsuccess, na.rm = TRUE)))*(devsuccess-min(devsuccess, na.rm = TRUE)) #max.h - ((max.h-min.h)/(max(photoperiod)-min(photoperiod)))*(photoperiod-min(photoperiod))
          resume = (1/(max(devsuccess, na.rm = TRUE)-min(devsuccess, na.rm = TRUE)))*(devsuccess-min(devsuccess, na.rm = TRUE))
          eggC <- approxfun(x = egg.correction, method = 'constant', rule = 2)
          dev1rate <- approxfun(x = dev.1, method = "linear", rule = 2)
          mu1rate <- approxfun(x = mu.1, method = "linear", rule = 2)
          mu2rate <- approxfun(x = mu.2, method = "linear", rule = 2)
          mu3rate <- approxfun(x = mu.3, method = "linear", rule = 2)
          mu4rate <- approxfun(x = mu.4, method = "linear", rule = 2)
          mu5rate <- approxfun(x = mu.5, method = "linear", rule = 2)
          hmigrate <- approxfun(x = h.mig, method = "linear", rule = 2)
          vmigrate <- approxfun(x = v.mig, method = "linear", rule = 2)
          hrate = approxfun(hypobiosis, method = "linear")
          rrate = approxfun(resume, method = 'linear')
          
          print('Cooperia parameters sourced')

          (sapply(ls(),function(x)get(x),simplify=F,USE.NAMES=T))

        }
      }
    }
  }
}
