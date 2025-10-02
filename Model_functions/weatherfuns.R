#' EOBS data extraction from point location
#'
#' Extracts a time series of daily interpolated weather data from the EOBS gridded dataset based on the given latitude and longitude.
#' EOBS gridded datasets must first be downloaded from ECAD: https://www.ecad.eu/download/ensembles/download.php
#' @param data The path to the .nc file downloaded from ECAD.
#' @param var variable name. See ECAD for details At the time of writing, the options were "tg", "rr", "pp", "tn", "tx", "fg", "hu" and "qq"
#' @param lat latitude in decimal degrees
#' @param lon longitude in decimal degrees
#' @param start start date in the format "YYYY-MM-DD"
#' @param end end date in the format "YYYY-MM-DD"
#' @return Vector of weather data
#' @examples 
#' temperature = eobspoint(data = "rr_ens_mean_0.1deg_reg_1995-2010_v29.0e.nc", var = "tg", start = "2000-01-01", end = "2000-01-31")
#' @import ncdf4
#' @import chron
#' @export

eobspoint = function(data, var, lat, lon, start, end) {
  
  dat = nc_open(data)
  
  latitude <- ncvar_get(dat,"latitude")
  longitude <- ncvar_get(dat,"longitude")  
  
  t <- ncvar_get(dat,"time")
  tunits <- ncatt_get(dat,"time","units")
  
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth = as.integer(unlist(tdstr)[2])
  tday = as.integer(unlist(tdstr)[3])
  tyear = as.integer(unlist(tdstr)[1])
  tchron = as.Date(chron(t, origin = c(tmonth, tday, tyear)))
  
  start. = c(which.min(longitude < lon)-1,which.min(latitude < lat)-1,which(tchron == start))
  count = c(1,1,which(tchron == end)+1-which(tchron == start))
  
  dat2 <- ncvar_get(dat, var, start., count)
  
  nc_close(dat)
  
  return(dat2)
}