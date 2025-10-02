
# Author: Hannah Vineer <hannah.vineer@liverpool.ac.uk>
# Last updated: 19/3/24
# Last updated by: Hannah Vineer <hannah.vineer@liverpool.ac.uk>

# 

#' Parasitological data interpolation function
#'
#' Function to linearly interpolate parasitological data to generate daily time series from data collected at less frequent intervals. 
#' 
#' Note that if you are using this function to produce input for the gloworm models, the parasitological data (para) and corresponding dates (dates) must include the start and end date of the model simulation or the model will not run (the parasitological data time series and the simulation will not be equal lengths). 
#' 
#' @param para vector of parasitological data e.g. faecal egg counts
#' @param dates dates data were collected, in the format 'YYYY-MM-DD'
#' @param method "linear" or "constant". See ?approxfun
#' @return Vector of interpolated daily parasitological data
#' @examples 
#' # NEEDS EXAMPLE HERE
#' @export

parainterp = function(para, dates, method="linear") {
  indices = which(seq(as.Date(dates[1]), as.Date(rev(dates)[1]), "days") %in% as.Date(dates))
  days = seq(1, length(seq(as.Date(dates[1]), as.Date(rev(dates)[1]), "days")), 1)
  interpfun = approxfun(y = para, x = indices, method = method)
  return(interpfun(days))
}