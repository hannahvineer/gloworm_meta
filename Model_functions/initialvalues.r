#' Define initial state variable conditions for gloworm models
#'
#' Defines the initial conditions for the state variables required by the gloworm models.
#' For gloworm(), use only the state variables for plotA and hostA, representing a single pasture and grazed by a single group of hosts.
#' For gloworm_meta(), which allows movements of groups of hosts between pastures, use as many host and plot groups as required, leaving all others at their default values.
#' See ?gloworm, ?gloworm_meta, and ?gloworm_ar. Defaults for all GI nematode state variables = 0. Thus users are only required to enter values for the state variables that are non-zero.
#'
#' Argument descriptions are only shown above for plotA and hostA, but the same definition applies to other plots with subscripts A-F.
#' @param Eggs_plotA Initial number of eggs per hectare on plot A
#' @param L1L2_plotA Initial number of preinfective freeliving larvae (L1 and L2) per hectare on plot A
#' @param L3faeces_plotA Initial number of infective freeliving larvae (L3) in faeces per hectare on plot A
#' @param L3herbage_plotA Initial number of L3 on herbage per hectare in plot A. Note that this is the number of herbage EXCLUDING soil. The number of L3 in soil is estimated from the number on herbage since L3 in soil is rarely measured in practice.
#' @param q frequency of resistance alleles in the population. Necessary only for gloworm::gloworm_ar() simulations.
#' @param Preadult_in_hostA Initial number of ingested L3, L4 and sexually immature adults per host in group A (only one host group currently configured. This is to allow for further development of functions incorporating multiple host groups)
#' @param Arrested_in_hostA Initial number of hypobiotic larvae per host in group A (only one host group currently configured. This is to allow for further development of functions incorporating multiple host groups)
#' @param Adult_in_hostA Initial number of adults per host in group A (only one host group currently configured. This is to allow for further development of functions incorporating multiple host groups)
#' @param immunity_hostA Initial relative level of immunity to the simulated species of GI nematode in host group A. Values between 0.0001 and 1 are accepted. 0.0001 (default) indicates no acquired immunity to infection (only one host group currently configured. This is to allow for further development of functions incorporating multiple host groups)
#' @return Dataframe assigning values to the initial state variables
#' @examples
#' init.vals() # using defaults
#' init.vals(Eggs_plotC = 10e6) # user input value for eggs on pasture C
#' @export
init.vals = function(Eggs_plotA = 0, L1L2_plotA = 0, L3faeces_plotA = 0, L3herbage_plotA = 0,
                     Eggs_plotB = 0, L1L2_plotB = 0, L3faeces_plotB = 0, L3herbage_plotB = 0,
                     Eggs_plotC = 0, L1L2_plotC = 0, L3faeces_plotC = 0, L3herbage_plotC = 0,
                     Eggs_plotD = 0, L1L2_plotD = 0, L3faeces_plotD = 0, L3herbage_plotD = 0,
                     Eggs_plotE = 0, L1L2_plotE = 0, L3faeces_plotE = 0, L3herbage_plotE = 0,
                     Eggs_plotF = 0, L1L2_plotF = 0, L3faeces_plotF = 0, L3herbage_plotF = 0,
                     q = 0.01,
                     Preadult_in_hostA = 0, Arrested_in_hostA = 0, Adult_in_hostA = 0, immunity_hostA = 0.0001,
                     L3pasture_plotA = 0, L3pasture_plotB = 0, L3pasture_plotC = 0, L3pasture_plotD = 0, L3pasture_plotE = 0, L3pasture_plotF = 0) {
  
  .init = cbind(Eggs_plotA, L1L2_plotA, L3faeces_plotA, L3herbage_plotA,
                Eggs_plotB, L1L2_plotB, L3faeces_plotB, L3herbage_plotB,
                Eggs_plotC, L1L2_plotC, L3faeces_plotC, L3herbage_plotC,
                Eggs_plotD, L1L2_plotD, L3faeces_plotD, L3herbage_plotD,
                Eggs_plotE, L1L2_plotE, L3faeces_plotE, L3herbage_plotE,
                Eggs_plotF, L1L2_plotF, L3faeces_plotF, L3herbage_plotF,
                q,
                Preadult_in_hostA, Arrested_in_hostA, Adult_in_hostA, immunity_hostA,
                L3pasture_plotA, L3pasture_plotB, L3pasture_plotC, L3pasture_plotD, L3pasture_plotE, L3pasture_plotF)
  return(as.data.frame(.init))
}
