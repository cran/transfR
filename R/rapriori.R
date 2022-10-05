#' Net rainfall a priori estimation
#'
#' A priori estimate of net rainfall as required for the inversion.
#' @name rapriori
#' @param Qobs vector of discharge value or object of class \code{transfR}. If no unit is provided,
#' \code{Qobs} is assumed to be in [m3/s].
#' @param area drainage area of the catchment. If no unit is provided, \code{area} is assumed to be in [km2].
#' @param lagtime lag time value of the catchment. If no unit is provided, \code{lagtime} is assumed to be in [h].
#' @param deltat time step of the time series. If no unit is provided, \code{deltat} is assumed to be in [min].
#' @param verbose boolean indicating if information messages should be written to the console
#' @param ... further arguments passed to or from other methods
#' @return An object of the same class of \code{Qobs}. If \code{Qobs} is a transfR object,
#' the same transfR object incremented by the new "RnAp" computed attributes.
#' @import sf stars
#' @importFrom units as_units drop_units
#' @details The function estimates an a priori of the net rainfall from Qobs. It converts Qobs to specific
#' discharge and removes the delay caused by transfer time in the river network (given by \code{lagtime}
#' and that could be estimated from the function \link{lagtime}). If an object of class \code{transfR} is provided,
#' \code{area} is estimated from its \code{st} attribute. Results are stored as a new space-time attribute,
#' called "RnAp", in the \code{transfR} object.
#' @examples
#' data(Oudon)
#' icatch <- 1
#' Qobs <- Oudon$obs[["Qobs"]][,icatch]
#' Qspec <- units::set_units(Qobs/st_area(st_geometry(Oudon$obs)[icatch]), "mm/h")
#' deltat <- units::set_units(1,"h")
#' uc <- velocity(hl = Oudon$hl[[icatch]])
#' uh <- uh(hl = Oudon$hl[[icatch]], uc = uc, deltat = deltat)$prob
#' RnAp <- rapriori(Qobs = Qspec, lagtime = lagtime(hl = Oudon$hl[[icatch]], uc = uc),
#' deltat = deltat)
#' @export
rapriori = function(Qobs, ...) UseMethod("rapriori")

#' @name rapriori
#' @export
rapriori.default = function(Qobs, area, lagtime, deltat, ...){
  warning(paste0("No units provided for the arguments of rapriori(). See the help of rapriori() for assumptions about the units."))
  #--- Assumed units
  Qobs <- units::set_units(Qobs,"m^3/s")
  area <- units::set_units(area,"km^2")
  lagtime <- units::set_units(lagtime,"h")
  deltat <- units::set_units(deltat,"min")
  #--- rapriori
  rapriori(Qobs, area, lagtime, deltat, ...)
}

#' @name rapriori
#' @export
rapriori.units = function(Qobs, area, lagtime, deltat, ...){
  lag  <- units::drop_units(round(lagtime/deltat))
  RnAp <- Qobs[(lag+1):length(Qobs)]
  RnAp[(length(RnAp)+1):(length(RnAp)+max(lag,1))] <- NA
  if(length(units(RnAp/units::set_units(1,"m^3/s"))$numerator)==0 & length(units(RnAp/units::set_units(1,"m^3/s"))$denominator)==0) RnAp <- RnAp/area
  if(deltat==units::set_units(1,"day")){out_units <- units::as_units("mm/d")}else{out_units <- units::as_units("mm/h")}
  return(units::set_units(RnAp[1:length(Qobs)],out_units,mode="standart"))
}

#' @name rapriori
#' @export
rapriori.transfR  <- function(Qobs, verbose=TRUE, ...){
  object <- Qobs
  if(!"lagtime"%in%names(object)){
    if(verbose) cat("Lag times attribute (lagtime) is missing in transfR object (see lagtime()).\nAttempting a default estimate:\n")
    object <- lagtime(object)}
  if(!"deltat"%in%names(object)) stop("Time step attribute (deltat) is missing in transfR object. See as_transfr().")
  if(!"Qobs"%in%names(object$st)) stop("Discharge observation (Qobs) is missing in the spatio-temporal arrays (st) of transfR object. See as_transfr().")
  object$st$RnAp <- units::set_units(NA,"mm/h")
  for(i in 1:dim(object$st)[2]){
    if(verbose) progress("Computing the a priori Rn for catchment ",i,dim(object$st)[2])
    RnAp <- rapriori(Qobs = object$st$Qobs[,i], area = st_area(st_geometry(object$st)[i]), lagtime = object$lagtime[i], deltat = object$deltat)
    object$st$RnAp[,i] <- RnAp
  }
  units(object$st[["RnAp"]]) <- units(RnAp) # could not find a way to keep the units provided by rapriori
  return(object)
}

