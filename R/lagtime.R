#' Lag time estimation
#'
#' Estimate the lag time of the catchment.
#' @name lagtime
#' @param hl hydraulic length of class \code{transfR} or \code{stars} or \code{matrix} or \code{vector}.
#' If no unit is provided, \code{hl} is assumed to be in [m].
#' @param uc streamflow velocity. If no unit is provided, \code{uc} is assumed to be in [m/s].
#' @param method integer describing the method to use for lag time estimation. Possible values: 1 (see details).
#' @param verbose boolean indicating if information messages should be written to the console
#' @param ... further arguments passed to or from other methods
#' @return A numeric value of class units, or if \code{hl} is a transfR object,
#' the same transfR object incremented by the "lagtime" attribute.
#' @importFrom units set_units
#' @details The function estimates the lag time of the catchment. It can be used to estimate one of the
#' inputs of the function \link{rapriori}. If \code{method} is 1, the lag time is estimated
#' from the ratio of the mean hydraulic length (\code{hl}) and the average streamflow velocity (\code{uc}).
#' @examples
#' data(Oudon)
#' icatch <- 1
#' lagtime(Oudon$hl[[icatch]],uc=units::set_units(0.5,"m/s"))
#' @export
lagtime = function(hl, ...) UseMethod("lagtime")

#' @name lagtime
#' @export
lagtime.default = function(hl, uc, ...){
  warning(paste0("No units provided for the arguments of lagtime(). See the help of lagtime() for assumptions about the units."))
  #--- Assumed units
  hl <- units::set_units(hl,"m")
  uc <- units::set_units(uc,"m/s")
  #--- lagtime
  lagtime(hl, uc, ...)
}

#' @name lagtime
#' @export
lagtime.units = function(hl, uc, method = 1, ...){
  if(method == 1) lagtime  = round(mean(hl, na.rm = TRUE)/uc, 3)
  return(units::set_units(lagtime,"h"))
}

#' @name lagtime
#' @export
lagtime.stars = function(hl, ...){
  lagtime(hl=hl[[1]], ...)
}

#' @name lagtime
#' @export
lagtime.transfR = function(hl, verbose=TRUE,  ...){
  object <- hl
  if(!"hl"%in%names(object)) stop("Hydraulic length attribute (hl) is missing in transfR argument 'hl'. See as_transfr().")
  if(!"uc"%in%names(object)) stop("Streamflow velocity attribute (uc) is missing in transfR argument 'hl'. See velocity().")
  for(i in 1:length(object$hl)){
    if(verbose) progress("Computing lag time for catchment ",i,length(object$hl))
    lti = lagtime(hl = object$hl[[i]], uc = object$uc[[i]], ...)
    if(i==1){lt = lti}else{lt = c(lt,lti)}
  }
  as_transfr(object = object, lagtime = lt)
}

