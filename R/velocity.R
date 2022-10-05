#' Streamflow velocity estimation
#'
#' Estimate streamflow velocity in average over the catchment.
#' @name velocity
#' @param hl hydraulic length of class \code{stars}, \code{matrix}, \code{vector} or \code{transfR}.
#' If no unit is provided, \code{hl} is assumed to be in [m].
#' @param lagtime lag time of the catchment. If no unit is provided, \code{lagtime} is assumed to be in [h].
#' @param method character string describing the method to estimate the velocity. One of "loire2016" (default), "brittany2013"
#' or "lagtime" (see details).
#' @param ... further arguments passed to or from other methods
#' @return A numeric value of class units, or if \code{hl} is a transfR object,
#' the same transfR object incremented by the "uc" attribute.
#' @details Estimate the average streamflow velocity of the catchment from three different
#' approaches. Method "lagtime" estimates the velocity from the ratio between the mean
#' hydraulic length and the lag time of the catchment. Method "loire2016" estimates the
#' velocity from a regression based on hydraulic length only:
#' \deqn{a \cdot hl^b}
#' where \eqn{a=4.38e-4} and \eqn{b=0.69} have been calibrated over the Loire river basin \insertCite{deLavenne2016}{transfR}.
#' Method "brittany2013" used a similar regression calibrated for the French Brittany region where \eqn{a=8.59e-4}
#' and \eqn{b=0.61} \insertCite{deLavenne2013}{transfR}.
#' @examples
#' data(Oudon)
#' velocity(Oudon$hl[[1]], method = "loire2016")
#'
#' object <- as_transfr(st = Oudon$obs, hl = Oudon$hl)
#' object <- velocity(object)
#' object$uc
#' @importFrom units set_units drop_units
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{deLavenne2013}{transfR}
#'
#' \insertRef{deLavenne2016}{transfR}
#' @export
velocity <- function(hl, ...) UseMethod("velocity")

#' @name velocity
#' @export
velocity.default = function(hl, lagtime, method = "loire2016", ...){
  warning(paste0("No units provided for the arguments of velocity(). See the help of velocity() for assumptions about the units."))
  hl <- units::set_units(hl,"m")
  if(method == "lagtime"){
    lagtime <- units::set_units(lagtime,"h")
    velocity.units(hl=hl, lagtime=lagtime, method=method)
  }
  if(method %in% c("loire2016","brittany2013")) velocity.units(hl=hl, method=method)
}

#' @name velocity
#' @export
velocity.units <- function(hl, lagtime = NULL, method = "loire2016", ...){
  method <- match.arg(method, c("lagtime","loire2016","brittany2013"), several.ok = FALSE)
  if(method == "lagtime"){
    if(is.null(lagtime)) stop("A lag time needs to be provided for this method.")
    uc  <- round(mean(hl, na.rm = TRUE)/lagtime, 3)}
  if(method == "loire2016"){
    hl <- units::set_units(hl,"m")
    uc <- round(4.38e-4 * mean(units::drop_units(hl), na.rm = TRUE)^0.69, 3)}
  if(method == "brittany2013"){
    hl <- units::set_units(hl,"m")
    uc <- round(8.59e-4 * mean(units::drop_units(hl), na.rm = TRUE)^0.61, 3)}
  uc <- units::set_units(uc,"m/s")
  return(uc)
}

#' @name velocity
#' @export
velocity.stars = function(hl, ...){
  hl = hl[[1]]
  velocity(hl=hl, ...)
}

#' @name velocity
#' @export
velocity.transfR = function(hl, ...){
  object <- hl
  if(!"hl"%in%names(object)) stop("Hydraulic length attribute (hl) is missing in transfR object. See as_transfr().")
  uc = sapply(object$hl,FUN=function(x){velocity(hl = x, ...)})
  as_transfr(object = object, uc = uc)
}

