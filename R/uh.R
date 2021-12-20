#' Unit hydrograph estimation
#'
#' Estimate the unit hydrograph from a sample of hydraulic lengths and a streamflow velocity.
#' @name uh
#' @param hl hydraulic length of class \code{stars}, \code{matrix}, \code{vector} or \code{transfR}.
#' If no unit is provided, \code{hl} is assumed to be in [m].
#' @param uc streamflow velocity. If no unit is provided, \code{uc} is assumed to be in [m/s].
#' @param deltat time step of the time series. If no unit is provided, \code{deltat} is assumed to be in [min].
#' @param verbose boolean indicating if information messages should be written to the console
#' @param ... further arguments passed to or from other methods
#' @return A data.frame with vectors of class units, or if \code{hl} is a transfR object,
#' the same transfR object incremented by the "uh" attribute.
#' @importFrom units set_units drop_units
#' @details The function estimates the unit hydrograph from geomorphometric information.
#' A travel time to the outlet is estimated by assuming an average streamflow velocity (\code{uc})
#' within the river network and  by applying \code{uc} over the sample of hydraulic lengths (\code{hl}).
#' The unit hydrograph is the probability distribution of this travel time to the outlet
#' given at each time step (\code{deltat}).
#' @examples
#' data(Oudon)
#' uh1 <- uh(hl=Oudon$hl[[1]], uc=units::set_units(0.5,"m/s"),
#' deltat=units::set_units(1,"h"))
#' plot(units::set_units(uh1$max_time,"h"), cumsum(uh1$prob), type = "b",
#' xlab = "Travel~time", ylab = "Probability~of~non-exceedance")
#'
#' object <- as_transfr(st=Oudon$obs,hl=Oudon$hl)
#' object <- velocity(object)
#' object <- uh(object)
#' plot(object,i=1,attribute=c("uh"))
#' @export
uh = function(hl, ...) UseMethod("uh")

#' @name uh
#' @export
uh.default = function(hl, uc, deltat, ...){
  warnings(paste0("No units provided for the arguments of uh(). See the help of uh() for assumptions about the units."))
  hl <- units::set_units(hl,"m")
  uc <- units::set_units(uc,"m/s")
  deltat <- units::set_units(deltat,"min")
  uh(hl=hl, uc=uc, deltat=deltat)
}

#' @name uh
#' @export
uh.units = function(hl, uc, deltat, ...){
  if(uc <= units::set_units(0,"m/s") | is.na(uc) | is.null(uc)) stop("Velocity must must be >= 0.")
  time_to_outlet <- hl/uc
  nb_classes <- ceiling(max(time_to_outlet,na.rm=TRUE) / deltat)
  lim_class_time <- (0:nb_classes) * deltat
  lim_class_dist <- uc * lim_class_time
  uh <- hist(units::drop_units(units::set_units(time_to_outlet,"min")), breaks = units::drop_units(units::set_units(lim_class_time,"min")), plot = FALSE)
  wf <- hist(units::drop_units(units::set_units(hl,"m")), breaks = units::drop_units(units::set_units(lim_class_dist,"m")), plot = FALSE)
  fdp <- data.frame("min_dist" = units::set_units(wf$breaks[1:nb_classes],"m"),
                    "max_dist" = units::set_units(wf$breaks[1:nb_classes+1],"m"),
                    "min_time" = units::set_units(uh$breaks[1:nb_classes],"min"),
                    "max_time" = units::set_units(uh$breaks[1:nb_classes+1],"min"),
                    "prob" = units::set_units(diff(uh$breaks) * uh$density,1))
  return(fdp)
}

#' @name uh
#' @export
uh.stars = function(hl, ...){
  uh(hl=hl[[1]], ...)
}

#' @name uh
#' @export
uh.transfR = function(hl, verbose=TRUE, ...){
  object <- hl
  if(!"hl"%in%names(object)) stop("Hydraulic length attribute (hl) is missing in transfR object. See as_transfr().")
  if(!"uc"%in%names(object)) stop("Streamflow velocity attribute (uc) is missing in transfR object. See velocity().")
  if(!"deltat"%in%names(object)) stop("Time step attribute (deltat) is missing in transfR object. See as_transfr().")
  fdp <- list()
  for(i in 1:length(object$hl)){
    if(verbose) progress("Computing unit hydrograph for catchment ",i,length(object$hl))
    fdp[[i]] <- uh(hl = object$hl[[i]], uc = object$uc[[i]], deltat = object$deltat)$prob
  }
  as_transfr(object = object, uh = fdp)
}



