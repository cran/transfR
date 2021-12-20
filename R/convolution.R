#' Convolution of net rainfall with unit hydrograph
#'
#' Simulate the discharge by a convolution between the unit hydrograph and the net rainfall.
#' @name convolution
#' @param Rn net rainfall vector or an object of class \code{transfR}
#' @param uh unit hydrograph vector
#' @param Rcol name of the space-time attribute for the discharge simulation in the \code{transfR} object
#' @param Qcol name of the space-time attribute for the net rainfall in the \code{transfR} object
#' @param continuous boolean indicating if, when one time step might be influenced by past or future rainfall
#' (according to the time span of the unit hydrograph), no simulated value is provided
#' @param save_donor boolean indicating if additional discharge simulations should be computed using the
#' net rainfall of each individual donor catchment instead of just the weighted average net rainfall. This
#' requires that \code{save_donor} was TRUE when using \link{mixr}
#' @param verbose boolean indicating if information messages should be written to the console
#' @param ... further arguments passed to or from other methods
#' @return An object of the same class of \code{Rn}. If \code{Rn} is a transfR object,
#' the same transfR object incremented by the new computed attributes.
#' @examples
#' data(Oudon)
#' icatch <- 1
#' uc <- velocity(hl = Oudon$hl[[icatch]])
#' uh <- uh(hl = Oudon$hl[[icatch]], uc = uc, deltat = units::set_units(1,"h"))$prob
#' Rn <- units::set_units(c(1,5,2),"mm/h")
#' Qsim <- convolution(Rn = Rn, uh = uh)
#' @import sf stars
#' @importFrom units set_units drop_units
#' @export
convolution = function(Rn, ...) UseMethod("convolution")

#' @name convolution
#' @export
convolution.default <- function(Rn, uh, continuous=FALSE, ...){
  Q <- rep(0,length(Rn)+length(uh))
  for(t in which(Rn>0)) Q[t:(t+length(uh)-1)] <- Rn[t]*uh + Q[t:(t+length(uh)-1)]
  for(t in which(is.na(Rn))) Q[t:(t+length(uh)-1)] <- NA
  if(continuous){
    Q[1:(length(uh)-1)] <- NA
    return(Q[1:length(Rn)])
  }else{return(Q)}
}

#' @name convolution
#' @export
convolution.units <- function(Rn, uh, ...){
  Q <- convolution.default(Rn=units::drop_units(Rn), uh=units::drop_units(uh), ...)
  units(Q) <- units(Rn)
  return(Q)
}


#' @name convolution
#' @export
convolution.transfR = function(Rn, Rcol="RnSim", Qcol="Qsim", save_donor=FALSE, verbose=TRUE, ...){
  object <- Rn
  if(!"st"%in%names(object)) stop("Spatio-temporal arrays (st) of class stars not found in stars object. See as_transfr().")
  if(!Rcol%in%names(object$st)) stop(paste0("Attribute '",Rcol,"' not found in stars object.\nUse mixr() function first to simulate net rainfall from neighbouring catchments or define an existing attribute using 'Rcol' argument."))
  object$st[[Qcol]] <- matrix(NA,nrow=dim(object$st)[1],ncol=dim(object$st)[2])
  dim(object$st[[Qcol]]) <- dim(object$st)
  for(i in 1:dim(object$st)[2]){
    if(verbose) progress("Computing simulated discharge for catchment ",i,dim(object$st)[2])
    Qsim <- convolution(Rn=object$st[[Rcol]][,i],uh=object$uh[[i]],continuous=TRUE)
    Qsim <- Qsim*st_area(st_geometry(object$st)[i])
    Qsim <- units::set_units(Qsim,"m^3/s")
    object$st[[Qcol]][,i] <- Qsim
  }
  if(save_donor){
    iRcol <- grep("RnDonor",names(object$st))
    if(length(iRcol)==0) stop("save_donor is TRUE but no 'RnDonor' attributes found. Use also save_donor=TRUE when using mixr() function to keep memory of donor's Rn.")
    if(verbose) cat("\n")
    for(i in iRcol){
      if(verbose) progress("Additional discharge simulation using the single donor ",which(i==iRcol),length(iRcol),"")
      object <- convolution(Rn = object, Rcol = names(object$st)[i], Qcol = gsub("RnDonor","QsimDonor",names(object$st)[i]),
                            continuous = TRUE, save_donor = FALSE, verbose = FALSE, ...)
    }
    #Re-order columns just for readability of numerous attributes
    n <- names(object$st)
    icol <- c(grep("Qsim",n),grep("RnSim",n),grep("RnDonor",n),grep("Idonor",n),grep("Wdonor",n))
    object$st <- object$st[c(icol,(1:length(n))[!(1:length(n))%in%icol])]
  }
  units(object$st[[Qcol]]) <- units(Qsim) # could not find a way to keep the units
  return(object)
}


