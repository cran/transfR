#' Create transfR object
#'
#' Create a transfR object or add new attributes to a transfR object.
#' @name as_transfr
#' @param object object of class \code{transfR}
#' @param st spatio-temporal arrays of class \code{stars}. Observed discharge must be described by the column name 'Qobs'.
#' Time should be the first dimension, space the second dimension. If no unit is provided, Qobs is assumed to be in [m3/s]
#' and RnInv is assumed to be in [mm/h] (or [mm/d] at daily time step).
#' @param uc vector of the streamflow velocities of the catchments. If no unit is provided, \code{uc} is assumed to be in [m/s].
#' @param lagtime vector of the lag times of the catchments. If no unit is provided, \code{lagtime} is assumed to be in [h].
#' @param surface vector of the surfaces of the catchments. If no unit is provided, \code{surface} is assumed to be in [km2].
#' @param delineation spatial layer of the boundary of the catchments of class \code{sfc_POLYGON}.
#' @param outlet spatial layer of the outlets of the catchments of class \code{sfc_POINT}.
#' @param centroid spatial layer of the centroids of the catchments of class \code{sfc_POINT}.
#' @param uh list of the unit hydrographs of the catchments.
#' @param hl hydraulic length of class \code{stars}, \code{matrix} or \code{vector}. If no unit is provided, \code{hl}
#' is assumed to be in [m]. See details below.
#' @return An object of class transfR.
#' @details This function creates an object of class \code{transfR} or increment an existing \code{transfR} object
#' with new attributes. It can be used to gather and organize most of the inputs and outputs of the other functions
#' like streamflow velocities, unit hydrograph, a priori on net rainfall, inversions and simulations of every catchments.
#'
#' This function can be used to organise the two user inputs required for a conventional use of the package, namely \code{st}
#' and \code{hl}. The hydraulic lengths are defined as the flow path length from each pixel to the outlet within
#' the river network \insertCite{Cudennec2004,Aouissi2013}{transfR}. Catchment delineations and hydraulic lengths
#' need to be prepared beforehand by the user. This package does not provide functions to create them.
#' However, several GIS software offer possibilities to extract them from a digital elevation model
#' such as GRASS toolkits \insertCite{Jasiewicz2011}{transfR},
#' Whitebox GAT (\insertCite{@see @Lindsay2016;textual}{transfR} or \href{https://github.com/jblindsay/whitebox-tools}{WhiteboxTools}),
#' TauDEM (D. Tarboton, Utah State University)
#' or online services (\insertCite{@see @Squividant2015;textual}{transfR} for catchment delineation in the Brittany French region).
#' @examples
#' data(Oudon)
#' object <- as_transfr(st=Oudon$obs,hl=Oudon$hl)
#' @import sf stars
#' @importFrom units set_units
#' @references
#' \insertRef{Aouissi2013}{transfR}
#'
#' \insertRef{Cudennec2004}{transfR}
#'
#' \insertRef{Jasiewicz2011}{transfR}
#'
#' \insertRef{Lindsay2016}{transfR}
#'
#' \insertRef{Squividant2015}{transfR}
#' @export
as_transfr <- function(object, st, uc, lagtime, surface, delineation, outlet, centroid, uh, hl){

  #--Create transfR object
  if(missing(object)){
    object <- list()
    class(object)=c("transfR","list")
  }

  #--Adding attributes to transfR object
  if(!missing(st)){

    #--Checking class
    if(!inherits(st, "stars")) stop("The class of 'st' must be stars.")
    if(!all(st_geometry_type(st) == c("POLYGON"))) stop("st geometry type has to be a POLYGON.")
    #--Checking attributes
    if(length(names(st))>0&!any(names(st)%in%c("Qobs","RnInv"))) stop("No column name is recognized in 'st'. Fill in 'Qobs' or 'RnInv' time series.")
    # if("Qobs"%in%names(st)) st[[which(names(st)=="Qobs")]] <- units::set_units(st[[which(names(st)=="Qobs")]],"m^3/s")

    #--Deducing time step
    deltat <- unique(difftime(st_get_dimension_values(st,1)[-1],st_get_dimension_values(st,1)[-dim(st)[1]],units="mins"))
    if(length(deltat)==1){
      deltat <- units::set_units(deltat,"min")
      object$deltat <- deltat
    }else{stop("Time step must be steady.")}

    #--Set up units
    if("Qobs"%in%names(st)) st[["Qobs"]] <- units::set_units(st[["Qobs"]],"m^3/s")
    if(object$deltat==units::set_units(1,"day")){rn_units <- units::as_units("mm/d")}else{rn_units <- units::as_units("mm/h")}
    if("RnInv"%in%names(st)) st[["RnInv"]] <- units::set_units(st[["RnInv"]],rn_units,mode="standart")
    object$st <- st

  }

  #--Checking new attributes
  if(!missing(uc)){
    if(!(inherits(uc, "numeric")|inherits(uc, "units"))) stop("The class of 'uc' must be numeric or units.")
    if(length(uc)!=dim(object$st)[2]) stop(paste0("Length of 'uc' is ",length(uc)," but there are ",dim(object$st)[2]," spatial features."))
    uc <- units::set_units(uc,"m/s")
    object$uc <- uc
  }
  if(!missing(lagtime)){
    if(!(inherits(lagtime, "numeric")|inherits(lagtime, "units"))) stop("The class of 'lagtime' must be numeric or units.")
    if(length(lagtime)!=dim(object$st)[2]) stop(paste0("Length of 'lagtime' is ",length(lagtime)," but there are ",dim(object$st)[2]," spatial features."))
    lagtime <- units::set_units(lagtime,"h")
    object$lagtime <- lagtime
  }
  if(!missing(surface)){
    if(!(inherits(surface, "numeric")|inherits(surface, "units"))) stop("The class of 'surface' must be numeric or units.")
    if(length(surface)!=dim(object$st)[2]) stop(paste0("Length of 'surface' is ",length(surface)," but there are ",dim(object$st)[2]," spatial features."))
    surface <- units::set_units(surface,"km2")
    object$surface <- surface
  }
  if(!missing(delineation)){
    if(!inherits(delineation, "sfc_POLYGON")) stop("The class of 'delineation' must be sfc_POLYGON")
    if(length(delineation)!=dim(object$st)[2]) stop(paste0("Length of 'delineation' is ",length(delineation)," but expecting ",dim(object$st)[2]," spatial features."))
    object$delineation <- delineation
  }
  if(!missing(outlet)){
    if(!inherits(outlet, "sfc_POINT")) stop("The class of 'outlet' must be sfc_POINT")
    if(length(outlet)!=dim(object$st)[2]) stop(paste0("Length of outlet is ",length(outlet)," but expecting ",dim(object$st)[2]," spatial features."))
    object$outlet <- outlet
  }
  if(!missing(centroid)){
    if(!inherits(centroid, "sfc_POINT")) stop("The class of 'centroid' must be sfc_POINT")
    if(length(centroid)!=dim(object$st)[2]) stop(paste0("Length of 'centroid' is ",length(centroid)," but expecting ",dim(object$st)[2]," spatial features."))
    object$centroid <- centroid
  }
  if(!missing(uh)){
    if(!inherits(uh, "list")) stop("The class of 'uh' must be a list")
    if(length(uh)!=dim(object$st)[2]) stop(paste0("Length of 'uh' is ",length(uh)," but there are ",dim(object$st)[2]," spatial features."))
    classes <- c("numeric","units")
    if(any(!unlist(lapply(uh,class))%in%classes)) stop(paste("Variables contained in the list of 'uh' can only be",paste(classes,collapse=" or ")))
    object$uh <- uh
  }
  if(!missing(hl)){
    if(!inherits(hl, "list")) stop("The class of 'hl' must be a list")
    if(length(hl)!=dim(object$st)[2]) stop(paste0("Length of 'hl' is ",length(hl)," but there are ",dim(object$st)[2]," spatial features."))
    classes <- c("stars","matrix","units")
    if(any(!unlist(lapply(hl,class))%in%classes)) stop(paste("Variables contained in the list of 'hl' can only be",paste(classes,collapse=" or ")))
    object$hl <- lapply(hl,FUN=function(x){x[[1]] <- units::set_units(x[[1]],"m");return(x)})
  }
  return(object)
}

