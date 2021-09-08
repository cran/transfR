#' Transfer of net rainfall and discharge simulation
#'
#' Combine net rainfall of gauged catchments to simulate the discharge
#' of an ungauged catchment.
#' @name mixr
#' @param obs "transfR" object of the gauged catchments
#' @param sim "transfR" object of the ungauged catchments
#' @param mdist the distance matrix between gauged and ungauged catchments as computed by
#' the function \link{hdist}
#' @param distance the method to use for computing distance matrix if \code{mdist} is not provided.
#' Possible values are "ghosh", "ghosh2", "rghosh", "points", "centroids", "combined" as available
#' in the function \link{hdist}
#' @param gres resolution of spatial discretisation (number of points by km²) for Ghosh
#' distance (see the function \link{hdist})
#' @param weightO weight given to the distance between outlets if \code{distance} is "combined" (see the function \link{hdist})
#' @param weightC weight given to the distance between centroids if \code{distance} is "combined" (see the function \link{hdist})
#' @param parallel logical indicating if the computation should be parallelised
#' @param cores the number of cores to use for parallel execution if \code{parallel} is TRUE.
#' If not specified, the number of cores is set to the value of \code{parallel::detectCores()}
#' @param power exponent applied in the inverse distance weighting strategy as defined by the function
#' \link{weightr}
#' @param ndonors maximum number of catchments to be used to simulate discharge of an ungauged
#' catchment as defined by the function \link{weightr}
#' @param donors vector of catchments id from which donors are selected. If empty, the \code{ndonors} closest
#' catchments are used
#' @param maxdist maximum distance between a gauged and an ungauged catchment to allow the net rainfall
#' to be transfered. This threshold is applied on the \code{mdist} distance matrix. If no units is provided,
#' \code{maxdist} is assumed to be in [m].
#' @param flexible_donor boolean indicating if the donor catchments can change during the simulation period
#' according to the availability of discharge observations. See \link{weightr} for more details
#' @param cv boolean indicating if cross validation evaluation should be done. If true, it will estimate
#' the net rainfall of every gauged catchments (\code{obs}) as if they were ungauged (leave-one-out evaluation)
#' @param save_donor boolean indicating if the net rainfall of each of the \code{ndonors} catchments
#' should be stored in the sim object for further analysis. If true, it is adding three new space-time attributes
#' in the \code{sim} object called "RnDonor", "Idonor" and "Wdonor" describing the net rainfall, the id and
#' the weight of the donor catchments respectively
#' @param verbose boolean indicating if information messages should be written to the console
#' @return The \code{sim} object incremented by the new computed attributes.
#' @details This function is a wrapper function for \link{hdist} and \link{weightr} to directly estimate
#' the net rainfall on a set of ungauged catchments (\code{sim}) from a set of gauged catchments (\code{obs}).
#' It returns the simulated net rainfall as a new space-time attribute in the \code{sim} object called "RnSim".
#' The simulated net rainfall of a given ungauged catchment \eqn{i} is a weighted average of the net rainfalls
#' of \code{ndonors} gauged catchments \eqn{j}:
#' \deqn{R_n^i =\Sigma_{j=1}^{ndonors} R_n^j \cdot \lambda_j}
#' where \eqn{\lambda_j} are defined by a inverse distance weighting function (see \link{weightr}).
#' @seealso \link{hdist}, \link{weightr}
#' @examples
#' \donttest{data(Oudon)
#' object <- as_transfr(st=Oudon$obs,hl=Oudon$hl)
#' object <- velocity(object)
#' object <- uh(object)
#' object <- lagtime(object)
#' object <- rapriori(object)
#' object <- inversion(object, parallel = TRUE, cores=2)
#' mdist  <- hdist(x=object, y=object, method="rghosh")
#' object <- mixr(obs = object, mdist=mdist, parallel = TRUE, cores=2,
#' cv=TRUE, flexible_donor=TRUE, save_donor = FALSE)
#' object <- convolution(object, save_donor = FALSE)
#' plot(object,i=1,attribute=c("Qobs","Qsim"))}
#' @export
mixr <- function(obs, sim, mdist, distance = "rghosh", gres = 5, weightO = 0.8, weightC = 0.2,
                    parallel = FALSE, cores = NULL, power = 1, ndonors = 5, donors = NULL, maxdist = 50e3,
                    flexible_donor = TRUE, cv = FALSE, save_donor = FALSE, verbose = TRUE){
  if(cv){
    if(missing(sim)){
      sim <- obs
    }else{
      stop("A cross validation can not be done on ungauged catchments ('sim' argument). Use only 'obs' argument or deactivate cross validation.")
    }
  }

  maxdist <- units::set_units(maxdist,"m")

  #--- Checking class
  if(!inherits(obs,"transfR")) stop("The class of 'obs' must be 'transfR'. See function 'as_transfr()'.")
  if(!inherits(sim,"transfR")) stop("The class of 'sim' must be 'transfR'. See function 'as_transfr()'.")

  #--- Checking needed inputs
  if(!"st"%in%names(sim)) stop("Spatio-temporal arrays (st) of class stars is missing in 'sim' argument. See as_transfr().")
  if(!"st"%in%names(obs)) stop("Spatio-temporal arrays (st) of class stars is missing in 'obs' argument. See as_transfr().")
  if(!"RnInv"%in%names(obs$st)) stop("Net rainfall simulated by inversion (RnInv) is missing in the spatio-temporal arrays (st) of 'obs' argument. See inversion().")

  #--- Checking donors
  if(!missing(donors)){
    if(length(donors)!= dim(sim$st)[2]) stop("Length of 'donors' list must be the same as the number of ungauged catchments.")
    if(any(sapply(donors,FUN = function(x) !all(x %in% 1:dim(obs$st)[2])))) stop("The value of 'donors' can not be matched to a gauged catchment.")
  }

  #--- Distance matrix beetween sim and obs
  if(!missing(mdist)){
    if(!(dim(mdist)[1]==dim(obs$st)[2] & dim(mdist)[2]==dim(sim$st)[2])) stop(paste0("The dimension of 'mdist' is (",dim(mdist)[1],",",dim(mdist)[2],") but (",dim(sim$st)[2],",",dim(obs$st)[2], ") was expected."))
  }else{
    if(verbose) cat("### Computing distance matrix ###\n")
    mdist <- hdist(x=obs, y=sim, method=distance, gres=gres, weightO=weightO, weightC=weightC,
                   parallel=parallel, cores=cores)
  }

  #--- Transfering Rn from obs to sim
  sim$st$RnSim <- NA
  for(i in 1:dim(sim$st)[2]){
    if(verbose) progress("Combining simulated Rn of neighbouring catchments for catchment ",i,dim(sim$st)[2])
    if(cv){
      RnInv <- obs$st$RnInv[,-i]
      distances <- mdist[-i,i]
    }else{
      RnInv <- obs$st$RnInv
      distances <- mdist[,i]}
    if(missing(donors)){idonors <- which(distances<=maxdist)}else{idonors <- donors[[i]][donors[[i]]%in%which(distances<=maxdist)]}
    if(length(idonors)==0){
      warning(paste0("No donor found for catchment ",i," below the distance of ",units::set_units(maxdist,"km")," km."))
      sim$st$RnSim[,i] <- rep(NA,dim(sim$st)[1])
    }else{
      weights <- weightr(Rn=RnInv, distances=distances, ndonors=ndonors, donors = idonors, power=power, flexible_donor=flexible_donor)
      sim$st$RnSim[,i] <- apply(RnInv*weights,MARGIN = 1,FUN = function(x){if(all(is.na(x))){NA}else{sum(x, na.rm = T)}})
      if(save_donor){
        for(idonor in 1:min(c(ndonors,length(idonors)))){
          #Create new columns
          if(!paste0("RnDonor",idonor)%in%names(sim$st)){
            sim$st[[paste0("RnDonor",idonor)]] <- matrix(NA,nrow=dim(sim$st)[1],ncol=dim(sim$st)[2])
            dim(sim$st[[paste0("RnDonor",idonor)]]) <- dim(sim$st)}
          if(!paste0("Idonor",idonor)%in%names(sim$st)){
            sim$st[[paste0("Idonor",idonor)]]  <- matrix(NA,nrow=dim(sim$st)[1],ncol=dim(sim$st)[2])
            dim(sim$st[[paste0("Idonor",idonor)]]) <- dim(sim$st)}
          if(!paste0("Wdonor",idonor)%in%names(sim$st)){
            sim$st[[paste0("Wdonor",idonor)]]  <- matrix(NA,nrow=dim(sim$st)[1],ncol=dim(sim$st)[2])
            dim(sim$st[[paste0("Wdonor",idonor)]]) <- dim(sim$st)}
          sim$st[[paste0("Idonor",idonor)]][,i]  <- apply(weights,MARGIN = 1, FUN = function(x){if(all(is.na(x))){NA}else{which(rank(1-x,ties.method = "first")==idonor)}})
          sim$st[[paste0("Wdonor",idonor)]][,i] <- apply(weights,MARGIN = 1, FUN = function(x){if(all(is.na(x))){NA}else{x[which(rank(1-x,ties.method = "first")==idonor)]}})
          for(t in 1:dim(sim$st)[1]) sim$st[[paste0("RnDonor",idonor)]][t,i] <- RnInv[t,sim$st[[paste0("Idonor",idonor)]][t,i]]
          # units(sim$st[[which(names(sim$st)==paste0("RnDonor",idonor))]]) <- units(RnInv) # could not find a way to keep the units provided by weighting
          units(sim$st[[paste0("RnDonor",idonor)]]) <- units(RnInv) # could not find a way to keep the units provided by weighting
        }
      }
    }
  }
  units(sim$st[["RnSim"]]) <- units(RnInv) # could not find a way to keep the units provided by weighting
  return(sim)
}

