#' Weights of donor catchments
#'
#' Estimate the weighting applied at each time step and to each gauged catchment (donors) for the calculation
#' of the average net rainfall of an ungauged catchment
#' @name weightr
#' @param Rn net rainfall matrix of donor catchments (rows for time index, and columns for donors index)
#' @param distances vector of the distances to each donor catchment (see \link{hdist})
#' @param ndonors maximum number of donor catchments to use
#' @param donors vector of catchments id from which donors are selected. If empty, the \code{ndonors} closest
#' catchments are used
#' @param power exponent applied in the inverse distance weighting function (see details)
#' @param flexible_donor boolean indicating if the donor catchments can change during the simulation period
#' according to the availability of discharge observations (see details)
#' @return A matrix with the same dimensions as \code{Rn}.
#' @importFrom units set_units as_units drop_units
#' @details This function returns a matrix of weights for each time steps (rows) and each gauged catchments
#' (columns) for the calculation of the average net rainfall of an ungauged catchment (see \link{mixr}).
#' The weights \eqn{\lambda} are estimated by an inverse distance weighting function \insertCite{deLavenne2016}{transfR}:
#' \deqn{\lambda_i=\frac{1}{d_i^p}}
#' \deqn{\Sigma_{i=1}^{ndonors}\lambda_i=1}
#' where \eqn{d} is the \code{distances} argument and \eqn{p} is the \code{power} argument. The weights are rescaled
#' so the sum is equal to 1.
#'
#' Two strategies to handle missing data in the \code{Rn} matrix are possible.
#' If \code{flexible_donor} is TRUE, donors catchments are redefined at each time steps, and are chosen among
#' the ones that are effectively gauged at this given time step. This aims to keep a constant number of donor
#' catchments throughout the simulation period.
#' If \code{flexible_donor} is FALSE, the donor catchments are chosen once within all the gauged catchments,
#' regardless of missing data and remain the same throughout the entire simulation period. This stability of
#' donor catchments might however leads to a reduced number of donors (below \code{ndonors}) during periods
#' of missing data.
#' @seealso \link{hdist}, \link{mixr}
#' @references
#' \insertRef{deLavenne2016}{transfR}
#' @export
weightr <- function(Rn, distances, ndonors = 5, donors, power = 1, flexible_donor = TRUE){
  if(inherits(distances, "units")) distances <- units::drop_units(distances)
  weights <- matrix(0,nrow=nrow(Rn),ncol=ncol(Rn))
  if(flexible_donor){ #Choose the best ndonor donors within available data at each time step.
    for(t in 1:nrow(Rn)){
      gaps    <- is.na(Rn[t,])
      if(all(gaps)){weights[t,]=NA}else{
        if(missing(donors)){tdonors  <- which(!gaps)}else{tdonors  <- donors[donors%in%which(!gaps)]}
        tdonors  <- tdonors[order(distances[tdonors])[1:min(c(length(tdonors),ndonors))]]
        if(any(distances[tdonors]==0)){
          weights[t,tdonors] <- 0
          weights[t,tdonors[which(distances[tdonors]==0)]] <- 1/sum(distances[tdonors]==0)
        }else{
          weights[t,tdonors]=1/distances[tdonors]^power
          weights[t,tdonors]=weights[t,tdonors]/sum(weights[t,tdonors])
        }
      }
    }
  }else{ #Keep always the same donors all along the time series. The number of donors might be reduced according to the available data.
    norn    <- apply(Rn,MARGIN=2,FUN=function(x){all(is.na(x))})
    stable_donors  <- which(!norn)
    if(missing(donors)){stable_donors  <- which(!norn)}else{stable_donors  <- donors[donors%in%which(!norn)]}
    stable_donors  <- stable_donors[order(distances[stable_donors])[1:min(c(length(stable_donors),ndonors))]]
    for(t in 1:nrow(Rn)){
      gaps    <- is.na(Rn[t,stable_donors])
      tdonors  <- stable_donors[!gaps]
      if(all(gaps)){weights[t,]=NA}else{
        if(any(distances[tdonors]==0)){
          weights[t,tdonors] <- 0
          weights[t,tdonors[which(distances[tdonors]==0)]] <- 1/sum(distances[tdonors]==0)
        }else{
          weights[t,tdonors]=1/distances[tdonors]^power
          weights[t,tdonors]=weights[t,tdonors]/sum(weights[t,tdonors])
        }
      }
    }
  }
  return(weights)
}
