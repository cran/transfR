#' Transfer of observed discharge from gauged catchments to ungauged catchments
#'
#' Wrap up all the modelling steps into one function for a quick implementation of this R package.
#' @name quick_transfr
#' @param obs "transfR" object of the gauged catchments
#' @param sim "transfR" object of the ungauged catchments
#' @param velocity character string describing the method to estimate the streamflow velocity.
#' See \link{velocity} for the available options (\code{method} argument)
#' @param distance character string describing the method to compute the distance between catchments.
#' See \link{hdist} for the available options (\code{method} argument)
#' @param gres resolution of spatial discretisation (number of points by km²) for Ghosh
#' distance. See \link{hdist} for more details
#' @param weightO weight given to the distance between outlets if distance method is "combined".
#' See \link{hdist} for more details
#' @param weightC weight given to the distance between centroids if distance method is "combined".
#' See \link{hdist} for more details
#' @param power exponent applied in the inverse distance weighting strategy. See \link{weightr} for more details
#' @param ndonors maximum number of catchments to be used to simulate discharge of an ungauged
#' catchment. See \link{weightr} for more details
#' @param maxdist maximum distance between a gauged and an ungauged catchment to allow the net rainfall
#' to be transfered. This threshold is applied on the \code{mdist} distance matrix. If no units is provided,
#' \code{maxdist} is assumed to be in [m]. See \link{mixr} for more details
#' @param flexible_donor boolean indicating if the donor catchments can change during the simulation period
#' according to the availability of discharge observations. See \link{weightr} for more details
#' @param cv boolean indicating if cross validation evaluation should be done. If true, it will estimate
#' the net rainfall of every gauged catchments (\code{obs}) as if they were ungauged (leave-one-out evaluation)
#' @param save_donor boolean indicating if the net rainfall of each of the \code{ndonors} catchments
#' should be stored in the sim object for further analysis. If true, it is adding three new space-time attributes
#' in the \code{sim} object called "RnDonor", "Idonor" and "Wdonor" describing the net rainfall, the id and
#' the weight of the donor catchments respectively. See \link{mixr} for more details
#' @param warmup length of the warmup period. If no unit is provided, \code{warmup} is assumed to be in [days].
#' See \link{inversion} for more details
#' @param cooldown length of the period removed at the end of the simulation. If no unit is provided,
#' \code{cooldown} is assumed to be in [days]. See \link{inversion} for more details
#' @param dosplit boolean, if true the inversion is performed by
#' subperiods of length defined by \code{split}. See \link{inversion} for more details
#' @param split length the subperiods if dosplit is true. If no unit is provided, \code{split} is assumed to be
#' in [days]. See \link{inversion} for more details
#' @param parallel logical indicating if the computation should be parallelised
#' @param cores the number of cores to use for parallel execution if \code{parallel} is TRUE.
#' If not specified, the number of cores is set to the value of \code{parallel::detectCores()}
#' @param verbose boolean indicating if information messages should be written to the console
#' @return The \code{sim} object incremented by the new computed attributes
#' @details The function applies sequentially the following functions: \link{velocity}, \link{uh},
#' \link{lagtime}, \link{rapriori}, \link{inversion}, \link{hdist}, \link{mixr} and \link{convolution}.
#' Please refer to the help of each of these functions and to \link{transfR} for a general description of the
#' modelling approach.
#' @seealso \link{velocity}, \link{uh}, \link{lagtime}, \link{rapriori}, \link{inversion},
#' \link{hdist}, \link{mixr}, \link{convolution}
#' @examples
#' \donttest{data(Oudon)
#' obs <- as_transfr(st = Oudon$obs[,,1:3], hl = Oudon$hl[1:3]) #gauged catchments
#' sim <- as_transfr(st = Oudon$obs[,,4:6], hl = Oudon$hl[4:6]) #catchments considered as ungauged
#' sim <- quick_transfr(obs, sim)}
#' @export
quick_transfr <- function(obs, sim, velocity = "loire2016", distance = "rghosh", gres = 5, weightO = 0.8, weightC = 0.2,
                 power = 1, ndonors = 5, maxdist = 50e3, flexible_donor = TRUE, cv = FALSE, save_donor = FALSE,
                 warmup = 10, cooldown = 8, dosplit = TRUE, split = 30,
                 parallel = FALSE, cores = NULL, verbose = TRUE){

  # Cross validation
  if(cv){
    if(missing(sim)){
      sim <- obs
    }else{
      stop("A cross validation can not be done on ungauged catchments ('sim' argument). Use only 'obs' argument or deactivate cross validation.")
    }
  }

  if(verbose) cat("### 1. BUILDING GEOMORPHOLOGY-BASED TRANSFER FUNCTION \n")

  # Streamflow velocity
  obs <- velocity(obs, method = velocity)
  if(!cv) sim <- velocity(sim, method = velocity)

  # A geomorphology based unit hydrograph
  obs <- uh(obs, verbose = verbose)
  if(!cv) sim <- uh(sim, verbose = verbose)

  if(verbose) cat("\n\n### 2. NET RAINFALL SIMULATION OF GAUGED CATCHMENTS BY DECONVOLUTION \n")

  # Net rainfall estimated a priori
  obs <- lagtime(obs, verbose = verbose)
  obs <- rapriori(obs, verbose = verbose)
  if(verbose) cat("\n")

  # Net rainfall estimated by inversion
  obs <- inversion(obs, warmup = warmup, cooldown = cooldown, dosplit = dosplit, split = split,
                   parallel = parallel, cores = cores, verbose = verbose)

  if(verbose) cat("\n\n### 3. NET RAINFALL SIMULATION OF UNGAUGED BASINS BY SPATIAL INTERPOLATION \n")

  # Compute distance matrix
  mdist <- hdist(x = obs, y = sim, method = distance, gres = gres, weightO = weightO, weightC = weightC,
                 parallel = parallel, cores = cores, verbose = verbose)

  # Estimate net rainfall at ungauged locations
  if(!cv){
    sim <- mixr(obs = obs, sim = sim, mdist = mdist, power = power, ndonors = ndonors, maxdist = maxdist,
                flexible_donor = flexible_donor, cv = cv, save_donor = save_donor,
                parallel = parallel, cores = cores, verbose = verbose)
  }else{
    obs <- mixr(obs = obs, mdist = mdist, power = power, ndonors = ndonors, maxdist = maxdist,
                 flexible_donor = flexible_donor, cv = cv, save_donor = save_donor,
                parallel = parallel, cores = cores, verbose = verbose)}

  if(verbose) cat("\n\n### 4. DISCHARGE SIMULATION OF UNGAUGED BASINS BY CONVOLUTION \n")

  # Simulate discharge at ungauged locations
  if(!cv){
    sim <- convolution(sim, verbose = verbose)
    return(sim)
  }else{
    obs <- convolution(obs, verbose = verbose)
    return(obs)
  }

}
