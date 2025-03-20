#' Transfer of net rainfall to ungauged catchments
#'
#' Combine the net rainfall of gauged catchments to simulate the net rainfall
#' of an ungauged catchment.
#' @name mixr
#' @param obs "transfR" object of the gauged catchments
#' @param sim "transfR" object of the ungauged catchments
#' @param mdist the geographical distance matrix between gauged and ungauged catchments, as computed by
#' the function \link{hdist}, used for donor selection
#' @param msim the similarity matrix of net rainfall time series for gauged catchments, as computed
#' by \link{rsimilarity}. It serves as the target variable in a multivariate model that predicts this
#' similarity at ungauged locations using several predictors (see the \code{similarity} argument) and
#' optimizes the donors' weighting strategy. Thus, \code{msim} is required only if \code{similarity}
#' and \code{similarity_obs} are lists containing several predictors
#' @param distance the method to use for computing distance matrix if \code{mdist} is not provided.
#' Possible values are "ghosh", "rghosh", "points", "centroids", "combined" as available
#' in the function \link{hdist}
#' @param gres resolution of spatial discretisation (number of points by km\eqn{^2}) for Ghosh
#' distance (see the function \link{hdist})
#' @param weightO weight given to the distance between outlets if \code{distance} is "combined"
#' (see the function \link{hdist})
#' @param weightC weight given to the distance between centroids if \code{distance} is "combined"
#' (see the function \link{hdist})
#' @param similarity a hydrological similarity matrix between gauged and ungauged catchments
#' used for donor weighting (\code{1/mdist^power[1]} is used if not provided). Gauged catchments
#' should be in rows, ungauged catchments in columns. Several predictors of hydrological similarity
#' can be provided in a list of matrices in order to be used by a multivariate model for predicting similarity
#' between a gauged and an ungauged catchment (see \code{predictors} of the \link{rsimilarity_model} function).
#' @param similarity_obs list of square matrices of hydrological similarity predictors between
#' gauged catchments in order to train a multivariate model for predicting similarity (\code{msim}).
#' Similar to \code{similarity} but between gauged catchments only. \code{similarity} will be used instead
#' of \code{similarity_obs} if \code{cv} is TRUE
#' @param mdist_obs the geographical distance matrix among gauged catchments, computed by
#' the function \link{hdist}, used for training a first model of hydrological similarity
#' @param FUN either a function or a character string specifying the name of a predifined function to quantify the similarity
#' of two net rainfall time series. See \link{rsimilarity} for more details
#' @param symmetrize a function to combine \code{FUN(Rn[,i],Rn[,j])} and \code{FUN(Rn[,j],Rn[,i])} into one similarity value
#' (typically if \code{FUN=KGE}, \code{KGE(Rn[,i],Rn[,j])!=KGE(Rn[,j],Rn[,i])}). Default is mean
#' @param parallel logical indicating if the computation should be parallelised
#' @param cores the number of cores to use for parallel execution if \code{parallel} is TRUE.
#' If not specified, the number of cores is set to the value of \code{parallel::detectCores()}
#' @param power exponents. The first vector value is used in the inverse distance weighting strategy
#' (as defined by the function \link{weightr}) and the second vector value is used during computation
#' of Rn similarities (as defined by the function \link{rsimilarity_model})
#' @param ndonors maximum number of catchments to be used to simulate discharge of an ungauged
#' catchment as defined by the function \link{weightr}
#' @param donors list of vectors indicating the catchment IDs from which donors are selected for each ungauged catchments.
#'  If empty, the \code{ndonors} closest catchments are used
#' @param maxdist maximum distance between a gauged and an ungauged catchment to allow the net rainfall
#' to be transferred. This threshold is applied on the \code{mdist} distance matrix. If no units is provided,
#' \code{maxdist} is assumed to be in [m].
#' @param min_ndonors minimum number of gauged catchments to start using \link{rsimilarity_model} to build a similarity model
#' of Rn time series for the weighting strategy. So this argument is only used when \code{similarity} is a list of matrices.
#' It cannot be less than 3. If less than \code{min_ndonors} donors are found inverse distance weighting
#' is applied using \code{mdist}
#' @param flexible_donor boolean indicating if the donor catchments can change during the simulation period
#' according to the availability of discharge observations. See \link{weightr} for more details
#' @param two_steps boolean indicating if a first model of hydrological similarity should be built first
#' on geographical distances. Residuals of this first model will be predicted from predictors provided by
#' \code{similarity} and \code{similarity_obs} in a second model (see \code{first_model} of \link{rsimilarity_model})
#' @param cv boolean indicating if cross validation evaluation should be done. If true, it will estimate
#' the net rainfall of every gauged catchments (\code{obs}) as if they were ungauged (leave-one-out evaluation)
#' @param save_donor boolean indicating if the net rainfall of each of the \code{ndonors} catchments
#' should be stored in the sim object for further analysis. If true, it is adding three new space-time attributes
#' in the \code{sim} object called "RnDonor", "Idonor" and "Wdonor" describing the net rainfall, the id and
#' the weight of the donor catchments respectively
#' @param save_model boolean indicating whether additional results of the local similarity model should be saved.
#' If true, then a list of data frames of observed and simulated similarities among gauged catchments
#' neighbouring each ungauged catchment will be available in the "model_training" attribute of the output
#' @param verbose boolean indicating if information messages should be written to the console
#' @param ... other arguments to be passed to \link{rsimilarity_model} if \code{similarity} is a list of matrices
#' @return The \code{sim} object incremented by the new computed attributes.
#' @details This function is a wrapper function for \link{hdist} and \link{weightr} to directly estimate
#' the net rainfall on a set of ungauged catchments (\code{sim}) from a set of gauged catchments (\code{obs}).
#' It returns the simulated net rainfall as a new space-time attribute in the \code{sim} object called "RnSim".
#' The simulated net rainfall of a given ungauged catchment \eqn{i} is a weighted average of the net rainfalls
#' of \code{ndonors} gauged catchments \eqn{j}:
#' \deqn{R_n^i =\Sigma_{j=1}^{ndonors} R_n^j \cdot \lambda_j}
#' where \eqn{\lambda_j} are defined by an inverse distance weighting function (see \link{weightr}).
#' If \code{similarity} is a list of predictors of hydrological similarity, then \eqn{\lambda_j}
#' are provided by \link{rsimilarity_model}.
#' @seealso \link{hdist}, \link{weightr}, \link{rsimilarity_model}, \link{rsimilarity}
#' @examples
#' \donttest{data(Oudon)
#' object <- as_transfr(st = Oudon$obs, hl = Oudon$hl)
#' object <- velocity(object)
#' object <- uh(object)
#' object <- lagtime(object)
#' object <- rapriori(object)
#' object <- inversion(object, parallel = TRUE, cores = 2)
#' mdist  <- hdist(x = object, y = object, method = "rghosh")
#' object <- mixr(obs = object, mdist = mdist, parallel = TRUE, cores=2,
#' cv = TRUE, flexible_donor = TRUE, save_donor = FALSE)
#' object <- convolution(object, save_donor = FALSE)
#' plot(object, i = 1, attribute = c("Qobs", "Qsim"))}
#' @importFrom stats var
#' @importFrom units drop_units set_units
#' @export
mixr <- function(obs, sim, mdist, msim = NULL, distance = "rghosh", gres = 5, weightO = 0.8, weightC = 0.2,
                 similarity, similarity_obs, mdist_obs, FUN="invKGE", symmetrize = mean,
                 parallel = FALSE, cores = NULL, power = c(1, 0.5), ndonors = 5, donors = NULL, maxdist = 50e3,
                 min_ndonors = 3, flexible_donor = TRUE, two_steps = FALSE, cv = FALSE,
                 save_donor = FALSE, save_model = FALSE, verbose = TRUE, ...){
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
  if(!is.null(donors)){
    if(!inherits(donors,"list")) stop(paste0("The class of 'donors' should be a list, with each elements being the catchment IDs of the donors for each ungauged catchment."))
    if(length(donors)!= dim(sim$st)[2]) stop(paste0("The length of 'donors' list is ",length(donors)," but the number of ungauged catchments is ",dim(sim$st)[2],"."))
    if(any(sapply(donors,FUN = function(x) !all(x %in% 1:dim(obs$st)[2])))) stop("The value of 'donors' can not be matched to a gauged catchment.")
  }

  #--- Parallelisation
  if(parallel & (missing(cores)|is.null(cores))) cores <- parallel::detectCores()
  if(!parallel) cores <- 1

  #--- Distance matrix between sim and obs
  if(!missing(mdist)){
    if(!(dim(mdist)[1]==dim(obs$st)[2] & dim(mdist)[2]==dim(sim$st)[2])) stop(paste0("The dimension of 'mdist' is (",dim(mdist)[1],",",dim(mdist)[2],") but (",dim(sim$st)[2],",",dim(obs$st)[2], ") was expected."))
  }else{
    if(verbose) cat("--- Computing distance matrix\n")
    mdist <- hdist(x=obs, y=sim, method=distance, gres=gres, weightO=weightO, weightC=weightC,
                   parallel=parallel, cores=cores, verbose=verbose)
  }
  if(missing(similarity)){similarity <- units::drop_units(1/mdist^power[1])}else{
    if(is.matrix(similarity) & any(dim(mdist)!=dim(similarity))) stop("The dimensions of 'similarity' should be similar to 'mdist'.")
    if(is.list(similarity) & any(sapply(lapply(similarity,FUN=dim),FUN=function(x) any(dim(x)!=dim(mdist))))) stop("The dimensions of each 'similarity' matrices should be similar to 'mdist'.")
  }
  if(cv){similarity_obs <- similarity}else{
    if(!missing(similarity_obs)){
      if(is.matrix(similarity_obs) & any(dim(similarity_obs)!=dim(obs$st)[2])) stop(paste0("Expected dimensions of 'similarity_obs': ",dim(obs$st)[2]," x ",dim(obs$st)[2]))
      if(is.list(similarity_obs) & any(sapply(lapply(similarity_obs,FUN=dim),FUN=function(x) any(dim(x)!=dim(obs$st)[2])))) stop(paste0("Expected dimensions of each matrices of 'similarity_obs': ",dim(obs$st)[2]," x ",dim(obs$st)[2]))
    }
  }

  #--- First core model of Rn similarity
  if(two_steps){
	if(verbose) cat("Building a first similarity model of Rn among all gauged catchments\n")
    if(cv){
      mdist_obs <- mdist
    }else{
      if(missing(mdist_obs)){
        if(verbose) cat("--- Computing distances between gauged catchments\n")
        mdist_obs <- hdist(x=obs, y=obs, method=distance, gres=gres, weightO=weightO, weightC=weightC,
                           parallel=parallel, cores=cores, verbose=verbose)

      }else{
        if(any(dim(mdist_obs)!=dim(obs$st)[2])) stop(paste0("Expected dimensions of 'mdist_obs': ",dim(obs$st)[2]," x ",dim(obs$st)[2]))
      }
    }
    if(inherits(mdist_obs,"units")) mdist_obs <- units::set_units(mdist_obs,"m") |>  units::drop_units()
    if(missing(msim) | is.null(msim)){
      if(verbose) cat("--- Computing similarity matrix of Rn between gauged catchments\n")
      Rn <- units::drop_units(obs$st$RnInv)
      msim <- rsimilarity(Rn=Rn^power[2], FUN=FUN, parallel= parallel, cores=cores)
    }
    if(verbose) cat("--- Prepare input data for global Rn similarity model\n")
    # Convert matrix to vector
    model_list <- list(mdist_obs, msim)
    nrow <- (nrow(mdist_obs)^2-nrow(mdist_obs))/2
    model_data <- data.frame(distance = vector(length=nrow, mode="numeric"), similarity = vector(length=nrow, mode="numeric"))
    for(i in 1:length(model_list)){
      m <- model_list[[i]]
      tm <- t(m)
      dat <- data.frame(ij=m[upper.tri(m, diag = FALSE)],ji=tm[upper.tri(tm, diag = FALSE)])
      model_data[,i] <- apply(dat, MARGIN=1, FUN=symmetrize)
    }
    # Remove NA or Inf value
    to_remove <- apply(model_data, MARGIN=1, FUN=function(x) any(is.infinite(x) | is.na(x)))
    model_data <- model_data[!to_remove,]
    # Summarize by classes
    nb_classes <- 20 #may become an argument of the function
    breaks <- seq(0, units::drop_units(maxdist), length.out = nb_classes + 1)
    model_data[,"dist_class"] <- as.character(cut(model_data[,"distance"], breaks = breaks, include.lowest = TRUE))
    model_data_summary <- stats::aggregate(cbind(distance, similarity) ~ dist_class, data = model_data,
              FUN = function(x) stats::median(x, na.rm = TRUE))
    # Fit a linear model with optimal power transformation
    if(verbose) cat("--- Fit a global Rn similarity model\n")
    powers <- seq(-5, 5, by = 0.01)
    r2 <- vector(length=length(powers), mode="numeric")
    for(ipower in seq_along(powers)) r2[ipower] <- summary(lm(similarity ~ distance_power, data = data.frame(model_data_summary,distance_power=model_data_summary$distance^powers[ipower])))$r.squared
    bestpower <- powers[which.max(r2)]
    first_model <- eval(bquote(lm(similarity ~ I(distance^(.(bestpower))), data = model_data_summary)))
    if(verbose) cat(paste0("R\u00B2: ",round(max(r2),4),"\nBest power transformation: ",bestpower,"\n\n"))
    model_data_summary$fitted.value <- first_model$fitted.values
    first_model$data <- model_data_summary
    sim$first_similarity_model <- first_model
  }else{
    first_model <- NULL
  }

  #--- Transferring Rn from obs to sim
  if(is.list(similarity)){
    if(is.null(names(similarity))) names(similarity) <- paste0("predictor",1:length(similarity))
	  cnames <- c(names(similarity), "intercept", "r2", "ncatch", "ndist", "n50", "n75", "n90")
	  similarity_models <- as.data.frame(matrix(ncol = length(cnames), nrow = dim(sim$st)[2]))
	  colnames(similarity_models) <- cnames
	  if(two_steps){
	    similarity <- c(list(first_predictor = mdist), similarity)
	    similarity_obs <- c(list(first_predictor = mdist_obs), similarity_obs)
	  }
	  if(save_model) model_training <- list()
  }
  sim$st$RnSim <- NA
  for(i in 1:dim(sim$st)[2]){
    if(verbose) progress("Combining simulated Rn of neighbouring catchments for catchment ",i,dim(sim$st)[2])
    if(verbose & is.list(similarity)) cat("\n--- Selection of potential donors based on distance\n")
    if(cv){distances <- mdist[-i,i]}else{distances <- mdist[,i]}
    if(is.null(donors)){idonors <- which(distances<=maxdist)}else{idonors <- donors[[i]][donors[[i]]%in%which(distances<=maxdist)]}
    if(length(idonors)==0){
      warning(paste0("No donor found for catchment ",i," within ",units::set_units(maxdist,"km")," km."))
      sim$st$RnSim[,i] <- rep(NA,dim(sim$st)[1])
	}else{
      if(cv){
        RnInv <- obs$st$RnInv[,-i][,idonors]
        if(!(missing(msim) | is.null(msim))){msim_local <- msim[-i,-i][idonors,idonors]}else{msim_local <- NULL}
      }else{
        RnInv <- obs$st$RnInv[,idonors]
        if(!(missing(msim) | is.null(msim))){msim_local <- msim[idonors,idonors]}else{msim_local <- NULL}
      }
	  if(length(idonors)==1){
		  warning(paste0("Only one donor found for catchment ",i," within ", units::set_units(maxdist,"km"),
						  " km.\nNo weighting strategy applied."))
		  sim$st$RnSim[,i] <- RnInv
		  weights <- matrix(NA, nrow=length(RnInv), ncol=1)
		  weights[!is.na(RnInv),1] <- 1
	  }else{
		  if(is.list(similarity)){
			  if(cv){
				  predictors <- lapply(similarity_obs, FUN = function(x) x[-i,-i][idonors,idonors])
				  newpredictors<- sapply(similarity, FUN = function(x) x[-i,i])[idonors,]
			  }else{
				  predictors <- lapply(similarity_obs, FUN = function(x) x[idonors,idonors])
				  newpredictors<- sapply(similarity, FUN = function(x) x[idonors,i])
			  }
		    if(two_steps){
		      no_variance <- sum(sapply(predictors[-1], FUN=function(x) stats::var(x[!is.infinite(x)], na.rm=TRUE)))==0
		    }else{
		      no_variance <- sum(sapply(predictors, FUN=function(x) stats::var(x[!is.infinite(x)], na.rm=TRUE)))==0
		    }
			  if(is.na(no_variance)) no_variance <- TRUE
			  nsimilarities <- max(apply(RnInv, MARGIN=1, FUN=function(x) sum(!is.na(x)))) # How many similarity metrics can be computed
			  if(length(idonors)>=min_ndonors & length(idonors)>=3 & !no_variance & nsimilarities>=min_ndonors){ #Need to have at least 3 catchments to have several Rn similarities
				  similarity_model <- rsimilarity_model(Rn=RnInv,
				      msim=msim_local,
						  FUN=FUN,
						  power=power[2],
						  predictors=predictors,
						  newpredictors=newpredictors,
						  first_model=first_model,
						  parallel=parallel,
						  cores=cores,
						  ...)
				  similarities <- similarity_model$similarity
				  # Save description of the similarity model
				  similarity_models[i,] <- data.frame(t(similarity_model$coef), intercept=similarity_model$intercept,
						  r2=similarity_model$r2, ncatch=similarity_model$ncatch, ndist=similarity_model$ndist,
						  n50=similarity_model$n50, n75=similarity_model$n75, n90=similarity_model$n90)
				  if(save_model) model_training[[i]] <- similarity_model$model_training
			  }else{
				  if(no_variance){
					  warning(paste0("Predictors have zero variance among the ",length(idonors),
									  " donors found for catchment ",i," within ", units::set_units(maxdist,"km"),
									  " km.\nInverse Distance Weighting using 'mdist' will be used by default."))
				  }else{
					  warning(paste0("Only ",length(idonors)," donors found for catchment ",i," within ", units::set_units(maxdist,"km"),
									  " km. Considering the gauging period, only ", nsimilarities," Rn similarities can be calculated.",
									  "\nInverse Distance Weighting using 'mdist' will be used by default."))
				  }
				  similarities <- units::drop_units(1/distances[idonors]^power[1])
			  }
		  }else{
			  if(cv){similarities <- similarity[-i,i][idonors]}else{similarities <- similarity[idonors,i]}
		  }
		  if(verbose & is.list(similarity)) cat("--- Weighting Rn of each donor catchments\n")
		  weights <- weightr(Rn=RnInv, distances=distances[idonors], similarities=similarities, ndonors=ndonors, power=power[1], flexible_donor=flexible_donor)
		  sim$st$RnSim[,i] <- apply(RnInv*weights,MARGIN = 1,FUN = function(x){if(all(is.na(x))){NA}else{sum(x, na.rm = T)}})
	  }
      if(save_donor){
		if(verbose) cat("--- Saving information about donors\n")
        if(cv){gdonors <- seq_len(dim(sim$st)[2])[-i][idonors]}else{gdonors <- seq_len(dim(sim$st)[2])[idonors]} #to come back to general spatial index
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
          sim$st[[paste0("Idonor",idonor)]][,i]  <- apply(weights,MARGIN = 1, FUN = function(x){if(all(is.na(x))){NA}else{gdonors[which(rank(1-x,ties.method = "first")==idonor)]}})
          sim$st[[paste0("Wdonor",idonor)]][,i] <- apply(weights,MARGIN = 1, FUN = function(x){if(all(is.na(x))){NA}else{x[which(rank(1-x,ties.method = "first")==idonor)]}})
          sim$st[[paste0("RnDonor", idonor)]][, i] <- sapply(seq_len(dim(sim$st)[1]), function(t) sim$st$RnInv[t, sim$st[[paste0("Idonor", idonor)]][t, i]])
          units(sim$st[[paste0("RnDonor",idonor)]]) <- units(RnInv) # could not find a way to keep the units provided by weighting
        }
      }
    }
  }
  if(is.list(similarity)){
    sim[["similarity_models"]] <- similarity_models # save similarity models summary
    if(save_model) sim[["model_training"]] <- model_training # save similarity models training data set and results
  }
  units(sim$st[["RnSim"]]) <- units(obs$st$RnInv) # could not find a way to keep the units provided by weighting
  return(sim)
}
