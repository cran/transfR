#' Prediction of hydrological similarity to drive transfer of hydrograph to ungauged catchments
#'
#' Quantify the similarity of net rainfall and model this hydrological similarity from catchment
#' similarity predictors
#' @name rsimilarity_model
#' @param Rn net rainfall matrix of gauged catchments (rows for time index, and columns for catchment index)
#' @param msim similarity matrix of net rainfall time series (as produced by \link{rsimilarity})
#' @param predictors a list of at least two squared matrices that could be used to predict similarity of net rainfall.
#' The matrices should be squared to describe the similarity of each pair of gauged catchments. High values should
#' indicate high probability of similarity in net rainfall
#' @param newpredictors a matrix or a data.frame of predictors, to be used as new input to the model being built,
#' and from which we want to predict net rainfall similarity. Typically, a similarity between an ungauged catchment
#' and several gauged catchments that we want to weight as donors. Each column should correspond to each \code{predictors},
#' and each row should correspond to each pair of catchments analysed.
#' @param FUN either a function or a character string specifying the name of a predifined function to quantify the similarity
#' of two net rainfall time series. See \link{rsimilarity} for more details
#' @param power exponent applied on \code{Rn} time series when computing \code{Rn} similarity (generally between -1 and 1,
#' a lower value emphasising the similarity of low flows). Default value is 0.5
#' @param symmetrize a function to combine \code{FUN(Rn[,i],Rn[,j])} and \code{FUN(Rn[,j],Rn[,i])} into one similarity value
#' (typically if \code{FUN=KGE}, \code{KGE(Rn[,i],Rn[,j])!=KGE(Rn[,j],Rn[,i])}). Default is mean
#' @param model the method to estimate similarity of \code{Rn} from \code{descriptors}. This must be one of "lm", "glmnet"
#' @param args_glmnet list of arguments to be passed to \code{glmnet::cv.glmnet()}
#' @param standardisation boolean indicating whether the predictors should be standardised
#' (in particular to make their coefficient comparable)
#' @param first_model a first model to use for predicting similarity (generally based on distance, see \code{two_steps} of
#' \link{mixr}). The actual model will be built on residuals of the first one
#' @param parallel logical indicating whether the computation of Rn similarities should be parallelised
#' @param cores the number of cores to use for parallel execution if \code{parallel} is TRUE.
#' If not specified, the number of cores is set to the value of \code{parallel::detectCores()}
#' @param verbose boolean indicating if warning messages should be written to the console
#' @param seed integer value to be used by \code{set.seed()} for reproducible results. Used only if \code{model="glmnet"}
#' @seealso \link{rsimilarity}, \link{mixr}
#' @examples
#' \donttest{data(Oudon)
#' obs <- as_transfr(st = Oudon$obs, hl = Oudon$hl)
#' obs <- velocity(obs, method = "loire2016")
#' obs <- uh(obs)
#' obs <- lagtime(obs)
#' obs <- rapriori(obs)
#' obs <- inversion(obs, parallel = TRUE, cores=2)
#' mdist1 <- hdist(x = obs, y = obs, method = "rghosh", parallel = c(FALSE,TRUE), cores=2)
#' mdist2 <- mdist1^2
#' rghosh1 <- seq(1000, 5000, by=100)
#' rghosh2 <- rghosh1^2
#' res <- rsimilarity_model(Rn = obs$st$RnInv,
#'                    predictors = list(rghosh1=1/mdist1, rghosh2=1/mdist2),
#'                    newpredictors = data.frame(rghosh1=1/rghosh1, rghosh2=1/rghosh2),
#'                    seed=1234)
#' plot(rghosh1, res$similarity, ylab = "Predicted Rn similarity")
#' plot(rghosh2, res$similarity, ylab = "Predicted Rn similarity")
#' # rsimilarity_model() is automatically called by mixr() if mdist is a list
#' obs <- mixr(obs = obs, mdist = mdist1,
#'             similarity = list(rghosh1=1/mdist1, rghosh2=1/mdist2),
#'             parallel = TRUE, cores=2, cv = TRUE, save_donor = TRUE)
#' obs$similarity_models
#' obs <- convolution(obs)
#' plot(obs, i = 1, attribute = c("Qobs", "Qsim"))}
#' @importFrom units drop_units
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom stats lm predict complete.cases cor as.formula
#' @export
rsimilarity_model <- function(Rn, msim=NULL, predictors, newpredictors, FUN = "invKGE", power = 0.5, symmetrize = mean,
                              model = "glmnet", args_glmnet=list(s="lambda.min", lower.limits=0), standardisation = FALSE,
                              first_model = NULL, parallel = FALSE, cores = NULL, verbose = TRUE, seed = NULL){

	# Checks inputs: predictors and newpredictors
	if(!inherits(predictors, "list")) predictors <- list(predictors)
	if(any(!sapply(predictors,FUN=is.matrix))) stop("predictors should be a list with 2 or more matrices")
	if(is.null(ncol(newpredictors))) stop("newpredictors should be a matrix or a data.frame with 2 or more columns")
	if(length(predictors)<=1) stop("predictors should be a list with 2 or more matrices")
	if(ncol(newpredictors)<=1) stop("newpredictors should be a matrix or a data.frame with 2 or more columns")
	if(length(predictors) != ncol(newpredictors)) stop("length(predictors) and ncol(newpredictors) should be the same")

	# Provide names of distances if not provided
	if(is.null(names(predictors))) names(predictors) <- paste0("predictor",1:length(predictors))
	if(is.null(names(newpredictors))) names(newpredictors) <- paste0("predictor",1:length(newpredictors))
	if(any(names(predictors) != colnames(newpredictors))) stop("predictors and newpredictors do not have the same names")

	# Compute similarity matrix index of Rn
	if(is.null(msim)){
		# Checks inputs: Rn
		if(missing(Rn)){
			stop("Need to provide either Rn or msim")
		}else{
			if(inherits(Rn, "units")) Rn <- units::drop_units(Rn)
			if(any(ncol(Rn)!=sapply(predictors,FUN=ncol)) | any(ncol(Rn)!=sapply(predictors,FUN=nrow))){
				stop(paste("All predictors elements should be a square matrix of",ncol(Rn),"x",ncol(Rn),"to be consistent with Rn"))}
		}
		# Parallelisation
		if(parallel & (missing(cores)|is.null(cores))) cores <- parallel::detectCores()
		if(!parallel) cores <- 1
		if(verbose) cat("--- Computing similarity matrix of Rn\n")
		msim <- rsimilarity(Rn=Rn^power, FUN=FUN, parallel= parallel, cores=cores)
	}else{
		# Checks inputs: msim
		if(!inherits(msim, "matrix")) stop("msim should be a matrix")
		if(nrow(msim) != ncol(msim)) stop("msim should be a square matrix")
		if(any(nrow(msim) != nrow(predictors))) stop("msim should be a square matrix")
		if(any(ncol(msim)!=sapply(predictors,FUN=ncol)) | any(ncol(msim)!=sapply(predictors,FUN=nrow))){
			stop(paste("msim should be a square matrix with dimensions similar to predictors"))}
	}

	# Convert matrix to vector for Rn similarity
	if(verbose) cat("--- Preparing the data for the local similarity model\n")
	tmsim <- t(msim)
	dat_similarity <- data.frame(ij=msim[upper.tri(msim, diag = FALSE)],ji=tmsim[upper.tri(tmsim, diag = FALSE)])
	similarity_gauged <- apply(dat_similarity, MARGIN=1, FUN=symmetrize)

	# Convert each matrix of predictors into vectors, to prepare data for linear model
	model_data <- matrix(nrow=length(similarity_gauged), ncol=length(predictors), dimnames=list(NULL,names(predictors)))
	for(i in 1:length(predictors)){
		mat_pred <- predictors[[i]]
		tmat_pred <- t(mat_pred)
		dat_pred <- data.frame(ij=mat_pred[upper.tri(mat_pred, diag = FALSE)],ji=tmat_pred[upper.tri(tmat_pred, diag = FALSE)])
		model_data[,i] <- apply(dat_pred, MARGIN=1, FUN=symmetrize)
	}

	# Remove NA or Inf value
	to_remove <- apply(model_data, MARGIN=1, FUN=function(x) any(is.infinite(x) | is.na(x)))
	to_remove <- to_remove | is.infinite(similarity_gauged) | is.na(similarity_gauged)
	model_data <- model_data[!to_remove,]
	similarity_gauged <- similarity_gauged[!to_remove]

	# Compute residuals of the first model, if provided
	if(!is.null(first_model)){
	  if(verbose) cat("--- Computing the residuals of the first global similarity model\n")
	  if(is.null(utils::getS3method("predict", class(first_model), optional = TRUE))) stop("No predict() method is defined for the class of 'first_model'.")
	  first_similarity_gauged <- stats::predict(first_model, newdata = data.frame(distance=model_data[,1]))
	  first_residuals_gauged <- first_similarity_gauged-similarity_gauged
	  target_variable <- asinh(first_residuals_gauged) # Now the new goal is to predict residuals of the first model
	  model_data <- model_data[,-1] # Remove variable that was carried out for the first model only
	  first_similarity_ungauged <- stats::predict(first_model, newdata = data.frame(distance=newpredictors[,1]))
	  newpredictors <- newpredictors[,-1] # Remove variable that was carried out for the first model only
	}else{
	  target_variable <- similarity_gauged
	}

	if(verbose & nrow(model_data) < 10 * ncol(model_data)) warning(paste0("Regression is done over a limited number of similarities (",
						nrow(model_data)," observations for ",ncol(model_data)," predictors)"))

	if(standardisation){
	  for(i in 1:ncol(model_data)){
	    imean <- mean(model_data[,i], na.rm=TRUE)
	    isd <- sd(model_data[,i], na.rm=TRUE)
	    model_data[,i] <- (model_data[,i]-imean)/isd
	    newpredictors[,i] <- (newpredictors[,i]-imean)/isd
	  }
	}

	if(verbose) cat("--- Building the similarity model")
	if(model == "glmnet"){
		if (!is.null(seed)) set.seed(seed)
		# Linear model to predict Rn similarity from predictors using a GLM with lasso or elasticnet regularization
		valid_args <- c(names(formals(glmnet::glmnet)),names(formals(glmnet::cv.glmnet)))
		glm_similarity <- do.call(glmnet::cv.glmnet, args=c(list(x=model_data, y=target_variable), args_glmnet[names(args_glmnet) %in% valid_args]))
		if(is.character(args_glmnet[["s"]])){s <- glm_similarity[[args_glmnet[["s"]]]]}else{s <- args_glmnet[["s"]]}
		prediction <- as.vector(stats::predict(glm_similarity, newx=as.matrix(newpredictors), s=s))
		coef <- as.matrix(coef(glm_similarity, s=s))
		intercept <- coef[1,]
		coef <- coef[-1,]
		# Infinite similarity of predictors leads to infinite predicted similarity (instead of NA)
		inf_prediction <- apply(newpredictors, MARGIN=1, FUN = function(x) any(is.infinite(x)))
		if(any(inf_prediction)) prediction[inf_prediction] <- Inf
		# Compute R2 manually
		prediction_gauged <- stats::predict(glm_similarity, newx = model_data, s = s)
		r2 <- 1 - sum((target_variable - prediction_gauged)^2) / sum((target_variable - mean(target_variable))^2)
		# Indicators of the distribution of the weights
		n50 <- which.min(cumsum(sort(prediction/sum(prediction), decreasing=TRUE)) <= 0.50)
		n75 <- which.min(cumsum(sort(prediction/sum(prediction), decreasing=TRUE)) <= 0.75)
		n90 <- which.min(cumsum(sort(prediction/sum(prediction), decreasing=TRUE)) <= 0.90)
		# Gather results
		if(!is.null(first_model)){ # get back to original value
			prediction <- sinh(prediction)
			prediction_gauged <- sinh(prediction_gauged)}
		results <- list(prediction=prediction, coef=coef, intercept=intercept, r2=r2, ncatch=as.integer(ncol(msim)), ndist=length(target_variable), n50=n50, n75=n75, n90=n90)
	}

	if(model == "lm"){
		# Linear model to predict Rn similarity from predictors using lm
		lm_data <- cbind(target_variable, model_data)
		formula <- stats::as.formula(paste(colnames(lm_data)[1],"~", paste(colnames(lm_data)[-1], collapse = " + ")))
		lm_similarity <- stats::lm(formula, data=as.data.frame(lm_data))
		prediction <- as.vector(stats::predict(lm_similarity, newdata=as.data.frame(newpredictors)))
		prediction_gauged <- lm_similarity$fitted.values
		coef <- as.matrix(coef(lm_similarity))
		intercept <- coef[1,]
		coef <- coef[-1,]
		r2 <- summary(lm_similarity)$r.squared
		pvalue <- summary(lm_similarity)$coefficients[-1, "Pr(>|t|)"]
		# Infinite similarity of predictors leads to infinite predicted similarity (instead of NA)
		inf_prediction <- apply(newpredictors, MARGIN=1, FUN = function(x) any(is.infinite(x)))
		if(any(inf_prediction)) prediction[inf_prediction] <- Inf
		# Indicators of the distribution of the weights
		n50 <- which.min(cumsum(sort(prediction/sum(prediction), decreasing=TRUE)) <= 0.50)
		n75 <- which.min(cumsum(sort(prediction/sum(prediction), decreasing=TRUE)) <= 0.75)
		n90 <- which.min(cumsum(sort(prediction/sum(prediction), decreasing=TRUE)) <= 0.90)
		# Gather results
		if(!is.null(first_model)){ # get back to original value
			prediction <- sinh(prediction)
			prediction_gauged <- sinh(prediction_gauged)}
		results <- list(prediction=prediction, coef=coef, intercept=intercept, r2=r2, ncatch=as.integer(ncol(msim)), ndist=length(target_variable), n50=n50, n75=n75, n90=n90, pvalue=pvalue)
	}
	if(verbose) cat(paste0(" (R\u00B2 = ",round(r2,2),")\n"))

	if(!is.null(first_model)){
	  results[["similarity"]] <- first_similarity_ungauged - prediction
	  results[["first_similarity"]] <- first_similarity_ungauged
	  results[["model_training"]] <- data.frame(similarity=similarity_gauged, similarity_pred=c(first_similarity_gauged-prediction_gauged), first_similarity_pred=first_similarity_gauged)
	}else{
	  results[["similarity"]] <- prediction
	  results[["model_training"]] <- data.frame(similarity=similarity_gauged, similarity_pred=prediction_gauged)
	}

	return(results)
}
