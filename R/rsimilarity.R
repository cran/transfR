#' Prediction of hydrological similarity to drive transfer of hydrograph to ungauged catchments
#'
#' Quantify the similarity of net rainfall and model this hydrological similarity from catchment
#' similarity predictors
#' @name rsimilarity
#' @param Rn net rainfall matrix of gauged catchments (rows for time index, and columns for catchment index)
#' @param predictors a list of at least two squared matrices that could be used to predict similarity of net rainfall.
#' The matrices should be squared to describe the similarity of each pair of gauged catchments. High values should
#' indicate high probability of similarity in net rainfall
#' @param newpredictors a matrix or a data.frame of predictors, to be used as new input to the model being built,
#' and from which we want to predict net rainfall similarity. Typically, a similarity between an ungauged catchment
#' and several gauged catchments that we want to weight as donors. Each column should correspond to each \code{predictors},
#' and each row should correspond to each pair of catchments analysed.
#' @param FUN a function that will be used to quantify the similarity of two net rainfall time series. High values should
#' indicate high probability of similarity in net rainfall. Default is the inverse of root mean square error \code{invRMSE},
#' but \code{KGE} (Kling–Gupta efficiency) can also be used
#' @param power exponent applied on \code{Rn} time series when computing \code{Rn} similarity (generally between -1 and 1,
#' a lower value emphasising the similarity of low flows). Default value is 0.5
#' @param symmetrize a function to combine \code{FUN(Rn[,i],Rn[,j])} and \code{FUN(Rn[,j],Rn[,i])} into one similarity value
#' (typically if \code{FUN=KGE}, \code{KGE(Rn[,i],Rn[,j])!=KGE(Rn[,j],Rn[,i])}). Default is mean
#' @param model the method to estimate similarity of \code{Rn} from \code{descriptors}. This must be one of "lm", "glmnet"
#' @param args_glmnet list of arguments to be passed to \code{glmnet::cv.glmnet()}
#' @param verbose boolean indicating if warning messages should be written to the console
#' @param seed integer value to be used by \code{set.seed()} for reproducible results. Used only if \code{model="glmnet"}
#' @seealso \link{mixr}
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
#' res <- rsimilarity(Rn = obs$st$RnInv,
#'                    predictors = list(rghosh1=1/mdist1, rghosh2=1/mdist2),
#'                    newpredictors = data.frame(rghosh1=1/rghosh1, rghosh2=1/rghosh2),
#'                    seed=1234)
#' plot(rghosh1, res$similarity, ylab = "Predicted Rn similarity")
#' plot(rghosh2, res$similarity, ylab = "Predicted Rn similarity")
#' # rsimilarity() is automatically called by mixr() if mdist is a list
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
rsimilarity <- function(Rn, predictors, newpredictors, FUN=invRMSE, power=0.5, symmetrize=mean, model="glmnet",
		args_glmnet=list(alpha=0.5, s="lambda.min", lower.limits=0), verbose = TRUE, seed=NULL){
	
	# Checks inputs
	if(inherits(Rn, "units")) Rn <- units::drop_units(Rn)
	if(!inherits(predictors, "list")) predictors <- list(predictors)
	if(any(!sapply(predictors,FUN=is.matrix))) stop("predictors should be a list with 2 or more matrices")
	if(is.null(ncol(newpredictors))) stop("newpredictors should be a matrix or a data.frame with 2 or more columns")
	if(length(predictors)<=1) stop("predictors should be a list with 2 or more matrices")
	if(ncol(newpredictors)<=1) stop("newpredictors should be a matrix or a data.frame with 2 or more columns")
	if(any(ncol(Rn)!=sapply(predictors,FUN=ncol)) | any(ncol(Rn)!=sapply(predictors,FUN=nrow))){
		stop(paste("All predictors elements should be a square matrix of",ncol(Rn),"x",ncol(Rn),"to be consistent with Rn"))}
	if(length(predictors) != ncol(newpredictors)) stop("length(predictors) and ncol(newpredictors) should be the same")
	
	# Provide names of distances if not provided
	if(is.null(names(predictors))) names(predictors) <- paste0("dist",1:length(predictors))
	if(is.null(names(newpredictors))) names(newpredictors) <- paste0("dist",1:length(newpredictors))
	if(any(names(predictors) != colnames(newpredictors))) stop("predictors and newpredictors do not have the same names")
	
	# Compute similarity matrix index of Rn
	if(verbose) cat("\n### Computing similarity matrix of Rn ###\n")
	n <- ncol(Rn)
	mat_similarity <- matrix(NA, n, n)
	mat_similarity <- sapply(1:n, function(i) {
				sapply(1:n, function(j) FUN(Rn[, i]^power, Rn[, j]^power))
			})
	
	# Convert matrix to vector for Rn similarity
	if(verbose) cat("### Preparing data for the similarity model ###\n")
	tmat_similarity <- t(mat_similarity)
	dat_similarity <- data.frame(ij=mat_similarity[upper.tri(mat_similarity, diag = FALSE)],ji=tmat_similarity[upper.tri(tmat_similarity, diag = FALSE)])
	vec_similarity <- apply(dat_similarity, MARGIN=1, FUN=symmetrize)
	
	# Convert each matrix of predictors into vectors, to prepare data for linear model
	model_data <- matrix(nrow=length(vec_similarity), ncol=length(predictors), dimnames=list(NULL,names(predictors)))
	for(i in 1:length(predictors)){
		mat_pred <- predictors[[i]]
		tmat_pred <- t(mat_pred)
		dat_pred <- data.frame(ij=mat_pred[upper.tri(mat_pred, diag = FALSE)],ji=tmat_pred[upper.tri(tmat_pred, diag = FALSE)])
		model_data[,i] <- apply(dat_pred, MARGIN=1, FUN=symmetrize)
	}
	
	# Remove NA or Inf value
	to_remove <- apply(model_data, MARGIN=1, FUN=function(x) any(is.infinite(x) | is.na(x)))
	to_remove <- to_remove | is.na(vec_similarity)
	model_data <- model_data[!to_remove,]
	vec_similarity <- vec_similarity[!to_remove]
	
	if(verbose & nrow(model_data) < 10 * ncol(model_data)) warning(paste0("Regression is done over a limited number of similarities (",
						nrow(model_data)," observations for ",ncol(model_data)," predictors)"))
	
	if(verbose) cat("### Building the similarity model ###\n")
	if(model == "glmnet"){
		if (!is.null(seed)) set.seed(seed)
		# Linear model to predict Rn similarity from predictors using a GLM with lasso or elasticnet regularization
		valid_args <- c(names(formals(glmnet::glmnet)),names(formals(glmnet::cv.glmnet)))
		glm_similarity <- do.call(glmnet::cv.glmnet, args=c(list(x=model_data, y=vec_similarity), args_glmnet[names(args_glmnet) %in% valid_args]))
		if(is.character(args_glmnet[["s"]])){s <- glm_similarity[[args_glmnet[["s"]]]]}else{s <- args_glmnet[["s"]]}
		similarity <- as.vector(stats::predict(glm_similarity, newx=as.matrix(newpredictors), s=s))
		coef <- as.matrix(coef(glm_similarity, s=s))[-1,]
		# Compute R2 manually
		y <- vec_similarity
		y_pred <- stats::predict(glm_similarity, newx = model_data, s = s)
		y_mean <- mean(y)
		ss_total <- sum((y - y_mean)^2)
		ss_residual <- sum((y - y_pred)^2)
		r2 <- 1 - (ss_residual / ss_total)
		# Indicators of the distribution of the weights
		n50 <- which.min(cumsum(sort(similarity/sum(similarity), decreasing=TRUE)) <= 0.50)
		n75 <- which.min(cumsum(sort(similarity/sum(similarity), decreasing=TRUE)) <= 0.75)
		n90 <- which.min(cumsum(sort(similarity/sum(similarity), decreasing=TRUE)) <= 0.90)
		# Gather results
		results <- list(similarity=similarity, coef=coef, r2=r2, ncatch=as.integer(ncol(Rn)), ndist=length(vec_similarity), n50=n50, n75=n75, n90=n90)
	}
	
	if(model == "lm"){
		# Linear model to predict Rn similarity from predictors using lm
		lm_data <- cbind(Rn_similarity=vec_similarity, model_data)
		formula <- stats::as.formula(paste(colnames(lm_data)[1],"~", paste(colnames(lm_data)[-1], collapse = " + ")))
		lm_similarity <- stats::lm(formula, data=as.data.frame(lm_data))
		similarity <- as.vector(stats::predict(lm_similarity, newdata=as.data.frame(newpredictors)))
		coef <- as.matrix(coef(lm_similarity))[-1,]
		r2 <- summary(lm_similarity)$r.squared
		pvalue <- summary(lm_similarity)$coefficients[-1, "Pr(>|t|)"]
		# Indicators of the distribution of the weights
		n50 <- which.min(cumsum(sort(similarity/sum(similarity), decreasing=TRUE)) <= 0.50)
		n75 <- which.min(cumsum(sort(similarity/sum(similarity), decreasing=TRUE)) <= 0.75)
		n90 <- which.min(cumsum(sort(similarity/sum(similarity), decreasing=TRUE)) <= 0.90)
		# Gather results
		results <- list(similarity=similarity, coef=coef, r2=r2, ncatch=ncol(Rn), ndist=length(vec_similarity), n50=n50, n75=n75, n90=n90, pvalue=pvalue)
	}
	
	return(results)
}

KGE <- function(obs, sim) {
	# Remove missing values
	valid <- !is.na(obs) & !is.na(sim)
	if (!any(valid)) return(NA)
	obs <- obs[valid]
	sim <- sim[valid]
	
	# KGE components
	r <- stats::cor(obs, sim, method = c("pearson"))
	beta <- mean(sim) / mean(obs)
	gamma <- (sd(sim) / mean(sim)) / (sd(obs) / mean(obs))
	
	# Final KGE
	kge <- 1 - sqrt((r - 1)^2 + (beta - 1)^2 + (gamma - 1)^2)
	return(kge)
}

invRMSE <- function(obs, sim) {
	# Remove missing values
	valid <- !is.na(obs) & !is.na(sim)
	if (!any(valid)) return(NA)
	obs <- obs[valid]
	sim <- sim[valid]
	
	# RMSE
	rmse <- sqrt(mean((sim - obs)^2))
	return(1 / rmse)
}
