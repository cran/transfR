#' Similarity of net rainfall time series
#'
#' Calculate the similarity of net rainfall time series
#' @name rsimilarity
#' @param Rn net rainfall matrix of gauged catchments (rows for time index, and columns for catchment index)
#' @param FUN either a function or a character string specifying the name of a predifined function to quantify the similarity of
#' two net rainfall time series. Higher values should indicate a higher probability of similarity. Predefined functions include:
#' "KGE" (Kling–Gupta efficiency), "invRMSE" (inverse of root mean square error), "invKGE" (inverse of (KGE-1))
#' and "RMSE" (root mean square error, should not be directly used as a similarity metric). The default is "invKGE"
#' @param parallel logical indicating whether the computation should be parallelised
#' @param cores the number of cores to use for parallel execution if \code{parallel} is TRUE.
#' If not specified, the number of cores is set to the value of \code{parallel::detectCores()}
#' @return A square matrix of the similarity metric between each pair of catchments
#' @seealso \link{rsimilarity_model}
#' @examples
#' \donttest{data(Oudon)
#' obs <- as_transfr(st = Oudon$obs, hl = Oudon$hl)
#' obs <- velocity(obs, method = "loire2016")
#' obs <- uh(obs)
#' obs <- lagtime(obs)
#' obs <- rapriori(obs)
#' obs <- inversion(obs, parallel = TRUE, cores=2)
#' msim <- rsimilarity(Rn = obs$st$RnInv, FUN="KGE", parallel = TRUE, cores=2)}
#' @importFrom units drop_units
#' @export
rsimilarity <- function(Rn, FUN="invKGE", parallel=FALSE, cores=NULL){

	# Checks
	if(inherits(Rn, "units")) Rn <- units::drop_units(Rn)
	if(!(inherits(Rn, "matrix") | inherits(Rn, "data.frame"))) stop("Rn should be a matrix or a data.frame")
	if(ncol(Rn) < 2) stop("Rn should have at least two columns")

	# Parallelisation
	if(parallel & (missing(cores)|is.null(cores))) cores <- parallel::detectCores()
	if(!parallel)	cores <- 1

	# Initialise msim
	n <- ncol(Rn)
	msim <- matrix(NA, n, n)

	# Compute similarity matrix
	if(inherits(FUN, "function")){
		msim <- sapply(1:n, function(i) {sapply(1:n, function(j) FUN(Rn[, j], Rn[, i]))})
	}else{
		valid_FUN <- c("KGE", "RMSE", "invKGE", "invRMSE")
		FUN <- match.arg(FUN, choices = valid_FUN)
		crit <- match(FUN, valid_FUN)
		msim <- .Call("c_similarity", Rn=Rn, crit=crit, nthreads=cores)
	}
	return(msim)
}
