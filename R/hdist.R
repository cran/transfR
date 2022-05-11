#' Geographical distance between catchments
#'
#' Calculate distances between two sets of catchments using their spatial support.
#' @name hdist
#' @param x sf, stars or transfR object of the first catchments
#' @param y sf, stars or transfR object of the second catchments
#' @param method the method to use for computing distance. This must be one of "ghosh",
#' "rghosh", "points", "centroids", "combined"
#' @param gres resolution of spatial discretisation (number of points by km²) for Ghosh
#' distance
#' @param weightO weight given to the distance between outlets if method is "combined"
#' @param weightC weight given to the distance between centroids if method is "combined"
#' @param ditself logical value indicating if the distance to itself should be computed.
#' It will add one row and one column in the distance matrix. Only used if method is "ghosh"
#' @param maxsample maximum size of sampling points for each catchments during spatial discretisation
#' @param proj logical indicating if spatial layer are using a projection. If TRUE, euclidean
#' distance is used. If FALSE, the great-circle distance is used
#' @param parallel logical indicating if the computation should be parallelised
#' @param cores the number of cores to use for parallel execution if \code{parallel} is TRUE.
#' If not specified, the number of cores is set to the value of \code{parallel::detectCores()}
#' @param verbose boolean indicating if information messages should be written to the console
#' @param ... further arguments passed to or from other methods
#' @return A matrix of class units with the catchments of \code{x} organised in rows
#' and the catchments of \code{y} organised in columns.
#' @details The \code{method} "ghosh" refers to a simplification of the distance defined
#' by \insertCite{Ghosh1951;textual}{transfR} as proposed by \insertCite{Gottschalk1993,Gottschalk2011;textual}{transfR}.
#' The rescaled Ghosh distance (\code{method} "rghosh") is calculted following \insertCite{deLavenne2016;textual}{transfR}.
#' @references
#' \insertRef{Ghosh1951}{transfR}
#'
#' \insertRef{Gottschalk1993}{transfR}
#'
#' \insertRef{Gottschalk2011}{transfR}
#'
#' \insertRef{deLavenne2016}{transfR}
#' @examples
#' data(Oudon)
#' catchments <- st_geometry(Oudon$obs)
#' hdist(x = catchments[1:2], y = catchments[3:5], gres = 5, method = "rghosh")
#' @import sf stars doParallel foreach
#' @importFrom units set_units drop_units
#' @useDynLib transfR, .registration = TRUE
#' @export
hdist <- function(x, y, ...) UseMethod("hdist")

#' @name hdist
#' @export
hdist.sfc <- function(x, y, method="rghosh", gres=5, ditself=FALSE, maxsample=2.5e4, proj=NULL, parallel=FALSE, cores=NULL, verbose=TRUE, ...){
  method <- match.arg(method,c("ghosh","rghosh","points","centroids","combined"),several.ok = FALSE)
  if(missing(proj)){
    if(!identical(st_crs(x),st_crs(y))) stop("x and y do not use the same coordinate reference system.")
    if(is.na(st_crs(x)$input)){stop("Don't know how to compute distances because no coordinate reference system was found. Define a coordinate reference system (using st_crs()) or define 'proj' argument.")}
    if(st_crs(x)$input=="EPSG:4326"){proj=FALSE}else{proj=TRUE}
  }
  if(method=="ghosh"){
    if(!(all(st_geometry_type(x) == c("POLYGON")) & all(st_geometry_type(y) == c("POLYGON")))) stop("Geometry type has to be a POLYGON if method 'rghosh' is used.")
    if(verbose) cat("Sampling the catchments at a resolution of about",gres,"pts/km2\n")
    size <- round(units::set_units(gres,"1/km^2")*units::set_units(st_area(st_union(c(x,y))),"km^2"))
    xydisc <- st_sample(st_union(c(x,y)),size=units::drop_units(size),type="regular")
    xdisc <- list()
    ydisc <- list()
    for(i in 1:length(x)) xdisc[[i]] <- xydisc[x[i]]
    for(i in 1:length(y)) ydisc[[i]] <- xydisc[y[i]]
    if(identical(x,y)){
      for(i in 1:length(xdisc)) if(length(xdisc[[i]])>maxsample){
        xdisc[[i]] <- ydisc[[i]] <- xdisc[[i]][sort(sample.int(length(xdisc[[i]]), size = maxsample, replace = FALSE))]
        if(verbose) cat("Random resampling of catchment",i,"with a maximum size sample of",maxsample,"points.\n")}
    }else{
      for(i in 1:length(xdisc)){
        if(length(xdisc[[i]])>maxsample){
          xdisc[[i]] <- xdisc[[i]][sort(sample.int(length(xdisc[[i]]), size = maxsample, replace = FALSE))]
          if(verbose) cat("Random resampling of catchment x[",i,"] with a maximum size sample of ",maxsample," points.\n",sep="")}
      }
      for(i in 1:length(ydisc)){
        if(length(ydisc[[i]])>maxsample){
          ydisc[[i]] <- ydisc[[i]][sort(sample.int(length(ydisc[[i]]), size = maxsample, replace = FALSE))]
          if(verbose) cat("Random resampling of catchment y[",i,"] with a maximum size sample of ",maxsample," points.\n",sep="")}
      }
    }
    if(verbose) cat("Computing Ghosh distance between catchments\n")
    if(parallel){
      if(missing(cores)|is.null(cores)) cores=parallel::detectCores()
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl=cl)
      on.exit(parallel::stopCluster(cl))
      if(!identical(x,y)){
        res <- foreach::"%dopar%"(foreach::foreach(i=1:length(x),.packages="sf"),
                                  sapply(ydisc,FUN=function(x){.Fortran("gdist",coord1=st_coordinates(x),coord2=st_coordinates(xdisc[[i]]),n1=length(x),n2=length(xdisc[[i]]),proj=proj,rescale=FALSE,diag=FALSE,mdist=0)$mdist}))
        gdist <- matrix(unlist(res), nrow = length(x), ncol = length(y), byrow = TRUE)
      }else{
        res <- foreach::"%dopar%"(foreach::foreach(i=1:length(x),.packages="sf"),{
          tmp <- sapply(ydisc[i:length(x)],FUN=function(x){.Fortran("gdist",coord1=st_coordinates(x),coord2=st_coordinates(xdisc[[i]]),n1=length(x),n2=length(xdisc[[i]]),proj=proj,rescale=FALSE,diag=FALSE,mdist=0)$mdist})
          c(rep(NA,i-1),tmp)
        })
        gdist <- matrix(unlist(res), nrow = length(x), ncol = length(y), byrow = TRUE)
        for(i in 1:length(x)) gdist[,i] <- gdist[i,]
      }
    }else{
      gdist <- matrix(NA,nrow=length(x),ncol=length(y))
      if(!identical(x,y)){
        for(i in 1:length(x)) gdist[i,] <- sapply(ydisc,FUN=function(x){.Fortran("gdist",coord1=st_coordinates(x),coord2=st_coordinates(xdisc[[i]]),n1=length(x),n2=length(xdisc[[i]]),proj=proj,rescale=FALSE,diag=FALSE,mdist=0)$mdist})
      }else{
        for(i in 1:length(x)) gdist[i,i:length(x)] <- gdist[i:length(x),i] <- sapply(ydisc[i:length(x)],FUN=function(x){.Fortran("gdist",coord1=st_coordinates(x),coord2=st_coordinates(xdisc[[i]]),n1=length(x),n2=length(xdisc[[i]]),proj=proj,rescale=FALSE,diag=FALSE,mdist=0)$mdist})
      }
    }
    if(ditself){
      if(verbose) cat(paste0("Computing Ghosh distance within catchments\n"))
      if(!identical(x,y)){
        gdist <- cbind(gdist,sapply(xdisc,FUN=function(x){.Fortran("gdist",coord1=st_coordinates(x),coord2=st_coordinates(x),n1=length(x),n2=length(x),proj=proj,rescale=FALSE,diag=TRUE,mdist=0)$mdist}))
        gdist <- rbind(gdist,c(sapply(ydisc,FUN=function(x){.Fortran("gdist",coord1=st_coordinates(x),coord2=st_coordinates(x),n1=length(x),n2=length(x),proj=proj,rescale=FALSE,diag=TRUE,mdist=0)$mdist}),NA))
      }else{
        gdist <- cbind(gdist,diag(gdist))
        gdist <- rbind(gdist,c(diag(gdist),NA))
      }
    }
    units(gdist) <- units(st_distance(xdisc[[1]][1],ydisc[[1]][1])) #units are lost by some operations
    return(gdist)
  }
  if(method=="rghosh"){
    gdist <- hdist(x=x, y=y, method="ghosh", gres=gres, ditself=TRUE, maxsample=maxsample, proj=proj, parallel=parallel, cores=cores, verbose=verbose)
    if(verbose) cat("Rescaling Ghosh distance\n")
    rgdist <- matrix(NA,nrow=length(x),ncol=length(y))
    if(identical(x,y)){
      for(i in 1:length(x)){
        for(j in 1:i){
          rgdist[i,j] <- gdist[i,j]-0.5*(gdist[i,length(y)+1]+gdist[length(x)+1,j])
          rgdist[j,i] <- rgdist[i,j]
        }
      }
    }else{
      for(i in 1:length(x)){
        for(j in 1:length(y)){
          rgdist[i,j] <- gdist[i,j]-0.5*(gdist[i,length(y)+1]+gdist[length(x)+1,j])
        }
      }
    }
    units(rgdist) <- units(gdist) #units are lost by some operations
    return(rgdist)
  }
  if(method=="points"){
    if(!all(st_geometry_type(x) == c("POINT") & st_geometry_type(y) == c("POINT"))) stop("Geometry type has to be a POINT if method 'points' is used.")
    mdist <- st_distance(x,y)
    return(mdist)
  }
  if(method=="centroids"){
    if(!all(st_geometry_type(x) == c("POLYGON") & st_geometry_type(y) == c("POLYGON"))) stop("Geometry type has to be a POLYGON if method 'centroids' is used.")
    mdist <- hdist(x=st_centroid(x), y=st_centroid(y), method="points")
    return(mdist)
  }
  if(method=="combined") stop("The class of 'x' and 'y' has to be 'transfR' if method 'combined' is used.")
}

#' @name hdist
#' @export
hdist.sf <- function(x, y, ...){
  mdist <- hdist(x=st_geometry(x), y=st_geometry(y), ...)
  return(mdist)
}

#' @name hdist
#' @export
hdist.stars <- function(x, y, ...){
  mdist <- hdist(x=st_geometry(x), y=st_geometry(y), ...)
  return(mdist)
}

#' @name hdist
#' @export
hdist.transfR <- function(x, y, method="rghosh", weightO=0.8, weightC=0.2, ...){
  if(method=="combined"){
    if(!"outlet"%in%names(x)) stop("Outlet attribute is missing in 'x' argument. See as_transfr().")
    if(!"outlet"%in%names(y)) stop("Outlet attribute is missing in 'y' argument. See as_transfr().")
    if(!"centroid"%in%names(x)) stop("Centroid attribute is missing in 'x' argument. See as_transfr().")
    if(!"centroid"%in%names(y)) stop("Centroid attribute is missing in 'y' argument. See as_transfr().")
    odist <- hdist(x$outlet, y$outlet, method = "points", ...)
    cdist <- hdist(x$centroid, y$centroid, method = "points", ...)
    mdist <- odist*weightO + cdist*weightC
  }else{
    if(!"st"%in%names(x)) stop("Spatio-temporal arrays (st) of class stars is missing in 'x' argument. See as_transfr().")
    if(!"st"%in%names(y)) stop("Spatio-temporal arrays (st) of class stars is missing in 'y' argument. See as_transfr().")
    mdist <- hdist(x=x$st, y=y$st, ...)
  }
  return(mdist)
}
