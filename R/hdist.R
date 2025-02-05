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
#' @param parallel logical indicating whether the computation should be parallelised. Could be a vector
#' of length two to distinguish between parallelisation of sampling and Ghosh distance
#' (because sampling over large areas can be memory intensive)
#' @param cores the number of cores to use for parallel execution if \code{parallel} is TRUE.
#' If not specified, the number of cores is set to the value of \code{parallel::detectCores()}.
#' Similarly to \code{parallel}, it could be a vector of length two to distinguish between
#' parallelisation of sampling and Ghosh distance
#' @param verbose boolean indicating if information messages should be written to the console
#' @param ... further arguments passed to or from other methods
#' @return A matrix of class units with the catchments of \code{x} organised in rows
#' and the catchments of \code{y} organised in columns.
#' @details The \code{method} "ghosh" refers to a simplification of the distance defined
#' by \insertCite{Ghosh1951;textual}{transfR} as proposed by \insertCite{Gottschalk1993,Gottschalk2011;textual}{transfR}.
#' The rescaled Ghosh distance (\code{method} "rghosh") is calculated following \insertCite{deLavenne2016;textual}{transfR}.
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
#' hdist(x = catchments[1:2], y = catchments[3:5], gres = 5, method = "rghosh",
#' parallel=c(FALSE, TRUE), cores=2)
#' @import sf stars doParallel foreach
#' @importFrom units set_units drop_units
#' @useDynLib transfR, .registration = TRUE
#' @export
hdist <- function(x, y, ...) UseMethod("hdist")

#' @name hdist
#' @export
hdist.sfc <- function(x, y, method="rghosh", gres=5, ditself=FALSE, maxsample=2.5e4, proj=NULL, parallel=FALSE, cores=NULL, verbose=TRUE, ...){
  method <- match.arg(method,c("ghosh","rghosh","rghosh2","points","centroids","combined"),several.ok = FALSE)
  if(missing(proj)){
    if(!identical(st_crs(x),st_crs(y))) stop("x and y do not use the same coordinate reference system.")
    if(is.na(st_crs(x)$input)){stop("Don't know how to compute distances because no coordinate reference system was found. Define a coordinate reference system (using st_crs()) or define 'proj' argument.")}
    if(st_crs(x)$input=="EPSG:4326"){proj=FALSE}else{proj=TRUE}
  }
  if(method%in%c("ghosh","rghosh2")){
    if(!(all(st_geometry_type(x) == c("POLYGON")) & all(st_geometry_type(y) == c("POLYGON")))) stop("Geometry type has to be a POLYGON if method 'rghosh' is used.")
    if(verbose) cat("Sampling the catchments at a resolution of about",gres,"pts/km2\n")
    if(length(parallel)==1) parallel[1:2] <- parallel
    if(any(parallel) & (missing(cores)|is.null(cores))) cores <- parallel::detectCores()
    if(length(cores)==1) cores[1:2] <- cores
    if(parallel[1]){
      cl <- parallel::makeCluster(cores[1])
      doParallel::registerDoParallel(cl=cl)
      on.exit(parallel::stopCluster(cl))
    }
    xy_union <- st_union(c(x,y))
    size <- round(units::set_units(gres,"1/km^2")*units::set_units(st_area(xy_union),"km^2"))
    xydisc <- st_sample(xy_union, size=units::drop_units(size), type="regular")
    if(parallel[1]){
      xdisc <- foreach::"%dopar%"(foreach::foreach(i=1:length(x),.packages="sf"), xydisc[x[i]])
    }else{
      xdisc <- list()
      for(i in 1:length(x)) xdisc[[i]] <- xydisc[x[i]]
    }
    if(identical(x,y)){
      ydisc <- xdisc
      for(i in 1:length(xdisc)) if(length(xdisc[[i]])>maxsample){
        xdisc[[i]] <- ydisc[[i]] <- xdisc[[i]][sort(sample.int(length(xdisc[[i]]), size = maxsample, replace = FALSE))]
        if(verbose) cat("Random resampling of catchment",i,"with a maximum size sample of",maxsample,"points.\n")}
    }else{
      if(parallel[1]){
        ydisc <- foreach::"%dopar%"(foreach::foreach(i=1:length(y),.packages="sf"), xydisc[y[i]])
      }else{
        ydisc <- list()
        for(i in 1:length(y)) ydisc[[i]] <- xydisc[y[i]]
      }
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
    gdist <- loop_gdist(xdisc=xdisc, ydisc=ydisc, proj=proj, intersect=FALSE, parallel=parallel[2], cores=cores[2])
    if(ditself){
      if(verbose) cat("Computing Ghosh distance within catchments\n")
      if(!identical(x,y)){
        gdist <- cbind(gdist,sapply(xdisc,FUN=function(x){call_gdist(pts1=x, pts2=x, proj=proj, intersect=FALSE, rescale=FALSE, diag=TRUE, parallel=parallel[2], cores=cores[2])}))
        gdist <- rbind(gdist,c(sapply(ydisc,FUN=function(x){call_gdist(pts1=x, pts2=x, proj=proj, intersect=FALSE, rescale=FALSE, diag=TRUE, parallel=parallel[2], cores=cores[2])}),NA))
      }else{
        gdist <- cbind(gdist,diag(gdist))
        gdist <- rbind(gdist,c(diag(gdist),NA))
      }
    units(gdist) <- units(st_distance(xdisc[[1]][1],ydisc[[1]][1])) #units are lost by some operations
    }
    if(method=="rghosh2"){
      gc() # Clean memory
      if(verbose) cat("Computing Ghosh distance within the shared areas\n")
      idist <- loop_gdist(xdisc=xdisc, ydisc=ydisc, proj=proj, intersect=TRUE, parallel=parallel[2], cores=cores[2])
      if(verbose) cat("Rescaling Ghosh distance\n")
      gdist <- gdist-idist
    }
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

call_gdist <- function(pts1, pts2, proj, intersect, rescale, diag, parallel, cores){
  if(!parallel) cores = 1
  if(intersect){
    # Compute the Ghosh distance matrix for common points only
    iseq <- st_equals(pts1, pts2, sparse = FALSE)
    pts_inter <- pts1[apply(iseq, MARGIN = 1, FUN = any)]
    if(length(pts_inter)>0){
      call_gdist(pts1=pts_inter, pts2=pts_inter, proj=proj, intersect=FALSE, rescale=FALSE, diag=FALSE, parallel=parallel, cores=cores)
    }else{
      return(0)
    }
  }else{
    # Compute the Ghosh distance matrix between all points
    .Fortran("gdist", coord1=st_coordinates(pts1), coord2=st_coordinates(pts2), n1=length(pts1), n2=length(pts2),
             proj=proj, rescale=rescale, diag=diag, nthreads=as.integer(cores), mdist=0)$mdist
  }
}

loop_gdist <- function(xdisc, ydisc, proj, intersect, parallel, cores){
  if(!parallel) cores = 1
  gdist <- matrix(NA,nrow=length(xdisc),ncol=length(ydisc))
  if(!identical(xdisc,ydisc)){
    for(i in 1:length(xdisc)) gdist[i,] <- sapply(ydisc,FUN=function(x){call_gdist(pts1=x, pts2=xdisc[[i]], proj=proj, intersect=intersect, rescale=FALSE, diag=FALSE, parallel=parallel, cores=cores)})
  }else{
    for(i in 1:length(xdisc)) gdist[i,i:length(xdisc)] <- gdist[i:length(xdisc),i] <- sapply(ydisc[i:length(xdisc)],FUN=function(x){call_gdist(pts1=x, pts2=xdisc[[i]], proj=proj, intersect=intersect, rescale=FALSE, diag=FALSE, parallel=parallel, cores=cores)})
  }
  units(gdist) <- units(st_distance(xdisc[[1]][1],ydisc[[1]][1])) #units are lost by some operations
  return(gdist)
}
