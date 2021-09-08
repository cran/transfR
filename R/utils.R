progress <- function(txt1 = "", i, imax, txt2 = ""){
  cat(txt1, i, "/", imax, txt2, paste(rep("\b", nchar(txt1)+nchar(i)+1+nchar(imax)+nchar(txt2)), collapse = ""),sep = "")
}

uh2mat <- function(uh,nrow,ncol){
  mat <- units::set_units(matrix(0,nrow,ncol),"1")
  for(c in 1:ncol){
    r <- min(c(nrow,c+length(uh)-1))
    mat[c:r,c] <- uh[1:(r-c+1)]
  }
  return(mat)
}

solve.units <- function(a, ...){
  x <- drop_units(a)
  y <- solve(x)
  units(y) <- units(1/a)
  return(y)
}

#' @method %*% transfR
'%*%' <- function(x,...) UseMethod('%*%',x)

'%*%.default' <- .Primitive('%*%')

'%*%.units' <-  function(e1, e2){
  x1 <- drop_units(e1)
  x2 <- drop_units(e2)
  y <- x1 %*% x2
  units(y) <- structure(list(numerator = c(units(e1)$numerator,units(e2)$numerator),
                             denominator = c(units(e1)$denominator,units(e2)$denominator)),
                        class = "symbolic_units")
  return(y)
}
