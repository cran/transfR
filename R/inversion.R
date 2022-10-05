#' Estimate net rainfall by inversion
#'
#' Estimate net rainfall by inverse modelling, where the model is a convolution between net rainfall
#' and a unit hydrograph in order to simulate discharge.
#' @name inversion
#' @param Qobs discharge vector or object of class \code{transfR}. If no unit is provided,
#' \code{Qobs} is assumed to be in [mm/h]
#' @param uh unit hydrograph vector
#' @param deltat time step of the time series. If no unit is provided, \code{deltat} is assumed to be in [min]
#' @param RnAp net rainfall a priori. If no unit is provided, \code{RnAp} is assumed to be in [mm/h]
#' @param Bd parameter used to maintain a minimum value of standart deviation for low discharge values.
#' If no unit is provided, \code{Bd} is assumed to be in [mm/h]
#' @param Dd decorrelation time of discharge errors. If no unit is provided, \code{Dd} is assumed to be
#' in [h]
#' @param Bp parameter used to maintain a minimum value of standart deviation for low net rainfall values.
#' If no unit is provided, \code{Bp} is assumed to be in [mm/h]
#' @param Tp decorrelation time of net rainfall errors. If no unit is provided, \code{Tp} is assumed to
#' be in [h]
#' @param Ad parameter equivalent to the coefficient of variation of the discharge measurement error. If
#' no unit is provided, \code{Ad} is assumed to be dimensionless
#' @param Ap parameter equivalent to the coefficient of variation of the net rainfall error. If no unit
#' is provided, \code{Ap} is assumed to be dimensionless
#' @param warmup length of the warmup period. If no unit is provided, \code{warmup} is assumed to be in [days]
#' @param cooldown length of the period removed at the end of the simulation. If no unit is provided,
#' \code{cooldown} is assumed to be in [days]
#' @param dosplit boolean, if true the inversion is performed by
#' subperiods of length defined by \code{split}
#' @param split length the subperiods if dosplit is true. If no unit is provided, \code{split} is assumed to be
#' in [days]
#' @param fixedpar boolean, if false Ap and Ad are calibrated dynamically according to the coefficient of variation of
#' RnAp and Qobs respectively (see details)
#' @param parallel boolean, if true the splitting of the inversion
#' by subperiods is parallelised
#' @param cores the number of cores to use for parallel execution if \code{parallel} is TRUE. If not specified, the number of cores is set to the value of \code{parallel::detectCores()}
#' @param verbose boolean indicating if information messages should be written to the console
#' @param ... further arguments passed to or from other methods
#' @return An object of the same class of \code{Qobs}. If \code{Qobs} is a transfR object,
#' the same transfR object incremented by the new computed attributes.
#' @import sf stars doParallel foreach
#' @importFrom units set_units as_units drop_units
#' @importFrom stats na.omit sd
#' @details In a convolution between the unit hydrograph (\code{uh}) and net rainfall that is simulating
#' streamflow at the outltet (\code{Qobs}), and where net rainfall is the only unknown variable, this function estimates
#' net rainfall by inversion \insertCite{Tarantola1982,Menke1989,Boudhraa2018}{transfR}. It requires an
#' a priori on this net rainfall (that could be estimated by the function \link{rapriori}), a description
#' of the errors on the discharge (\code{Ad}, \code{Bd}, \code{Dd}) and on the net rainfall (\code{Ap},
#' \code{Bp}, \code{Tp}) that are assumed to be Gaussian and unbiased. Default values of these parameters
#' are taken from \insertCite{deLavenne2016;textual}{transfR}. If \code{fixedpar} is deactivated, \code{Ap}
#' is estimated at 20% of the coefficient of variation of RnAp, and \code{Ad} estimated at 5% of the coefficient
#' of variation of Qobs.
#'
#' It is recommanded to use \code{warmup} and \code{cooldown} periods in order to reduce the problem of oscillations
#' created by inversion.
#'
#' If \code{object} is provided, results are stored as a new space-time attribute in the \code{object}
#' called "RnAp".
#' @seealso
#' \link{rapriori}
#' @references
#' \insertRef{Boudhraa2018}{transfR}
#'
#' \insertRef{deLavenne2016}{transfR}
#'
#' \insertRef{Menke1989}{transfR}
#'
#' \insertRef{Tarantola1982}{transfR}
#' @examples
#' \donttest{data(Oudon)
#' icatch <- 1 # Catchment index
#' itime <- 1:1000 # Using the first values for a quicker example
#' Qobs <- Oudon$obs[["Qobs"]][itime,icatch]
#' Qspec <- units::set_units(Qobs/st_area(st_geometry(Oudon$obs)[icatch]), "mm/h")
#' deltat <- units::set_units(1, "h")
#' uc <- velocity(hl = Oudon$hl[[icatch]])
#' uh <- uh(hl = Oudon$hl[[icatch]], uc = uc, deltat = units::set_units(1,"h"))$prob
#' RnAp <- rapriori(Qobs = Qspec, lagtime = lagtime(hl = Oudon$hl[[icatch]], uc = uc),
#' deltat = deltat)
#' RnInv <- inversion(Qobs = Qspec, RnAp = RnAp, uh = uh, deltat = deltat)}
#' @export
inversion <- function(Qobs,...) UseMethod("inversion")

#' @name inversion
#' @export
inversion.default <- function(Qobs, uh, RnAp, deltat, ...){
  #--- Assumed units
  Qobs <- units::set_units(Qobs,"mm/h")
  uh <- units::set_units(uh,1)
  RnAp <- units::set_units(RnAp,"mm/h")
  deltat <- units::set_units(deltat,"min")
  #--- Inversion
  inversion.units(Qobs = Qobs, uh = uh, RnAp = RnAp, deltat = deltat, ...)
}

#' @name inversion
#' @export
inversion.units <- function(Qobs, uh, RnAp, deltat, Bd = 0.01, Dd = 1, Bp = 0.001, Tp = 20, Ad = 0.01,
                            Ap = 0.9, warmup = 10, cooldown = 8, dosplit = TRUE, split = 30,
                            fixedpar = TRUE, parallel = FALSE, cores = NULL, ...){

  #--- Define some parameter values according to inputs characteristics
  if(!fixedpar){
    if(missing(Ad)) Ad <- units::drop_units(stats::sd(Qobs,na.rm = T)/mean(Qobs,na.rm = T)*0.05)
    if(missing(Ap)) Ap <- units::drop_units(stats::sd(RnAp,na.rm = T)/mean(RnAp,na.rm = T)*0.2)
  }

  #--- Inversion parameters
  Ad <- units::set_units(Ad,1)
  Bd <- units::set_units(Bd,"mm/h")
  Ap <- units::set_units(Ap,1)
  Bp <- units::set_units(Bp,"mm/h")
  Dd <- units::set_units(Dd,"h")
  Tp <- units::set_units(Tp,"h")

  #--- Warmup and cooldown periods
  deltat   <- units::set_units(deltat,"min")
  warmup   <- units::set_units(warmup,"days")
  cooldown <- units::set_units(cooldown,"days")
  split    <- units::set_units(split,"days")
  npdt_warmup   <- ceiling(units::drop_units(warmup/deltat))
  npdt_cooldown <- ceiling(units::drop_units(cooldown/deltat))
  npdt_split    <- ceiling(units::drop_units(split/deltat))

  #--- Check units
  if(length(units(Qobs/units::set_units(1,"m"))$numerator)!=0) stop("Qobs should be a specific discharge. Check units.")
  if(length(units(RnAp/units::set_units(1,"m"))$numerator)!=0) stop("RnAp should be a specific discharge. Check units.")
  if(!(length(units(Ad)$numerator)==0 & length(units(Ad)$denominator)==0)) stop("Ad should be without units.")
  if(!(length(units(Bd/units::set_units(1,"mm/h"))$numerator)==0 & length(units(Bd/units::set_units(1,"mm/h"))$denominator)==0)) stop("Bd should be a specific discharge. Check units.")
  if(!(length(units(Ap)$numerator)==0 & length(units(Ap)$denominator)==0)) stop("Ap should be without units.")
  if(!(length(units(Bp/units::set_units(1,"mm/h"))$numerator)==0 & length(units(Bp/units::set_units(1,"mm/h"))$denominator)==0)) stop("Bp should be a specific discharge. Check units.")
  if(!(length(units(Dd/units::set_units(1,"h"))$numerator)==0 & length(units(Dd/units::set_units(1,"h"))$denominator)==0)) stop("Dd should be a time. Check units.")
  if(!(length(units(Tp/units::set_units(1,"h"))$numerator)==0 & length(units(Tp/units::set_units(1,"h"))$denominator)==0)) stop("Tp should be a time. Check units.")

  #--- Managing periods with NA
  bna <- !is.na(RnAp)&!is.na(Qobs)
  nna <- which(bna)
  if(length(nna)<=(npdt_warmup+npdt_cooldown)) return(rep(NA,length(Qobs)))

  #--- Detect if several periods without NA and run individual inversions for each period
  period <- as.numeric(bna[1])
  for(i in 1:(length(bna)-1)) period[i+1] <- period[i]+abs(bna[i]-bna[i+1])
  if(max(period)>2){
    Rn <- NA
    for(p in unique(period[bna])) Rn[period==p] <- inversion.units(Qobs = Qobs[period==p], uh = uh, RnAp = RnAp[period==p], deltat = deltat, Bd = Bd, Dd = Dd, Bp = Bp, Tp = Tp, Ad = Ad, Ap = Ap, warmup = warmup, cooldown = cooldown, dosplit = dosplit, split = split, fixedpar = TRUE, parallel = parallel)
    return(Rn)
  }

  #--- Checking
  if(any((nna[2:length(nna)]-nna[2:length(nna)-1])>1)) stop("Times series must be continuous. NA values can be handle only at the beginning or the end of the time series.")
  # if(sum(!bna)>0) warning(paste0(sum(!bna)," NA values at the beginning or end of Qobs or RnAp times series. This part of the time series will not be used."))
  if(warmup<units::set_units(10,"days")) warning("Warmup period is short (< 10 days), you might observe oscillations on Rn time series.")
  if(cooldown<units::set_units(8,"days")) warning("Cool down period is short (< 8 days), you might observe oscillations on Rn time series.")
  if(npdt_warmup<(5*length(uh))) warning("Warmup period might be too short, you might observe oscillations on Rn time series.")
  if(npdt_cooldown<(4*length(uh))) warning("Cool down period might be too short, you might observe oscillations on Rn time series.")

  #--- Vector dimension
  Q <- Qobs[nna]
  R <- RnAp[nna]
  npdtR<-length(R)
  npdtQ<-length(Q)
  ndec <- max(sapply(strsplit(as.character(Q), "\\."),FUN=function(x)(nchar(x[2]))),na.rm=T)

  #--- Default values
  ZRn <- 0
  NRn <- NA
  units(ZRn) <- units(NRn) <- units(R)

  if(dosplit&(length(Q)>npdt_split)){
    if((npdt_warmup+npdt_cooldown)>=npdt_split) stop("All the subperiod is used for initialisation only: (warmup+cooldown)>=split. Increase the value of 'split' or decrease the values of 'warmup' and 'cooldown'.")
    cuts_start <- seq(from=1,to=(length(Q)+npdt_split-npdt_warmup-npdt_cooldown),by=(npdt_split-npdt_warmup-npdt_cooldown))
    cuts_end   <- cuts_start+npdt_split-1
    cuts_end[cuts_end>length(Q)] <- length(Q)
    cuts_start <- cuts_start[!duplicated(cuts_end)]
    cuts_end   <- cuts_end[!duplicated(cuts_end)]
    if(parallel){
      if(missing(cores)|is.null(cores)) cores <- parallel::detectCores()
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl=cl)
      on.exit(parallel::stopCluster(cl))
      tmp <- foreach::"%dopar%"(foreach::foreach(i = 1:length(cuts_start), .combine='c'),
                                na.omit(inversion.units(Qobs = Q[cuts_start[i]:cuts_end[i]], uh = uh, RnAp = R[cuts_start[i]:cuts_end[i]], deltat = deltat, Bd = Bd, Dd = Dd, Bp = Bp, Tp = Tp, Ad = Ad, Ap = Ap, warmup = warmup, cooldown = cooldown, dosplit=FALSE, split = split, fixedpar = TRUE, parallel = FALSE)))
      Rn <- c(rep(NRn,npdt_warmup),tmp,rep(NRn,npdt_cooldown))
    }else{
      Rn <- rep(NRn,length(Q))
      for(i in 1:length(cuts_start)){
        tmp <- inversion.units(Qobs = Q[cuts_start[i]:cuts_end[i]], uh = uh, RnAp = R[cuts_start[i]:cuts_end[i]], deltat = deltat, Bd = Bd, Dd = Dd, Bp = Bp, Tp = Tp, Ad = Ad, Ap = Ap, warmup = warmup, cooldown = cooldown, dosplit=FALSE, split = split, fixedpar = TRUE, parallel = FALSE)
        Rn[cuts_start[i]:cuts_end[i]][!is.na(tmp)] <- tmp[!is.na(tmp)]
      }
    }
    # If NA in Qobs or RnAp
    FullRn <- rep(NRn,length(Qobs))
    FullRn[nna] <- Rn
    return(FullRn)
  }

  #--- Transfer function matrix
  M <- uh2mat(uh,nrow=npdtQ,ncol=npdtR)

  #------------------------------------
  # Covariance matrix on runoff
  #------------------------------------
  #--- Standart deviation (confidence to measurement)
  TQobs = Ad * Q + Bd
  #--- Matrix computation
  MTQobs1<-matrix(rep(TQobs,npdtQ),npdtQ,npdtQ)
  MTQobs2<-matrix(0,npdtQ,npdtQ)
  units(MTQobs1) <- units(MTQobs2) <- units(TQobs)
  for(i in 1:npdtQ){MTQobs2[i,i]<-TQobs[i]}
  u<-npdtQ-1
  v<-0:u
  Mtmp1<-matrix(rep(0:u,npdtQ),npdtQ,npdtQ)
  Mtmp2<-matrix(rep(0:-u,npdtQ),npdtQ,npdtQ,byrow=TRUE)
  MTQobs3<-Mtmp1+Mtmp2

  CovQ_tmp <- -0.5*((abs(MTQobs3)*deltat)/Dd)^2
  if(as.character(units(CovQ_tmp))=="1") CovQ_tmp <- units::drop_units(CovQ_tmp) else stop("Units are not consistent. Check Qobs, deltat and Dd.")
  CovQ<-MTQobs1%*%MTQobs2*exp(CovQ_tmp)
  #------------------------------------
  # Covariance matrix on net rainfall
  #------------------------------------
  #--- Standart deviation (confidence to apriori rainfall)
  TRnAp = Ap * R + Bp
  #--- Matrix computation
  MTRnAp1<-matrix(rep(TRnAp,npdtR),npdtR,npdtR)
  MTRnAp2<-matrix(0,npdtR,npdtR)
  units(MTRnAp1) <- units(MTRnAp2) <- units(TRnAp)
  for(i in 1:npdtR){MTRnAp2[i,i]<-TRnAp[i]}
  u<-npdtR-1
  v<-0:u
  Mtmp1<-matrix(rep(0:u,npdtR),npdtR,npdtR)
  Mtmp2<-matrix(rep(0:-u,npdtR),npdtR,npdtR,byrow=TRUE)
  MTRnAp3<-Mtmp1+Mtmp2

  CovAp_tmp <- -0.5*((abs(MTRnAp3)*deltat)/Tp)^2
  if(as.character(units(CovAp_tmp))=="1") CovAp_tmp <- units::drop_units(CovAp_tmp) else stop("Units are not consistent. Check RnAp, deltat and Tp.")
  CovAp<-MTRnAp1%*%MTRnAp2*exp(CovAp_tmp)

  #------------------------------------
  # Solution: net rainfall estimation
  #------------------------------------
  Rn <- R + CovAp%*%t(M) %*% solve((M %*% CovAp %*% t(M) + CovQ)) %*% (Q - M%*%R)
  #--- Removing negative Rn, warmup period and cooldown period
  Rn[Rn<ZRn] <- ZRn
  if(npdt_warmup>0)   Rn[1:npdt_warmup] <- NRn
  if(npdt_cooldown>0) Rn[(length(Rn)-npdt_cooldown+1):length(Rn)] <- NRn
  # If NA in Qobs or RnAp
  FullRn <- rep(NRn,length(Qobs))
  FullRn[nna] <- round(Rn,ndec)
  return(FullRn)
}


#' @name inversion
#' @export
inversion.transfR <- function(Qobs, verbose=TRUE, ...){
  object <- Qobs
  if(!"deltat"%in%names(object)) stop("Time step attribute (deltat) is missing in transfR object. See as_transfr().")
  if(!"Qobs"%in%names(object$st)) stop("Discharge observation (Qobs) is missing in the spatio-temporal arrays (st) of in transfR object. See as_transfr().")
  if(!"RnAp"%in%names(object$st)) stop("A priori on net rainfall (RnAp) is missing in the spatio-temporal arrays (st) of in transfR object. See rapriori().")
  if(object$deltat==units::set_units(1,"day")){out_units <- units::as_units("mm/d")}else{out_units <- units::as_units("mm/h")}
  object$st$RnInv <- NA
  for(i in 1:dim(object$st)[2]){
    if(verbose) progress("Computing inversion for catchment ",i,dim(object$st)[2])
    Qspec <- object$st$Qobs[,i]/st_area(st_geometry(object$st)[i])
    Qspec <- units::set_units(Qspec,out_units,mode="standart")
    RnInv <- inversion(Qobs = Qspec, uh = object$uh[[i]], RnAp = object$st$RnAp[,i], deltat = object$deltat, ...)
    object$st$RnInv[,i] <- RnInv
  }
  object$st[["RnInv"]] <- units::set_units(object$st[["RnInv"]],out_units,mode="standart") # could not find a way to keep the units provided by inversion
  return(object)
}



