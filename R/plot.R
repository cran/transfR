#' Plot transfR object
#'
#' Plot transfR object.
#' @param x transfR object
#' @param y ignored
#' @param i spatial index to plot
#' @param attribute attribute of the transfR object to plot
#' @param format format for labels of time series on x axis
#' @param at a date-time or date object for ticks on x axis
#' @param nticks number of ticks on x axis
#' @param main a main title for the plot, see also \link[graphics]{title}
#' @param xlab a label for the x axis, defaults to a description of x
#' @param ylab a label for the y axis, defaults to a description of y
#' @param type 1-character string giving the type of plot desired (for details,
#' see \link[graphics]{plot})
#' @param lwd the line width (for details, see \link[graphics]{par})
#' @param las the style of axis labels (for details, see \link[graphics]{par})
#' @param cex.names expansion factor for axis names (for details, see \link[graphics]{barplot})
#' @param col a specification for the default plotting color (for details, see \link[graphics]{par})
#' @param keeplocal boolean to preserve local graphical parameters
#' @param ... further specifications, see \link[graphics]{plot}
#' @method plot transfR
#' @import sf stars graphics
#' @examples
#' data(Oudon)
#' object <- as_transfr(st=Oudon$obs,hl=Oudon$hl)
#' plot(object,attribute="Qobs")

#' @name plot
#' @export
plot.transfR <- function(x, y, i, attribute,
                         main = sprintf("Catchment %i",i),
                         xlab,
                         ylab,
                         format,
                         at,
                         nticks = 5,
                         type = "l",
                         lwd = 2,
                         las = 1,
                         cex.names = 1,
                         col = c("#045a8d","#fb8072","#bebada","#ffffb3","#8dd3c7"),
                         keeplocal = TRUE,
                         ...){

  # save default graphical & time parameters and resetting on exit
  if(keeplocal){
    opar <- par(no.readonly = TRUE)
    olctime <- Sys.getlocale("LC_TIME")
    if(Sys.info()[['sysname']]=="Windows"){
      Sys.setlocale("LC_TIME", "English")
    }else{
      Sys.setlocale("LC_TIME", "en_US.UTF-8")
    }
    on.exit({
      par(opar)
      Sys.setlocale("LC_TIME", olctime)
    })
  }
  if(any(!attribute %in% c(names(x),names(x$st)))) stop("Attribute not found in the transfR object.")
  if(missing(i)) i <- 1:dim(x$st)[2]
  if(any(i>dim(x$st)[2]|i<1)) stop(paste0("'i' can not be greater than the number of catchments (",dim(x$st)[2],")."))
  if(length(i)>1){
    for(j in i){
      if(j!=i[1]) par(ask=T)
      if(missing(format)){
        plot(x = x, i = j, attribute = attribute, xlab = xlab, ylab = ylab,
             type = type, col = col, keeplocal = FALSE, ...)
      }else{
        plot(x = x, i = j, attribute = attribute, xlab = xlab, ylab = ylab,
             format = format, type = type, col = col, keeplocal = FALSE, ...)
      }
    }
    return()
  }
  if(all(attribute %in% names(x$st))){
    st <- x$st
    if(missing(xlab)) xlab <- ""
    if(missing(ylab)){
      unit <- as.character(units(x$st[[attribute[1]]]))
      ylab <- parse(text=paste0(attribute[1],"~(",unit,")"))
    }
    ymax <- 0
    for(att in attribute) ymax <- max(c(ymax,max(st[[att]][,i],na.rm=T)))
    par(mar=c(3,5,3,1))
    plot(x = st_get_dimension_values(st,1), y = st[[attribute[1]]][,i],
         ylim = c(0,ymax),
         main = main, xlab = xlab, ylab = ylab,
         bty = "n", xaxt = "n", type = "n", ...)
    grid(lwd = 1)
    for(att in attribute) lines(x = st_get_dimension_values(st,1), y = st[[att]][,i],
                                col = col[which(att==attribute)], lwd = lwd)
    if(missing(at)) at <- seq(from = min(st_get_dimension_values(st,1)),
                                  to = max(st_get_dimension_values(st,1)),
                                  length = nticks)
    axis.POSIXct(side = 1, x = st_get_dimension_values(st,1), format = format, at = at)
    legend("topright", legend = c(attribute), col = col, lty = 1, lwd = lwd, bty = "n")
  }
  if(all(attribute %in% "uh")){
    uh <- x$uh[[i]]
    if(missing(main)) main <- paste("Unit hydrograph of catchment",i)
    if(missing(xlab)) xlab <- "Travel time to the outlet [h]"
    if(missing(ylab)) ylab <- "Probability"
    if(length(uh)>10){
      if(missing(las)) las <- 2
      if(missing(cex.names)) cex.names <- 0.75
    }
    tt_class <- paste0(units::set_units((1:length(uh)-1)*x$deltat,"h"),"-",units::set_units((1:length(uh))*x$deltat,"h"))
    barplot(height = matrix(uh,nrow = 1), space = 0, names.arg = tt_class,
            ylab = ylab, xlab = xlab, main = main, col = "#045a8d", border = "gray",
            las = las, cex.names = cex.names, ...)
  }
}
