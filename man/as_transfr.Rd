% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_transfr.R
\name{as_transfr}
\alias{as_transfr}
\title{Create transfR object}
\usage{
as_transfr(
  object,
  st,
  uc,
  lagtime,
  surface,
  delineation,
  outlet,
  centroid,
  uh,
  hl
)
}
\arguments{
\item{object}{object of class \code{transfR}}

\item{st}{spatio-temporal arrays of class \code{stars}. Observed discharge must be described by the column name 'Qobs'.
Time should be the first dimension, space the second dimension. If no unit is provided, Qobs is assumed to be in [m3/s]
and RnInv is assumed to be in [mm/h] (or [mm/d] at daily time step).}

\item{uc}{vector of the streamflow velocities of the catchments. If no unit is provided, \code{uc} is assumed to be in [m/s].}

\item{lagtime}{vector of the lag times of the catchments. If no unit is provided, \code{lagtime} is assumed to be in [h].}

\item{surface}{vector of the surfaces of the catchments. If no unit is provided, \code{surface} is assumed to be in [km2].}

\item{delineation}{spatial layer of the boundary of the catchments of class \code{sfc_POLYGON}.}

\item{outlet}{spatial layer of the outlets of the catchments of class \code{sfc_POINT}.}

\item{centroid}{spatial layer of the centroids of the catchments of class \code{sfc_POINT}.}

\item{uh}{list of the unit hydrographs of the catchments.}

\item{hl}{hydraulic length of class \code{stars}, \code{matrix} or \code{vector}. If no unit is provided, \code{hl}
is assumed to be in [m]. See details below.}
}
\value{
An object of class transfR.
}
\description{
Create a transfR object or add new attributes to a transfR object.
}
\details{
This function creates an object of class \code{transfR} or increment an existing \code{transfR} object
with new attributes. It can be used to gather and organize most of the inputs and outputs of the other functions
like streamflow velocities, unit hydrograph, a priori on net rainfall, inversions and simulations of every catchments.

This function can be used to organise the two user inputs required for a conventional use of the package, namely \code{st}
and \code{hl}. The hydraulic lengths are defined as the flow path length from each pixel to the outlet within
the river network \insertCite{Cudennec2004,Aouissi2013}{transfR}. Catchment delineations and hydraulic lengths
need to be prepared beforehand by the user. This package does not provide functions to create them.
However, several GIS software offer possibilities to extract them from a digital elevation model
such as GRASS toolkits \insertCite{Jasiewicz2011}{transfR},
Whitebox GAT (\insertCite{@see @Lindsay2016;textual}{transfR} or \href{https://github.com/jblindsay/whitebox-tools}{WhiteboxTools}),
TauDEM (D. Tarboton, Utah State University)
or online services (\insertCite{@see @Squividant2015;textual}{transfR} for catchment delineation in the Brittany French region).
}
\examples{
data(Oudon)
object <- as_transfr(st = Oudon$obs, hl = Oudon$hl)
}
\references{
\insertRef{Aouissi2013}{transfR}

\insertRef{Cudennec2004}{transfR}

\insertRef{Jasiewicz2011}{transfR}

\insertRef{Lindsay2016}{transfR}

\insertRef{Squividant2015}{transfR}
}
