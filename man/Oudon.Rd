\docType{data}
\encoding{UTF-8}


\name{Oudon}
\alias{Oudon}


\title{Oudon French river dataset}


\description{
'Oudon' is a dataset of the Oudon French river, part of the wider Loire Catchment. It contains all the necessary inputs to test the package and perform discharge prediction at the outlet of six sub-catchments:
\itemize{
  \item M3771810: Oudon at Châtelais (734 km²)
  \item M3774010: Chéran at la Boissière (85 km²)
  \item M3823010: Verzée at Bourg-d'Iré (205 km²)
  \item M3834030: Argos at Sainte-Gemmes-d'Andigné (153 km²)
  \item M3851810: Oudon at Segré (1310 km²)
  \item M3711810: Oudon at Cossé-le-Vivien (133 km²)
}
Hourly discharge observations of the six sub-catchments (Oudon French river) are provided from 2019-12-01 to 2020-03-01, and extracted from the French HYDRO database (http://www.hydro.eaufrance.fr). Catchment delineations and respective maps of hydraulic length have been extracted from a digital elevation model of 100 m resolution.
}


\format{'Oudon' is a list of two objects:
\itemize{
  \item hl: A list of stars objects containing the six rasters maps of hydraulic length.
  \item obs: A stars object with two dimensions (time and space, with catchment delineations as spatial support) and one attribute (discharge observations).
}}


\source{http://www.hydro.eaufrance.fr}


\keyword{datasets}
