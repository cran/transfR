\docType{data}
\encoding{UTF-8}


\name{Blavet}
\alias{Blavet}


\title{Blavet French river dataset}


\description{
'Blavet' is a dataset of the Blavet French river in Brittany peninsula and two neighouring rivers (Claie and Coët-Organ). It contains all the necessary inputs to test the package and perform discharge prediction at the outlet of six catchments:
\itemize{
  \item {J5613010}  {Evel at Guénin (316 km²)}
  \item {J5618310}  {Fremeur et Guénin (15.1 km²)}
  \item {J5618320}  {Fremeur et Pluméliau (5.88 km²)}
  \item {J5704810}  {Coët-Organ at Quistinic (47.7 km²)}
  \item {J8433020}  {Claie at Saint-Jean-Brévelay (135 km²)}
  \item {AgrHys}  {Coët-Dan at Naizin (4.9 km²)}
}
Hourly discharge observations of the six catchments are provided for one hydrological year, from 2013-10-01 to 2014-10-01. It has been extracted from the French HYDRO database (http://www.hydro.eaufrance.fr). Discharge observations for the Coët-Dan river is provided by the AgrHys Environment Research Observatories \insertCite{Fovet2018}{transfR} managed by INRAE (https://www6.inrae.fr/ore_agrhys_eng). Catchment delineations and respective maps of hydraulic length have been extracted from a digital elevation model of 100 m resolution.
}


\format{'Blavet' is a list of two objects:
\itemize{
  \item {hl}  {A list of stars objects containing the six rasters maps of hydraulic length.}
  \item {obs}  {A stars object with two dimensions (time and space, with catchment delineations as spatial support) and one attribute (discharge observations).}
}}


\source{
http://www.hydro.eaufrance.fr

https://www6.inrae.fr/ore_agrhys_eng
}


\references{
\insertRef{Fovet2018}{transfR}
}


\keyword{datasets}
