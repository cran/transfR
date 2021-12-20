\name{transfR}
\alias{transfR}
\docType{package}
\encoding{UTF-8}
\title{\packageTitle{transfR}}


\description{
This R package aims to propose a geomorphology-based hydrological modelling to transfer streamflow measurements from gauged catchments to ungauged catchments (where there is no stations monitoring the streamflow). It follows a runoff-runoff approach, i.e. it directly combines the observed streamflow series available at monitoring stations to estimate the streamflow series anywhere else in the surroundings rivers and without the need to implement a full rainfall-runoff model. \cr\cr



## --- Short description of the modelling approach

The hydrological modelling is based on a description of the hydro-geomorphometry of the river network which can be easily observed for any given outlet. An inversion of this model for a gauged catchment allows the observed streamflow series to be deconvoluted in order to estimate an almost scale-independent signal, namely the net rainfall \insertCite{Boudhraa2018}{transfR}. Transferring this estimate of the net rainfall series to a targeted ungauged catchment then allows to simulate the streamflow there. The use of streamflow observations from several gauged catchments of the neighbourhood increases the robustness of the simulation \insertCite{deLavenne2016}{transfR}. The methodology has first been implemented on a few catchments in semiarid Tunisia at the event time scale \insertCite{Boudhraa2009}{transfR}, then in dense configurations of neighbouring and nesting catchments in France with mainly temperate oceanic climate \insertCite{deLavenne2015,deLavenne2016,deLavenne2019}{transfR} and in snow-influenced Québec, Canada \insertCite{Ecrepont2019}{transfR}. \cr\cr



## --- Functions and objects

To implement the method, it is advised to explore the following functions in this order:

\itemize{
\item \code{\link{as_transfr}} create a ``transfR'' database from a ``stars'' object and morphometric description of the catchments (hydraulic lengths)
\item \code{\link{velocity}} estimates the main model parameter, i.e. the streamflow velocity, from different regionalisation strategies
\item \code{\link{uh}} estimates a simple linear model, i.e. the unit hydrograph, based on the analysis of catchment geomorphology and streamflow velocity
\item \code{\link{rapriori}} provides an a priori on the net rainfall, as needed for the model's inversion
\item \code{\link{inversion}} estimates the net rainfall by an inverse modelling
\item \code{\link{hdist}} computes hydrological distances between catchments, such as the rescaled Ghosh distances
\item \code{\link{mixr}} estimates the net rainfall of one catchment by averaging the net rainfall of neighbouring gauged catchments and according to hydrological distances
\item \code{\link{convolution}} computes the convolution of the net rainfall by the unit hydrograph to estimate streamflow
}


## --- How to get started

This package comes with two datasets (\link{Blavet} and \link{Oudon}) that contains all the necessary inputs to test the package and perform discharge prediction. Users are advised to check the 'Get started with transfR' vignette (\code{vignette("V01_get_started", package = "transfR")}) that provides a complete implementation of the method with the Oudon dataset. In addition, each function comes with different examples.

For the French region of Brittany, a \href{https://geosas.fr/simfen/}{web service} using this package was developed to facilitate the implementation of the method without the need for the user to have programming skills in R or to collect the necessary input data \insertCite{Dallery2020}{transfR}.

}

\references{
\insertRef{Boudhraa2009}{transfR}

\insertRef{Boudhraa2018}{transfR}

\insertRef{Ecrepont2019}{transfR}

\insertRef{Dallery2020}{transfR}

\insertRef{deLavenne2015}{transfR}

\insertRef{deLavenne2016}{transfR}

\insertRef{deLavenne2019}{transfR}
}

\concept{model}
\concept{hydrology}
\concept{geomorphology}
\concept{ungauged catchment}
