# transfR: Transfer of Hydrograph from Gauged to Ungauged Catchments


## Overview

This R package aims to propose a geomorphology-based hydrological modelling to transfer streamflow measurements from gauged catchments to ungauged catchments (where there is no stations monitoring the streamflow). It follows a runoff-runoff approach, i.e. it directly combines the observed streamflow series available at monitoring stations to estimate the streamflow series anywhere else in the surroundings rivers and without the need to implement a full rainfall-runoff model. 

## Short description of the modelling approach

The hydrological modelling is based on a description of the hydro-geomorphometry of the river network which can be easily observed for any given outlet. An inversion of this model for a gauged catchment allows the observed streamflow series to be deconvoluted in order to estimate an almost scale-independent signal, namely the net rainfall (Boudhraâ et al. 2018). Transferring this estimate of the net rainfall series to a targeted ungauged catchment then allows to simulate the streamflow there. The use of streamflow observations from several gauged catchments of the neighbourhood increases the robustness of the simulation (de Lavenne et al. 2016). The methodology has first been implemented on a few catchments in semiarid Tunisia at the event time scale (Boudhraâ et al. 2009), then in dense configurations of neighbouring and nesting catchments in France with mainly temperate oceanic climate (de Lavenne et al. 2015; de Lavenne et al. 2016; de Lavenne and Cudennec 2019) and in snow-influenced Québec, Canada (Ecrepont et al. 2019).

## Installation

``` r
install.packages("transfR")
```

## Functions and objects

To implement the method, it is advised to explore the following functions in this order:

- `as_transfr` create a “transfR” database from a “stars” object and morphometric description of the catchments (hydraulic lengths)
- `velocity` estimates the main model parameter, i.e. the streamflow velocity, from different regionalisation strategies
- `uh` estimates a simple linear model, i.e. the unit hydrograph, based on the analysis of catchment geomorphology and streamflow velocity
- `rapriori` provides an a priori on the net rainfall, as needed for the model's inversion
- `inversion` estimates the net rainfall by an inverse modelling
- `hdist` computes hydrological distances between catchments, such as the rescaled Ghosh distances
- `mixr` estimates the net rainfall of one catchment by averaging the net rainfall of neighbouring gauged catchments and according to hydrological distances
- `convolution` computes the convolution of the net rainfall by the unit hydrograph to estimate streamflow

## How to get started

This package comes with two datasets (Blavet and Oudon) that contains all the necessary inputs to test the package and perform discharge prediction. Users are advised to check the 'Get started with transfR' vignette (`vignette("V01_get_started", package = "transfR")`) that provides a complete implementation of the method with the Oudon dataset. In addition, each function comes with different examples.

For the French region of Brittany, a [web service (SIMFEN)](https://geosas.fr/simfen/) using this package was developed to facilitate the implementation of the method without the need for the user to have programming skills in R or to collect the necessary input data (Dallery et al. 2020).

## References

Boudhraâ H, Cudennec C, Slimani M, Andrieu H (2009). “Hydrograph transposition between basins through a geomorphology-based deconvolution-reconvolution approach.” IAHS publication, 333, 76.

Boudhraâ H, Cudennec C, Andrieu H, Slimani M (2018). “Net rainfall estimation by the inversion of a geomorphology-based transfer function and discharge deconvolution.” Hydrological Sciences Journal, 63(2), 285–301. doi: [10.1080/02626667.2018.1425801](https://doi.org/10.1080/02626667.2018.1425801).

Dallery D, Squividant H, de Lavenne A, Launay J, Cudennec C (2020). “An end-user-friendly hydrological Web Service for hydrograph prediction in ungauged basins.” Hydrological Sciences Journal, 1–9. doi: [10.1080/02626667.2020.1797045](https://doi.org/10.1080/02626667.2020.1797045).

Ecrepont S, Cudennec C, Anctil F, Jaffrézic A (2019). “PUB in Québec: A robust geomorphology-based deconvolution-reconvolution framework for the spatial transposition of hydrographs.” Journal of Hydrology, 570, 378–392. doi: [10.1016/j.jhydrol.2018.12.052](https://doi.org/10.1016/j.jhydrol.2018.12.052).

de Lavenne A, Boudhraâ H, Cudennec C (2015). “Streamflow prediction in ungauged basins through geomorphology-based hydrograph transposition.” Hydrology Research, 46(2), 291–302. doi: [10.2166/nh.2013.099](https://doi.org/10.2166/nh.2013.099).

de Lavenne A, Skøien JO, Cudennec C, Curie F, Moatar F (2016). “Transferring measured discharge time series: Large-scale comparison of Top-kriging to geomorphology-based inverse modeling.” Water Resources Research, 52(7), 5555–5576. doi: [10.1002/2016WR018716](https://doi.org/10.1002/2016WR018716).

de Lavenne A, Cudennec C (2019). “Assessment of freshwater discharge into a coastal bay through multi-basin ensemble hydrological modelling.” Science of The Total Environment, 669, 812 - 820. ISSN 0048-9697, doi: [10.1016/j.scitotenv.2019.02.387](https://doi.org/10.1016/j.scitotenv.2019.02.387).

