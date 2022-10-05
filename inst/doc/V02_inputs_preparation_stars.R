## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=TRUE, message=FALSE, results='hide', eval=TRUE---------------------
library(transfR)
data(Oudon)

wd <- tempdir(check = TRUE)
st_write(st_sf(ID = paste0("ID", 1:6), geom = st_geometry(Oudon$obs)), 
         dsn = file.path(wd, "catchments.shp"), delete_layer = TRUE)
write.table(data.frame(DateTime = format(st_get_dimension_values(Oudon$obs,1),
                                         "%Y-%m-%d %H:%M:%S"), 
                       ID1 = Oudon$obs$Qobs[,1], 
                       ID2 = Oudon$obs$Qobs[,2], 
                       ID3 = Oudon$obs$Qobs[,3], 
                       ID4 = Oudon$obs$Qobs[,4], 
                       ID5 = Oudon$obs$Qobs[,5], 
                       ID6 = Oudon$obs$Qobs[,6]), 
            file = file.path(wd, "discharge.txt"), 
            col.names = TRUE, row.names = FALSE, sep = ";", quote = FALSE)

## ---- echo=TRUE, message=FALSE, results='hide', eval=TRUE---------------------
library(sf)
catchments <- st_read(file.path(wd, "catchments.shp"), "catchments", stringsAsFactors = FALSE)
obs_sf <- catchments[1:5,] # Gauged catchments
sim_sf <- catchments[6,]   # Ungauged catchments

## ---- echo=TRUE, message=FALSE, results='hide', eval=TRUE---------------------
library(units)
Q <- read.table(file.path(wd, "discharge.txt"), header = TRUE, sep = ";", 
                colClasses = c("character", rep("numeric", 6)))
Qmatrix  <- as.matrix(Q[,-1])
Qmatrix  <- set_units(Qmatrix, "m^3/s")

## ---- echo=TRUE, message=FALSE, results='hide', eval=TRUE---------------------
library(stars)
Qmatrix  <- Qmatrix[,obs_sf$ID] #to have the same order as in the spacial data layer
obs_st   <- st_as_stars(list(Qobs = Qmatrix), 
                            dimensions = st_dimensions(time = as.POSIXct(Q$DateTime, tz="UTC"), 
                                                       space = obs_sf$geometry))
sim_st   <- st_as_stars(dimensions = st_dimensions(time = as.POSIXct(Q$DateTime, tz="UTC"), 
                                                       space = sim_sf$geometry))

## ---- echo=TRUE, message=FALSE, results='hide', eval=TRUE---------------------
obs <- as_transfr(st = obs_st, hl = Oudon$hl[1:5])
sim <- as_transfr(st = sim_st, hl = Oudon$hl[6])

## ---- echo=TRUE, message=FALSE, results='hide', eval=TRUE---------------------
sim <- quick_transfr(obs, sim, parallel = TRUE, cores = 2)

## ---- echo=TRUE, message=TRUE, eval=TRUE--------------------------------------
sim$st

## ---- echo=FALSE, message=FALSE, warning=FALSE, eval=TRUE, results='hide'-----
# Cleaning temporary directory
unlink(wd, recursive = TRUE)

