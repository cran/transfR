## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=TRUE, message=FALSE, results='hide', eval=FALSE--------------------
#  library(transfR)
#  data(Oudon)
#  
#  st_write(st_sf(ID = paste0("ID",1:6), geom = st_geometry(Oudon$obs)),
#           dsn = "catchments.shp", delete_layer = T)
#  write.table(data.frame(DateTime = st_get_dimension_values(Oudon$obs,1),
#                         ID1 = Oudon$obs$Qobs[,1],
#                         ID2 = Oudon$obs$Qobs[,2],
#                         ID3 = Oudon$obs$Qobs[,3],
#                         ID4 = Oudon$obs$Qobs[,4],
#                         ID5 = Oudon$obs$Qobs[,5],
#                         ID6 = Oudon$obs$Qobs[,6] ),
#              file = "discharge.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## ---- echo=TRUE, message=FALSE, results='hide', eval=FALSE--------------------
#  library(sf)
#  catchments <- st_read("catchments.shp", "catchments", stringsAsFactors = F)
#  obs_sf <- catchments[1:5,] # Gauged catchments
#  sim_sf <- catchments[6,]   # Ungauged catchments

## ---- echo=TRUE, message=FALSE, results='hide', eval=FALSE--------------------
#  library(units)
#  Q <- read.table("discharge.txt", header = T, sep = "\t")
#  Qmatrix  <- as.matrix(Q[,-1])
#  Qmatrix  <- set_units(Qmatrix,"m^3/s")

## ---- echo=TRUE, message=FALSE, results='hide', eval=FALSE--------------------
#  library(stars)
#  Qmatrix  <- Qmatrix[,obs_sf$ID] #to have the same order as in the spacial data layer
#  obs_st   <- st_as_stars(list(Qobs = Qmatrix),
#                              dimensions = st_dimensions(time = as.POSIXct(Q$DateTime, tz="UTC"),
#                                                         space = obs_sf$geometry))
#  sim_st   <- st_as_stars(dimensions = st_dimensions(time = as.POSIXct(Q$DateTime, tz="UTC"),
#                                                         space = sim_sf$geometry))

## ---- echo=TRUE, message=FALSE, results='hide', eval=FALSE--------------------
#  obs <- as_transfr(st = obs_st, hl = Oudon$hl[1:5])
#  sim <- as_transfr(st = sim_st, hl = Oudon$hl[6])

## ---- echo=TRUE, message=FALSE, results='hide', eval=FALSE--------------------
#  sim <- quick_transfr(obs, sim)
#  sim$st

