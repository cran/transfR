## ----setup, include = FALSE---------------------------------------------------
library(transfR)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=TRUE, message=FALSE, warning=FALSE, results='hide'-----------------
wbt_wd <- tempdir(check = TRUE)

## ---- download_dem, echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE, results='hide'----
#  library(elevatr)
#  library(rgdal) # Still often needed by elevatr
#  library(progress) # Still often needed by elevatr
#  
#  # Set up a projection (French Lambert-93 projection)
#  EPSG <- 2154
#  
#  # Define a bbox that will encompass the catchments of the study area
#  blavet_bbox <- st_bbox(c(xmin = -3.3, xmax = -2.7, ymax = 48.11, ymin = 47.77),
#                             crs = st_crs(4326))
#  
#  # Retrieve elevation data as raster
#  dem_raw <- elevatr::get_elev_raster(st_as_sfc(blavet_bbox), z = 10) # ~76m resolution
#  
#  # Project and define spatial resolution:
#  dem_100m  <- st_warp(st_as_stars(dem_raw), cellsize = 100, crs = st_crs(EPSG))
#  names(dem_100m) <- "warp"
#  
#  # Set negative values (ocean) to NA
#  dem_100m[dem_100m < 0] <- NA
#  
#  # Write to file
#  write_stars(dem_100m["warp"], file.path(wbt_wd,"dem_100m.tif"))

## ---- echo=FALSE, message=FALSE, warning=TRUE, eval=TRUE, results='hide'------
try_chunk <- try({
library(elevatr)
library(rgdal) # Still often needed by elevatr
library(progress) # Still often needed by elevatr

# Set up a projection (French Lambert-93 projection)
EPSG <- 2154 

# Define a bbox that will encompass the catchments of the study area
blavet_bbox <- st_bbox(c(xmin = -3.3, xmax = -2.7, ymax = 48.11, ymin = 47.77), 
                           crs = st_crs(4326))

# Retrieve elevation data as raster
dem_raw <- elevatr::get_elev_raster(st_as_sfc(blavet_bbox), z = 10) # ~76m resolution

# Project and define spatial resolution: 
dem_100m  <- st_warp(st_as_stars(dem_raw), cellsize = 100, crs = st_crs(EPSG))
names(dem_100m) <- "warp"

# Set negative values (ocean) to NA
dem_100m[dem_100m < 0] <- NA

# Write to file
write_stars(dem_100m["warp"], file.path(wbt_wd,"dem_100m.tif"))
}, silent = TRUE)
if(inherits(try_chunk, "try-error")){
  warning("\nIssue when downloading elevation data. \nThe vignette will not be fully built.")
  running <- FALSE
}else{
  running <- TRUE
}

## ---- install_whitebox, echo=FALSE, message=FALSE, warning=TRUE, eval=running, results='hide'----
# If WhiteboxTools executable are not present, install it in the temporary directory
library(whitebox)
if(!wbt_init()){
  wbt_inst <- try(install_whitebox(pkg_dir = wbt_wd), silent = TRUE)
  if(!inherits(wbt_inst, "try-error")){
    exe_path <- file.path(wbt_wd, "WBT", "whitebox_tools") # Unix
    if(!file.exists(exe_path)) exe_path <- paste0(exe_path,".exe") # Windows
    if(!file.exists(exe_path)){
      warning("WhiteboxTools executable not found")
      running <- FALSE
    }else{
        wbt_options(exe_path = exe_path,
                max_procs = 2,
                wd = wbt_wd,
                verbose = TRUE)
      }
    }else{running <- FALSE}
  }else{
    wbt_options(max_procs = 2,
                wd = wbt_wd,
                verbose = TRUE)
  }

## ---- burn_stream, echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE, results='hide'----
#  library(transfR)
#  library(whitebox)
#  
#  # Get the French Topage river network from the Blavet dataset
#  data(Blavet)
#  CoursEau_Topage2019 <- Blavet$network
#  
#  # Change projection and write files
#  network_topage <- st_transform(CoursEau_Topage2019, EPSG)
#  st_write(network_topage, file.path(wbt_wd, "network_topage.shp"),
#               delete_layer = TRUE, quiet = TRUE)
#  whitebox::wbt_rasterize_streams("network_topage.shp",
#                                  base = "dem_100m.tif",
#                                  output = "network_topage.tif",
#                                  nodata = 0,
#                                  wd = wbt_wd)
#  
#  # Burn this river network on the DEM
#  # We will neglect the effect of the road embankments at this DEM resolution of 100m
#  # by creating an empty shapefile for roads
#  st_write(st_sfc(st_multilinestring(),crs = EPSG), file.path(wbt_wd,"roads.shp"),
#               delete_layer = TRUE, quiet = TRUE)
#  whitebox::wbt_burn_streams_at_roads(dem = "dem_100m.tif",
#                          streams = "network_topage.shp",
#                          roads = "roads.shp",
#                          output = "dem_100m_burn.tif",
#                          wd = wbt_wd)

## ---- echo=FALSE, message=FALSE, warning=TRUE, eval=running, results='hide'----
try_chunk <- try({
library(transfR)
library(whitebox)

# Get the French Topage river network from the Blavet dataset
data(Blavet)
CoursEau_Topage2019 <- Blavet$network

# Change projection and write files
network_topage <- st_transform(CoursEau_Topage2019, EPSG)
st_write(network_topage, file.path(wbt_wd, "network_topage.shp"), 
             delete_layer = TRUE, quiet = TRUE)
whitebox::wbt_rasterize_streams("network_topage.shp", 
                                base = "dem_100m.tif", 
                                output = "network_topage.tif", 
                                nodata = 0, 
                                wd = wbt_wd)

# Burn this river network on the DEM
# We will neglect the effect of the road embankments at this DEM resolution of 100m 
# by creating an empty shapefile for roads
st_write(st_sfc(st_multilinestring(),crs = EPSG), file.path(wbt_wd,"roads.shp"), 
             delete_layer = TRUE, quiet = TRUE)
whitebox::wbt_burn_streams_at_roads(dem = "dem_100m.tif", 
                        streams = "network_topage.shp", 
                        roads = "roads.shp",
                        output = "dem_100m_burn.tif", 
                        wd = wbt_wd)
}, silent = TRUE)
if(inherits(try_chunk, "try-error")){
  warning("\nIssue when burning the river network on the DEM. \nThe vignette will not be fully built.")
  running <- FALSE
}

## ---- extract_stream, echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE, results='hide'----
#  # Remove the depressions on the DEM
#  whitebox::wbt_fill_depressions(dem = "dem_100m_burn.tif",
#                                 output = "dem_fill.tif",
#                                 wd = wbt_wd)
#  
#  # Flow direction raster
#  whitebox::wbt_d8_pointer(dem = "dem_fill.tif",
#                           output = "d8.tif",
#                           wd = wbt_wd)
#  
#  # Compute flow accumulation
#  whitebox::wbt_d8_flow_accumulation(input = "d8.tif",
#                                     pntr = TRUE,
#                                     output ="facc.tif",
#                                     wd = wbt_wd)
#  
#  # Extract a stream network (threshold = 1 km2) consistent with flow direction
#  whitebox::wbt_extract_streams(flow_accum =  "facc.tif",
#                                threshold = 100, # 100 cells for 1 km2
#                                output = "network_1km2.tif",
#                                zero_background = TRUE,
#                                wd = wbt_wd)
#  whitebox::wbt_remove_short_streams(d8_pntr = "d8.tif",
#                                     streams = "network_1km2.tif",
#                                     output = "network_d8.tif",
#                                     min_length= 200,
#                                     wd = wbt_wd)

## ---- echo=FALSE, message=FALSE, warning=TRUE, eval=running, results='hide'----
try_chunk <- try({
# Remove the depressions on the DEM
whitebox::wbt_fill_depressions(dem = "dem_100m_burn.tif",
                               output = "dem_fill.tif",
                               wd = wbt_wd)

# Flow direction raster
whitebox::wbt_d8_pointer(dem = "dem_fill.tif",
                         output = "d8.tif",
                         wd = wbt_wd)

# Compute flow accumulation
whitebox::wbt_d8_flow_accumulation(input = "d8.tif",
                                   pntr = TRUE,
                                   output ="facc.tif",
                                   wd = wbt_wd)

# Extract a stream network (threshold = 1 km2) consistent with flow direction
whitebox::wbt_extract_streams(flow_accum =  "facc.tif",
                              threshold = 100, # 100 cells for 1 km2
                              output = "network_1km2.tif",
                              zero_background = TRUE,
                              wd = wbt_wd)
whitebox::wbt_remove_short_streams(d8_pntr = "d8.tif",
                                   streams = "network_1km2.tif",
                                   output = "network_d8.tif",
                                   min_length= 200,
                                   wd = wbt_wd)
}, silent = TRUE)
if(inherits(try_chunk, "try-error")){
  warning("\nIssue when extracting the river network from the DEM. \nThe vignette will not be fully built.")
  running <- FALSE
}

## ---- delineate_catchments, echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE, results='hide'----
#  # Localize the outlets of the studied catchments (with manual adjustments to help snapping)
#  outlets_coordinates <- data.frame(id = names(Blavet$hl),
#                               X = c(254010.612,255940-100,255903,237201,273672,265550),
#                               Y = c(6772515.474,6776418-200,6776495,6774304-200,6762681,6783313))
#  outlets <- st_as_sf(outlets_coordinates, coords = c("X", "Y"), crs=2154)
#  st_write(outlets, dsn = file.path(wbt_wd, "outlets.shp"),
#               delete_layer = TRUE, quiet = TRUE)
#  
#  # Snap the outlets on the stream raster
#  whitebox::wbt_jenson_snap_pour_points(pour_pts = "outlets.shp",
#                                        streams = "network_d8.tif",
#                                        output = "outlets_snapped.shp",
#                                        snap_dist = 200,
#                                        wd = wbt_wd)
#  outlets_snapped <- st_read(file.path(wbt_wd, "outlets_snapped.shp"), quiet = TRUE)
#  
#  # Delineate catchments
#  catchments <- st_sfc(crs = EPSG)
#  for(id in outlets_snapped$id){
#    st_write(outlets_snapped[outlets_snapped$id==id,],
#                 file.path(wbt_wd, paste0(id, "_outlet.shp")), delete_layer = TRUE, quiet = TRUE)
#    whitebox::wbt_watershed(d8_pntr = "d8.tif",
#                            pour_pts = paste0(id, "_outlet.shp"),
#                            output = paste0(id, "_catchment.tif"),
#                            wd = wbt_wd)
#    # Vectorize catchments
#    drainage <- read_stars(file.path(wbt_wd, paste0(id, "_catchment.tif")))
#    contours <- st_contour(drainage, breaks = 1) |> st_geometry() |> st_cast("POLYGON")
#    contours <- contours[which.max(st_area(contours))]
#    catchments <- rbind(catchments, st_sf(data.frame(id, geom = contours)))
#  }

## ---- echo=FALSE, message=FALSE, warning=TRUE, eval=running, results='hide'----
try_chunk <- try({
# Localize the outlets of the studied catchments (with manual adjustments to help snapping)
outlets_coordinates <- data.frame(id = names(Blavet$hl),
                             X = c(254010.612,255940-100,255903,237201,273672,265550),
                             Y = c(6772515.474,6776418-200,6776495,6774304-200,6762681,6783313))
outlets <- st_as_sf(outlets_coordinates, coords = c("X", "Y"), crs=2154)
st_write(outlets, dsn = file.path(wbt_wd, "outlets.shp"), 
             delete_layer = TRUE, quiet = TRUE)

# Snap the outlets on the stream raster
whitebox::wbt_jenson_snap_pour_points(pour_pts = "outlets.shp",
                                      streams = "network_d8.tif",
                                      output = "outlets_snapped.shp",
                                      snap_dist = 200,
                                      wd = wbt_wd) 
outlets_snapped <- st_read(file.path(wbt_wd, "outlets_snapped.shp"), quiet = TRUE)

# Delineate catchments 
catchments <- st_sfc(crs = EPSG)
for(id in outlets_snapped$id){
  st_write(outlets_snapped[outlets_snapped$id==id,], 
               file.path(wbt_wd, paste0(id, "_outlet.shp")), delete_layer = TRUE, quiet = TRUE)
  whitebox::wbt_watershed(d8_pntr = "d8.tif",
                          pour_pts = paste0(id, "_outlet.shp"),
                          output = paste0(id, "_catchment.tif"),
                          wd = wbt_wd)
  # Vectorize catchments
  drainage <- read_stars(file.path(wbt_wd, paste0(id, "_catchment.tif")))
  contours <- st_contour(drainage, breaks = 1) |> st_geometry() |> st_cast("POLYGON")
  contours <- contours[which.max(st_area(contours))]
  catchments <- rbind(catchments, st_sf(data.frame(id, geom = contours)))
}
}, silent = TRUE)
if(inherits(try_chunk, "try-error")){
  warning("\nIssue when delineating catchments. \nThe vignette will not be fully built.")
  running <- FALSE
}

## ---- echo=TRUE, message=FALSE, warning=FALSE, eval=running-------------------
# Compare drainage areas to the dataset provided with transfR
compare_areas <- data.frame(
  name = names(Blavet$hl),
  expected_area =  st_area(st_geometry(Blavet$obs)) |> units::set_units("km^2") |> round(1),
  computed_area = st_area(catchments) |> units::set_units("km^2") |> round(1)
)

print(compare_areas)

## ---- echo=TRUE, message=FALSE, warning=FALSE, eval=running, results='hide', fig.width=7, fig.height=4----
# Plot catchment delineation
par(oma = c(0, 0, 0, 6))
plot(catchments[,"id"], main = "Catchments delineation", key.pos = 4, reset = FALSE)
plot(st_geometry(st_intersection(network_topage,catchments)), 
     col = "white", lwd = 1.5, add = TRUE)
plot(outlets_snapped, col = "black", pch = 16, add = TRUE)

## ---- hydraulic_length, echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE, results='hide'----
#  # Compute hydraulic length
#  whitebox::wbt_downslope_flowpath_length(d8_pntr = "d8.tif",
#                                          output = "fpl.tif",
#                                          wd = wbt_wd)
#  whitebox::wbt_downslope_distance_to_stream(dem = "dem_fill.tif",
#                                             streams = "network_topage.tif",
#                                             output = "d2s.tif",
#                                             wd = wbt_wd)
#  fpl <- read_stars(file.path(wbt_wd, "fpl.tif"))
#  d2s <- read_stars(file.path(wbt_wd, "d2s.tif"))
#  hl_region <- fpl-d2s
#  names(hl_region) <- "hl"
#  
#  # Crop hydraulic length for each catchment
#  hl <- list()
#  for(id in catchments$id){
#    crop <- st_crop(hl_region, catchments[catchments$id==id,])
#    crop <- crop-min(crop$hl, na.rm = TRUE)
#    crop$hl <- units::set_units(crop$hl, "m")
#    hl[[id]] <- crop
#  }
#  

## ---- echo=FALSE, message=FALSE, warning=TRUE, eval=running, results='hide'----
try_chunk <- try({
# Compute hydraulic length
whitebox::wbt_downslope_flowpath_length(d8_pntr = "d8.tif",
                                        output = "fpl.tif",
                                        wd = wbt_wd)
whitebox::wbt_downslope_distance_to_stream(dem = "dem_fill.tif",
                                           streams = "network_topage.tif",
                                           output = "d2s.tif",
                                           wd = wbt_wd)
fpl <- read_stars(file.path(wbt_wd, "fpl.tif"))
d2s <- read_stars(file.path(wbt_wd, "d2s.tif"))
hl_region <- fpl-d2s
names(hl_region) <- "hl"

# Crop hydraulic length for each catchment
hl <- list()
for(id in catchments$id){
  crop <- st_crop(hl_region, catchments[catchments$id==id,])
  crop <- crop-min(crop$hl, na.rm = TRUE)
  crop$hl <- units::set_units(crop$hl, "m")
  hl[[id]] <- crop
}

}, silent = TRUE)
if(inherits(try_chunk, "try-error")){
  warning("\nIssue when computing hydraulic length. \nThe vignette will not be fully built.")
  running <- FALSE
}

## ---- echo=TRUE, message=FALSE, warning=FALSE, eval=running, results='hide', fig.width=7, fig.height=4----
i <- 1
network <- st_geometry(st_intersection(network_topage,catchments[i,]))
plot(hl[[i]], main = paste("Hydraulic length of catchment", i,"[m]"), 
     col = hcl.colors(20, palette = "Teal"), key.pos = 1, reset = FALSE)
plot(network, col = "white", lwd = 1.5, add = TRUE)

## ---- echo=TRUE, message=FALSE, warning=FALSE, eval=running, results='hide'----
obs_st <- st_as_stars(list(Qobs = Blavet$obs$Qobs), 
                            dimensions = st_dimensions(
                              time = st_get_dimension_values(Blavet$obs,1),
                              space = st_geometry(catchments)))

## ---- echo=TRUE, message=FALSE, warning=FALSE, eval=running, results='hide'----
obs <- as_transfr(st = obs_st, hl = hl)

## ---- echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE---------------------
#  obs <- quick_transfr(obs, velocity = "brittany2013", parallel = TRUE, cores = 2, cv = TRUE)

## ---- echo=TRUE, message=FALSE, warning=FALSE, eval=FALSE---------------------
#  obs$st

## ---- echo=FALSE, message=FALSE, warning=FALSE, eval=TRUE, results='hide'-----
# Cleaning temporary directory
unlink(wbt_wd, recursive = TRUE)

