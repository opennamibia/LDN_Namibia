## ----echo=FALSE----------------------------------------------------------
# NDVI data preprocessing
# Author: John M
# International Center for Tropical Agriculture (CIAT)
# Last modified: March 1 2018


setwd("C:/Users/computation/Desktop/Final_NDVI_27_Mar_18")
## ------------------------------------------------------------------------
# clear your work space
rm(list = ls(all = TRUE))

## ---- message=FALSE, warning=FALSE---------------------------------------
# load packages
library(rgdal)
library(raster)
library(snow)

## ---- eval=FALSE---------------------------------------------------------
## # load packages
## .packages = c("rgdal","raster","snow")
## .inst <- .packages %in% installed.packages()
## if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
## lapply(.packages, require, character.only=TRUE)

## ---- message=FALSE, warning=FALSE---------------------------------------
# set variables USE THIS INSTEAD OF SETTING THE WORKING DIRECTORY
iDir <- "C:/Users/computation/Desktop/Final_NDVI_27_Mar_18/data"
oDir <- "C:/Users/computation/Desktop/Final_NDVI_27_Mar_18"
#aoi <- "Omusati" #this would be use in case the .shp was named Omusati

# define projections
utm33s <- CRS('+proj=utm +zone=33 +south +datum=WGS84 +units=m +no_defs 
              +ellps=WGS84 +towgs84=0,0,0')

# read in polygon file for use in clipping
aoi = shapefile(paste0(iDir, "/", "Omus_33s" ,".shp"))




## ----message=FALSE, warning=FALSE----------------------------------------
# create output folder
if (!dir.exists(paste0(oDir, "/", "processed_data")))
  dir.create(paste0(oDir, "/", "processed_data"))

## ------------------------------------------------------------------------
# read, list and stack files
r.list <- list.files(paste0(iDir, "/"), pattern = ".tif$", full.names = T)
raw.stack <- stack(r.list)

#raw.stack <- stack(r.list[1:87]) #if we are to slipt it

## ------------------------------------------------------------------------
# project stack to utm 33s
beginCluster()
r.stack = projectRaster(raw.stack, crs = utm33s)

## ------------------------------------------------------------------------
# mask using polygon
r.stack <- mask(r.stack, aoi)

## ------------------------------------------------------------------------
# create function to calculate float values
fun <- function(x) {
  x = (x - 50) / 200
  return(x)
}

# apply function on the stack
r.stack<-calc(r.stack, fun=fun)
endCluster()

## ------------------------------------------------------------------------
# save as GeoTIFF
writeRaster(r.stack, filename = paste0(oDir, "/", "processed_data", "/"), names(raw.stack), bylayer=TRUE, format = "GTiff", overwrite=TRUE)

