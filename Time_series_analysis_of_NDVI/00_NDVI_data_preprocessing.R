# NDVI data preprocessing
# Author: John Mutua
# International Center for Tropical Agriculture (CIAT)
# Last modified: March 16 2018

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
# set variables
iDir <- "E:/Omusati_covariates/ndvi/data"
oDir <- "E:/Omusati_covariates/ndvi"
aoi <- "Omusati"

# define projections
utm33s <- CRS('+proj=utm +zone=33 +south +datum=WGS84 +units=m +no_defs 
              +ellps=WGS84 +towgs84=0,0,0')

# read in polygon file for use in clipping
aoi = shapefile(paste0(iDir, "/", aoi, ".shp"))
              

## ----message=FALSE, warning=FALSE----------------------------------------
# create output folder
if (!dir.exists(paste0(oDir, "/", "processed_data")))
  dir.create(paste0(oDir, "/", "processed_data"))

## ------------------------------------------------------------------------
# read, list and stack files
r.list <- list.files(paste0(iDir, "/"), pattern = ".tif$", full.names = T)
raw.stack <- stack(r.list)

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
writeRaster(r.stack, filename = paste0(oDir, "/", "processed_data", "/"), 
            names(raw.stack), bylayer=TRUE, format = "GTiff", overwrite=TRUE)

