#########################################################################################################
## Downloading TIF files from an FTP server and reprojection
## Version: 1.0, 4 April2018
## Author: Dr. B. Kempen
## ISRIC - World Soil Information
#########################################################################################################

###########################################################################
##### 1 Initialization

# load libraries
library(gdalUtils)
library(RCurl)

# set working directory
setwd("D:/_temp")

# define projections
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
utm33s <- "+proj=utm +zone=33 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# create directory in the working directory to store the covariate layers
dir.create(paste0(getwd(), "/covariates"))


###########################################################################
##### 2 Define FTP location and retrieve file names

# define URL
url <- paste0("ftp://isric.org/Namibia/",  "covs", "/")
url

# define username and password
userpwd <- "gsp:gspisric"

# retrieve filenames
filenames <- getURL(url, userpwd = userpwd, ftp.use.epsv = FALSE, dirlistonly = TRUE) 
filenames

# clean-up filenames
filenames <- strsplit(filenames, "\r*\n")[[1]]
filenames

# select tif files
sel <- grep(pattern=".tif$" , filenames)
filenames <- filenames[sel]
filenames


###########################################################################
##### 3 Download covariate layers

# download covariates
for(j in 1:length(filenames)){
  # define input file
  inname <- paste0(getwd(), "/covariates/",  filenames[j])
  
  # download
  download.file(paste0("ftp://gsp:gspisric@isric.org/Namibia/",  "covs/", filenames[j]), destfile = inname, method="auto",mode="wb")
}


###########################################################################
##### 4 Reproject (here from WGS84 to UTM33s)

# create a directory to store reprojected rasters
dir.create(paste0(getwd(), "/covariates/reproject"))

### example how to reproject one layer
# define input and output files
inname <- paste0(getwd(), "/covariates/",  filenames[1])
inname

outname <- paste0(getwd(), "/covariates/reproject/",  filenames[1])
outname

# reproject (note this takes might while; be patient)
gdalwarp(srcfile = inname, dstfile = outname,  s_srs=wgs84, t_srs=utm33s , r="near", tr=c(1000, 1000))

### reproject all layers in the stack: TAKES TIME!!!!
for(j in 1:length(filenames)){
  # define input file
  inname <- paste0(getwd(), "/covariates/",  filenames[j])
  
  # define output file
  outname <- paste0(getwd(), "/covariates/reproject/",  filenames[j])
  
  # reproject
  gdalwarp(srcfile = inname, dstfile = outname,  s_srs=wgs84, t_srs=aeqd , r="near", tr=c(1000, 1000))
}

# end of script;