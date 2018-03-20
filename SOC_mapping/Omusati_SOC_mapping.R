#########################################################################################################
## Mapping SOC concentrations and stocks for the Omusati Region, Nambia
## Version: 1.0, 20 March 2018
## Author: Dr. B. Kempen
## ISRIC - World Soil Information
#########################################################################################################

# CHANGE PATH to your LOCAL PATH where the folder with training materials is saved
setwd("D:/ISRIC/DSM_Namibia/Workingdir/Omusati")


###########################################################################
##### 1 Initialization

# load libraries
require(sp)
require(raster)
require(rgdal)
require(e1071)
require(caret)
require(randomForest)
require(leaflet)
require(htmlwidgets)
require(plotKML)

# set random seed
set.seed(2018)

# define projections
wgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
utm33s <- "+proj=utm +zone=33 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#aeqd <- "+proj=aeqd +lat_0=8.5 +lon_0=21.5 +x_0=5621452.01998 +y_0=5990638.42298 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"


###########################################################################
##### 2 Read and process soil point data

## SAMPLE DATA ------------------------------------------------------------------------
# read sampling data
d <- read.csv("Combined SOC - BD.csv", header = TRUE)

# compute SOC stock
d$ocs1 <- d$soc1/100 * d$bd1*1000 * 0.3


## STATISTICS ------------------------------------------------------------------------
# summary statistics; check for suspicious or unrealistic values
summary(d)

# skewness statistic; na.rm=T ignores NA values when computing the statistic
skewness(d$ocs1, na.rm=T);skewness(d$soc1, na.rm=T);skewness(d$soc2, na.rm=T);skewness(d$bd1, na.rm=T)

# plot histograms
hist(d$soc1, col = "lightblue", xlab="%", main="SOC, 0-30 cm")
hist(d$soc2, col = "lightblue", xlab="%", main="SOC, 30-100 cm")
hist(d$bd1, col = "lightblue", xlab="g/cm3", main="Bulk density, 0-30 cm")
hist(d$ocs1, col = "lightblue", xlab="kg/m2", main="SOC stock, 0-30 cm")

## PLOTTING ------------------------------------------------------------------------
# read boundary shapefile and set projection
admin <- readOGR(dsn = "./administrative", layer = "Omusati_admin_wgs84")
proj4string(admin) <- CRS(wgs84)

# prepare boundary file for plotting with spplot
admin.layout <- list(list("sp.polygons", admin))

# define color ramps for plotting
ramp.ocs1 <- c(0,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,4,5)
ramp.soc1 <- seq(0,1.2,0.2)

# spatially plot the SOC stock and SOC concentration (first convert to SpatialPointsDataFrame)
coordinates(d) <- ~lon+lat
spplot(d, zcol = "ocs1", cex = 1, main = "SOC stock(kg/m2) - Omusati", cuts = ramp.ocs1, key.space = "right", sp.layout = admin.layout, xlim = bbox(admin)[1, ], ylim = bbox(admin)[2, ], col.regions = bpy.colors(10), scales = list(draw=TRUE))
spplot(d, zcol = "soc1", cex = 1, main = "SOC (%) - Omusati", cuts = ramp.soc1, key.space = "right", sp.layout = admin.layout, xlim = bbox(admin)[1, ], ylim = bbox(admin)[2, ], col.regions = bpy.colors(10), scales = list(draw=TRUE))

# convert back to data.frame
d <- as.data.frame(d)


###########################################################################
##### 3 Reading covariate layers

## ------------------------------------------------------------------------
# read mask
m <- raster(paste0(getwd(),"/covariates/1km/mask.tif"))

# explore
str(m)
summary(m)

# plot
plot(m)


## ------------------------------------------------------------------------
# list file names of covariate layers
cov.lst <- list.files(path = "./covariates/1km/stack1", pattern =".tif")

# explore the set
head(cov.lst)
length(cov.lst)


## ------------------------------------------------------------------------
# create stack
r <- stack(paste0(getwd(), "/covariates/1km/stack1/", cov.lst)) 

# check object class
class(r)

# show dimensions
dim(r)

# check presence of projection
proj4string(r)

## ------------------------------------------------------------------------
# mask the covariate layers
r2 <- mask(r,m)


###########################################################################
##### 4 Processing covariate layers

## ------------------------------------------------------------------------
# convert RasterLayer with masked stack to data.frame (first convert to matrix)
p <- as(r2, Class = "SpatialGridDataFrame")
p <- as(p, Class = "data.frame")

# exclude pixels that do not have complete covariate data
p <- p[complete.cases(p),]

# list names of categorical covariates in stack (note this is Omusati specific, list can be different for other areas)
cat.list <- c("africasoil","hypsclass","riclass", "slopeclass","soterdomsoil","soterlandform","soterlitho","thermo")

# check data structure
str(p[,cat.list])

# convert to factor in a loop
for(i in 1:length(cat.list)){
  p[,cat.list[i]] <- as.factor(p[,cat.list[i]])
}

# check data structure; note the variables have been converted to factor
str(p[,cat.list])


## ------------------------------------------------------------------------
# convert categorical covariates to binary layers
for(i in 1:length(cat.list)){
  dum <- model.matrix(as.formula(paste0("~",cat.list[i],"+0")), p)
  dum <- as.data.frame(dum)
  dum <- lapply(dum, FUN = factor)
  dum <- as.data.frame(dum)
  p <- cbind(p, dum)
}

# exclude original categorical covariates
p <- p[,!names(p) %in% cat.list]

# clean-up
rm(dum, i)


## ------------------------------------------------------------------------
# check for zero or near zero variance covariates
nzv <- nearZeroVar(p, saveMetrics = TRUE)
nzv

# drop covariate layers that are zero or near zero variance based on the nzv column
p2 <- p[,!nzv[,4]]


## ------------------------------------------------------------------------
# statistical summary of covariate layers to check for suspicious/unrealistic values
summary(p2)


###########################################################################
##### 5 Creating the Regression Matrix

# convert data.frame with sampling sites to SpatialPointDataFrame
coordinates(d) <- ~lon+lat

# check if there is a projection
proj4string(d)

# project the data
proj4string(d) <- CRS(wgs84)
proj4string(d)

# identify sampling sites with identical coordinates
zerodist(d, zero = 0.0, unique.ID = FALSE) # no sampling sites with identical coordinates

# remove duplicated sampling sites
d <- remove.duplicates(d, zero = 0.0, remove.second = TRUE)

# reproject to azimuthal 
d <- spTransform(d, CRS(utm33s))

# convert covariate stack to SpatialPointsDataFrame object
coordinates(p2) <- ~s1+s2

# define projection of the covariate stack
proj4string(p2) <- CRS(utm33s)
proj4string(p2) # check projection

# convert covariate SpatialPointDataFrame to SpatialPixelsDataFrame
gridded(p2) <- TRUE
#fullgrid(p2) <- TRUE # convert SpatialPixelsDataFrame to SpatialGridDataFrame

# overlay points with covariate stack; output is a data.frame
dum <- over(x = d, y = p2)

# convert SpatialPointDF with sampling sites to data.frame object
d <- as(d, Class="data.frame")

# append the covariate data to the sampling sites
d <- cbind(d,dum)

# clean-up
rm(dum)


###########################################################################
##### 6 Calibrating a prediction model

# copy regression matrix to a new object
rm <- d

# exclude sampling sites with missing OCS/SOC data
rm <- rm[!is.na(rm$ocs1),]
#rm <- rm[!is.na(rm$soc1),]

# summary statsistics: check categorical covariates
summary(rm)

# drop categorical covariate that have zero variation
nzv3 <- nearZeroVar(rm, saveMetrics = TRUE)
rm <- rm[,!nzv3[,3]]

# identify the columns in which the covariates are stored
names(rm) # here: columns 10 to 104

# remove observations without covariate data
rm <- rm[complete.cases(rm[,10:ncol(rm)]),]

# generate random forest input objects 1: vector with carbon data
z <- rm$ocs1 # or: z <- rm[,"ocs1"]
#z <- rm$soc1 # or: z <- rm[,"soc1"]

# generate random forest input objects 2: data.frame with covariate data
covar <- rm[,10:ncol(rm)]

# fit a random forest model
rf <- randomForest(x = covar, y = z, ntree = 500)

# explore model results
str(rf)

# variable importance plot
varImpPlot(rf)

# r-square
rf$rsq[500] # number should be equal to the total number of trees, or use rf$ntree to be more generic

# correlation plot
plot(rf$predicted,rf$y)
abline(0,1)


###########################################################################
##### 7 Mapping

# convert covariate stack from SpatialPixel/GridDF to data.frame
p2 <- as.data.frame(p2)
#p2 <- as(p2, Class = "data.frame")

# apply model to predict SOC/OCS at each pixel
pred <- predict(rf, newdata = p2)

# combine predictions with coordinates
p.map <- data.frame(
  lon = p2$s1,
  lat = p2$s2,
  pred = pred
)

# convert data.frame with predictions to a SpatialGridDataFrame
coordinates(p.map) <- ~lon+lat
gridded(p.map) <- TRUE
fullgrid(p.map) <- TRUE
proj4string(p.map) <- CRS(utm33s)

# plot
spplot(p.map, zcol="pred", main = "Organic Carbon Stock (kg/m2)")

# save plot
png("OCS1_Omusatie.png", width = 500, height = 700)
spplot(p.map, zcol="pred", main = "Organic Carbon Stock (kg/m2)")
dev.off()

# export map to GeoTiff
writeGDAL(p.map["pred"], fname = "OCS1_Omusati.tif", drivername = "GTiff", type = "Float32")


###########################################################################
##### 8 Plotting

## KML ------------------------------------------------------------------------
# write maps to KML (google earth)
kml(p.map["pred"], file.name = "OCS1_Omusati.kml", colour_scale = SAGA_pal[[1]])
#kml(p.map["pred"], file.name = "SOC1_Omusati.kml", colour_scale = SAGA_pal[[1]])

# write sample data to KML (google earth)
coordinates(d) <- ~lon+lat
proj4string(d) <- CRS(utm33s)
shape <- "http://maps.google.com/mapfiles/kml/pal2/icon18.png"
kml(d, file.name = "OCS1_Omusati_sample.kml", colour = d@data[,"ocs1"], balloon = TRUE, shape = shape,
    colour_scale = SAGA_pal[[1]], points_names = round(d$ocs1, 1))


## LEAFLET ------------------------------------------------------------------------
# convert the map to RasterLayer class
p.map.r <- raster(p.map)

# define a color palette
pal <- colorNumeric(SAGA_pal[[1]], values(p.map.r), na.color = "transparent")

# generate leaflet map
(ll <- leaflet() %>% 
  addTiles() %>%
  addProviderTiles("OpenStreetMap") %>%
  addRasterImage(p.map.r, colors=pal, opacity=0.5) %>%
  addLegend(pal = pal, values=values(p.map.r), title="OCS (kg/m2) (0-30 cm)")
)

# save leaflet as an htmlWidget
saveWidget(ll, file = paste0("OCS1_Omusati.html",sep=""))

# end of script;