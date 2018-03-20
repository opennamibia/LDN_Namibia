# Bush density mapping
# Author: John M
# International Center for Tropical Agriculture (CIAT)
# Last modified: March 1 2017

## ------------------------------------------------------------------------
# clear your work space
rm(list = ls(all = TRUE))

## ------------------------------------------------------------------------
# set the random seed
set.seed(123)

## ----eval=TRUE-----------------------------------------------------------
# set variables
iDir <- "D:/OneDrive - CGIAR/ToBackup/_GitHub/LDN_Namibia/Bush_density_mapping/data"
oDir <- "D:/OneDrive - CGIAR/ToBackup/_GitHub/LDN_Namibia/Bush_density_mapping"
aoi <- "Otjiwarongo"

## ---- message=FALSE, warning=FALSE---------------------------------------
# load packages
library(sp)
library(rgdal)
library(raster)
library(randomForest)
library(plyr)
library(dplyr)
library(caret)
library(dygraphs)
library(car)
library(e1071)
library(snow)
library(plotKML)

## ---- eval=FALSE---------------------------------------------------------
## # load packages
## .packages = c("sp","rgdal","raster","randomForest","plyr",
##               "dplyr","caret","car", "e1071","snow","plotKML")
## .inst <- .packages %in% installed.packages()
## if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
## lapply(.packages, require, character.only=TRUE)

## ------------------------------------------------------------------------
# read in data
raw.d <- read.csv(paste0(iDir, "/", aoi, "_bd_sampling_points", ".csv"), header=TRUE)

## ------------------------------------------------------------------------
# calculate values
raw.d$shrubs_less_1.5 <- apply(raw.d[,8:11], 1, sum, na.rm=TRUE)
raw.d$shrubs_more_1.5_no_stem <- apply(raw.d[,16:19], 1, sum, na.rm=TRUE)
raw.d$shrubs_more_1.5_stem <- apply(raw.d[,24:27], 1, sum, na.rm=TRUE)

## ------------------------------------------------------------------------
# create new dataframe with columns you need
raw.d<-raw.d[,c(1,2,3,34,35,36)]

## ------------------------------------------------------------------------
# add two new columns of shrubs > 1.5 and all shrubs in general
raw.d$shrubs_more_1.5 <- apply(raw.d[,5:6], 1, sum, na.rm=TRUE)

## ------------------------------------------------------------------------
# select the columns you need
raw.d<-raw.d[,c(1,2,3,4,7)]

## ------------------------------------------------------------------------
# compute shrubs per hectare
raw.d$shrubs_less_1.5 <- raw.d$shrubs_less_1.5*25
raw.d$shrubs_more_1.5 <- raw.d$shrubs_more_1.5*25

## ------------------------------------------------------------------------
# remove all NAs
raw.d<-raw.d[complete.cases(raw.d),]

## ------------------------------------------------------------------------
# plot histograms of the three variables
hist(raw.d$shrubs_less_1.5, col = "lightblue", xlab="Count", main="Shrubs, [<1.5]")
hist(raw.d$shrubs_more_1.5, col = "lightblue", xlab="Count", main="Shrubs, [>1.5]")

## ---- message=FALSE, warning=FALSE---------------------------------------
# create output folder
if (!dir.exists(paste0(oDir, "/", "outputs")))
  dir.create(paste0(oDir, "/", "outputs"))

## ------------------------------------------------------------------------
# export data to .csv
write.csv(raw.d, file = paste0(oDir, "/", "outputs", "/", aoi, "_OBS_Data", ".csv"),
          row.names=FALSE)

## ------------------------------------------------------------------------
# make shapefiles
xy <- raw.d[,c(3,2)]
trainDatageo <- SpatialPointsDataFrame(coords = xy, data = raw.d,
                                    proj4string = CRS("+proj=longlat 
                                                      +datum=WGS84"))
trainData <- spTransform(trainDatageo, CRS('+proj=utm +zone=33 +south 
                                           +datum=WGS84'))

## ------------------------------------------------------------------------
# check column names
names(trainData)

## ------------------------------------------------------------------------
# import the rest of input data and stack
r.list<-list.files(paste0(iDir, "/"), pattern = ".tif$", full.names = TRUE)
r.stack <- stack(r.list)
names(r.stack) <- c("asp","b1","b15","b16","b17","b18","cc","elev",
                    "l8b3","ndvi","soil")

## ------------------------------------------------------------------------
# import the bush area mask
o.mask <- raster(paste0(iDir, "/", "other_data", "/", "Otji_BushArea_2016", ".tif"))

## ------------------------------------------------------------------------
# set extent of the training data
trainData@bbox <- bbox(o.mask)

## ------------------------------------------------------------------------
# mask and remove NAs in the covariates
covs <- mask(r.stack, o.mask )
covs <- na.omit(covs)

## ------------------------------------------------------------------------
# assign raster values to the training data
v<-as.data.frame(extract(covs,trainData))
trainData@data=data.frame(trainData@data, v[match(rownames(trainData@data),
                                                  rownames(v)),])

## ------------------------------------------------------------------------
# rename fields in the training dataset, remove NAs, write the dataset 
names(trainData) <- c("waypoint.no","lat","lon","shrubs.l1.5","shrubs.g1.5",
                      "asp","b1","b15","b16","b17","b18","cc", "elev","l8b3",
                      "ndvi","soil")
trainData@data<-trainData@data[complete.cases(trainData@data),]
write.csv(trainData@data, file = paste0(oDir, "/", "outputs", "/", aoi, "_trainingData.csv"), row.names=FALSE)

## ------------------------------------------------------------------------
# convert trainData object to data.frame
d <- as.data.frame(trainData@data, na.rm=TRUE)

# skewness statistics for shrubs <1.5m
skewness(d$shrubs.l1.5, na.rm = T)

# skewness statistics for shrubs >1.5m
skewness(d$shrubs.g1.5, na.rm = T)

## ---- eval=FALSE---------------------------------------------------------
## # compute correlation coefficients and plot correlations
## cor(d[,4:16])
## pairs(d[,4:16])

## ---- eval=FALSE---------------------------------------------------------
## # correlate count of shrubs with NDVI and Landsat 8 band 2-7
## cor(d$shrubs.l1.5,d$cc)
## cor(d$shrubs.g1.5,d$cc)

## ---- eval=FALSE---------------------------------------------------------
## # histogram before transformation
## hist(d$shrubs.l1.5)
## plot(shrubs.l1.5~b1,data=d)
## 
## # histogram after transformation
## hist(log(d$shrubs.l1.5))
## plot(log(shrubs.l1.5)~log(b1),data=d)
## 
## # regression model for log tranformed data for shrubs<1.5
## reg1 <- lm(log(shrubs.l1.5)~log(b1),data=d)
## summary(reg1)
## 
## # regression model for log tranformed data for shrubs>1.5
## reg2 <- lm(log(shrubs.g1.5)~log(b1),data=d)
## summary(reg2)

## ------------------------------------------------------------------------
# log tranformed lm1
x1 <- log(shrubs.l1.5)~b1+b15+cc+ndvi+soil
lm1 <- lm(x1,data=d)

# log tranformed lm2
x2 <- log(shrubs.g1.5)~b1+b15+cc+ndvi+soil
lm2 <- lm(x2,data=d)

## ---- eval=FALSE---------------------------------------------------------
## # model 1 summary
## summary(lm1)
## 
## # model 2 summary
## summary(lm2)

## ------------------------------------------------------------------------
# select covariate values
s.cov <- d[,c(6:16)]

# select response variable [shrubs.l1.5]
s.l1.5<-d[,4]

# select response variable [shrubs.g1.5]
s.g1.5<-d[,5]

## ----message=FALSE, warning=FALSE----------------------------------------
# fit random forest model for shrubs <1.5m
model1 <- randomForest(x = s.cov, y = s.l1.5, ntree=5000, 
                       mtry=5, importance=T)

# fit random forest model for shrubs >1.5m
model2 <- randomForest(x = s.cov, y = s.g1.5, ntree=5000, 
                       mtry=5, importance=T)

# predicted versus observed
cor(model1$predicted,model1$y)**2
cor(model2$predicted,model2$y)**2

## ---- eval=FALSE---------------------------------------------------------
## # get trees
## getTree(model1,k=2)
## getTree(model2,k=2)
## 
## # explore the models
## str(model1, max.level=2)
## str(model2, max.level=2)

## ------------------------------------------------------------------------
# plot model 1
plot(model1, main = "Error rate across decision trees [model1]")

# plot model 2
plot(model2, main = "Error rate across decision trees [model2]")

## ------------------------------------------------------------------------
# variable importance plot for model 1
varImpPlot(model1, type = 2, main="Variable Importance [Model 1]")

# variable importance plot for model 2
varImpPlot(model2, type = 2, main="Variable Importance [Model 2]")

## ------------------------------------------------------------------------
# predict models
beginCluster()
prediction1 <- clusterR(covs, raster::predict, args = list(model = model1))
prediction2 <- clusterR(covs, raster::predict, args = list(model = model2))
endCluster()

## ------------------------------------------------------------------------
# round the numbers
prediction1<-round(prediction1, digits = 0)
prediction2<-round(prediction2, digits = 0)
writeRaster(prediction1, filename = paste0(oDir, "/", 
                                           "outputs", "/", aoi, "_bd1", 
                                           ".tif"), overwrite=TRUE)
writeRaster(prediction2, filename = paste0(oDir, "/", 
                                           "outputs", "/", aoi, "_bd2", 
                                           ".tif"), overwrite=TRUE)

## ------------------------------------------------------------------------
# plot the two maps
plot(prediction1, main="Density for shrubs <1.5m", axes=FALSE)
plot(prediction2, main="Density for shrubs >1.5m", axes=FALSE)

