# Image Classification part 2
# Author: John Mutua
# International Center for Tropical Agriculture (CIAT)
# Last modified: 28-02-2018

## ------------------------------------------------------------------------
# clear your workspace
rm(list = ls(all = TRUE))

## ----eval=TRUE-----------------------------------------------------------
# set variables
iDir <- "E:/Omusati_LDN/Image_classification/data"
oDir <- "E:/Omusati_LDN/Image_classification"
aoi <- "Omusati"
responseCol <- "LC_Code"

## ---- message=FALSE, warning=FALSE---------------------------------------
# load packages
library(rgdal)
library(raster)
library(caret)
library(randomForest)
library(e1071)
library(snow)
library(leaflet)
library(htmlwidgets)
library(sp)
library(RSAGA)
library(plyr)
library(plotKML)

## ---- eval=FALSE---------------------------------------------------------
# .packages = c("rgdal","raster","caret", "randomForest","e1071","snow",
#                "leaflet", "htmlwidgets","sp","RSAGA","plyr", "plotKML")
# .inst <- .packages %in% installed.packages()
#  if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
#  lapply(.packages, require, character.only=TRUE)

## ----help, eval=FALSE----------------------------------------------------
## # this is how you read more on functions
## help(calc)
## ?calc

## ------------------------------------------------------------------------
# set random seed
set.seed(1124)

## ------------------------------------------------------------------------
# import the image into R
img <- brick(paste0(iDir, "/", "TOA", "/", aoi, "_", "TOA", ".dat"))
names(img) <- c(paste0("B", 2:7, coll = ""))

## ---- message=FALSE, warning=FALSE---------------------------------------
# create output folder
if (!dir.exists(paste0(oDir, "/", "outputs")))
  dir.create(paste0(oDir, "/", "outputs"))

## ------------------------------------------------------------------------
# define the function to calculate NDVI 
NDVI.Over <- function(x, y) {
    ndvi <- (y - x) / (x + y)
    return(ndvi)
}

# calculate NDVI and save the output
NDVI <- overlay(x=img[[3]], y=img[[4]], fun=NDVI.Over)
writeRaster(NDVI, filename = paste0(oDir, "/", "outputs", "/", aoi, "_", 
                                    "NDVI", ".tif", sep=""), 
            format = "GTiff", overwrite=TRUE)
plot(NDVI, main = "NDVI map")

## ---- eval=FALSE---------------------------------------------------------
## # make a natural colour visualization of the Landsat image
## plotRGB(img * (img >= 0), r = 4, g = 3, b = 2, scale = 10000)

## ------------------------------------------------------------------------
# import the training data into R
trainData <- shapefile(paste0(iDir, "/", "Training_data", "/", aoi, "_", 
                              "trainingData", ".shp"))

## ------------------------------------------------------------------------
names(trainData)
head(trainData,4)

## ------------------------------------------------------------------------
nsamples <- 300
trainData <- subset(trainData[sample(1:nrow(trainData), nsamples), ])

## ------------------------------------------------------------------------
# add land cover categories based on codes
index <- c(1, 2, 31, 32, 52, 71)
values <- c("Forest", "Woodland", "Bushland", "Grassland", "Water body", 
                "Bare land")
trainData$lulc_name <- values[match(trainData$LC_Code, index)]

## ---- eval=FALSE---------------------------------------------------------
## # inspect the data slot of the trainData object
## trainData@data
## 
## # the 'lulc_name' column is character
## trainData@data$lulc_name
## 
## # explore the data slot of the trainData object
## str(trainData@data$lulc_name)

## ------------------------------------------------------------------------
# plot the training data.
plot(trainData, main="Distribution of training data", axes=FALSE)

## ---- warning=FALSE, message=FALSE---------------------------------------
# extract the pixel values in the training areas
trainSet = data.frame(matrix(vector(), nrow = 0, ncol = length(names(img)) + 1))
for (i in 1:length(unique(trainData[[responseCol]]))){
  category <- unique(trainData[[responseCol]])[i]
  categorymap <- trainData[trainData[[responseCol]] == category,]
  dataSet <- extract(img, categorymap)
  
  if(is(trainData, "SpatialPointsDataFrame")){
    dataSet <- cbind(dataSet, class = as.numeric(category))
    trainSet <- rbind(trainSet, dataSet)
  }
  if(is(trainData, "SpatialPolygonsDataFrame")){
    dataSet <- lapply(dataSet, function(x){cbind(x, class = 
                                                   as.numeric(rep(category, 
                                                                  nrow(x))))})
    df <- do.call("rbind", dataSet)
    trainSet <- rbind(trainSet, df)
  }
}

## ------------------------------------------------------------------------
# partition the data into training and testing
inData <- createDataPartition(y = trainSet$class, p = 0.7, list = FALSE)
training <- trainSet[inData,]
testing <- trainSet[-inData,]

## ------------------------------------------------------------------------
# how does the class look like?
table(training$class)
table(testing$class)

## ------------------------------------------------------------------------
# fit the trandom forest model
rf <- randomForest(as.factor(class) ~ B3 + B4 + B5, data=training,
                    importance=TRUE,
                    ntree=2000)

## ------------------------------------------------------------------------
varImpPlot(rf)

## ------------------------------------------------------------------------
# cluster the predictions
beginCluster()
rf.pred <- clusterR(img, raster::predict, args = list(model = rf))
endCluster()

## ------------------------------------------------------------------------
# predict using testing dataset
rf.predT <- predict(rf, testing)

# overall accuracy
confusionMatrix(rf.predT , testing$class)$overall[1]

# accuracy by class
confusionMatrix(rf.predT , testing$class)$byClass[, 1]

## ------------------------------------------------------------------------
# save as GeoTIFF
writeRaster(rf.pred, filename = paste0(oDir, "/", "outputs", "/", aoi, "_", 
                                       "classified", ".tif", sep=""), 
            format = "GTiff", overwrite=TRUE)

## ------------------------------------------------------------------------
# load study area
studyarea <- shapefile(paste0(iDir, "/", aoi, ".shp"))

# save 'rf.pred' into another object just incase
r <- rf.pred

header <- "Otjiwarongo Land cover 2016"
pal <- colorFactor(c("darkgreen", "forestgreen", "yellow", "magenta", "blue",
                      "red"), values(r), na.color = "transparent")

lb <- labels(c("Forest", "Woodland", "Bushland", "Grassland", "Water body",
                "Bare land"))

# create leaflet
lulc.ll <- leaflet() %>% 
    addProviderTiles("Esri.WorldImagery") %>%
    
    # add data
    addPolygons(data =studyarea, fill = FALSE, stroke = TRUE, 
              color = "#f93", group = "Study area") %>%
  
    addRasterImage(r, colors=pal, opacity=0.5, project = FALSE,
                   group = "Land cover", maxBytes = 4 * 1024 * 1024) %>%

    # add legend
    addLegend("bottomleft", pal=pal, values=values(r), na.label = "NA", 
              labels = lb, opacity = 1, title=header) %>%
  
    # add layers control
    addLayersControl(
      overlayGroups = c("Study area","Land cover"),
      options = layersControlOptions(collapsed = FALSE)
      
    )

# save as widget
saveWidget(lulc.ll, file = paste0(oDir, "/", "outputs", "/", aoi, "_", "LULC", 
                                  ".html", sep = ""))

