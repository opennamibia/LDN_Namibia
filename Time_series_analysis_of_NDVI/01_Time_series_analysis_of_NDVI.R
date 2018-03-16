# Time series analysis of NDVI
# Author: John M
# International Center for Tropical Agriculture (CIAT)
# Last modified: March 16 2018

## ------------------------------------------------------------------------
# clear your work space
rm(list = ls(all = TRUE))

## ----eval=TRUE-----------------------------------------------------------
# set variables
iDir <- "E:/Omusati_LDN/Time_series_analysis_of_NDVI/data"
oDir <- "E:/Omusati_LDN/Time_series_analysis_of_NDVI"
aoi <- "Omusati"
sYr <- "2005"
eYr <- "2015"

## ---- message=FALSE, warning=FALSE---------------------------------------
# load packages
library(rgdal)
library(raster)
library(ggplot2)
library(dplyr)
library(lubridate)
library(grid)
library(gridExtra)
library(dygraphs)
library(signal)
library(leaflet)
library(htmlwidgets)

## ---- eval=FALSE---------------------------------------------------------
## # load packages
## .packages = c("rgdal","raster","ggplot2","dplyr","lubridate","grid",
##               "gridExtra","dygraphs","signal","leaflet","htmlwidgets")
## .inst <- .packages %in% installed.packages()
## if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
## lapply(.packages, require, character.only=TRUE)

## ----message=FALSE, warning=FALSE----------------------------------------
# create output folder
if (!dir.exists(paste0(oDir, "/", "outputs")))
  dir.create(paste0(oDir, "/", "outputs"))

## ------------------------------------------------------------------------
# import rasters
r.list <- list.files(paste0(iDir, "/"), pattern = ".tif$", full.names = T)
r.stack <- stack(r.list)

## ------------------------------------------------------------------------
# insert temporal ID in the dataframe
oldnames<-names(r.stack)
head(oldnames)
year<-substr(oldnames,20,23)
table(year)
julianDay<-substr(oldnames,24,26)

## ------------------------------------------------------------------------
# calculate mean NDVI
avg.NDVI <- as.data.frame(cellStats(r.stack,mean))

## ------------------------------------------------------------------------
# rename NDVI column
names(avg.NDVI) <- "NDVI"

## ------------------------------------------------------------------------
# add julian day as a column
avg.NDVI$julianDay <- julianDay

## ------------------------------------------------------------------------
# check out the new column
class(avg.NDVI$julianDay)

## ------------------------------------------------------------------------
# create time vector and convert to date
tVector<-paste(year,julianDay,sep="-")
timeNDVI<-as.Date(tVector,format = "%Y-%j")

## ------------------------------------------------------------------------
# add date values as a column in the data frame
avg.NDVI$date <- timeNDVI

## ------------------------------------------------------------------------
# check out 'class' of the two columns now
class(avg.NDVI$date)
class(avg.NDVI$julianDay)

## ------------------------------------------------------------------------
# add site and year columns to the data frame
avg.NDVI$site <- paste0(aoi, " ", "Region")
avg.NDVI$year <- year

## ------------------------------------------------------------------------
# plot NDVI by year although this doesn't make sense
ggplot(avg.NDVI, aes(year, NDVI), na.rm=TRUE) +
  geom_point(size=2,colour = "PeachPuff4") + 
  ggtitle(paste0("MODIS Derived NDVI", " - ", aoi, " ", "Region")) +
  xlab("Year") + ylab("Mean NDVI") +
  theme(text = element_text(size=16))

## ------------------------------------------------------------------------
# write NDVI data to .csv
NDVI.Values<-avg.NDVI
row.names(NDVI.Values)<-NULL
write.csv(NDVI.Values, file=paste0(oDir, "/", "outputs", "/", aoi, "_", 
                                   "NDVI_", sYr, "-", eYr, ".csv"))

## ------------------------------------------------------------------------
# add month to data frame
avg.NDVI$month  <- month(avg.NDVI$date)

## ------------------------------------------------------------------------
# subset data by season
avg.NDVI <- avg.NDVI %>% 
  mutate(season = 
           ifelse(month %in% c(12, 1, 2, 3, 4, 5), "Hot-Wet",
                  ifelse(month %in% c(6, 7, 8), "Cool-Dry",
                         ifelse(month %in% c(9, 10, 11), "Hot-Dry", "Error"))))

## ------------------------------------------------------------------------
# did the above work?
head(avg.NDVI$month)
head(avg.NDVI$season)
tail(avg.NDVI$month)
tail(avg.NDVI$season)

## ------------------------------------------------------------------------
# aggregate data by month and calculate standard deviation and error
month.summary<-avg.NDVI %>%
  group_by(month) %>%
  summarise(mean.NDVI=mean(NDVI, na.rm=TRUE),
            sd_NDVI = sd(NDVI), # standard deviation of each group
            n_NDVI = n(),  # sample size per group
            SE_NDVI = sd(NDVI)/sqrt(n())) # standard error of each group

## ------------------------------------------------------------------------
# convert month numeric to month abbreviation
month.summary$month_name <- month.abb[month.summary$month]

## ------------------------------------------------------------------------
# reassign the 'month_name' field to a factor
month.summary$month_name = factor(month.summary$month_name,
                              levels=c('Jan','Feb','Mar',
                                       'Apr','May','Jun','Jul',
                                       'Aug','Sep','Oct',
                                       'Nov','Dec'))

## ------------------------------------------------------------------------
# plot of months vs mean NDVI with standard errors as the error bar
monthly.NDVI.plot<-ggplot(month.summary, aes(month_name, mean.NDVI, group=4)) +
  theme_update(plot.title = element_text(hjust = 0.5)) + # centre title
  geom_line(colour="red") +
  geom_errorbar(aes(ymin=mean.NDVI-sd_NDVI, ymax=mean.NDVI+sd_NDVI),width=0.2) +
    ggtitle(paste0("Average Monthly NDVI", " - ", aoi, " ", "Region")) +
  xlab("Month") + ylab("Mean NDVI") +
  theme(plot.title = element_text(lineheight=.8, face="bold",
                                  size = 16)) +
  theme(text = element_text(size=14)) + geom_point()
monthly.NDVI.plot
ggsave(file=paste0(oDir, "/", "outputs", "/", aoi, "_Monthly_NDVI", ".pdf" ), 
       width = 297, height = 210, units = "mm")

## ------------------------------------------------------------------------
# plot HOT-WET season and save as PDF
Hot.Wet <- subset(month.summary, month >= 12 | month <= 5)
target <- c("12", "5", "4", "3", "2", "1")
Hot.Wet<-Hot.Wet[match(target, Hot.Wet$month),]
Hot.Wet$month_name <- factor(Hot.Wet$month_name, c("Dec", "Jan", "Feb", "Mar", 
                                                   "Apr", "May"))
Hot.Wet.NDVI.Plot<-ggplot(Hot.Wet, aes(month_name, mean.NDVI, group=1)) +
  geom_line(colour="red") +
  ggtitle(paste0("Average NDVI", " - ", "Hot Wet Season")) +
  xlab("Month") + ylab("Mean NDVI") +
  theme(plot.title = element_text(lineheight=.8, face="bold",
                                  size = 16)) +
  theme(text = element_text(size=14)) + geom_point()
Hot.Wet.NDVI.Plot
ggsave(file=paste0(oDir, "/", "outputs", "/", aoi, "_Hot_Wet_NDVI", ".pdf"), 
       width = 297, height = 210, units = "mm")

## ------------------------------------------------------------------------
# plot COOL-DRY season and save as PDF
Cool_Dry <- subset(month.summary, month <= 8 & month >= 6)
Cool_Dry_NDVI_Plot<-ggplot(Cool_Dry, aes(month_name, mean.NDVI, group=1)) +
  geom_line(colour="red") +
  ggtitle(paste0("Average NDVI", " - ", "Cool Dry Season")) +
  xlab("Month") + ylab("Mean NDVI") +
  theme(plot.title = element_text(lineheight=.8, face="bold",
                                  size = 16)) +
  theme(text = element_text(size=14)) + geom_point()
Cool_Dry_NDVI_Plot
ggsave(file=paste0(oDir, "/", "outputs", "/", aoi, "_Cool_Dry_NDVI", ".pdf"), 
       width = 297, height = 210, units = "mm")

## ------------------------------------------------------------------------
# plot HOT DRY season
Hot.Dry <- subset(month.summary, month <= 11 & month >= 9)
Hot.Dry.NDVI.Plot<-ggplot(Hot.Dry, aes(month_name, mean.NDVI, group=1)) +
  geom_line(colour="red") +
  ggtitle(paste0("Average NDVI", " - ", "Hot Dry Season")) +
  xlab("Month") + ylab("Mean NDVI") +
  theme(plot.title = element_text(lineheight=.8, face="bold",
                                  size = 16)) +
  theme(text = element_text(size=14)) + geom_point()
Hot.Dry.NDVI.Plot
ggsave(file=paste0(oDir, "/", "outputs", "/", aoi, "_Hot_Dry_NDVI", ".pdf"), 
       width = 297, height = 210, units = "mm")

## ------------------------------------------------------------------------
# aggregate data by season
season.NDVI<-avg.NDVI %>%
  group_by(year, season) %>%
  summarise(season.NDVI=mean(NDVI, na.rm=TRUE),
            sd_NDVI = sd(NDVI), # standard deviation of each group
            n_NDVI = n(),  # sample size per group
            SE_NDVI = sd(NDVI)/sqrt(n())) # standard error of each group)

## ------------------------------------------------------------------------
# export .csv file of seasonal NDVI values
write.csv(season.NDVI, file=paste0(oDir, "/", "outputs", "/", aoi, 
                                   "_Seasonal_NDVI", "_", 
                                   sYr, "-", eYr, ".csv"))

## ------------------------------------------------------------------------
# plot data by season and save as PDF
Seasonal.NDVI.Plot<-ggplot(season.NDVI, aes(year, season.NDVI, group=3)) +
  geom_line(colour="red") +
  geom_errorbar(aes(ymin=season.NDVI-sd_NDVI, ymax=season.NDVI+sd_NDVI),width=0.2) +
  ggtitle(paste0("Average Seasonal NDVI")) +
  xlab("Year") + ylab("Mean NDVI") +
  theme(plot.title = element_text(lineheight=.8, face="bold",
                                  size = 16)) +
  theme(text = element_text(size=14)) + geom_point()
Seasonal.NDVI.Plot + facet_grid(season ~ .)
ggsave(file=paste0(oDir, "/", "outputs", "/", aoi, "_Seasonal_NDVI", ".pdf"), 
       width = 297, height = 210, units = "mm")

## ------------------------------------------------------------------------
# end of exercise

