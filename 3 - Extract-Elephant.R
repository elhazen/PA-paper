#' Code to extract environmental covariates for elephant movement data in Etosha National Park, Namibia
#' Elephant and environmental data from Tsalyuk et al. 2018 Ecol. Mon.
#' Acknowledgements for movement and environmental layers go to 
#' Wayne Getz, Mirium Tsalyuk, Dana Seidel, and associated funding sources
#' Data files saved to Dropbox in PseudoAbsence-MS folder
#' Code by Briana Abrahms

#1. ------------> Load libraries and global objects ####
library(sp) # For creating spatial objects and reprojecting them
library(rgdal) # For shapefiles
library(raster) # For rasters
library(rgeos) #For calculating distances

shapefilesFolder <- './Data/Etosha shapefiles/'
eleFolder <- './Data/'
crwFolder <- './Data/CRW_output/'


#2. ------------> Load environmental layer shapefiles and reproject to UTM ####

## Reprojecting to UTM so distance units are in meters

##Roads
roadlayer <- rgdal::readOGR(paste0(shapefilesFolder, "enp roads.shp"))
proj4string(roadlayer) <- CRS('+proj=longlat +datum=WGS84')
roadlayer <-sp::spTransform(roadlayer, CRS("+proj=utm +south +zone=33 +ellps=WGS84"))

##NDVI
NDVIlayer <- raster::raster(paste0(shapefilesFolder, "ndvi_mean_utm.tif"))

##Water
waterlayer <- rgdal::readOGR(paste0(shapefilesFolder, "functional water.shp"))
proj4string(waterlayer) <- CRS('+proj=longlat +datum=WGS84')
waterlayer <-sp::spTransform(waterlayer, CRS("+proj=utm +south +zone=33 +ellps=WGS84"))


#3. ------------> Function to extract covariates ####

Etosha_extract <- function(data=data){
  
  #1.------------> Convert points to SpatialPoints and reproject to UTM
  pts <- sp::SpatialPoints(cbind(data$long, data$lat), proj4string = CRS('+proj=longlat +datum=WGS84'))
  pts <- spTransform(pts, CRS("+proj=utm +south +zone=33 +ellps=WGS84")) # reproject to UTM
  
  #1.------------> Extract
  print('EXTRACTING DISTANCE TO ROADS')
  data$roaddist <- apply(rgeos::gDistance(pts,roadlayer,byid=TRUE),2,min)
  
  print('EXTRACTING NDVI')
  data$NDVI <- raster::extract(NDVIlayer, pts)
  
  print('EXTRACTING DISTANCE TO WATER')
  data$waterdist <- apply(rgeos::gDistance(pts,waterlayer,byid=TRUE),2,min)
  
  return(data)
  
}

#3b. -----------> Create rasterstack 

cs = res(NDVIlayer) # Raster cell size

# Create extent and coerce to SpatialPolygons 
e <- as( raster::extent(roadlayer), "SpatialPolygons")
proj4string(e) <- CRS("+proj=utm +south +zone=33 +ellps=WGS84")
#e <- sp::spTransform(e, CRS('+proj=longlat +datum=WGS84'))
class(e)
plot(e)

# Create raster from defined extent
r <- raster(e, resolution = cs)

r[] <- 1:ncell(r)
plot(r)
cat("\n", "Number of cells in raster: ", ncell(r), "\n")  

r_spdf <- as(r, "SpatialPoints")
r_spdf <- sp::spTransform(r_spdf, CRS('+proj=longlat +datum=WGS84'))

r_spdf$lat<-r_spdf$y
r_spdf$long<-r_spdf$x

class(r_spdf)

d_gridded <- Etosha_extract(data=r_spdf)
d_gridded_utm<-sp::spTransform(d_gridded, CRS("+proj=utm +south +zone=33 +ellps=WGS84"))

d_NDVI <- rasterize(d_gridded_utm,r,"NDVI", fun=mean)
d_waterdist <- rasterize(d_gridded_utm,r,"waterdist", fun=mean)
d_roaddist <- rasterize(d_gridded_utm,r,"roaddist", fun=mean)
d_lat <- rasterize(d_gridded_utm,r,"lat", fun=mean)
d_long <- rasterize(d_gridded_utm,r,"long", fun=mean)
  
ephant_rasterstack<-raster::stack(d_NDVI,d_waterdist,d_roaddist,d_lat,d_long)
names(ephant_rasterstack)<-c("NDVI","waterdist","roaddist","lat","long")  

saveRDS(ephant_rasterstack,file=paste0(shapefilesFolder,"ephant_rasterdata.RDS"))

#4. ------------> Extract empirical elephant data ####

## Load data
d <- read.csv(paste0(eleFolder, "Elephants_EcolMonographs2018_cleaned.csv"))

## Add lat/long because data must contain 'long' ,'lat' to match crw files
xy <- sp::SpatialPoints(cbind(d$x, d$y), proj4string = CRS("+proj=utm +south +zone=33 +ellps=WGS84"))
latlong <- spTransform(xy, CRS('+proj=longlat +datum=WGS84')) # reproject to lat/long
d$long <- latlong@coords[,1]; d$lat <- latlong@coords[,2]

## Extract
d <- Etosha_extract(data=d)

## Save
write.csv(d,paste0(eleFolder,'XtractedData/Extracted_Elephants.csv'), row.names = FALSE)


#5. ------------> Extract PseudoAbsence data ####
files <- list.files(crwFolder,pattern=glob2rx("*AG*csv*"))

#Buffer and Background
for (f in 1:30){
  print(paste0("Extracting absence file: " ,f, " of 60"))
  d <- read.csv(paste0(crwFolder,files[f]))
  d <- Etosha_extract(data=d)
    write.csv(d,paste0(eleFolder,'XtractedData/Extracted_',files[f]), row.names = FALSE)
} 

#CRW and Reverse CRW
for (f in 31:60){
  print(paste0("Extracting absence file: " ,f, " of 60"))
  d <- read.csv(paste0(crwFolder,files[f]))
  d$long <- d$x; d$lat <- d$y
  d <- Etosha_extract(data=d)
  write.csv(d,paste0(eleFolder,'XtractedData/Extracted_',files[f]), row.names = FALSE)
} 

