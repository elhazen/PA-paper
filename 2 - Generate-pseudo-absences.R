# 2) generate pseudo-absences for tracking data

#' Code to generate pseudo-absences for tracking data for blue whales in California Current
#' and elephant movement data in Etosha National Park, Namibia
#' Elephant and environmental data from Tsalyuk et al. 2018 Ecol. Mon.
#' Blue whale and environmental data from Abrahms et al. 2019 Div. Dist.
#' Acknowledgements for movement and environmental layers for elephants go to 
#' Wayne Getz, Mirium Tsalyuk, Dana Seidel, and associated funding sources
#' Acknowledgements for movement and environmental layers for blue whales go to
#' Briana Abrahms, Heather Welch, Stephanie Brodie, Michael Jacox, Elizabeth Becker, 
#' Steven Bograd, LAdd Irvine, Daniel Palacios, Bruce Mate, Elliott Hazen
#' Code saved to GitHub, Data files saved to Dropbox 
#' 

### Start with Blue Whales
# read in tag data
tags = read.csv(in.csv)

tags$date <- paste(gsub(" ","",paste(tags$Month,tags$Day,tags$Year,sep="/"), fixed=TRUE), " 12:00:00 GMT")
tags$dTime = as.POSIXct(strptime(as.character(tags$date), "%m/%d/%Y %H:%M", tz="GMT"))

# clear sim.alltags if exists in workspace
if (exists('sim.alltags')) rm(sim.alltags)
if (exists('sim.allbufftags')) rm(sim.allbufftags)
if (exists('sim.allbacktags')) rm(sim.allbacktags)

#tags to be run:
#tagsleft<-unique(tags$eventid)[7:64]

for (tagid in unique(tags$tag)[unique(tags$tag)>0]){    
  #Simulate CRW
  sim.alldata <- createCRW(tags, tagid, n.sim=100)
  #Simulate reverse CRW
  sim.alldata <- createCRW(tags, tagid, n.sim=100, reverse = TRUE)
  graphics.off()
  #Simulate buffers
  sim.buffer.alldata <- createbufferabsence(tags, tagid, n.sim=100, buffdist = 1000.000)
  #Simulate background
  sim.background.alldata <- createbackgroundabsence(tags, tagid, n.sim=20)  
}  # end for (tagid in unique(tags$ptt)){

tags$id<-tags$tag
for (tagid in unique(tags$tag)[unique(tags$tag)>0]){  
  sim.buffer.alldata <- createbufferabsence(tags, tagid, n.sim=100, buffdist = 1000.000)
  graphics.off()
}

#Next run Elephant pseudo-tracks
tags$date <- as.character(tags$date)
tags$dTime = as.POSIXct(strptime(as.character(tags$date), "%Y-%m-%dT%H:%M:%S", tz="GMT"))
if (exists('sim.alltags')) rm(sim.alltags)

sputm <- SpatialPoints(cbind(tags$x,tags$y), proj4string=CRS("+proj=utm +south +zone=33K +datum=WGS84 "))  
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))

tags$lat<-coordinates(spgeo)[,2]
tags$long<-coordinates(spgeo)[,1]

tags.hr<-tags[seq(1,dim(tags)[1],by=3),]

saveRDS(tags.hr,file="ElephantData_v1.RDS")
tags.hr<-readRDS(file="ElephantData_v1.RDS")

#tags to be run:
tagsleft<-unique(tags.hr$id)[10:15]

for (tagid in tagsleft){    #CAN SUBSET TAGS USING: metad[metad[,1]>510241600,1]  [unique(as.character(tags$id))>0]
  #Simulate CRW
  sim.alldata <- create.land.CRW(tags.hr, tagid, n.sim=100)
  #Simulate reverse CRW
  sim.alldata <- create.land.CRW(tags.hr, tagid, n.sim=100, reverse = TRUE)
  graphics.off()
  #Simulate buffer absences
  sim.buffer.alldata <- createbufferabsence(tags.hr, tagid, n.sim=100, buffdist = 4.4)
  #Simulate background absences
  sim.background.alldata <- createbackgroundabsence(tags.hr, tagid, n.sim=100)  
}  # end for (tagid in unique(tags.hr$ptt)){