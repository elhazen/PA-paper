# $Id$
# build GAM and BRT models (and GLM?) for elephants and blue whales


##### test\
### read data start models

datawd<-('./Data/XtractedData')

#source('ImportDatasets.R') # to rerun data-org code OR

bwhaleabsbuffer <- readRDS("./Data/bwhaleabsbuffer.RDS")
bwhaleabsbackgrnd <- readRDS("./Data/bwhaleabsbackgrnd.RDS")
bwhaleabsCRW <- readRDS("./Data/bwhaleabsCRW.RDS")
bwhaleabsCRWrev <- readRDS("./Data/bwhaleabsCRWrev.RDS")
ephantabsbuffer <- readRDS("./Data/ephantabsbuffer.RDS")
ephantabsbackgrnd <- readRDS("./Data/ephantabsbackgrnd.RDS")
ephantabsCRW <- readRDS("./Data/ephantabsCRW.RDS")
ephantabsCRWrev <- readRDS("./Data/ephantabsCRWrev.RDS")
bwhalepresence <- readRDS("./Data/bwhalepresence.RDS")
elephantdata <- readRDS("./Data/elephantdata.RDS")

#### RUN TO RE-DO DATA SAMPLES
createWhaleSamples(bwhalepresence,bwhaleabsbuffer,bwhaleabsbackgrnd,bwhaleabsCRW,bwhaleabsCRWrev,tag="reduce_space")

tag="final"
bwhaledata_buff <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_buff_",tag,".RDS")))
bwhaledata_back <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_back_",tag,".RDS")))
bwhaledata_CRW <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_CRW_",tag,".RDS")))
bwhaledata_rev <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_rev_",tag,".RDS")))

extent1<-minmaxextent(bwhaledata_buff,bwhaledata_back,bwhaledata_CRW,bwhaledata_rev)
extent1<-quartileextent(bwhaledata_buff,bwhaledata_back,bwhaledata_CRW,bwhaledata_rev)
bwhaledata_buff<-bwhaledata_buff[bwhaledata_buff$lat<extent1[4]&bwhaledata_buff$lat>extent1[3]&bwhaledata_buff$lon<extent1[2]&bwhaledata_buff$lon>extent1[1],]
bwhaledata_back<-bwhaledata_back[bwhaledata_back$lat<extent1[4]&bwhaledata_back$lat>extent1[3]&bwhaledata_back$lon<extent1[2]&bwhaledata_back$lon>extent1[1],]
bwhaledata_CRW<-bwhaledata_CRW[bwhaledata_CRW$lat<extent1[4]&bwhaledata_CRW$lat>extent1[3]&bwhaledata_CRW$lon<extent1[2]&bwhaledata_CRW$lon>extent1[1],]
bwhaledata_rev<-bwhaledata_rev[bwhaledata_rev$lat<extent1[4]&bwhaledata_rev$lat>extent1[3]&bwhaledata_rev$lon<extent1[2]&bwhaledata_rev$lon>extent1[1],]

tag="final"
bwhaledata_buff <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_buff_",tag,".RDS")))
bwhaledata_back <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_back_",tag,".RDS")))
bwhaledata_CRW <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_CRW_",tag,".RDS")))
bwhaledata_rev <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_rev_",tag,".RDS")))

### predictors > 1%
bwhaleGAMM_CRW<-gamm(presabs~s(SST,bs="ts", k=5)+s(SST_SD,bs="ts", k=5)+s(Oxygen_100m,bs="ts", k=5)+s(log10(Chla_4km_8day+0.001),bs="ts", k=5)+s(log10(Chla_25km_monthly+0.001),bs="ts", k=5)+s(SLA,bs="ts", k=5)+s(Bathymetry,bs="ts", k=5), random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(bwhaledata_CRW)) #get(paste("GAMdataRun",l,sep=''))
bwhaleGLMM_CRW<-gamm(presabs~SST + SST_SD + Oxygen_100m + log10(Chla_4km_8day+0.001) + log10(Chla_25km_monthly+0.001) + SLA + Bathymetry, random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(bwhaledata_CRW))
try(bwhaleBRT.lr005.CRW<-gbm.fixed(data=BRTtransformDataFrame_BW(bwhaledata_CRW),gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, family="bernoulli",tree.complexity=5, learning.rate = 0.005, n.trees = 2000, bag.fraction=0.75))

bwhaleGAMM_buff<-gamm(presabs~s(SST,bs="ts", k=5)+s(SST_SD,bs="ts", k=5)+s(Oxygen_100m,bs="ts", k=5)+s(log10(Chla_4km_8day+0.001),bs="ts", k=5)+s(log10(Chla_25km_monthly+0.001),bs="ts", k=5)+s(SLA,bs="ts", k=5)+s(Bathymetry,bs="ts", k=5), random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(bwhaledata_buff)) #get(paste("GAMdataRun",l,sep=''))
bwhaleGLMM_buff<-gamm(presabs~SST + SST_SD + Oxygen_100m + log10(Chla_4km_8day+0.001) + log10(Chla_25km_monthly+0.001) + SLA + Bathymetry, random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(bwhaledata_buff))
try(bwhaleBRT.lr005.buff<-gbm.fixed(data=BRTtransformDataFrame_BW(bwhaledata_buff),gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, family="bernoulli",tree.complexity=5, learning.rate = 0.005, n.trees = 2000, bag.fraction=0.75))

bwhaleGAMM_back<-gamm(presabs~s(SST,bs="ts", k=5)+s(SST_SD,bs="ts", k=5)+s(Oxygen_100m,bs="ts", k=5)+s(log10(Chla_4km_8day+0.001),bs="ts", k=5)+s(log10(Chla_25km_monthly+0.001),bs="ts", k=5)+s(SLA,bs="ts", k=5)+s(Bathymetry,bs="ts", k=5), random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(bwhaledata_back)) #get(paste("GAMdataRun",l,sep=''))
bwhaleGLMM_back<-gamm(presabs~SST + SST_SD + Oxygen_100m + log10(Chla_4km_8day+0.001) + log10(Chla_25km_monthly+0.001) + SLA + Bathymetry, random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(bwhaledata_back))
try(bwhaleBRT.lr005.back<-gbm.fixed(data=BRTtransformDataFrame_BW(bwhaledata_back),gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, family="bernoulli",tree.complexity=5, learning.rate = 0.005, n.trees = 2000, bag.fraction=0.75))

bwhaleGAMM_rev<-gamm(presabs~s(SST,bs="ts", k=5)+s(SST_SD,bs="ts", k=5)+s(Oxygen_100m,bs="ts", k=5)+s(log10(Chla_4km_8day+0.001),bs="ts", k=5)+s(log10(Chla_25km_monthly+0.001),bs="ts", k=5)+s(SLA,bs="ts", k=5)+s(Bathymetry,bs="ts", k=5), random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(bwhaledata_rev)) #get(paste("GAMdataRun",l,sep=''))
bwhaleGLMM_rev<-gamm(presabs~SST + SST_SD + Oxygen_100m + log10(Chla_4km_8day+0.001) + log10(Chla_25km_monthly+0.001) + SLA + Bathymetry, random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(bwhaledata_rev))
try(bwhaleBRT.lr005.rev<-gbm.fixed(data=BRTtransformDataFrame_BW(bwhaledata_rev),gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, family="bernoulli",tree.complexity=5, learning.rate = 0.005, n.trees = 2000, bag.fraction=0.75))

tag="final"
save(bwhaleGAMM_CRW,bwhaleGLMM_CRW,bwhaleBRT.lr005.CRW,bwhaleGAMM_rev,bwhaleGLMM_rev,bwhaleBRT.lr005.rev,bwhaleGAMM_back,bwhaleGLMM_back,bwhaleBRT.lr005.back,bwhaleGAMM_buff,bwhaleGLMM_buff,bwhaleBRT.lr005.buff,file=paste0("bwhalemodels_",tag,".RData"))
save(bwhaledata_buff,bwhaledata_back,bwhaledata_CRW,bwhaledata_rev, file=paste0("bwhaledata_",tag,".RData"))
  
### Make an elephant dataset for buffer scenario

createElephantSamples(ephantpresence,ephantabsbuffer,ephantabsbackgrnd,ephantabsCRW,ephantabsCRWrev,tag="test2")

tag="test"
ephantdata_buff <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/ephantdata_buff_",tag,".RDS")))
ephantdata_back <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/ephantdata_back_",tag,".RDS")))
ephantdata_CRW <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/ephantdata_CRW_",tag,".RDS")))
ephantdata_rev <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/ephantdata_rev_",tag,".RDS")))

ephantdata_buff$dTime <- fixdTime(ephantdata_buff$dTime)
ephantdata_back$dTime <- fixdTime(ephantdata_back$dTime)
ephantdata_CRW$dTime <- fixdTime(ephantdata_CRW$dTime)
ephantdata_rev$dTime <- fixdTime(ephantdata_rev$dTime)

ephantdata_buff[ephantdata_buff$iteration==0,]<-UTM2WGS(ephantdata_buff[ephantdata_buff$iteration==0,])
ephantdata_back[ephantdata_back$iteration==0,]<-UTM2WGS(ephantdata_back[ephantdata_back$iteration==0,])
ephantdata_CRW[ephantdata_CRW$iteration==0,]<-UTM2WGS(ephantdata_CRW[ephantdata_CRW$iteration==0,])
ephantdata_rev[ephantdata_rev$iteration==0,]<-UTM2WGS(ephantdata_rev[ephantdata_rev$iteration==0,])

##cut PA's by total extent
tag="final"
extent1<-c(min(ephantdata_back$long),max(ephantdata_back$long),min(ephantdata_back$lat),max(ephantdata_back$lat))
ephantdata_buff<-ephantdata_buff[ephantdata_buff$lat<extent1[4]&ephantdata_buff$lat>extent1[3]&ephantdata_buff$lon<extent1[2]&ephantdata_buff$lon>extent1[1],]
ephantdata_back<-ephantdata_back[ephantdata_back$lat<extent1[4]&ephantdata_back$lat>extent1[3]&ephantdata_back$lon<extent1[2]&ephantdata_back$lon>extent1[1],]
ephantdata_CRW<-ephantdata_CRW[ephantdata_CRW$lat<extent1[4]&ephantdata_CRW$lat>extent1[3]&ephantdata_CRW$lon<extent1[2]&ephantdata_CRW$lon>extent1[1],]
ephantdata_rev<-ephantdata_rev[ephantdata_rev$lat<extent1[4]&ephantdata_rev$lat>extent1[3]&ephantdata_rev$lon<extent1[2]&ephantdata_rev$lon>extent1[1],]

map = st_read("./Data/Etosha shapefiles", "enp fence poly")
pnts_sf <- st_as_sf(ephantdata_buff, coords = c('long', 'lat'), crs = st_crs(map))
pnts <- pnts_sf %>% mutate(
  intersection = as.integer(st_intersects(geometry, map))) 
ephantdata_buff<-ephantdata_buff[!is.na(pnts$intersection),]

pnts_sf <- st_as_sf(ephantdata_back, coords = c('long', 'lat'), crs = st_crs(map))
pnts <- pnts_sf %>% mutate(
  intersection = as.integer(st_intersects(geometry, map))) 
ephantdata_back<-ephantdata_back[!is.na(pnts$intersection),]

pnts_sf <- st_as_sf(ephantdata_CRW, coords = c('long', 'lat'), crs = st_crs(map))
pnts <- pnts_sf %>% mutate(
  intersection = as.integer(st_intersects(geometry, map))) 
ephantdata_CRW<-ephantdata_CRW[!is.na(pnts$intersection),]

pnts_sf <- st_as_sf(ephantdata_rev, coords = c('long', 'lat'), crs = st_crs(map))
pnts <- pnts_sf %>% mutate(
  intersection = as.integer(st_intersects(geometry, map))) 
ephantdata_rev<-ephantdata_rev[!is.na(pnts$intersection),]

save(ephantdata_buff,ephantdata_back,ephantdata_CRW,ephantdata_rev,file=paste0("ephantdata",tag,".RData"))
load(paste0("ephantmodels",tag,".RData"))

ephantGAMM_CRW<-gamm(presabs~s(roaddist,bs="ts",k=5)+s(NDVI,bs="ts",k=5)+s(waterdist,bs="ts", k=5), random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(ephant_downsample_6hrs(ephantdata_CRW)))
ephantGLMM_CRW<-gamm(presabs~roaddist+NDVI+waterdist, random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(ephant_downsample_6hrs(ephantdata_CRW)))
try(ephantBRT.lr005.CRW<-gbm.fixed(data=ephant_downsample_6hrs(ephantdata_CRW),gbm.x=6:8, gbm.y=9, family="bernoulli",tree.complexity=5, learning.rate = 0.005, n.trees = 2000, bag.fraction=0.75))

ephantGAMM_buff<-gamm(presabs~s(roaddist,bs="ts",k=5)+s(NDVI,bs="ts",k=5)+s(waterdist,bs="ts", k=5), random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(ephant_downsample_6hrs(ephantdata_buff)))
ephantGLMM_buff<-gamm(presabs~roaddist+NDVI+waterdist, random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(ephant_downsample_6hrs(ephantdata_buff)))
try(ephantBRT.lr005.buff<-gbm.fixed(data=ephant_downsample_6hrs(ephantdata_buff),gbm.x=6:8, gbm.y=9, family="bernoulli",tree.complexity=5, learning.rate = 0.005, n.trees = 2000, bag.fraction=0.75))

ephantGAMM_back<-gamm(presabs~s(roaddist,bs="ts",k=5)+s(NDVI,bs="ts",k=5)+s(waterdist,bs="ts", k=5), random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(ephant_downsample_6hrs(ephantdata_back)))
ephantGLMM_back<-gamm(presabs~roaddist+NDVI+waterdist, random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(ephant_downsample_6hrs(ephantdata_back)))
try(ephantBRT.lr005.back<-gbm.fixed(data=ephant_downsample_6hrs(ephantdata_back),gbm.x=6:8, gbm.y=9, family="bernoulli",tree.complexity=5, learning.rate = 0.005, n.trees = 2000, bag.fraction=0.75))

ephantGAMM_rev<-gamm(presabs~s(roaddist,bs="ts",k=5)+s(NDVI,bs="ts",k=5)+s(waterdist,bs="ts", k=5), random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(ephant_downsample_6hrs(ephantdata_rev)))
ephantGLMM_rev<-gamm(presabs~roaddist+NDVI+waterdist, random=list(tag=~1),family=binomial, niterPQL=50, data=na.omit(ephant_downsample_6hrs(ephantdata_rev)))
try(ephantBRT.lr005.rev<-gbm.fixed(data=ephant_downsample_6hrs(ephantdata_rev),gbm.x=6:8, gbm.y=9, family="bernoulli",tree.complexity=5, learning.rate = 0.005, n.trees = 2000, bag.fraction=0.75))

save(ephantGAMM_CRW,ephantGLMM_CRW,ephantBRT.lr005.CRW,ephantGAMM_rev,ephantGLMM_rev,ephantBRT.lr005.rev,ephantGAMM_back,ephantGLMM_back,ephantBRT.lr005.back,ephantGAMM_buff,ephantGLMM_buff,ephantBRT.lr005.buff,file=paste0("ephantmodels",tag,".RData"))
load(paste0("ephantmodels",tag,".RData"))

