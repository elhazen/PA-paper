### LOAD PACKAGES
require(tidyverse)
require(mapdata)
require(maps)
require(ggplot2)
require(ggmap)
require(lubridate)
require(RColorBrewer)
require("viridis")
require("gridExtra")
require("scales")
require(mgcv)
require(rnaturalearth)
require(rgeos)
require(sf)
#require("rgdal") # requires sp, will use proj.4 if installed

source("LOOModels.R")
source("ModelFunctions.R")

### LOAD FUNCTIONS

predGAM <- function(model, savefilename = "GAMplot.pdf", rasterdata, rasterextent = newdata.r, logtransform=FALSE, rescale0to1=FALSE, dt = "unknown date", makeplot=FALSE, ggplottitle="Map of predicted blue whale habitat"){
  
  if (mean((rasterdata$SST), na.rm=TRUE) > 100) rasterdata$SST <- rasterdata$SST-273.15
  newdata<-predict.gam(model, rasterdata, se.fit=TRUE, type="response")
  #newdata<-raster::predict(POPA_CPUE_gam,rasterstack, fun=predict.gam, na.rm=TRUE, type="response")
  #newdata.r[]<-newdata$fit
  if (logtransform) {
    temp<-exp(setValues(rasterextent,as.vector(newdata$fit)))-0.001# %>% range01(.)
  } else {
    temp<-(setValues(rasterextent,as.vector(newdata$fit)))
  }
  
  if (rescale0to1) temp <- temp %>% range01(.)
  
  ggplottitle <- paste0("Blue whale habitat likelihood predicted on ",as.Date(dt))
  
  ### plot the data, eventually add fishing locations?
  if (makeplot){
    pdf (file=savefilename, height=8, width=8)
    par(mar=c(5,4,4,2))
    print(AzoresPlot(temp, xvis=c(-130,-110), yvis=c(30,50), addNarrow = FALSE, addkmscale = FALSE, plottitle = ggplottitle))
    dev.off()
  }
  
  return(temp)
  ### c(-50,-8,8,66)
}

predBRT <- function(model, savefilename = "BRTplot.pdf", rasterdata, rasterextent = newdata.r, logtransform=FALSE, rescale0to1=FALSE, dt = "unknown date", makeplot=FALSE, ggplottitle="Map of predicted blue whale habitat"){
  if (mean(rasterdata$SST, na.rm=TRUE) > 100) rasterdata$SST <- rasterdata$SST-273.15
  temp<-predict.gbm(model, newdata=rasterdata,
                    n.trees=model$gbm.call$best.trees, type="response")
  if (logtransform) {
    temp<-exp(setValues(newdata.r,as.vector(temp)))-0.001 # %>% range01(.)
  } else {
    temp<-(setValues(newdata.r,as.vector(temp))) # %>% range01(.)
  }
  
  if (rescale0to1) temp <- temp %>% range01(.)
  
  if (makeplot){
    pdf (file=savefilename, height=8, width=8)
    par(mar=c(5,4,4,2))
    print(AzoresPlot(temp, xvis=c(-130,-110), yvis=c(30,50), addNarrow = FALSE, addkmscale = FALSE, plottitle = ggplottitle))
    dev.off()
  }
  return(temp)
}

automatePredictions <- function(dts, GAMmodel=bestGAMmodel, BRTmodel=bestBRTmodel, savefilenms = c("",""), xvis=c(-50,-10), yvis=c(25,45), zvisBRT=NA, zvisGAM=NA, tunapts=NA, saveCPUEsummary=FALSE, adddate2plot=FALSE){
  #dts<-readRDS("./RasterData/ACastdates.RDS")
  dt <- dts[1]
  #dt <-as.POSIXct("2006-08-01 12:00:00")
  outDir<-paste(getwd(),"RasterData",as_date(dt),sep='/')
  rasterdata <- readRDS(file=paste0(outDir,"/rasterdata.RDS"))
  #rasterdata <- readRDS(file="/Users/elliotthazen/Documents/R/github/POPA/RasterData/2003-05-15/rasterdata.RDS")
  newdata.r<-rasterdata[[1]]
  rasterdata<-addLayer(rasterdata,log(rasterdata$Chla_4km_8day + 0.001))
  names(rasterdata)[20]<-"logChl8D"
  rasterdata<-addLayer(rasterdata,log(rasterdata$Chla_25km_monthly + 0.001))
  names(rasterdata)[21]<-"logChlM"
  rasterdata<-addLayer(rasterdata,log(rasterdata$Chla_4km_8day_lag1mo + 0.001))
  names(rasterdata)[22]<-"logChl8D1mo"
  rasterdata<-addLayer(rasterdata,log(rasterdata$Chla_25km_monthly_lag1mo + 0.001))
  names(rasterdata)[23]<-"logChlM1mo"
  
  fishingCat <- newdata.r
  fishingCat[]<-as.factor("Isca")
  names(fishingCat)<-"fishingCategory"
  BaitType<-fishingCat
  BaitType[]<-as.factor("sardinha")
  names(BaitType)<-"BaitType"
  fishingtime<-newdata.r
  fishingtime[]<-50
  names(fishingtime)<-"timefishing"
  
  # if skipjack grab this
  if (!is.na(str_extract(savefilenms[1],"KP"))) {
    abundata<-read_csv("SkipjackTunaAbundance.csv")
  } else  if (!is.na(str_extract(savefilenms[1],"TO"))) { # if bigeye grab this
    abundata<-read_csv("BigeyeTunaAbundance.csv")
  }
  
  abunindex<-newdata.r
  abunindex[]<-Extract_abundance(as.data.frame(dt),abundata)$abunindex
  names(abunindex)<-"abunindex"
  
  rasterdata<-addLayer(rasterdata,fishingCat,BaitType,fishingtime,abunindex)
  
  #names(rasterdata)[20]<-"fishingCategory"
  
  #rasterdata<-addLayer(rasterdata,"Isca")
  #names(rasterdata)[21]<-"fishingCategory"
  
  rasterstack=rasterdata
  rasterdata=as.data.frame(rasterdata,stringsAsFactors=F) ## predicts on data.frame, not on a stack
  #rasterdata=rasterToPoints(rasterstack[[1]],spatial=TRUE) ## predicts on data.frame, not on a stack
  names(rasterdata)[24]<-"fishingCategory"
  names(rasterdata)[25]<-"BaitType"
  names(rasterdata)[26]<-"timefishing"
  if (mean(rasterdata$SST, na.rm=TRUE) > 100) rasterdata$SST <- rasterdata$SST-273.15
  
  if (length(tunapts)>1){
    tunaptsbyyear<-tunapts[as.numeric(tunapts$Year) == year(dt) & as.numeric(tunapts$Month) == month(dt),]
    if (!is.numeric(dim(tunaptsbyyear)[1]) | (dim(tunaptsbyyear)[1]==0)) tunaptsbyyear<-NA
  } else {
    tunaptsbyyear<-NA
  }
  
  ### Predict GAM for CPUE
  temp<-predGAM(GAMmodel, savefilename = paste0("GAMplot_",savefilenms,as.Date(dt),".pdf"), rasterdata, rasterextent = newdata.r, logtransform=TRUE, rescale0to1=FALSE, dt = dt)
  
  pdf (file=paste0("./predmaps/GAMplot_",savefilenms[1],"_",as.Date(dt),".pdf"), height = 8, width = 8)
  par(mar=c(5,4,4,2))
  print(AzoresPlot(temp, xvis=xvis, yvis=yvis, addNarrow = FALSE, addkmscale = FALSE, plottitle = paste0("Tuna catch predicted on ",as.Date(dt)), labellegend="CPUE", tunapts=tunaptsbyyear,adddate2plot=TRUE, zvis=zvisGAM))
  dev.off()
  
  saveRDS(temp,file=paste0("./preddata/GAMprediction_",savefilenms[1],"_",as.Date(dt),".RDS"))
  
  ### Predict BRT for # fish
  temp2<-predBRT(BRTmodel, savefilename = paste0("BRTplot_",as.Date(dt),".pdf"), rasterdata, rasterextent = newdata.r, logtransform=TRUE, rescale0to1=FALSE, dt = dt)
  
  pdf (file=paste0("./predmaps/BRTplot_",savefilenms[2],"_",as.Date(dt),".pdf"), height=8, width=8)
  par(mar=c(5,4,4,2))
  print(AzoresPlot(temp2, xvis=xvis, yvis=yvis, addNarrow = FALSE, addkmscale = FALSE, plottitle = paste0("Tuna catch predicted on ",as.Date(dt)), labellegend="# caught",tunapts=tunaptsbyyear,adddate2plot=TRUE, zvis=zvisBRT))
  dev.off()
  
  saveRDS(temp2,file=paste0("./preddata/BRTprediction_",savefilenms[1],"_",as.Date(dt),".RDS"))
  
  if (saveCPUEsummary){
    ### read in shapefiles
    #    princess.alice.shp<-st_read("./shapefiles/Box_Pricesa_Alice.shp")
    EEZ.shp<-st_read("./shapefiles/eez/eez.shp")
    POPA.shp<-st_read("./shapefiles/POPA95KD/homerange95.shp")
    #    princess.alice.poly<-readShapePoly("./shapefiles/Box_Pricesa_Alice.shp")
    #    EEZ.poly <- readShapePoly("./shapefiles/eez/eez.shp")
    
    #extract raster cell count (sum) within each polygon area (poly)
    ex <- raster::extract(temp, EEZ.shp, fun=sum, na.rm=TRUE, df=TRUE)
    ex2 <- raster::extract(temp2, EEZ.shp, fun=sum, na.rm=TRUE, df=TRUE)
    #      ex.PA <- raster::extract(temp, princess.alice.shp, fun=sum, na.rm=TRUE, df=TRUE)
    #      ex.PA2 <- raster::extract(temp2, princess.alice.shp, fun=sum, na.rm=TRUE, df=TRUE)
    ex.POPA <- raster::extract(temp, POPA.shp, fun=sum, na.rm=TRUE, df=TRUE)
    ex2.POPA <- raster::extract(temp2, POPA.shp, fun=sum, na.rm=TRUE, df=TRUE)
    
    #write to a data frame
    df <- data.frame(ex) ## in CPUE so multiply by mean fishing time of 50 mins
    df2 <- data.frame(ex2) ## in # catch
    df3 <- data.frame(ex.POPA) ## in CPUE so multiply by mean fishing time of 50 mins
    df4 <- data.frame(ex2.POPA) ## in # catch
    #df5 <- data.frame(ex.PA) ## in CPUE so multiply by mean fishing time of 50 mins
    #df6 <- data.frame(ex.PA2) ## in # catch
    
    df$layer <- df$layer *50
    df3$layer <- df3$layer *50
    
    #df<-rbind(df,df2)
    
    #write to a CSV file
    write.csv(df, file = paste0("./predmaps/CPUE_EEZ",savefilenms[1],"_",as.Date(dt),".csv"))
    write.csv(df2, file = paste0("./predmaps/CPUE_EEZ",savefilenms[2],"_",as.Date(dt),".csv"))
    write.csv(df3, file = paste0("./predmaps/CPUE_POPA",savefilenms[1],"_",as.Date(dt),".csv"))
    write.csv(df4, file = paste0("./predmaps/CPUE_POPA",savefilenms[2],"_",as.Date(dt),".csv"))
    #write.csv(df5, file = paste0("./predmaps/CPUE_PA",savefilenms[1],"_",as.Date(dt),".csv"))
    #write.csv(df6, file = paste0("./predmaps/CPUE_PA",savefilenms[2],"_",as.Date(dt),".csv"))
    
    #saveRDS(a1,"./predmaps/CPUE_EEZ_",savefilenms[1],"_",as.Date(dt),".pdf")
    #saveRDS(a2,"./predmaps/CPUE_EEZ_",savefilenms[2],"_",as.Date(dt),".pdf")
    #saveRDS(PA1,"./predmaps/CPUE_EEZ_",savefilenms[1],"_",as.Date(dt),".pdf")
    #saveRDS(PA2,"./predmaps/CPUE_EEZ_",savefilenms[2],"_",as.Date(dt),".pdf")
    
    
  }
}

range01 <- function(r){
  r.min = cellStats(r, "min")
  r.max = cellStats(r, "max")
  r.scale <- ((r - r.min) / (r.max - r.min))
  return(r.scale) #(r-rmin)/(rmax-rmin)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


### LOAD DATA
tag="final"

bwhaledata_buff <- (readRDS(paste0("./Data/bwhaledata_buff_",tag,".RDS")))
bwhaledata_back <- (readRDS(paste0("./Data/bwhaledata_back_",tag,".RDS")))
bwhaledata_CRW <- (readRDS(paste0("~./Data/bwhaledata_CRW_",tag,".RDS")))
bwhaledata_rev <- (readRDS(paste0("~./Data/bwhaledata_rev_",tag,".RDS")))

tag="BRTeqGAMM"
load(paste0("bwhalemodels_",tag,".RData"))
tag="reduce_space_1Chl"
load(paste0("bwhaledata_",tag,".RData"))
newdata.r<-raster("RasterExtent.grd")

dt<-"2012-09-01"

rasterdata<-readRDS(paste0("./RasterData/",dt,"/rasterdata.RDS"))
pdf(paste0("./RasterData/",dt,"/rasterdata.pdf"))
  plot(rasterdata)
dev.off()

rasterstack<-rasterdata
rasterdata<-as.data.frame(rasterstack)
rasterdata$eke<-rasterdata$EKE
rasterdata$RN<-0



predGAM(bwhaleGAMM_CRW$gam, savefilename= paste0("GAMplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_back$gam, savefilename= paste0("GAMplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_buff$gam, savefilename= paste0("GAMplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_rev$gam, savefilename= paste0("GAMplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

predBRT(bwhaleBRT.lr005.CRW, savefilename= paste0("BRTplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.rev, savefilename= paste0("BRTplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.buff, savefilename= paste0("BRTplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.back, savefilename= paste0("BRTplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

dt<-"2012-06-01"

rasterdata<-readRDS(paste0("./RasterData/",dt,"/rasterdata.RDS"))
pdf(paste0("./RasterData/",dt,"/rasterdata.pdf"))
plot(rasterdata)
dev.off()

#isthere<-which(names(rasterdata) %in% all.vars(formula(bwhaleGAMM_CRW$gam)[-2]))
#rasterstack<-raster::subset(rasterdata, isthere)
rasterstack<-rasterdata
rasterdata<-as.data.frame(rasterstack)
rasterdata$eke<-rasterdata$EKE
rasterdata$RN<-0

predGAM(bwhaleGAMM_CRW$gam, savefilename= paste0("GAMplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_back$gam, savefilename= paste0("GAMplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_buff$gam, savefilename= paste0("GAMplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_rev$gam, savefilename= paste0("GAMplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

predBRT(bwhaleBRT.lr005.CRW, savefilename= paste0("BRTplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.rev, savefilename= paste0("BRTplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.buff, savefilename= paste0("BRTplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.back, savefilename= paste0("BRTplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

dt<-"2012-03-01"

rasterdata<-readRDS(paste0("./RasterData/",dt,"/rasterdata.RDS"))
pdf(paste0("./RasterData/",dt,"/rasterdata.pdf"))
plot(rasterdata)
dev.off()

rasterstack<-rasterdata
rasterdata<-as.data.frame(rasterstack)
rasterdata$eke<-rasterdata$EKE
rasterdata$RN<-0

predGAM(bwhaleGAMM_CRW$gam, savefilename= paste0("GAMplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_back$gam, savefilename= paste0("GAMplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_buff$gam, savefilename= paste0("GAMplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_rev$gam, savefilename= paste0("GAMplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

predBRT(bwhaleBRT.lr005.CRW, savefilename= paste0("BRTplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.rev, savefilename= paste0("BRTplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.buff, savefilename= paste0("BRTplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.back, savefilename= paste0("BRTplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

dt<-"2012-12-01"

rasterdata<-readRDS(paste0("./RasterData/",dt,"/rasterdata.RDS"))
pdf(paste0("./RasterData/",dt,"/rasterdata.pdf"))
plot(rasterdata)
dev.off()

#isthere<-which(names(rasterdata) %in% all.vars(formula(bwhaleGAMM_CRW$gam)[-2]))
#rasterstack<-raster::subset(rasterdata, isthere)
rasterstack<-rasterdata
rasterdata<-as.data.frame(rasterstack)
rasterdata$eke<-rasterdata$EKE
rasterdata$RN<-0

predGAM(bwhaleGAMM_CRW$gam, savefilename= paste0("GAMplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_back$gam, savefilename= paste0("GAMplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_buff$gam, savefilename= paste0("GAMplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_rev$gam, savefilename= paste0("GAMplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

predBRT(bwhaleBRT.lr005.CRW, savefilename= paste0("BRTplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.rev, savefilename= paste0("BRTplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.buff, savefilename= paste0("BRTplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.back, savefilename= paste0("BRTplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

dt<-"2006-08-01"

rasterdata<-readRDS(paste0("./RasterData/",dt,"/rasterdata.RDS"))
pdf(paste0("./RasterData/",dt,"/rasterdata.pdf"))
plot(rasterdata)
dev.off()

rasterstack<-rasterdata
rasterdata<-as.data.frame(rasterstack)
rasterdata$eke<-rasterdata$EKE
rasterdata$RN<-0

predGAM(bwhaleGAMM_CRW$gam, savefilename= paste0("GAMplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_back$gam, savefilename= paste0("GAMplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_buff$gam, savefilename= paste0("GAMplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predGAM(bwhaleGAMM_rev$gam, savefilename= paste0("GAMplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

predBRT(bwhaleBRT.lr005.CRW, savefilename= paste0("BRTplot","_CRW_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.rev, savefilename= paste0("BRTplot","_rev_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.buff, savefilename= paste0("BRTplot","_buff_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)
predBRT(bwhaleBRT.lr005.back, savefilename= paste0("BRTplot","_back_",tag,"_",dt,".pdf"), rasterdata=rasterdata, dt = dt, makeplot=TRUE)

dt='2006-08-01'
rasterdata<-as.data.frame(rasterstack)
names(rasterdata)

alldata<-rbind(bwhaledata_back,bwhaledata_buff,bwhaledata_CRW,bwhaledata_rev)
allextent<-extent(c(min(bwhaledata_back$long),max(bwhaledata_back$long),min(bwhaledata_back$lat),max(bwhaledata_back$lat)))
## to predict on original data, run: rasterdata <- alldata

newdata<-predict.gam(bwhaleGAMM_CRW$gam, rasterdata, se.fit=TRUE, type="response")
newGLMM<-predict.gam(bwhaleGLMM_CRW$gam, rasterdata, se.fit=TRUE, type="response")
model<-bwhaleBRT.lr005.CRW
newdata_b<-predict.gbm(model, newdata=rasterdata,
                       n.trees=model$gbm.call$best.trees, type="response")
newdata<-cbind(rasterdata,newdata$fit,newGLMM$fit,newdata_b)
t<-cbind(newdata$lon,newdata$lat,newdata$`newdata$fit`,newdata$`newGLMM$fit`,newdata_b)
e <- extent(c(min(newdata$lon),max(newdata$lon),min(newdata$lat),max(newdata$lat)))
r <- raster(allextent, ncol=500, nrow=500)
x_CRW <- rasterize(t[,1:2], r, t[,3], fun=mean)
x_CRW_GLMM <- rasterize(t[,1:2], r, t[,4], fun=mean)
x_CRW_brt <- rasterize(t[,1:2], r, t[,5], fun=mean)

newdata<-predict.gam(bwhaleGAMM_rev$gam, rasterdata, se.fit=TRUE, type="response")
newGLMM<-predict.gam(bwhaleGLMM_rev$gam, rasterdata, se.fit=TRUE, type="response")
model<-bwhaleBRT.lr005.rev
newdata_b<-predict.gbm(model, newdata=rasterdata,
                       n.trees=model$gbm.call$best.trees, type="response")
newdata<-cbind(rasterdata,newdata$fit,newGLMM$fit,newdata_b)
t<-cbind(newdata$lon,newdata$lat,newdata$`newdata$fit`,newdata$`newGLMM$fit`,newdata_b)
e <- extent(c(min(newdata$lon),max(newdata$lon),min(newdata$lat),max(newdata$lat)))
r <- raster(allextent, ncol=500, nrow=500)
x_rev<- rasterize(t[,1:2], r, t[,3], fun=mean)
x_rev_GLMM <- rasterize(t[,1:2], r, t[,4], fun=mean)
x_rev_brt <- rasterize(t[,1:2], r, t[,5], fun=mean)

newdata<-predict.gam(bwhaleGAMM_back$gam, rasterdata, se.fit=TRUE, type="response")
newGLMM<-predict.gam(bwhaleGLMM_back$gam, rasterdata, se.fit=TRUE, type="response")
model<-bwhaleBRT.lr005.back
newdata_b<-predict.gbm(model, newdata=rasterdata,
                       n.trees=model$gbm.call$best.trees, type="response")
newdata<-cbind(rasterdata,newdata$fit,newGLMM$fit,newdata_b)
t<-cbind(newdata$lon,newdata$lat,newdata$`newdata$fit`,newdata$`newGLMM$fit`,newdata_b)
e <- extent(c(min(newdata$lon),max(newdata$lon),min(newdata$lat),max(newdata$lat)))
r <- raster(allextent, ncol=500, nrow=500)
x_back <- rasterize(t[,1:2], r, t[,3], fun=mean)
x_back_GLMM <- rasterize(t[,1:2], r, t[,4], fun=mean)
x_back_brt <- rasterize(t[,1:2], r, t[,5], fun=mean)

newdata<-predict.gam(bwhaleGAMM_buff$gam, rasterdata, se.fit=TRUE, type="response")
newGLMM<-predict.gam(bwhaleGLMM_buff$gam, rasterdata, se.fit=TRUE, type="response")
model<-bwhaleBRT.lr005.buff
newdata_b<-predict.gbm(model, newdata=rasterdata,
                       n.trees=model$gbm.call$best.trees, type="response")
newdata<-cbind(rasterdata,newdata$fit,newGLMM$fit,newdata_b)
t<-cbind(newdata$lon,newdata$lat,newdata$`newdata$fit`,newdata$`newGLMM$fit`,newdata_b)
e <- extent(c(min(newdata$lon),max(newdata$lon),min(newdata$lat),max(newdata$lat)))
r <- raster(allextent, ncol=500, nrow=500)
x_buff <- rasterize(t[,1:2], r, t[,3], fun=mean)
x_buff_GLMM <- rasterize(t[,1:2], r, t[,4], fun=mean)
x_buff_brt <- rasterize(t[,1:2], r, t[,5], fun=mean)

save(x_CRW,x_CRW_brt,x_CRW_GLMM,x_rev,x_rev_brt,x_rev_GLMM,x_back,x_back_brt,x_back_GLMM,x_buff,x_buff_brt,x_buff_GLMM, file=paste0("bwhalePredict_",Sys.Date(),".RData"))

pdf(file=paste0("bwhalePredict_",Sys.Date(),".pdf"))
par(mfrow=c(4,3))
plot(x_CRW)
plot(x_CRW_brt)
plot(x_CRW_GLMM)
plot(x_rev)
plot(x_rev_brt)
plot(x_rev_GLMM)
plot(x_back)
plot(x_back_brt)
plot(x_back_GLMM)
plot(x_buff)
plot(x_buff_brt)
plot(x_buff_GLMM)
dev.off()





#### Evaluate LOO by month
tag="final"

bwhaledata_buff <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_buff_",tag,".RDS")))
bwhaledata_back <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_back_",tag,".RDS")))
bwhaledata_CRW <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_CRW_",tag,".RDS")))
bwhaledata_rev <- (readRDS(paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_rev_",tag,".RDS")))

load(paste0("bwhalemodels_",tag,".RData"))

formula=formula(bwhaleGAMM_CRW$gam)
bwhaleGAMM_CRW.LOOm<-eval_LOOm_GAMM(bwhaledata_CRW, formula, months=c(1:12))
formula=formula(bwhaleGAMM_back$gam)
bwhaleGAMM_back.LOOm<-eval_LOOm_GAMM(bwhaledata_back, formula, months=c(1:12))
formula=formula(bwhaleGAMM_buff$gam)
bwhaleGAMM_buff.LOOm<-eval_LOOm_GAMM(bwhaledata_buff, formula, months=c(1:12))
formula=formula(bwhaleGAMM_rev$gam)
bwhaleGAMM_rev.LOOm<-eval_LOOm_GAMM(bwhaledata_rev, formula, months=c(1:12))

bwhaleBRT_CRW.LOOm<-eval_LOOm_BRT(BRTtransformDataFrame_BW(bwhaledata_CRW), gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, lr=0.005, months=c(1:12))
bwhaleBRT_back.LOOm<-eval_LOOm_BRT(BRTtransformDataFrame_BW(bwhaledata_back), gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, lr=0.005, months=c(1:12))
bwhaleBRT_buff.LOOm<-eval_LOOm_BRT(BRTtransformDataFrame_BW(bwhaledata_buff), gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, lr=0.005, months=c(1:12))
bwhaleBRT_rev.LOOm<-eval_LOOm_BRT(BRTtransformDataFrame_BW(bwhaledata_rev), gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, lr=0.005, months=c(1:12))

bwhaleGAMM_CRW.LOOm$type<-paste0(tag,"_CRW")
bwhaleGAMM_back.LOOm$type<-paste0(tag,"_back")
bwhaleGAMM_buff.LOOm$type<-paste0(tag,"_buff")
bwhaleGAMM_rev.LOOm$type<-paste0(tag,"_rev")
write.csv(rbind(bwhaleGAMM_CRW.LOOm,bwhaleGAMM_back.LOOm,bwhaleGAMM_buff.LOOm,bwhaleGAMM_rev.LOOm), file=paste0("GAMM_LOOm_",tag,".csv"))


bwhaleBRT_CRW.LOOm$type<-paste0(tag,"_CRW")
bwhaleBRT_back.LOOm$type<-paste0(tag,"_back")
bwhaleBRT_buff.LOOm$type<-paste0(tag,"_buff")
bwhaleBRT_rev.LOOm$type<-paste0(tag,"_rev")
write.csv(rbind(bwhaleBRT_CRW.LOOm,bwhaleBRT_back.LOOm,bwhaleBRT_buff.LOOm,bwhaleBRT_rev.LOOm), file=paste0("BRT_LOOm_",tag,".csv"))


#### Evaluate LOO by space, for whales shallower than 400 vs deeper than 400

load(paste0("bwhalemodels_",tag,".RData"))

traindata<-bwhaledata_buff[which(bwhaledata_buff$Bathymetry<(-400)),]
testdata<-bwhaledata_buff[which(bwhaledata_buff$Bathymetry>(-400)),]
formula=formula(bwhaleGAMM_buff$gam)
bwhaleGAMM_buff.LOOsp<-eval_LOOspace_GAMM(bwhaledata_buff, formula, traindata=traindata, testdata=testdata)
bwhaleBRT_buff.LOOsp<-eval_LOOspace_BRT(BRTtransformDataFrame_BW(bwhaledata_buff), gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, lr=0.005, traindata=BRTtransformDataFrame_BW(traindata), testdata=BRTtransformDataFrame_BW(testdata))

traindata<-bwhaledata_CRW[which(bwhaledata_CRW$Bathymetry<(-400)),]
testdata<-bwhaledata_CRW[which(bwhaledata_CRW$Bathymetry>(-400)),]
formula=formula(bwhaleGAMM_CRW$gam)
bwhaleGAMM_CRW.LOOsp<-eval_LOOspace_GAMM(bwhaledata_CRW, formula, traindata=traindata, testdata=testdata)
bwhaleBRT_CRW.LOOsp<-eval_LOOspace_BRT(BRTtransformDataFrame_BW(bwhaledata_CRW), gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, lr=0.005, traindata=BRTtransformDataFrame_BW(traindata), testdata=BRTtransformDataFrame_BW(testdata))

traindata<-bwhaledata_back[which(bwhaledata_back$Bathymetry<(-400)),]
testdata<-bwhaledata_back[which(bwhaledata_back$Bathymetry>(-400)),]
formula=formula(bwhaleGAMM_back$gam)
bwhaleGAMM_back.LOOsp<-eval_LOOspace_GAMM(bwhaledata_back, formula, traindata=traindata, testdata=testdata)
bwhaleBRT_back.LOOsp<-eval_LOOspace_BRT(BRTtransformDataFrame_BW(bwhaledata_back), gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, lr=0.005, traindata=BRTtransformDataFrame_BW(traindata), testdata=BRTtransformDataFrame_BW(testdata))

traindata<-bwhaledata_rev[which(bwhaledata_rev$Bathymetry<(-400)),]
testdata<-bwhaledata_rev[which(bwhaledata_rev$Bathymetry>(-400)),]
formula=formula(bwhaleGAMM_rev$gam)
bwhaleGAMM_rev.LOOsp<-eval_LOOspace_GAMM(bwhaledata_rev, formula, traindata=traindata, testdata=testdata)
bwhaleBRT_rev.LOOsp<-eval_LOOspace_BRT(BRTtransformDataFrame_BW(bwhaledata_rev), gbm.x=c(9,12,14,17,18,21,22), gbm.y=1, lr=0.005, traindata=BRTtransformDataFrame_BW(traindata), testdata=BRTtransformDataFrame_BW(testdata))

bwhaleGAMM_CRW.LOOsp$type<-paste0(tag,"_CRW")
bwhaleGAMM_back.LOOsp$type<-paste0(tag,"_back")
bwhaleGAMM_buff.LOOsp$type<-paste0(tag,"_buff")
bwhaleGAMM_rev.LOOsp$type<-paste0(tag,"_rev")
write.csv(rbind(bwhaleGAMM_CRW.LOOsp,bwhaleGAMM_back.LOOsp,bwhaleGAMM_buff.LOOsp,bwhaleGAMM_rev.LOOsp), file=paste0("GAMM_LOOsp_",tag,".csv"))


bwhaleBRT_CRW.LOOsp$type<-paste0(tag,"_CRW")
bwhaleBRT_back.LOOsp$type<-paste0(tag,"_back")
bwhaleBRT_buff.LOOsp$type<-paste0(tag,"_buff")
bwhaleBRT_rev.LOOsp$type<-paste0(tag,"_rev")
#write.csv(rbind(bwhaleGAMM_CRW.LOOm,bwhaleGAMM_back.LOOm,bwhaleGAMM_buff.LOOm,bwhaleGAMM_rev.LOOm), file=paste0("GAMM_LOOm_",tag,".csv"))
write.csv(rbind(bwhaleBRT_CRW.LOOsp,bwhaleBRT_back.LOOsp,bwhaleBRT_buff.LOOsp,bwhaleBRT_rev.LOOsp), file=paste0("BRT_LOOsp_",tag,".csv"))



#### Plot elephant background & then evaluate LOO by space for elephants

tag="final"

ephantdata_buff <- (readRDS(paste0("./Data/ephantdata_buff_",tag,".RDS")))
ephantdata_back <- (readRDS(paste0("./Data/ephantdata_back_",tag,".RDS")))
ephantdata_CRW <- (readRDS(paste0("./Data/ephantdata_CRW_",tag,".RDS")))
ephantdata_rev <- (readRDS(paste0("./Data/ephantdata_rev_",tag,".RDS")))

ephantdata_buff$dTime <- fixdTime(ephantdata_buff$dTime)
ephantdata_back$dTime <- fixdTime(ephantdata_back$dTime)
ephantdata_CRW$dTime <- fixdTime(ephantdata_CRW$dTime)
ephantdata_rev$dTime <- fixdTime(ephantdata_rev$dTime)

ephantdata_buff[ephantdata_buff$iteration==0,]<-UTM2WGS(ephantdata_buff[ephantdata_buff$iteration==0,])
ephantdata_back[ephantdata_back$iteration==0,]<-UTM2WGS(ephantdata_back[ephantdata_back$iteration==0,])
ephantdata_CRW[ephantdata_CRW$iteration==0,]<-UTM2WGS(ephantdata_CRW[ephantdata_CRW$iteration==0,])
ephantdata_rev[ephantdata_rev$iteration==0,]<-UTM2WGS(ephantdata_rev[ephantdata_rev$iteration==0,])

tag="etosha_only"
load(paste0("ephantmodels",tag,".RData"))
load(paste0("ephantdata",tag,".RData"))

rasterstack<-readRDS(paste0("./RasterData/ephant/ephant_rasterdata.RDS"))
pdf(paste0("./RasterData/ephant/ephant_rasterdata.pdf"))
  plot(rasterstack)
dev.off()


dt='ephant'
rasterdata<-as.data.frame(rasterstack)
names(rasterdata)

alldata<-rbind(ephantdata_back,ephantdata_buff,ephantdata_CRW,ephantdata_rev)
allextent<-extent(c(min(ephantdata_back$long),max(ephantdata_back$long),min(ephantdata_back$lat),max(ephantdata_back$lat)))
## to predict on original data, run: rasterdata <- alldata

newdata<-predict.gam(ephantGAMM_CRW$gam, rasterdata, se.fit=TRUE, type="response")
newGLMM<-predict.gam(ephantGLMM_CRW$gam, rasterdata, se.fit=TRUE, type="response")
model<-ephantBRT.lr005.CRW
newdata_b<-predict.gbm(model, newdata=rasterdata,
                  n.trees=model$gbm.call$best.trees, type="response")
newdata<-cbind(rasterdata,newdata$fit,newGLMM$fit,newdata_b)
t<-cbind(newdata$long,newdata$lat,newdata$`newdata$fit`,newdata$`newGLMM$fit`,newdata_b)
e <- extent(c(min(newdata$long),max(newdata$long),min(newdata$lat),max(newdata$lat)))
r <- raster(allextent, ncol=500, nrow=500)
x_CRW <- rasterize(t[,1:2], r, t[,3], fun=mean)
x_CRW_GLMM <- rasterize(t[,1:2], r, t[,4], fun=mean)
x_CRW_brt <- rasterize(t[,1:2], r, t[,5], fun=mean)

newdata<-predict.gam(ephantGAMM_rev$gam, rasterdata, se.fit=TRUE, type="response")
newGLMM<-predict.gam(ephantGLMM_rev$gam, rasterdata, se.fit=TRUE, type="response")
model<-ephantBRT.lr005.rev
newdata_b<-predict.gbm(model, newdata=rasterdata,
                       n.trees=model$gbm.call$best.trees, type="response")
newdata<-cbind(rasterdata,newdata$fit,newGLMM$fit,newdata_b)
t<-cbind(newdata$long,newdata$lat,newdata$`newdata$fit`,newdata$`newGLMM$fit`,newdata_b)
e <- extent(c(min(newdata$long),max(newdata$long),min(newdata$lat),max(newdata$lat)))
r <- raster(allextent, ncol=500, nrow=500)
x_rev<- rasterize(t[,1:2], r, t[,3], fun=mean)
x_rev_GLMM <- rasterize(t[,1:2], r, t[,4], fun=mean)
x_rev_brt <- rasterize(t[,1:2], r, t[,5], fun=mean)

newdata<-predict.gam(ephantGAMM_back$gam, rasterdata, se.fit=TRUE, type="response")
newGLMM<-predict.gam(ephantGLMM_back$gam, rasterdata, se.fit=TRUE, type="response")
model<-ephantBRT.lr005.back
newdata_b<-predict.gbm(model, newdata=rasterdata,
                       n.trees=model$gbm.call$best.trees, type="response")
newdata<-cbind(rasterdata,newdata$fit,newGLMM$fit,newdata_b)
t<-cbind(newdata$long,newdata$lat,newdata$`newdata$fit`,newdata$`newGLMM$fit`,newdata_b)
e <- extent(c(min(newdata$long),max(newdata$long),min(newdata$lat),max(newdata$lat)))
r <- raster(allextent, ncol=500, nrow=500)
x_back <- rasterize(t[,1:2], r, t[,3], fun=mean)
x_back_GLMM <- rasterize(t[,1:2], r, t[,4], fun=mean)
x_back_brt <- rasterize(t[,1:2], r, t[,5], fun=mean)

newdata<-predict.gam(ephantGAMM_buff$gam, rasterdata, se.fit=TRUE, type="response")
newGLMM<-predict.gam(ephantGLMM_buff$gam, rasterdata, se.fit=TRUE, type="response")
model<-ephantBRT.lr005.buff
newdata_b<-predict.gbm(model, newdata=rasterdata,
                       n.trees=model$gbm.call$best.trees, type="response")
newdata<-cbind(rasterdata,newdata$fit,newGLMM$fit,newdata_b)
t<-cbind(newdata$long,newdata$lat,newdata$`newdata$fit`,newdata$`newGLMM$fit`,newdata_b)
e <- extent(c(min(newdata$long),max(newdata$long),min(newdata$lat),max(newdata$lat)))
r <- raster(allextent, ncol=500, nrow=500)
x_buff <- rasterize(t[,1:2], r, t[,3], fun=mean)
x_buff_GLMM <- rasterize(t[,1:2], r, t[,4], fun=mean)
x_buff_brt <- rasterize(t[,1:2], r, t[,5], fun=mean)

save(x_CRW,x_CRW_brt,x_CRW_GLMM,x_rev,x_rev_brt,x_rev_GLMM,x_back,x_back_brt,x_back_GLMM,x_buff,x_buff_brt,x_buff_GLMM, file=paste0("EphantPredict_",Sys.Date(),".RData"))

pdf(file=paste0("EphantPredict_",Sys.Date(),".pdf"))
  par(mfrow=c(4,3))
    plot(x_CRW)
    plot(x_CRW_brt)
    plot(x_CRW_GLMM)
    plot(x_rev)
    plot(x_rev_brt)
    plot(x_rev_GLMM)
    plot(x_back)
    plot(x_back_brt)
    plot(x_back_GLMM)
    plot(x_buff)
    plot(x_buff_brt)
    plot(x_buff_GLMM)
  dev.off()

### Eval LOO by space
trainextent <- c(221252.9,602480.1,3922484,7897285)
### in GPS now
trainextent <- c(15.5,16.1,-18.45,-19.25)

testdata <- ephantdata_buff[which(ephantdata_buff$lat<trainextent[4]),]
traindata <- ephantdata_buff[which(ephantdata_buff$lat>trainextent[4]),]
formula=formula(ephantGAMM_buff$gam)
ephantGAMM_buff.LOOsp<-eval_LOOspace_GAMM(ephantdata_buff, formula, traindata=traindata, testdata=testdata)
ephantBRT_buff.LOOsp<-eval_LOOspace_BRT(ephantdata_buff, gbm.x=c(6,7,8), gbm.y=9, lr=0.005, traindata=traindata, testdata=testdata)

traindata<-ephantdata_CRW[which(ephantdata_CRW$lat<trainextent[4]),]
testdata<-ephantdata_CRW[which(ephantdata_CRW$lat>trainextent[4]),]
formula=formula(ephantGAMM_CRW$gam)
ephantGAMM_CRW.LOOsp<-eval_LOOspace_GAMM(ephantdata_CRW, formula, traindata=traindata, testdata=testdata)
ephantBRT_CRW.LOOsp<-eval_LOOspace_BRT((ephantdata_CRW), gbm.x=c(6,7,8), gbm.y=9, lr=0.005, traindata=(traindata), testdata=(testdata))

traindata<-ephantdata_back[which(ephantdata_back$lat<trainextent[4]),]
testdata<-ephantdata_back[which(ephantdata_back$lat>trainextent[4]),]
formula=formula(ephantGAMM_back$gam)
ephantGAMM_back.LOOsp<-eval_LOOspace_GAMM(ephantdata_back, formula, traindata=traindata, testdata=testdata)
ephantBRT_back.LOOsp<-eval_LOOspace_BRT((ephantdata_back), gbm.x=c(6,7,8), gbm.y=9, lr=0.005, traindata=(traindata), testdata=(testdata))

traindata<-ephantdata_rev[which(ephantdata_rev$lat<trainextent[4]),]
testdata<-ephantdata_rev[which(ephantdata_rev$lat>trainextent[4]),]
formula=formula(ephantGAMM_rev$gam)
ephantGAMM_rev.LOOsp<-eval_LOOspace_GAMM(ephantdata_rev, formula, traindata=traindata, testdata=testdata)
ephantBRT_rev.LOOsp<-eval_LOOspace_BRT((ephantdata_rev), gbm.x=c(6,7,8), gbm.y=9, lr=0.005, traindata=(traindata), testdata=(testdata))

ephantGAMM_CRW.LOOsp$type<-paste0(tag,"_CRW")
ephantGAMM_back.LOOsp$type<-paste0(tag,"_back")
ephantGAMM_buff.LOOsp$type<-paste0(tag,"_buff")
ephantGAMM_rev.LOOsp$type<-paste0(tag,"_rev")
write.csv(rbind(ephantGAMM_CRW.LOOsp,ephantGAMM_back.LOOsp,ephantGAMM_buff.LOOsp,ephantGAMM_rev.LOOsp), file=paste0("GAMM_ephant_LOOsp_",tag,".csv"))


ephantBRT_CRW.LOOsp$type<-paste0(tag,"_CRW")
ephantBRT_back.LOOsp$type<-paste0(tag,"_back")
ephantBRT_buff.LOOsp$type<-paste0(tag,"_buff")
ephantBRT_rev.LOOsp$type<-paste0(tag,"_rev")
write.csv(rbind(ephantBRT_CRW.LOOsp,ephantBRT_back.LOOsp,ephantBRT_buff.LOOsp,ephantBRT_rev.LOOsp), file=paste0("BRT_ephant_LOOsp_",tag,".csv"))

