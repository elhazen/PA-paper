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






