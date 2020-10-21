#Data extraction Function
#Function to extract satellite data for the pseudo-absence paper, intitially for Blue whales. 
#Extraction resolution set at 0.25 for blue whales. 
#Written by Steph


Build_Raster_Stack <- function(dt, rasterextent=newextent, dirname='', savePieces=T){
  #------load packages and library------
  library("ncdf4")
  library("raster")
  library("tidyverse")
  library("lubridate")
  
  outDir=paste0("./RasterData/",dirname,as.character(as.Date(dt)))
  bathyDir="./RasterData/"
  dir.create(file.path(outDir), showWarnings = FALSE)
  desired.resolution = 0.25
  desired.resolution = desired.resolution/2
  sd.resolution = 1
  mth <- month(dt)
  yr <- year(dt)
  
  print(dt)
  
  #----BEGIN EXTRACTION----
  #Code differs for each satellite product, all downloaded from CMEMS
  #For more details on each product see Steph who has a product summary table
  
  if (file.exists(paste0(bathyDir,"bathymetry.RDS"))){
    b<-readRDS(paste0(bathyDir,"bathymetry.RDS"))
    b_sd<-readRDS(paste0(bathyDir,"rugosity.RDS"))
    lat.r<-readRDS(paste0(bathyDir,"lat.RDS"))
    lon.r<-readRDS(paste0(bathyDir,"lon.RDS"))
    print('LOADING EXISTING BATHYMETRY')    
  } else {
    
    #----Bathymetry & Rugosity----
    print('RASTER STACK BATHYMETRY')
    #data["Bathymetry"] <- NA
    #data["Rugosity"] <- NA
    dir <- "/Volumes/HazWht_8TB/Data/GlobalData/"
    nc <- paste0(dir,'ETOPO180_Global/etopo180_4558_b013_12d8.nc')
    nc.data <- nc_open(nc, write=FALSE)
    lat <- ncvar_get(nc.data,'latitude')
    lon <- ncvar_get(nc.data,'longitude')
    nrows <- length(lat); ncols <- length(lon)
    nc_close(nc.data)
    
    b<-raster(nc)
    n.sd.cells<-(sd.resolution/res(b))/2-1
    
    b2<-b
    
    rtemp<-rasterextent
    extent(rtemp)<-extent(rasterextent) + c(-5,5,-5,5)
    res(rtemp) <- res(b)
    b2 <- resample(b,rtemp, method="bilinear")
    
    #extent(b2) <- extent(rasterextent) + c(-5,5,-5,5)
    #res(b2) <- res(b)
    n.sd.cells<-(sd.resolution/res(b2))/2-1
    
    #b_sd<-focal(b2, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
    
    beginCluster()
    b_sd <- clusterR(b2, focal, args=list(w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=TRUE, pad=T))
    endCluster()
    
    b<-resample(b,rasterextent, method="bilinear")
    b_sd<-resample(b_sd,rasterextent, method="bilinear")
    
    b_sd[b > 0] <- NA
    b[b > 0] <- NA

    b.pts <- rasterToPoints(b, spatial=TRUE)
    b.pts@data <- data.frame(b.pts@data, lon=coordinates(b.pts)[,1],
                             lat=coordinates(b.pts)[,2])
    lat.r <- init(b,'y')
    lon.r <- init(b,'x')
    
    
    #head(b.pts@data)
        
    if (savePieces){
      saveRDS(b,paste(bathyDir,"bathymetry.RDS",sep=""))
      saveRDS(b_sd,paste(bathyDir,"rugosity.RDS",sep=""))
      saveRDS(lat.r,paste(bathyDir,"lat.RDS",sep=""))
      saveRDS(lon.r,paste(bathyDir,"lon.RDS",sep=""))
    }
  }
  
  if (file.exists(paste0(outDir,"/oxygen_100.RDS"))){
    o100<-readRDS(paste0(outDir,"/oxygen_100.RDS"))
    osurf<-readRDS(paste0(outDir,"/oxygen_surf.RDS"))
    print('LOADING EXISTING OXYGEN')    
  } else {
    print('RASTER STACK OXYGEN')
    #1 degree resolution so just extraction point value
    #data["Oxygen_Surface"] <- NA
    #data["Oxygen_100m"] <- NA
    dir <- "/Volumes/HazWht_8TB/Data/GlobalData/WOA 2013 oxygen monthly clim/"
    if (nchar(month(dt))==1) mth <- paste0("0",month(dt))
    nc <- paste0(dir,'woa13_all_o',mth,'_01.nc')
    
    nc.data <- nc_open(nc, write=FALSE)
    o100 <- ncvar_get(nc.data,varid="o_an")
    osurf<-o100[,,1]
    o100<-o100[,,21]
    #o100=brick(nc,level=21,varname="o_an") %>% raster()
    dims <- nc.data$dim
    nc_close(nc.data)
    e <- extent(c(range(dims$lon$vals), range(dims$lat$vals)))
    o100 <- raster((o100))
    extent(o100) <- e
    osurf <- raster((osurf))
    extent(osurf) <- e
    o100<-flip(t(o100),direction='y')
    osurf<-flip(t(osurf),direction='y')
    extent(o100) <- e
    extent(osurf) <- e
    
    
    #### Need to get the 21st layer here!!!
    #o100 <- stack(nc, varname="o_an")
    #o100 <- o100[[,,21]] #raster(nc, varname="o_an")   # need to grab the 21st layer
    #n.sd.cells<-(sd.resolution/res(b))/2-1
    #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
    o100<-resample(o100,rasterextent, method="bilinear")
    osurf<-resample(osurf,rasterextent, method="bilinear")
    
    if (savePieces){
      saveRDS(o100,paste(outDir,"oxygen_100.RDS",sep="/"))
      saveRDS(osurf,paste(outDir,"oxygen_surf.RDS",sep="/"))
    }
  }   
  #----SSH, SLA, U, and V----
  if (file.exists(paste0(outDir,"/eke.RDS"))){
    ssh.r<-readRDS(paste(outDir,"SSHa.RDS", sep="/"))
    ssh.sd.r<-readRDS(paste(outDir,"SSHa_sd.RDS", sep="/"))
    sla.r<-readRDS(paste(outDir,"SLA.RDS", sep="/"))
    U.r<-readRDS(paste(outDir,"U.RDS", sep="/"))
    V.r<-readRDS(paste(outDir,"V.RDS",sep="/"))
    eke.r<-readRDS(paste(outDir,"eke.RDS", sep="/"))
    print('LOADING EXISTING SSH PRODUCTS')    
  } else {
    print('RASTER STACK SSH')
    dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_SSH_uv/"
    
    if (yr >=1998 & yr < 2018){
      
      files <- list.files(paste0(dir,yr))
      parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
      dtname <- as.Date(parsename[,6], format="%Y%m%d")
      file_index <- which.min(abs(as.Date(dt) - dtname))
      
      nc <- paste0(dir,yr,'/',files[file_index])
      ssh.r <- raster(nc, varname="sla")
      sla.r <- raster(nc, varname="adt")
      U.r <- raster(nc, varname="ugos")
      V.r <- raster(nc, varname="vgos")
      eke.r <- (U.r^2+V.r^2)/2
      
      #   n.sd.cells<-(sd.resolution/res(ssh.r))/2-1
      #   ssh.r2<-ssh.r
      #   extent(ssh.r2) <- extent(rasterextent) + c(-5,5,-5,5)
      
      #ssh.sd.r <- focal(ssh.r2, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
      rtemp<-rasterextent
      extent(rtemp)<-extent(rasterextent) + c(-5,5,-5,5)
      res(rtemp) <- res(ssh.r)
      ssh.r<-rotate(ssh.r)
      ssh.r2 <- resample(ssh.r,rtemp, method="bilinear")
      
      n.sd.cells<-(sd.resolution/res(ssh.r2))-1#/2-1
      
      beginCluster()
      ssh.sd.r <- clusterR(ssh.r2, focal, args=list(w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=TRUE, pad=T))
      endCluster()
      
      ssh.r<-resample(ssh.r,rasterextent, method="bilinear")
      sla.r <- sla.r %>% rotate() %>% resample(., rasterextent, method="bilinear")
      U.r <- U.r %>% rotate() %>% resample(., rasterextent, method="bilinear")
      V.r <- V.r %>% rotate() %>% resample(., rasterextent, method="bilinear")
      eke.r <- eke.r %>% rotate() %>% resample(., rasterextent, method="bilinear")
      ssh.sd.r <- ssh.sd.r %>% resample(., rasterextent, method="bilinear")
      
    } 
    
    if (savePieces){
      saveRDS(ssh.r,paste(outDir,"SSHa.RDS",sep="/"))
      saveRDS(ssh.sd.r,paste(outDir,"SSHa_sd.RDS",sep="/"))
      saveRDS(sla.r,paste(outDir,"SLA.RDS",sep="/"))
      saveRDS(U.r,paste(outDir,"U.RDS",sep="/"))
      saveRDS(V.r,paste(outDir,"V.RDS",sep="/"))
      saveRDS(eke.r,paste(outDir,"eke.RDS",sep="/"))
    }
  }  
  
  #---- LAGGED Chla historical and NRT: 4km and 8day----
  dtorig<-dt
  outdir1mo<-paste0("./RasterData/",as.character(as.Date(as.Date(dt) %m-% months(1))))
  if (file.exists(paste0(outDir,"/Chla-4km-lag1mo.RDS"))){
    chl_25k_lag.r<-readRDS(paste(outDir,"Chla-25km-lag1mo.RDS",sep="/"))
    chl_4k_lag.r<-readRDS(paste(outDir,"Chla-4km-lag1mo.RDS", sep="/"))
    print('LOADING EXISTING LAGGED CHL-A PRODUCTS')    
  } else if (file.exists(paste0(outdir1mo,"/Chla-25km.RDS"))) {
    chl_25k.r<-readRDS(paste(outdir1mo,"Chla-25km.RDS",sep="/"))
    chl_4k.r<-readRDS(paste(outdir1mo,"Chla-4km.RDS", sep="/"))
    print('LOADING LAST MONTH CHL-A PRODUCTS')    
  } else {
    print('EXTRACTING CHLOROPHYLL 4KM - 1 MO lagged')
    #Chlorophyll-a at 4km over 8day
    if (yr>=1998){
      dt <- as.Date(dt) %m-% months(1)
      if (dt < "2017-01-01"){
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_Chla_4km_8day_Historical/"
        files <- list.files(paste0(dir,yr),pattern=".nc")
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        nc <- paste0(dir,yr,'/',files[file_index])
        #data$Chla_4km_8day[i] <- mean(data.var1,na.rm=T)
        chl_4k_lag.r<-raster(nc)
        #n.sd.cells<-(sd.resolution/res(b))/2-1
        #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        chl_4k_lag.r<-resample(chl_4k_lag.r,rasterextent, method="bilinear")
        
      }
      
      if (dt >= "2017-01-01" & dt<"2018-01-01"){
        nc <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_Chla_4km_8day_NRT/dataset-oc-glo-chl-multi-l4-gsm_4km_8days-rt-v02_1520637318265.nc"
        chl_4k_lag.r<-raster(nc)
        #n.sd.cells<-(sd.resolution/res(b))/2-1
        #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        chl_4k_lag.r<-resample(chl_4k_lag.r,rasterextent, method="bilinear")
      }
    }
    if (savePieces){
      saveRDS(chl_4k_lag.r,paste(outDir,"Chla-4km-lag1mo.RDS",sep="/"))
    }
    
    
    print('EXTRACTING CHLOROPHYLL 25KM - 1 MO lagged')
    if (yr>=1998){
      if (dt < "2017-01-01"){
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_Chla_25km_monthly_Historical/"
        files <- list.files(paste0(dir,yr),pattern=".nc")
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        nc <- paste0(dir,yr,'/',files[file_index])
        nc.data <- nc_open(nc, write=FALSE)
        chl_25k_lag.r<-raster(nc)
        #n.sd.cells<-(sd.resolution/res(b))/2-1
        #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        chl_25k_lag.r<-resample(chl_25k_lag.r,rasterextent, method="bilinear")
      } 
      
      if (dt >= "2017-01-01" & dt < "2018-01-01"){
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_Chla_25km_monthly_NRT/"
        files <- list.files(paste0(dir,yr),pattern=".nc")
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        nc <- paste0(dir,yr,'/',files[file_index])
        chl_25k_lag.r<-raster(nc)
        #n.sd.cells<-(sd.resolution/res(b))/2-1
        #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        chl_25k_lag.r<-resample(chl_25k_lag.r,rasterextent, method="bilinear")
      }
    }
    if (savePieces){
      saveRDS(chl_25k_lag.r,paste(outDir,"Chla-25km-lag1mo.RDS",sep="/"))
    }
  } #else
  
  dt<-dtorig
  
  #----Chla historical and NRT: 4km and 8day----
  if (file.exists(paste0(outDir,"/Chla-25km.RDS"))){
    chl_25k.r<-readRDS(paste(outDir,"Chla-25km.RDS",sep="/"))
    chl_4k.r<-readRDS(paste(outDir,"Chla-4km.RDS", sep="/"))
    print('LOADING EXISTING CHL-A PRODUCTS')    
  } else {
    print('EXTRACTING CHLOROPHYLL 4KM')
    #Chlorophyll-a at 4km over 8day
    if (yr>=1998){
      if (dt < "2017-01-01"){
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_Chla_4km_8day_Historical/"
        files <- list.files(paste0(dir,yr),pattern=".nc")
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        nc <- paste0(dir,yr,'/',files[file_index])
        #data$Chla_4km_8day[i] <- mean(data.var1,na.rm=T)
        chl_4k.r<-raster(nc)
        #n.sd.cells<-(sd.resolution/res(b))/2-1
        #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        chl_4k.r<-resample(chl_4k.r,rasterextent, method="bilinear")
        
      }
      
      if (dt >= "2017-01-01" & dt<"2018-01-01"){
        nc <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_Chla_4km_8day_NRT/dataset-oc-glo-chl-multi-l4-gsm_4km_8days-rt-v02_1520637318265.nc"
        chl_4k.r<-raster(nc)
        #n.sd.cells<-(sd.resolution/res(b))/2-1
        #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        chl_4k.r<-resample(chl_4k.r,rasterextent, method="bilinear")
      }
    }
    if (savePieces){
      saveRDS(chl_4k.r,paste(outDir,"Chla-4km.RDS",sep="/"))
    }
    
    print('EXTRACTING CHLOROPHYLL 25KM')
    if (yr>=1998){
      if (dt < "2017-01-01"){
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_Chla_25km_monthly_Historical/"
        files <- list.files(paste0(dir,yr),pattern=".nc")
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        nc <- paste0(dir,yr,'/',files[file_index])
        nc.data <- nc_open(nc, write=FALSE)
        chl_25k.r<-raster(nc)
        #n.sd.cells<-(sd.resolution/res(b))/2-1
        #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        chl_25k.r<-resample(chl_25k.r,rasterextent, method="bilinear")
      } 
      
      if (dt >= "2017-01-01" & dt < "2018-01-01"){
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_Chla_25km_monthly_NRT/"
        files <- list.files(paste0(dir,yr),pattern=".nc")
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        nc <- paste0(dir,yr,'/',files[file_index])
        chl_25k.r<-raster(nc)
        #n.sd.cells<-(sd.resolution/res(b))/2-1
        #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        chl_25k.r<-resample(chl_25k.r,rasterextent, method="bilinear")
      }
    }
    if (savePieces){
      saveRDS(chl_25k.r,paste(outDir,"Chla-25km.RDS",sep="/"))
    }
  } #else
  
  
  #----AVISO FSLE----
  if (file.exists(paste0(outDir,"/FSLE.RDS"))){
    FSLE.r<-readRDS(paste(outDir,"FSLE.RDS",sep="/"))
    Theta.r<-readRDS(paste(outDir,"Theta.RDS", sep="/"))
    print('LOADING EXISTING FSLE PRODUCTS')    
  } else {
    print('EXTRACTING FSLE')
    #Native resolution is 0.25 so extracting point value
    #Looks like I have a corrupt file (dt_global_allsat_madt_fsle_20040908_20141027.nc) and have code below to skip it
    #data["FSLE_max"] <- NA
    #data["Theta_max"] <- NA
    dir <- "/Volumes/HazWht_8TB/Data/GlobalData/AVISO_FSLE/"
    
    if (yr>=1998 & yr < 2018){
      #Figure out which file to grab
      files <- list.files(paste0(dir,yr))
      parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
      dtname <- as.Date(parsename[,6], format="%Y%m%d")
      file_index <- which.min(abs(as.Date(dt) - dtname))
      nc <- paste0(dir,yr,'/',files[file_index])
      if (files[file_index] %in% 'dt_global_allsat_madt_fsle_20040908_20141027.nc'){
      } else { #Skip corrupt file
        #nc.data <- nc_open(nc, write=FALSE)
        #data.var1  <-  ncvar_get(nc.data,'fsle_max',start=c(c,r,length(time)), #time is indexed this way as I want it to throw an error if there is >1 time field
        #                         count=c(1,1,1),verbose=FALSE)
        #data.var2  <-  ncvar_get(nc.data,'theta_max',start=c(c,r,length(time)),
        #                         count=c(1,1,1),verbose=FALSE)
        FSLE.r<-(raster(nc, varname="fsle_max"))
        Theta.r<-(raster(nc, varname="theta_max"))
        
        FSLE.r <- rotate(FSLE.r)
        Theta.r <- rotate(Theta.r)
        #n.sd.cells<-(sd.resolution/res(b))/2-1
        #b_sd<-focal(b, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        FSLE.r<-resample(FSLE.r,rasterextent, method="bilinear")
        Theta.r<-resample(Theta.r,rasterextent, method="bilinear")
        #nc_close(nc)
      }
    } 
    
    if (savePieces){
      saveRDS(FSLE.r,paste(outDir,"FSLE.RDS",sep="/"))
      saveRDS(Theta.r,paste(outDir,"Theta.RDS",sep="/"))
    }
  } #else
  
  #----SST----
  if (file.exists(paste0(outDir,"/SST_SD.RDS"))){
    SST.r<-readRDS(paste(outDir,"SST.RDS",sep="/"))
    SST_sd.r<-readRDS(paste(outDir,"SST_SD.RDS", sep="/"))
    print('LOADING EXISTING SST PRODUCTS')    
  } else {
    print('EXTRACTING SST')
    #data["SST"] <- NA
    #data["SST_SD"] <- NA
    
    # print(i)
    #year <- format(data$dt[i], "%Y")
    if (yr >=1998 & yr < 2018){
      if (dt < "2007-01-01"){
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_SST/Historical product/"
        files <- list.files(paste0(dir,yr))
        parsename <- matrix(unlist(strsplit(files,'-')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        nc <- paste0(dir,yr,'/',files[file_index])
        SST.r<-raster(nc, varname="analysed_sst")
        SST.r <-SST.r-273.15
        #Theta.r<-raster(nc, varname="theta_max")
        rtemp<-rasterextent
        extent(rtemp)<-extent(rasterextent) + c(-5,5,-5,5)
        res(rtemp) <- res(SST.r)
        SST.r2 <- resample(SST.r,rtemp, method="bilinear")
        
        n.sd.cells<-(sd.resolution/res(SST.r2))/2-1
        
        beginCluster()
        SST_sd.r <- clusterR(SST.r2, focal, args=list(w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[2]), fun=sd, na.rm=TRUE, pad=TRUE))
        endCluster()
        
        #SST_sd.r<-focal(SST.r2, w=matrix(1, nc=n.sd.cells[2], nr=n.sd.cells[1]), fun=sd, na.rm=T)
        SST.r <- resample(SST.r,rasterextent, method="bilinear")
        SST_sd.r<-resample(SST_sd.r,rasterextent, method="bilinear")        
        #data.var1  <-  ncvar_get(nc.data,'analysed_sst',start=c(c_low,r_low,1), count=c(numcols,numrows,1),verbose=FALSE)
        #data.var1 <- data.var1 - 273.15
        #data$SST[i] <- mean(data.var1,na.rm=T)
        #data$SST_SD[i] <- sd(data.var1,na.rm=T)
        
        #nc_close(nc.data)
      }
      
      if (dt >= "2007-01-01"){
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_SST/NRT product/"
        #year <- format(dt, "%Y")
        files <- list.files(paste0(dir,yr))
        parsename <- matrix(unlist(strsplit(files,'120000')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        #Open nc and extract
        nc <- paste0(dir,yr,'/',files[file_index])
        SST.r<-raster(nc, varname="analysed_sst")
        #Theta.r<-raster(nc, varname="theta_max")
        rtemp<-rasterextent
        extent(rtemp)<-extent(rasterextent) + c(-5,5,-5,5)
        res(rtemp) <- res(SST.r)
        SST.r2 <- resample(SST.r,rtemp, method="bilinear")
        
        n.sd.cells<-(sd.resolution/res(SST.r2))/2-1
        cellsize<-ceiling(n.sd.cells[1])
        
        beginCluster()
        SST_sd.r <- clusterR(SST.r2, focal, args=list(w=matrix(1, nc=cellsize, nr=cellsize), fun=sd, na.rm=TRUE, pad=T))
        endCluster()
        
        SST.r <- resample(SST.r,rasterextent, method="bilinear")
        SST_sd.r<-resample(SST_sd.r,rasterextent, method="bilinear")
        
        #data.var1  <-  ncvar_get(nc.data,'analysed_sst',start=c(c_low,r_low,1), count=c(numcols,numrows,1),verbose=FALSE)
        #data.var1 <- data.var1 - 273.15
        #data$SST[i] <- mean(data.var1, na.rm=T)
        #data$SST_SD[i] <- sd(data.var1, na.rm=T)
        
        #nc_close(nc.data)
      }
    } 
    
    if (savePieces){
      saveRDS(SST.r,paste(outDir,"SST.RDS",sep="/"))
      saveRDS(SST_sd.r,paste(outDir,"SST_SD.RDS",sep="/"))
    }
  } #else 
  
  #----MLD----
  if (file.exists(paste0(outDir,"/MLD.RDS"))){
    MLD.r<-readRDS(paste(outDir,"MLD.RDS",sep="/"))
    print('LOADING EXISTING MLD PRODUCTS')    
  } else {
    
    #Native resolution is 0.25 so extracting point value
    print('EXTRACTING MLD')
    #data["MLD"] <- NA
    
    # print(i)
    #year <- data$Year[i]
    if (yr>=1998){
      if (dt < "2016-01-01"){
        
        #Figure out which file to grab
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS MLD/"
        files <- list.files(paste0(dir))
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,4], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        nc <- paste0(dir,'/',files[file_index])
        MLD.r<-raster(nc, varname="mlp")
        MLD.r<-resample(MLD.r,rasterextent, method="bilinear")
        #data.var1  <-  ncvar_get(nc.data,'mlp',start=c(c,r,1), #time is indexed this way as I want it to throw an error if there is >1 time field
        #                         count=c(1,1,1),verbose=FALSE)
        #data$MLD[i] <- data.var1
        #nc_close(nc.data)
      }
      
      if (dt >= "2016-01-01" & dt < "2019-01-01"){
        dir <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS MLD Forecast/"
        files <- list.files(paste0(dir))
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        parsename <- matrix(unlist(strsplit(parsename[,6],'b')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,2], format="%Y%m%d")
        file_index <- which.min(abs(as.Date(dt) - dtname))
        nc <- paste0(dir,'/',files[file_index])
        MLD.r<-raster(nc, varname="mlotst")
        MLD.r<-resample(rotate(MLD.r),rasterextent, method="bilinear")
        #data.var1  <-  ncvar_get(nc.data,'mlotst',start=c(c,r,1), #time is indexed this way as I want it to throw an error if there is >1 time field
        #                         count=c(1,1,1),verbose=FALSE)
        #data$MLD[i] <- data.var1
        #nc_close(nc.data)
      }
    }
    ###    MLD.r    
    if (savePieces){
      saveRDS(MLD.r,paste(outDir,"MLD.RDS",sep="/"))
    }    
  } #else
  namesrasterstack<-c("Bathymetry","Rugosity","lat","lon","Oxygen_Surface","Oxygen_100m","SSH","SLA","U","V","Chla_4km_8day","Chla_25km_monthly","FSLE_max","Theta_max","SST","SST_SD","MLD","EKE","SSH_SD","Chla_4km_8day_lag1mo","Chla_25km_monthly_lag1mo")
  rasterdata<-stack(b,b_sd,lat.r,lon.r,osurf,o100,ssh.r,sla.r,U.r,V.r,chl_4k.r,chl_25k.r,FSLE.r,Theta.r,SST.r,SST_sd.r,MLD.r,eke.r,ssh.sd.r,chl_4k_lag.r,chl_25k_lag.r)
  names(rasterdata)<-namesrasterstack
  saveRDS(rasterdata,paste(outDir,"rasterdata.RDS",sep="/"))
  
  return(rasterdata)
  
}

Extract_Satellite <- function(data=data){
  
  #------load packages and library------
  library("ncdf4")
  
  #----BEGIN EXTRACTION----
  #Code differs for each satellite product, all downloaded from CMEMS
  #For more details on each product see Steph who has a product summary table
  
  #----Bathymetry & Rugosity----
  print('EXTRACTING BATHYMETRY')
  data["Bathymetry"] <- NA
  data["Rugosity"] <- NA
  dir <- "/Volumes/Harrable_Starage/Data/GlobalData/"
  nc <- paste0(dir,'ETOPO180_Global/etopo180_4558_b013_12d8.nc')
  nc.data <- nc_open(nc, write=FALSE)
  lat <- ncvar_get(nc.data,'latitude')
  lon <- ncvar_get(nc.data,'longitude')
  nrows <- length(lat); ncols <- length(lon)
  desired.resolution = 0.25
  desired.resolution = desired.resolution/2
  for (i in 1:nrow(data)){
    # print(i)
    # #Compute bathymetry: point extraction
    # c <- which.min(abs(lon-data$lon[i]))
    # r <- which.min(abs(lat-data$lat[i]))
    # data.var  <-  ncvar_get(nc.data,'altitude',start=c(c,r),count=c(1,1),verbose=FALSE)
    
    #Compute Bathymetry & rugosity
    c <- which.min(abs(lon-data$lon[i]))
    c_low <- which.min(abs(lon-(data$lon[i]-desired.resolution)))
    c_up <- which.min(abs(lon-(data$lon[i]+desired.resolution)))
    r <- which.min(abs(lat-data$lat[i]))
    r_low <- which.min(abs(lat-(data$lat[i]-desired.resolution)))
    r_up <- which.min(abs(lat-(data$lat[i]+desired.resolution)))
    numcols=abs(c_up-c_low)+1; numrows=abs(r_up-r_low)+1
    data.var  <-  ncvar_get(nc.data,'altitude',start=c(c_low,r_low),count=c(numcols,numrows),verbose=FALSE)
    data.var <- data.var[data.var<0]
    data$Bathymetry[i] <- mean(data.var, na.rm=T)
    data$Rugosity[i] <- sd(data.var, na.rm=T)
  }
  nc_close(nc.data)
  #----Oxygen: Surface & 100m----
  print('EXTRACTING OXYGEN')
  #1 degree resolution so just extraction point value
  data["Oxygen_Surface"] <- NA
  data["Oxygen_100m"] <- NA
  dir <- "/Volumes/Harrable_Starage/Data/GlobalData/WOA 2013 oxygen monthly clim/"
  for (i in 1:nrow(data)){
    # print(i)
    nc <- paste0(dir,'woa13_all_o',data$Month[i],'_01.nc')
    nc.data <- nc_open(nc, write=FALSE)
    lat <- ncvar_get(nc.data,'lat')
    lon <- ncvar_get(nc.data,'lon')
    depth <- ncvar_get(nc.data,'depth')
    nrows <- length(lat); ncols <- length(lon)
    c <- which.min(abs(lon-data$lon[i]))
    r <- which.min(abs(lat-data$lat[i]))
    data.var.surface  <-  ncvar_get(nc.data,'o_an',start=c(c,r,1,1),count=c(1,1,1,1),verbose=FALSE)
    data$Oxygen_Surface[i] <- data.var.surface

    data.var.depth  <-  ncvar_get(nc.data,'o_an',start=c(c,r,21,1),count=c(1,1,1,1),verbose=FALSE)
    data$Oxygen_100m[i] <- data.var.depth
    nc_close(nc.data)
  }

  #----SSH, SLA, U, and V----
  print('EXTRACTING SSH')
  #Native resolution is 0.25 so extracting point value

  data["SSH"] <- NA
  data["SLA"] <- NA
  data["U"] <- NA
  data["V"] <- NA

  dir <- "/Volumes/Harrable_Starage/Data/GlobalData/CMEMS_SSH_uv/"
  for (i in 1:nrow(data)){
    # print(i)
    #Figure out which file to grab
    year <- data$Year[i]
    # if (year >=2000){
    if (year >=1998){
      files <- list.files(paste0(dir,year))
      parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
      dtname <- as.Date(parsename[,6], format="%Y%m%d")
      file_index <- which.min(abs(data$dt[i] - dtname))

      nc <- paste0(dir,year,'/',files[file_index])
      nc.data <- nc_open(nc, write=FALSE)
      lat <- ncvar_get(nc.data,'latitude')
      lon <- ncvar_get(nc.data,'longitude') #360 degrees
      nrows <- length(lat); ncols <- length(lon)
      data_lon_conv <- ifelse(data$lon[i]<0,data$lon[i]+360,data$lon[i]) #convert lon to 360 degree scale
      c <- which.min(abs(lon-data_lon_conv))
      r <- which.min(abs(lat-data$lat[i]))
      #Variable 1: SLA
      data.var1 <-  ncvar_get(nc.data,'sla',start=c(c,r,1),  count=c(1,1,1),verbose=FALSE)
      data$SSH[i] <- data.var1

      #Variable 2: SSH
      data.var2 <-  ncvar_get(nc.data,'adt',start=c(c,r,1),count=c(1,1,1),verbose=FALSE)
      data$SLA[i] <- data.var2

      #Variable 3: U
      data.var3  <-  ncvar_get(nc.data,'ugos',start=c(c,r,1),count=c(1,1,1),verbose=FALSE)
      data$U[i] <- data.var3

      #Variable 4: V
      data.var4  <-  ncvar_get(nc.data,'vgos',start=c(c,r,1),count=c(1,1,1),verbose=FALSE)
      data$V[i] <- data.var4

      nc_close(nc.data)
    } else{
      data$SSH[i] <- NA
      data$SLA[i] <- NA
      data$U[i] <- NA
      data$V[i] <- NA
    }
  }

  #----Chla historical and NRT: 4km and 8day----
  print('EXTRACTING CHLOROPHYLL 4KM')
  #Chlorophyll-a at 4km over 8day
  data["Chla_4km_8day"] <- NA
  for (i in 1:nrow(data)){
    # print(i)
    year <- data$Year[i]
    # if (year>=2000){
    if (year>=1998){
      if (data$dt[i] < "2017-01-01"){
        dir <- "/Volumes/Harrable_Starage/Data/GlobalData/CMEMS_Chla_4km_8day_Historical/"
        files <- list.files(paste0(dir,year),pattern=".nc")
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(data$dt[i] - dtname))
        nc <- paste0(dir,year,'/',files[file_index])
        nc.data <- nc_open(nc, write=FALSE)
        lat <- ncvar_get(nc.data,'lat')
        lon <- ncvar_get(nc.data,'lon') #180 degrees
        nrows <- length(lat); ncols <- length(lon)
        c <- which.min(abs(lon-data$lon[i]))
        c_low <- which.min(abs(lon-(data$lon[i]-desired.resolution)))
        c_up <- which.min(abs(lon-(data$lon[i]+desired.resolution)))
        r <- which.min(abs(lat-data$lat[i]))
        r_low <- which.min(abs(lat-(data$lat[i]-desired.resolution)))
        r_up <- which.min(abs(lat-(data$lat[i]+desired.resolution)))
        numcols=abs(c_up-c_low)+1; numrows=abs(r_up-r_low)+1

        #Variable 1: Chla
        data.var1 <-  ncvar_get(nc.data,'CHL',start=c(c_low,r_low,1), count=c(numcols,numrows,1),verbose=FALSE)
        data$Chla_4km_8day[i] <- mean(data.var1,na.rm=T)
        nc_close(nc.data)
      }
    } else {
      data$Chla_4km_8day[i]  <- NA
    }

    # if (data$dt[i] >= "2017-01-01"){
    #   nc <- "/Volumes/LaCie1TB/CMEMS_Chla_4km_8day_NRT/dataset-oc-glo-chl-multi-l4-gsm_4km_8days-rt-v02_1520637318265.nc"
    #   nc.data <- nc_open(nc, write=FALSE)
    #   lat <- ncvar_get(nc.data,'lat')
    #   lon <- ncvar_get(nc.data,'lon') #180 degrees
    #   time <- ncvar_get(nc.data,'time') #days since 1900-01-01
    #   time_real <- as.Date(time, origin="1900-01-01")
    #   time_index <- which.min(abs(data$dt[i] - time_real))
    #   nrows <- length(lat); ncols <- length(lon)
    #   c <- which.min(abs(lon-data$lon[i]))
    #   r <- which.min(abs(lat-data$lat[i]))
    #   #Variable 1: SST
    #   data.var1 <-  ncvar_get(nc.data,'CHL',start=c(c,r,time_index),
    #                           count=c(1,1,1),verbose=FALSE)
    #   data$Chla[i] <- data.var1
    #   nc_close(nc.data)
    # }
  }

  #----Chla historical and NRT: 25km and monthly----
  #Native resolution is 0.25 so extracting at the point
  print('EXTRACTING CHLOROPHYLL 25KM')
  data["Chla_25km_monthly"] <- NA
  for (i in 1:nrow(data)){
    # print(i)
    year <- data$Year[i]
    if (year>=1998){
      if (data$dt[i] < "2017-01-01"){
        dir <- "/Volumes/Harrable_Starage/Data/GlobalData/CMEMS_Chla_25km_monthly_Historical/"
        files <- list.files(paste0(dir,year),pattern=".nc")
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(data$dt[i] - dtname))
        nc <- paste0(dir,year,'/',files[file_index])
        nc.data <- nc_open(nc, write=FALSE)
        lat <- ncvar_get(nc.data,'lat')
        lon <- ncvar_get(nc.data,'lon') #180 degrees
        nrows <- length(lat); ncols <- length(lon)
        c <- which.min(abs(lon-data$lon[i]))
        r <- which.min(abs(lat-data$lat[i]))
        #Variable 1: SST
        data.var1 <-  ncvar_get(nc.data,'CHL',start=c(c,r,1), count=c(1,1,1),verbose=FALSE)
        data$Chla_25km_monthly[i] <- data.var1
        nc_close(nc.data)
      } else{
        data$Chla_25km_monthly[i] <- NA
      }
    }
    # if (data$dt[i] >= "2017-01-01"){
    # dir <- "/Volumes/LaCie1TB/CMEMS_Chla_25km_monthly_NRT//"
    #   year <- data$Year[i]
    #   files <- list.files(paste0(dir,year),pattern=".nc")
    #   parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
    #   dtname <- as.Date(parsename[,1], format="%Y%m%d")
    #   file_index <- which.min(abs(data$dt[i] - dtname))
    #   nc <- paste0(dir,year,'/',files[file_index])
    #   nc.data <- nc_open(nc, write=FALSE)
    #   lat <- ncvar_get(nc.data,'lat')
    #   lon <- ncvar_get(nc.data,'lon') #180 degrees
    #   nrows <- length(lat); ncols <- length(lon)
    #   c <- which.min(abs(lon-data$lon[i]))
    #   r <- which.min(abs(lat-data$lat[i]))
    #   #Variable 1: SST
    #   data.var1 <-  ncvar_get(nc.data,'CHL',start=c(c,r,1),
    #                           count=c(1,1,1),verbose=FALSE)
    #   data$Chla_25km_monthly[i] <- data.var1
    #   nc_close(nc.data)
    # }
  }

  #----AVISO FSLE----
  print('EXTRACTING FSLE')
  #Native resolution is 0.25 so extracting point value
  #Looks like I have a corrupt file (dt_global_allsat_madt_fsle_20040908_20141027.nc) and have code below to skip it
  data["FSLE_max"] <- NA
  data["Theta_max"] <- NA
  dir <- "/Volumes/Harrable_Starage/Data/GlobalData/AVISO_FSLE/"
  for (i in 1:nrow(data)){
    # print(i)
    year <- format(data$dt[i], "%Y")
    if (year>=1998){
      #Figure out which file to grab
      files <- list.files(paste0(dir,year))
      parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
      dtname <- as.Date(parsename[,6], format="%Y%m%d")
      file_index <- which.min(abs(data$dt[i] - dtname))
      nc <- paste0(dir,year,'/',files[file_index])
      if (files[file_index] %in% 'dt_global_allsat_madt_fsle_20040908_20141027.nc'){
        data$FSLE_max[i] <- NA
        data$Theta_max[i] <- NA
      } else { #Skip corrupt file
        nc.data <- nc_open(nc, write=FALSE)
        lat <- ncvar_get(nc.data,'lat')
        lon <- ncvar_get(nc.data,'lon')
        time <- ncvar_get(nc.data,'time') #days since 1950-01-01
        time_real <- as.Date(time,origin="1950-01-01")
        nrows <- length(lat); ncols <- length(lon)
        data_lon_conv <- ifelse(data$lon[i]<0,data$lon[i]+360,data$lon[i]) #convert lon to 360 degree scale
        c <- which.min(abs(lon-data_lon_conv))
        r <- which.min(abs(lat-data$lat[i]))
        data.var1  <-  ncvar_get(nc.data,'fsle_max',start=c(c,r,length(time)), #time is indexed this way as I want it to throw an error if there is >1 time field
                                 count=c(1,1,1),verbose=FALSE)
        data$FSLE_max[i] <- data.var1

        data.var2  <-  ncvar_get(nc.data,'theta_max',start=c(c,r,length(time)),
                                 count=c(1,1,1),verbose=FALSE)
        data$Theta_max[i] <- data.var2
        nc_close(nc.data)
      }
    } else {
      data$FSLE_max[i] <- NA
      data$Theta_max[i] <- NA
    }
  }



  #----SST----
  print('EXTRACTING SST')
  data["SST"] <- NA
  data["SST_SD"] <- NA

  for (i in 1:nrow(data)){
    # print(i)
    year <- format(data$dt[i], "%Y")
    if (year >=1998){
      if (data$dt[i] < "2007-01-01"){
        dir <- "/Volumes/Harrable_Starage/Data/GlobalData/CMEMS_SST/Historical product/"
        files <- list.files(paste0(dir,year))
        parsename <- matrix(unlist(strsplit(files,'-')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(data$dt[i] - dtname))
        nc <- paste0(dir,year,'/',files[file_index])
        nc.data <- nc_open(nc, write=FALSE)
        lat <- ncvar_get(nc.data,'lat')
        lon <- ncvar_get(nc.data,'lon') #in 180 degrees
        nrows <- length(lat); ncols <- length(lon)
        c <- which.min(abs(lon-data$lon[i]))
        c_low <- which.min(abs(lon-(data$lon[i]-desired.resolution)))
        c_up <- which.min(abs(lon-(data$lon[i]+desired.resolution)))
        r <- which.min(abs(lat-data$lat[i]))
        r_low <- which.min(abs(lat-(data$lat[i]-desired.resolution)))
        r_up <- which.min(abs(lat-(data$lat[i]+desired.resolution)))
        numcols=abs(c_up-c_low)+1; numrows=abs(r_up-r_low)+1

        data.var1  <-  ncvar_get(nc.data,'analysed_sst',start=c(c_low,r_low,1), count=c(numcols,numrows,1),verbose=FALSE)
        data.var1 <- data.var1 - 273.15

        data$SST[i] <- mean(data.var1,na.rm=T)
        data$SST_SD[i] <- sd(data.var1,na.rm=T)

        nc_close(nc.data)
      }

      if (data$dt[i] >= "2007-01-01"){
        dir <- "/Volumes/Harrable_Starage/Data/GlobalData/CMEMS_SST/NRT product/"
        year <- format(data$dt[i], "%Y")
        files <- list.files(paste0(dir,year))
        parsename <- matrix(unlist(strsplit(files,'120000')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,1], format="%Y%m%d")
        file_index <- which.min(abs(data$dt[i] - dtname))
        #Open nc and extract
        nc <- paste0(dir,year,'/',files[file_index])
        nc.data <- nc_open(nc, write=FALSE)
        lat <- ncvar_get(nc.data,'lat')
        lon <- ncvar_get(nc.data,'lon') #in 180 degrees
        nrows <- length(lat); ncols <- length(lon)
        c <- which.min(abs(lon-data$lon[i]))
        c_low <- which.min(abs(lon-(data$lon[i]-desired.resolution)))
        c_up <- which.min(abs(lon-(data$lon[i]+desired.resolution)))
        r <- which.min(abs(lat-data$lat[i]))
        r_low <- which.min(abs(lat-(data$lat[i]-desired.resolution)))
        r_up <- which.min(abs(lat-(data$lat[i]+desired.resolution)))
        numcols=abs(c_up-c_low)+1; numrows=abs(r_up-r_low)+1

        data.var1  <-  ncvar_get(nc.data,'analysed_sst',start=c(c_low,r_low,1), count=c(numcols,numrows,1),verbose=FALSE)
        data.var1 <- data.var1 - 273.15
        data$SST[i] <- mean(data.var1, na.rm=T)
        data$SST_SD[i] <- sd(data.var1, na.rm=T)

        nc_close(nc.data)
      }
    } else{
      data$SST[i] <- NA
      data$SST_SD[i] <- NA
    }
  }



  #----MLD----
  #Native resolution is 0.25 so extracting point value
  print('EXTRACTING MLD')
  data["MLD"] <- NA

  for (i in 1:nrow(data)){
    # print(i)
    year <- data$Year[i]
    if (year>=1998){
      if (data$dt[i] < "2016-01-01"){

        #Figure out which file to grab
        dir <- "/Volumes/Harrable_Starage/Data/GlobalData/CMEMS MLD/"
        files <- list.files(paste0(dir))
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,4], format="%Y%m%d")
        file_index <- which.min(abs(data$dt[i] - dtname))
        nc <- paste0(dir,'/',files[file_index])
        nc.data <- nc_open(nc, write=FALSE)
        lat <- ncvar_get(nc.data,'latitude')
        lon <- ncvar_get(nc.data,'longitude') #in 180 degrees
        nrows <- length(lat); ncols <- length(lon)
        c <- which.min(abs(lon-data$lon[i]))
        r <- which.min(abs(lat-data$lat[i]))
        data.var1  <-  ncvar_get(nc.data,'mlp',start=c(c,r,1), #time is indexed this way as I want it to throw an error if there is >1 time field
                                 count=c(1,1,1),verbose=FALSE)
        data$MLD[i] <- data.var1
        nc_close(nc.data)
      }

      if (data$dt[i] >= "2016-01-01"){
        dir <- "/Volumes/Harrable_Starage/Data/GlobalData/CMEMS MLD Forecast/"
        files <- list.files(paste0(dir))
        parsename <- matrix(unlist(strsplit(files,'_')),nrow=length(files), byrow=T)
        parsename <- matrix(unlist(strsplit(parsename[,6],'b')),nrow=length(files), byrow=T)
        dtname <- as.Date(parsename[,2], format="%Y%m%d")
        file_index <- which.min(abs(data$dt[i] - dtname))
        nc <- paste0(dir,'/',files[file_index])
        nc.data <- nc_open(nc, write=FALSE)
        lat <- ncvar_get(nc.data,'lat')
        lon <- ncvar_get(nc.data,'lon') #in 360 degrees
        nrows <- length(lat); ncols <- length(lon)
        data_lon_conv <- ifelse(data$lon[i]<0,data$lon[i]+360,data$lon[i]) #convert lon to 360 degree scale
        c <- which.min(abs(lon-data_lon_conv))
        r <- which.min(abs(lat-data$lat[i]))
        data.var1  <-  ncvar_get(nc.data,'mlotst',start=c(c,r,1), #time is indexed this way as I want it to throw an error if there is >1 time field
                                 count=c(1,1,1),verbose=FALSE)
        data$MLD[i] <- data.var1
        nc_close(nc.data)
      }
    }else{
      data$MLD[i] <- NA
    }
  }
  
  return(data)
  
}



