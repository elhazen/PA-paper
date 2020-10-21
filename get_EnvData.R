# This is an R function to download all environmental data used in the GAMM models. It names the variables properly to run in the GAMM prediction.
#This currently has CHLA, SST, and Y-Wind operational.

# IMPORTANT - The script assumes that it is in the correct parent directory, and then sets relative paths from there.

#Code validated to get env variables 5/18/2015 by EAH

get_EnvData <- function() {

  ##############HEADER########################
  # Define functions used in script.
  
  # Function to generalize grabbing data from ERDDAP. 
  #Download chl and change dataframe format for GMT to regrid chlorophyll to match SST
  #Use [(last)] when obtaining real-time data from ERDDAP, e.g. [(2014-08-16T00:00:00Z):1:(2014-08-16T00:00:00Z)] = [(last)]
  curldap <- function(dapurl, outfile) {
  	f = CFILE(outfile,mode="wb")
  	id = curlPerform(url=dapurl,writedata=f@ref) 
  	close(f)
  	return(id)
  }
  
  #Load required libraries. Function pkgTest is in the file Code/load_Functions.R and should have been loaded. If not test here and load file.
  
  if(exists("pkgTest")==FALSE) {
     logprint("Function pkgTest not found, loading file Code/load_Functions.R...")
     source("Code/load_Functions.R")
     }
  
  logprint("Loading required libraries")   
  pkgTest("gmt")
  pkgTest("SDMTools")
  pkgTest("ncdf")
  pkgTest("RCurl")
  pkgTest("raster")
  pkgTest("reshape")
  
  #load master raster grid for resampling (0.25x0.25 degrees)
  load("Data/MasterGrid.RData")
  logprint("Getting the current month and year...")
  todaydate = format(Sys.time(), "%Y-%m-%d")
  year = as.numeric(format(Sys.time(), "%Y"))
  month = as.numeric(format(Sys.time(), "%m"))
  
  # Create temporary directory to hold files while creating product. Folder is named tmp.[seconds since 1970]. You can convert to date with as.Date([seconds since 1970]/86400,origin="1970-01-01")
  tmpdir = paste("tmp.",as.integer(Sys.time()),sep="")
  logprint(paste("Creating temp directory",tmpdir))
  dir.create(tmpdir, showWarnings = TRUE)
  
  #http://coastwatch.pfeg.noaa.gov/erddap/convert/time.txt?n=473472000&units=seconds%20since%201970-01-01T00:00:00Z #Debugging
  
  #Get Environmental variables. 
  # First is to grab Chlorophyll data as this is the limiting factor. First step is to try from ERDDAP. If this fails can grab from other sources.
  initialtime = "last" #use this to get the last available ERDDAP file
  
  #For debugging you can enter in a specific ERDDAP time in the format 2015-04-16T00:00:00Z
  #initialtime = "2015-10-16T00:00:00Z" #use this to get the last available ERDDAP file
  
  #Call ERDDAP to get the date of "latest" chlorophyll file
  dapurl="http://coastwatch.pfeg.noaa.gov/erddap/tabledap/allDatasets.csv0?maxTime&datasetID=%22erdMBchlamday%22"
  h = basicTextGatherer()
  id=curlPerform(url=dapurl, writefunction=h$update)
  testdate=as.Date(h$value())
  factorfile = paste("Data/WhaleWatchFactors_",format(testdate,"%Y_%m"),".csv",sep="")
  
  # Do first sanity check - Make sure we actually need the file before moving on
  logprint(paste("Checking for existence of factor file",factorfile))
  
  if (file.exists(factorfile)){
    logprint(paste("Factor file ", factorfile, " present, moving to modeling phase",sep="")) #Problem here
    logprint(paste("Cleaning up temp directory", tmpdir))
    unlink(tmpdir,recursive=TRUE)
    return(factorfile)
  }
  
  logprint("Factorfile not found, grabbing data")
  
  #Now get latest file
  dapurl=paste("http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMBchlamday.nc?chlorophyll[(",initialtime,")][(0.0):1:(0.0)][(29):1:(49)][(224):1:(245)]",sep="")
  logprint(paste("Attempting to grab Chlorophyll using",dapurl))
  chlfile = paste(tmpdir,"/chl.nc",sep="")
  
  idchl = curldap(dapurl, chlfile)
  
  #Check id to make sure it exits cleanly (e.g., id=0) #Maybe add try command here?
  if(idchl!=0) {
  	logprint(paste("Download failed, exit code =",id,"Trying secondary source"))
    #Need to actually put secondary source here...
  } else {
    logprint("Chlorophyll data file download successful!!!")
  }
  
  #Read in chl.nc file to get month and year to load other data files
  logprint(paste("Reading in chlorophyll data file",chlfile))
  chlnc = open.ncdf(chlfile,write=FALSE)
  chldate = as.Date(get.var.ncdf(chlnc,"time")/(60*60*24),origin="1970-01-01")
  chlmonth = as.integer(format(chldate,"%m"))
  chlyear = as.integer(format(chldate,"%Y"))
  
  chllon=get.var.ncdf(chlnc,"longitude")
  chllat=get.var.ncdf(chlnc,"latitude")
  
  erdtime = format(chldate,"%Y-%m-%dT00:00:00Z")
  
  #Get Chlorophyll data and regrid to mastergrid
  chlraw = get.var.ncdf( chlnc, "chlorophyll")
  rownames(chlraw) = chllon #Need to have lon and lat for melt reshape
  colnames(chlraw) = chllat
  chlxyz = melt(chlraw, c("longitude","latitude"))
  coordinates(chlxyz) <- ~ longitude + latitude     #### check order here of V1 and V2
  gridded(chlxyz) <- TRUE
  # coerce to raster
  rasterchl <- raster(chlxyz)
  rasterchlmedian <- focal(rasterchl, w=matrix(1/25,nrow=5,ncol=5))
  chl<-resample(rasterchlmedian, mastergrid, method="bilinear")
  
  #Download SST and change dataframe format for GMT to regrid 
  # erdMBsstdmday - ERDDAP PFEG Pacific Ocean at 0.0125 (PREFERRED AND SAME RES AS CHL)
  sstfile = paste(tmpdir,"/sst.nc",sep="")
  dapurl = paste("http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMBsstdmday.nc?sst[(",erdtime,"):1:(",erdtime,")][(0.0):1:(0.0)][(29):1:(49)][(224):1:(245)]",sep="")
  
  logprint(paste("Attempting SST download", dapurl))
  idsst = curldap(dapurl, sstfile)
  
  #Check id to make sure it exits cleanly (e.g., id=0)
  if(idsst!=0) {
  	print(paste("Download failed, exit code =",id,"Trying secondary source"))
  }
  
  sstnc = open.ncdf(sstfile,write=FALSE)
  sstdate = as.Date(get.var.ncdf(sstnc,"time")/(60*60*24),origin="1970-01-01")
  sstlon=get.var.ncdf(sstnc,"longitude")
  sstlat=get.var.ncdf(sstnc,"latitude")
  
  sstraw = get.var.ncdf( sstnc, "sst")
  rownames(sstraw) = sstlon #Need to have lon and lat for melt reshape
  colnames(sstraw) = sstlat
  sstxyz = melt(sstraw, c("longitude","latitude"))
  coordinates(sstxyz) <- ~ longitude + latitude     #### check order here of V1 and V2
  gridded(sstxyz) <- TRUE
  # coerce to raster
  rastersst <- raster(sstxyz)
  sst<-resample(rastersst, mastergrid, method="bilinear")
  
  #Do some sanity checks on dates and locations
  if(sstdate-chldate!=0) { logprint("Problem: Dates not identical")}
  if(sum(sstlon-chllon)!=0) { logprint("Problem: Longitudes not identical")}
  if(sum(sstlat-chllat)!=0) { logprint("Problem: Latitudes not identical")}
  
  
  #Everything OK, close Chla and SST NetCDF files
  close.ncdf(sstnc)
  close.ncdf(chlnc)
  
  #Download SSH and change dataframe format for GMT to regrid the same way to match other variables
  # The primary data place is AVISO, and we will pull 16 days centered on the mid-point day of the latest Chla file. 
  #
  #>>>>>>>MAY NEED TO CHANGE THIS AS THERE IS DELAY<<<<<<<<
  #
  #
  # Normally will pull from near-real-time data, unless running particular month where we need delayed time data
  # First try near real time data
  #IF THIS FAILS!!! Try the username and pass aviso-users:grid2010 as this may have changed?
  dapurl = "http://aviso-users:grid2010@opendap.aviso.altimetry.fr/thredds/dodsC/dataset-duacs-nrt-over30d-global-allsat-msla-h"
  
  #Load SSH OpenDAP file from AVISO
  logprint(paste("Attempting SSH download", dapurl))
  sshnc = open.ncdf(dapurl,write=FALSE)
  
  #Get all dates in NetCDF file
  sshdates = as.Date(get.var.ncdf(sshnc,"time"),origin="1950-01-01")
  
  #Find initial indices to pull data chunk. Pixels are centered and SST and Chl are not so subtracting 0.125 to get left edge of pixel. May want to look at this but shouldn"t affect model.
  LonStartIdx <- which( sshnc$dim$lon$vals-0.125 == 225)
  LatStartIdx <- which( sshnc$dim$lat$vals-0.125 == 30)
  TimeStartIdx <- which(sshdates==as.character(sstdate-16)) #16 days before SST and Chla
  SLAIdx <- which(names(sshnc$var)=="sla") # find variable index for SLA
  
  #This subroutine not exactly foolproof, may fail...
  if (length(TimeStartIdx) == 0) {
    logprint("Cannot find date in NRT, need to using delayed time product")
    close(sshnc)
    dapurl = "http://aviso-users:grid2010@opendap.aviso.altimetry.fr/thredds/dodsC/dataset-duacs-dt-global-twosat-msla-h"
    sshnc = open.ncdf(dapurl,write=FALSE)
    sshdates = as.Date(get.var.ncdf(sshnc,"time"),origin="1950-01-01")
    LonStartIdx <- which( sshnc$dim$lon$vals-0.125 == 225)
    LatStartIdx <- which( sshnc$dim$lat$vals-0.125 == 30)
    TimeStartIdx <- which(sshdates==as.character(sstdate-16)) #16 days before SST and Chla
    SLAIdx <- which(names(sshnc$var)=="sla") # find variable index for SLA
  }
  
  #Get data slice
  #If you change grid you will have to change counts for lon and lat
  sshslice = get.var.ncdf( sshnc, "sla", start=c(LonStartIdx,LatStartIdx,TimeStartIdx), count=c(81,77,16))
  
  #Get lat and long, need to subtract 0.125 to fix pixel center vs left corner mismatch with SST and Chla
  sshlon=get.var.ncdf(sshnc,"lon",start=c(LonStartIdx),count=c(81))-0.125 
  sshlat=get.var.ncdf(sshnc,"lat",start=c(LatStartIdx),count=c(77))-0.125
  
  close.ncdf(sshnc)
  
  #Do average in third (time) dimension
  sshraw = apply(sshslice, c(1,2), mean)
  rownames(sshraw) = sshlon #Need to have lon and lat for melt reshape
  colnames(sshraw) = sshlat
  
   #Turn matrix into xyz to regrid, then make raster
  logprint("Regridding SSH to mastergrid")
  sshxyz = melt(sshraw, c("longitude","latitude"))
  coordinates(sshxyz) <- ~ longitude + latitude     #### check order here of V1 and V2
  gridded(sshxyz) <- TRUE
  # coerce to raster
  rasterssh <- raster(sshxyz)
  ssh<-resample(rasterssh, mastergrid, method="bilinear")
  
  #calculate SD for ssh
  logprint("Calculating SD for SSH raster")
  sshsdraw<-focal(ssh, w=matrix(1,nrow=5,ncol=5), fun=sd,na.rm=TRUE)
  sshrms<-resample(sshsdraw, mastergrid, method="bilinear") #Do we need this if in same grid?
   
  #Load in reshaped SST, Chla, and bathy data
  logprint("Loading in reshaped SST and Chl-a data")
  sst = as.data.frame(sst,xy=TRUE) #turn SST raster into xyz file
  chl = as.data.frame(chl,xy=TRUE) #turn Chl raster into xyz file
  ssh = as.data.frame(ssh,xy=TRUE) #turn SSH raster into xyz file
  sshrms = as.data.frame(sshrms,xy=TRUE) #turn SSHSD raster into xyz file
  bathy = read.table("Data/bathy.txt")
  bathy = bathy[,c(5,6,1,2)] #Only take bathy, bathyrms, long, and lat
  
  colnames(sst) = c("longitude","latitude","sst")
  colnames(chl) = c("longitude","latitude","chl")
  colnames(ssh) = c("longitude","latitude","ssh")
  colnames(sshrms) = c("longitude","latitude","sshrms")
  colnames(bathy) = c("longitude","latitude","bathy","bathyrms")
  
  #Can we turn this into a function do have a dynamic number of merges?
  modelin = merge(merge(merge(merge(bathy,sst,by=c("longitude","latitude")),chl,by=c("longitude","latitude")),ssh,by=c("longitude","latitude")),sshrms,by=c("longitude","latitude"))
  colnames(modelin) = c("longitude","latitude","bathy","bathyrms","sst","chl","ssh","sshrms")
  
  modelin$month = chlmonth #month derived from chlorophyll file
  modelin$year = chlyear #month derived from chlorophyll file
  
  logprint(paste("Writing factorfile to file",factorfile))
  write.csv(modelin,factorfile,row.names=FALSE)
  
  #Do checks on files and then remove temp directory
  logprint(paste("Cleaning up temp directory", tmpdir))
  unlink(tmpdir,recursive=TRUE)
  
  return(factorfile)
}
