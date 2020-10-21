#Data MUST contain: a date file called 'dt', 'lon' ,'lat', "Year", and "Month" 
data <- read.csv('./Data/BlueWhale_SSM.csv')
data$dt <- as.Date(paste(data$Year,data$Month,data$Day,sep="-")) #make dt
data$lon <- data$long #easy name change
data$Month <- format(data$dt, "%m") #important for oxygen files later
t1 <- Sys.time()
data <- Extract_Satellite(data=data)
t2 <- Sys.time() #Takes 1.7 hours for the SSM file
write.csv(data,'./Data/XtractedData/BlueWhale_SSM.csv', row.names = FALSE)

dir <- '~/Dropbox/PseudoAbsence-MS/Data/CRW_output/'
files <- list.files(dir,full.names = FALSE,pattern=".csv")
files <- files[grep('AG',files,invert = TRUE)]

#Process Background and Buffer point files   
for (f in 1:416){
  print(paste0("Extracting absence file: " ,f, " of 416"))
  d <- read.csv(paste0(dir,files[f]))
  d$dt <- as.Date(d$dTime, format="%Y-%m-%d")
  d$lon <- d$long
  d$Month <- format(d$dt, "%m")
  d$Year <- format(d$dt, "%Y")
  if (d$Year[1] >=1998 & d$Year[1] <=1999){
    d <- Extract_Satellite(data=d)
    write.csv(d,paste0('~/Dropbox/PseudoAbsence-MS/Data/XtractedData/Extracted_',files[f]), row.names = FALSE)
  } else {
    print ("File not in correct time period")
  }
}
