
#library(lattice)

plotROC <- function(truth, predicted, ...){
  pred <- prediction(as.vector(abs(predicted)), as.vector(truth))    
  roc <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  plot(roc, ...)
  abline(a=0, b= 1)
  print(auc@y.values)
}

UTM2WGS <- function(d){
  xy <- sp::SpatialPoints(cbind(d$long, d$lat), proj4string = CRS("+proj=utm +south +zone=33 +ellps=WGS84"))
  
  latlong <- spTransform(xy, CRS('+proj=longlat +datum=WGS84')) # reproject to lat/long
  d$long[d$long>1000] <- latlong@coords[d$long>1000,1]; d$lat[d$lat>1000] <- latlong@coords[d$lat>1000,2]
  return(d)
}

saveAUC <- function(truth, predicted, ...){
  pred <- prediction(as.vector(abs(predicted)), as.vector(truth))    
  roc <- performance(pred,"tpr","fpr")
  auc <- performance(pred,"auc")
  #plot(roc, ...)
  #abline(a=0, b= 1)
  return(auc@y.values)
}

saveTSS <- function(truth, predicted, ...){
  pred <- prediction(as.vector(abs(predicted)), as.vector(truth))    
  TSS <- performance(pred, "sens", "spec")
  TSSvals <- max(unlist(TSS@y.values) + unlist(TSS@x.values) - 1)
  return(TSSvals)
}

kfolds_eval <- function(dataInput, gbm.x, gbm.y, lr=lr){
  DataInput <- dataInput
  DataInput$Kset <- kfold(DataInput,5) #randomly allocate k groups
  Evaluations_kfold <- as.data.frame(matrix(data=0,nrow=5,ncol=4))
  colnames(Evaluations_kfold) <- c("k","Deviance","AUC","TSS")
  counter=1
  for (k in 1:5){
    print(k)
    DataInput_train <- DataInput[DataInput$Kset!=k,]
    DataInput_test <- DataInput[DataInput$Kset==k,]
    DataInput.kfolds <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                 family="bernoulli", tree.complexity=3,
                                 learning.rate = lr, bag.fraction = 0.6)
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=DataInput.kfolds$gbm.call$best.trees, type="response")
    dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$presabs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    e <- dismo::evaluate(p=pres, a=abs)
    Evaluations_kfold[counter,1] <- k
    Evaluations_kfold[counter,2] <- dev
    Evaluations_kfold[counter,3] <- e@auc
    Evaluations_kfold[counter,4] <- max(e@TPR + e@TNR-1)
    counter=counter+1 
  }
  return(Evaluations_kfold)}

pseudoR2.BRT <- function(x){
  d2 <- 1-(x$self.statistics$resid.deviance/x$self.statistics$null.deviance)
  return(d2)
}

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
addslopeaspect<-function(datafile,slope=NA,aspect=NA)
{
  datafile$slope[]<-NA
  datafile$aspect[]<-NA
  for (i in 1:length(datafile$lat)){
    nearpt<-spDistsN1(cbind(slope$lon,slope$lat),c(datafile$lon[i],datafile$lat[i]), longlat = TRUE)
    slopept<-slope[which(nearpt==min(nearpt)),]
    aspectpt<-aspect[which(nearpt==min(nearpt)),]
    datafile$slope[i]<-slopept
    datafile$aspect[i]<-aspectpt
  }
  
  return(datafile)
}

matchpts<-function(masterfile,tempfile){
  coordinates(masterfile) <- c("lon", "lat")
  coordinates(tempfile) <- c("lon", "lat")      
  
  closestSiteVecS <- vector(mode = "numeric",length = nrow(masterfile))
  minDistVecS     <- vector(mode = "numeric",length = nrow(masterfile))
  
  for (i in 1 : nrow(masterfile))
  {
    distVec <- spDistsN1(tempfile,masterfile[i,],longlat = TRUE)
    minDistVecS[i] <- min(distVec)
    closestSiteVecS[i] <- which.min(distVec)
  }
  return(closestSiteVecS)
}

makedataset<-function(presence,absencedata){
  #numsims<-max(absencedata$iteration)
  datasample<-sample(absencedata$iteration)[1:length(unique(presence$tag))]
  presence$presabs<-1; absencedata$presabs<-0
  if(exists("spdata") && is.data.frame(get("spdata"))) rm(spdata)
  for (x in 1:length(unique(presence$tag))){
    inddata<-rbind(presence[presence$tag==unique(presence$tag)[x],],absencedata[absencedata$tag==unique(presence$tag)[x]&absencedata$iteration==datasample[x],])
    if(exists("spdata") && is.data.frame(get("spdata"))) spdata<-rbind(spdata,inddata) else spdata<-inddata
  }
  
  return(spdata)
}

bhattacharyya.stat.ephant<-function(data){
  ### test Bhattacharyya's coefficient ranges from 0 to 1, with 0 being no overlap in distributions and 1 being perfect overlap. 
  source("http://tguillerme.github.io/R/bhatt.coef.R")
  data.pres<-(data$NDVI[data$presabs==1]); data.pres<-data.pres[!is.na(data.pres)]
  data.abs<-(data$NDVI[data$presabs==0]); data.abs<-data.abs[!is.na(data.abs)]
  bh.NDVI<-bhatt.coeff(data.pres,data.abs)
  data.pres<-(data$roaddist[data$presabs==1]); data.pres<-data.pres[!is.na(data.pres)]
  data.abs<-(data$roaddist[data$presabs==0]); data.abs<-data.abs[!is.na(data.abs)]
  bh.road<-bhatt.coeff(data.pres,data.abs)
  data.pres<-(data$lat[data$presabs==1]); data.pres<-data.pres[!is.na(data.pres)]
  data.abs<-(data$lat[data$presabs==0]); data.abs<-data.abs[!is.na(data.abs)]
  bh.lat<-bhatt.coeff(data.pres,data.abs)
  data.pres<-(data$waterdist[data$presabs==1]); data.pres<-data.pres[!is.na(data.pres)]
  data.abs<-(data$waterdist[data$presabs==0]); data.abs<-data.abs[!is.na(data.abs)]
  bh.water<-bhatt.coeff(data.pres,data.abs)

  
  #bhatt.coeff(BRTData[[i]]$Oxygen_100m[BRTData[[i]]$presabs==1], BRTData[[i]]$Oxygen_100m[BRTData[[i]]$presabs==0])
  #bhatt.coeff(BRTData[[i]]$lat[BRTData[[i]]$presabs==1], BRTData[[i]]$lat[BRTData[[i]]$presabs==0])
  #bhatt.coeff(BRTData[[i]]$SST[BRTData[[i]]$presabs==1], BRTData[[i]]$SST[BRTData[[i]]$presabs==0])
  
  outvect<-c(bh.NDVI,bh.road,bh.lat,bh.water)
  names(outvect)<-c("bh.NDVI","bh.road","bh.lat","bh.water")
  return(outvect)
}

bhattacharyya.stat.whales<-function(data){
  ### test Bhattacharyya's coefficient ranges from 0 to 1, with 0 being no overlap in distributions and 1 being perfect overlap. 
  source("http://tguillerme.github.io/R/bhatt.coef.R")
  data.pres<-(data$Chla_25km_monthly[data$presabs==1]); data.pres<-data.pres[!is.na(data.pres)]
  data.abs<-(data$Chla_25km_monthly[data$presabs==0]); data.abs<-data.abs[!is.na(data.abs)]
  bh.chl<-bhatt.coeff(data.pres,data.abs)
  data.pres<-(data$SST[data$presabs==1]); data.pres<-data.pres[!is.na(data.pres)]
  data.abs<-(data$SST[data$presabs==0]); data.abs<-data.abs[!is.na(data.abs)]
  bh.SST<-bhatt.coeff(data.pres,data.abs)
  data.pres<-(data$lat[data$presabs==1]); data.pres<-data.pres[!is.na(data.pres)]
  data.abs<-(data$lat[data$presabs==0]); data.abs<-data.abs[!is.na(data.abs)]
  bh.lat<-bhatt.coeff(data.pres,data.abs)
  data.pres<-(data$Oxygen_100m[data$presabs==1]); data.pres<-data.pres[!is.na(data.pres)]
  data.abs<-(data$Oxygen_100m[data$presabs==0]); data.abs<-data.abs[!is.na(data.abs)]
  bh.oxy<-bhatt.coeff(data.pres,data.abs)
  
  #bhatt.coeff(BRTData[[i]]$Oxygen_100m[BRTData[[i]]$presabs==1], BRTData[[i]]$Oxygen_100m[BRTData[[i]]$presabs==0])
  #bhatt.coeff(BRTData[[i]]$lat[BRTData[[i]]$presabs==1], BRTData[[i]]$lat[BRTData[[i]]$presabs==0])
  #bhatt.coeff(BRTData[[i]]$SST[BRTData[[i]]$presabs==1], BRTData[[i]]$SST[BRTData[[i]]$presabs==0])
  
  outvect<-c(bh.chl,bh.SST,bh.lat,bh.oxy)
  names(outvect)<-c("bh.chl","bh.SST","bh.lat","bh.oxy")
  return(outvect)
}

createWhaleSamples<-function(bwhalepresence,bwhaleabsbuffer,bwhaleabsbackgrnd,bwhaleabsCRW,bwhaleabsCRWrev,tag=Sys.Date()){
  bwhalemerge <- bwhalepresence %>% select(dTime, tag, iteration, long, lat, Month, Year, Bathymetry, Rugosity, Oxygen_Surface, Oxygen_100m, SSH, SLA, U, V, Chla_4km_8day, Chla_25km_monthly, FSLE_max, Theta_max, SST, SST_SD, MLD)
  
  bwhaleabsmerge <- bwhaleabsbuffer %>% select(dTime, ID, iteration, long, lat, Month, Year, Bathymetry, Rugosity, Oxygen_Surface, Oxygen_100m, SSH, SLA, U, V, Chla_4km_8day, Chla_25km_monthly, FSLE_max, Theta_max, SST, SST_SD, MLD)
  names(bwhaleabsmerge)[2] <- "tag"
  bwhaledata_buff<-makedataset(bwhalemerge,bwhaleabsmerge)
  bwhaledata_buff$eke<-(bwhaledata_buff$U^2+bwhaledata_buff$V^2)/2
  
  bwhaleabsmerge <- bwhaleabsbackgrnd %>% select(dTime, ID, iteration, long, lat, Month, Year, Bathymetry, Rugosity, Oxygen_Surface, Oxygen_100m, SSH, SLA, U, V, Chla_4km_8day, Chla_25km_monthly, FSLE_max, Theta_max, SST, SST_SD, MLD)
  names(bwhaleabsmerge)[2] <- "tag"
  bwhaledata_back<-makedataset(bwhalemerge,bwhaleabsmerge)
  bwhaledata_back$eke<-(bwhaledata_back$U^2+bwhaledata_back$V^2)/2
  
  bwhaleabsmerge <- bwhaleabsCRW %>% select(t, ID, iteration, x, y, Month, Year, Bathymetry, Rugosity, Oxygen_Surface, Oxygen_100m, SSH, SLA, U, V, Chla_4km_8day, Chla_25km_monthly, FSLE_max, Theta_max, SST, SST_SD, MLD)
  names(bwhaleabsmerge)[c(1,2,4,5)] <- c("dTime","tag","long","lat")
  bwhaledata_CRW<-makedataset(bwhalemerge,bwhaleabsmerge)
  bwhaledata_CRW$eke<-(bwhaledata_CRW$U^2+bwhaledata_CRW$V^2)/2
  
  bwhaleabsmerge <- bwhaleabsCRWrev %>% select(t, ID, iteration, x, y, Month, Year, Bathymetry, Rugosity, Oxygen_Surface, Oxygen_100m, SSH, SLA, U, V, Chla_4km_8day, Chla_25km_monthly, FSLE_max, Theta_max, SST, SST_SD, MLD)
  names(bwhaleabsmerge)[c(1,2,4,5)] <- c("dTime","tag","long","lat")
  bwhaledata_rev<-makedataset(bwhalemerge,bwhaleabsmerge)
  bwhaledata_rev$eke<-(bwhaledata_rev$U^2+bwhaledata_rev$V^2)/2
  
  saveRDS(bwhaledata_buff, file=paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_buff_",tag,".RDS"))
  saveRDS(bwhaledata_back, file=paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_back_",tag,".RDS"))
  saveRDS(bwhaledata_CRW, file=paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_CRW_",tag,".RDS"))
  saveRDS(bwhaledata_rev, file=paste0("~/Documents/R/github/PseudoAbsence/Data/bwhaledata_rev_",tag,".RDS"))
}

createElephantSamples<-function(ephantpresence,ephantabsbuffer,ephantabsbackgrnd,ephantabsCRW,ephantabsCRWrev,tag=Sys.Date()){
  ephantmerge <- elephantdata %>% select(date, id, iteration, x, y, roaddist, NDVI, waterdist)
  names(ephantmerge)[c(1,2,4,5)] <- c("dTime","tag","long","lat")
  
  ephantabsmerge <- ephantabsbuffer %>% select(dTime, ID, iteration, long, lat, roaddist, NDVI, waterdist)
  names(ephantabsmerge)[2] <- "tag"
  
  ephantdata_buff<-makedataset(ephantmerge,ephantabsmerge)
  
  ephantabsmerge <- ephantabsbackgrnd %>% select(dTime, ID, iteration, long, lat, roaddist, NDVI, waterdist)
  names(ephantabsmerge)[2] <- "tag"
  
  ephantdata_back<-makedataset(ephantmerge,ephantabsmerge)
  
  ephantabsmerge <- ephantabsCRW %>% select(t, ID, iteration, long, lat, roaddist, NDVI, waterdist)
  names(ephantabsmerge)[c(1:2)] <- c("dTime","tag")
  
  ephantdata_CRW<-makedataset(ephantmerge,ephantabsmerge)
  
  ephantabsmerge <- ephantabsCRWrev %>% select(t, ID, iteration, long, lat, roaddist, NDVI, waterdist)
  names(ephantabsmerge)[c(1:2)] <- c("dTime","tag")
  
  ephantdata_rev<-makedataset(ephantmerge,ephantabsmerge)
  
  saveRDS(ephantdata_buff, file=paste0("~/Documents/R/github/PseudoAbsence/Data/ephantdata_buff_",tag,".RDS"))
  saveRDS(ephantdata_back, file=paste0("~/Documents/R/github/PseudoAbsence/Data/ephantdata_back_",tag,".RDS"))
  saveRDS(ephantdata_CRW, file=paste0("~/Documents/R/github/PseudoAbsence/Data/ephantdata_CRW_",tag,".RDS"))
  saveRDS(ephantdata_rev, file=paste0("~/Documents/R/github/PseudoAbsence/Data/ephantdata_rev_",tag,".RDS"))
}


BRTtransformDataFrame_BW<-function(datafr){
  tempdata <- datafr
  tempdata <- tempdata[c(23,1:22,24)]
  names(tempdata)[5]<-"lon"
  tempdata$Chla_4km_8day <- log10(tempdata$Chla_4km_8day + 0.001)
  tempdata$Chla_25km_monthly <- log10(tempdata$Chla_25km_monthly + 0.001)
  tempdata$eke <- log10(tempdata$eke + 0.001)
  #names(tempdata)[c(17,18,24)]<-c("log10(Chla_4km_8day)","log10(Chla_25km_monthly)","log10(eke)")
  tempdata$RN<-sample(1000000,size=dim(tempdata)[1])
  return(tempdata)
}

fixdTime<-function(timedata){
  #temptime<-str_match(elephantdata$date, "T(.*?)Z")  #str_match(CRWfiles[-grep("AG", CRWfiles)], "sim_(.*?).csv")
  tempdate<-str_split(timedata,"T")
  tempdate<-lapply(tempdate,function(x) gsub("Z","",x))
  tempdate<-lapply(tempdate,function(x) paste(x[1],x[2],sep=' '))
  return(as.POSIXct(unlist(tempdate), format = "%Y-%m-%d %H:%M:%S"))
}

ephant_downsample_6hrs<-function(datafr){
  times <- format(datafr$dTime, "%H%M")
  dsample <- which (times=="0600" | times=="1200" | times=="1800" | times=="0000")
  return(datafr[dsample,])
}

saveModelStats<-function(listofGAMMs,listofGLMMs,listofBRTs,BRTData,GAMMData,listnames,run.bh=FALSE){
  library(caret)
  nummodels = length(listofGAMMs)*3
  summarystats<-data.frame("Model_Name" = rep(NA, nummodels), "R_squared" = rep(NA, nummodels), "AIC" = rep(NA, nummodels), "AUC" = rep(NA, nummodels), "TSS" = rep(NA,nummodels), "Sensitivity" = rep(NA,nummodels), "Specificity" = rep(NA,nummodels), "FalsePos" = rep(NA,nummodels), "FalseNeg" = rep(NA,nummodels), "Accuracy" = rep(NA,nummodels))
  # FOR DEBUGGING 
  #listofGAMMs<-listGAMMs;listofGLMMs<-listGLMMs;listofBRTs<-listBRTs;BRTData<-listofBRTData;GAMMData<-listofGAMMData
  
  items<-seq(1,length(listofGAMMs),by=1)
  times<-seq(0,6,by=2)
  for(i in items){
    summarystats$Model_Name[i+times[i]]<-paste(listnames[i],"GAMM",sep="_")
    summarystats$Model_Name[i+1+times[i]]<-paste(listnames[i],"GLMM",sep="_")
    summarystats$Model_Name[i+2+times[i]]<-paste(listnames[i],"BRT",sep="_")
    
    GAMMmodel<-listofGAMMs[[i]]
    GLMMmodel<-listofGLMMs[[i]]
    BRTmodel<-listofBRTs[[i]]
    BRTdata<-BRTData[[i]]
    GAMMdata<-GAMMData[[i]]
    
    pred=predict.gbm(BRTmodel,newdata=BRTdata,type="response", n.trees=BRTmodel$gbm.call$best.trees,na.rm=F)
    summarystats$AIC[i+times[i]]<-AIC(GAMMmodel$lme)
    summarystats$AIC[i+1+times[i]]<-AIC(GLMMmodel$lme)
    summarystats$AIC[i+2+times[i]]<-NA
    
    CV_RSq <- (cor(pred, BRTdata$presabs))^2
    Train_RSq <- (cor(pred, BRTdata$presabs))^2
    
    summarystats$R_squared[i+times[i]]<-summary(GAMMmodel$gam)[10]
    summarystats$R_squared[i+1+times[i]]<-summary(GLMMmodel$gam)[10]
    summarystats$R_squared[i+2+times[i]]<-pseudoR2.BRT(BRTmodel)

    newdata<-predict.gam(GAMMmodel$gam, GAMMdata,se=TRUE, type="response"); newdata$presabs <- GAMMdata$presabs; 
    AUCval<-saveAUC(newdata$presabs,newdata$fit)
    TSSval<-saveTSS(newdata$presabs,newdata$fit)
    specval<-specificity(factor(round(newdata$presabs)),factor(round(newdata$fit)))
    sensval<-sensitivity(factor(round(newdata$presabs)),factor(round(newdata$fit)))
    cm<-confusionMatrix(factor(round(newdata$presabs)),factor(round(newdata$fit)))
      
    summarystats$Sensitivity[i+times[i]]<-sensval
    summarystats$Specificity[i+times[i]]<-specval
    summarystats$FalsePos[i+times[i]]<-mean(newdata$fit[newdata$presabs==0])
    summarystats$FalseNeg[i+times[i]]<-mean(newdata$fit[newdata$presabs==1])
    summarystats$Accuracy[i+times[i]]<-cm$overall['Accuracy']
    summarystats$AUC[i+times[i]]<-AUCval
    summarystats$TSS[i+times[i]]<-TSSval
    pdf(paste("GAMM_ROCR_",listnames[i],".pdf",sep=''))
    plotROC(newdata$presabs,newdata$fit, colorize = TRUE)
    dev.off()
    
    newdata<-predict.gam(GLMMmodel$gam, GAMMdata,se=TRUE, type="response"); newdata$presabs <- GAMMdata$presabs; 
    AUCval<-saveAUC(newdata$presabs,newdata$fit)
    TSSval<-saveTSS(newdata$presabs,newdata$fit)
    specval<-specificity(factor(round(newdata$presabs)),factor(round(newdata$fit)))
    sensval<-sensitivity(factor(round(newdata$presabs)),factor(round(newdata$fit)))
    cm<-confusionMatrix(factor(round(newdata$presabs)),factor(round(newdata$fit)))
    
    summarystats$FalsePos[i+1+times[i]]<-mean(newdata$fit[newdata$presabs==0])
    summarystats$FalseNeg[i+1+times[i]]<-mean(newdata$fit[newdata$presabs==1])
    summarystats$Accuracy[i+1+times[i]]<-cm$overall['Accuracy']
    summarystats$Sensitivity[i+1+times[i]]<-sensval
    summarystats$Specificity[i+1+times[i]]<-specval
    summarystats$AUC[i+1+times[i]]<-AUCval
    summarystats$TSS[i+1+times[i]]<-TSSval
    pdf(paste("GLMM_ROCR_",listnames[i],".pdf",sep=''))
    plotROC(newdata$presabs,newdata$fit, colorize = TRUE)
    dev.off()
    
    #newdata<-predict.gam(GAMMmodel$gam, GAMMdata,se=TRUE, type="response"); newdata$presabs <- GAMMdata$presabs; 
    AUCval<-saveAUC(BRTdata$presabs,pred)
    TSSval<-saveTSS(BRTdata$presabs,pred)
    specval<-specificity(factor((BRTdata$presabs)),factor(round(pred)))
    sensval<-sensitivity(factor((BRTdata$presabs)),factor(round(pred)))
    cm<-confusionMatrix(factor((BRTdata$presabs)),factor(round(pred)))
    
    summarystats$FalsePos[i+2+times[i]]<-mean(pred[BRTdata$presabs==0])
    summarystats$FalseNeg[i+2+times[i]]<-mean(pred[BRTdata$presabs==1])
    summarystats$Accuracy[i+2+times[i]]<-cm$overall['Accuracy']
    summarystats$Sensitivity[i+2+times[i]]<-sensval
    summarystats$Specificity[i+2+times[i]]<-specval
    summarystats$AUC[i+2+times[i]]<-AUCval
    summarystats$TSS[i+2+times[i]]<-TSSval
    
    pdf(paste("BRT_ROCR_",listnames[i],".pdf",sep=''))
    plotROC(BRTdata$presabs,pred, colorize = TRUE)
    dev.off()
    
    if (run.bh){
      if (grepl("Blue",listnames[1])) bhvals<-bhattacharyya.stat.whales(BRTdata)
      else (bhvals<-bhattacharyya.stat.ephant(BRTdata))
      if (i==1) bhvect<-bhvals else bhvect<-rbind(bhvect,bhvals) 
      #summarystats[i+times[i],]<-c(summarystats[i+times[i],],bhvals)
      #summarystats[i+1+times[i],]<-summarystats[i+1+times[i],]
      #summarystats[i+2+times[i],]<-summarystats[i+2+times[i],]
      
    }
    
    ## Pseudo-R2
     #pseudoR2 <- pseudoR2.BRT(fittedModel)
     
     ## Explained deviance
     DN <- BRTmodel$self.statistics$null.deviance 
     DR <- BRTmodel$cv.statistics$resid.deviance 
     DE <- (DN-DR)*100/DN; DE
    
    pdf(paste("BRT_responses_",listnames[i],".pdf",sep=''))
    ### gbm.plot is failing below:
    #Error in if ((min(i.var) < 1) || (max(i.var) > length(x$var.names))) { : 
    #    missing value where TRUE/FALSE needed
    par(mfrow=c(4,4))
    gbm.plot(BRTmodel, write.title = FALSE)
    dev.off()
    
    pdf(paste("GAMM_responses_",listnames[i],".pdf",sep=''))
    par(mfrow=c(4,3))
    plot(GAMMmodel$gam, shade=TRUE)
    dev.off()
    
    pdf(paste("GAMM_check_",listnames[i],".pdf",sep=''))
    par(mfrow=c(2,2))
    gam.check(GAMMmodel$gam, shade=TRUE)
    dev.off()
    
    pdf(paste("BRT_summary_",listnames[i],".pdf",sep=''))
    summary(listofBRTs[[i]])
    dev.off()
    
  }
  
  if (run.bh) {
    bhtibble<-as_tibble(bhvect) %>% 
    slice(rep(1:n(), each = 3))
    summarystats<-cbind(summarystats,bhtibble)
  } 
  return(summarystats)
}

plotPresAbsPts<-function(species_name="Blue_Whale", species_data, xlimits=c(-130,-110), ylimits=c(20,50)){
  library(maps)
  if(mean(species_data$lat) > 1000) {
    species_data<-UTM2WGS(species_data)
  } else if (species_name == "Blue_Whale") species_data <- species_data[species_data$Bathymetry<(-25),]
    ### Approach to test on land
    #points <- expand.grid(species_data$long, species_data$lat)  # Note that I reversed OP's ordering of lat/long
    #pts <- SpatialPoints(points, proj4string=CRS(proj4string(wrld_simpl)))
    #ii <- !is.na(over(pts, wrld_simpl)$FIPS)
  
  
#    if (mean(species_data$Bathymetry,na.rm=TRUE)<(-100)) species_data <- species_data[species_data$Bathymetry<(-25),]
    pdf(paste(species_name,"_presabs_map_",Sys.Date(),".pdf",sep=''))
    maps::map(xlim=xlimits, ylim=ylimits)         # pacific-centric projection
    map.axes()
    points(species_data$long[species_data$presabs==0],species_data$lat[species_data$presabs==0], col=makeTransparent('red',60), pch=16, cex=.5)
    points(species_data$long[species_data$presabs==1],species_data$lat[species_data$presabs==1], col=makeTransparent('blue',60), pch=16, cex=.5)
  
    dev.off()

  # pdf(paste(species_name,"_presabs_map.pdf",sep=''))
  # map('worldHires', xlim=c(-140,-110), ylim=c(20,60))         # pacific-centric projection
  # map.axes()
  # points(GAMdataRunX$lon[GAMdataRunCCS$presabs==0],GAMdataRunX$lat[GAMdataRunCCS$presabs==0], col=makeTransparent('red',85), pch=16, cex=.5)
  # points(GAMdataRunX$lon[GAMdataRunCCS$presabs==1],GAMdataRunX$lat[GAMdataRunCCS$presabs==1], col=makeTransparent('blue',85), pch=16, cex=.5)
  # rect(-135, 30, -115, 49, density = NULL, angle = 45, col = NA, border = "black")
  # dev.off()
  
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



predictsurfaces<-function(BRTmodels,GAMMmodels,BRTData,GAMMdata,listnames,stacksurfaces){
  
}

create_env_rasters=function(day_nc){  # fast verion, use this
  ### Test code
  #nc <- "/Volumes/HazWht_8TB/Data/GlobalData/CMEMS_SST/NRT\ product/2012/20120801120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc"
  #varnames=c("analysed_sst","lon","lat")
  #day_nc <- "2001-08-01 12:00:00 UTC"
  #dname <- varnames[1]
  dir <- "/Volumes/HazWht_8TB/Data/GlobalData/"
  outdir <- "/Users/elliotthazen/Documents/R/github/PseudoAbsence/Data/"
  
  ### BATHY AND RUGOSITY
  nc <- paste0(dir,'ETOPO180_Global/etopo180_4558_b013_12d8.nc')
  nc.data <- nc_open(nc, write=FALSE)
  lat <- ncvar_get(nc.data,'latitude')
  lon <- ncvar_get(nc.data,'longitude')
  nrows <- length(lat); ncols <- length(lon)
  dat1=list()
  dat1$x=lon #- 360
  dat1$y=lat
  dat1$z=ncvar_get(nc.data, 'altitude',verbose=FALSE)
  nc_close(nc.data)
  
  r <-raster(
    dat1$z,
    xmn=range(dat1$x)[1], xmx=range(dat1$x)[2],
    ymn=range(dat1$y)[1], ymx=range(dat1$y)[2], 
    crs=CRS("+proj=longlat +datum=WGS84")
  )
  r2=flip(r,2)
  r2sd=flip(focal(r2, w=matrix(1,nrow=21,ncol=21), fun=sd,na.rm=TRUE,pad=T),2)
  
  r2 <- crop(r2,c(-150,-11,20,60))
  template <- raster("/Users/elliotthazen/Dropbox/shared folders/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd")
  Bathy_r <- raster::resample(r2, template, method="bilinear")  
  Rugos_r <- raster::resample(r2sd, template, method="bilinear")  
  
  x##SST  
  nc <- paste0(dir,'CMEMS_SST/NRT\ product/2012/20120801120000-UKMO-L4_GHRSST-SSTfnd-OSTIA-GLOB-v02.0-fv02.0.nc')
  nc.data=nc_open(nc)
#  lat <- ncvar_get(nc.data,'latitude')
#  lon <- ncvar_get(nc.data,'longitude')
  lat <- ncvar_get(nc.data,"lat")
  lon <- ncvar_get(nc.data,"lon")
  nrows <- length(lat); ncols <- length(lon)
  tim <- ncvar_get(nc.data,'time')
  day <- as.POSIXlt(tim,origin='1970-01-01',tz= "UTC")
  indx <- which.min(abs(day - as.POSIXlt(day_nc)))
  startdt <- as.Date(day[indx])
  tmp.array <- ncvar_get(nc.data, dname)
  fillvalue <- ncatt_get(nc.data, dname, "_FillValue")
  # lat.index <- which(lat>30 & lat < 48)
  # lon.index <- which(lon>(-134+360) & lon < (-115.5+360))
  tmp.array.day=tmp.array#[,,indx]
  tmp.array.day[tmp.array.day==fillvalue$value]=NA #setting fill value
  
  dat1=list()
  dat1$x=lon #- 360
  dat1$y=lat
  dat1$z=log(t(tmp.array.day))
  
  r <-raster(
    dat1$z,
    xmn=range(dat1$x)[1], xmx=range(dat1$x)[2],
    ymn=range(dat1$y)[1], ymx=range(dat1$y)[2], 
    crs=CRS("+proj=longlat +datum=WGS84")
  )
  r2=flip(r,2)
  r2 <- crop(r2,c(-150,-11,20,60))
  r2sd=focal(r2, w=matrix(1,nrow=5,ncol=5), fun=sd,na.rm=TRUE,pad=T)
  
  #Force to template res
  template <- raster("/Users/elliotthazen/Dropbox/shared folders/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd")
  SST_r <- raster::resample(r2, template, method="bilinear")  
  SSTsd_r <- raster::resample(r2sd, template, method="bilinear")  
  
  print('EXTRACTING OXYGEN')
  #1 degree resolution so just extraction point value
  datafile <- "/Volumes/Harrable_Starage/Data/GlobalData/WOA 2013 oxygen monthly clim/"
  nc <- paste0(dir,'/GlobalData/WOA 2013 oxygen monthly clim/woa13_all_o08_01.nc')
  nc.data <- nc_open(nc, write=FALSE)
  lat <- ncvar_get(nc.data,'lat')
  lon <- ncvar_get(nc.data,'lon')
  depth <- ncvar_get(nc.data,'depth')
  nrows <- length(lat); ncols <- length(lon)
  c <- which.min(abs(lon-data$lon[i]))
  r <- which.min(abs(lat-data$lat[i]))
  dat1=list()
  dat1$x=lon #- 360
  dat1$y=lat
  dat1$z=ncvar_get(nc.data,'o_an',start=c(c,r,1,1),count=c(ncols,nrows,1,1),verbose=FALSE)
  dat1$z100  <-  ncvar_get(nc.data,'o_an',start=c(c,r,21,1),count=c(ncols,nrows,1,1),verbose=FALSE)
  nc_close(nc.data)
  r <-raster(
    dat1$z,
    xmn=range(dat1$x)[1], xmx=range(dat1$x)[2],
    ymn=range(dat1$y)[1], ymx=range(dat1$y)[2], 
    crs=CRS("+proj=longlat +datum=WGS84")
  )
  r2=flip(r,2)
  r2 <- crop(r2,c(-150,-11,20,60))

  #Force to template res
  template <- raster("/Users/elliotthazen/Dropbox/shared folders/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd")
  O2_r <- raster::resample(r2, template, method="bilinear")  
  
  r <-raster(
    dat1$z100,
    xmn=range(dat1$x)[1], xmx=range(dat1$x)[2],
    ymn=range(dat1$y)[1], ymx=range(dat1$y)[2], 
    crs=CRS("+proj=longlat +datum=WGS84")
  )
  r2=flip(r,2)
  r2 <- crop(r2,c(-150,-11,20,60))
  
  #Force to template res
  template <- raster("/Users/elliotthazen/Dropbox/shared folders/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd")
  O2_100_r <- raster::resample(r2, template, method="bilinear")  
  
  ### PICK UP HERE #----SSH, SLA, U, and V----
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
  
  #writeRaster(SST_r, paste0(outdir,"SST_",day_nc,".grd"))
  #writeRaster(SSTsd_r, paste0(outdir,"SSTsd_",day_nc,".grd"))
  
  }
}

stackephantpredictors<-function(){
  
}

histogramplots<-function(BRTData,GAMMData,listnames,sp="BW"){
  items<-seq(1,length(listnames))
  for(i in items){
    species_name<-listnames[i]

    if (sp == "BW"){    
      p1 <- ggplot(na.omit(BRTData[[i]])) + geom_density(aes(x = SST, fill = as.factor(presabs)), alpha = 0.2) + theme_classic() + guides(fill=guide_legend(title=species_name))
      p2 <- ggplot(na.omit(BRTData[[i]])) + geom_density(aes(x = Oxygen_100m, fill = as.factor(presabs)), alpha = 0.2) + theme_classic() + guides(fill=guide_legend(title=species_name))
      p3 <- ggplot(na.omit(BRTData[[i]])) + geom_density(aes(x = lat, fill = as.factor(presabs)), alpha = 0.2) + theme_classic() + guides(fill=guide_legend(title=species_name))
      p4 <- ggplot(na.omit(BRTData[[i]])) + geom_density(aes(x = Chla_25km_monthly, fill = as.factor(presabs)), alpha = 0.2) + theme_classic() + guides(fill=guide_legend(title=species_name))
      pdf(paste(species_name,"_envVars_densityplot_",Sys.Date(),".pdf",sep=''))
        multiplot(p1, p2, p3, p4, cols=2)
        #hist(GAMMData$SST[GAMMData$presabs==1,],breaks=30)
      dev.off()
      
      
      
      # if(exists("BRTData$long")) {
      #     xlims=c(summary(BRTData[[i]]$long)[2]-3,summary(BRTData[[i]]$long)[5]+3)
      #   } else {
      #     xlims=c(summary(BRTData[[i]]$lon)[2]-3,summary(BRTData[[i]]$lon)[5]+3)
      #   }
      if(!exists("BRTData[[i]]$long")) BRTData[[i]]$long<-BRTData[[i]]$lon
      xlims=c(summary(BRTData[[i]]$long)[2]-3,summary(BRTData[[i]]$long)[5]+3)
      ylims=c(summary(BRTData[[i]]$lat)[2]-3,summary(BRTData[[i]]$lat)[5]+3)
      plotPresAbsPts(species_name=listnames[i], BRTData[[i]], xlimits=c(-130,-110), ylimits=c(20,50))
      
    } else {
      p1 <- ggplot(na.omit(BRTData[[i]])) + geom_density(aes(x = NDVI, fill = as.factor(presabs)), alpha = 0.2) + theme_classic() + guides(fill=guide_legend(title=species_name))
      p2 <- ggplot(na.omit(BRTData[[i]])) + geom_density(aes(x = roaddist, fill = as.factor(presabs)), alpha = 0.2) + theme_classic() + guides(fill=guide_legend(title=species_name))
      p3 <- ggplot(na.omit(BRTData[[i]])) + geom_density(aes(x = waterdist, fill = as.factor(presabs)), alpha = 0.2) + theme_classic() + guides(fill=guide_legend(title=species_name))
      pdf(paste(species_name,"_SST_densityplot_",Sys.Date(),".pdf",sep=''))
        multiplot(p1, p2, p3, cols=1)
      dev.off()
      
      #xlims=c(summary(UTM2WGS(BRTData[[i]])$long)[2]-3,summary(UTM2WGS(BRTData[[i]])$long)[5]+3)
      #ylims=c(summary(UTM2WGS(BRTData[[i]])$lat)[2]-3,summary(UTM2WGS(BRTData[[i]])$lat)[5]+3)
      xlims=c(summary((BRTData[[i]])$long)[2]-3,summary((BRTData[[i]])$long)[5]+3)
      ylims=c(summary((BRTData[[i]])$lat)[2]-3,summary((BRTData[[i]])$lat)[5]+3)
      plotPresAbsPts(species_name=listnames[i], BRTData[[i]], xlimits=xlims, ylimits=ylims)
    }
    
  }
  
  
  
}

histogram_bytype<-function(BRTData,GAMMData,listnames,sp="BW",leg="none", colortypes=c("red","blue","green","yellow",NA)){
  items<-seq(1,length(listnames))

  #BRTDataAllPres<-BRTDataAll[BRTDataAll$presabs==1]

  #ggplot(na.omit(BRTDataAll)) + geom_density(aes(x = SST, fill = as.factor(abstype), color = as.factor(abstype)), alpha = 0.3) + theme_classic() + guides(fill=guide_legend(title="Blue whale")) + scale_fill_manual(values = c("green","blue","red","yellow",NA)) + scale_color_manual(values = c(NA,NA,NA,NA,"black"),guide = 'none') 
  
  BRTData[[1]]$abstype=listnames[1]
  BRTData[[2]]$abstype=listnames[2]
  BRTData[[3]]$abstype=listnames[3]
  BRTData[[4]]$abstype=listnames[4]
  
  BRTDataAll<-rbind(BRTData[[1]],BRTData[[2]],BRTData[[3]],BRTData[[4]])
  BRTDataAll$abstype[BRTDataAll$presabs==1]<-"Presence"
  
  
  if (sp == "BW"){        #Chla_25km_monthly

    p1 <- ggplot(na.omit(BRTDataAll)) + geom_density(aes(x = SST, fill = as.factor(abstype), color = as.factor(abstype)), alpha = 0.3) + theme_classic() + guides(fill=guide_legend(title="Blue whale")) + scale_fill_manual(values = colortypes) + scale_color_manual(values = c(NA,NA,NA,NA,"black"),guide = 'none') + theme(legend.position = "none")
    p2 <- ggplot(na.omit(BRTDataAll)) + geom_density(aes(x = Oxygen_100m, fill = as.factor(abstype), color = as.factor(abstype)), alpha = 0.3) + theme_classic() + guides(fill=guide_legend(title="Blue whale")) + scale_fill_manual(values = colortypes) + scale_color_manual(values = c(NA,NA,NA,NA,"black"),guide = 'none') + theme(legend.position = "none")
    p3 <- ggplot(na.omit(BRTDataAll)) + geom_density(aes(x = Bathymetry, fill = as.factor(abstype), color = as.factor(abstype)), alpha = 0.3) + theme_classic() + guides(fill=guide_legend(title="Blue whale")) + scale_fill_manual(values = colortypes) + scale_color_manual(values = c(NA,NA,NA,NA,"black"),guide = 'none')  + theme(legend.position = "none")
    p4 <- ggplot(na.omit(BRTDataAll)) + geom_density(aes(x = Chla_25km_monthly, fill = as.factor(abstype), color = as.factor(abstype)), alpha = 0.3) + theme_classic() + guides(fill=guide_legend(title="Blue whale")) + scale_fill_manual(values = colortypes) + scale_color_manual(values = c(NA,NA,NA,NA,"black"),guide = 'none')  + theme(legend.position = leg)
    pdf(paste("BlueWhale_envVars_densityplot_",Sys.Date(),".pdf",sep=''))
    multiplot(p1, p2, p3, p4, cols=2)
    dev.off()
    
  } else {
 
    p1 <- ggplot(na.omit(BRTDataAll)) + geom_density(aes(x = NDVI, fill = as.factor(abstype), color = as.factor(abstype)), alpha = 0.3) + theme_classic() + guides(fill=guide_legend(title="Elephant")) + scale_fill_manual(values = colortypes) + scale_color_manual(values = c(NA,NA,NA,NA,"black"),guide = 'none')  + theme(legend.position = "none")
    p2 <- ggplot(na.omit(BRTDataAll)) + geom_density(aes(x = roaddist, fill = as.factor(abstype), color = as.factor(abstype)), alpha = 0.3) + theme_classic() + guides(fill=guide_legend(title="Elephant")) + scale_fill_manual(values = colortypes) + scale_color_manual(values = c(NA,NA,NA,NA,"black"),guide = 'none')  + theme(legend.position = "none")
    p3 <- ggplot(na.omit(BRTDataAll)) + geom_density(aes(x = waterdist, fill = as.factor(abstype), color = as.factor(abstype)), alpha = 0.3) + theme_classic() + guides(fill=guide_legend(title="Elephant")) + scale_fill_manual(values = colortypes) + scale_color_manual(values = c(NA,NA,NA,NA,"black"),guide = 'none')  + theme(legend.position = leg)
    
    pdf(paste("Elephant_envVars_densityplot_",Sys.Date(),".pdf",sep=''))
      multiplot(p1, p2, p3, cols=1)
    dev.off()
      
    }
}

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

