### running model eval functions
### written by HW, SB, BA 2019

#### Summary Statistics for fit models
tag="final"
load(file=paste0("bwhalemodels_",tag,".RData"))
load(file=paste0("bwhaledata_",tag,".RData"))

listGAMMs<-list(bwhaleGAMM_CRW,bwhaleGAMM_back,bwhaleGAMM_buff,bwhaleGAMM_rev)
listGLMMs<-list(bwhaleGLMM_CRW,bwhaleGLMM_back,bwhaleGLMM_buff,bwhaleGLMM_rev)
listBRTs<-list(bwhaleBRT.lr005.CRW,bwhaleBRT.lr005.back,bwhaleBRT.lr005.buff,bwhaleBRT.lr005.rev)
listofBRTData<-list(BRTtransformDataFrame_BW(bwhaledata_CRW),BRTtransformDataFrame_BW(bwhaledata_back),BRTtransformDataFrame_BW(bwhaledata_buff),BRTtransformDataFrame_BW(bwhaledata_rev))
listofGAMMData<-list(na.omit(bwhaledata_CRW),na.omit(bwhaledata_back),na.omit(bwhaledata_buff),na.omit(bwhaledata_rev))
listnames<-c("Blue_Whale_CRW","Blue_Whale_background","Blue_Whale_buffer","Blue_Whale_reverseCRW")

summarystatistics <- saveModelStats(listGAMMs,listGLMMs,listBRTs,listofBRTData,listofGAMMData,listnames,run.bh=TRUE)
summarystatistics$R_squared<-unlist(summarystatistics$R_squared)
summarystatistics$AUC<-unlist(summarystatistics$AUC)
summarystatistics$TSS<-unlist(summarystatistics$TSS)
write.csv(summarystatistics, file=paste0("Output/SummaryStatistics_bwhale_",tag,".csv"))
summarystatistics<-read.csv(file=paste0("Output/SummaryStatistics_bwhale_",tag,".csv"))
summarystatistics<-summarystatistics %>% subset(select = -c(X))

histogramplots(listofBRTData,listofGAMMData,listnames)

listGAMMs<-list(ephantGAMM_CRW,ephantGAMM_back,ephantGAMM_buff,ephantGAMM_rev)
listGLMMs<-list(ephantGLMM_CRW,ephantGLMM_back,ephantGLMM_buff,ephantGLMM_rev)
listBRTs<-list(ephantBRT.lr005.CRW,ephantBRT.lr005.back,ephantBRT.lr005.buff,ephantBRT.lr005.rev)
listnames<-c("elephant_CRW","elephant_background","elephant_buffer","elephant_reverseCRW")
listofBRTData<-list(ephant_downsample_6hrs(ephantdata_CRW),ephant_downsample_6hrs(ephantdata_back),ephant_downsample_6hrs(ephantdata_buff),ephant_downsample_6hrs(ephantdata_rev))
listofGAMMData<-list(na.omit(ephant_downsample_6hrs(ephantdata_CRW)),na.omit(ephant_downsample_6hrs(ephantdata_back)),na.omit(ephant_downsample_6hrs(ephantdata_buff)),na.omit(ephant_downsample_6hrs(ephantdata_rev)))

summarystatistics <- saveModelStats(listGAMMs,listGLMMs,listBRTs,listofBRTData,listofGAMMData,listnames,run.bh=TRUE)
summarystatistics$R_squared<-unlist(summarystatistics$R_squared)
summarystatistics$AUC<-unlist(summarystatistics$AUC)
summarystatistics$TSS<-unlist(summarystatistics$TSS)
write.csv(summarystatistics, file=paste0("Output/SummaryStatistics_ephant_",tag,".csv"))

histogramplots(listofBRTData,listofGAMMData,listnames,sp="ele")


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

