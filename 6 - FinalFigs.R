###### BLAH BLAH BLAH BLAH BLAH

world <- ne_countries(scale = "medium", returnclass = "sf")
#class(world)

etosha <- readOGR("./Data/Etosha shapefiles", "ENP_Pans_clipped")
# Next the shapefile has to be converted to a dataframe for use in ggplot2
shapefile_df <- fortify(etosha)

# Now the shapefile can be plotted as either a geom_path or a geom_polygon.
# Paths handle clipping better. Polygons can be filled.
# You need the aesthetics long, lat, and group.
map <- ggplot() +
  geom_path(data = shapefile_df, 
            aes(x = long, y = lat, group = group),
            color = 'gray', fill = 'white', size = .2)

print(map) 


tag="final"
load(paste0("bwhaledata_",tag,".RData"))
load(paste0("bwhalemodels_",tag,".RData"))

### Read in ELEPHANTS
tag="final"

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

###  match with ephant_downsample_6hrs
ephantdata_CRW <- ephant_downsample_6hrs(ephantdata_CRW)
ephantdata_back <- ephant_downsample_6hrs(ephantdata_back)
ephantdata_buff <- ephant_downsample_6hrs(ephantdata_buff)

load(paste0("ephantdata",tag,".RData"))
load(paste0("ephantmodels",tag,".RData"))

bwhaledata_CRW<-bwhaledata_CRW[order(bwhaledata_CRW$presabs),]
bwhaledata_rev<-bwhaledata_rev[order(bwhaledata_rev$presabs),]
bwhaledata_back<-bwhaledata_back[order(bwhaledata_back$presabs),]
bwhaledata_buff<-bwhaledata_buff[order(bwhaledata_buff$presabs),]
ephantdata_CRW<-ephantdata_CRW[order(ephantdata_CRW$presabs),]
ephantdata_rev<-ephantdata_rev[order(ephantdata_rev$presabs),]
ephantdata_back<-ephantdata_back[order(ephantdata_back$presabs),]
ephantdata_buff<-ephantdata_buff[order(ephantdata_buff$presabs),]

sz=0.05
leg="none" #=c(0.9, 0.9)
savefilename=paste0("Absence_plots_",as.Date(Sys.Date(),format='%m/%d/%Y'),".pdf")

### First figure - tracks v. pseudotracks each approach   
#https://www.rdocumentation.org/packages/ggExtra/versions/0.9/topics/ggMarginal
#(sites <- data.frame(long = c(-140, -115), lat = c(32, 45)))
bw_CRW_plot <- ggplot(data = world) +
  geom_point(data = bwhaledata_CRW, aes(x = long, y = lat, colour = factor(presabs, levels=c(1,0))), size = sz, 
             shape = 16) +
  geom_sf() +
  coord_sf(xlim = c(-140, -115), ylim = c(32, 45), expand = FALSE) +
  theme_classic()  + 
  scale_color_manual(values=c("blue","red"),name = "Blue Whales", labels = c("CRW Absences","Presences")) +
  theme(legend.position = leg)

bw_rev_plot <- ggplot(data = world) +
  geom_point(data = bwhaledata_rev, aes(x = long, y = lat, colour = factor(presabs, levels=c(1,0))), size = sz, 
             shape = 16) +
  geom_sf() +
  coord_sf(xlim = c(-140, -115), ylim = c(32, 45), expand = FALSE) +
  theme_classic()  + 
  scale_color_manual(values=c("blue","red"),name = "Blue Whales", labels = c("Reverse CRW Absences","Presences")) +
  theme(legend.position = leg)

bw_buffer_plot <- ggplot(data = world) +
  geom_point(data = bwhaledata_buff, aes(x = long, y = lat, colour = factor(presabs, levels=c(1,0))), size = sz, 
             shape = 16) +
  geom_sf() +
  coord_sf(xlim = c(-140, -115), ylim = c(32, 45), expand = FALSE) +
  theme_classic()  + 
  scale_color_manual(values=c("blue","red"),name = "Blue Whales", labels = c("Buffer Absences","Presences")) +
  theme(legend.position = leg) #c(0.9, 0.9)

bw_back_plot <- ggplot(data = world) +
  geom_point(data = bwhaledata_back, aes(x = long, y = lat, colour = factor(presabs, levels=c(1,0))), size = sz, 
             shape = 16) +
  geom_sf() +
  coord_sf(xlim = c(-140, -115), ylim = c(32, 45), expand = FALSE) +
  theme_classic()  + 
  scale_color_manual(values=c("blue","red"),name = "Blue Whales", labels = c("Background Absences","Presences")) +
  theme(legend.position = leg) #c(0.9, 0.9)

ele_CRW_plot <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = ephantdata_CRW, aes(x = long, y = lat, colour = factor(presabs, levels=c(1,0))), size = sz, 
             shape = 16) +
  coord_sf(xlim = c(14, 17), ylim = c(-20, -18), expand = FALSE) +
  theme_classic()  + 
  scale_color_manual(values=c("blue","red"),name = "Elephants", labels = c("CRW Absences","Presences")) +
  geom_path(data = shapefile_df, aes(x = long, y = lat, group = group)) +
  theme(legend.position = leg) #c(0.9, 0.9)

ele_rev_plot <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = ephantdata_rev, aes(x = long, y = lat, colour = factor(presabs, levels=c(1,0))), size = sz, 
             shape = 16) +
  coord_sf(xlim = c(14, 17), ylim = c(-20, -18), expand = FALSE) +
  theme_classic()  + 
  scale_color_manual(values=c("blue","red"),name = "Elephants", labels = c("Reverse CRW Absences","Presences")) +
  geom_path(data = shapefile_df, aes(x = long, y = lat, group = group)) +
  theme(legend.position = leg) #c(0.9, 0.9)

ele_buffer_plot <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = ephantdata_buff, aes(x = long, y = lat, colour = factor(presabs, levels=c(1,0))), size = sz, 
             shape = 16) +
  coord_sf(xlim = c(14, 17), ylim = c(-20, -18), expand = FALSE) +
  theme_classic()  + 
  scale_color_manual(values=c("blue","red"),name = "Elephants", labels = c("Buffer Absences","Presences")) +
  geom_path(data = shapefile_df, aes(x = long, y = lat, group = group)) +
  theme(legend.position = leg) #c(0.9, 0.9)

ele_back_plot <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = ephantdata_back, aes(x = long, y = lat, colour = factor(presabs, levels=c(1,0))), size = sz, 
             shape = 16) +
  coord_sf(xlim = c(14, 17), ylim = c(-20, -18), expand = FALSE) +
  theme_classic()  + 
  scale_color_manual(values=c("blue","red"),name = "Elephants", labels = c("Background Absences","Presences")) +
  geom_path(data = shapefile_df, aes(x = long, y = lat, group = group)) +
  theme(legend.position = leg) #c(0.9, 0.9)

pdf (file=savefilename, height=8, width=11)
par(mar=c(5,4,4,2))
multiplot(bw_CRW_plot,bw_rev_plot,bw_back_plot,bw_buffer_plot,ele_CRW_plot,ele_rev_plot,ele_back_plot,ele_buffer_plot, cols=2)
dev.off()

pdf (file="studyarea.pdf", height=8, width=11)
par(mar=c(5,4,4,2))
ggplot(data = world) +
  geom_sf() +
  geom_point(data = ephantdata_buff[ephantdata_buff$presabs==1,], aes(x = long, y = lat, colour = "green"), size = sz, 
             shape = 16) +
  theme_classic()  + 
  theme(legend.position = "none") #c(0.9, 0.9)
dev.off()


### Second figure - env. space density plots
listGAMMs<-list(bwhaleGAMM_CRW,bwhaleGAMM_back,bwhaleGAMM_buff,bwhaleGAMM_rev)
listGLMMs<-list(bwhaleGLMM_CRW,bwhaleGLMM_back,bwhaleGLMM_buff,bwhaleGLMM_rev)
listBRTs<-list(bwhaleBRT.lr005.CRW,bwhaleBRT.lr005.back,bwhaleBRT.lr005.buff,bwhaleBRT.lr005.rev)
listofBRTData<-list(BRTtransformDataFrame_BW(bwhaledata_CRW),BRTtransformDataFrame_BW(bwhaledata_back),BRTtransformDataFrame_BW(bwhaledata_buff),BRTtransformDataFrame_BW(bwhaledata_rev))
listofGAMMData<-list(na.omit(bwhaledata_CRW),na.omit(bwhaledata_back),na.omit(bwhaledata_buff),na.omit(bwhaledata_rev))
listnames<-c("CRW","background","buffer","CRWrev")

histogram_bytype(listofBRTData,listofGAMMData,listnames,sp="BW",leg=c(0.8,0.9), colortypes = c("orange","red","green","blue",NA)) #plot

listGAMMs<-list(ephantGAMM_CRW,ephantGAMM_back,ephantGAMM_buff,ephantGAMM_rev)
listGLMMs<-list(ephantGLMM_CRW,ephantGLMM_back,ephantGLMM_buff,ephantGLMM_rev)
listBRTs<-list(ephantBRT.lr005.CRW,ephantBRT.lr005.back,ephantBRT.lr005.buff,ephantBRT.lr005.rev)
listnames<-c("elephant_CRW","elephant_background","elephant_buffer","elephant_reverseCRW")
listofBRTData<-list(ephant_downsample_6hrs(ephantdata_CRW),ephant_downsample_6hrs(ephantdata_back),ephant_downsample_6hrs(ephantdata_buff),ephant_downsample_6hrs(ephantdata_rev))
listofGAMMData<-list(na.omit(ephant_downsample_6hrs(ephantdata_CRW)),na.omit(ephant_downsample_6hrs(ephantdata_back)),na.omit(ephant_downsample_6hrs(ephantdata_buff)),na.omit(ephant_downsample_6hrs(ephantdata_rev)))

histogram_bytype(listofBRTData,listofGAMMData,listnames,sp="elephant",leg=c(0.8,0.7), colortypes = c("orange","red","green","blue",NA))   #c("red","blue","cadetblue3","green",NA)


### Third figure - model response curves (highlight SST differences GAM & BRT across PA types, roaddist for elephants?). Full models in supplement
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
m1<-bwhaleGAMM_back
m2<-bwhaleGAMM_buff
m3<-bwhaleGAMM_CRW
m4<-bwhaleGAMM_rev

pdf(file = paste0("F3-BWhale_PartialPlots_",Sys.Date(),".pdf"))
par(mfrow=c(2,2))

plot(m1$gam,select=1, shade=F, shade.col=add.alpha("red",0.2),col=add.alpha("red",0.4),xlim=c(5,30),ylim=c(-6,4),ylab="Partial Response",xlab="Sea Surface Temperature")
par(new=T)
plot(m2$gam,select=1, shade=F, shade.col=add.alpha("cyan",0.2),col=add.alpha("cyan",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")
par(new=T)
plot(m3$gam,select=1, shade=F, shade.col=add.alpha("chartreuse4",0.2),col=add.alpha("chartreuse4",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")
par(new=T)
plot(m4$gam,select=1, shade=F, shade.col=add.alpha("blue",0.2),col=add.alpha("blue",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")

plot(m1$gam,select=6, shade=F, shade.col=add.alpha("red",0.2),col=add.alpha("red",0.4),xlim=c(-6000,0),ylim=c(-5,5),ylab="Partial Response",xlab="Bathymetry")
par(new=T)
plot(m2$gam,select=6, shade=F, shade.col=add.alpha("cyan",0.2),col=add.alpha("cyan",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")
par(new=T)
plot(m3$gam,select=6, shade=F, shade.col=add.alpha("chartreuse4",0.2),col=add.alpha("chartreuse4",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")
par(new=T)
plot(m4$gam,select=6, shade=F, shade.col=add.alpha("blue",0.2),col=add.alpha("blue",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")

plotdata4 <- plot(bwhaleBRT.lr005.rev, return.grid=TRUE, i.var=6)
plotdata3 <- plot(bwhaleBRT.lr005.CRW, return.grid=TRUE, i.var=6)
plotdata2 <- plot(bwhaleBRT.lr005.buff, return.grid=TRUE, i.var=6)
plotdata1 <- plot(bwhaleBRT.lr005.back, return.grid=TRUE, i.var=6)

plot(plotdata4,type="l",lwd=2,col=add.alpha("blue",0.4),xlim=c(5,30),ylim=c(-6,4),ylab="Partial Response",xlab="Sea Surface Temperature")
par(new=T)
plot(plotdata3,type="l",lwd=2,col=add.alpha("chartreuse4",0.4),xlim=c(5,30),ylim=c(-6,4),ylab="",xlab="")
par(new=T)
plot(plotdata2,type="l",lwd=2,col=add.alpha("cyan",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata1,type="l",lwd=2,col=add.alpha("red",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")

plotdata4 <- plot(bwhaleBRT.lr005.rev, return.grid=TRUE, i.var=1)
plotdata3 <- plot(bwhaleBRT.lr005.CRW, return.grid=TRUE, i.var=1)
plotdata2 <- plot(bwhaleBRT.lr005.buff, return.grid=TRUE, i.var=1)
plotdata1 <- plot(bwhaleBRT.lr005.back, return.grid=TRUE, i.var=1)

plot(plotdata4,type="l",lwd=2,col=add.alpha("blue",0.4),xlim=c(-6000,0),ylim=c(-5,5),ylab="",xlab="")
par(new=T)
plot(plotdata3,type="l",lwd=2,col=add.alpha("chartreuse4",0.4),xlim=c(-6000,0),ylim=c(-5,5),ylab="Partial Response",xlab="Bathymetry")
par(new=T)
plot(plotdata2,type="l",lwd=2,col=add.alpha("cyan",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata1,type="l",lwd=2,col=add.alpha("red",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")
dev.off()

### Elephant partial plots!
m1<-ephantGAMM_back
m2<-ephantGAMM_buff
m3<-ephantGAMM_CRW
m4<-ephantGAMM_rev

pdf(file = paste0("F3-Ephant_PartialPlots_",Sys.Date(),".pdf"))
par(mfrow=c(2,2))

plot(m1$gam,select=1, shade=F, shade.col=add.alpha("red",0.2),col=add.alpha("red",0.4),xlim=c(0,30000),ylim=c(-100,0),ylab="Partial Response",xlab="Distance to Road")
legend("bottomleft", legend=c("Background", "Buffer", "CRW","Reverse CRW"), col=c(add.alpha("red",0.4), add.alpha("cyan",0.4),add.alpha("chartreuse4",0.4),add.alpha("blue",0.4)), lty=1, cex=0.8)
par(new=T)
plot(m2$gam,select=1, shade=F, shade.col=add.alpha("cyan",0.2),col=add.alpha("cyan",0.4),xlim=c(0,30000),ylim=c(-100,0),axis=F,ylab="",xlab="")
par(new=T)
plot(m3$gam,select=1, shade=F, shade.col=add.alpha("chartreuse4",0.2),col=add.alpha("chartreuse4",0.4),xlim=c(0,30000),ylim=c(-100,0),axis=F,ylab="",xlab="")
par(new=T)
plot(m4$gam,select=1, shade=F, shade.col=add.alpha("blue",0.2),col=add.alpha("blue",0.4),xlim=c(0,30000),ylim=c(-100,0),axis=F,ylab="",xlab="")

plot(m1$gam,select=2, shade=F, shade.col=add.alpha("red",0.2),col=add.alpha("red",0.4),xlim=c(0,5000),ylim=c(-20,2),ylab="Partial Response",xlab="NDVI")
par(new=T)
plot(m2$gam,select=2, shade=F, shade.col=add.alpha("cyan",0.2),col=add.alpha("cyan",0.4),xlim=c(0,5000),ylim=c(-20,2),axis=F,ylab="",xlab="")
par(new=T)
plot(m3$gam,select=2, shade=F, shade.col=add.alpha("chartreuse4",0.2),col=add.alpha("chartreuse4",0.4),xlim=c(0,5000),ylim=c(-20,2),axis=F,ylab="",xlab="")
par(new=T)
plot(m4$gam,select=2, shade=F, shade.col=add.alpha("blue",0.2),col=add.alpha("blue",0.4),xlim=c(0,5000),ylim=c(-20,2),axis=F,ylab="",xlab="")

plotdata4 <- plot(ephantBRT.lr005.rev, return.grid=TRUE, i.var=1)
plotdata3 <- plot(ephantBRT.lr005.CRW, return.grid=TRUE, i.var=1)
plotdata2 <- plot(ephantBRT.lr005.buff, return.grid=TRUE, i.var=1)
plotdata1 <- plot(ephantBRT.lr005.back, return.grid=TRUE, i.var=1)

plot(plotdata4,type="l",lwd=2,col=add.alpha("blue",0.4),xlim=c(0,30000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata3,type="l",lwd=2,col=add.alpha("chartreuse4",0.4),xlim=c(0,30000),ylim=c(-3.5,3),ylab="Partial Response",xlab="Distance to Road")
par(new=T)
plot(plotdata2,type="l",lwd=2,col=add.alpha("cyan",0.4),xlim=c(0,30000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata1,type="l",lwd=2,col=add.alpha("red",0.4),xlim=c(0,30000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")

plotdata3 <- plot(ephantBRT.lr005.CRW, return.grid=TRUE, i.var=2)
plotdata2 <- plot(ephantBRT.lr005.buff, return.grid=TRUE, i.var=2)
plotdata1 <- plot(ephantBRT.lr005.back, return.grid=TRUE, i.var=2)

plot(plotdata2,type="l",lwd=2,col=add.alpha("blue",0.4),xlim=c(0,5000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata3,type="l",lwd=2,col=add.alpha("chartreuse4",0.4),xlim=c(0,5000),ylim=c(-3.5,3),ylab="Partial Response",xlab="NDVI")
par(new=T)
plot(plotdata2,type="l",lwd=2,col=add.alpha("cyan",0.4),xlim=c(0,5000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata1,type="l",lwd=2,col=add.alpha("red",0.4),xlim=c(0,5000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")
dev.off()

pdf(file = paste0("F3-8-panel_PartialPlots_",Sys.Date(),".pdf"))
par(mfrow=c(2,4))

m1<-bwhaleGAMM_back
m2<-bwhaleGAMM_buff
m3<-bwhaleGAMM_CRW
m4<-bwhaleGAMM_rev

plot(m1$gam,select=1, shade=F, shade.col=add.alpha("red",0.2),col=add.alpha("red",0.4),xlim=c(5,30),ylim=c(-6,4),ylab="Partial Response",xlab="Sea Surface Temperature")
par(new=T)
plot(m2$gam,select=1, shade=F, shade.col=add.alpha("cyan",0.2),col=add.alpha("cyan",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")
par(new=T)
plot(m3$gam,select=1, shade=F, shade.col=add.alpha("chartreuse4",0.2),col=add.alpha("chartreuse4",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")
par(new=T)
plot(m4$gam,select=1, shade=F, shade.col=add.alpha("blue",0.2),col=add.alpha("blue",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")

plot(m1$gam,select=6, shade=F, shade.col=add.alpha("red",0.2),col=add.alpha("red",0.4),xlim=c(-6000,0),ylim=c(-5,5),ylab="Partial Response",xlab="Bathymetry")
par(new=T)
plot(m2$gam,select=6, shade=F, shade.col=add.alpha("cyan",0.2),col=add.alpha("cyan",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")
par(new=T)
plot(m3$gam,select=6, shade=F, shade.col=add.alpha("chartreuse4",0.2),col=add.alpha("chartreuse4",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")
par(new=T)
plot(m4$gam,select=6, shade=F, shade.col=add.alpha("blue",0.2),col=add.alpha("blue",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")

plotdata4 <- plot(bwhaleBRT.lr005.rev, return.grid=TRUE, i.var=5)
plotdata3 <- plot(bwhaleBRT.lr005.CRW, return.grid=TRUE, i.var=5)
plotdata2 <- plot(bwhaleBRT.lr005.buff, return.grid=TRUE, i.var=5)
plotdata1 <- plot(bwhaleBRT.lr005.back, return.grid=TRUE, i.var=5)

plot(plotdata4,type="l",lwd=2,col=add.alpha("blue",0.4),xlim=c(5,30),ylim=c(-6,4),ylab="Partial Response",xlab="Sea Surface Temperature")
par(new=T)
plot(plotdata3,type="l",lwd=2,col=add.alpha("chartreuse4",0.4),xlim=c(5,30),ylim=c(-6,4),ylab="",xlab="")
par(new=T)
plot(plotdata2,type="l",lwd=2,col=add.alpha("cyan",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata1,type="l",lwd=2,col=add.alpha("red",0.4),xlim=c(5,30),ylim=c(-6,4),axis=F,ylab="",xlab="")

plotdata4 <- plot(bwhaleBRT.lr005.rev, return.grid=TRUE, i.var=1)
plotdata3 <- plot(bwhaleBRT.lr005.CRW, return.grid=TRUE, i.var=1)
plotdata2 <- plot(bwhaleBRT.lr005.buff, return.grid=TRUE, i.var=1)
plotdata1 <- plot(bwhaleBRT.lr005.back, return.grid=TRUE, i.var=1)

plot(plotdata4,type="l",lwd=2,col=add.alpha("blue",0.4),xlim=c(-6000,0),ylim=c(-5,5),ylab="",xlab="")
legend("topright", legend=c("Background", "Buffer", "CRW","Reverse CRW"), col=c(add.alpha("red",0.4), add.alpha("cyan",0.4),add.alpha("chartreuse4",0.4),add.alpha("blue",0.4)), lty=1, cex=0.8)
par(new=T)
plot(plotdata3,type="l",lwd=2,col=add.alpha("chartreuse4",0.4),xlim=c(-6000,0),ylim=c(-5,5),ylab="Partial Response",xlab="Bathymetry")
par(new=T)
plot(plotdata2,type="l",lwd=2,col=add.alpha("cyan",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata1,type="l",lwd=2,col=add.alpha("red",0.4),xlim=c(-6000,0),ylim=c(-5,5),axis=F,ylab="",xlab="")


### Now include elephant partial plots!
m1<-ephantGAMM_back
m2<-ephantGAMM_buff
m3<-ephantGAMM_CRW
m4<-ephantGAMM_rev

plot(m1$gam,select=1, shade=F, shade.col=add.alpha("red",0.2),col=add.alpha("red",0.4),xlim=c(0,30000),ylim=c(-100,0),ylab="Partial Response",xlab="Distance to Road")
par(new=T)
plot(m2$gam,select=1, shade=F, shade.col=add.alpha("cyan",0.2),col=add.alpha("cyan",0.4),xlim=c(0,30000),ylim=c(-100,0),axis=F,ylab="",xlab="")
par(new=T)
plot(m3$gam,select=1, shade=F, shade.col=add.alpha("chartreuse4",0.2),col=add.alpha("chartreuse4",0.4),xlim=c(0,30000),ylim=c(-100,0),axis=F,ylab="",xlab="")
par(new=T)
plot(m4$gam,select=1, shade=F, shade.col=add.alpha("blue",0.2),col=add.alpha("blue",0.4),xlim=c(0,30000),ylim=c(-100,0),axis=F,ylab="",xlab="")

plot(m1$gam,select=2, shade=F, shade.col=add.alpha("red",0.2),col=add.alpha("red",0.4),xlim=c(0,5000),ylim=c(-20,2),ylab="Partial Response",xlab="NDVI")
par(new=T)
plot(m2$gam,select=2, shade=F, shade.col=add.alpha("cyan",0.2),col=add.alpha("cyan",0.4),xlim=c(0,5000),ylim=c(-20,2),axis=F,ylab="",xlab="")
par(new=T)
plot(m3$gam,select=2, shade=F, shade.col=add.alpha("chartreuse4",0.2),col=add.alpha("chartreuse4",0.4),xlim=c(0,5000),ylim=c(-20,2),axis=F,ylab="",xlab="")
par(new=T)
plot(m4$gam,select=2, shade=F, shade.col=add.alpha("blue",0.2),col=add.alpha("blue",0.4),xlim=c(0,5000),ylim=c(-20,2),axis=F,ylab="",xlab="")

plotdata4 <- plot(ephantBRT.lr005.rev, return.grid=TRUE, i.var=1)
plotdata3 <- plot(ephantBRT.lr005.CRW, return.grid=TRUE, i.var=1)
plotdata2 <- plot(ephantBRT.lr005.buff, return.grid=TRUE, i.var=1)
plotdata1 <- plot(ephantBRT.lr005.back, return.grid=TRUE, i.var=1)

plot(plotdata4,type="l",lwd=2,col=add.alpha("blue",0.4),xlim=c(0,30000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata3,type="l",lwd=2,col=add.alpha("chartreuse4",0.4),xlim=c(0,30000),ylim=c(-3.5,3),ylab="Partial Response",xlab="Distance to Road")
par(new=T)
plot(plotdata2,type="l",lwd=2,col=add.alpha("cyan",0.4),xlim=c(0,30000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata1,type="l",lwd=2,col=add.alpha("red",0.4),xlim=c(0,30000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")

plotdata3 <- plot(ephantBRT.lr005.CRW, return.grid=TRUE, i.var=2)
plotdata2 <- plot(ephantBRT.lr005.buff, return.grid=TRUE, i.var=2)
plotdata1 <- plot(ephantBRT.lr005.back, return.grid=TRUE, i.var=2)

plot(plotdata2,type="l",lwd=2,col=add.alpha("blue",0.4),xlim=c(0,5000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata3,type="l",lwd=2,col=add.alpha("chartreuse4",0.4),xlim=c(0,5000),ylim=c(-3.5,3),ylab="Partial Response",xlab="NDVI")
par(new=T)
plot(plotdata2,type="l",lwd=2,col=add.alpha("cyan",0.4),xlim=c(0,5000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")
par(new=T)
plot(plotdata1,type="l",lwd=2,col=add.alpha("red",0.4),xlim=c(0,5000),ylim=c(-3.5,3),axis=F,ylab="",xlab="")

dev.off()  


### Fourth figure - model predictions 
## blue whales ---------------------------------------------------------------------------------
eez=st_read(paste0(dropboxdir,"Data/CA_shapefiles/World_EEZ_v10_20180221/eez_boundaries_v10.shp"))%>% st_simplify(preserveTopology=TRUE, dTolerance = .2) 
map.world = map_data(map="world")
testt=map.world %>% filter(long<=180)
extent1<-c(-140,-115,32,45)

load(paste0(dropboxdir,"Data/bwhalemodels_reduce_space_1Chl.RData"))
load(paste0(dropboxdir,"Data/bwhalemodels_BRTeqGAMM.RData"))
rasterDir=paste0(dropboxdir,"Data/RasterData")
r2006.08.01=readRDS(paste0(rasterDir,"/2006-08-01/rasterdata.RDS"))# %>% 
#as.data.frame()%>% mutate(RN=0) #%>% rename("EKE"="eke") 
r2012.03.01=readRDS(paste0(rasterDir,"/2012-03-01/rasterdata.RDS"))
r2012.06.01=readRDS(paste0(rasterDir,"/2012-06-01/rasterdata.RDS"))
r2012.09.01=readRDS(paste0(rasterDir,"/2012-09-01/rasterdata.RDS"))
r2012.12.01=readRDS(paste0(rasterDir,"/2012-12-01/rasterdata.RDS"))

# function
predict_models=function(modtype,patype,rasStack,model,date){
  
  if(modtype=="BRT"){
    pred=predict(rasStack,model,type="response", n.trees=model$gbm.call$best.trees,na.rm=F)
    pred=crop(pred,extent1)
  } else if(modtype=="GAMM"){
    pred=predict(rasStack,model$gam,type="response")
    pred=crop(pred,extent1)
  } else if(modtype=="GLMM"){
    pred=predict(rasStack,model$gam,type="response")
    pred=crop(pred,extent1)
  }
  
  df_map=rasterToPoints(pred)%>% as.data.frame()
  colnames(df_map)=c("rows","cols","value")    
  a=df_map %>% 
    ggplot() + 
    # geom_tile(aes(x = rows, y = cols, fill = ntile(value,100))) + coord_equal()+
    geom_tile(aes(x = rows, y = cols, fill = value),color=NA) + coord_equal()+
    scale_fill_gradientn("Habitat suitability",colours = pals::parula(100),limits = c(0,1),breaks=c(0,.5,1),labels=c("0",".5","1"))+
    scale_x_continuous(expand=c(0,0),limits = c(-140,-115)) +
    scale_y_continuous(expand=c(0,0),limits = c(32,45))+
    geom_map(data=testt,map=testt,aes(map_id=region,x=long,y=lat),fill="darkgrey",color="black")+
    geom_sf(data = eez, color = "black",fill=NA,size=1)+
    theme(legend.position = "bottom",
          # legend.key.size = unit(.9,'lines'),
          # legend.background = NULL,
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          plot.margin = unit(c(0,0,3,0), "mm"))+theme(axis.ticks=element_blank(), 
                                                      panel.background=element_rect(fill = "white"), 
                                                      axis.text.x=element_blank(), axis.text.y=element_blank(),  
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(),
                                                      axis.title.x=element_blank(), axis.title.y=element_blank())+
    ggtitle(glue("{modtype} {patype} {date}"))
  
  pdf(file=glue("{dropboxdir}Output/predictions/{modtype}_{patype}_{date}.pdf"), height = 5, width = 6)
  par(mar=c(0,0,0,0))
  # plot_grid(backplot,buffplot,crwplot,revplot, ncol=2, align="h",axis = "bt")
  grid.arrange(a,ncol=1 )
  dev.off()
}

# brt ####
# brt background
predict_models(modtype = "BRT",patype = "background",rasStack = r2006.08.01,model= bwhaleBRT.lr005.back,date="2006.08.01")
predict_models(modtype = "BRT",patype = "background",rasStack = r2012.06.01,model= bwhaleBRT.lr005.back,date="2012.06.01")
predict_models(modtype = "BRT",patype = "background",rasStack = r2012.03.01,model= bwhaleBRT.lr005.back,date="2012.03.01")
predict_models(modtype = "BRT",patype = "background",rasStack = r2012.09.01,model= bwhaleBRT.lr005.back,date="2012.09.01")
predict_models(modtype = "BRT",patype = "background",rasStack = r2012.12.01,model= bwhaleBRT.lr005.back,date="2012.12.01")

# brt crw
predict_models(modtype = "BRT",patype = "CRW",rasStack = r2006.08.01,model= bwhaleBRT.lr005.CRW,date="2006.08.01")
predict_models(modtype = "BRT",patype = "CRW",rasStack = r2012.06.01,model= bwhaleBRT.lr005.CRW,date="2012.06.01")
predict_models(modtype = "BRT",patype = "CRW",rasStack = r2012.03.01,model= bwhaleBRT.lr005.CRW,date="2012.03.01")
predict_models(modtype = "BRT",patype = "CRW",rasStack = r2012.09.01,model= bwhaleBRT.lr005.CRW,date="2012.09.01")
predict_models(modtype = "BRT",patype = "CRW",rasStack = r2012.12.01,model= bwhaleBRT.lr005.CRW,date="2012.12.01")

# brt rev
predict_models(modtype = "BRT",patype = "rev.CRW",rasStack = r2006.08.01,model= bwhaleBRT.lr005.rev,date="2006.08.01")
predict_models(modtype = "BRT",patype = "rev.CRW",rasStack = r2012.06.01,model= bwhaleBRT.lr005.rev,date="2012.06.01")
predict_models(modtype = "BRT",patype = "rev.CRW",rasStack = r2012.03.01,model= bwhaleBRT.lr005.rev,date="2012.03.01")
predict_models(modtype = "BRT",patype = "rev.CRW",rasStack = r2012.09.01,model= bwhaleBRT.lr005.rev,date="2012.09.01")
predict_models(modtype = "BRT",patype = "rev.CRW",rasStack = r2012.12.01,model= bwhaleBRT.lr005.rev,date="2012.12.01")

# brt buff
predict_models(modtype = "BRT",patype = "buffer",rasStack = r2006.08.01,model= bwhaleBRT.lr005.buff,date="2006.08.01")
predict_models(modtype = "BRT",patype = "buffer",rasStack = r2012.06.01,model= bwhaleBRT.lr005.buff,date="2012.06.01")
predict_models(modtype = "BRT",patype = "buffer",rasStack = r2012.03.01,model= bwhaleBRT.lr005.buff,date="2012.03.01")
predict_models(modtype = "BRT",patype = "buffer",rasStack = r2012.09.01,model= bwhaleBRT.lr005.buff,date="2012.09.01")
predict_models(modtype = "BRT",patype = "buffer",rasStack = r2012.12.01,model= bwhaleBRT.lr005.buff,date="2012.12.01")

# GAMM ####
# GAMM background
predict_models(modtype = "GAMM",patype = "background",rasStack = r2006.08.01,model= bwhaleGAMM_back,date="2006.08.01")
predict_models(modtype = "GAMM",patype = "background",rasStack = r2012.06.01,model= bwhaleGAMM_back,date="2012.06.01")
predict_models(modtype = "GAMM",patype = "background",rasStack = r2012.03.01,model= bwhaleGAMM_back,date="2012.03.01")
predict_models(modtype = "GAMM",patype = "background",rasStack = r2012.09.01,model= bwhaleGAMM_back,date="2012.09.01")
predict_models(modtype = "GAMM",patype = "background",rasStack = r2012.12.01,model= bwhaleGAMM_back,date="2012.12.01")

# GAMM crw
predict_models(modtype = "GAMM",patype = "CRW",rasStack = r2006.08.01,model= bwhaleGAMM_CRW,date="2006.08.01")
predict_models(modtype = "GAMM",patype = "CRW",rasStack = r2012.06.01,model= bwhaleGAMM_CRW,date="2012.06.01")
predict_models(modtype = "GAMM",patype = "CRW",rasStack = r2012.03.01,model= bwhaleGAMM_CRW,date="2012.03.01")
predict_models(modtype = "GAMM",patype = "CRW",rasStack = r2012.09.01,model= bwhaleGAMM_CRW,date="2012.09.01")
predict_models(modtype = "GAMM",patype = "CRW",rasStack = r2012.12.01,model= bwhaleGAMM_CRW,date="2012.12.01")

# GAMM rev
predict_models(modtype = "GAMM",patype = "rev.CRW",rasStack = r2006.08.01,model= bwhaleGAMM_rev,date="2006.08.01")
predict_models(modtype = "GAMM",patype = "rev.CRW",rasStack = r2012.06.01,model= bwhaleGAMM_rev,date="2012.06.01")
predict_models(modtype = "GAMM",patype = "rev.CRW",rasStack = r2012.03.01,model= bwhaleGAMM_rev,date="2012.03.01")
predict_models(modtype = "GAMM",patype = "rev.CRW",rasStack = r2012.09.01,model= bwhaleGAMM_rev,date="2012.09.01")
predict_models(modtype = "GAMM",patype = "rev.CRW",rasStack = r2012.12.01,model= bwhaleGAMM_rev,date="2012.12.01")

# GAMM buff
predict_models(modtype = "GAMM",patype = "buffer",rasStack = r2006.08.01,model= bwhaleGAMM_buff,date="2006.08.01")
predict_models(modtype = "GAMM",patype = "buffer",rasStack = r2012.06.01,model= bwhaleGAMM_buff,date="2012.06.01")
predict_models(modtype = "GAMM",patype = "buffer",rasStack = r2012.03.01,model= bwhaleGAMM_buff,date="2012.03.01")
predict_models(modtype = "GAMM",patype = "buffer",rasStack = r2012.09.01,model= bwhaleGAMM_buff,date="2012.09.01")
predict_models(modtype = "GAMM",patype = "buffer",rasStack = r2012.12.01,model= bwhaleGAMM_buff,date="2012.12.01")

# GLMM ####
# GLMM background
predict_models(modtype = "GLMM",patype = "background",rasStack = r2006.08.01,model= bwhaleGLMM_back,date="2006.08.01")
predict_models(modtype = "GLMM",patype = "background",rasStack = r2012.06.01,model= bwhaleGLMM_back,date="2012.06.01")
predict_models(modtype = "GLMM",patype = "background",rasStack = r2012.03.01,model= bwhaleGLMM_back,date="2012.03.01")
predict_models(modtype = "GLMM",patype = "background",rasStack = r2012.09.01,model= bwhaleGLMM_back,date="2012.09.01")
predict_models(modtype = "GLMM",patype = "background",rasStack = r2012.12.01,model= bwhaleGLMM_back,date="2012.12.01")

# GLMM crw
predict_models(modtype = "GLMM",patype = "CRW",rasStack = r2006.08.01,model= bwhaleGLMM_CRW,date="2006.08.01")
predict_models(modtype = "GLMM",patype = "CRW",rasStack = r2012.06.01,model= bwhaleGLMM_CRW,date="2012.06.01")
predict_models(modtype = "GLMM",patype = "CRW",rasStack = r2012.03.01,model= bwhaleGLMM_CRW,date="2012.03.01")
predict_models(modtype = "GLMM",patype = "CRW",rasStack = r2012.09.01,model= bwhaleGLMM_CRW,date="2012.09.01")
predict_models(modtype = "GLMM",patype = "CRW",rasStack = r2012.12.01,model= bwhaleGLMM_CRW,date="2012.12.01")

# GLMM rev
predict_models(modtype = "GLMM",patype = "rev.CRW",rasStack = r2006.08.01,model= bwhaleGLMM_rev,date="2006.08.01")
predict_models(modtype = "GLMM",patype = "rev.CRW",rasStack = r2012.06.01,model= bwhaleGLMM_rev,date="2012.06.01")
predict_models(modtype = "GLMM",patype = "rev.CRW",rasStack = r2012.03.01,model= bwhaleGLMM_rev,date="2012.03.01")
predict_models(modtype = "GLMM",patype = "rev.CRW",rasStack = r2012.09.01,model= bwhaleGLMM_rev,date="2012.09.01")
predict_models(modtype = "GLMM",patype = "rev.CRW",rasStack = r2012.12.01,model= bwhaleGLMM_rev,date="2012.12.01")

# GLMM buff
predict_models(modtype = "GLMM",patype = "buffer",rasStack = r2006.08.01,model= bwhaleGLMM_buff,date="2006.08.01")
predict_models(modtype = "GLMM",patype = "buffer",rasStack = r2012.06.01,model= bwhaleGLMM_buff,date="2012.06.01")
predict_models(modtype = "GLMM",patype = "buffer",rasStack = r2012.03.01,model= bwhaleGLMM_buff,date="2012.03.01")
predict_models(modtype = "GLMM",patype = "buffer",rasStack = r2012.09.01,model= bwhaleGLMM_buff,date="2012.09.01")
predict_models(modtype = "GLMM",patype = "buffer",rasStack = r2012.12.01,model= bwhaleGLMM_buff,date="2012.12.01")

## elephants ---------------------------------------------------------------------------------
load("./Data/Predict/EphantPredict_2020-08-13.RData")
fence=st_read("./Data/Etosha\ shapefiles/enp\ fence\ poly.shp")%>% st_simplify(preserveTopology=TRUE, dTolerance = .2) 
fence2=fence %>% summarise(AREA=sum(AREA))

predict_models_elephants=function(modtype,patype,raster){
  pred=mask(raster,fence2)
  
  df_map=rasterToPoints(pred)%>% as.data.frame()
  colnames(df_map)=c("rows","cols","value")    
  a=df_map %>% 
    ggplot() + 
    geom_tile(aes(x = rows, y = cols, fill = value),color=NA) + coord_equal()+
    scale_fill_gradientn("Habitat suitability",colours = pals::parula(100),limits = c(0,1),breaks=c(0,.5,1),labels=c("0",".5","1"))+
    theme(legend.position = "bottom",
          # legend.key.size = unit(.9,'lines'),
          # legend.background = NULL,
          legend.text = element_text(size=10),
          legend.title = element_text(size=12),
          plot.margin = unit(c(0,0,3,0), "mm"))+theme(axis.ticks=element_blank(), 
                                                      panel.background=element_rect(fill = "white"), 
                                                      axis.text.x=element_blank(), axis.text.y=element_blank(),  
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(),
                                                      axis.title.x=element_blank(), axis.title.y=element_blank())+
    ggtitle(glue("{modtype} {patype}"))
  
  pdf(file=glue("./Output/predictions/{modtype}_{patype}_elephant.pdf"), height = 3, width = 5)
  par(mar=c(0,0,0,0))
  grid.arrange(a,ncol=1 )
  dev.off()
}

predict_models_elephants(modtype="GAMM",patype = "background",raster=x_back)
predict_models_elephants(modtype="BRT",patype = "background",raster=x_back_brt)
predict_models_elephants(modtype="GLMM",patype = "background",raster=x_back_GLMM)
predict_models_elephants(modtype="GAMM",patype = "buffer",raster=x_buff)
predict_models_elephants(modtype="BRT",patype = "buffer",raster=x_buff_brt)
predict_models_elephants(modtype="GLMM",patype = "buffer",raster=x_buff_GLMM)
predict_models_elephants(modtype="GAMM",patype = "CRW",raster=x_CRW)
predict_models_elephants(modtype="BRT",patype = "CRW",raster=x_CRW_brt)
predict_models_elephants(modtype="GLMM",patype = "CRW",raster=x_CRW_GLMM)
predict_models_elephants(modtype="GAMM",patype = "rev.CRW",raster=x_rev)
predict_models_elephants(modtype="BRT",patype = "rev.CRW",raster=x_rev_brt)
predict_models_elephants(modtype="GLMM",patype = "rev.CRW",raster=x_rev_GLMM)

### Fifth figure - model distance comparison
#Blue whales
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

tag="final"
summarystatistics<-read_csv(paste0("Output/SummaryStatistics_bwhale_",tag,".csv"))
summarystatistics <- dplyr::select(summarystatistics, -X1)
#write.csv(summarystatistics, file=paste0("Output/SummaryStatistics_bwhale_",tag,".csv"))

### plot bh
summarystatistics$samptype<-c(rep("CRW",3),rep("back",3),rep("buff",3),rep("rev",3))
summarystatistics$modtype<-rep(c("GAMM","GLMM","BRT"),4)
plot_data <- summarystatistics %>% 
  melt(id.vars=names(summarystatistics)[c(1:10,13,15,16)])   
BW_bh_facet<-ggplot(plot_data, aes(x = value,y = AUC)) +
  geom_point(aes(colour = (samptype), shape = (variable))) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  xlab("Bhattacharyya's coefficient") +
  ggtitle("Blue whale model skill as a function of overlap") +
  theme_minimal()

lm.1<-(lm(value~AUC, data=plot_data[which(plot_data$modtype=="BRT"),]))
lm.2<-(lm(value~AUC, data=plot_data[which(plot_data$modtype=="GAMM"),]))
lm.3<-(lm(value~AUC, data=plot_data[which(plot_data$modtype=="GLMM"),]))

R2.BRT<-round(summary(lm.1)$r.squared,3)
pval.BRT <- round(lmp(lm.1),3)
R2.GAMM<-round(summary(lm.2)$r.squared,3)
pval.GAMM <- round(lmp(lm.2),3)
R2.GLMM<-round(summary(lm.3)$r.squared, 3)
pval.GLMM <- round(lmp(lm.3), 3)

BW_bh_facet +  
  facet_grid(. ~ modtype) #+

pdf(file=paste0("AUC v. Bhatt - Blue Whale facet_",Sys.Date(),".pdf"))
BW_bh_facet + facet_grid(. ~ modtype) 
dev.off()

#Elephants
tag="etosha_only"
summarystatistics<-read_csv(paste0("Output/SummaryStatistics_ephant_",tag,".csv"))  #write.csv(summarystatistics, file=paste0("Output/SummaryStatistics_ephant_",tag,".csv"))
summarystatistics <- dplyr::select(summarystatistics, -X1)

summarystatistics$samptype<-c(rep("CRW",3),rep("back",3),rep("buff",3),rep("rev",3))
summarystatistics$modtype<-rep(c("GAMM","GLMM","BRT"),4)
plot_data <- summarystatistics %>% 
  melt(id.vars=names(summarystatistics)[c(1:10,13,15,16)]) 

ephant_bh_facet<-ggplot(plot_data, aes(x = value,y = AUC)) +
  geom_point(aes(colour = (samptype), shape = (variable))) +
  #  geom_smooth(method='lm') +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  xlab("Bhattacharyya's coefficient") +
  ggtitle("Elephant model skill as a function of overlap") +
  theme_minimal()

lm.1<-(lm(value~AUC, data=plot_data[which(plot_data$modtype=="BRT"),]))
lm.2<-(lm(value~AUC, data=plot_data[which(plot_data$modtype=="GAMM"),]))
lm.3<-(lm(value~AUC, data=plot_data[which(plot_data$modtype=="GLMM"),]))

R2.BRT<-round(summary(lm.1)$r.squared,5)
pval.BRT <- round(lmp(lm.1),5)
R2.GAMM<-round(summary(lm.2)$r.squared,5)
pval.GAMM <- round(lmp(lm.2),5)
R2.GLMM<-round(summary(lm.3)$r.squared, 5)
pval.GLMM <- round(lmp(lm.3), 5)
ephant_bh_facet + facet_grid(. ~ modtype)


pdf(file=paste0("AUC v. Bhatt - Elephant facet_",Sys.Date(),".pdf"))
ephant_bh_facet + facet_grid(. ~ modtype)
dev.off()


### Fifth figure - model distance comparison -HEATHER WELCH EDITS
#Blue whales
tag="reduce_space_1Chl"
summarystatistics<-read.csv(paste0("/Users/heatherwelch/Dropbox/PseudoAbsence-MS/Output/SummaryStatistics_bwhale_",tag,".csv"))
summarystatistics <- select(summarystatistics, -X)

#write.csv(summarystatistics, file=paste0("Output/SummaryStatistics_bwhale_",tag,".csv"))

### plot bh
summarystatistics$samptype<-c(rep("CRW",3),rep("back",3),rep("buff",3),rep("rev",3))
summarystatistics$modtype<-rep(c("GAMM","GLMM","BRT"),4)
plot_data <- summarystatistics %>% 
  melt(id.vars=names(summarystatistics)[c(1:4,9,10)])  

## auc
BW_bh_facet<-ggplot(plot_data, aes(x = value,y = AUC)) +
  geom_point(aes(colour = (samptype), shape = (variable))) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  xlab("Bhattacharyya's coefficient") +
  ggtitle("Blue whale model skill as a function of overlap") +
  theme_minimal()

BW_bh_facet + facet_grid(. ~ modtype)


pdf(file="AUC v. Bhatt - Blue Whale facet.pdf")
BW_bh_facet + facet_grid(. ~ modtype)
dev.off()

## r-sq
BW_bh_facet<-ggplot(plot_data, aes(x = value,y = R_squared)) +
  geom_point(aes(colour = (samptype), shape = (variable))) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  xlab("Bhattacharyya's coefficient") +
  ggtitle("Blue whale model skill as a function of overlap") +
  theme_minimal()

BW_bh_facet + facet_grid(. ~ modtype)


pdf(file="r2 v. Bhatt - Blue Whale facet.pdf")
BW_bh_facet + facet_grid(. ~ modtype)
dev.off()

## AIC
BW_bh_facet<-ggplot(plot_data, aes(x = value,y = AIC)) +
  geom_point(aes(colour = (samptype), shape = (variable))) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  xlab("Bhattacharyya's coefficient") +
  ggtitle("Blue whale model skill as a function of overlap") +
  theme_minimal()

BW_bh_facet + facet_grid(. ~ modtype)

pdf(file="AIC v. Bhatt - Blue Whale facet.pdf")
BW_bh_facet + facet_grid(. ~ modtype)
dev.off()


#Elephants
tag="reduced_space"
summarystatistics<-read.csv(paste0("/Users/heatherwelch/Dropbox/PseudoAbsence-MS/Output/SummaryStatistics_ephant_",tag,".csv"))  #write.csv(summarystatistics, file=paste0("Output/SummaryStatistics_ephant_",tag,".csv"))
summarystatistics <- select(summarystatistics, -X)

summarystatistics$samptype<-c(rep("CRW",3),rep("back",3),rep("buff",3),rep("rev",3))
summarystatistics$modtype<-rep(c("GAMM","GLMM","BRT"),4)
plot_data <- summarystatistics %>% 
  melt(id.vars=names(summarystatistics)[c(1:4,9,10)]) 

## AUC
ephant_bh_facet<-ggplot(plot_data, aes(x = value,y = AUC)) +
  geom_point(aes(colour = (samptype), shape = (variable))) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  xlab("Bhattacharyya's coefficient") +
  ggtitle("Elephant model skill as a function of overlap") +
  theme_minimal()+ facet_grid( ~ modtype)

pdf(file="AUC v. Bhatt - Elephant facet.pdf")
ephant_bh_facet
dev.off()

a=plot_data %>% rename(Bhatt=value) %>% gather(metric,value,-c(Model_Name,R_squared,AIC,samptype,modtype,variable))
ggplot(a,aes(x=samptype,y=value,fill=metric))+geom_boxplot()+facet_wrap(~modtype)
ggplot(plot_data)+geom_boxplot(aes(x=samptype,y=value))+
  geom_point(aes(x=samptype,y=AUC))+
  facet_wrap(~modtype)

## R2
ephant_bh_facet<-ggplot(plot_data, aes(x = value,y = R_squared)) +
  geom_point(aes(colour = (samptype), shape = (variable))) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  xlab("Bhattacharyya's coefficient") +
  ggtitle("Elephant model skill as a function of overlap") +
  theme_minimal()+ facet_grid( ~ modtype)

pdf(file="r2 v. Bhatt - Elephant facet.pdf")
ephant_bh_facet
dev.off()

## AIC
ephant_bh_facet<-ggplot(plot_data, aes(x = value,y = AIC)) +
  geom_point(aes(colour = (samptype), shape = (variable))) +
  geom_smooth(method='lm', se=FALSE, colour="black") +
  xlab("Bhattacharyya's coefficient") +
  ggtitle("Elephant model skill as a function of overlap") +
  theme_minimal()+ facet_grid( ~ modtype)

pdf(file="AIC v. Bhatt - Elephant facet.pdf")
ephant_bh_facet
dev.off()