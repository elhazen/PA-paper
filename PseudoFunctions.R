require(adehabitatLT) # help(package='adehabitat') # help.search('angle',package='adehabitat')
require(maps)       # for map.where
require(mapdata)    # for worldHires
require(sp)
require(maptools)
require(raster)

createCRW<-function(tags, tagid=tagid, n.sim = 100, reverse = FALSE){
  
  #for (tagid in unique(tags$ptt)[14:17]){    #starting with 66884
  # clear sim.alltags if exists in workspace
  #if (exists('sim.alltags')) rm(sim.alltags)
  
   
  print(tagid)
  tag = tags[which(tags$tag==tagid),c('long', 'lat','dTime')]
  tag = arrange(tag, dTime)
  
  out.alltags.csv = sprintf('%s/crw_sim_%s.csv', out.dir,tagid)
  
  if (reverse) {
    out.alltags.csv = sprintf('%s/crw_sim_reverse_%d.csv', out.dir,tagid)
    tag = tag[seq(dim(tag)[1],1),]  
  }
  
  # remove non-unique dates
  dupes = which(duplicated(tag$dTime))
  if (length(dupes)!=0){
    tag = tag[-dupes,]
  }
  
  # remove NAs
  tag = tag[!is.na(tag$lat),]
  
  # Creation of an object of class "ltraj"
  tr = as.ltraj(cbind(tag$lon, tag$lat), date=tag$dTime, id=tagid)
  tr1 = tr[[1]]
  
  head(tr1)
  x11()
  plot(tr)
  
  # get indices to rows of trajectory that are not NA for sampling in simulation
  i.tr1.nona = which(!is.na(tr1[['dist']]) & !is.na(tr1[['rel.angle']]))
  #tr1<-tr1[i.tr1.nona,]
  # setup map per tag
  #mex.map = map('worldHires', 'Gabon', fill=T, col='transparent', plot=F)
  CC.map = maps::map('worldHires', fill=T, col='transparent', plot=F, ylim = c(0,50))
  
  CC.IDs = sapply(strsplit(CC.map$names, ":"), function(x) x[1])
  CC.sp = map2SpatialPolygons(CC.map, IDs=CC.IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # setup data frame for simulation
  n.tag = nrow(tag)
  #n.sim = n.tag # or change this to 1000, or n.tag+1000
  sim = data.frame(x=numeric(n.sim), y=numeric(n.sim), t=as.POSIXct(rep(NA,n.sim)), tagid=rep(tagid,n.sim), iteration=1:n.sim)
  #sim = data.frame(x=numeric(n.sim), y=numeric(n.sim), t=as.POSIXct(rep(NA,n.sim)),iteration=1:n.sim)
  
  if(n.tag > 0){
    for (k in 1:n.sim){
      
      print(sprintf("  k'th simulation: %d",k))
      sim = data.frame(x=numeric(n.tag), y=numeric(n.tag), t=as.POSIXct(rep(NA,n.tag)), flag=rep(NA,n.tag))     # new location
      
      # create column to keep track of sim iteration number
      sim$iteration = rep(k, n.tag)
      
      # populate initial location
      sim[1,'x'] = tag$lon[1]
      sim[1,'y'] = tag$lat[1]
      sim[1,'t'] = tr1$date[1]
      angle      = tr1[2,'abs.angle'] # starting angle
      dtime      = tr1[1,'date']      # starting time
      
      #debug plot check
      #par(mar=c(3.5, 2, 1, 1)) # margin = c(bottom, left, top, right)
      maps::map('worldHires', xlim=c(-140,-100), ylim=c(15,50))         # Bluewhale-centric projection
      #  map('worldHires', xlim=c(-100,5), ylim=c(18,50))         # atlantic-centric projection
      map.axes()
      lines(tag$lon, tag$lat,col='grey')                          # plot grey lines for original track
      points(sim[1,'x'],sim[1,'y'], col='blue', pch=2, cex=2)     # plot initial point like ltraj symbology
      title(main = sprintf("%d 'th CRW Simulation for tagid %s", k,tagid))
      
      
      # for (j in 2:n.sim){ # j = 2
      for (j in 2:n.tag){ # j = 2
        on.land = T # start with assumption on land to force getting new location
        while (on.land == T){
          # get random index of trajectory row, outside NA values
          i = sample(i.tr1.nona, 1)
          
          # get distance and angle from original data
          dist  = tr1[i,'dist']
          angle = angle + tr1[i,'rel.angle']
          
          # force westward if outside of gulf
          #	  if (tr1[i,"x"]<-80) {angle<-runif(1,3.14,5.05)}
          
          # calculate a new x and y
          x = sim[j-1,'x'] + dist * cos(angle)
          y = sim[j-1,'y'] + dist * sin(angle)
          
          # check if on land
          pt = SpatialPoints(matrix(c(x,y),nrow=1), proj4string=CRS("+proj=longlat +datum=WGS84"))
          place = over(pt, CC.sp)
          
          #debug plot check
          if (is.na(place)){
            # assign valid on.water location to sim data frame
            sim[j,'x'] = x
            sim[j,'y'] = y
            sim[j,'t'] = sim[j-1,'t'] + tr1[i,'dt']
            
            # update plot with points
            lines(sim[(j-1):j, c('x','y')], col='black')
            points(x, y, col='blue', pch=20)
            
            # set to F so bumps out of while loop and to next simulated location
            on.land = F
          } else {
            points(x, y, col='red', pch=20)
          }
        } # end while (on.land == T){
        
      } # end for (j in 2:n.sim)
      
      #### Calculate flag for quality of pseudo-track
      tagstart<-rbind(tag[1,],tag[2,],tag[n.tag,])
      tagsim<-rbind(sim[1,],sim[2,],sim[n.tag,])
      tagstart$lat[2]<-tagstart$lat[1]-0.01
      tagsim$y[2]<-tagsim$y[1]-0.01
      trstart = as.ltraj(cbind(tagstart$lon, tagstart$lat), date=tagstart$dTime,id=tagid)
      trsim = as.ltraj(cbind(tagsim$x, tagsim$y), date=tagsim$t,id=tagid)
      # trstart[[1]]$rel.angle[2]*180/pi
      distdiff<-abs(trstart[[1]]$dist[2]-trsim[[1]]$dist[2])
      angdiff<-abs(trstart[[1]]$rel.angle[2]*180/pi-trsim[[1]]$rel.angle[2]*180/pi)
      sim$flag<-distdiff/trstart[[1]]$dist[2]*3+angdiff/45
      #if (distdiff/trstart[[1]]$dist[2]<(1/3)) sim$flag<-1  
      #if (distdiff/trstart[[1]]$dist[2]<(2/3)) sim$flag<-2
      #if (distdiff/trstart[[1]]$dist[2]<(3/3)) sim$flag<-3
      #if (angdiff < 45 | angdiff >315) sim$flag<-sim$flag
      #if (angdiff < 90 & angdiff >45) sim$flag<-sim$flag+1
      #if (angdiff > 90 & angdiff <135) sim$flag<-sim$flag+2
      #if (angdiff < 180 & angdiff >135) sim$flag<-sim$flag+3
      #if (angdiff > 180 & angdiff <225) sim$flag<-sim$flag+3
      #if (angdiff < 270 & angdiff >225) sim$flag<-sim$flag+2
      #if (angdiff < 315 & angdiff >270) sim$flag<-sim$flag+1
      
      #  head(trstart[[1]])
      #  head(trsim[[1]])
      
      
      # append full tag simulation to sim.alltags
      if (exists('sim.alltags')){
        sim.alltags = rbind(sim.alltags, sim)
      } else {
        sim.alltags = sim
      }
      
    } # end for (k in 2:n.sim)
  } # end if n.tag > 0
  # write out csv per tag
  #out.csv = sprintf('%s/crw_sim_%d.csv', out.dir, tagid)
  #write.csv(sim, file=out.csv, row.names=F)
  
  # write out
  write.csv(sim.alltags, out.alltags.csv, row.names=F)      #OLD (commented out 4ssm)
  
  
  ################################################################
  # output map
  out.png = sprintf('%s/crw_sim_%s.png', out.dir, tagid); res = 72 # resolution in dots per inch
  dev.print(device=png, file=out.png, width=8*res, height=8*res, res=res)

  return(sim.alltags)
} #end function

create.land.CRW<-function(tags, tagid=tagid, n.sim = 100, reverse = FALSE){
  
  #for (tagid in unique(tags$ptt)[14:17]){    #starting with 66884
  # clear sim.alltags if exists in workspace
  #if (exists('sim.alltags')) rm(sim.alltags)
  
  
  print(tagid)
  tag = tags[which(tags$id==tagid),c('long', 'lat','dTime')]
  
  out.alltags.csv = sprintf('%s/crw_sim_%s.csv', out.dir,tagid)
  
    if (reverse) {
    out.alltags.csv = sprintf('%s/crw_sim_reverse_%s.csv', out.dir,tagid)
    tag = tag[seq(dim(tag)[1],1),]  
  }
  
  # remove non-unique dates
  dupes = which(duplicated(tag$dTime))
  if (length(dupes)!=0){
    tag = tag[-dupes,]
  }
  
  # remove NAs
  tag = tag[!is.na(tag$lat),]
  
  # Creation of an object of class "ltraj"
  tr = as.ltraj(cbind(tag$lon, tag$lat), date=tag$dTime, id=tagid)
  tr1 = tr[[1]]
  
  head(tr1)
  x11()
  plot(tr)
  
  # get indices to rows of trajectory that are not NA for sampling in simulation
  i.tr1.nona = which(!is.na(tr1[['dist']]) & !is.na(tr1[['rel.angle']]))
  #tr1<-tr1[i.tr1.nona,]
  # setup map per tag
  #mex.map = map('worldHires', 'Gabon', fill=T, col='transparent', plot=F)
  CC.map = maps::map('worldHires', fill=T, col='transparent', plot=F, ylim = c(0,50))
  
  CC.IDs = sapply(strsplit(CC.map$names, ":"), function(x) x[1])
  CC.sp = map2SpatialPolygons(CC.map, IDs=CC.IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # setup data frame for simulation
  n.tag = nrow(tag)
  #n.sim = n.tag # or change this to 1000, or n.tag+1000
  sim = data.frame(x=numeric(n.sim), y=numeric(n.sim), t=as.POSIXct(rep(NA,n.sim)), tagid=rep(tagid,n.sim), iteration=1:n.sim)
  #sim = data.frame(x=numeric(n.sim), y=numeric(n.sim), t=as.POSIXct(rep(NA,n.sim)),iteration=1:n.sim)
  
  if(n.tag > 0){
    for (k in 1:n.sim){
      
      print(sprintf("  k'th simulation: %d",k))
      sim = data.frame(x=numeric(n.tag), y=numeric(n.tag), t=as.POSIXct(rep(NA,n.tag)), flag=rep(NA,n.tag))     # new location
      
      # create column to keep track of sim iteration number
      sim$iteration = rep(k, n.tag)
      
      # populate initial location
      sim[1,'x'] = tag$lon[1]
      sim[1,'y'] = tag$lat[1]
      sim[1,'t'] = tr1$date[1]
      angle      = tr1[2,'abs.angle'] # starting angle
      dtime      = tr1[1,'date']      # starting time
      
      #debug plot check
      #par(mar=c(3.5, 2, 1, 1)) # margin = c(bottom, left, top, right)
      maps::map('worldHires', xlim=c(14,18), ylim=c(-20,-15))         # Bluewhale-centric projection
      #  map('worldHires', xlim=c(-100,5), ylim=c(18,50))         # atlantic-centric projection
      map.axes()
      lines(tag$lon, tag$lat,col='grey')                          # plot grey lines for original track
      points(sim[1,'x'],sim[1,'y'], col='blue', pch=2, cex=2)     # plot initial point like ltraj symbology
      title(main = sprintf("%d 'th CRW Simulation for tagid %s", k,tagid))
      
      
      # for (j in 2:n.sim){ # j = 2
      for (j in 2:n.tag){ # j = 2
        on.land = T # start with assumption on land to force getting new location
        while (on.land == T){
          # get random index of trajectory row, outside NA values
          i = sample(i.tr1.nona, 1)
          
          # get distance and angle from original data
          dist  = tr1[i,'dist']
          angle = angle + tr1[i,'rel.angle']
          
          # force westward if outside of gulf
          #	  if (tr1[i,"x"]<-80) {angle<-runif(1,3.14,5.05)}
          
          # calculate a new x and y
          x = sim[j-1,'x'] + dist * cos(angle)
          y = sim[j-1,'y'] + dist * sin(angle)
          
          # check if on land
          pt = SpatialPoints(matrix(c(x,y),nrow=1))
          projection(pt)<-projection(CC.sp)
          place = over(pt, CC.sp)
          
          #debug plot check
            # assign valid on.water location to sim data frame
            sim[j,'x'] = x
            sim[j,'y'] = y
            sim[j,'t'] = sim[j-1,'t'] + tr1[i,'dt']
            
            # update plot with points
            lines(sim[(j-1):j, c('x','y')], col='black')
            points(x, y, col='blue', pch=20)
            
            # set to F so bumps out of while loop and to next simulated location
            on.land = F
        } # end while (on.land == T){
        
      } # end for (j in 2:n.sim)
      
      #### Calculate flag for quality of pseudo-track
      tagstart<-rbind(tag[1,],tag[2,],tag[n.tag,])
      tagsim<-rbind(sim[1,],sim[2,],sim[n.tag,])
      tagstart$lat[2]<-tagstart$lat[1]-0.01
      tagsim$y[2]<-tagsim$y[1]-0.01
      trstart = as.ltraj(cbind(tagstart$lon, tagstart$lat), date=tagstart$dTime,id=tagid)
      trsim = as.ltraj(cbind(tagsim$x, tagsim$y), date=tagsim$t,id=tagid)
      # trstart[[1]]$rel.angle[2]*180/pi
      distdiff<-abs(trstart[[1]]$dist[2]-trsim[[1]]$dist[2])
      angdiff<-abs(trstart[[1]]$rel.angle[2]*180/pi-trsim[[1]]$rel.angle[2]*180/pi)
      sim$flag<-distdiff/trstart[[1]]$dist[2]*3+angdiff/45
      #if (distdiff/trstart[[1]]$dist[2]<(1/3)) sim$flag<-1  
      #if (distdiff/trstart[[1]]$dist[2]<(2/3)) sim$flag<-2
      #if (distdiff/trstart[[1]]$dist[2]<(3/3)) sim$flag<-3
      #if (angdiff < 45 | angdiff >315) sim$flag<-sim$flag
      #if (angdiff < 90 & angdiff >45) sim$flag<-sim$flag+1
      #if (angdiff > 90 & angdiff <135) sim$flag<-sim$flag+2
      #if (angdiff < 180 & angdiff >135) sim$flag<-sim$flag+3
      #if (angdiff > 180 & angdiff <225) sim$flag<-sim$flag+3
      #if (angdiff < 270 & angdiff >225) sim$flag<-sim$flag+2
      #if (angdiff < 315 & angdiff >270) sim$flag<-sim$flag+1
      
      #  head(trstart[[1]])
      #  head(trsim[[1]])
      
      
      # append full tag simulation to sim.alltags
      if (exists('sim.alltags')){
        sim.alltags = rbind(sim.alltags, sim)
      } else {
        sim.alltags = sim
      }
      
    } # end for (k in 2:n.sim)
  } # end if n.tag > 0
  # write out csv per tag
  #out.csv = sprintf('%s/crw_sim_%d.csv', out.dir, tagid)
  #write.csv(sim, file=out.csv, row.names=F)
  
  # write out
  write.csv(sim.alltags, out.alltags.csv, row.names=F)      #OLD (commented out 4ssm)
  
  ################################################################
  # output map
  out.png = sprintf('%s/crw_sim_%s.png', out.dir, tagid); res = 72 # resolution in dots per inch
  dev.print(device=png, file=out.png, width=8*res, height=8*res, res=res)
  
  return(sim.alltags)
} #end function



createbufferabsence<-function(tags, tagid=tagid, n.sim = 100, buffdist = 100.000){
  
  print(tagid)
  tag = tags[which(tags$id==tagid),c('long', 'lat','dTime')]
  
  out.alltags.csv = sprintf('%s/bufferpts_sim_%s.csv', out.dir,tagid)
  
  # remove non-unique dates
  dupes = which(duplicated(tag$dTime))
  if (length(dupes)!=0){
    tag = tag[-dupes,]
  }
  
  # remove NAs
  tag = tag[!is.na(tag$lat),]
  tag.sp = tag
  
  # Creation of an object of class "ltraj"
  tr = as.ltraj(cbind(tag$lon, tag$lat), date=tag$dTime, id=tagid)
  tr1 = tr[[1]]
  
  head(tr1)
#  x11()
#  plot(tr)

  coordinates(tag.sp) <- ~long+lat
  dmat <- spDists(tag.sp)
  buffdist.deg<-buffdist * 0.0089982311916
  
  randomsamp <- function(tag, limit = buffdist.deg, n.samp = n.sim){

    coordinates(tag) <- ~long+lat
    bufferlims<-bbox(tag)
    bufferlims[1,1] <- bufferlims[1,1] - limit
    bufferlims[1,2] <- bufferlims[1,2] + limit
    bufferlims[2,1] <- bufferlims[2,1] - limit
    bufferlims[2,2] <- bufferlims[2,2] + limit
    
    return(spsample(tag, n.samp, type="random", bb=bufferlims))
    
  }
  
  #test code - THIS WORKS
  #randomsamp(tag)
  
  ### WHY ISN"T THIS WORKING???
  #buffer.list<-lapply(as.matrix(tag[,1:2]),randomsamp)
  
  for (i in 1:dim(tag)[1]){
    sim <- randomsamp(tag[i,])
    sim$dTime<-rep(tag$dTime[i],n.sim)
    sim$step<-rep(i,n.sim)
    sim$ID<-rep(tagid,n.sim)
    sim$iteration <- seq(1,n.sim)
    if (exists('sim.allbufftags')){
      sim.allbufftags = rbind(sim.allbufftags, sim)
    } else {
      sim.allbufftags = sim
    }
    
    write.csv(sim.allbufftags, out.alltags.csv, row.names=F)      #OLD (commented out 4ssm)
    
  }
  
  return(sim.allbufftags)
  
  #min.dist <- buffdist 
  #dmat[dmat <= min.dist] <- NA
  
}

createbackgroundabsence<-function(tags, tagid=tagid, n.sim = 100){
  
  
  print(tagid)
  tag = tags[which(tags$id==tagid),c('long', 'lat','dTime')]
  
  out.alltags.csv = sprintf('%s/backgroundpts_sim_%s.csv', out.dir,tagid)
  
  # remove non-unique dates
  dupes = which(duplicated(tag$dTime))
  if (length(dupes)!=0){
    tag = tag[-dupes,]
  }
  
  # remove NAs
  tag = tag[!is.na(tag$lat),]
  tag.sp = tag
  tags.sp = tags

  
  # Creation of an object of class "ltraj"
  tr = as.ltraj(cbind(tag$lon, tag$lat), date=tag$dTime, id=tagid)
  tr1 = tr[[1]]
  
  head(tr1)
  #  x11()
  #  plot(tr)
  
  coordinates(tag.sp) <- ~long+lat
  coordinates(tags.sp) <- ~long+lat
  
  
  studyarea=bbox(tags.sp)
  
  for (i in 1:dim(tag)[1]){
    tagi<-tag[i,]
    coordinates(tagi) <- ~long+lat
    sim <- spsample(tagi, n.sim, type="random", bb=studyarea)
    sim$dTime<-rep(tag$dTime[i],n.sim)
    sim$step<-rep(i,n.sim)
    sim$ID<-rep(tagid,n.sim)
    sim$iteration <- seq(1,n.sim)
    if (exists('sim.allbacktags')){
      sim.allbacktags = rbind(sim.allbacktags, sim)
    } else {
      sim.allbacktags = sim
    }
    
    write.csv(sim.allbacktags, out.alltags.csv, row.names=F)      #OLD (commented out 4ssm)
    
  }
  
  return(sim.allbacktags)
  
    
  
}

ResourceSelectionFunction<-function(tagid=tagid, n.sim = 100){
  
  #need code from Briana here?
  
}