### model eval functions
### written by HW 2019; edited by EH 2020; modified from code written by HW, SB, BA

library(caret)
library(mlbench)
library(gbm)
library(dismo)

### BRTS ####
## These functions call the gbm.fixed() function - differs from gbm.step in that is fits a predefined number of trees.
# dataInput: complete.cases input data with predictors. y term needs to be "presabs"
# gbm.x = character list of predictor names
# gbm. y = response term, must be "presabs"
# lr = learning rate
# tc = tree complexity
# bf = bag fraction
# nt = number of trees

kfolds_eval_brt <- function(dataInput, gbm.x, gbm.y, lr, tc, bf,nt){
  DataInput <- dataInput
  DataInput$Kset <- kfold(DataInput,10) #randomly allocate k groups
  Evaluations_kfold <- as.data.frame(matrix(data=0,nrow=10,ncol=4))
  colnames(Evaluations_kfold) <- c("k","AUC","TSS","TPR")
  counter=1
  for (k in 1:10){
    print(k)
    DataInput_train <- DataInput[DataInput$Kset!=k,]
    DataInput_test <- DataInput[DataInput$Kset==k,]
    DataInput.kfolds <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                                 family="bernoulli", tree.complexity=tc,
                                 learning.rate = lr, bag.fraction = bf,n.trees = nt)
    preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                         n.trees=nt, type="response")
    dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$presabs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    e <- evaluate(p=pres, a=abs)
    Evaluations_kfold[counter,1] <- k
    Evaluations_kfold[counter,2] <- e@auc
    Evaluations_kfold[counter,3] <- max(e@TPR + e@TNR-1)
    Evaluations_kfold[counter,4] <- mean(e@TPR)
    counter=counter+1 
  }
  return(Evaluations_kfold)
}

LOO_eval_brt <- function(DataInput, gbm.x, gbm.y, lr, tc, bf,nt){
  Evaluations_LOO <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_LOO) <- c("k","AUC","TSS","TPR")
  counter=1
  for (y in min(DataInput$Year):max(DataInput$Year)){
    print(y)
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    DataInput.loo <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                              family="bernoulli", tree.complexity=tc,
                              learning.rate = lr, bag.fraction = bf,n.trees = nt)
    preds <- predict.gbm(DataInput.loo, DataInput_test,
                         n.trees=nt, type="response")
    dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$presabs, preds)
    pres <- d[d[,1]==1,2]
    abs <- d[d[,1]==0,2]
    if(length(pres)>0 & length(abs)>0){
      e <- evaluate(p=pres, a=abs)
      
      Evaluations_LOO[counter,1] <- y
      Evaluations_LOO[counter,2] <- e@auc
      Evaluations_LOO[counter,3] <- max(e@TPR + e@TNR-1)
      Evaluations_LOO[counter,4] <- mean(e@TPR)
      counter=counter+1 
    }
  }
  return(Evaluations_LOO)}

eval_7525_brt <- function(dataInput, gbm.x, gbm.y, lr, tc, bf,nt){
  DataInput <- dataInput
  Evaluations_7525 <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_7525) <- c("AUC","TSS","TPR")
  DataInput_bound <- floor((nrow(DataInput)/4)*3)         #define % of training and test set
  DataInput_train<- DataInput[sample(nrow(DataInput),DataInput_bound),]
  DataInput_test<- DataInput[sample(nrow(DataInput),nrow(DataInput)-DataInput_bound),]
  DataInput.kfolds <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction = bf,n.trees = nt)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=nt, type="response")
  dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$presabs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_7525[1,1] <- e@auc
  Evaluations_7525[1,2] <- max(e@TPR + e@TNR-1)
  Evaluations_7535[1,3] <- mean(e@TPR)
  
  return(Evaluations_7525)}

eval_100_brt <- function(dataInput, gbm.x, gbm.y, lr, tc, bf,nt){
  DataInput <- dataInput
  Evaluations_100_percent <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_100_percent) <- c("AUC","TSS","TPR")
  DataInput_train<- DataInput
  DataInput_test<- DataInput
  DataInput.kfolds <- gbm.fixed(data=DataInput_train, gbm.x= gbm.x, gbm.y = gbm.y, 
                               family="bernoulli", tree.complexity=tc,
                               learning.rate = lr, bag.fraction =bf,n.trees = nt)
  preds <- predict.gbm(DataInput.kfolds, DataInput_test,
                       n.trees=nt, type="response")
  dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$presabs, preds)
  pres <- d[d[,1]==1,2]
  abs <- d[d[,1]==0,2]
  e <- evaluate(p=pres, a=abs)
  Evaluations_100_percent[1,1] <- e@auc
  Evaluations_100_percent[1,2] <- max(e@TPR + e@TNR-1)
  Evaluations_100_percent[1,2] <- mean(e@TPR)
  
  return(Evaluations_100_percent)}

dev_eval_brt=function(model_object){
  null <- model_object$self.statistics$null.deviance
  res <- model_object$self.statistics$resid.deviance
  dev=((null - res)/null)*100 
  return(dev)
}

ratio_brt <- function(dataInput,model_object,nt){
  ## Predict on model data using the best tree for predicting
  dataInput$RN=1:nrow(dataInput)
  BRTpred <- predict.gbm(model_object, dataInput, n.trees = nt, "response")
  # calculate ratio of observed to predicted values for study area
  ratio.BRTpred <- sum(dataInput$presabs)/sum(BRTpred)
  return(ratio.BRTpred)
}


### GAMM and GLMM ####

eval_100_gamm_glmm <- function(dataInput, formula) {
  DataInput=dataInput
  Evaluations_100_gamm <- as.data.frame(matrix(data=0,nrow=1,ncol=2))
  colnames(Evaluations_100_gamm) <- c("AUC","TSS")
  DataInput_train<- DataInput
  DataInput_test<- DataInput
  DataInput.100 <- mgcv::gamm(formula=formula,data=DataInput_train,random=list(tag=~1), family=binomial,niterPQL=50)
  preds<-as.data.frame(mgcv::predict.gam(DataInput.100$gam, DataInput_test, se=TRUE, type="response"))
  d <- cbind(DataInput_test$presabs, preds)
  pres <- as.numeric(d[d[,1]==1,2])
  abs <- as.numeric(d[d[,1]==0,2])
  e <- evaluate(p=pres, a=abs)
  Evaluations_100_gamm[1,1] <- e@auc
  Evaluations_100_gamm[1,2] <- max(e@TPR + e@TNR-1)
  return(Evaluations_100_gamm)
} 

eval_kfold_gamm_glmm <- function(dataInput, formula){
  DataInput=dataInput
  DataInput$Kset <- kfold(DataInput,10) #randomly allocate k groups
  Evaluations_kfold_gamm <- as.data.frame(matrix(data=0,nrow=10,ncol=3))
  colnames(Evaluations_kfold_gamm) <- c("k","AUC","TSS")
  counter=1
  for (k in 1:10){
    print(k)
    DataInput_train <- DataInput[DataInput$Kset!=k,]
    DataInput_test <- DataInput[DataInput$Kset==k,]
    DataInput.kfolds <- mgcv::gamm(formula=formula,data=DataInput_train,random=list(tag=~1), family=binomial,niterPQL=50)
    preds <- as.data.frame(mgcv::predict.gam(DataInput.kfolds$gam, DataInput_test,se=TRUE, type="response"))
    d <- cbind(DataInput_test$presabs, preds)
    pres <- as.numeric(d[d[,1]==1,2])
    abs <- as.numeric(d[d[,1]==0,2])
    e <- evaluate(p=pres, a=abs)
    Evaluations_kfold_gamm[counter,1] <- k
    Evaluations_kfold_gamm[counter,2] <- e@auc
    Evaluations_kfold_gamm[counter,3] <- max(e@TPR + e@TNR-1)
    counter=counter+1 
  }
  return(Evaluations_kfold_gamm)}

eval_LOO_gamm_glmm <- function(dataInput, formula){
  DataInput=dataInput
  Evaluations_LOO_gamm <- as.data.frame(matrix(data=0,nrow=1,ncol=3))
  colnames(Evaluations_LOO_gamm) <- c("k","AUC","TSS")
  counter=1
  for (y in unique(DataInput$Year)){
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    print(y)
    DataInput.loo <- mgcv::gamm(formula=formula,data=DataInput_train,random=list(tag=~1), family=binomial,niterPQL=50)
    preds <- as.data.frame(mgcv::predict.gam(DataInput.loo$gam, DataInput_test,se=TRUE, type="response"))
    d <- cbind(DataInput_test$presabs, preds)
    pres <- as.numeric(d[d[,1]==1,2])
    abs <- as.numeric(d[d[,1]==0,2])
    e <- evaluate(p=pres, a=abs)
    Evaluations_LOO_gamm[counter,1] <- y
    Evaluations_LOO_gamm[counter,2] <- e@auc
    Evaluations_LOO_gamm[counter,3] <- max(e@TPR + e@TNR-1)
    counter=counter+1 
  }
  return(Evaluations_LOO_gamm)}

eval_7525_gamm_glmm <- function(dataInput, formula){
  Evaluations_7525 <- as.data.frame(matrix(data=0,nrow=1,ncol=2))
  colnames(Evaluations_7525) <- c("AUC","TSS")
  dataInput_bound <- floor((nrow(dataInput)/4)*3)         #define % of training and test set
  dataInput_train<- dataInput[sample(nrow(dataInput),dataInput_bound),]
  dataInput_test<- dataInput[sample(nrow(dataInput),nrow(dataInput)-dataInput_bound),]
  dataInput.7525 <- mgcv::gamm(formula=formula,data=dataInput_train,random=list(tag=~1), family=binomial,niterPQL=50)
  preds<-as.data.frame(mgcv::predict.gam(dataInput.7525$gam, dataInput_test, se=TRUE, type="response"))
  dev <- calc.deviance(obs=dataInput_test$presabs, pred=preds$fit, calc.mean=TRUE)
  d <- cbind(dataInput_test$presabs, preds$fit)
  pres <- as.numeric(d[d[,1]==1,2])
  abs <- as.numeric(d[d[,1]==0,2])
  e <- evaluate(p=pres, a=abs)
  Evaluations_7525[1,1] <- e@auc
  Evaluations_7525[1,2] <- max(e@TPR + e@TNR-1)
  
  return(Evaluations_7525)}

rsq_gamm_glmm=function(model_object){ ### still lost
   rsq=summary(glm$gam)$r.sq
  return(rsq)
}

ratio_gamm_glmm <- function(dataInput,model_object){
  ## Predict on model data using the best tree for predicting
  dataInput$RN=1:nrow(dataInput)
  GMpred <- as.data.frame(mgcv::predict.gam(dataInput.7525$gam, dataInput_test, se=TRUE, type="response"))
  # calculate ratio of observed to predicted values for study area
  ratio.GMpred <- sum(dataInput$presabs)/sum(GMpred$fit)
  return(ratio.GMpred)
}

pseudoR2.BRT <- function(x){
  d2 <- 1-(x$self.statistics$mean.resid/x$self.statistics$mean.null)
  return(d2)
}

eval_LOO_BRT <- function(dataInput, gbm.x, gbm.y, lr=0.05, months=c(1:12)){
  DataInput <- dataInput[dataInput$month %in% months,]
  DataInput$Year <- lubridate::year(DataInput$dt)
  Evaluations_LOO_BRT <- as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(Evaluations_LOO_BRT) <- c("k","Deviance","AUC","TSS","TPR")
  counter=1
  for (y in unique(DataInput$Year)){
    if(any(months %in% DataInput[DataInput$Year==y & DataInput$presence==1,]$month)==FALSE) next #skip year if no months in dataset are in months vector
    print(y)
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    DataInput.loo <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = c("presence"), 
                              family="bernoulli", tree.complexity=3,
                              learning.rate = lr, bag.fraction = 0.6)
    preds <- predict.gbm(DataInput.loo, DataInput_test,
                         n.trees=DataInput.loo$gbm.call$best.trees, type="response")
    dev <- calc.deviance(obs=DataInput_test$presence, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$presence, preds)
    pres <- as.numeric(d[d[,1]==1,2])
    abs <- as.numeric(d[d[,1]==0,2])
    e <- dismo::evaluate(p=pres, a=abs)
    
    Evaluations_LOO_BRT[counter,1] <- y
    Evaluations_LOO_BRT[counter,2] <- dev
    Evaluations_LOO_BRT[counter,3] <- e@auc
    Evaluations_LOO_BRT[counter,4] <- max(e@TPR + e@TNR-1)
    Evaluations_LOO_BRT[counter,5] <- mean(e@TPR)
    counter=counter+1 
  }
  return(Evaluations_LOO_BRT)}

#demo BRT: eval_LOO_BRT(DataInput_Fit, gbm.x=c("log_epi_mnk_pb","log_epi_mnk_pp", "log_meso_mnk_pp", "log_mmeso_mnk_pp","log_pk_pb", "log_pk_pp", "temperature", "sqrt_distcoast", "log_EKE","sqrt_SLPd","sqrt_dFronts","sqrt_dEddies", "Dpth"), "presence")

eval_LOO_GAMM <- function(dataInput, formula, months=c(1:12)){
  DataInput <- dataInput[dataInput$month %in% months,]
  DataInput$ptt <- as.numeric(DataInput$ptt)
  DataInput$Year <- lubridate::year(DataInput$dt)
  Evaluations_LOO_GAMM <- as.data.frame(matrix(data=0,nrow=1,ncol=4))
  colnames(Evaluations_LOO_GAMM) <- c("k","AUC","TSS","TPR")
  counter=1
  for (y in unique(DataInput$Year)){
    #skip year if no months in dataset are in months vector
    if(any(months %in% DataInput[DataInput$Year==y & DataInput$presence==1,]$month)==FALSE |
       any(months %in% DataInput[DataInput$Year==y & DataInput$presence==0,]$month)==FALSE) 
      next 
    DataInput_train <- DataInput[DataInput$Year!=y,]
    DataInput_test <- DataInput[DataInput$Year==y,]
    print(y)
    DataInput.loo <- mgcv::gam(formula=formula,data=DataInput_train, family=binomial,method="REML")
    preds <- as.data.frame(predict.gam(DataInput.loo, DataInput_test,se=TRUE, type="response"))
    d <- cbind(DataInput_test$presence, preds)
    pres <- as.numeric(d[d[,1]==1,2])
    abs <- as.numeric(d[d[,1]==0,2])
    e <- evaluate(p=pres, a=abs)
    Evaluations_LOO_GAMM[counter,1] <- y
    Evaluations_LOO_GAMM[counter,2] <- e@auc
    Evaluations_LOO_GAMM[counter,3] <- max(e@TPR + e@TNR-1)
    Evaluations_LOO_GAMM[counter,4] <- mean(e@TPR)
    counter=counter+1 
  }
  return(Evaluations_LOO_GAMM)}

eval_LOOm_BRT <- function(dataInput, gbm.x, gbm.y, lr=0.05, months=c(1:12)){
  DataInput <- dataInput[dataInput$Month %in% months,]
  #DataInput$Year <- lubridate::year(DataInput$dt)
  Evaluations_LOO_BRT <- as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(Evaluations_LOO_BRT) <- c("k","Deviance","AUC","TSS","TPR")
  counter=1
  for (m in unique(months)){
    if(any(unique(DataInput$Year) %in% DataInput[DataInput$Month==m & DataInput$presabs==1,]$Year)==FALSE) next #skip year if no months in dataset are in months vector
    print(paste0("Month: ",m))
    DataInput_train <- DataInput[DataInput$Month!=m,]
    DataInput_test <- DataInput[DataInput$Month==m,]
    DataInput.loo <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = c("presabs"), 
                              family="bernoulli", tree.complexity=3,
                              learning.rate = lr, bag.fraction = 0.6)
    preds <- predict.gbm(DataInput.loo, DataInput_test,
                         n.trees=DataInput.loo$gbm.call$best.trees, type="response")
    dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
    d <- cbind(DataInput_test$presabs, preds)
    pres <- as.numeric(d[d[,1]==1,2])
    abs <- as.numeric(d[d[,1]==0,2])
    e <- dismo::evaluate(p=pres, a=abs)
    
    Evaluations_LOO_BRT[counter,1] <- m
    Evaluations_LOO_BRT[counter,2] <- dev
    Evaluations_LOO_BRT[counter,3] <- e@auc
    Evaluations_LOO_BRT[counter,4] <- max(e@TPR + e@TNR-1)
    Evaluations_LOO_BRT[counter,5] <- mean(e@TPR)
    counter=counter+1 
  }
  return(Evaluations_LOO_BRT)}

#demo BRT: eval_LOO_BRT(DataInput_Fit, gbm.x=c("log_epi_mnk_pb","log_epi_mnk_pp", "log_meso_mnk_pp", "log_mmeso_mnk_pp","log_pk_pb", "log_pk_pp", "temperature", "sqrt_distcoast", "log_EKE","sqrt_SLPd","sqrt_dFronts","sqrt_dEddies", "Dpth"), "presence")

eval_LOOm_GAMM <- function(dataInput, formula, months=c(1:12)){
  DataInput <- dataInput[dataInput$Month %in% months,]
  years=unique(DataInput$Year)
  DataInput$tag <- as.numeric(DataInput$tag)
  #DataInput$Year <- lubridate::year(DataInput$dt)
  Evaluations_LOO_GAMM <- as.data.frame(matrix(data=0,nrow=1,ncol=4))
  colnames(Evaluations_LOO_GAMM) <- c("k","AUC","TSS","TPR")
  counter=1
  for (m in unique(DataInput$Month)){
    #skip year if no months in dataset are in months vector
    if(any(years %in% DataInput[DataInput$Month==m & DataInput$presabs==1,]$Year)==FALSE |
       any(years %in% DataInput[DataInput$Month==m & DataInput$presabs==0,]$Year)==FALSE) 
      next 
    DataInput_train <- DataInput[DataInput$Month!=m,]
    DataInput_test <- DataInput[DataInput$Month==m,]
    print(paste0("Month: ",m))
    DataInput.loo <- mgcv::gam(formula=formula,data=DataInput_train, family=binomial,method="REML")
    preds <- as.data.frame(predict.gam(DataInput.loo, DataInput_test,se=TRUE, type="response"))
    d <- cbind(DataInput_test$presabs, preds)
    pres <- as.numeric(d[d[,1]==1,2])
    abs <- as.numeric(d[d[,1]==0,2])
    e <- evaluate(p=pres, a=abs)
    Evaluations_LOO_GAMM[counter,1] <- m
    Evaluations_LOO_GAMM[counter,2] <- e@auc
    Evaluations_LOO_GAMM[counter,3] <- max(e@TPR + e@TNR-1)
    Evaluations_LOO_GAMM[counter,4] <- mean(e@TPR)
    counter=counter+1 
  }
  return(Evaluations_LOO_GAMM)}

eval_LOOspace_BRT <- function(dataInput, gbm.x, gbm.y, lr=0.05, trainextent=c(20,60,-140,-110), testextent=c(20,60,-140,-110), traindata=NA, testdata=NA){
  
  DataInput<-dataInput
  DataInput_train <- traindata
  DataInput_test <- testdata
  
  if (is.null(dim(traindata))){
    DataInput_train <- dataInput[dataInput$lat>trainextent[1]&dataInput$lat<trainextent[2]&dataInput$lon>trainextent[3]&dataInput$lon<trainextent[4]]
    DataInput_test <- dataInput[dataInput$lat>testextent[1]&dataInput$lat<testextent[2]&dataInput$lon>testextent[3]&dataInput$lon<testextent[4]]
  }
  
  #DataInput$Year <- lubridate::year(DataInput$dt)
  Evaluations_LOO_BRT <- as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(Evaluations_LOO_BRT) <- c("k","Deviance","AUC","TSS","TPR")
  DataInput.loo <- gbm.step(data=DataInput_train, gbm.x= gbm.x, gbm.y = c("presabs"), 
                            family="bernoulli", tree.complexity=3,
                            learning.rate = lr, bag.fraction = 0.6)
  preds <- predict.gbm(DataInput.loo, DataInput_test,
                       n.trees=DataInput.loo$gbm.call$best.trees, type="response")
  dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds, calc.mean=TRUE)
  d <- cbind(DataInput_test$presabs, preds)
  pres <- as.numeric(d[d[,1]==1,2])
  abs <- as.numeric(d[d[,1]==0,2])
  e <- dismo::evaluate(p=pres, a=abs)
  
  Evaluations_LOO_BRT[1,1] <- dim(DataInput_train)[1]
  Evaluations_LOO_BRT[1,2] <- dev
  Evaluations_LOO_BRT[1,3] <- e@auc
  Evaluations_LOO_BRT[1,4] <- max(e@TPR + e@TNR-1)
  Evaluations_LOO_BRT[1,4] <- mean(e@TPR)
  
  return(Evaluations_LOO_BRT)}

#demo BRT: eval_LOO_BRT(DataInput_Fit, gbm.x=c("log_epi_mnk_pb","log_epi_mnk_pp", "log_meso_mnk_pp", "log_mmeso_mnk_pp","log_pk_pb", "log_pk_pp", "temperature", "sqrt_distcoast", "log_EKE","sqrt_SLPd","sqrt_dFronts","sqrt_dEddies", "Dpth"), "presence")

eval_LOOspace_GAMM <- function(dataInput, formula, trainextent=c(20,60,-140,-110), testextent=c(20,60,-140,-110), traindata=NA, testdata=NA){
  
  DataInput<-dataInput
  DataInput_train <- traindata
  DataInput_test <- testdata
  
  if (is.null(dim(traindata))){
    DataInput_train <- dataInput[dataInput$lat>trainextent[1]&dataInput$lat<trainextent[2]&dataInput$lon>trainextent[3]&dataInput$lon<trainextent[4]]
    DataInput_test <- dataInput[dataInput$lat>testextent[1]&dataInput$lat<testextent[2]&dataInput$lon>testextent[3]&dataInput$lon<testextent[4]]
  }
  
  DataInput$tag <- as.numeric(DataInput$tag)
  #DataInput$Year <- lubridate::year(DataInput$dt)
  Evaluations_LOO_GAMM <- as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(Evaluations_LOO_GAMM) <- c("k","Deviance","AUC","TSS","TPR")
  
  DataInput.loo <- mgcv::gam(formula=formula,data=DataInput_train, family=binomial,method="REML")
  preds <- as.data.frame(predict.gam(DataInput.loo, DataInput_test,se=TRUE, type="response"))
  d <- cbind(DataInput_test$presabs, preds)
  dev <- calc.deviance(obs=DataInput_test$presabs, pred=preds$fit, calc.mean=TRUE)
  pres <- as.numeric(d[d[,1]==1,2])
  abs <- as.numeric(d[d[,1]==0,2])
  e <- evaluate(p=pres, a=abs)
  Evaluations_LOO_GAMM[1,1] <- dim(DataInput_train)[1]
  Evaluations_LOO_GAMM[1,2] <- dev
  Evaluations_LOO_GAMM[1,3] <- e@auc
  Evaluations_LOO_GAMM[1,4] <- max(e@TPR + e@TNR-1)
  Evaluations_LOO_GAMM[1,5] <- mean(e@TPR)
  
  return(Evaluations_LOO_GAMM)
}
