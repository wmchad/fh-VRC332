require(glmnet)
require(glmnetUtils)
require(dplyr)
require(survival)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetDeltaData.r"))
source(file.path(fnFolder, "GetVariableSetData.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))
source(file.path(fnFolder, "SummarizePredictionResults.r"))

vlData <- GetTimepointData(fcData, 0:8, predSummary,
                           response="NoChallenges")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

event <- rep(1, length(vlData$y))
event[vlData$y==13] <- 0

vlData$surv <- Surv(vlData$y, event)

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$surv, family="cox")
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$surv, family="cox")

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1$i]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=as.numeric(coef1$x))
coefSummary %>% arrange(desc(abs(coef)))
##   tp             re                      ag shortVarName         coef
## 1   5          R3A.3                    G119      var1014 -0.258161862
## 2   0  aRhIgG.PE.low      SIVcpz.EK505.gp120        var32 -0.168982380
## 3   7  aRhIgG.PE.low          SIV.E543.gp140      var1235  0.149624970
## 4   8            C1q          SIV.E543.gp140      var1428  0.148341996
## 5   7          R3A.1         SIVsmH4.p55.Gag      var1354  0.143970423
## 6   1     R2A.4.high J08.V1V2.mac239.AVI.His       var280  0.113516567
## 7   5          R3A.3                   C1.TR      var1013 -0.089949544
## 8   7 aRhIgG.PE.high       SIVmac251.BK.PR55      var1220 -0.055079717
## 9   5          R3A.1      SIVcpz.EK505.gp120      var1003 -0.047728104
## 10  8 aRhIgG.PE.high       SIVmac251.BK.PR55      var1392 -0.046834865
## 11  3          R3A.3 SIVmac239.gp140.AVI.His       var683  0.039980181
## 12  5          R3A.1           SIVsm.E660.84      var1009 -0.035324286
## 13  1          R2A.3         SIVmac239.gp120       var263  0.028843243
## 14  5          R3A.3                     G73      var1017 -0.026091262
## 15  3          R3A.3      SIVcpz.EK505.gp120       var679  0.012822499
## 16  2          R2A.2         SIVmac239.gp140       var422  0.006472271
## 17  6     R2A.4.high      SIVcpz.EK505.gp120      var1143  0.004822233
## 18  1          R3A.1         SIVsmH4.p55.Gag       var322 -0.001865904



coxIter <- 100

## Maybe how to do prediction?
vars <- coefSummary$shortVarName
subData <- GetVariableSetData(fcData, vars, TRUE, "NoChallenges")
coxModel <- coxph(vlData$surv~as.matrix(subData$x[,-(1:4)]), init=coefSummary$coef, iter=coxIter)
coxFit <- survfit(coxModel, newData=as.matrix(subData$x[,-(1:4)]))
s0 <- exp(-coxFit$cumhaz)
preds <- sapply(1:nrow(subData$x), function(i) {
    sum(s0^exp(sum(subData$x[i,-(1:4)] * coefSummary$coef)))
})
preds2 <- sapply(1:nrow(subData$x), function(i) {
    sum(s0^predict(glmnetRes,
                   newx=as.matrix(vlData$x[i,]),
                   s=glmnetRes$lambda.min,
                   type="response"))
})
rmse(preds, vlData$y)
## 3.762668
rmse(preds2, vlData$y)
## 3.310601

pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Projects/VRC332/Data/PredictNoChallenges/CoxPh/Original/All")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
rmseFile <- "coxph-CompRmse.txt"
fitsFile <- "coxph-Fits.rdata"

vars <- (predSummary %>%
         filter(shortVarName %in% names(vlData$x)[coef1$i]) %>%
         select(shortVarName))$shortVarName
subData <- GetVariableSetData(fcData, vars, TRUE, "NoChallenges")
event <- rep(1, length(subData$y))
event[subData$y==13] <- 0
subData$surv <- Surv(subData$y, event)

vtimestamp <- function(){timestamp()}
vprint <- function(..., sep=""){print(paste(..., sep=sep))}
progressEvery <- 5

curdir <- getwd()
compRmse <- matrix(nrow=nRand, ncol=4)
colnames(compRmse) <- c("baseline", "group", "predicted", "predicted2")
fits <- NULL


vtimestamp()
vprint("Starting...")
for ( i in 1:nRand ) {
    inTrain = createDataPartition(vlData$groups, p=0.75, list=FALSE)
    train.x <- vlData$x[inTrain,]
    train.y <- vlData$y[inTrain]
    train.surv <- vlData$surv[inTrain]
    test.x <- vlData$x[-inTrain,]
    test.y <- vlData$y[-inTrain]
    test.surv <- vlData$surv[-inTrain]

    res <- cv.glmnet(x=as.matrix(train.x), y=train.surv, family="cox")
    coef <- as.data.frame(summary(predict(res, type="coef", s=res$lambda.min)))
    vars <- (predSummary %>% filter(shortVarName %in% names(vlData$x)[coef$i]))$shortVarName

    coxModel <- coxph(train.surv~as.matrix(train.x[,coef$i]), init=coef$x, iter=coxIter)
    coxFit <- survfit(coxModel, newData=as.matrix(train.x[,coef$i]))
    s0 <- exp(-coxFit$cumhaz)
    preds <- sapply(1:nrow(test.x), function(i) {
        sum(s0^exp(sum(test.x[i,coef$i] * coef$x)))
    })
    preds2 <- sapply(1:nrow(test.x), function(i) {
        sum(s0^predict(res,
                       newx=as.matrix(test.x[i,]),
                       s=res$lambda.min,
                       type="response"))
    })

    fit <- list(fit=res, vars=vars,
                predicted=preds, rmse=rmse(preds, test.y),
                predicted2=preds2, rmse2=rmse(preds2, test.y))

    fits[[i]] <- list(fit=fit$fit, vars=vars, predicted=fit$predicted,
                      predicted2=fit$predicted2, actual=test.y,
                      testAnimals=vldata$animalIds[-inTrain])
    compRmse[i,] <- c(rmse(mean(train.y), test.y),
                      GroupRmse(train.x, train.y, test.x, test.y),
                      fit$rmse, fit$rmse2)
    if ( i %% progressEvery == 0 ) {
        vtimestamp()
        vprint(paste(i, "iterations finished"))
        setwd(outdir)
        write.table(compRmse, rmseFile,
                    quote=FALSE, row.names=FALSE, col.names=FALSE)
        save(fits, file=fitsFile)
    }
}



setwd(outdir)
write.table(compRmse, rmseFile,
            quote=FALSE, row.names=FALSE, col.names=FALSE)
save(fits, file=fitsFile)



allVars <- NULL
allVars2 <- NULL
for ( i in 1:nRand ) {
    ## res <- fits[[i]]$fit
    ## coef <- as.data.frame(summary(predict(res, type="coef", s=res$lambda.min)))
    ## vars <- (predSummary %>% filter(shortVarName %in% names(vlData$x)[coef$i]))$shortVarName

    coef <- as.data.frame(summary(predict(fits[[i]]$fit$glmnet.fit,
                                          type="coef",
                                          s=fits[[i]]$fit$glmnet.fit$lambda.min)))
    vars <- (predSummary %>% filter(shortVarName %in% names(vlData$x)[coef$i]))$shortVarName
    allVars <- c(allVars, vars)
    allVars2 <- c(allVars2, fits[[i]]$vars)
}

GetVars <- function(vars) {
    vtbl <- table(vars)
    names(vtbl)[vtbl >= nRand/2]
}

varset1 <- GetVars(allVars)
varset2 <- GetVars(allVars2)

subData1 <- GetVariableSetData(fcData, varset1, TRUE, "NoChallenges")
event <- rep(1, length(subData1$y))
event[subData1$y==13] <- 0
subData1$surv <- Surv(subData1$y, event)
subData2 <- GetVariableSetData(fcData, varset1, TRUE, "NoChallenges")
event <- rep(1, length(subData2$y))
event[subData2$y==13] <- 0
subData2$surv <- Surv(subData2$y, event)



compRmse <- matrix(nrow=nRand, ncol=6)
colnames(compRmse) <- c("baseline", "group",
                        "predicted1a", "predicted1b",
                        "predicted2a", "predicted2b")
fits1 <- NULL
fits2 <- NULL

rmseFile <- "coxph-BestCompRmse.txt"
fits1File <- "coxph-BestFits1.rdata"
fits2File <- "coxph-BestFits2.rdata"


vtimestamp()
vprint("Starting...")
for ( i in 1:nRand ) {
    inTrain = createDataPartition(subData1$groups, p=0.75, list=FALSE)
    train1.x <- subData1$x[inTrain,]
    train1.y <- subData1$y[inTrain]
    train1.surv <- subData1$surv[inTrain]
    test1.x <- subData1$x[-inTrain,]
    test1.y <- subData1$y[-inTrain]
    test1.surv <- subData1$surv[-inTrain]
    train2.x <- subData2$x[inTrain,]
    train2.y <- subData2$y[inTrain]
    train2.surv <- subData2$surv[inTrain]
    test2.x <- subData2$x[-inTrain,]
    test2.y <- subData2$y[-inTrain]
    test2.surv <- subData2$surv[-inTrain]

    res1 <- cv.glmnet(x=as.matrix(train1.x), y=train1.surv, family="cox")
    coef1 <- as.data.frame(summary(predict(res1, type="coef", s=res1$lambda.min)))
    vars1 <- (predSummary %>% filter(shortVarName %in% names(subData1$x)[coef1$i]))$shortVarName

    res2 <- cv.glmnet(x=as.matrix(train2.x), y=train2.surv, family="cox")
    coef2 <- as.data.frame(summary(predict(res2, type="coef", s=res2$lambda.min)))
    vars2 <- (predSummary %>% filter(shortVarName %in% names(subData2$x)[coef2$i]))$shortVarName

    coxModel1 <- coxph(train1.surv~as.matrix(train1.x[,coef1$i]), init=coef1$x, iter=coxIter)
    coxFit1 <- survfit(coxModel1, newData=as.matrix(train1.x[,coef1$i]))
    s01 <- exp(-coxFit1$cumhaz)
    preds1a <- sapply(1:nrow(test1.x), function(i) {
        sum(s01^exp(sum(test1.x[i,coef1$i] * coef1$x)))
    })
    preds1b <- sapply(1:nrow(test1.x), function(i) {
        sum(s01^predict(res1,
                       newx=as.matrix(test1.x[i,]),
                       s=res1$lambda.min,
                       type="response"))
    })
    
    coxModel2 <- coxph(train2.surv~as.matrix(train2.x[,coef2$i]), init=coef2$x, iter=coxIter)
    coxFit2 <- survfit(coxModel2, newData=as.matrix(train2.x[,coef2$i]))
    s02 <- exp(-coxFit2$cumhaz)
    preds2a <- sapply(1:nrow(test2.x), function(i) {
        sum(s02^exp(sum(test2.x[i,coef2$i] * coef2$x)))
    })
    preds2b <- sapply(1:nrow(test2.x), function(i) {
        sum(s02^predict(res2,
                       newx=as.matrix(test2.x[i,]),
                       s=res2$lambda.min,
                       type="response"))
    })

    fit <- list(fit=res, vars=vars,
                predicted=preds, rmse=rmse(preds, test.y),
                predicted2=preds2, rmse2=rmse(preds2, test.y))

    fits1[[i]] <- list(fit=res1, vars=vars1, predicted=preds1a,
                      predicted2=preds1b, actual=test1.y,
                      testAnimals=subData$animalIds[-inTrain])
    fits2[[i]] <- list(fit=res2, vars=vars2, predicted=preds2a,
                      predicted2=preds2b, actual=test2.y,
                      testAnimals=subData$animalIds[-inTrain])
    compRmse[i,] <- c(rmse(mean(train.y), test.y),
                      GroupRmse(train.x, train.y, test.x, test.y),
                      rmse(preds1a, test1.y), rmse(preds1b, test1.y),
                      rmse(preds2a, test2.y), rmse(preds2b, test2.y))
    if ( i %% progressEvery == 0 ) {
        vtimestamp()
        vprint(paste(i, "iterations finished"))
        setwd(outdir)
        write.table(compRmse, rmseFile,
                    quote=FALSE, row.names=FALSE, col.names=FALSE)
        save(fits1, file=fits1File)
        save(fits1, file=fits2File)
    }
}



predResults <- list(compRmse=compRmse,
                    fits=fits)


predResults <- RunRandomPartitionPredictions(
    subData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-CompRmse.txt",
    fitsFile="glmnet-Fits.rdata",
    verbose=TRUE, progressEvery=10)

vars2 <- (predSummary %>%
          filter(shortVarName %in% names(vlData$x)[coef2[-1,"i"]-1]) %>%
          select(shortVarName))$shortVarName
subData2 <- GetVariableSetData(deltaData, vars2, TRUE, "NoChallenges")
predResults2 <- RunRandomPartitionPredictions(
    subData2, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-alpha-CompRmse.txt",
    fitsFile="glmnet-alpha-Fits.rdata",
    verbose=TRUE, progressEvery=10)

resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
##  3.699131  3.579199  2.244707 

## Model sizes range from 15 to 20 

## PickedCoeffs
##        SIV_Gag        var1233        var1237        var1245         var204         var434 
##             50             50             50             50             50             50 
##         var679         var684         var690         var732         var863        var1239 
##             50             50             50             50             50             49 
##         var263         var505        var1270        var1242 SIV_Mosaic_Env        SIV_Env 
##             49             48             45             44             42             33 

resSummary2 <- SummarizePredictionResults(predResults2)
## Average RMSE: 
##  baseline     group predicted 
##  3.749298  3.692931  3.052328 

## Model sizes range from 5 to 12 

## PickedCoeffs
##         var505        var1239        var1237        var1231        var1242        var1245 
##             50             49             48             46             46             40 
##         x_PARI SIV_Mosaic_Env 
##             36             26 
