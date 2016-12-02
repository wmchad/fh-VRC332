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
## 1   5          R3A.3                    G119      var1014 -0.240197502
## 2   7          R3A.1         SIVsmH4.p55.Gag      var1354  0.153688069
## 3   0  aRhIgG.PE.low      SIVcpz.EK505.gp120        var32 -0.144013992
## 4   8            C1q          SIV.E543.gp140      var1428  0.126951015
## 5   7  aRhIgG.PE.low          SIV.E543.gp140      var1235  0.124963595
## 6   1     R2A.4.high J08.V1V2.mac239.AVI.His       var280  0.097555202
## 7   5          R3A.3                   C1.TR      var1013 -0.070783505
## 8   5          R3A.1      SIVcpz.EK505.gp120      var1003 -0.054112953
## 9   7 aRhIgG.PE.high       SIVmac251.BK.PR55      var1220 -0.044278760
## 10  8 aRhIgG.PE.high       SIVmac251.BK.PR55      var1392 -0.032838486
## 11  5          R3A.3                     G73      var1017 -0.021182020
## 12  3          R3A.3 SIVmac239.gp140.AVI.His       var683  0.016889695
## 13  1          R2A.3         SIVmac239.gp120       var263  0.013744691
## 14  5          R3A.1           SIVsm.E660.84      var1009 -0.013229371
## 15  2          R2A.2         SIVmac239.gp140       var422  0.005258732
## 16  6     R2A.4.high      SIVcpz.EK505.gp120      var1143  0.004377974
## 17  3          R3A.3      SIVcpz.EK505.gp120       var679  0.002951191


## Maybe how to do prediction?
vars <- coefSummary$shortVarName
subData <- GetVariableSetData(fcData, vars, TRUE, "NoChallenges")
coxModel <- coxph(vlData$surv~as.matrix(subData$x[,-(1:4)]), init=coefSummary$coef, iter=0)
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
## 3.847226
rmse(preds2, vlData$y)
## 2.991508

pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictNoChallenges/CoxPh")
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
    inTrain = createDataPartition(subData$groups, p=0.75, list=FALSE)
    train.x <- subData$x[inTrain,]
    train.y <- subData$y[inTrain]
    train.surv <- subData$surv[inTrain]
    test.x <- subData$x[-inTrain,]
    test.y <- subData$y[-inTrain]
    test.surv <- subData$surv[-inTrain]

    res <- cv.glmnet(x=as.matrix(train.x), y=train.surv, family="cox")
    coef <- as.data.frame(summary(predict(res, type="coef", s=glmnetRes$lambda.min)))
    vars <- (predSummary %>% filter(shortVarName %in% names(subData$x)[coef$i]))$shortVarName

    coxModel <- coxph(train.surv~as.matrix(train.x[,coef$i]), init=coef$x, iter=0)
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

    fit <- list(fit=res, predicted=preds, rmse=rmse(preds, test.y),
                predicted2=preds2, rmse2=rmse(preds2, test.y))

    fits[[i]] <- list(fit=fit$fit, predicted=fit$predicted,
                      predicted2=fit$predicted2, actual=test.y,
                      testAnimals=subData$animalIds[-inTrain])
    compRmse[i,] <- c(rmse(mean(train.y), test.y),
                      GroupRmse(train.x, train.y, test.x, test.y),
                      fit$rmse, fit$rmse2)
    if ( i %% progressEvery == 0 ) {
        vtimestamp()
        vprint(paste(i, "iterations finished"))
        ## setwd(outdir)
        ## write.table(compRmse, rmseFile,
        ##             quote=FALSE, row.names=FALSE, col.names=FALSE)
        ## save(fits, file=fitsFile)
    }
}

setwd(outdir)
write.table(compRmse, rmseFile,
            quote=FALSE, row.names=FALSE, col.names=FALSE)
save(fits, file=fitsFile)

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
