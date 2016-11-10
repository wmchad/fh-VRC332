require(glmnet)
require(glmnetUtils)
require(dplyr)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetVariableSetData.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))
source(file.path(fnFolder, "SummarizePredictionResults.r"))

vlData <- GetTimepointData(fcData, 0:6, predSummary,
                           response="LogPeakVL")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-(1:2),"i"]-1]) %>%
    select(tp:ag, corPeak, shortVarName) %>% mutate(coef=coef1$x[-(1:2)])
##    tp         re                        ag    corPeak shortVarName         coef
## 1   4        MBL                       G73  0.3922100       var755  0.102162362
## 2   5      R2A.3                       V1a -0.2042609       var959 -0.065107218
## 3   4      R2A.3                     C1.Ak  0.2672538       var767  0.063294329
## 4   1      R3A.3            SIVsm.E660.2A5 -0.3631306       var341 -0.042306352
## 5   2 R2A.4.high           SIVsmH4.p55.Gag -0.3672093       var463 -0.030127008
## 6   0  R2A.4.low                     C1.Ak -0.2734684       var121 -0.029760898
## 7   5      R2A.3                       G49 -0.2731291       var943 -0.017510833
## 8   0      R2A.3             SIVsm.E660.84 -0.2992283        var97 -0.010690596
## 9   5        C1q J08.V1V2.E660.084.AVI.His -0.1897875       var908 -0.007747502
## 10  5        C1q                       G73 -0.2694557       var907 -0.001460580

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
dim(coef2)
## 712 variables
coefSummary2 <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef2[-(1:4),"i"]-1]) %>%
    select(tp:ag, corPeak, shortVarName) %>% mutate(coef=coef2$x[-(1:4)])

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 0.9062212

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)]) -
           vlData$y)^2))
## 0.8889812


pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictPeakVL/StraightGlmnet/Scaled")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

subData <- GetVariableSetData(fcData, coefSummary$shortVarName, TRUE, "LogPeakVL")
predResults <- RunRandomPartitionPredictions(
    subData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="PreChallenge-glmnet-CompRmse.txt",
    fitsFile="PreChallenge-glmnet-Fits.rdata",
    verbose=TRUE, progressEvery=10)


subData2 <- GetVariableSetData(fcData, coefSummary2$shortVarName, TRUE, "LogPeakVL")
predResults2 <- RunRandomPartitionPredictions(
    subData2, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="PreChallenge-glmnet-alpha-CompRmse.txt",
    fitsFile="PreChallenge-glmnet-alpha-Fits.rdata",
    verbose=TRUE, progressEvery=10)

resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
## 1.0295027 0.9813392 0.9187757 

## Model sizes range from 9 to 15 

## PickedCoeffs
##         var121         var341         var755         var767         var959 SIV_Mosaic_Env 
##             50             50             50             50             50             49 
##          var97         var907         var463         var908 
##             48             43             38             35 

resSummary2 <- SummarizePredictionResults(predResults2)
## Average RMSE: 
##  baseline     group predicted 
## 1.0107474 0.9739792 1.0303023 

## Model sizes range from 5 to 172 

## PickedCoeffs
## SIV_Mosaic_Env         var767         var755         var341         var897         var880 
##             43             43             42             35             34             29 
##         var907         var463         var959 
##             26             25             25 
