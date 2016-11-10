require(glmnet)
require(glmnetUtils)
require(dplyr)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetVariableSetData.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))
source(file.path(fnFolder, "SummarizePredictionResults.r"))

vlData <- GetTimepointData(fcData, 0:8, predSummary,
                           response="LogPeakVL")

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1])
##   tp             re                        ag numObs numObs2    corPeak  shortVarName
## 1   0          R2A.3             SIVsm.E660.84     69      69 -0.2992283        var97
## 2   0      R2A.4.low                       G73     69      69 -0.2738896       var126
## 3   0          R3A.3            SIV.E543.gp140     69      69 -0.3133685       var162
## 4   0          R3A.3            SIVsm.E660.2A5     69      69 -0.2801948       var169
## 5   1          R2A.2 J08.V1V2.E660.2A5.AVI.His     70      70 -0.2849406       var248
## 6   1          R3A.3            SIVsm.E660.2A5     70      70 -0.3631306       var341
## 7   2     R2A.4.high            SIVsm.E660.2A5     70      70  0.2334368       var461
## 8   4            MBL                       G73     70      70  0.3922100       var755
## 9   4          R2A.3                     C1.Ak     70      70  0.2672538       var767
## 10  4          R3A.1        SIVcpz.EK505.gp120     70      70 -0.3189028       var831
## 11  5  aRhIgG.PE.low         SIVmac251.BK.PR55     87      87  0.1665472       var897
## 12  5            C1q J08.V1V2.E660.084.AVI.His     87      87 -0.1897875       var908
## 13  5          R2A.3                       G49     87      87 -0.2731291       var943
## 14  5          R2A.3                       V1a     87      87 -0.2042609       var959
## 15  7 aRhIgG.PE.high           SIVsmH4.p55.Gag     87      87  0.1157801      var1223
## 16  7            C1q                  G145.146     87      87  0.1533543      var1249
## 17  7            C1q                       G73     87      87  0.2101395      var1251
## 18  8            C1q                       G49     87      87  0.1911714      var1422
## 19  8            C1q           SIVsmH4.p55.Gag     87      87 -0.1975664      var1436
## 20  8      R2A.4.low            SIV.E543.gp140     87      87  0.1922405      var1506
## 21  8      R2A.4.low           SIVmac239.gp140     87      87  0.3730738      var1507

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=0)))
dim(coef2)
## 1553 x 3

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 0.7844463

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=0) -
           vlData$y)^2))
## 0.9256206


pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictPeakVL/StraightGlmnet")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

vars <- (predSummary %>%
         filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
         select(shortVarName))$shortVarName
subData <- GetVariableSetData(fcData, vars, TRUE, "LogPeakVL")
predResults <- RunRandomPartitionPredictions(
    subData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-CompRmse.txt",
    fitsFile="glmnet-Fits.rdata",
    verbose=TRUE, progressEvery=10)

resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
## 0.9914205 0.9446177 0.7191968 

## Model sizes range from 20 to 26 

## PickedCoeffs
##        var1223        var1251         var126        var1436        var1506        var1507 
##             50             50             50             50             50             50 
##         var169         var248         var341         var461         var755         var959 
##             50             50             50             50             50             50 
##        var1422         var162         var767         var897          var97         var831 
##             49             49             49             49             49             48 
## SIV_Mosaic_Env        SIV_Gag         var908        var1249        SIV_Env         var943 
##             47             45             43             40             35             32 
##         x_PARI 
##             31 
