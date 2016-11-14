require(glmnet)
require(glmnetUtils)
require(dplyr)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetDeltaData.r"))
source(file.path(fnFolder, "GetVariableSetData.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))
source(file.path(fnFolder, "SummarizePredictionResults.r"))

deltaData <- GetDeltaData(fcData.na, predSummary)

vlData <- GetTimepointData(fcData, 1:8, predSummary,
                           response="LogPeakVL")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-(1:2),"i"]-1]) %>%
    select(tp:ag, corPeak, shortVarName) %>% mutate(coef=coef1$x[-(1:2)])
##    tp             re                        ag     corPeak shortVarName         coef
## 1   5 aRhIgG.PE.high                       V1a -0.18655099       var880 -0.160609859
## 2   8            C1q           SIVsmH4.p55.Gag -0.19756641      var1436 -0.159731221
## 3   7            C1q                       G73  0.21013952      var1251  0.149689975
## 4   8      R2A.4.low           SIVmac239.gp140  0.37307378      var1507  0.140831945
## 5   4          R2A.3                     C1.Ak  0.26725377       var767  0.116861625
## 6   4            MBL                       G73  0.39220997       var755  0.111395009
## 7   7            C1q         SIVmac251.BK.PR55 -0.05783829      var1261 -0.105362216
## 8   7 aRhIgG.PE.high           SIVsmH4.p55.Gag  0.11578011      var1223  0.083780858
## 9   8      R2A.4.low            SIV.E543.gp140  0.19224055      var1506  0.082492396
## 10  2     R2A.4.high           SIVsmH4.p55.Gag -0.36720929       var463 -0.074290147
## 11  2     R2A.4.high            SIVsm.E660.2A5  0.23343676       var461  0.071703126
## 12  8            C1q                       G49  0.19117143      var1422  0.069988347
## 13  8          R3A.1                      G119  0.17186295      var1511  0.069789891
## 14  1          R3A.3            SIVsm.E660.2A5 -0.36313061       var341 -0.068621967
## 15  4          R2A.3        SIVcpz.EK505.gp120  0.21461710       var778  0.067829931
## 16  5  aRhIgG.PE.low         SIVmac251.BK.PR55  0.16654715       var897  0.057731794
## 17  5            C1q J08.V1V2.E660.084.AVI.His -0.18978747       var908 -0.054472123
## 18  5          R2A.3                       V1a -0.20426090       var959 -0.052873473
## 19  1          R2A.2 J08.V1V2.E660.2A5.AVI.His -0.28494057       var248 -0.033700292
## 20  7          R3A.3                  G145.146  0.15032071      var1359  0.032904699
## 21  1          R3A.1                       V1a -0.17683157       var323 -0.030778706
## 22  5            MBL                     C1.TR -0.10967465       var923 -0.030343566
## 23  3          R3A.3        SIVcpz.EK505.gp120  0.23145542       var679  0.016795923
## 24  8            MBL                       G49  0.14589246      var1442  0.014158564
## 25  7            MBL                     C1.Ak -0.11688720      var1266 -0.013559093
## 26  7     R2A.4.high         SIVmac251.BK.PR55 -0.24623255      var1320 -0.011855091
## 27  8            C1q                     C1.TR  0.23169116      var1419  0.011204269
## 28  8          R3A.1                       G49  0.19528004      var1513  0.011138790
## 29  8          R3A.3           SIVsmH4.p55.Gag -0.17138748      var1547 -0.009033663
## 30  4          R3A.1        SIVcpz.EK505.gp120 -0.31890276       var831 -0.008899570
## 31  3          R2A.3        SIVcpz.EK505.gp120  0.26767900       var606  0.008694290
## 32  1          R2A.3           SIVmac239.gp120 -0.14379717       var263 -0.006218442
## 33  4            MBL   J08.V1V2.mac239.AVI.His  0.34018032       var758  0.005650447
## 34  3            MBL                       G49  0.10863155       var582  0.003525945

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
dim(coef2)
## 397 x 3

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 0.6044718

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)]) -
           vlData$y)^2))
## 0.975958


pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictPeakVL/StraightGlmnet/Delta")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

vars <- (predSummary %>%
         filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
         select(shortVarName))$shortVarName
subData <- GetVariableSetData(deltaData, vars, TRUE, "LogPeakVL")
predResults <- RunRandomPartitionPredictions(
    subData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-CompRmse.txt",
    fitsFile="glmnet-Fits.rdata",
    verbose=TRUE, progressEvery=10)

vars2 <- (predSummary %>%
          filter(shortVarName %in% names(vlData$x)[coef2[-1,"i"]-1]) %>%
          select(shortVarName))$shortVarName
subData2 <- GetVariableSetData(deltaData, vars2, TRUE, "LogPeakVL")
predResults2 <- RunRandomPartitionPredictions(
    subData2, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-alpha-CompRmse.txt",
    fitsFile="glmnet-alpha-Fits.rdata",
    verbose=TRUE, progressEvery=10)

resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
## 0.8861933 0.8929262 0.8768514 

## Model sizes range from 16 to 38 

## PickedCoeffs
## var1223  var679  var923  var461 var1251 var1507  var755  var463 var1422  var248  var959  var263 
##      50      50      49      48      47      47      47      45      42      42      41      40 
## SIV_Gag  var582 var1359  var897 var1513 var1547  var323  var341  var778 var1506  var606 var1320 
##      38      37      36      35      34      34      32      31      31      30      29      28 
## var1442 var1261 
##      28      26 

resSummary2 <- SummarizePredictionResults(predResults2)
## Average RMSE: 
##  baseline     group predicted 
## 0.8668670 0.8724611 0.8982137 

## Model sizes range from 6 to 119 

## PickedCoeffs
##  var943  var944  var931  var595 var1518 var1246  var291  var212  var752  var805  var926  var182 
##      39      38      37      35      33      32      31      29      29      29      27      25 
