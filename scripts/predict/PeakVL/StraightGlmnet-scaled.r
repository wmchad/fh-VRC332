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

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-(1:2),"i"]-1]) %>%
    select(tp:ag, corPeak, shortVarName) %>% mutate(coef=coef1$x[-(1:2)])
##    tp             re                        ag    corPeak shortVarName          coef
## 1   8      R2A.4.low           SIVmac239.gp140  0.3730738      var1507  0.1719910968
## 2   5          R2A.3                       V1a -0.2042609       var959 -0.1037377559
## 3   4            MBL                       G73  0.3922100       var755  0.1016400163
## 4   8            C1q           SIVsmH4.p55.Gag -0.1975664      var1436 -0.0660874552
## 5   7            C1q                       G73  0.2101395      var1251  0.0598798340
## 6   0      R2A.4.low                       G73 -0.2738896       var126 -0.0433253957
## 7   4          R2A.3                     C1.Ak  0.2672538       var767  0.0408172021
## 8   1          R3A.3            SIVsm.E660.2A5 -0.3631306       var341 -0.0352135586
## 9   0          R2A.3             SIVsm.E660.84 -0.2992283        var97 -0.0327974906
## 10  5          R2A.3                       G49 -0.2731291       var943 -0.0237251981
## 11  4          R3A.1        SIVcpz.EK505.gp120 -0.3189028       var831 -0.0185997963
## 12  5            C1q J08.V1V2.E660.084.AVI.His -0.1897875       var908 -0.0177452902
## 13  8            C1q                       G49  0.1911714      var1422  0.0153716515
## 14  8      R2A.4.low            SIV.E543.gp140  0.1922405      var1506  0.0137558982
## 15  7 aRhIgG.PE.high           SIVsmH4.p55.Gag  0.1157801      var1223  0.0130093885
## 16  2     R2A.4.high            SIVsm.E660.2A5  0.2334368       var461  0.0122441658
## 17  0          R3A.3            SIVsm.E660.2A5 -0.2801948       var169 -0.0113801834
## 18  5  aRhIgG.PE.low         SIVmac251.BK.PR55  0.1665472       var897  0.0066147958
## 19  0          R3A.3            SIV.E543.gp140 -0.3133685       var162 -0.0037170293
## 20  7            C1q                  G145.146  0.1533543      var1249  0.0012534543
## 21  1          R2A.2 J08.V1V2.E660.2A5.AVI.His -0.2849406       var248 -0.0008899876

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
dim(coef2)
## 1 x 3

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 0.8027164

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)]) -
           vlData$y)^2))
## 1.017055


pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictPeakVL/StraightGlmnet/Scaled")
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
## 1.0114538 0.9578107 0.7495702 

## Model sizes range from 19 to 26 

## PickedCoeffs
##        var1223        var1251         var126        var1436        var1506         var755 
##             50             50             50             50             50             50 
##         var959          var97        var1507         var169         var341         var461 
##             50             50             49             49             49             49 
##         var767         var897 SIV_Mosaic_Env        var1422         var248        SIV_Gag 
##             49             49             48             48             48             47 
##         var162         var831         var908        var1249         var943        SIV_Env 
##             46             45             45             38             38             30 
##         x_PARI 
##             25 
