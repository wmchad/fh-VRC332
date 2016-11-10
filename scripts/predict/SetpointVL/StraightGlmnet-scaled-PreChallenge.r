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
                           response="LogSetpointVL")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
    select(tp:ag, corSetpoint, shortVarName) %>% mutate(coef=coef1$x[-1])
##    tp             re                        ag corSetpoint shortVarName         coef
## 1   0  aRhIgG.PE.low             SIVsm.E660.84  0.32756399        var39  0.259346391
## 2   1 aRhIgG.PE.high            SIV.E543.gp140 -0.16850068       var182 -0.193650642
## 3   3          R2A.3                       G49  0.35430960       var599  0.189258364
## 4   0          R3A.3           SIVmac239.gp120 -0.36506237       var164 -0.163301695
## 5   4     R2A.4.high           SIVsmH4.p55.Gag -0.35612630       var807 -0.163179718
## 6   6  aRhIgG.PE.low                     C1.TR  0.41909923      var1054  0.131315107
## 7   3          R3A.3        SIVcpz.EK505.gp120  0.29138827       var679  0.117305209
## 8   2  aRhIgG.PE.low         SIVmac251.BK.PR55  0.36046966       var381  0.110575490
## 9   0      R2A.4.low            SIV.E543.gp140 -0.31307479       var130 -0.102131610
## 10  0          R3A.1           SIVmac239.gp130  0.21806357       var145  0.081708763
## 11  5          R2A.3        SIVcpz.EK505.gp120 -0.08457092       var950 -0.071004416
## 12  1      R2A.4.low                       G73 -0.27886344       var298 -0.066273789
## 13  0            C1q                       G73 -0.25642545        var47 -0.048460870
## 14  0          R2A.3            SIV.E543.gp140 -0.15810981        var89 -0.039156885
## 15  1          R3A.1 J08.V1V2.E660.2A5.AVI.His -0.21768334       var311 -0.038376913
## 16  2 aRhIgG.PE.high         SIVmac251.BK.PR55  0.32642465       var360  0.036170151
## 17  4          R2A.3             SIVsm.E660.84  0.25051885       var785  0.034116339
## 18  0      R2A.4.low   J08.V1V2.mac239.AVI.His -0.35496690       var129 -0.032109522
## 19  1          R3A.1            SIV.E543.gp140 -0.20186190       var314 -0.031542990
## 20  1          R3A.3                     C1.Ak -0.28032764       var324 -0.025323832
## 21  6            C1q            SIV.1A11.gp140 -0.17589498      var1083 -0.022776320
## 22  4            MBL   J08.V1V2.mac239.AVI.His  0.35978886       var758  0.015709812
## 23  1          R3A.1             SIVsm.E660.84  0.13616918       var321  0.008601751
## 24  0 aRhIgG.PE.high           SIVmac239.gp130  0.24117758        var13  0.006035165
## 25  3          R3A.1        SIVcpz.EK505.gp120 -0.21760688       var659 -0.002213295
## 26  6     R2A.4.high                       G73 -0.25505974      var1137 -0.001053021

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
dim(coef2)
## no variables

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 0.8354788


pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictSetpointVL/StraightGlmnet/Scaled")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

subData <- GetVariableSetData(fcData, coefSummary$shortVarName, TRUE, "LogSetpointVL")
predResults <- RunRandomPartitionPredictions(
    subData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="PreChallenge-glmnet-CompRmse.txt",
    fitsFile="PreChallenge-glmnet-Fits.rdata",
    verbose=TRUE, progressEvery=10)


resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
## 1.3696874 1.3583008 0.8730802 

## Model sizes range from 22 to 31 

## PickedCoeffs
##        var1054         var130         var145         var164         var182         var314 
##             50             50             50             50             50             50 
##         var321         var360          var39          var47         var599         var679 
##             50             50             50             50             50             50 
##         var807        var1083         var298         var311         var324         var785 
##             50             49             49             49             49             49 
## SIV_Mosaic_Env         var129          var89         var950         x_PARI         var381 
##             48             48             48             48             47             44 
##        SIV_Env         var659         var758        var1137          var13        SIV_Gag 
##             42             41             37             36             36             32 
