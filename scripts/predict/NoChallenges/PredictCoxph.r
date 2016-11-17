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
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1$i-1]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=as.numeric(coef1$x))
coefSummary %>% arrange(desc(abs(coef)))
##    tp             re                        ag shortVarName         coef
## 1   5          R3A.3                     C1.TR      var1013 -0.249159481
## 2   0  aRhIgG.PE.low            SIV.E543.gp140        var31 -0.156534056
## 3   7          R3A.1             SIVsm.E660.84      var1353  0.148413418
## 4   8            C1q            SIV.1A11.gp140      var1427  0.138364310
## 5   7  aRhIgG.PE.low            SIV.1A11.gp140      var1234  0.137292435
## 6   1     R2A.4.high J08.V1V2.E660.2A5.AVI.His       var279  0.105786372
## 7   5          R3A.3                     C1.Ak      var1012 -0.080875241
## 8   5          R3A.1            SIV.E543.gp140      var1002 -0.050433619
## 9   7 aRhIgG.PE.high   SIVmac239.gp140.AVI.His      var1219 -0.049992260
## 10  8 aRhIgG.PE.high   SIVmac239.gp140.AVI.His      var1391 -0.039623242
## 11  3          R3A.3           SIVmac239.gp140       var682  0.029190090
## 12  5          R3A.1            SIVsm.E660.2A5      var1008 -0.024624666
## 13  5          R3A.3                       G49      var1016 -0.023656457
## 14  1          R2A.3        SIVcpz.EK505.gp120       var262  0.021476501
## 15  3          R3A.3            SIV.E543.gp140       var678  0.008011964
## 16  2          R2A.2   J08.V1V2.mac239.AVI.His       var421  0.006003471
## 17  6     R2A.4.high            SIV.E543.gp140      var1142  0.004563079

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
coefSummary2 <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef2$i-1]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=as.numeric(coef2$x))
coefSummary2 %>% arrange(desc(abs(coef)))
##   tp            re                        ag shortVarName        coef
## 1   7         R3A.1             SIVsm.E660.84      var1353  0.15559746
## 2   5         R3A.3                     C1.TR      var1013 -0.15523835
## 3   0 aRhIgG.PE.low            SIV.E543.gp140        var31 -0.06410573
## 4   5         R3A.1            SIV.E543.gp140      var1002 -0.04923751
## 5   7 aRhIgG.PE.low            SIV.1A11.gp140      var1234  0.04599086
## 6   8           C1q            SIV.1A11.gp140      var1427  0.04562859
## 7   1    R2A.4.high J08.V1V2.E660.2A5.AVI.His       var279  0.03680156
## 8   5         R3A.3                     C1.Ak      var1012 -0.02133747
## 9   8         R3A.1            SIVsm.E660.2A5      var1524  0.01666741
## 10  5         R3A.3                       G49      var1016 -0.01472086
## 11  5     R2A.4.low            SIV.E543.gp140       var990 -0.01171111
## 12  8    R2A.4.high            SIV.1A11.gp140      var1485  0.00222758

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min, type="response") -
           vlData$y)^2))
## 5.687537

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)],
                   type="response") -
           vlData$y)^2))
## 5.606923

pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictNoChallenges/StraightGlmnet/Delta")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

vars <- (predSummary %>%
         filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
         select(shortVarName))$shortVarName
subData <- GetVariableSetData(deltaData, vars, TRUE, "NoChallenges")
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
