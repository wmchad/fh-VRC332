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

vlData <- GetTimepointData(deltaData, 1:8, predSummary,
                           response="NoChallenges")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=coef1$x[-1])
##    tp             re                      ag shortVarName        coef
## 1   1  aRhIgG.PE.low      SIVcpz.EK505.gp120       var204 -0.84993236
## 2   7  aRhIgG.PE.low          SIVsm.E660.2A5      var1242 -0.83853641
## 3   7  aRhIgG.PE.low J08.V1V2.mac239.AVI.His      var1233 -0.46184020
## 4   3          R3A.3      SIVcpz.EK505.gp120       var679 -0.41836333
## 5   4            C1q                    G119       var732  0.30219396
## 6   6 aRhIgG.PE.high                    G119      var1035 -0.28592024
## 7   3          R3A.3       SIVmac251.BK.PR55       var684 -0.28129118
## 8   5 aRhIgG.PE.high                    G119       var863  0.27927984
## 9   1          R2A.3         SIVmac239.gp120       var263 -0.22266643
## 10  4 aRhIgG.PE.high                   C1.TR       var690  0.20929366
## 11  2          R2A.3      SIVcpz.EK505.gp120       var434  0.18733617
## 12  7  aRhIgG.PE.low                     V1a      var1245 -0.17924497
## 13  2            C1q         SIVsmH4.p55.Gag       var404 -0.17640890
## 14  6     R2A.4.high      SIVcpz.EK505.gp120      var1143 -0.15199260
## 15  2          R3A.3          SIV.1A11.gp140       var505  0.14763619
## 16  4     R2A.4.high          SIV.1A11.gp140       var797  0.12986969
## 17  4            C1q                     G73       var735  0.11936014
## 18  1     R2A.4.high       SIVmac251.BK.PR55       var288  0.11924560
## 19  3          R3A.1          SIV.E543.gp140       var658  0.07043638
## 20  4      R2A.4.low          SIV.E543.gp140       var818 -0.05845978
## 21  7  aRhIgG.PE.low         SIVmac239.gp120      var1237 -0.05626338
## 22  3          R3A.3         SIVsmH4.p55.Gag       var687 -0.05121061
## 23  7            MBL                     G49      var1270 -0.04580186
## 24  6  aRhIgG.PE.low         SIVmac239.gp140      var1067  0.04495614
## 25  3          R3A.3          SIV.1A11.gp140       var677 -0.03167274
## 26  7     R2A.4.high                    G119      var1306 -0.02064295
## 27  2            C1q                     G73       var391 -0.01851624
## 28  3  aRhIgG.PE.low                   C1.TR       var538  0.01322479
## 29  1          R3A.1          SIV.E543.gp140       var314 -0.01035192

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
coefSummary2 <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef2[-1,"i"]-1]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=coef2$x[-1])
##   tp            re             ag shortVarName       coef
## 1  7 aRhIgG.PE.low SIVsm.E660.2A5      var1242 -0.8155072

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 2.009955

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)]) -
           vlData$y)^2))
## 3.389466

pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictNoChallenges/StraightGlmnet/Delta")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

vars <- (predSummary %>%
         filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
         select(shortVarName))$shortVarName
subData <- GetVariableSetData.delta(deltaData, vars, TRUE, "NoChallenges")
predResults <- RunRandomPartitionPredictions(
    subData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-CompRmse.txt",
    fitsFile="glmnet-Fits.rdata",
    verbose=TRUE, progressEvery=10)

vars2 <- (predSummary %>%
          filter(shortVarName %in% names(vlData$x)[coef2[-1,"i"]-1]) %>%
          select(shortVarName))$shortVarName
subData2 <- GetVariableSetData.delta(deltaData, vars2, TRUE, "NoChallenges")
predResults2 <- RunRandomPartitionPredictions(
    subData2, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-alpha-CompRmse.txt",
    fitsFile="glmnet-alpha-Fits.rdata",
    verbose=TRUE, progressEvery=10)

resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
##  3.760994  3.684869  2.098374 

## Model sizes range from 20 to 33 

## PickedCoeffs
##        var1233         var204         var434         var538         var684         var690 
##             50             50             50             50             50             50 
##         var735         var797         var818        var1035        var1306         var263 
##             50             50             50             49             49             49 
##         var288         var404         var679         var863         var677         var687 
##             49             49             49             49             48             48 
##        var1143         var658         var732        var1242        var1067         var314 
##             47             47             47             45             44             44 
##        var1245         var391 SIV_Mosaic_Env        var1237        SIV_Gag        SIV_Env 
##             43             43             42             38             37             36 
##         var505        var1270 
##             33             31 


resSummary2 <- SummarizePredictionResults(predResults2)
## Average RMSE: 
##  baseline     group predicted 
##  3.841224  3.708755  3.212188 

## Model sizes range from 2 to 6 

## PickedCoeffs
## var1242  x_PARI 
##      50      49 
