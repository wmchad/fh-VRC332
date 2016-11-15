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

vlData <- GetTimepointData(deltaData, 1:6, predSummary,
                           response="NoChallenges")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=coef1$x[-1])
coefSummary %>% arrange(desc(abs(coef)))
##    tp             re                        ag shortVarName         coef
## 1   1  aRhIgG.PE.low        SIVcpz.EK505.gp120       var204 -1.022442353
## 2   5 aRhIgG.PE.high                      G119       var863  0.699424678
## 3   4     R2A.4.high            SIV.1A11.gp140       var797  0.562566709
## 4   1          R2A.3           SIVsmH4.p55.Gag       var270  0.449945873
## 5   2            C1q           SIVsmH4.p55.Gag       var404 -0.400003971
## 6   3          R3A.1            SIV.E543.gp140       var658  0.390211961
## 7   3          R3A.3            SIV.1A11.gp140       var677 -0.381108885
## 8   4          R3A.3                       G73       var845  0.366183011
## 9   4            C1q                       G73       var735  0.324929638
## 10  2            C1q                       G73       var391 -0.322577165
## 11  3          R3A.3         SIVmac251.BK.PR55       var684 -0.306921775
## 12  1     R2A.4.high         SIVmac251.BK.PR55       var288  0.300917883
## 13  1          R2A.3           SIVmac239.gp120       var263 -0.296723536
## 14  3          R3A.3           SIVsmH4.p55.Gag       var687 -0.288304309
## 15  6          R3A.1        SIVcpz.EK505.gp120      var1175 -0.215707058
## 16  3            MBL                  G145.146       var581  0.213435171
## 17  5            C1q           SIVmac239.gp120       var913  0.153245282
## 18  6  aRhIgG.PE.low           SIVmac239.gp140      var1067  0.148314464
## 19  6          R3A.3                      G119      var1186 -0.094644322
## 20  6          R3A.1                      G119      var1167 -0.084368306
## 21  6            MBL   J08.V1V2.mac239.AVI.His      var1102  0.066723599
## 22  2          R3A.1                      G119       var479 -0.036892109
## 23  4      R2A.4.low            SIV.E543.gp140       var818 -0.030263989
## 24  3  aRhIgG.PE.low                     C1.TR       var538  0.029682388
## 25  3     R2A.4.high            SIV.E543.gp140       var626  0.029089381
## 26  3          R3A.3        SIVcpz.EK505.gp120       var679 -0.018447824
## 27  3     R2A.4.high   SIVmac239.gp140.AVI.His       var631  0.012387577
## 28  3          R2A.3   J08.V1V2.mac239.AVI.His       var603  0.010247846
## 29  3  aRhIgG.PE.low J08.V1V2.E660.2A5.AVI.His       var544  0.007398183
## 30  2          R3A.3            SIV.1A11.gp140       var505  0.003766673

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
## 1037 x 3

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 1.851807

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)]) -
           vlData$y)^2))
## 3.78923

pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictNoChallenges",
                    "StraightGlmnet/Delta/PreChallenge")
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

resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
##  3.778215  3.643187  2.153396 

## Model sizes range from 28 to 33 

## PickedCoeffs
##        var1186         var204         var270         var288         var391         var581 
##             50             50             50             50             50             50 
##         var626         var658         var677         var684         var687         var735 
##             50             50             50             50             50             50 
##         var797         var845         var863        var1167         var263         var404 
##             50             50             50             49             49             49 
##         var818         var538         var913        var1102         var679        var1175 
##             49             48             48             45             45             44 
##         var479        SIV_Env SIV_Mosaic_Env        var1067         var603         var544 
##             44             43             39             36             34             30 
##         var631         var505 
##             28             25 
