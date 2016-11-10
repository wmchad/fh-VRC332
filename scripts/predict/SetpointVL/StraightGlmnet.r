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
                           response="LogSetpointVL")

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1])
##    tp             re                      ag numObs numObs2 corSetpoint shortVarName
## 1   0  aRhIgG.PE.low          SIV.E543.gp140     69      69  0.22307893        var31
## 2   0  aRhIgG.PE.low           SIVsm.E660.84     69      69  0.32756399        var39
## 3   0            C1q                     G73     69      69 -0.25642545        var47
## 4   0          R2A.3                   C1.Ak     69      69 -0.06372115        var79
## 5   0          R2A.3          SIV.E543.gp140     69      69 -0.15810981        var89
## 6   0      R2A.4.low          SIV.E543.gp140     69      69 -0.31307479       var130
## 7   0          R3A.1         SIVmac239.gp130     69      69  0.21806357       var145
## 8   0          R3A.3         SIVmac239.gp120     69      69 -0.36506237       var164
## 9   1 aRhIgG.PE.high          SIV.E543.gp140     70      70 -0.16850068       var182
## 10  1     R2A.4.high          SIV.1A11.gp140     69      69 -0.13917344       var281
## 11  1          R3A.1                    G119     70      70 -0.28853666       var307
## 12  1          R3A.1                     G49     70      70 -0.28191082       var309
## 13  1          R3A.3                   C1.Ak     70      70 -0.28032764       var324
## 14  1          R3A.3          SIVsm.E660.2A5     70      70 -0.25016767       var341
## 15  2  aRhIgG.PE.low       SIVmac251.BK.PR55     70      70  0.36046966       var381
## 16  2     R2A.4.high         SIVsmH4.p55.Gag     70      70 -0.30699397       var463
## 17  3 aRhIgG.PE.high                     G73     70      70  0.03991239       var521
## 18  3          R2A.3                     G49     70      70  0.35430960       var599
## 19  3          R3A.1      SIVcpz.EK505.gp120     70      70 -0.21760688       var659
## 20  3          R3A.3                     G73     70      70  0.02315775       var673
## 21  4            C1q         SIVsmH4.p55.Gag     70      70 -0.21014862       var748
## 22  4          R2A.3                   C1.Ak     70      70  0.20631931       var767
## 23  4          R2A.3          SIVsm.E660.2A5     70      70  0.19552949       var784
## 24  4      R2A.4.low J08.V1V2.mac239.AVI.His     70      70  0.15336267       var817
## 25  4          R3A.1                G145.146     70      70 -0.41336217       var824
## 26  4          R3A.1          SIV.1A11.gp140     70      70  0.04715944       var829
## 27  4          R3A.1      SIVcpz.EK505.gp120     70      70 -0.38429515       var831
## 28  5          R2A.3      SIVcpz.EK505.gp120     87      87 -0.08457092       var950
## 29  6 aRhIgG.PE.high                     G49     70      70  0.20607246      var1036
## 30  6  aRhIgG.PE.low                   C1.TR     70      70  0.41909923      var1054
## 31  6            C1q       SIVmac251.BK.PR55     70      70 -0.17388383      var1089
## 32  7 aRhIgG.PE.high         SIVsmH4.p55.Gag     87      87 -0.09530865      var1223
## 33  7            C1q       SIVmac251.BK.PR55     87      87 -0.07777936      var1261
## 34  7     R2A.4.high                    G119     87      87 -0.08356839      var1306
## 35  7     R2A.4.high                     V1a     87      87 -0.02661804      var1324
## 36  7      R2A.4.low                   C1.TR     87      87 -0.04115498      var1326
## 37  7          R3A.1                G145.146     87      87  0.23309302      var1340
## 38  8 aRhIgG.PE.high       SIVmac251.BK.PR55     87      87 -0.31031129      var1392
## 39  8  aRhIgG.PE.low                G145.146     87      87  0.30397207      var1400
## 40  8  aRhIgG.PE.low      SIVcpz.EK505.gp120     87      87  0.30231496      var1408
## 41  8            C1q          SIV.E543.gp140     87      87  0.32260044      var1428
## 42  8            C1q         SIVsmH4.p55.Gag     87      87 -0.31273679      var1436
## 43  8          R2A.2                G145.146     87      87  0.20901549      var1450
## 44  8          R2A.2 J08.V1V2.mac239.AVI.His     87      87  0.21596848      var1453
## 45  8          R2A.2         SIVmac239.gp140     87      87  0.24911016      var1454
## 46  8          R2A.3       SIVmac251.BK.PR55     87      87 -0.30267812      var1471
## 47  8     R2A.4.high      SIVcpz.EK505.gp120     87      87  0.24968574      var1487
## 48  8      R2A.4.low          SIV.E543.gp140     87      87  0.22154030      var1506
## 49  8      R2A.4.low         SIVmac239.gp140     87      87  0.35455809      var1507
## 50  8          R3A.3         SIVsmH4.p55.Gag     87      87 -0.37932999      var1547

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs==min(minCvs)])))
dim(coef2)
## 165 x 3

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 0.382113

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)]) -
           vlData$y)^2))
## 0.5455156

pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictSetpointVL/StraightGlmnet")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

vars <- (predSummary %>%
         filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
         select(shortVarName))$shortVarName
subData <- GetVariableSetData(fcData, vars, TRUE, "LogSetpointVL")
predResults <- RunRandomPartitionPredictions(
    subData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-CompRmse.txt",
    fitsFile="glmnet-Fits.rdata",
    verbose=TRUE, progressEvery=10)

vars2 <- (predSummary %>%
         filter(shortVarName %in% names(vlData$x)[coef2[-1,"i"]-1]) %>%
         select(shortVarName))$shortVarName
subData2 <- GetVariableSetData(fcData, vars2, TRUE, "LogSetpointVL")
predResults2 <- RunRandomPartitionPredictions(
    subData2, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-alpha-CompRmse.txt",
    fitsFile="glmnet-alpha-Fits.rdata",
    verbose=TRUE, progressEvery=10)

resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
## 1.3466849 1.3322638 0.5264561 

## Model sizes range from 40 to 55 

## PickedCoeffs
##        var1036        var1054        var1223         var130        var1306        var1326 
##             50             50             50             50             50             50 
##        var1428        var1454        var1547         var164         var182         var281 
##             50             50             50             50             50             50 
##         var307         var463          var47         var599          var79         var831 
##             50             50             50             50             50             50 
##        var1340        var1392        var1400        var1408        var1436        var1487 
##             49             49             49             49             49             49 
##        var1506         var381         var659         var817          var89        var1324 
##             49             49             49             49             49             48 
##        var1453          var39        var1471          var31         var767        var1089 
##             48             48             46             46             46             45 
##        var1261        var1450         var341         var145         var748         var950 
##             44             44             44             43             43             43 
##        var1507        SIV_Gag         var784         var829         var673         var324 
##             42             41             39             36             35             29 
##         var521        SIV_Env         var309 SIV_Mosaic_Env 
##             29             28             28             25 

resSummary2 <- SummarizePredictionResults(predResults2)
## Average RMSE: 
##  baseline     group predicted 
## 1.3276296 1.3129018 0.6873014 

## Model sizes range from 101 to 125 

## PickedCoeffs
## var1036 var1223 var1306 var1375 var1400 var1408 var1416 var1436 var1454 var1487 var1506 var1519 
##      50      50      50      50      50      50      50      50      50      50      50      50 
## var1526 var1547  var164  var599 var1326 var1392 var1474  var307   var39 var1309 var1450  var182 
##      50      50      50      50      49      49      49      49      49      48      48      48 
## var1054  var130 var1341 var1359 var1428  var324   var47  var472 var1135  var129 var1397 var1471 
##      47      47      47      47      47      47      47      47      46      46      46      46 
## var1495 var1534  var463   var89 var1286 var1297  var145 var1535  var818 var1205 var1370 var1407 
##      46      46      46      46      45      45      45      45      45      44      44      44 
## var1507   var78  var784  var817 var1172   var27   var31  var341  var807   var93  var281  var635 
##      44      44      44      44      43      43      43      43      43      43      42      42 
##   var79 var1280 var1420  var169 var1089 var1431  var309  var679  var778  var831  var855  var917 
##      42      41      41      41      40      40      40      40      40      40      40      40 
##  var261  var311  var598  var824  var950 var1402  var381   var68  var806 var1278  var342  var993 
##      39      39      39      39      39      38      38      38      38      36      36      36 
## var1083 var1134 var1451  var798   var95 var1034  var814  var314  var473  var919 var1461 SIV_Gag 
##      35      35      35      35      35      34      34      33      33      33      32      31 
##  var748  var785  var809 var1350 var1395 var1399 var1035  var581 var1033 var1453  var595  var659 
##      31      31      31      30      30      30      29      29      28      28      28      28 
##  var758 var1149  var142  var789  var840 var1316 var1502  var293 var1047  var191  var521  var710 
##      28      27      27      27      27      26      26      26      25      25      25      25 
