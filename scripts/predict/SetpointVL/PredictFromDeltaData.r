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
                           response="LogSetpointVL")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
    select(tp:ag, corSetpoint, shortVarName) %>% mutate(coef=coef1$x[-1])
##    tp             re                      ag corSetpoint shortVarName         coef
## 1   8          R3A.3         SIVsmH4.p55.Gag -0.37932999      var1547 -0.420924787
## 2   3          R2A.3                     G49  0.35430960       var599  0.192001163
## 3   8            C1q          SIV.E543.gp140  0.32260044      var1428  0.188641406
## 4   8      R2A.4.low          SIV.E543.gp140  0.22154030      var1506  0.182093415
## 5   8            C1q         SIVsmH4.p55.Gag -0.31273679      var1436 -0.175287987
## 6   8      R2A.4.low         SIVmac239.gp140  0.35455809      var1507  0.156443953
## 7   6  aRhIgG.PE.low                   C1.TR  0.41909923      var1054  0.150445479
## 8   1 aRhIgG.PE.high          SIV.E543.gp140 -0.16850068       var182 -0.146244015
## 9   4          R3A.1      SIVcpz.EK505.gp120 -0.38429515       var831 -0.130772532
## 10  7      R2A.4.low                   C1.TR -0.04115498      var1326 -0.122110015
## 11  2  aRhIgG.PE.low       SIVmac251.BK.PR55  0.36046966       var381  0.118211395
## 12  7     R2A.4.high                    G119 -0.08356839      var1306 -0.115325705
## 13  8          R2A.2         SIVmac239.gp140  0.24911016      var1454  0.115140979
## 14  8  aRhIgG.PE.low                G145.146  0.30397207      var1400  0.093352449
## 15  4     R2A.4.high         SIVsmH4.p55.Gag -0.35612630       var807 -0.090498716
## 16  4      R2A.4.low J08.V1V2.mac239.AVI.His  0.15336267       var817  0.088641583
## 17  1          R3A.1                    G119 -0.28853666       var307 -0.083136287
## 18  4            C1q         SIVsmH4.p55.Gag -0.21014862       var748 -0.074564548
## 19  8          R2A.3       SIVmac251.BK.PR55 -0.30267812      var1471 -0.072616485
## 20  8          R2A.2                G145.146  0.20901549      var1450  0.061459724
## 21  5          R2A.3      SIVcpz.EK505.gp120 -0.08457092       var950 -0.061171291
## 22  7          R3A.1                G145.146  0.23309302      var1340  0.055907770
## 23  1          R3A.3                   C1.Ak -0.28032764       var324 -0.054338895
## 24  8          R2A.2 J08.V1V2.mac239.AVI.His  0.21596848      var1453  0.052880714
## 25  3          R3A.1      SIVcpz.EK505.gp120 -0.21760688       var659 -0.043034859
## 26  3 aRhIgG.PE.high                     G73  0.03991239       var521 -0.042630681
## 27  1          R3A.1                     G49 -0.28191082       var309 -0.036348493
## 28  4          R3A.1                G145.146 -0.41336217       var824 -0.034050623
## 29  8     R2A.4.high      SIVcpz.EK505.gp120  0.24968574      var1487  0.032358293
## 30  7            MBL                G145.146 -0.13009980      var1269 -0.017929234
## 31  1          R3A.1          SIV.E543.gp140 -0.20186190       var314 -0.017014690
## 32  1     R2A.4.high          SIV.1A11.gp140 -0.13917344       var281 -0.016150580
## 33  3          R2A.3                G145.146  0.23090366       var598  0.012987061
## 34  7 aRhIgG.PE.high         SIVsmH4.p55.Gag -0.09530865      var1223  0.010872248
## 35  6            C1q          SIV.E543.gp140 -0.17013389      var1084 -0.006388640
## 36  3          R3A.3      SIVcpz.EK505.gp120  0.29138827       var679  0.004721512
## 37  6            C1q       SIVmac251.BK.PR55 -0.17388383      var1089 -0.004277747
## 38  4          R2A.3                   C1.Ak  0.20631931       var767  0.000776197

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
coefSummary2 <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef2[-1,"i"]-1]) %>%
    select(tp:ag, corSetpoint, shortVarName) %>% mutate(coef=coef2$x[-1])
##    tp             re                      ag corSetpoint shortVarName          coef
## 1   8      R2A.4.low          SIV.E543.gp140   0.2215403      var1506  9.180849e-02
## 2   8            C1q          SIV.E543.gp140   0.3226004      var1428  8.517200e-02
## 3   8      R2A.4.low         SIVmac239.gp140   0.3545581      var1507  7.765092e-02
## 4   8          R3A.3         SIVsmH4.p55.Gag  -0.3793300      var1547 -7.392201e-02
## 5   8          R2A.3         SIVsmH4.p55.Gag  -0.4017975      var1474 -6.584112e-02
## 6   8          R3A.1         SIVsmH4.p55.Gag  -0.3623140      var1526 -6.012964e-02
## 7   6  aRhIgG.PE.low                   C1.TR   0.4190992      var1054  5.965096e-02
## 8   8  aRhIgG.PE.low      SIVcpz.EK505.gp120   0.3023150      var1408  5.671478e-02
## 9   4          R3A.1                G145.146  -0.4133622       var824 -5.664255e-02
## 10  8  aRhIgG.PE.low         SIVsmH4.p55.Gag  -0.2500835      var1416 -5.501605e-02
## 11  2  aRhIgG.PE.low       SIVmac251.BK.PR55   0.3604697       var381  5.284031e-02
## 12  3          R2A.3                     G49   0.3543096       var599  4.951381e-02
## 13  8          R2A.2         SIVmac239.gp140   0.2491102      var1454  4.933403e-02
## 14  4     R2A.4.high         SIVsmH4.p55.Gag  -0.3561263       var807 -4.373330e-02
## 15  8            C1q         SIVsmH4.p55.Gag  -0.3127368      var1436 -3.339417e-02
## 16  4          R3A.1      SIVcpz.EK505.gp120  -0.3842951       var831 -3.226406e-02
## 17  3  aRhIgG.PE.low       SIVmac251.BK.PR55   0.2955088       var553  1.461241e-02
## 18  6  aRhIgG.PE.low                   C1.Ak   0.3375813      var1053  1.440895e-02
## 19  3          R2A.3                    G119   0.3333186       var597  1.251230e-02
## 20  4            MBL J08.V1V2.mac239.AVI.His   0.3597889       var758  1.007775e-02
## 21  8 aRhIgG.PE.high       SIVmac251.BK.PR55  -0.3103113      var1392 -8.102323e-03
## 22  6 aRhIgG.PE.high                   C1.Ak   0.3334741      var1033  7.278978e-03
## 23  8            C1q         SIVmac239.gp140   0.3407227      var1431  7.190711e-03
## 24  8  aRhIgG.PE.low          SIV.E543.gp140   0.4204668      var1407  4.229972e-03
## 25  6  aRhIgG.PE.low SIVmac239.gp140.AVI.His   0.3618778      var1068  1.908959e-03
## 26  3          R2A.3                G145.146   0.2309037       var598  1.825875e-03
## 27  6     R2A.4.high                     G73  -0.2550597      var1137 -7.093354e-04
## 28  1          R3A.3                   C1.Ak  -0.2803276       var324 -9.868165e-05

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 0.4681463

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)]) -
           vlData$y)^2))
## 0.9790285


pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictSetpointVL/StraightGlmnet/Delta")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

vars <- (predSummary %>%
         filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
         select(shortVarName))$shortVarName
subData <- GetVariableSetData(deltaData, vars, TRUE, "LogSetpointVL")
predResults <- RunRandomPartitionPredictions(
    subData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-CompRmse.txt",
    fitsFile="glmnet-Fits.rdata",
    verbose=TRUE, progressEvery=10)

vars2 <- (predSummary %>%
          filter(shortVarName %in% names(vlData$x)[coef2[-1,"i"]-1]) %>%
          select(shortVarName))$shortVarName
subData2 <- GetVariableSetData(deltaData, vars2, TRUE, "LogSetpointVL")
predResults2 <- RunRandomPartitionPredictions(
    subData2, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-alpha-CompRmse.txt",
    fitsFile="glmnet-alpha-Fits.rdata",
    verbose=TRUE, progressEvery=10)

resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
##  1.415096  1.421355  1.077553 

## Model sizes range from 21 to 43 

## PickedCoeffs
##        SIV_Gag        var1084        var1340        var1450        var1471         var314 
##             50             50             50             50             50             50 
##         var679         var824         var659         var817        var1428        var1547 
##             50             50             49             49             48             48 
## SIV_Mosaic_Env        var1054        var1487         var598         var767         x_PARI 
##             47             46             46             46             46             46 
##        var1506         var182         var324         var381         var748        var1306 
##             45             45             45             45             45             44 
##         var599        var1436        var1454         var521        var1269         var307 
##             44             43             42             42             41             41 
##        var1223        var1400        var1453         var309        var1326         var281 
##             40             37             37             35             34             34 
##         var807         var831        var1507         var950        var1089        SIV_Env 
##             33             32             31             29             27             26 

resSummary2 <- SummarizePredictionResults(predResults2)
## Average RMSE: 
##  baseline     group predicted 
##  1.407875  1.414864  1.291528 

## Model sizes range from 11 to 33 

## PickedCoeffs
##        SIV_Gag        var1407        var1416         x_PARI        var1068        var1474 
##             50             50             50             50             49             49 
##        var1428         var599         var807         var381 SIV_Mosaic_Env        var1506 
##             48             48             48             47             46             46 
##         var324         var598         var758         var824        var1392        var1431 
##             46             46             43             43             42             42 
##         var831        var1507         var597         var553        var1054        var1454 
##             42             41             41             40             37             37 
##        var1137        var1547        var1526        var1436        var1033        SIV_Env 
##             36             34             30             29             28             26 
##        var1408 
##             25 
