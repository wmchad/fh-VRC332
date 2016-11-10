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

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
    select(tp:ag, corSetpoint, shortVarName) %>% mutate(coef=coef1$x[-1])

##    tp             re                      ag corSetpoint shortVarName         coef
## 1   8          R3A.3         SIVsmH4.p55.Gag -0.37932999      var1547 -0.435897004
## 2   8            C1q          SIV.E543.gp140  0.32260044      var1428  0.243614247
## 3   1 aRhIgG.PE.high          SIV.E543.gp140 -0.16850068       var182 -0.181053082
## 4   3          R2A.3                     G49  0.35430960       var599  0.178024611
## 5   8            C1q         SIVsmH4.p55.Gag -0.31273679      var1436 -0.159852626
## 6   8      R2A.4.low          SIV.E543.gp140  0.22154030      var1506  0.143315463
## 7   7     R2A.4.high                    G119 -0.08356839      var1306 -0.137448772
## 8   7 aRhIgG.PE.high         SIVsmH4.p55.Gag -0.09530865      var1223  0.133368420
## 9   2     R2A.4.high         SIVsmH4.p55.Gag -0.30699397       var463 -0.123689256
## 10  4          R3A.1      SIVcpz.EK505.gp120 -0.38429515       var831 -0.122801919
## 11  8  aRhIgG.PE.low                G145.146  0.30397207      var1400  0.118495825
## 12  7      R2A.4.low                   C1.TR -0.04115498      var1326 -0.116368158
## 13  0      R2A.4.low          SIV.E543.gp140 -0.31307479       var130 -0.103623683
## 14  0            C1q                     G73 -0.25642545        var47 -0.103037145
## 15  2  aRhIgG.PE.low       SIVmac251.BK.PR55  0.36046966       var381  0.100994422
## 16  0          R2A.3                   C1.Ak -0.06372115        var79 -0.096433105
## 17  6  aRhIgG.PE.low                   C1.TR  0.41909923      var1054  0.095095724
## 18  8          R2A.2         SIVmac239.gp140  0.24911016      var1454  0.091852404
## 19  0          R3A.3         SIVmac239.gp120 -0.36506237       var164 -0.089439651
## 20  1          R3A.1                    G119 -0.28853666       var307 -0.087582752
## 21  8          R2A.2 J08.V1V2.mac239.AVI.His  0.21596848      var1453  0.087186906
## 22  4      R2A.4.low J08.V1V2.mac239.AVI.His  0.15336267       var817  0.085347257
## 23  8          R2A.3       SIVmac251.BK.PR55 -0.30267812      var1471 -0.084804526
## 24  8 aRhIgG.PE.high       SIVmac251.BK.PR55 -0.31031129      var1392 -0.077174440
## 25  0          R2A.3          SIV.E543.gp140 -0.15810981        var89 -0.066999258
## 26  8     R2A.4.high      SIVcpz.EK505.gp120  0.24968574      var1487  0.066334365
## 27  7          R3A.1                G145.146  0.23309302      var1340  0.064847790
## 28  1     R2A.4.high          SIV.1A11.gp140 -0.13917344       var281 -0.061925951
## 29  8  aRhIgG.PE.low      SIVcpz.EK505.gp120  0.30231496      var1408  0.060376072
## 30  4            C1q         SIVsmH4.p55.Gag -0.21014862       var748 -0.058500907
## 31  3          R3A.1      SIVcpz.EK505.gp120 -0.21760688       var659 -0.054196893
## 32  8      R2A.4.low         SIVmac239.gp140  0.35455809      var1507  0.051630789
## 33  1          R3A.3                   C1.Ak -0.28032764       var324 -0.048044289
## 34  0  aRhIgG.PE.low           SIVsm.E660.84  0.32756399        var39  0.044566671
## 35  6 aRhIgG.PE.high                     G49  0.20607246      var1036  0.043276885
## 36  4          R2A.3                   C1.Ak  0.20631931       var767  0.042536643
## 37  5          R2A.3      SIVcpz.EK505.gp120 -0.08457092       var950 -0.042428631
## 38  0  aRhIgG.PE.low          SIV.E543.gp140  0.22307893        var31  0.037375740
## 39  8          R2A.2                G145.146  0.20901549      var1450  0.037175236
## 40  0          R3A.1         SIVmac239.gp130  0.21806357       var145  0.028759173
## 41  6            C1q       SIVmac251.BK.PR55 -0.17388383      var1089 -0.026658990
## 42  7            C1q                G145.146  0.16696670      var1249  0.024848138
## 43  4          R3A.3 SIVmac239.gp140.AVI.His -0.05507961       var855  0.019578933
## 44  7            C1q       SIVmac251.BK.PR55 -0.07777936      var1261 -0.019450978
## 45  7     R2A.4.high                     V1a -0.02661804      var1324 -0.019294021
## 46  1          R2A.3          SIV.E543.gp140 -0.24351392       var261 -0.018973061
## 47  6 aRhIgG.PE.high                    G119  0.26121178      var1035  0.018556494
## 48  7          R3A.3         SIVmac239.gp140  0.08710605      var1370 -0.016097569
## 49  8            C1q                    G119  0.14649234      var1420  0.015159862
## 50  4      R2A.4.low                   C1.Ak  0.11802566       var809  0.014824108
## 51  5 aRhIgG.PE.high                     G73  0.09588807       var865 -0.012588191
## 52  0          R2A.3       SIVmac251.BK.PR55  0.30735824        var95  0.011879879
## 53  4          R2A.3          SIVsm.E660.2A5  0.19552949       var784  0.011141477
## 54  3  aRhIgG.PE.low       SIVmac251.BK.PR55  0.29550883       var553  0.006463124
## 55  7 aRhIgG.PE.high                    G119 -0.03735070      var1207 -0.005386678
## 56  1          R3A.1                     G49 -0.28191082       var309 -0.004673842
## 57  4          R3A.1          SIV.1A11.gp140  0.04715944       var829  0.003084343

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs==min(minCvs)])))
dim(coef2)
## 54 x 3
coefSummary2 <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef2[-1,"i"]-1]) %>%
    select(tp:ag, corSetpoint, shortVarName) %>% mutate(coef=coef2$x[-1])

##    tp             re                      ag corSetpoint shortVarName          coef
## 1   8      R2A.4.low          SIV.E543.gp140   0.2215403      var1506  0.0743648741
## 2   8            C1q          SIV.E543.gp140   0.3226004      var1428  0.0684379888
## 3   8      R2A.4.low         SIVmac239.gp140   0.3545581      var1507  0.0617485758
## 4   8          R3A.3         SIVsmH4.p55.Gag  -0.3793300      var1547 -0.0578488497
## 5   8          R2A.3         SIVsmH4.p55.Gag  -0.4017975      var1474 -0.0538874326
## 6   8          R3A.1         SIVsmH4.p55.Gag  -0.3623140      var1526 -0.0507055569
## 7   8  aRhIgG.PE.low         SIVsmH4.p55.Gag  -0.2500835      var1416 -0.0483499442
## 8   8  aRhIgG.PE.low      SIVcpz.EK505.gp120   0.3023150      var1408  0.0468515362
## 9   8          R2A.2         SIVmac239.gp140   0.2491102      var1454  0.0445301713
## 10  6  aRhIgG.PE.low                   C1.TR   0.4190992      var1054  0.0425621059
## 11  4          R3A.1                G145.146  -0.4133622       var824 -0.0408941117
## 12  3          R2A.3                     G49   0.3543096       var599  0.0382716947
## 13  8            C1q         SIVsmH4.p55.Gag  -0.3127368      var1436 -0.0376491276
## 14  4     R2A.4.high         SIVsmH4.p55.Gag  -0.3561263       var807 -0.0367218956
## 15  2  aRhIgG.PE.low       SIVmac251.BK.PR55   0.3604697       var381  0.0366104457
## 16  4          R3A.1      SIVcpz.EK505.gp120  -0.3842951       var831 -0.0270909056
## 17  0          R3A.3         SIVmac239.gp120  -0.3650624       var164 -0.0240445857
## 18  0  aRhIgG.PE.low           SIVsm.E660.84   0.3275640        var39  0.0236469649
## 19  0      R2A.4.low J08.V1V2.mac239.AVI.His  -0.3549669       var129 -0.0228916199
## 20  0          R2A.3       SIVmac251.BK.PR55   0.3073582        var95  0.0224282342
## 21  8 aRhIgG.PE.high       SIVmac251.BK.PR55  -0.3103113      var1392 -0.0202108770
## 22  8  aRhIgG.PE.low          SIV.E543.gp140   0.4204668      var1407  0.0184964256
## 23  0      R2A.4.low          SIV.E543.gp140  -0.3130748       var130 -0.0171555413
## 24  8            C1q         SIVmac239.gp140   0.3407227      var1431  0.0165403838
## 25  3  aRhIgG.PE.low       SIVmac251.BK.PR55   0.2955088       var553  0.0159432785
## 26  6  aRhIgG.PE.low                   C1.Ak   0.3375813      var1053  0.0156895151
## 27  4            MBL J08.V1V2.mac239.AVI.His   0.3597889       var758  0.0140429726
## 28  3          R2A.3                    G119   0.3333186       var597  0.0132300759
## 29  8 aRhIgG.PE.high         SIVsmH4.p55.Gag  -0.3554162      var1395 -0.0124814903
## 30  8  aRhIgG.PE.low                G145.146   0.3039721      var1400  0.0116036366
## 31  6  aRhIgG.PE.low SIVmac239.gp140.AVI.His   0.3618778      var1068  0.0112488044
## 32  3          R2A.3                G145.146   0.2309037       var598  0.0104973084
## 33  6 aRhIgG.PE.high                   C1.Ak   0.3334741      var1033  0.0093933644
## 34  8     R2A.4.high      SIVcpz.EK505.gp120   0.2496857      var1487  0.0085220876
## 35  1          R3A.3                   C1.Ak  -0.2803276       var324 -0.0077933103
## 36  6     R2A.4.high                     G73  -0.2550597      var1137 -0.0058721422
## 37  0            C1q                     G73  -0.2564255        var47 -0.0054516269
## 38  1          R3A.1                     G49  -0.2819108       var309 -0.0054118486
## 39  1          R3A.1                    G119  -0.2885367       var307 -0.0048091964
## 40  6 aRhIgG.PE.high SIVmac239.gp140.AVI.His   0.3890780      var1047  0.0047294015
## 41  6     R2A.4.high                    G119  -0.2717097      var1134 -0.0043367555
## 42  0 aRhIgG.PE.high         SIVmac239.gp130   0.2411776        var13  0.0042631483
## 43  8          R2A.3       SIVmac251.BK.PR55  -0.3026781      var1471 -0.0041786882
## 44  2     R2A.4.high         SIVsmH4.p55.Gag  -0.3069940       var463 -0.0038893907
## 45  4          R3A.1                    G119  -0.3387945       var823 -0.0038658598
## 46  6     R2A.4.high                   C1.TR  -0.2961030      var1133 -0.0033709392
## 47  8  aRhIgG.PE.low                    G119   0.2821330      var1399  0.0027773519
## 48  8  aRhIgG.PE.low                   C1.Ak   0.2952612      var1397  0.0024649469
## 49  8  aRhIgG.PE.low                     G73   0.2483427      var1402  0.0022496115
## 50  2  aRhIgG.PE.low                   C1.TR   0.3531313       var366  0.0021533675
## 51  2  aRhIgG.PE.low          SIV.1A11.gp140   0.3328245       var374  0.0011153880
## 52  8  aRhIgG.PE.low         SIVmac239.gp130   0.3301959      var1410  0.0004218929
## 53  8  aRhIgG.PE.low         SIVmac239.gp140   0.5029202      var1411  0.0001154916

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 0.2911641

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)]) -
           vlData$y)^2))
## 0.9606868

pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictSetpointVL/StraightGlmnet/Scaled")
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
## 1.3630084 1.3376722 0.5148773 

## Model sizes range from 47 to 62 

## PickedCoeffs
##        var1036        var1223        var1261         var130        var1306        var1326 
##             50             50             50             50             50             50 
##        var1428        var1436        var1454        var1487        var1506        var1547 
##             50             50             50             50             50             50 
##         var164         var182         var261         var281         var381          var47 
##             50             50             50             50             50             50 
##         var599         var659         var817         var831          var89        var1054 
##             50             50             50             50             50             49 
##        var1089        var1392        var1453         var307         var463         var748 
##             49             49             49             49             49             49 
##         var855        var1207        var1249        var1370         var145        var1471 
##             49             48             48             48             48             48 
##        var1507          var31        var1450          var39          var79        var1035 
##             48             48             47             46             46             45 
##        var1340        var1400        var1408        var1420         var784         var553 
##             45             45             45             45             45             42 
##         var809         var324        SIV_Gag         var950         var865 SIV_Mosaic_Env 
##             42             41             40             40             38             35 
##          var95        SIV_Env        var1324         var829         var309         var767 
##             35             33             33             32             31             30 
##         x_PARI 
##             28 


resSummary2 <- SummarizePredictionResults(predResults2)
## Average RMSE: 
##  baseline     group predicted 
## 1.3342643 1.3183362 0.9182014 

## Model sizes range from 23 to 50 

## PickedCoeffs
##         var130        var1428        var1436        var1454         var324          var47 
##             50             50             50             50             50             50 
##        var1506        var1547         var129        var1487          var39         var599 
##             49             49             46             46             46             46 
##        var1471        var1054          var13        var1392        var1474         var758 
##             44             43             42             42             42             42 
##        var1507         var553         var309         var381         var831        var1400 
##             41             41             35             35             35             34 
##         var598         var463        var1526 SIV_Mosaic_Env         var164          var95 
##             34             32             30             29             29             29 
##        var1416         var824        var1407 
##             28             27             26 
