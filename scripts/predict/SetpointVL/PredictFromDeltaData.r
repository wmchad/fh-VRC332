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
                           response="LogSetpointVL")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=coef1$x[-1])
##   tp            re             ag shortVarName        coef
## 1  7           C1q SIV.1A11.gp140      var1255 0.073701752
## 2  3         R2A.3          C1.Ak       var595 0.073530380
## 3  2 aRhIgG.PE.low          C1.TR       var366 0.022384058
## 4  2 aRhIgG.PE.low       G145.146       var368 0.009767604

minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
## 1 x 3

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 1.336225

pTrain <- 0.75
nRand <- 50
ctrl <- trainControl(method="repeatedcv", repeats=20)
outdir <- file.path("~/Documents/Projects/MonkeySIV/Data/PredictSetpointVL/StraightGlmnet/Delta")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

vars <- (predSummary %>%
         filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
         select(shortVarName))$shortVarName
subData <- GetVariableSetData.delta(deltaData, vars, TRUE, "LogSetpointVL")
predResults <- RunRandomPartitionPredictions(
    subData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
    outdir=outdir,
    rmseFile="glmnet-CompRmse.txt",
    fitsFile="glmnet-Fits.rdata",
    verbose=TRUE, progressEvery=10)

resSummary <- SummarizePredictionResults(predResults)
## Average RMSE: 
##  baseline     group predicted 
##  1.412480  1.376615  1.255242 

## Model sizes range from 5 to 9 

## PickedCoeffs
## var1255  var368  var595  var366 SIV_Gag SIV_Env 
##      50      50      50      49      28      26 
