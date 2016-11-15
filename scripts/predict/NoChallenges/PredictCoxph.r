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

vlData <- GetTimepointData(fcData, 0:8, predSummary,
                           response="NoChallenges")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$y, family="cox")
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$y, family="cox")

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1[-1,"i"]-1]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=coef1$x[-1])


minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })
coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=glmnetRes2$alpha[minCvs == min(minCvs)])))
coefSummary2 <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef2[-1,"i"]-1]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=coef2$x[-1])

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(vlData$x),
                   s=glmnetRes$lambda.min) -
           vlData$y)^2))
## 2.499679

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(vlData$x),
                   alpha=glmnetRes2$alpha[minCvs==min(minCvs)]) -
           vlData$y)^2))
## 3.583404

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
