require(glmnet)
require(glmnetUtils)
require(dplyr)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))

fcData <- fcData %>% filter(NoChallenges >0)

noChallengesData <- GetTimepointData(fcData, 0:8, predSummary,
                                     response="NoChallenges")

glmnetRes <- cv.glmnet(x=as.matrix(noChallengesData$x), y=noChallengesData$y)
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(noChallengesData$x), y=noChallengesData$y)

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
predSummary %>% filter(shortVarName %in% c("var32", "var863"))
##   tp             re                 ag numObs numObs2     corPeak corSetpoint shortVarName
## 1  0  aRhIgG.PE.low SIVcpz.EK505.gp120     69      69 0.001888941   0.1506623        var32
## 2  5 aRhIgG.PE.high               G119     87      87 0.070262290   0.1851554       var863


minCvs <- sapply(glmnetRes2$modlist, function(x) { min(x$cvm) })

coef2 <- as.data.frame(summary(coef(glmnetRes2, alpha=1)))

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(noChallengesData$x),
                   s=glmnetRes$lambda.min) -
           noChallengesData$y)^2))
## 2.517042

sqrt(mean((predict(glmnetRes2,
                   newx=as.matrix(noChallengesData$x),
                   alpha=1) -
           noChallengesData$y)^2))
## 2.734363

sqrt(mean((predict(glmnetRes,
                   newx=as.matrix(noChallengesData$x),
                   s=glmnetRes$lambda.1se) -
           noChallengesData$y)^2))
## 2.881858
