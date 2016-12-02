brrequire(glmnet)
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

deltaData <- GetDeltaData(fcData.na, predSummary)

vlData <- GetTimepointData(deltaData, 1:6, predSummary,
                           response="NoChallenges")

vlData$x[,-(1:4)] <- scale(vlData$x[,-(1:4)])

event <- rep(1, length(vlData$y))
event[vlData$y==13] <- 0

vlData$surv <- Surv(vlData$y, event)

glmnetRes <- cv.glmnet(x=as.matrix(vlData$x), y=vlData$surv, family="cox")
glmnetRes2 <- cvAlpha.glmnet(x=as.matrix(vlData$x), y=vlData$surv, family="cox")

coef1 <- as.data.frame(summary(predict(glmnetRes, type="coef", s=glmnetRes$lambda.min)))
coefSummary <- predSummary %>% filter(shortVarName %in% names(vlData$x)[coef1$i]) %>%
    select(tp:ag, shortVarName) %>% mutate(coef=as.numeric(coef1$x[-1]))
coefSummary %>% arrange(desc(abs(coef)))
##  tp            re              ag shortVarName        coef
## 1  6 aRhIgG.PE.low SIVmac239.gp140      var1067 -0.05523819
## 2  2         R3A.3  SIV.1A11.gp140       var505 -0.03026863



## Maybe how to do prediction?
vars <- coefSummary$shortVarName
subData <- GetVariableSetData(deltaData, vars, TRUE, "NoChallenges")
event <- rep(1, length(subData$y))
event[subData$y==13] <- 0
subData$surv <- Surv(subData$y, event)
coxModel <- coxph(subData$surv~as.matrix(subData$x[,-(1:3)]), init=coef1$x, iter=0)
coxFit <- survfit(coxModel, newData=as.matrix(subData$x[,-(1:3)]))
s0 <- exp(-coxFit$cumhaz)
preds <- sapply(1:nrow(subData$x), function(i) {
    sum(s0^exp(sum(subData$x[i,-(1:3)] * coef1$x)))
})
preds2 <- sapply(1:nrow(subData$x), function(i) {
    sum(s0^predict(glmnetRes,
                   newx=as.matrix(vlData$x[i,]),
                   s=glmnetRes$lambda.min,
                   type="response"))
})
rmse(preds, subData$y)
## 4.11666
rmse(preds2, subData$y)
## 4.150214
