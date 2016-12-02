require(pryr)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "Modeling/LinearModels-Support.r"))

origData <- fcData.na[apply(fcData.na[,-(1:7)], 1, function(x){!any(x<0)}),]
adjData <- GetBaselineAdjustedData(origData, predSummary)
lagData <- GetDeltaData(origData, predSummary)

longData.orig <- ConvertToLongData(origData, predSummary) %>%
    filter(value >= 0) %>%
    select(AnimalId, GroupNm, tp, re, ag, value)

longData.adj <- ConvertToLongData(adjData, predSummary) %>%
    select(AnimalId, GroupNm, tp, re, ag, value)

longData.lag <- ConvertToLongData(lagData, predSummary) %>%
    select(AnimalId, GroupNm, tp, re, ag, value)

rm(origData)
rm(adjData)
rm(lagData)

model1 <- as.formula(value ~ GroupNm * as.factor(tp) + AnimalId * re + AnimalId * ag)
model2 <- as.formula(value ~ GroupNm * as.factor(tp) +
                         GroupNm * re + GroupNm * ag +
                         AnimalId * re + AnimalId * ag)
tps <- 1:8
grps <- unique(longData.orig$GroupNm)
animalIds <- unique(longData.orig$AnimalId)
ags <- unique(longData.orig$ag)
res <- unique(longData.orig$re)

##################################################
## First model
## See FirstModelSummaries.txt for results
##################################################

##------------------------------------------------
## Original Data
##------------------------------------------------
setwd("~/Projects/VRC332/Data/Modeling/Exploratory")
saveRDS(lm(model1, dat=longData.orig, model=FALSE),
        "model1.fit.orig.rds")

setwd("~/Projects/VRC332/Data/Modeling/Exploratory")
fit.orig <- readRDS("model1.fit.orig.rds")

resid.orig <- resid(fit.orig)
predict.orig <- predict(fit.orig)
coef.orig <- BuildFitCoefData(fit.orig, tps, grps, animalIds, ags, res)
rm(fit.orig)

longData.orig <- longData.orig %>% mutate(predicted=predict.orig,
                                          resid=resid.orig)
rm(predict.orig)
rm(resid.orig)

plotDir.orig <- "~/Projects/VRC332/Plots/Modeling/LinearModels/Orig"
PlotResids(longData.orig, ags, res, plotDir.orig)
PlotVsTruth(longData.orig, plotDir.orig)
setwd(plotDir.orig)
png("qqplot.png", height=600, width=600)
qqnorm(longData.orig$resid)
abline(0, 1, lty=2)
dev.off()



