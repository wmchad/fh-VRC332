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
model3 <- as.formula(value ~ GroupNm * as.factor(tp) + GroupNm * re + GroupNm * ag)

tps <- 1:8
grps <- unique(longData.orig$GroupNm)
animalIds <- unique(longData.orig$AnimalId)
ags <- unique(longData.orig$ag)
res <- unique(longData.orig$re)

dataDir <- "~/Projects/VRC332/Data/Modeling/Exploratory"
plotDir <- "~/Projects/VRC332/Plots/Modeling/LinearModels"

##################################################
## First model
## See FirstModelSummaries.txt for results
##################################################


##------------------------------------------------
## Original Data
##------------------------------------------------

results.orig <- RunAndSummarizeModel(longData.orig, "orig", model1, "model1",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)

ModelSummaryPlots(results.orig$data, file.path(plotDir, "Model1", "Orig"), ags, res)


##------------------------------------------------
## Baseline-Adjusted Data
##------------------------------------------------

results.adj <- RunAndSummarizeModel(longData.adj, "adj", model1, "model1",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)

ModelSummaryPlots(results.adj$data, file.path(plotDir, "Model1", "Adj"), ags, res)


##------------------------------------------------
## Lagged Data
##------------------------------------------------

results.lag <- RunAndSummarizeModel(longData.lag, "lag", model1, "model1",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)

ModelSummaryPlots(results.lag$data, file.path(plotDir, "Model1", "Lag"), ags, res)

##################################################
## Second model
## See SecondModelSummaries.txt for results
##################################################


##------------------------------------------------
## Original Data
##------------------------------------------------

results.orig2 <- RunAndSummarizeModel(longData.orig, "orig", model2, "model2",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)

ModelSummaryPlots(results.orig$data, file.path(plotDir, "Model2", "Orig"), ags, res)


##------------------------------------------------
## Baseline-Adjusted Data
##------------------------------------------------

longData.adj$tp <- as.factor(longData.adj$tp)
results.adj2 <- RunAndSummarizeModel(longData.adj, "adj", model2, "model2",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)

ModelSummaryPlots(results.adj2$data, file.path(plotDir, "Model2", "Adj"), ags, res)


##------------------------------------------------
## Lagged Data
##------------------------------------------------

longData.lag$tp <- as.factor(longData.adj$tp)
results.lag2 <- RunAndSummarizeModel(longData.lag, "lag", model2, "model2",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)

ModelSummaryPlots(results.lag2$data, file.path(plotDir, "Model2", "Lag"), ags, res)

##################################################
## Third model
## See ThirdModelSummaries.txt for results
##################################################


##------------------------------------------------
## Original Data
##------------------------------------------------

results.orig3 <- RunAndSummarizeModel(longData.orig, "orig", model3, "model3",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)

ModelSummaryPlots(results.orig3$data, file.path(plotDir, "Model3", "Orig"), ags, res)


##------------------------------------------------
## Baseline-Adjusted Data
##------------------------------------------------

longData.adj$tp <- as.factor(longData.adj$tp)
results.adj3 <- RunAndSummarizeModel(longData.adj, "adj", model3, "model3",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)

ModelSummaryPlots(results.adj3$data, file.path(plotDir, "Model3", "Adj"), ags, res)


##------------------------------------------------
## Lagged Data
##------------------------------------------------

longData.lag$tp <- as.factor(longData.adj$tp)
results.lag3 <- RunAndSummarizeModel(longData.lag, "lag", model3, "model3",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)

ModelSummaryPlots(results.lag3$data, file.path(plotDir, "Model3", "Lag"), ags, res)
