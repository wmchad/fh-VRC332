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
model4 <- as.formula(value ~ GroupNm * as.factor(tp) +
                         GroupNm * re + GroupNm * ag +
                         AnimalId) 
model5 <- as.formula(value ~ re * ag)
model6 <- as.formula(value ~ re * ag * as.factor(tp))
model7 <- as.formula(value ~ re * ag * as.factor(tp) + AnimalId * as.factor(tp))
model8 <- as.formula(value ~ re * ag * as.factor(tp) + GroupNm * as.factor(tp))
model9 <- as.formula(value ~ re * ag * as.factor(tp) * GroupNm)

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

##################################################
## Fourth model
## See FourthModelSummaries.txt for results
##################################################


##------------------------------------------------
## Original Data
##------------------------------------------------

results.orig4 <- RunAndSummarizeModel(longData.orig, "orig", model4, "model4",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.orig4$data, file.path(plotDir, "Model4", "Orig"), ags, res)

##------------------------------------------------
## Baseline-Adjusted Data
##------------------------------------------------

results.adj4 <- RunAndSummarizeModel(longData.adj, "adj", model4, "model4",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.adj4$data, file.path(plotDir, "Model4", "Adj"), ags, res)

##------------------------------------------------
## Lagged Data
##------------------------------------------------

results.lag4 <- RunAndSummarizeModel(longData.lag, "lag", model4, "model4",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.lag4$data, file.path(plotDir, "Model4", "Lag"), ags, res)


##################################################
## Fifth model
## See FifthModelSummaries.txt for results
##################################################


##------------------------------------------------
## Original Data
##------------------------------------------------

results.orig5 <- RunAndSummarizeModel(longData.orig, "orig", model5, "model5",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.orig5$data, file.path(plotDir, "Model5", "Orig"), ags, res)

##------------------------------------------------
## Baseline-Adjusted Data
##------------------------------------------------

results.adj5 <- RunAndSummarizeModel(longData.adj, "adj", model5, "model5",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.adj5$data, file.path(plotDir, "Model5", "Adj"), ags, res)

##------------------------------------------------
## Lagged Data
##------------------------------------------------

results.lag5 <- RunAndSummarizeModel(longData.lag, "lag", model5, "model5",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.lag5$data, file.path(plotDir, "Model5", "Lag"), ags, res)

##################################################
## Sixth model
## See SixthModelSummaries.txt for results
##################################################

##------------------------------------------------
## Original Data
##------------------------------------------------

results.orig6 <- RunAndSummarizeModel(longData.orig, "orig", model6, "model6",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.orig6$data, file.path(plotDir, "Model6", "Orig"), ags, res)

##------------------------------------------------
## Baseline-Adjusted Data
##------------------------------------------------

results.adj6 <- RunAndSummarizeModel(longData.adj, "adj", model6, "model6",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.adj6$data, file.path(plotDir, "Model6", "Adj"), ags, res)

##------------------------------------------------
## Lagged Data
##------------------------------------------------

results.lag6 <- RunAndSummarizeModel(longData.lag, "lag", model6, "model6",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.lag6$data, file.path(plotDir, "Model6", "Lag"), ags, res)


##################################################
## Seventh model
## See SeventhModelSummaries.txt for results
##################################################


##------------------------------------------------
## Original Data
##------------------------------------------------

results.orig7 <- RunAndSummarizeModel(longData.orig, "orig", model7, "model7",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.orig7$data, file.path(plotDir, "Model7", "Orig"), ags, res)

##------------------------------------------------
## Baseline-Adjusted Data
##------------------------------------------------

results.adj7 <- RunAndSummarizeModel(longData.adj, "adj", model7, "model7",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.adj7$data, file.path(plotDir, "Model7", "Adj"), ags, res)

##------------------------------------------------
## Lagged Data
##------------------------------------------------

results.lag7 <- RunAndSummarizeModel(longData.lag, "lag", model7, "model7",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.lag7$data, file.path(plotDir, "Model7", "Lag"), ags, res)


##################################################
## Eighth model
## See EighthModelSummaries.txt for results
##################################################


##------------------------------------------------
## Original Data
##------------------------------------------------

results.orig8 <- RunAndSummarizeModel(longData.orig, "orig", model8, "model8",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.orig8$data, file.path(plotDir, "Model8", "Orig"), ags, res)

##------------------------------------------------
## Baseline-Adjusted Data
##------------------------------------------------

results.adj8 <- RunAndSummarizeModel(longData.adj, "adj", model8, "model8",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.adj8$data, file.path(plotDir, "Model8", "Adj"), ags, res)

##------------------------------------------------
## Lagged Data
##------------------------------------------------

results.lag8 <- RunAndSummarizeModel(longData.lag, "lag", model8, "model8",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.lag8$data, file.path(plotDir, "Model8", "Lag"), ags, res)


##################################################
## Ninth model
## See NinthModelSummaries.txt for results
##################################################


##------------------------------------------------
## Original Data
##------------------------------------------------

results.orig9 <- RunAndSummarizeModel(longData.orig, "orig", model9, "model9",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.orig9$data, file.path(plotDir, "Model9", "Orig"), ags, res)

##------------------------------------------------
## Baseline-Adjusted Data
##------------------------------------------------

results.adj9 <- RunAndSummarizeModel(longData.adj, "adj", model9, "model9",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.adj9$data, file.path(plotDir, "Model9", "Adj"), ags, res)

##------------------------------------------------
## Lagged Data
##------------------------------------------------

results.lag9 <- RunAndSummarizeModel(longData.lag, "lag", model9, "model9",
                                     dataDir,
                                     tps, grps, animalIds, ags, res)
ModelSummaryPlots(results.lag9$data, file.path(plotDir, "Model9", "Lag"), ags, res)

