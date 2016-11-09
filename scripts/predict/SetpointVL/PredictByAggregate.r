require(caret)
require(glmnet)
require(reshape2)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetBestPredictorsData.r"))
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))
source(file.path(fnFolder, "ConvertToLongData.r"))

longData <- ConvertToLongData(fcData, predSummary)

aggAgData <- longData %>%
    mutate(tp=paste("tp", tp, sep="")) %>%
    unite(aggName, tp, ag) %>%
    group_by(AnimalId, aggName) %>%
    summarize(avg=mean(value)) %>%
    spread(aggName, avg)
fcAgData <- fcData[,1:7] %>% left_join(aggAgData)

aggReData <- longData %>%
    mutate(tp=paste("tp", tp, sep="")) %>%
    unite(aggName, tp, re) %>%
    group_by(AnimalId, aggName) %>%
    summarize(avg=mean(value)) %>%
    spread(aggName, avg)
fcReData <- fcData[,1:7] %>% left_join(aggReData)

fcAggData <- fcAgData %>% left_join(aggReData)

aggPredSummary <- data.frame(varName=names(fcAggData)[-(1:7)],
                             shortVarName=names(fcAggData)[-(1:7)])
aggPredSummary <- aggPredSummary %>%
    separate(varName, c("tp", "pred"), "_") %>%
    mutate(tp=substr(tp, 3, 3))

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
outputFolder <- "~/Projects/VRC332/Data/PredictSetpointVL/AggregatePredictors"
for ( tp in 0:8 ) {
  cat("    Starting Timepoint", tp, "\n")
  tpData <- GetTimepointData(fcAggData, tp, aggPredSummary, TRUE, "LogSetpointVL")
  outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  predResults <- RunRandomPartitionPredictions(
                   tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-CompRmse.txt", sep=""),
                   fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)

  bestPredData <- GetBestPredictorsData(fcAggData, predResults, nRand/2,
                                        TRUE, "LogSetpointVL", "tp")
  predResults <- RunRandomPartitionPredictions(
                   bestPredData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-CompRmse.txt", sep=""),
                   fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)
}

tpLevels <- c("AllPreChallenge", "AllPreTerminal", "All")
tpSets <- list(apc=0:6, apt=0:7, all=0:8)

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
for ( i in 1:length(tpLevels) ) {
  cat("    Starting", tpLevels[i], "\n")
  tpData <- GetTimepointData(fcAggData, tpSets[[i]], aggPredSummary, TRUE, "LogSetpointVL")
  outdir <- file.path("~/Projects/VRC332/Data/PredictSetpointVL/AggregatePredictors", tpLevels[i])
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  predResults <- RunRandomPartitionPredictions(
                   tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste(tpLevels[i], "-CompRmse.txt", sep=""),
                   fitsFile=paste(tpLevels[i], "-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)

  bestPredData <- GetBestPredictorsData(fcAggData, predResults, nRand/2,
                                        TRUE, "LogSetpointVL", "tp")
  predResults <- RunRandomPartitionPredictions(
                   bestPredData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste(tpLevels[i], "-BestPred-CompRmse.txt", sep=""),
                   fitsFile=paste(tpLevels[i], "-BestPred-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)
}
