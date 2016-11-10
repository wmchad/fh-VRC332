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

aggData <- AggregateData(fcData, predSummary)
fcAggData <- aggData$aggData
aggPredSummary <- aggData$predSummary

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
outputFolder <- "~/Projects/VRC332/Data/PredictPeakVL/AggregatePredictors"
for ( tp in 0:8 ) {
  cat("    Starting Timepoint", tp, "\n")
  tpData <- GetTimepointData(fcAggData, tp, aggPredSummary, TRUE, "LogPeakVL")
  outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  predResults <- RunRandomPartitionPredictions(
                   tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-CompRmse.txt", sep=""),
                   fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)

  bestPredData <- GetBestPredictorsData(fcAggData, predResults, nRand/2,
                                        TRUE, "LogPeakVL", "tp")
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
  tpData <- GetTimepointData(fcAggData, tpSets[[i]], aggPredSummary, TRUE, "LogPeakVL")
  outdir <- file.path("~/Projects/VRC332/Data/PredictPeakVL/AggregatePredictors", tpLevels[i])
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  predResults <- RunRandomPartitionPredictions(
                   tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste(tpLevels[i], "-CompRmse.txt", sep=""),
                   fitsFile=paste(tpLevels[i], "-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)

  bestPredData <- GetBestPredictorsData(fcAggData, predResults, nRand/2,
                                        TRUE, "LogPeakVL", "tp")
  predResults <- RunRandomPartitionPredictions(
                   bestPredData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste(tpLevels[i], "-BestPred-CompRmse.txt", sep=""),
                   fitsFile=paste(tpLevels[i], "-BestPred-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)
}


outputFolder <- "~/Projects/VRC332/Data/PredictPeakVL/AggregatePredictors"
for ( tp in 0:8 ) {
    print("==============================")
    print(paste("tp ", tp, ": ", tpLabels[tp+1], sep=""))
    print("==============================")
    outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
    resSummary <- GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-CompRmse.txt", sep=""),
        fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Fits.rdata", sep="")) %>%
        SummarizePredictionResults()
    resSummary <- GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-CompRmse.txt", sep=""),
        fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-Fits.rdata", sep="")) %>%
        SummarizePredictionResults()
}

for ( i in 1:length(tpLevels) ) {
    print("==============================")
    print(tpLevels[i])
    print("==============================")
    outdir <- file.path("~/Projects/VRC332/Data/PredictPeakVL/AggregatePredictors", tpLevels[i])
    resSummary <- GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste(tpLevels[i], "-CompRmse.txt", sep=""),
        fitsFile=paste(tpLevels[i], "-Fits.rdata", sep="")) %>%
        SummarizePredictionResults()
    resSummary <- GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste(tpLevels[i], "-BestPred-CompRmse.txt", sep=""),
        fitsFile=paste(tpLevels[i], "-BestPred-Fits.rdata", sep="")) %>%
        SummarizePredictionResults()
}
