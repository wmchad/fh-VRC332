require(caret)
require(glmnet)
require(reshape2)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetBestPredictorsData.r"))
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))

tpLevels <- c("AllPreChallenge", "AllPreTerminal", "All")
tpSets <- list(apc=0:6, apt=0:7, all=0:8)

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
for ( i in 1:length(tpLevels) ) {
  cat("    Starting", tpLevels[i], "\n")
  tpData <- GetTimepointData(fcData, tpSets[[i]], predSummary, TRUE, "NoChallenges")
  tpData$x <- ReplaceMissingWithMean(tpData$x)
  outdir <- file.path("~/Projects/VRC332/Data/PredictNoChallenges", tpLevels[i])
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  predResults <- RunRandomPartitionPredictions(
                   tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste(tpLevels[i], "-CompRmse.txt", sep=""),
                   fitsFile=paste(tpLevels[i], "-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)

  bestPredData <- GetBestPredictorsData(fcData, predResults, nRand/2, TRUE, "NoChallenges")
  bestPredData$x <- ReplaceMissingWithMean(bestPredData$x)
  predResults <- RunRandomPartitionPredictions(
                   bestPredData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste(tpLevels[i], "-BestPred-CompRmse.txt", sep=""),
                   fitsFile=paste(tpLevels[i], "-BestPred-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)
}
