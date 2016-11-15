require(caret)
require(glmnet)
require(reshape2)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetBestPredictorsData.r"))
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
outputFolder <- "~/Projects/VRC332/Data/PredictNoChallenges"
for ( tp in 0:8 ) {
  cat("    Starting Timepoint", tp, "\n")
  tpData <- GetTimepointData(fcData, tp, predSummary, TRUE, "NoChallenges")
  tpData$x <- ReplaceMissingWithMean(tpData$x)
  outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  predResults <- RunRandomPartitionPredictions(
                   tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-CompRmse.txt", sep=""),
                   fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)

  bestPredData <- GetBestPredictorsData(fcData, predResults, nRand/2, TRUE, "NoChallenges")
  bestPredData$x <- ReplaceMissingWithMean(bestPredData$x)
  predResults <- RunRandomPartitionPredictions(
                   bestPredData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
                   outdir=outdir,
                   rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-CompRmse.txt", sep=""),
                   fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-Fits.rdata", sep=""),
                   verbose=TRUE, progressEvery=10)
}
