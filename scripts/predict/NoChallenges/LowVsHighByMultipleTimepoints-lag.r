require(caret)
require(glmnet)
require(reshape2)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetBestPredictorsData.r"))
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions_Classification.r"))
source(file.path(fnFolder, "GetDeltaData.r"))

deltaData <- GetDeltaData(fcData.na, predSummary)

tpLevels <- c("AllPreChallenge", "AllPreTerminal", "All")
tpSets <- list(apc=0:6, apt=0:7, all=0:8)

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
outputFolder <- "~/Projects/VRC332/Data/PredictNoChallenges/LowVsHigh/lag"
maxLow <- 6
for ( i in 1:length(tpLevels) ) {
  cat("    Starting ", tpLevels[i], "\n")
  tpData <- GetTimepointData(deltaData, tpSets[[i]], predSummary, TRUE, "NoChallenges")
  tpData$isHigh <- rep(FALSE, length(tpData$y))
  tpData$isHigh[tpData$y > maxLow] <- TRUE
  tpData$y <- tpData$isHigh
  outdir <- file.path(outputFolder, tpLevels[i])
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  predResults <- RunRandomPartitionPredictions_Classification(
      tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
      outdir=outdir,
      resultsFile=paste(tpLevels[i], "-Results.rdata", sep=""),
      verbose=TRUE, progressEvery=10)
}

for ( i in 1:length(tpLevels) ) {
    outdir <- file.path(outputFolder, tpLevels[i])
    predResults <- GetRandomPartitionPredictions_Classification(
        outdir,
        resultsFile=paste(tpLevels[i], "-Results.rdata", sep=""))
    cat("============================================================\n")
    cat(tpLevels[i], " Data:\n", sep="")
    cat("Training Accuracy: ", mean(predResults$accuracy$train), "\n", sep="")
    cat("Testing Accuracy: ", mean(predResults$accuracy$test), "\n", sep="")
    cat("Number of true predictions: ", sum(as.logical(predResults$testResults$pred)), "\n", sep="")
    cat("Number of false predictions: ", sum(!as.logical(predResults$testResults$pred)), "\n", sep="")
}
