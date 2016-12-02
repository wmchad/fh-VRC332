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

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
outputFolder <- "~/Projects/VRC332/Data/PredictNoChallenges/LowVsHigh/lag"
maxLow <- 6
for ( tp in 1:8 ) {
  cat("    Starting Timepoint", tp, "\n")
  tpData <- GetTimepointData(deltaData, tp, predSummary, TRUE, "NoChallenges")
  tpData$isHigh <- rep(FALSE, length(tpData$y))
  tpData$isHigh[tpData$y > maxLow] <- TRUE
  tpData$y <- tpData$isHigh
  outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  predResults <- RunRandomPartitionPredictions_Classification(
      tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
      outdir=outdir,
      resultsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Results.rdata", sep=""),
      verbose=TRUE, progressEvery=10)
}

for ( tp in 1:8 ) {
    outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
    predResults <- GetRandomPartitionPredictions_Classification(
        outdir,
        resultsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Results.rdata", sep=""))
    cat("============================================================\n")
    cat("Timepoint ", tp, " - ", tpLabels[tp+1], ":\n", sep="")
    cat("Training Accuracy: ", mean(predResults$accuracy$train), "\n", sep="")
    cat("Testing Accuracy: ", mean(predResults$accuracy$test), "\n", sep="")
    cat("Number of true predictions: ", sum(as.logical(predResults$testResults$pred)), "\n", sep="")
    cat("Number of false predictions: ", sum(!as.logical(predResults$testResults$pred)), "\n", sep="")
}

## ============================================================
## Timepoint 1 - 1xDna:
## Training Accuracy: 0.8835714
## Testing Accuracy: 0.7506667
## Number of true predictions: 42
## Number of false predictions: 708
## ============================================================
## Timepoint 2 - 2xDna:
## Training Accuracy: 0.8646429
## Testing Accuracy: 0.74
## Number of true predictions: 49
## Number of false predictions: 701
## ============================================================
## Timepoint 3 - PeakPostPrime:
## Training Accuracy: 0.8721429
## Testing Accuracy: 0.7573333
## Number of true predictions: 27
## Number of false predictions: 723
## ============================================================
## Timepoint 4 - PreBoost:
## Training Accuracy: 0.8735714
## Testing Accuracy: 0.7426667
## Number of true predictions: 59
## Number of false predictions: 691
## ============================================================
## Timepoint 5 - PeakPostBoost:
## Training Accuracy: 0.81
## Testing Accuracy: 0.7666667
## Number of true predictions: 22
## Number of false predictions: 728
## ============================================================
## Timepoint 6 - PreChallenge:
## Training Accuracy: 0.8053571
## Testing Accuracy: 0.7733333
## Number of true predictions: 18
## Number of false predictions: 732
## ============================================================
## Timepoint 7 - PostChallenge:
## Training Accuracy: 0.8767857
## Testing Accuracy: 0.8186667
## Number of true predictions: 57
## Number of false predictions: 693
## ============================================================
## Timepoint 8 - Terminal:
## Training Accuracy: 0.8782143
## Testing Accuracy: 0.7613333
## Number of true predictions: 44
## Number of false predictions: 706
