require(caret)
require(glmnet)
require(reshape2)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetBestPredictorsData.r"))
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions_Classification.r"))
source(file.path(fnFolder, "GetBaselineAdjustedData.r"))

adjData <- GetBaselineAdjustedData(fcData.na, predSummary)

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
outputFolder <- "~/Projects/VRC332/Data/PredictNoChallenges/LowVsHigh/adj"
maxLow <- 6
for ( tp in 7:8 ) {
  cat("    Starting Timepoint", tp, "\n")
  tpData <- GetTimepointData(adjData, tp, predSummary, TRUE, "NoChallenges")
  tpData$x <- ReplaceMissingWithMean(tpData$x)
  tpData$x <- tpData$x[,!apply(tpData$x, 2,
                              function(x){sum(is.nan(x))}) == nrow(tpData$x)]
  tpData$isHigh <- rep(FALSE, length(tpData$y))
  tpData$isHigh[tpData$y > maxLow] <- TRUE
  tpData$y <- tpData$isHigh
  outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  predResults <- RunRandomPartitionPredictions_Classification(
      tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
      outdir=outdir,
      resultsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Results.txt", sep=""),
      verbose=TRUE, progressEvery=10)
}

for ( tp in 1:8 ) {
    outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
    predResults <- GetRandomPartitionPredictions_Classification(
        outdir,
        resultsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Results.txt", sep=""))
    cat("============================================================\n")
    cat("Timepoint ", tp, " - ", tpLabels[tp+1], ":\n", sep="")
    cat("Training Accuracy: ", mean(predResults$accuracy$train), "\n", sep="")
    cat("Testing Accuracy: ", mean(predResults$accuracy$test), "\n", sep="")
    cat("Number of true predictions: ", sum(as.logical(predResults$testResults$pred)), "\n", sep="")
    cat("Number of false predictions: ", sum(!as.logical(predResults$testResults$pred)), "\n", sep="")
}

## ============================================================
## Timepoint 1 - 1xDna:
## Training Accuracy: 0.9471186
## Testing Accuracy: 0.8
## Number of true predictions: 33
## Number of false predictions: 717
## ============================================================
## Timepoint 2 - 2xDna:
## Training Accuracy: 0.9423729
## Testing Accuracy: 0.8106667
## Number of true predictions: 60
## Number of false predictions: 690
## ============================================================
## Timepoint 3 - PeakPostPrime:
## Training Accuracy: 0.8342373
## Testing Accuracy: 0.784
## Number of true predictions: 18
## Number of false predictions: 732
## ============================================================
## Timepoint 4 - PreBoost:
## Training Accuracy: 0.860678
## Testing Accuracy: 0.748
## Number of true predictions: 37
## Number of false predictions: 713
## ============================================================
## Timepoint 5 - PeakPostBoost:
## Training Accuracy: 0.8722034
## Testing Accuracy: 0.7706667
## Number of true predictions: 53
## Number of false predictions: 697
## ============================================================
## Timepoint 6 - PreChallenge:
## Training Accuracy: 0.9345763
## Testing Accuracy: 0.7786667
## Number of true predictions: 83
## Number of false predictions: 667
## ============================================================
## Timepoint 7 - PostChallenge:
## Training Accuracy: 0.9264407
## Testing Accuracy: 0.812
## Number of true predictions: 62
## Number of false predictions: 688
## ============================================================
## Timepoint 8 - Terminal:
## Training Accuracy: 0.8637288
## Testing Accuracy: 0.8173333
## Number of true predictions: 45
## Number of false predictions: 705
