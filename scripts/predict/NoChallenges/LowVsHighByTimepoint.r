require(caret)
require(glmnet)
require(reshape2)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetBestPredictorsData.r"))
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions_Classification.r"))

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
outputFolder <- "~/Projects/VRC332/Data/PredictNoChallenges/LowVsHigh"
maxLow <- 6
for ( tp in 0:8 ) {
  cat("    Starting Timepoint", tp, "\n")
  tpData <- GetTimepointData(fcData, tp, predSummary, TRUE, "NoChallenges")
  tpData$x <- ReplaceMissingWithMean(tpData$x)
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

for ( tp in 0:8 ) {
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
## Timepoint 0 - Baseline:
## Training Accuracy: 0.899322
## Testing Accuracy: 0.788
## Number of true predictions: 46
## Number of false predictions: 704
## ============================================================
## Timepoint 1 - 1xDna:
## Training Accuracy: 0.8589831
## Testing Accuracy: 0.76375
## Number of true predictions: 43
## Number of false predictions: 757
## ============================================================
## Timepoint 2 - 2xDna:
## Training Accuracy: 0.8535593
## Testing Accuracy: 0.775
## Number of true predictions: 41
## Number of false predictions: 759
## ============================================================
## Timepoint 3 - PeakPostPrime:
## Training Accuracy: 0.9027119
## Testing Accuracy: 0.77875
## Number of true predictions: 51
## Number of false predictions: 749
## ============================================================
## Timepoint 4 - PreBoost:
## Training Accuracy: 0.8755932
## Testing Accuracy: 0.75625
## Number of true predictions: 51
## Number of false predictions: 749
## ============================================================
## Timepoint 5 - PeakPostBoost:
## Training Accuracy: 0.8429333
## Testing Accuracy: 0.7209091
## Number of true predictions: 127
## Number of false predictions: 973
## ============================================================
## Timepoint 6 - PreChallenge:
## Training Accuracy: 0.8277966
## Testing Accuracy: 0.78875
## Number of true predictions: 21
## Number of false predictions: 779
## ============================================================
## Timepoint 7 - PostChallenge:
## Training Accuracy: 0.9210667
## Testing Accuracy: 0.7918182
## Number of true predictions: 168
## Number of false predictions: 932
## ============================================================
## Timepoint 8 - Terminal:
## Training Accuracy: 0.8925333
## Testing Accuracy: 0.7572727
## Number of true predictions: 133
## Number of false predictions: 967
