require(caret)
require(glmnet)
require(reshape2)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetBestPredictorsData.r"))
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))
source(file.path(fnFolder, "GetBaselineAdjustedData.r"))
source(file.path(fnFolder, "GetDeltaData.r"))

origData <- fcData.na[apply(fcData.na[,-(1:7)], 1, function(x){!any(x<0)}),]
adjData <- GetBaselineAdjustedData(origData, predSummary)
lagData <- GetDeltaData(origData, predSummary)

tpLevels <- c("AllPreChallenge", "AllPreTerminal", "All")
tpSets <- list(apc=0:6, apt=0:7, all=0:8)

nRand <- 50
pTrain <- 0.75
ctrl <- trainControl(method="repeatedcv", repeats=20)
dataFolder <- "~/Projects/VRC332/Data"
noChallengesFolder <- "PredictNoChallenges/CommonAnimals"
peakVLFolder <- "PredictPeakVL/CommonAnimals"
setpointVLFolder <- "PredictSetpointVL/CommonAnimals"

RunPrediction <- function(data, target, outputFolder, tps=0:8) {
    for ( tp in tps ) {
        cat("    Starting Timepoint", tp, "\n")
        tpData <- GetTimepointData(data, tp, predSummary, TRUE, target)
        outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
        dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
        predResults <- RunRandomPartitionPredictions(
            tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
            outdir=outdir,
            rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-CompRmse.txt", sep=""),
            fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Fits.rdata", sep=""),
            verbose=TRUE, progressEvery=10)

        bestPredData <- GetBestPredictorsData(data, predResults, nRand/2, TRUE, target)
        predResults <- RunRandomPartitionPredictions(
            bestPredData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
            outdir=outdir,
            rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-CompRmse.txt", sep=""),
            fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-Fits.rdata", sep=""),
            verbose=TRUE, progressEvery=10)
    }

    for ( i in 1:length(tpLevels) ) {
        cat("    Starting", tpLevels[i], "\n")
        tpData <- GetTimepointData(data, tpSets[[i]], predSummary, TRUE, target)
        outdir <- file.path(outputFolder, tpLevels[i])
        dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
        predResults <- RunRandomPartitionPredictions(
            tpData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
            outdir=outdir,
            rmseFile=paste(tpLevels[i], "-CompRmse.txt", sep=""),
            fitsFile=paste(tpLevels[i], "-Fits.rdata", sep=""),
            verbose=TRUE, progressEvery=10)

        bestPredData <- GetBestPredictorsData(data, predResults, nRand/2, TRUE, target)
        predResults <- RunRandomPartitionPredictions(
            bestPredData, pTrain=pTrain, nRand=nRand, method="glmnet", ctrl=ctrl,
            outdir=outdir,
            rmseFile=paste(tpLevels[i], "-BestPred-CompRmse.txt", sep=""),
            fitsFile=paste(tpLevels[i], "-BestPred-Fits.rdata", sep=""),
            verbose=TRUE, progressEvery=10)
    }
}

cat("============================================================\n")
cat("Running # Challenges, Original Data\n")
cat("============================================================\n")
RunPrediction(origData, "NoChallenges",
              file.path(dataFolder, noChallengesFolder, "Original"))
cat("============================================================\n")
cat("Running # Challenges, Adjusted Data\n")
cat("============================================================\n")
RunPrediction(adjData, "NoChallenges",
              file.path(dataFolder, noChallengesFolder, "Adjusted"), 1:8)
cat("============================================================\n")
cat("Running # Challenges, Lag Data\n")
cat("============================================================\n")
RunPrediction(lagData, "NoChallenges",
              file.path(dataFolder, noChallengesFolder, "Lag"), 1:8)

cat("============================================================\n")
cat("Running Peak VL, Original Data\n")
cat("============================================================\n")
RunPrediction(origData, "LogPeakVL",
              file.path(dataFolder, peakVLFolder, "Original"))
cat("============================================================\n")
cat("Running Peak VL, Adjusted Data\n")
cat("============================================================\n")
RunPrediction(adjData, "LogPeakVL",
              file.path(dataFolder, peakVLFolder, "Adjusted"), 1:8)
cat("============================================================\n")
cat("Running Peak VL, Lag Data\n")
cat("============================================================\n")
RunPrediction(lagData, "LogPeakVL",
              file.path(dataFolder, peakVLFolder, "Lag"), 1:8)

cat("============================================================\n")
cat("Running Setpoint VL, Original Data\n")
cat("============================================================\n")
RunPrediction(origData, "LogSetpointVL",
              file.path(dataFolder, setpointVLFolder, "Original"))
cat("============================================================\n")
cat("Running Setpoint VL, Adjusted Data\n")
cat("============================================================\n")
RunPrediction(adjData, "LogSetpointVL",
              file.path(dataFolder, setpointVLFolder, "Adjusted"), 1:8)
cat("============================================================\n")
cat("Running Setpoint VL, Lag Data\n")
cat("============================================================\n")
RunPrediction(lagData, "LogSetpointVL",
              file.path(dataFolder, setpointVLFolder, "Lag"), 1:8)

