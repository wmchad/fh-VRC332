require(caret)
require(glmnet)
require(glmnetUtils)
require(dplyr)
require(survival)
require(reshape2)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "GetBestPredictorsData.r"))
source(file.path(fnFolder, "GetVariableSetData.r"))
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions_CoxPh.r"))
source(file.path(fnFolder, "GetBaselineAdjustedData.r"))
source(file.path(fnFolder, "GetDeltaData.r"))
source(file.path(fnFolder, "Prediction/CoxPhHelper.r"))

origData <- fcData.na[apply(fcData.na[,-(1:7)], 1, function(x){!any(x<0)}),]
adjData <- GetBaselineAdjustedData(origData, predSummary)
lagData <- GetDeltaData(origData, predSummary)
origData[,-(1:7)] <- scale(origData[,-(1:7)])
adjData[,-(1:7)] <- scale(adjData[,-(1:7)])
lagData[,-(1:7)] <- scale(lagData[,-(1:7)])

tpLevels <- c("AllPreChallenge", "AllPreTerminal", "All")
tpSets <- list(apc=0:6, apt=0:7, all=0:8)

nRand <- 200
pTrain <- 0.75
coxIter <- 100
dataFolder <- "~/Projects/VRC332/Data"
noChallengesFolder <- "PredictNoChallenges/CommonAnimals/CoxPh"

BuildSurv <- function(y) {
    event <- rep(1, length(y))
    event[y==13] <- 0
    Surv(y, event)
}

RunPrediction_CoxPh<- function(data, target, outputFolder, tps=0:8, tpLabs=tpLabels,
                               tpLevs=tpLevels, nRand=200, pTrain=0.75,
                               coxIter=100) {
    for ( tp in tps ) {
        cat("    Starting Timepoint", tp, "\n")
        tpData <- GetTimepointData(data, tp, predSummary, TRUE, target)
        tpData$surv <- BuildSurv(tpData$y)
        outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabs[tp+1], sep=""))
        dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
        predResults <- RunRandomPartitionPredictions_CoxPh(
            tpData, pTrain=pTrain, nRand=nRand, coxIter=coxIter,
            outdir=outdir,
            rmseFile=paste("tp", tp, "-", tpLabs[tp+1], "-CompRmse.txt", sep=""),
            fitsFile=paste("tp", tp, "-", tpLabs[tp+1], "-Fits.rdata", sep=""),
            verbose=TRUE, progressEvery=50)
        z <- apply(predResults$compRmse, 2, mean)
        cat("Full data improvement:\n")
        cat("  cox:        ", round(100*(z[2]-z[3])/z[2], 2), "%\n", sep="")
        cat("  cox shrunk: ", round(100*(z[2]-z[4])/z[2], 2), "%\n", sep="")

        bestPredData <- GetBestPredictorsData_CoxPh(data, predResults, nRand/2, TRUE, target)
        if ( ncol(bestPredData$x) > 4 ) {
            bestPredData$surv <- BuildSurv(bestPredData$y)
            predResults <- RunRandomPartitionPredictions_CoxPh(
                bestPredData, pTrain=pTrain, nRand=nRand, coxIter=coxIter,
                outdir=outdir,
                rmseFile=paste("tp", tp, "-", tpLabs[tp+1], "-BestPred-CompRmse.txt", sep=""),
                fitsFile=paste("tp", tp, "-", tpLabs[tp+1], "-BestPred-Fits.rdata", sep=""),
                verbose=TRUE, progressEvery=50)
            z <- apply(predResults$compRmse, 2, mean)
            cat("Best data improvement:\n")
            cat("  cox:        ", round(100*(z[2]-z[3])/z[2], 2), "%\n", sep="")
            cat("  cox shrunk: ", round(100*(z[2]-z[4])/z[2], 2), "%\n", sep="")
        }
    }

    for ( i in 1:length(tpLevs) ) {
        cat("    Starting", tpLevs[i], "\n")
        tpData <- GetTimepointData(data, tpSets[[i]], predSummary, TRUE, target)
        tpData$surv <- BuildSurv(tpData$y)
        outdir <- file.path(outputFolder, tpLevs[i])
        dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
        predResults <- RunRandomPartitionPredictions_CoxPh(
            tpData, pTrain=pTrain, nRand=nRand, coxIter=coxIter,
            outdir=outdir,
            rmseFile=paste(tpLevs[i], "-CompRmse.txt", sep=""),
            fitsFile=paste(tpLevs[i], "-Fits.rdata", sep=""),
            verbose=TRUE, progressEvery=50)
        z <- apply(predResults$compRmse, 2, mean)
        cat("Full data improvement:\n")
        cat("  cox:        ", round(100*(z[2]-z[3])/z[2], 2), "%\n", sep="")
        cat("  cox shrunk: ", round(100*(z[2]-z[4])/z[2], 2), "%\n", sep="")

        bestPredData <- GetBestPredictorsData_CoxPh(data, predResults, nRand/2, TRUE, target)
        if ( ncol(bestPredData$x) > 4 ) {
            bestPredData$surv <- BuildSurv(bestPredData$y)
            predResults <- RunRandomPartitionPredictions_CoxPh(
                bestPredData, pTrain=pTrain, nRand=nRand, coxIter=coxIter,
                outdir=outdir,
                rmseFile=paste(tpLevs[i], "-BestPred-CompRmse.txt", sep=""),
                fitsFile=paste(tpLevs[i], "-BestPred-Fits.rdata", sep=""),
                verbose=TRUE, progressEvery=50)
            z <- apply(predResults$compRmse, 2, mean)
            cat("Full data improvement:\n")
            cat("  cox:        ", round(100*(z[2]-z[3])/z[2], 2), "%\n", sep="")
            cat("  cox shrunk: ", round(100*(z[2]-z[4])/z[2], 2), "%\n", sep="")
        }
    }
}

cat("============================================================\n")
cat("Running # Challenges, Original Data\n")
cat("============================================================\n")
RunPrediction_CoxPh(origData, "NoChallenges",
                    file.path(dataFolder, noChallengesFolder, "Original"), 3:8)
cat("============================================================\n")
cat("Running # Challenges, Adjusted Data\n")
cat("============================================================\n")
RunPrediction_CoxPh(adjData, "NoChallenges",
                    file.path(dataFolder, noChallengesFolder, "Adjusted"), NULL)
cat("============================================================\n")
cat("Running # Challenges, Lag Data\n")
cat("============================================================\n")
RunPrediction_CoxPh(lagData, "NoChallenges",
                    file.path(dataFolder, noChallengesFolder, "Lag"), 3:8)





