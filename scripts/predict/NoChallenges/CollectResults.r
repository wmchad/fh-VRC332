fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetTimepointData.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))
source(file.path(fnFolder, "SummarizePredictionResults.r"))

resultsTable <- data.frame(timepoint=character(0), subset=logical(0),
                           nAnimals=numeric(0), nPred=numeric(0),
                           baselineRmse=numeric(0), groupRmse=numeric(0),
                           predictionRmse=numeric(0),
                           stringsAsFactors=FALSE)

outputFolder <- "~/Projects/VRC332/Data/PredictNoChallenges"
for ( tp in 0:8 ) {
    tpData <- GetTimepointData(fcData, tp, predSummary, TRUE, "NoChallenges")
    nAnimals <- nrow(tpData$x)
    nPred <- ncol(tpData$x) - 4
    outdir <- file.path(outputFolder, paste("tp", tp, "-", tpLabels[tp+1], sep=""))
    predResults <- GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-CompRmse.txt", sep=""),
        fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-Fits.rdata", sep=""))
    resSum <- SummarizePredictionResults(predResults, FALSE)

    resultsTable[nrow(resultsTable)+1,] <- c(tpLabels[tp+1], FALSE, nAnimals,
                                             nPred, resSum$rmses)

    predResults2 <- GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-CompRmse.txt", sep=""),
        fitsFile=paste("tp", tp, "-", tpLabels[tp+1], "-BestPred-Fits.rdata", sep=""))
    resSum2 <- SummarizePredictionResults(predResults2, FALSE)
    resultsTable[nrow(resultsTable)+1,] <- c(tpLabels[tp+1], TRUE, nAnimals,
                                             resSum$nPicked, resSum2$rmses)
}

tpLevels <- c("AllPreChallenge", "AllPreTerminal", "All")
tpSets <- list(apc=0:6, apt=0:7, all=0:8)

for ( i in 1:length(tpLevels) ) {
    tpData <- GetTimepointData(fcData, tpSets[[i]], predSummary, TRUE, "NoChallenges")
    nAnimals <- nrow(tpData$x)
    nPred <- ncol(tpData$x) - 4
    outdir <- file.path("~/Projects/VRC332/Data/PredictNoChallenges", tpLevels[i])
    predResults <- GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste(tpLevels[i], "-CompRmse.txt", sep=""),
        fitsFile=paste(tpLevels[i], "-Fits.rdata", sep=""))
    resSum <- SummarizePredictionResults(predResults, FALSE)
    resultsTable[nrow(resultsTable)+1,] <- c(tpLevels[i], FALSE, nAnimals,
                                             nPred, resSum$rmses)

    predResults2 <- GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste(tpLevels[i], "-BestPred-CompRmse.txt", sep=""),
        fitsFile=paste(tpLevels[i], "-BestPred-Fits.rdata", sep=""))
    resSum2 <- SummarizePredictionResults(predResults2, FALSE)
    resultsTable[nrow(resultsTable)+1,] <- c(tpLevels[i], TRUE, nAnimals,
                                             resSum$nPicked, resSum2$rmses)
}

setwd("~/Projects/VRC332/Reports/201610-NoChallengesPrediction")
write.table(resultsTable, "TimepointAnalysisResults.txt", quote=F, row.names=F)
