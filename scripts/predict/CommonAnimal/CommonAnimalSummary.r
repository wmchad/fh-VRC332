require(glmnet)
fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))
source(file.path(fnFolder, "SummarizePredictionResults.r"))

dataNames <- c("Original", "Adjusted", "Lag")
tps <- 0:8
tpLevels <- c("AllPreChallenge", "AllPreTerminal", "All")
tpSets <- list(apc=0:6, apt=0:7, all=0:8)

dataFolder <- "~/Projects/VRC332/Data"
noChallengesFolder <- "PredictNoChallenges/CommonAnimals"
peakVLFolder <- "PredictPeakVL/CommonAnimals"
setpointVLFolder <- "PredictSetpointVL/CommonAnimals"

PctImprove <- function(rmses) {
    c((rmses[1]-rmses[3])/rmses[1], (rmses[2]-rmses[3])/rmses[2])
}

PickedVariables <- function(summary, predData) {
    pickedDf <- as.data.frame(summary$pickedVars)
    names(pickedDf) <- c("shortVarName", "nPicked")
    predData %>% filter(shortVarName %in% pickedDf$shortVarName[pickedDf$nPicked >= 25]) %>%
        left_join(pickedDf) %>% arrange(desc(nPicked))
}

CollectResults <- function(target, outputFolder, predData=predSummary, dataNms=dataNames,
                           tpset=tps, tpLabs=tpLabels, tpLevs=tpLevels) {
    
    rmses <- data.frame(target=character(0), data=character(0),
                        predictors=character(0),
                        timepoints=character(0), nvar=numeric(0),
                        averageRmse=numeric(0), groupRmse=numeric(0),
                        predictedRmse=numeric(0), pctImproveAverage=numeric(0),
                        pctImproveGroup=numeric(0), stringsAsFactors=FALSE)
    for ( dataName in dataNms ) {
        for ( tp in tpset ) {
            outdir <- file.path(outputFolder, dataName,
                                paste("tp", tp, "-", tpLabs[tp+1], sep=""))
            if (file.exists(outdir)) {
                resSummary <- SummarizePredictionResults(
                    GetRandomPartitionPredictions(
                        outdir=outdir,
                        rmseFile=paste("tp", tp, "-", tpLabs[tp+1], "-CompRmse.txt", sep=""),
                        fitsFile=paste("tp", tp, "-", tpLabs[tp+1], "-Fits.rdata", sep="")),
                    FALSE)
                rmses[nrow(rmses)+1,] <- c(target, dataName, "All", tpLabs[tp+1],
                                           resSummary$nvar, resSummary$rmses,
                                           PctImprove(resSummary$rmses))
                write.table(PickedVariables(resSummary, predData),
                            file.path(outdir, paste("tp", tp, "-", tpLabs[tp+1],
                                                    "-PickedVariables.txt", sep="")),
                            quote=FALSE, row.names=FALSE)

                resSummary <- SummarizePredictionResults(
                    GetRandomPartitionPredictions(
                        outdir=outdir,
                        rmseFile=paste("tp", tp, "-", tpLabs[tp+1],
                                       "-BestPred-CompRmse.txt", sep=""),
                        fitsFile=paste("tp", tp, "-", tpLabs[tp+1],
                                       "-BestPred-Fits.rdata", sep="")),
                    FALSE)
                rmses[nrow(rmses)+1,] <- c(target, dataName, "Best", tpLabs[tp+1],
                                              resSummary$nvar, resSummary$rmses,
                                              PctImprove(resSummary$rmses))
            }
        }

        for ( i in 1:length(tpLevs) ) {
            outdir <- file.path(outputFolder, dataName, tpLevs[i])
            resSummary <- SummarizePredictionResults(
                GetRandomPartitionPredictions(
                    outdir=outdir,
                    rmseFile=paste(tpLevs[i], "-CompRmse.txt", sep=""),
                    fitsFile=paste(tpLevs[i], "-Fits.rdata", sep="")),
                FALSE)
            rmses[nrow(rmses)+1,] <- c(target, dataName, "All", tpLevs[i],
                                       resSummary$nvar, resSummary$rmses,
                                       PctImprove(resSummary$rmses))
            write.table(PickedVariables(resSummary, predData),
                        file.path(outdir, paste(tpLevs[i], "-PickedVariables.txt", sep="")),
                        quote=FALSE, row.names=FALSE)


            resSummary <- SummarizePredictionResults(
                GetRandomPartitionPredictions(
                    outdir=outdir,
                    rmseFile=paste(tpLevs[i], "-BestPred-CompRmse.txt", sep=""),
                    fitsFile=paste(tpLevs[i], "-BestPred-Fits.rdata", sep="")),
                FALSE)
            rmses[nrow(rmses)+1,] <- c(target, dataName, "Best", tpLevs[i],
                                       resSummary$nvar, resSummary$rmses,
                                       PctImprove(resSummary$rmses))
        }
    }
    rmses
}

rmses.noCh <- CollectResults("NoChallenges", file.path(dataFolder, noChallengesFolder))
rmses.peakVL <- CollectResults("PeakVL", file.path(dataFolder, peakVLFolder))
rmses.spVL <- CollectResults("SetpointVL", file.path(dataFolder, setpointVLFolder))

for ( i in 5:10 ) {
    rmses.noCh[,i] <- as.numeric(rmses.noCh[,i])
    rmses.peakVL[,i] <- as.numeric(rmses.peakVL[,i])
    rmses.spVL[,i] <- as.numeric(rmses.spVL[,i])
}

write.table(rmses.noCh, file.path(dataFolder, noChallengesFolder, "rmseSummary.txt"),
            quote=FALSE, row.names=FALSE)
write.table(rmses.peakVL, file.path(dataFolder, peakVLFolder, "rmseSummary.txt"),
            quote=FALSE, row.names=FALSE)
write.table(rmses.spVL, file.path(dataFolder, setpointVLFolder, "rmseSummary.txt"),
            quote=FALSE, row.names=FALSE)


GetSelectedVars <- function(dataFolder, predictionFolder,
                            dataName, tpSetName) {
    read.table(file.path(dataFolder, predictionFolder, dataName, tpSetName,
                         paste(tpSetName, "-PickedVariables.txt", sep="")),
               header=T, stringsAsFactors=FALSE)
}

tbl1 <- GetSelectedVars(dataFolder, noChallengesFolder, "Lag", "AllPreChallenge") ## 29.7%
tbl2 <- GetSelectedVars(dataFolder, noChallengesFolder, "Adjusted", "AllPreChallenge") ## 28.3%
tbl3 <- GetSelectedVars(dataFolder, noChallengesFolder, "Original", "AllPreChallenge") ## 23.3%

## Common between Lag and Adjusted
## tp            re                 ag shortVarName nPicked
##  1 aRhIgG.PE.low SIVcpz.EK505.gp120       var204      45
##  1         R2A.3    SIVmac239.gp120       var263      35
##  1         R2A.3    SIVsmH4.p55.Gag       var270      35
## In Original
## tp            re                      ag shortVarName nPicked
##  0 aRhIgG.PE.low      SIVcpz.EK505.gp120        var32      42

## Common between Original and Adjusted
## tp       re                ag shortVarName nPicked
##  2    R2A.2   SIVmac239.gp140       var422      36
##  2    R3A.3    SIV.1A11.gp140       var505      33

## Similar
## Lag:
## tp            re                  ag shortVarName nPicked
##  3          R3A.3  SIVmac251.BK.PR55       var684      28
##  4     R2A.4.high     SIV.1A11.gp140       var797      30
## Adjusted:
## tp            re                 ag shortVarName nPicked
##  2         R3A.3  SIVmac251.BK.PR55       var512      37
##  3    R2A.4.high     SIV.1A11.gp140       var625      31

tbl1 <- GetSelectedVars(dataFolder, peakVLFolder, "Lag", "AllPreChallenge") ## 12.8%
tbl2 <- GetSelectedVars(dataFolder, peakVLFolder, "Adjusted", "AllPreChallenge") ## 14.8%
tbl3 <- GetSelectedVars(dataFolder, peakVLFolder, "Original", "AllPreChallenge") ## 10.3%

## Common between Lag and Original
## tp      re                ag shortVarName nPicked
##  5   R2A.3               G49       var943      38

## Common between Lag and Adjusted
## tp            re                      ag shortVarName nPicked
##  1 aRhIgG.PE.low J08.V1V2.mac239.AVI.His       var201      26

## Common between Adjusted and Original
## tp         re              ag shortVarName nPicked
##  2 R2A.4.high SIVsmH4.p55.Gag       var463      29
## Lag:
## tp         re              ag shortVarName nPicked
## 1  R2A.4.high SIVsmH4.p55.Gag       var291      28

## Similar:
## Lag:
## tp         re              ag shortVarName nPicked
##  5        MBL             G49       var926      25
##  3      R2A.3           C1.Ak       var595      25
## Adjusted:
## tp         re              ag shortVarName nPicked
##  4        MBL             G49       var754      39
## Original:
## tp         re              ag shortVarName nPicked
##  4      R2A.3           C1.Ak       var767      46

tbl1 <- GetSelectedVars(dataFolder, setpointVLFolder, "Original", "tp0-Baseline") ## 27.3%
tbl2 <- GetSelectedVars(dataFolder, setpointVLFolder, "Original", "AllPreChallenge") ## 26.2%
tbl3 <- GetSelectedVars(dataFolder, setpointVLFolder, "Lag", "tp2-2xDna") ## 15.3%

