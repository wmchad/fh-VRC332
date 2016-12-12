require(glmnet)
require(dplyr)
fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))
source(file.path(fnFolder, "RunRandomPartitionPredictions_CoxPh.r"))
source(file.path(fnFolder, "SummarizePredictionResults.r"))
## source(file.path(fnFolder, "Prediction/GetEffectDirections.r"))

dataNames <- c("Original", "Adjusted", "Lag")
tps <- 0:8
tpLevels <- c("AllPreChallenge", "AllPreTerminal", "All")
tpSets <- list(apc=0:6, apt=0:7, all=0:8)

dataFolder <- "~/Projects/VRC332/Data"
noChallengesFolder <- "PredictNoChallenges/CommonAnimals/CoxPh"

PctImprove_CoxPh <- function(rmses) {
    c((rmses[1]-rmses[3])/rmses[1], (rmses[2]-rmses[3])/rmses[2],
      (rmses[1]-rmses[4])/rmses[1], (rmses[2]-rmses[4])/rmses[2])
}

BuildRmseRow_CoxPh <- function(res) {
    rmses <- apply(res$compRmse, 2, mean)
    pctImp <- PctImprove_CoxPh(rmses)
    c(rmses[1:3], pctImp[1:2], rmses[4], pctImp[3:4])
}

PickedVariables_CoxPh <- function(res, predData, threshold) {
    pickedVars <- NULL
    for ( fit in res$fits ) {
        pickedVars <- c(pickedVars, fit$vars)
    }
    pickedDf <- as.data.frame(table(pickedVars))
    names(pickedDf) <- c("shortVarName", "nPicked")
    predData %>%
        filter(shortVarName %in% pickedDf$shortVarName[pickedDf$nPicked >= threshold]) %>%
        left_join(pickedDf) %>% arrange(desc(nPicked))
}

CollectResults_CoxPh<- function(target, outputFolder, pickThreshold=100,
                                predData=predSummary, dataNms=dataNames,
                                tpset=tps, tpLabs=tpLabels, tpLevs=tpLevels) {
    
    rmses <- data.frame(target=character(0), data=character(0),
                        predictors=character(0),
                        timepoints=character(0), nvar=numeric(0),
                        averageRmse=numeric(0), groupRmse=numeric(0),
                        predictedRmse=numeric(0), pctImproveAverage=numeric(0),
                        pctImproveGroup=numeric(0), predictedRmse.glm=numeric(0),
                        pctImproveAverage.glm=numeric(0), pctImproveGroup.glm=numeric(0),
                        stringsAsFactors=FALSE)
    for ( dataName in dataNms ) {
        for ( tp in tpset ) {
            outdir <- file.path(outputFolder, dataName,
                                paste("tp", tp, "-", tpLabs[tp+1], sep=""))
            if (file.exists(outdir)) {
                res <- GetRandomPartitionPredictions(
                    outdir=outdir,
                    rmseFile=paste("tp", tp, "-", tpLabs[tp+1], "-CompRmse.txt", sep=""),
                    fitsFile=paste("tp", tp, "-", tpLabs[tp+1], "-Fits.rdata", sep=""))
                rmses[nrow(rmses)+1,] <- c(target, dataName, "All", tpLabs[tp+1],
                                           nrow(res$fits[[1]]$fit$glmnet.fit$beta) - 4,
                                           BuildRmseRow_CoxPh(res))
                pickedVars <- PickedVariables_CoxPh(res, predData, pickThreshold)
                write.table(pickedVars,
                            file.path(outdir, paste("tp", tp, "-", tpLabs[tp+1],
                                                    "-PickedVariables.txt", sep="")),
                            quote=FALSE, row.names=FALSE)

                if ( nrow(pickedVars) > 0 ) {
                    res <- GetRandomPartitionPredictions(
                        outdir=outdir,
                        rmseFile=paste("tp", tp, "-", tpLabs[tp+1],
                                       "-BestPred-CompRmse.txt", sep=""),
                        fitsFile=paste("tp", tp, "-", tpLabs[tp+1],
                                       "-BestPred-Fits.rdata", sep=""))
                    rmses[nrow(rmses)+1,] <- c(target, dataName, "Best", tpLabs[tp+1],
                                               nrow(pickedVars), BuildRmseRow_CoxPh(res))
                }
            }
        }

        for ( i in 1:length(tpLevs) ) {
            outdir <- file.path(outputFolder, dataName, tpLevs[i])
            if (file.exists(outdir)) {
                res <- GetRandomPartitionPredictions(
                    outdir=outdir,
                    rmseFile=paste(tpLevs[i], "-CompRmse.txt", sep=""),
                    fitsFile=paste(tpLevs[i], "-Fits.rdata", sep=""))
                rmses[nrow(rmses)+1,] <- c(target, dataName, "All", tpLevs[i],
                                           nrow(res$fits[[1]]$fit$glmnet.fit$beta) - 4,
                                           BuildRmseRow_CoxPh(res))
                pickedVars <- PickedVariables_CoxPh(res, predData, pickThreshold)
                write.table(pickedVars,
                            file.path(outdir, paste(tpLevs[i],
                                                    "-PickedVariables.txt", sep="")),
                            quote=FALSE, row.names=FALSE)

                if ( nrow(pickedVars) > 0 ) {
                    res <- GetRandomPartitionPredictions(
                        outdir=outdir,
                        rmseFile=paste(tpLevs[i],
                                       "-BestPred-CompRmse.txt", sep=""),
                        fitsFile=paste(tpLevs[i],
                                       "-BestPred-Fits.rdata", sep=""))
                    rmses[nrow(rmses)+1,] <- c(target, dataName, "Best", tpLevs[i],
                                               nrow(pickedVars), BuildRmseRow_CoxPh(res))
                }
            }
        }
    }
    rmses
}

rmses.noCh <- CollectResults_CoxPh("NoChallenges", file.path(dataFolder, noChallengesFolder))
for ( i in 5:13 ) {
    rmses.noCh[,i] <- as.numeric(rmses.noCh[,i])
}
write.table(rmses.noCh, file.path(dataFolder, noChallengesFolder, "rmseSummary.txt"),
            quote=FALSE, row.names=FALSE)

rmses.noCh.summary <- rmses.noCh %>%
    filter(predictors=="Best") %>%
    select(data, timepoints, nvar,
           pctImproveAverage, pctImproveGroup,
           pctImproveAverage.glm, pctImproveGroup.glm) %>%
    mutate(pctImproveAverage=round(100*pctImproveAverage, 1),
           pctImproveGroup=round(100*pctImproveGroup, 1),
           pctImproveAverage.glm=round(100*pctImproveAverage.glm, 1),
           pctImproveGroup.glm=round(100*pctImproveGroup.glm, 1))

ggplot(rmses.noCh.summary) +
    geom_point(aes(pctImproveGroup, pctImproveGroup.glm, color=timepoints, shape=data)) +
    geom_abline(slope=1, intercept=0, lty=2)



GetSelectedVars <- function(dataFolder, predictionFolder,
                            dataName, tpSetName) {
    read.table(file.path(dataFolder, predictionFolder, dataName, tpSetName,
                         paste(tpSetName, "-PickedVariables.txt", sep="")),
               header=T, stringsAsFactors=FALSE)
}

## 12.3% improvement over group prediction
GetSelectedVars(dataFolder, noChallengesFolder, "Lag", "AllPreTerminal")
##   tp            re             ag shortVarName nPicked
## 1  7 aRhIgG.PE.low SIVsm.E660.2A5      var1242     185

## 9.6% improvement
GetSelectedVars(dataFolder, noChallengesFolder, "Lag", "tp7-PeakVL")
##   tp            re              ag shortVarName nPicked
## 1  7 aRhIgG.PE.low SIVmac239.gp140      var1239     185
## 2  7 aRhIgG.PE.low  SIVsm.E660.2A5      var1242     177

## 6.7% improvement (also adjusted - 6.6% improvement)
GetSelectedVars(dataFolder, noChallengesFolder, "Lag", "tp1-1xDna")
##   tp            re                 ag shortVarName nPicked
## 1  1 aRhIgG.PE.low SIVcpz.EK505.gp120       var204     155
## 2  1         R2A.3    SIVmac239.gp120       var263     144
## 3  1    R2A.4.high  SIVmac251.BK.PR55       var288     119

## 9.6% improvement
GetSelectedVars(dataFolder, noChallengesFolder, "Adjusted", "tp4-PreBoost")
##   tp             re                 ag shortVarName nPicked
## 1  4          R3A.1           G145.146       var824     161
## 2  4  aRhIgG.PE.low           G145.146       var712     149
## 3  4  aRhIgG.PE.low SIVcpz.EK505.gp120       var720     149
## 4  4 aRhIgG.PE.high              C1.TR       var690     130
## 5  4  aRhIgG.PE.low    SIVmac239.gp130       var722     115


plot1 <- ggplot(lagData) + geom_point(aes(var1242, NoChallenges, color=GroupNm)) +
    labs(x="aRhIgG.PE.low, SIVsm.E660.2A5, lagged, Peak VL Timepoint (scaled)",
         y="Number of Challenges") +
    guides(color=guide_legend(title="Group"))

plot2 <- ggplot(adjData) + geom_point(aes(var824, NoChallenges, color=GroupNm)) +
    labs(x="R3A.1, G145.146, adjusted, Pre-Boost Timepoint (scaled)",
         y="Number of Challenges") +
    guides(color=guide_legend(title="Group"))

plot3 <- ggplot(adjData) + geom_point(aes(var712, NoChallenges, color=GroupNm)) +
    labs(x="aRhIgG.PE.low, G145.146, adjusted, Pre-Boost Timepoint (scaled)",
         y="Number of Challenges") +
    guides(color=guide_legend(title="Group"))

plot4 <- ggplot(adjData) + geom_point(aes(var720, NoChallenges, color=GroupNm)) +
    labs(x="aRhIgG.PE.low, SIVcpz.EK505.gp120, adjusted, Pre-Boost Timepoint (scaled)",
         y="Number of Challenges") +
    guides(color=guide_legend(title="Group"))

plot5 <- ggplot(adjData) + geom_point(aes(var690, NoChallenges, color=GroupNm)) +
    labs(x="aRhIgG.PE.high, C1.TR, adjusted, Pre-Boost Timepoint (scaled)",
         y="Number of Challenges") +
    guides(color=guide_legend(title="Group"))

plot6 <- ggplot(adjData) + geom_point(aes(var722, NoChallenges, color=GroupNm)) +
    labs(x="aRhIgG.PE.low, SIVmac239.gp130, adjusted, Pre-Boost Timepoint (scaled)",
         y="Number of Challenges") +
    guides(color=guide_legend(title="Group"))

setwd("~/Projects/VRC332/Plots/PredictNoChallenges/CommonAnimals/CoxPh")
ggsave("BestVar-lag-tp7.png", plot1, height=6, width=8)
ggsave("BestVar-adj-tp4.png", plot2, height=6, width=8)
