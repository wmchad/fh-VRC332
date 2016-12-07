require(ggplot2)

CompLinePlots <- function(compRmse) {
  compRmse <- as.data.frame(compRmse)
  baseComp <- ggplot(compRmse) + geom_point(aes(baseline, predicted)) +
    geom_abline(intercept=0, slope=1, lty=2) +
    geom_segment(aes(x=baseline, y=baseline, xend=baseline, yend=predicted,
                     color=factor(sign(baseline-predicted)))) +
    guides(color=FALSE) +
    labs(x="Baseline RMSE", y="Predicted RMSE")

  groupComp <- ggplot(compRmse) + geom_point(aes(group, predicted)) +
    geom_abline(intercept=0, slope=1, lty=2) +
    geom_segment(aes(x=group, y=group, xend=group, yend=predicted,
                     color=factor(sign(group-predicted)))) +
    guides(color=FALSE) +
    labs(x="Group RMSE", y="Predicted RMSE")

  list(baselineComp=baseComp,
       groupComp=groupComp)
}








target <- "SetpointVL"
targetFolder <- setpointVLFolder
dataName <- "Original"
i <- 3

outdir <- file.path(dataFolder, targetFolder, dataName, tpLevs[i])
resSummary <- SummarizePredictionResults(
    GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste(tpLevs[i], "-CompRmse.txt", sep=""),
        fitsFile=paste(tpLevs[i], "-Fits.rdata", sep="")),
    FALSE)
resSummary2 <- SummarizePredictionResults(
    GetRandomPartitionPredictions(
        outdir=outdir,
        rmseFile=paste(tpLevs[i], "-BestPred-CompRmse.txt", sep=""),
        fitsFile=paste(tpLevs[i], "-BestPred-Fits.rdata", sep="")),
    FALSE)
