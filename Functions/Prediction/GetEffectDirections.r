require(caret)
require(glmnet)
require(caret)
require(dplyr)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "RunRandomPartitionPredictions.r"))

GetResultBetas <- function(results) {
    varNames <- c("Intercept", rownames(results$fits[[1]]$fit$finalModel$beta))
    nFits <- length(results$fits)
    dfBetas <- as.data.frame(matrix(ncol=length(varNames), nrow=nFits))
    names(dfBetas) <- varNames
    for ( i in 1:nFits ) {
        model <- results$fits[[i]]$fit$finalModel
        coefs <- as.data.frame(summary(predict(model, type="coef",
                                               s=model$lambdaOpt)))[,c(1,3)]
        dfBetas[i, coefs$i] <- coefs$x
    }
    apply(dfBetas, 2, mean, na.rm=T)
}

GetEffectDirections <- function(data, target, outputFolder, tps=0:8) {
    for ( tp in tps ) {
        outdir <- file.path(outputFolder,
                            paste("tp", tp, "-", tpLabels[tp+1], sep=""))
        rmseFile <- paste("tp", tp, "-", tpLabels[tp+1],
                          "-BestPred-CompRmse.txt", sep="")
        fitsFile <- paste("tp", tp, "-", tpLabels[tp+1],
                          "-BestPred-Fits.rdata", sep="")
        GetRandomPartitionPredictions(outdir, rmseFile, fitsFile) %>%
            GetResultBetas()
    }

    for ( i in 1:length(tpLevels) ) {
        outdir <- file.path(outputFolder, tpLevels[i])
        rmseFile <- paste(tpLevels[i], "-BestPred-CompRmse.txt", sep="")
        fitsFile <- paste(tpLevels[i], "-BestPred-Fits.rdata", sep="")
        GetRandomPartitionPredictions(outdir, rmseFile, fitsFile) %>%
            GetResultBetas()
    }
}
