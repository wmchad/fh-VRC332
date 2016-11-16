SummarizePredictionResults <- function(results, verbose=TRUE) {
    catn <- function(..., sep=" ") {}
    vprint <- function(...) {}
    if ( verbose ) {
        catn <- function(..., sep=" ") {cat(..., "\n", sep=sep)}
        vprint <- function(...) {print(...)}
    }
    rmses <- apply(results$compRmse, 2, mean)
    catn("Average RMSE:")
    vprint(rmses)
    catn()

    PickedCoeffs <- NULL
    nPicked <- NULL
    nRep <- length(results$fits)
    for ( i in 1:nRep ) {
        model <- results$fits[[i]]$fit$finalModel
        coeffs <- as.data.frame(summary(predict(model, type="coef", s=model$lambdaOpt)))[,c(1,3)]
        coeffNames <- c("(Intercept)", rownames(model$beta))[coeffs$i]
        PickedCoeffs <- c(PickedCoeffs, coeffNames)
        nPicked <- c(nPicked, length(coeffNames))
    }
    catn(paste("Model sizes range from", min(nPicked), "to", max(nPicked)))
    catn()
    
    pickedTbl <- table(PickedCoeffs)
    pickedTbl <- sort(pickedTbl, decreasing=TRUE)
    pickedTbl <- pickedTbl[substr(names(pickedTbl), 1, 3) == "var"]
    vprint(pickedTbl[pickedTbl >= nRep/2])
    list(rmses=rmses,
         pickedVars=pickedTbl,
         nPicked=sum(pickedTbl >= nRep/2))
}
