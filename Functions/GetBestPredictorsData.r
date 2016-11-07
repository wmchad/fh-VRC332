source("~/Projects/VRC332/Code/fh-vrc332/Functions/GetVariableSetData.r")

GetBestPredictorsData <- function(vlData, prevResults, nCutoff,
                                  includeGroups=TRUE, response="LogPeakVL") {
  allPickedCoeffs <- NULL
  nPicked <- NULL
  for ( i in 1:length(prevResults$fits) ) {
    model <- prevResults$fits[[i]]$fit$finalModel

    coeffs <- as.data.frame(summary(predict(model, type="coef", s=model$lambdaOpt)))[,c(1,3)]
    coeffNames <- c("(Intercept)", rownames(model$beta))[coeffs$i]
    allPickedCoeffs <- c(allPickedCoeffs, coeffNames)
    nPicked <- c(nPicked, length(coeffNames))
  }
  pickedTbl <- table(allPickedCoeffs)
  pickedTbl <- sort(pickedTbl, decreasing=TRUE)

  vars <- names(pickedTbl)[pickedTbl >= nCutoff]
  vars <- vars[grepl("var", vars)]

  GetVariableSetData(vlData, vars, includeGroups, response)
}
