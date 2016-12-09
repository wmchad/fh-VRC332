source("~/Projects/VRC332/Code/fh-vrc332/Functions/GetVariableSetData.r")

GetBestPredictorsData_CoxPh <- function(data, prevResults, nCutoff,
                                        includeGroups=TRUE, response="LogPeakVL",
                                        varIndicator="var") {
  allPickedCoeffs <- NULL
  nPicked <- NULL
  for ( i in 1:length(prevResults$fits) ) {
    allPickedCoeffs <- c(allPickedCoeffs, prevResults$fits[[i]]$vars)
    nPicked <- c(nPicked, length(prevResults$fits[[i]]$vars))
  }
  pickedTbl <- table(allPickedCoeffs)
  pickedTbl <- sort(pickedTbl, decreasing=TRUE)

  vars <- names(pickedTbl)[pickedTbl >= nCutoff]
  vars <- vars[grepl(varIndicator, vars)]

  GetVariableSetData.delta(data, vars, includeGroups, response)
}
