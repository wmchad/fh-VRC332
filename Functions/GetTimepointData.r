source("~/Projects/VRC332/Code/fh-vrc332/Functions/BuildGroupData.r")

## Does not do anything with NA/negative values
GetTimepointData <- function(vlData, tps, predSummary,
                             includeGroups=TRUE, response="LogPeakVL") {
  predSummary <- predSummary %>%
    filter(substr(varName, 1, 1)=="X") %>%
    filter(as.numeric(substr(varName, 2, 2)) %in% tps)

  x <- NULL
  if ( includeGroups ) {
    x <- BuildGroupData(vlData)
  }
  x <- cbind(x, vlData %>% select(one_of(as.character(predSummary$shortVarName))))
  y <- vlData[,response]

  goodAnimalIndices <- apply(x, 1, function(x) { any(x[-(1:4)]>0) }) & y >= 0
  
  list(x=x[goodAnimalIndices,],
       y=y[goodAnimalIndices],
       groups=vlData$GroupNm[goodAnimalIndices],
       animalIds=vlData$AnimalID[goodAnimalIndices])
}
