source("~/Projects/VRC332/Code/fh-vrc332/Functions/BuildGroupData.r")

## Does not do anything with NA/negative values
GetVariableSetData <- function(vlData, vars, includeGroups=TRUE, response="LogPeakVL") {

  x <- NULL
  if ( includeGroups ) {
    x <- BuildGroupData(vlData)
  }
  x <- cbind(x, vlData[,vars])
  y <- vlData[,response]

  goodAnimalIndices <- apply(x, 1, function(x) { any(x[-(1:4)]>0) })
  
  list(x=x[goodAnimalIndices,],
       y=y[goodAnimalIndices],
       groups=vlData$GroupNm[goodAnimalIndices],
       animalIds=vlData$AnimalID[goodAnimalIndices])
}
