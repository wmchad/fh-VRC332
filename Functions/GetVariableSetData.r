source("~/Projects/VRC332/Code/fh-vrc332/Functions/BuildGroupData.r")

## Does not do anything with NA/negative values
GetVariableSetData <- function(fcData, vars, includeGroups=TRUE, response="LogPeakVL") {

  x <- NULL
  if ( includeGroups ) {
    x <- BuildGroupData(fcData)
  }
  x <- cbind(x, fcData[,vars])
  y <- fcData[,response]

  goodAnimalIndices <- apply(x, 1, function(x) { any(x[-(1:4)]>0) })
  
  list(x=x[goodAnimalIndices,],
       y=y[goodAnimalIndices],
       groups=fcData$GroupNm[goodAnimalIndices],
       animalIds=fcData$AnimalID[goodAnimalIndices])
}
