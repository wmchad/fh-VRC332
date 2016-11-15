source("~/Projects/VRC332/Code/fh-vrc332/Functions/BuildGroupData.r")

## Does not do anything with NA/negative values
GetVariableSetData <- function(fcData, vars, includeGroups=TRUE, response="LogPeakVL") {

    x <- fcData %>% select(one_of(as.character(vars)))
    y <- fcData[,response]

    varMeans <- apply(x, 2, mean)
    goodAnimalIndices <- apply(x, 1, function(x) { any(x>0) & any(x!=varMeans) }) & y >= 0

    if ( includeGroups ) {
        x <- cbind(BuildGroupData(fcData), x)
    }
    
    list(x=x[goodAnimalIndices,],
         y=y[goodAnimalIndices],
         groups=fcData$GroupNm[goodAnimalIndices],
         animalIds=fcData$AnimalID[goodAnimalIndices])
}

GetVariableSetData.delta <- function(fcData, vars, includeGroups=TRUE, response="LogPeakVL") {

    x <- fcData %>% select(one_of(as.character(vars)))
    y <- fcData[,response]

    varMeans <- apply(x, 2, mean)
    goodAnimalIndices <- apply(x, 1, function(x) { any(x!=varMeans) }) & y >= 0

    if ( includeGroups ) {
        x <- cbind(BuildGroupData(fcData), x)
    }
    
    list(x=x[goodAnimalIndices,],
         y=y[goodAnimalIndices],
         groups=fcData$GroupNm[goodAnimalIndices],
         animalIds=fcData$AnimalID[goodAnimalIndices])
}
