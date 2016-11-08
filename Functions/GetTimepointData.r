require(dplyr)
source("~/Projects/VRC332/Code/fh-vrc332/Functions/BuildGroupData.r")

## Does not do anything with NA/negative values
GetTimepointData <- function(fcData, tps, predSummary,
                             includeGroups=TRUE, response="LogPeakVL") {
    predSub <- predSummary %>% filter(tp %in% tps)

    x <- fcData %>% select(one_of(as.character(predSub$shortVarName)))
    y <- fcData[,response]

    varMeans <- apply(x, 2, mean)
    goodAnimalIndices <- apply(x, 1, function(x) { any(x>0) & any(x!=varMeans) }) & y >= 0

    if ( includeGroups ) {
        x <- cbind(BuildGroupData(fcData), x)
    }
    
    list(x=x[goodAnimalIndices,],
         y=y[goodAnimalIndices],
         groups=fcData$GroupNm[goodAnimalIndices],
         animalIds=as.character(fcData$AnimalId[goodAnimalIndices]))
}
