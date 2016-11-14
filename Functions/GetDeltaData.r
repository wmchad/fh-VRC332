require(dplyr)
source("~/Projects/VRC332/Code/fh-vrc332/Functions/GetTimepointData.r")

## Does not do anything with NA/negative values
GetDeltaData <- function(fcData.na, predSummary) {
    ## baseline data
    tpPrevData <- fcData.na %>%
        select(one_of(as.character(predSummary$shortVarName[predSummary$tp %in% 0:7])))
    tpSubsData <- fcData.na %>%
        select(one_of(as.character(predSummary$shortVarName[predSummary$tp %in% 1:8])))
    allDataIndices <- apply(fcData.na[,-(1:7)], 1, function(x) { all(x >= 0) })
    tpPrevData <- tpPrevData[allDataIndices,]
    tpSubsData <- tpSubsData[allDataIndices,]
    deltaData <- tpSubsData - tpPrevData
    cbind(fcData.na[allDataIndices,1:7], deltaData)
}
