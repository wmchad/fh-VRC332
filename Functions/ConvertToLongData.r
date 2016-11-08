require(dplyr)
require(reshape2)

ConvertToLongData <- function(fcData, predSummary) {
    longData <- melt(fcData[,-(3:4)],
                     id.vars=c("AnimalId", "GroupNm", "NoChallenges",
                               "LogPeakVL", "LogSetpointVL"))
    longData <- left_join(longData,
                      predSummary %>% select(tp:ag, shortVarName),
                      by=c("variable"="shortVarName"))
    longData
}
