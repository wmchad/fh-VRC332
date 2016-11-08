require(dplyr)
require(reshape2)

ConvertToLongData <- function(fcData, predSummary) {
    longData <- melt(fcData[,-(3:4)],
                     id.vars=c("AnimalID2", "GroupNm", "NoChallenges",
                               "LogPeakVL", "LogSetpointVL"))
    longData <- merge(longData,
                      predSummary[,c("tp", "re", "ag", "shortVarName")],
                      by.x="variable", by.y="shortVarName", all.x=TRUE)
    longData[,-1]
}
