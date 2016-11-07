require(dplyr)
require(tidyr)
setwd("~/Projects/VRC332/Data")
predSummary = read.table("predSummary.txt", header=TRUE)
predSummary <- predSummary %>%
    filter(substr(varName, 1, 1)=="X") %>%
    separate(varName, c("tp", "re", "ag"), "_") %>%
    mutate(tp=substr(tp, 2, 2))
fcData <- read.table("longAnimalData2.narm.txt", header=TRUE)
fcData <- cbind(fcData[,1:7], fcData %>%
                              select(one_of(as.character(predSummary$shortVarName))))
fcData[,-(1:7)] <- ReplaceMissingWithMean(fcData[,-(1:7)])

tpLabels <- c("Baseline", "1xDna", "2xDna", "PeakPostPrime",
              "PreBoost", "PeakPostBoost", "PreChallenge",
              "PostChallenge", "Terminal")
