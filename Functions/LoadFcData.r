require(dplyr)
require(tidyr)
setwd("~/Projects/VRC332/Data")
fcData = read.table("longAnimalData2.narm.txt", header=TRUE)
predSummary = read.table("predSummary.txt", header=TRUE)
predSummary <- predSummary %>%
    filter(substr(varName, 1, 1)=="X") %>%
    separate(varName, c("tp", "re", "ag"), "_") %>%
    mutate(tp=substr(tp, 2, 2))


tpLabels <- c("Baseline", "1xDna", "2xDna", "PeakPostPrime",
              "PreBoost", "PeakPostBoost", "PreChallenge",
              "PostChallenge", "Terminal")
