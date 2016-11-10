require(dplyr)
require(tidyr)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))

setwd("~/Projects/VRC332/Data")
predSummary = read.table("predSummary.txt", header=TRUE, stringsAsFactors=FALSE)
predSummary <- predSummary %>%
    filter(substr(varName, 1, 1)=="X") %>%
    separate(varName, c("tp", "re", "ag"), "_") %>%
    mutate(tp=substr(tp, 2, 2))
fcData <- read.table("longAnimalData2.narm.txt", header=TRUE, stringsAsFactors=FALSE)
names(fcData)[1] <- "AnimalId"
fcData <- cbind(fcData[,1:7], fcData %>%
                              select(one_of(as.character(predSummary$shortVarName))))
fcData.na <- fcData # missing values are set to -1
fcData[,-(1:7)] <- ReplaceMissingWithMean(fcData[,-(1:7)])

tpLabels <- c("Baseline", "1xDna", "2xDna", "PeakPostPrime",
              "PreBoost", "PeakPostBoost", "PreChallenge",
              "PostChallenge", "Terminal")
