require(dplyr)
require(tidyr)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "ReplaceMissingWithMean.r"))

setwd("~/Projects/VRC332/Data/Processed")
predSummary = read.table("predSummary.txt", header=TRUE, stringsAsFactors=FALSE)
fcData <- read.table("longAnimalData.log.narm.txt", header=TRUE, stringsAsFactors=FALSE)
## Remove MBL
mblVars <- (predSummary %>%
    filter(re=="MBL") %>%
    select(shortVarName))$shortVarName
fcData <- fcData %>% select(-one_of(mblVars))
predSummary <- predSummary %>% filter(re!="MBL")
fcData.na <- fcData # missing values are set to -1
fcData[,-(1:7)] <- ReplaceMissingWithMean(fcData[,-(1:7)])

tpLabels <- c("Baseline", "1xDna", "2xDna", "PeakPostPrime",
              "PreBoost", "PeakPostBoost", "PreChallenge",
              "PeakVL", "Terminal")
