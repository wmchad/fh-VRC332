library(VRC332)
library(tidyr)
library(dplyr)
library(stringr)

## Load the base datasets
data("FcCoreAssay")
names(FcCoreAssay)[1] <- "AnimalId"
names(FcCoreAssay)[10] <- "NoChallenges"

## Get the animal-specific information
Animals = FcCoreAssay %>%
    select(AnimalId, GroupNm, PeakVL, SetpointVL, NoChallenges, Ack_Group) %>%
    distinct()
Animals$GroupNm <- str_replace_all(Animals$Ack_Group, " ", "_")

Animals %>% filter(is.na(GroupNm))
##      AnimalId GroupNm   PeakVL SetpointVL NoChallenges Ack_Group
## 1         8-2    <NA>     <NA>         NA           NA      <NA>
## 2        8-20    <NA>     <NA>         NA           NA      <NA>
## 3        8-24    <NA> 12665771      57958           NA      <NA>
## 4 T-4276_blip    <NA>     NULL         37           NA      <NA>

Animals <- Animals %>% mutate(PeakVL=as.numeric(paste(PeakVL))) %>%
    mutate(LogPeakVL=log10(PeakVL), LogSetpointVL=log10(SetpointVL)) %>%
    filter(!is.na(GroupNm)) %>% select(AnimalId:NoChallenges, LogPeakVL, LogSetpointVL)

## Collect all the data in long format
FcCoreAssay <- FcCoreAssay %>% mutate(logreadout = log(readout))

longData.Fc <- FcCoreAssay %>%
    filter(AnimalId %in% Animals$AnimalId) %>%
    mutate(Key=paste("X", Timepoint, "_", DetectionReagent, "_", Antigen, sep="")) %>%
    select(Key, AnimalId, readout) %>%
    spread(Key, readout) %>%
    select(1:1549)

longData.Fclog <- FcCoreAssay %>%
    filter(AnimalId %in% Animals$AnimalId) %>%
    mutate(Key=paste("X", Timepoint, "_", DetectionReagent, "_", Antigen, sep="")) %>%
    select(Key, AnimalId, logreadout) %>%
    spread(Key, logreadout) %>%
    select(1:1549)

longData.Fc <- longData.Fc[apply(longData.Fc[,-1], 1, function(x){any(!is.na(x))}),]
longData.Fclog <- longData.Fclog[apply(longData.Fclog[,-1], 1, function(x){any(!is.na(x))}),]
Animals <- Animals %>% filter(AnimalId %in% longData.Fc$AnimalId)

## Change predictor names and create summary table
predSummary <- data.frame(varName=names(longData.Fc)[-1],
                          shortVarName=paste("var", 1:(ncol(longData.Fc)-1), sep=""),
                          stringsAsFactors=FALSE)

predSummary <- predSummary %>%
    filter(substr(varName, 1, 1)=="X") %>%
    separate(varName, c("tp", "re", "ag"), "_") %>%
    mutate(tp=substr(tp, 2, 2))

names(longData.Fc)[-1] <- predSummary$shortVarName
names(longData.Fclog)[-1] <- predSummary$shortVarName

## Combine all the data
longAnimalData = Animals %>% left_join(longData.Fc, by="AnimalId")
longAnimalData.log = Animals %>% left_join(longData.Fclog, by="AnimalId")

## Save the data
setwd("~/Projects/VRC332/Data/Processed")
write.table(longAnimalData, "longAnimalData.txt", row.names=FALSE)
write.table(longAnimalData.log, "longAnimalData.log.txt", row.names=FALSE)
write.table(Animals, "animalData.txt", row.names=FALSE)
write.table(longData.Fc, "longFcData.txt", row.names=FALSE)
write.table(longData.Fclog, "longFcData.log.txt", row.names=FALSE)
write.table(predSummary, "predSummary.txt", row.names=FALSE)

longAnimalData.narm = longAnimalData;
longAnimalData.narm[is.na(longAnimalData.narm)] = -1
write.table(longAnimalData.narm, "longAnimalData.narm.txt", row.names=FALSE)

longAnimalData.log.narm = longAnimalData.log;
longAnimalData.log.narm[is.na(longAnimalData.log.narm)] = -1
write.table(longAnimalData.log.narm, "longAnimalData.log.narm.txt", row.names=FALSE)
