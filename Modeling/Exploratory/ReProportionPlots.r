fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetBaselineAdjustedData.r"))
source(file.path(fnFolder, "GetDeltaData.r"))
source(file.path(fnFolder, "ConvertToLongData.r"))

origData <- fcData.na[apply(fcData.na[,-(1:7)], 1, function(x){!any(x<0)}),]
adjData <- GetBaselineAdjustedData(origData, predSummary)
lagData <- GetDeltaData(origData, predSummary)

longData.orig <- ConvertToLongData(origData, predSummary) %>%
    filter(value >= 0) %>%
    select(AnimalId, GroupNm, tp, re, ag, value)

longData.adj <- ConvertToLongData(adjData, predSummary) %>%
    select(AnimalId, GroupNm, tp, re, ag, value)

longData.lag <- ConvertToLongData(lagData, predSummary) %>%
    select(AnimalId, GroupNm, tp, re, ag, value)

rm(origData)
rm(adjData)
rm(lagData)

animalIds <- unique(longData.orig$AnimalId)
ags <- unique(longData.orig$ag)
res <- unique(longData.orig$re)

propData <- longData.orig %>%
    group_by(AnimalId, tp, ag) %>%
    mutate(groupSum=sum(value)) %>%
    mutate(reProp=value/groupSum)

setwd("~/Projects/VRC332/Plots/Modeling/Exploratory/ReProportions")
for ( anId in animalIds ) {
    plotData <- propData %>% filter(AnimalId==anId)
    grp <- first(plotData$GroupNm)
    rePropPlot <- ggplot(plotData) +
        geom_line(aes(tp, reProp, group=paste(ag, re), color=re), alpha=0.6) +
        labs(x="Timepoint", y="Proportion",
             title=paste("Reagent proportion plot, animal ", anId, ", ", grp, sep=""))
    ggsave(paste("rePropPlot-", grp, "-", anId, ".png", sep=""), rePropPlot,
           height=6, width=12)
}



testData <- propData %>% filter(AnimalId==anId, re==res[2])
testGroup <- first(testData$GroupNm)

ggplot(testData) +
    geom_line(aes(tp, reProp, color=ag), alpha=0.6) +
    labs(x="Timepoint", y="Proportion",
         title=paste("Reagent proportion plot, animal ", anId, ", ", testGroup, sep=""))

testData <- propData %>% filter(AnimalId==anId, ag==ags[14])
ggplot(testData) + geom_area(aes(tp, reProp, color=re, fill=re), position="stack")
