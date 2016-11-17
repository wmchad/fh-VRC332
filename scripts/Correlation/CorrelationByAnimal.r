require(ggplot2)
library(corrplot)
library(RColorBrewer)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetBaselineAdjustedData.r"))

animalData <- fcData.na %>% select(AnimalId, GroupNm, var1:var1548) %>%
    arrange(GroupNm)

animalData[animalData < 0] <- NA

animalCor <- cor(t(animalData[,-(1:2)]), use="complete.obs")

grpCols <- c(brewer.pal(5, "Dark2"))
names(grpCols) <- unique(animalData$GroupNm)
groupColors <- grpCols[animalData$GroupNm]

setwd("~/Projects/VRC332/Plots/FcCorrelation/ByAnimal")
png("FcCorrelationsByAnimal-alltp.png", width=2200, height=2000)
corrplot(animalCor, "square", tl.cex=1.5, tl.col=groupColors, cl.cex=3)
dev.off()

for ( t in 0:8 ) {
    tpVars <- (predSummary %>% filter(tp==t) %>% select(shortVarName))$shortVarName
    tpData <- animalData %>% select(one_of(tpVars))
    goodIndices <- apply(tpData, 1, function(x){any(!is.na(x))})
    tpData <- tpData[goodIndices,]
    tpCor <- cor(t(tpData), use="complete.obs")
    png(paste("FcCorrelationsByAnimal-tp", t, "-", tpLabels[t+1], ".png", sep=""),
        width=2200, height=2000)
    corrplot(tpCor, "square", tl.cex=1.5, tl.col=groupColors[goodIndices], cl.cex=3)
    dev.off();
}


adjData <- GetBaselineAdjustedData(fcData.na, predSummary) %>% arrange(GroupNm)

adjCor <- cor(t(adjData[,-(1:5)]), use="complete.obs")

grpCols <- c(brewer.pal(5, "Dark2"))
names(grpCols) <- unique(adjData$GroupNm)
groupColors <- grpCols[adjData$GroupNm]

setwd("~/Projects/VRC332/Plots/FcCorrelation/ByAnimalAdj")
png("FcCorrelationsByAnimal-alltp.png", width=2200, height=2000)
corrplot(adjCor, "square", tl.cex=1.5, tl.col=groupColors, cl.cex=3)
dev.off()

for ( t in 1:8 ) {
    tpVars <- (predSummary %>% filter(tp==t) %>% select(shortVarName))$shortVarName
    tpData <- adjData %>% select(one_of(tpVars))
    goodIndices <- apply(tpData, 1, function(x){any(!is.na(x))})
    tpData <- tpData[goodIndices,]
    tpCor <- cor(t(tpData), use="complete.obs")
    png(paste("FcCorrelationsByAnimal-tp", t, "-", tpLabels[t+1], ".png", sep=""),
        width=2200, height=2000)
    corrplot(tpCor, "square", tl.cex=1.5, tl.col=groupColors[goodIndices], cl.cex=3)
    dev.off();
}
