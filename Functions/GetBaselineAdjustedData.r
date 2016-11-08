require(dplyr)
source("~/Projects/VRC332/Code/fh-vrc332/Functions/GetTimepointData.r")

## Does not do anything with NA/negative values
GetBaselineAdjustedData <- function(fcData.na, predSummary) {
    ## baseline data
    tp0Data <- fcData.na %>%
        select(one_of(as.character(predSummary$shortVarName[predSummary$tp==0])))
    someDataIndices <- apply(tp0Data, 1, function(x) { any(x >= 0) })
    tp0Data <- tp0Data[someDataIndices,]
    
    ## post-baseline data
    pbData <- fcData.na %>%
        select(one_of(as.character(predSummary$shortVarName[predSummary$tp %in% 1:8])))
    pbData <- pbData[someDataIndices,]

    for ( varName in names(pbData) ) {
        agRe <- predSummary %>% filter(shortVarName==varName) %>% select(ag, re)
        var0 <- as.character(predSummary %>%
                             filter(tp==0, ag==agRe$ag, re==agRe$re) %>%
                             select(shortVarName))
        missingIndices <- pbData[,varName] < 0 | tp0Data[,var0] < 0
        pbData[,varName] <- pbData[,varName] - tp0Data[,var0]
        pbData[missingIndices, varName] <- 0
    }
    cbind(fcData.na[someDataIndices,1:7], pbData  )
}
