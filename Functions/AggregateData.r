fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"

source(file.path(fnFolder, "ConvertToLongData.r"))

AggregateData <- function(fcData, predSummary) {
    longData <- ConvertToLongData(fcData, predSummary)

    aggAgData <- longData %>%
        mutate(tp=paste("tp", tp, sep="")) %>%
        unite(aggName, tp, ag) %>%
        group_by(AnimalId, aggName) %>%
        summarize(avg=mean(value)) %>%
        spread(aggName, avg)
    fcAgData <- fcData[,1:7] %>% left_join(aggAgData)

    aggReData <- longData %>%
        mutate(tp=paste("tp", tp, sep="")) %>%
        unite(aggName, tp, re) %>%
        group_by(AnimalId, aggName) %>%
        summarize(avg=mean(value)) %>%
        spread(aggName, avg)
    fcReData <- fcData[,1:7] %>% left_join(aggReData)

    fcAggData <- fcAgData %>% left_join(aggReData)

    aggPredSummary <- data.frame(varName=names(fcAggData)[-(1:7)],
                                 shortVarName=names(fcAggData)[-(1:7)])
    aggPredSummary <- aggPredSummary %>%
        separate(varName, c("tp", "pred"), "_") %>%
        mutate(tp=substr(tp, 3, 3))
    list(aggData=fcAggData,
         predSummary=aggPredSummary)
}
