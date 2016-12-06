require(dplyr)
require(stringr)
require(ggplot2)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "LoadFcData.r"))
source(file.path(fnFolder, "GetBaselineAdjustedData.r"))
source(file.path(fnFolder, "GetDeltaData.r"))
source(file.path(fnFolder, "ConvertToLongData.r"))

BuildFitCoefData <- function(fit, tps, grps, animalIds, ags, res) {
    coefDf <- as.data.frame(coef(summary(fit)))
    names(coefDf) <- c("beta", "se", "tval", "pval")
    coefDf$var <- rownames(coefDf)
    coefDf$tp <- NA
    for ( tp in tps ) {
        coefDf$tp[grep(paste("as.factor\\(tp\\)", tp, sep=""), coefDf$var)] <- tp
        coefDf$tp[grep(paste("tp", tp, sep=""), coefDf$var)] <- tp
    }
    coefDf$group <- NA
    for ( grp in grps ) {
        coefDf$group[grep(grp, coefDf$var)] <- grp
    }
    coefDf$animal <- NA
    for ( anId in animalIds ) {
        coefDf$animal[grep(anId, coefDf$var)] <- anId
    }
    coefDf$ag <- NA
    for ( ag in ags ) {
        coefDf$ag[grep(ag, coefDf$var)] <- ag
    }
    coefDf$re <- NA
    for ( re in res ) {
        coefDf$re[grep(re, coefDf$var)] <- re
    }
    coefDf$var <- apply(coefDf %>% select(tp:re), 1, function(x)
        paste(na.omit(x), collapse=":"))
    predGrpDf <- data.frame(tp=ifelse(is.na(coefDf$tp), NA, "Timepoint"),
                            grp=ifelse(is.na(coefDf$group), NA, "Group"),
                            animal=ifelse(is.na(coefDf$animal), NA, "Animal"),
                            ag=ifelse(is.na(coefDf$ag), NA, "Antigen"),
                            re=ifelse(is.na(coefDf$re), NA, "Reagent"))
    coefDf$predictorGroup <- apply(predGrpDf, 1, function(x)
        paste(na.omit(x), collapse=":"))
    coefDf
}

PlotResids <- function(data, ags, res, outdir) {
    basePlot <- ggplot(data)
    setwd(outdir)
    ggsave("ResidDensity-tp.png",
           basePlot +
           geom_density(aes(resid, color=factor(tp))) +
           labs(title="Residual density per timepoint"),
           height=6, width=10)
    ggsave("ResidDensity-ag.png",
           basePlot +
           geom_density(aes(resid, color=ag)) +
           labs(title="Residual density per antigen"),
           height=6, width=10)
    ggsave("ResidDensity-re.png",
           basePlot +
           geom_density(aes(resid, color=re)) +
           labs(title="Residual density per detection reagent"),
           height=6, width=10)
    ggsave("ResidDensity-group.png",
           basePlot +
           geom_density(aes(resid, color=GroupNm)) +
           labs(title="Residual density per treatment group"),
           height=6, width=10)

    for ( agTest in ags ) {
        agPlot <- ggplot(data %>% filter(ag==agTest)) +
            geom_density(aes(resid, color=re)) +
            labs(title=paste("Residual density for antigen", agTest))
        ggsave(paste("ResidDensity-ag-", agTest, ".png", sep=""), agPlot,
               height=6, width=10)
    }
    for ( reTest in res ) {
        rePlot <- ggplot(data %>% filter(re==reTest)) +
            geom_density(aes(resid, color=ag)) +
            labs(title=paste("Residual density for reagent", reTest))
        ggsave(paste("ResidDensity-re-", reTest, ".png", sep=""), rePlot,
               height=6, width=10)
    }
}

PlotVsTruth <- function(data, outdir) {
    basePlot <- ggplot(data)
    setwd(outdir)
    ggsave("ResidVsTruth.png",
           basePlot +
           geom_point(aes(value, resid), alpha=0.1) +
           labs(title="Residuals Vs Response",
                x="log(readout)",
                y="Residual"),
           height=6, width=10)
    ggsave("PredictedVsTruth.png",
           basePlot +
           geom_point(aes(value, predicted), alpha=0.1) +
           geom_abline(slope=1, intercept=0, lty=2, col="green") +
           labs(title="Predicted Vs Response",
                x="log(readout)",
                y="Predicted Value"),
           height=6, width=10)
}

FitName <- function(dataName, modelName) {
    paste(modelName, "fit", dataName, "rds", sep=".")
}

RunAndSummarizeModel <- function(data, dataName,
                                 modelFormula, modelName,
                                 dataDir, tps,
                                 grps, animalIds,
                                 ags, res) {
    fit <- lm(modelFormula, data=data, model=FALSE)
    setwd(dataDir)
    saveRDS(fit, FitName(dataName, modelName))
    coefs <- BuildFitCoefData(fit, tps, grps, animalIds, ags, res)
    data <- data %>% mutate(predicted=predict(fit, newdata=data),
                            resid=resid(fit))
    rm(fit)
    list(data=data, coefs=coefs)
}

ModelSummaryPlots <- function(data, plotDir, ags, res) {
    PlotResids(data, ags, res, plotDir)
    PlotVsTruth(data, plotDir)
    setwd(plotDir)
    png("qqplot.png", height=600, width=600)
    qqnorm(data$resid)
    abline(0, 1, lty=2)
    dev.off()
}

GetModelFit <- function(dataName, modelName, dataDir) {
    readRDS(file.path(dataDir, FitName(dataName, modelName)))
}

GetModelSummary <- function(data, dataName, modelName, dataDir,
                            tps, grps, animalIds, ags, res) {
    fit <- GetModelFit(dataName, modelName, dataDir)
    coefs <- BuildFitCoefData(fit, tps, grps, animalIds, ags, res)
    data <- data %>% mutate(predicted=predict(fit, newdata=data),
                            resid=resid(fit))
    rm(fit)
    list(data=data, coefs=coefs)
}
