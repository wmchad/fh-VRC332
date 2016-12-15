ggplot(results.orig$coef %>% filter(predictorGroup=="Reagent")) + geom_point(aes(re, beta))

ggplot(results.orig$coef %>% filter(predictorGroup=="Antigen")) + geom_point(aes(ag, beta))

betaHeat <- ggplot(results.orig5$coef) + geom_tile(aes(ag, re, fill=beta)) +
    scale_fill_gradient2() + theme(axis.text=element_text(size=4))

setwd("~/Projects/VRC332/Plots/Modeling/LinearModels/Model5/Orig/Exploratory")
ggsave("betaHeat.png", betaHeat, width=12, height=6)

simpleTbl <- results.orig5$coef %>% filter(predictorGroup=="Antigen:Reagent") %>%
    group_by(ag, re) %>% summarize(val=beta)

ags <-unique(results.orig5$coef$ag)
res <- unique(results.orig5$coef$re)

betaMatrix <- matrix(nrow=length(ags),
                     ncol=length(res),
                     data=0)
rownames(betaMatrix) <- ags
colnames(betaMatrix) <- res

for ( i in 1:length(ags) ) {
    for ( j in 1:length(res) ) {
        subtbl <- simpleTbl %>% filter(ag==ags[i], re==res[j])
        if ( nrow(subtbl) == 1 ) {
            betaMatrix[i,j] <- subtbl$val
        }
    }
}

betaMatrix2 <- betaMatrix
for ( i in 1:length(ags) ) {
    agMean <- mean(betaMatrix2[i,betaMatrix2[i,] != 0])
    betaMatrix2[i,betaMatrix2[i,] == 0] <- agMean
}

betaMatrix <- betaMatrix[-1,]
betaMatrix2 <- betaMatrix2[-1,]


agDist <- dist(betaMatrix)
agClust <- hclust(agDist)

agDist2 <- dist(betaMatrix2)
agClust2 <- hclust(agDist2)

setwd("~/Projects/VRC332/Plots/Modeling/LinearModels/Model5/Orig/Exploratory")

png("agClust1.png", width=1200, height=800)
plot(agClust)
dev.off()

png("agClust2.png", width=1200, height=800)
plot(agClust2)
dev.off()



library(formula.tools)

modelComp <- data.frame(data      = character(0),
                        modelName = character(0),
                        model     = character(0),
                        rse       = numeric(0),
                        r2        = numeric(0),
                        nSig      = numeric(0),
                        stringsAsFactors = FALSE)

for ( modelName in paste("model", 1:9, sep="") ) {
    for ( dataName in c("orig", "adj", "lag") ) {
        fit <- summary(GetModelFit(dataName, modelName, dataDir))
        modelComp[nrow(modelComp)+1,1:3] <- c(dataName, modelName,
                                              as.character(get(modelName)))
        modelComp[nrow(modelComp),4:6] <- c(fit$sigma, fit$r.squared,
                                            sum(fit$coefficients[,4] <= 0.05))
    }
}

## time-series clustering
modelComp[nrow(modelComp)+1,1:3] <- c("adj", "hclust", "Time series hclust")
modelComp[nrow(modelComp),4:6] <- c(0.903, 0.802, 0)
modelComp[nrow(modelComp)+1,1:3] <- c("lag", "hclust", "Time series hclust")
modelComp[nrow(modelComp),4:6] <- c(0.942, 0.627, 0)


## null model
modelComp[nrow(modelComp)+1,1:3] <- c("adj", "null", "value ~ 1")
modelComp[nrow(modelComp),4:6] <- c(2.036, 0, 0)
modelComp[nrow(modelComp)+1,1:3] <- c("lag", "null", "value ~ 1")
modelComp[nrow(modelComp),4:6] <- c(1.541, 0, 0)
modelComp[nrow(modelComp)+1,1:3] <- c("orig", "null", "value ~ 1")
modelComp[nrow(modelComp),4:6] <- c(2.245, 0, 0)

setwd(dataDir)
write.table(modelComp, "ModelComparison.txt", row.names=FALSE)



modelComp <- data.frame(data  = sort(rep(c("Original", "Adjusted", "Lagged"), 7)),
                        model = c("value ~ 1",
                                  "value ~ Group * tp + Animal * (ag + re)",
                                  "value ~ Group * (tp + ag + re) + Animal * (ag + re)",
                                  "value ~ Group * (tp + ag + re)",
                                  "value ~ Group * (tp + ag + re) + Animal",
                                  "value ~ ag * re",
                                  "value ~ ag * re * tp"),
                        rse   = 0,
                        r2    = 0,
                        nSig  = 0,
                        stringsAsFactors = FALSE)

## Baseline-Adjusted
modelComp[1,3:5] <- c(2.036, 0,      0)
modelComp[2,3:5] <- c(1.308, 0.5966, 668)
modelComp[3,3:5] <- c(1.308, 0.5966, 402)
modelComp[4,3:5] <- c(1.347, 0.563,  122)
modelComp[5,3:5] <- c(1.333, 0.5725, 158)
modelComp[6,3:5] <- c(1.739, 0.2715, 103)
modelComp[7,3:5] <- c(1.219, 0.6469, 450)

## Lagged
modelComp[8,3:5] <-  c(1.541, 0,      0)
modelComp[9,3:5] <-  c(1.279, 0.3277, 37)
modelComp[10,3:5] <-  c(1.279, 0.3277, 33)
modelComp[11,3:5] <- c(1.272, 0.3204, 47)
modelComp[12,3:5] <- c(1.269, 0.3244, 66)
modelComp[13,3:5] <- c(1.514, 0.0376, 37)
modelComp[14,3:5] <- c(1.163, 0.4391, 304)

## Original
modelComp[15,3:5] <- c(2.245, 0,      0)
modelComp[16,3:5] <- c(1.466, 0.5821, 500)
modelComp[17,3:5] <- c(1.466, 0.5821, 189)
modelComp[18,3:5] <- c(1.478, 0.5674, 124)
modelComp[19,3:5] <- c(1.47,  0.5724, 164)
modelComp[20,3:5] <- c(1.678, 0.4423, 136)
modelComp[21,3:5] <- c(1.107, 0.7604, 449)

## time-series clustering
modelComp[22,1:2] <- c("Adjusted", "Time series hclust")
modelComp[22,3:5] <- c(0.903, 0.802, 0)
modelComp[23,1:2] <- c("Lagged", "Time series hclust")
modelComp[23,3:5] <- c(0.942, 0.627, 0)
