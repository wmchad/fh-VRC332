require(rstan)
require(RColorBrewer)

options(mc.cores = parallel::detectCores())

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "Modeling/LinearModels-Support.r"))

log_sum_exp <- function(x) {
    maxx <- max(x)
    maxx + log(sum(exp(x-maxx)))
}

logsoftmax <- function(x) {
    x - log_sum_exp(x)
}

clust_assign <- function(x) {
    lsx <- logsoftmax(x)
    (1:length(lsx))[lsx == max(lsx)]
}

origData <- fcData.na[apply(fcData.na[,-(1:7)], 1, function(x){!any(x<0)}),]
adjData <- GetBaselineAdjustedData(origData, predSummary)

longData.orig <- ConvertToLongData(origData, predSummary) %>%
    filter(value >= 0) %>%
    select(AnimalId, GroupNm, tp, re, ag, value)

longData.adj <- ConvertToLongData(adjData, predSummary) %>%
    select(AnimalId, GroupNm, tp, re, ag, value)


rm(origData)
rm(adjData)

agretp <- longData.adj %>% select(ag, re, tp) %>% distinct()
agretp$groupId <- 1:nrow(agretp)

longData.adj <- longData.adj %>%
    left_join(agretp, by=c("ag" = "ag",
                            "re" = "re",
                            "tp" = "tp"))

N <- nrow(longData.adj)
N_subgroups <- length(unique(longData.adj$groupId))
obs_to_subgroup <- longData.adj$groupId
y <- longData.adj$value

model_dat <- list(N = N,
                  N_subgroups = N_subgroups,
                  obs_to_subgroup = obs_to_subgroup,
                  y = y)

setwd("~/Projects/VRC332/Code/fh-vrc332/Modeling/Stan/")
fit <- stan("model-beta_art.stan", data=model_dat, iter=1000, chains=4)

fit_summary <- as.data.frame(summary(fit)$summary)

agretp$mean <- fit_summary$mean[1:nrow(agretp)]

ggplot(agretp %>% filter(tp==5)) +
    geom_density(aes(mean, color=ag))

agHighlight <- "SIVmac239.gp140"

ggplot(agretp %>% group_by(ag, re)) +
    geom_line(aes(tp, mean, group=paste(ag, re)), alpha = 0.1) +
    geom_line(aes(tp, mean, group=paste(ag, re), color=re),
              agretp %>% filter(ag==agHighlight) %>% group_by(ag, re),
              size=1.2)

maxtp <- 8
agretp2 <- agretp %>% filter(tp <= maxtp)
longData.adj2 <- longData.adj %>% filter(tp <= maxtp)

N2 <- nrow(longData.adj2)
N_subgroups2 <- length(unique(longData.adj2$groupId))
obs_to_subgroup2 <- longData.adj2$groupId
y2 <- longData.adj2$value

model_dat2 <- list(N = N2,
                   N_subgroups = N_subgroups2,
                   obs_to_subgroup = obs_to_subgroup2,
                   N_t = maxtp,
                   obs_to_t = longData.adj2$tp,
                   y = y2)

setwd("~/Projects/VRC332/Code/fh-vrc332/Modeling/Stan/")
fit2 <- stan("model-beta_art_mix.stan", data=model_dat2, iter=1000, chains=4)
saveRDS(fit2, "model-beta_art_mix.rds")

fit2_summary <- as.data.frame(summary(fit2)$summary)
fit2_summary$var <- rownames(fit2_summary)
agretp$mean2 <- (fit2_summary %>% filter(str_sub(var, 1, 5) == "beta["))$mean
agretp$omega <- (fit2_summary %>% filter(str_sub(var, 1, 5) == "omega"))$mean

ggplot(agretp) + geom_point(aes(mean2, omega, color=factor(tp)))

ggplot(agretp) + geom_point(aes(mean, mean2, color=factor(tp), alpha=omega))

ggplot(agretp) + geom_density(aes(omega, color=ag))
ggplot(agretp) + geom_density(aes(omega, color=re))


sigBetas <- (fit2_summary %>%
             filter(sign(`2.5%`) == sign(`97.5%`),
                    str_sub(var, 1, 5) == "beta["))$var %>%
                                                  str_sub(6, -2)
sigBetas <- as.numeric(sigBetas)
sigBetas <- sigBetas[!is.na(sigBetas)]
agretp$isSig2 <- FALSE
agretp$isSig2[sigBetas] <- TRUE

betat <- (fit2_summary  %>% filter(str_sub(var, 1, 6)=="beta_t"))$mean
agretp$betat <- betat[agretp$tp]

ggplot(agretp) + geom_point(aes(mean, mean2+betat, color=factor(tp), shape=isSig2))

ggplot(agretp) + geom_density(aes(mean2, color=isSig2))

ggplot(agretp) + geom_density(aes(mean, color=isSig2))



agHighlight <- "SIVmac239.gp140"

ggplot(agretp %>% group_by(ag, re)) +
    geom_line(aes(tp, mean2+betat, group=paste(ag, re)), alpha = 0.1) +
    geom_line(aes(tp, mean2+betat, group=paste(ag, re), color=re),
              agretp %>% filter(isSig2, ag==agHighlight) %>% group_by(ag, re),
              size=1.2) +
    geom_line(aes(tp, mean, group=paste(ag, re), color=re),
              agretp %>% filter(ag==agHighlight) %>% group_by(ag, re),
              size=1, lty=2) +
  geom_line(aes(tp, meanval),
            longData.adj %>% group_by(tp) %>% summarize(meanval=mean(value)),
            size=2, lty=3) +
  geom_line(aes(tp, betat),
            data_frame(tp=1:8, betat=betat),
            size=2, lty=3, col="red")

ggplot(agretp) + geom_density(aes(omega, color=isSig2))


table(agretp_sig$ag, agretp_sig$re) / table(agretp$ag, agretp$re)
##                           aRhIgG.PE.high aRhIgG.PE.low   C1q R2A.2 R2A.3 R2A.4.high R2A.4.low R3A.1 R3A.3
## C1.Ak                              0.000         0.000 0.250       0.125      0.000     0.125       0.125
## C1.TR                              0.125         0.125 0.250       0.125      0.000     0.125       0.000
## G119                               0.125         0.000 0.125 0.000 0.000      0.000     0.000 0.125 0.250
## G145.146                                         0.000 0.125 0.125 0.000      0.125     0.000 0.000 0.000
## G49                                0.000         0.000 0.125       0.125      0.000     0.000 0.000 0.125
## G73                                0.000         0.000 0.375 0.000 0.125      0.000     0.000       0.000
## J08.V1V2.E660.084.AVI.His          0.125         0.000 0.125       0.750      0.125     0.000 0.000 0.625
## J08.V1V2.E660.2A5.AVI.His          0.250         0.125 0.000 0.125 0.750      0.125     0.125 0.000 0.625
## J08.V1V2.mac239.AVI.His            0.250         0.250 0.125 0.000 0.625      0.000     0.000 0.250 0.750
## SIV.1A11.gp140                     0.250         0.125 0.000       0.000      0.125           0.250 0.125
## SIV.E543.gp140                     0.875         0.375 0.000       0.500      0.125     0.000 0.500 0.750
## SIVcpz.EK505.gp120                 0.125         0.250             0.000      0.000           0.000 0.000
## SIVmac239.gp120                    0.125         0.000 0.000       0.750      0.250           0.750 0.750
## SIVmac239.gp130                    0.750         0.000 0.000       0.750      0.000           0.750 0.750
## SIVmac239.gp140                    0.875         0.875 0.250 0.000 0.875      0.750     0.375 1.000 1.000
## SIVmac239.gp140.AVI.His            0.000         0.125 0.125       0.000      0.000     0.125       0.125
## SIVmac251.BK.PR55                  0.000         0.000 0.000       0.000      0.125           0.125 0.000
## SIVsm.E660.2A5                     0.125         0.500 0.125       0.750      0.000           0.500 0.500
## SIVsm.E660.84                      0.125         0.125 0.125       0.125      0.000           0.375 0.250
## SIVsmH4.p55.Gag                    0.125         0.125 0.125       0.000      0.000     0.000 0.125 0.250
## V1a                                0.000         0.000 0.125       0.375      0.250     0.000 0.000 0.375

table(agretp_sig$ag, agretp_sig$tp) / table(agretp$ag, agretp$tp)
##                                   1         2         3         4         5         6         7         8
## C1.Ak                     0.2857143 0.1428571 0.1428571 0.0000000 0.1428571 0.0000000 0.0000000 0.0000000
## C1.TR                     0.2857143 0.0000000 0.0000000 0.2857143 0.2857143 0.0000000 0.0000000 0.0000000
## G119                      0.1111111 0.2222222 0.0000000 0.0000000 0.2222222 0.0000000 0.0000000 0.0000000
## G145.146                  0.1250000 0.0000000 0.0000000 0.1250000 0.0000000 0.0000000 0.0000000 0.1250000
## G49                       0.3750000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## G73                       0.1250000 0.0000000 0.1250000 0.1250000 0.1250000 0.0000000 0.0000000 0.0000000
## J08.V1V2.E660.084.AVI.His 0.1250000 0.3750000 0.2500000 0.0000000 0.3750000 0.1250000 0.2500000 0.2500000
## J08.V1V2.E660.2A5.AVI.His 0.1111111 0.2222222 0.3333333 0.1111111 0.2222222 0.2222222 0.4444444 0.2222222
## J08.V1V2.mac239.AVI.His   0.0000000 0.4444444 0.4444444 0.0000000 0.3333333 0.2222222 0.4444444 0.1111111
## SIV.1A11.gp140            0.0000000 0.0000000 0.1428571 0.0000000 0.2857143 0.1428571 0.1428571 0.2857143
## SIV.E543.gp140            0.0000000 0.2500000 0.2500000 0.1250000 0.6250000 0.5000000 0.7500000 0.6250000
## SIVcpz.EK505.gp120        0.0000000 0.0000000 0.0000000 0.1666667 0.0000000 0.1666667 0.0000000 0.1666667
## SIVmac239.gp120           0.0000000 0.4285714 0.5714286 0.0000000 0.4285714 0.4285714 0.5714286 0.5714286
## SIVmac239.gp130           0.0000000 0.5714286 0.5714286 0.1428571 0.5714286 0.5714286 0.5714286 0.4285714
## SIVmac239.gp140           0.3333333 0.6666667 0.7777778 0.6666667 0.7777778 0.6666667 0.7777778 0.6666667
## SIVmac239.gp140.AVI.His   0.2857143 0.2857143 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
## SIVmac251.BK.PR55         0.0000000 0.0000000 0.0000000 0.1428571 0.1428571 0.0000000 0.0000000 0.0000000
## SIVsm.E660.2A5            0.1428571 0.0000000 0.2857143 0.0000000 0.5714286 0.4285714 0.7142857 0.7142857
## SIVsm.E660.84             0.0000000 0.1428571 0.0000000 0.0000000 0.1428571 0.1428571 0.2857143 0.5714286
## SIVsmH4.p55.Gag           0.0000000 0.0000000 0.1250000 0.0000000 0.2500000 0.1250000 0.0000000 0.2500000
## V1a                       0.1250000 0.2500000 0.3750000 0.0000000 0.1250000 0.0000000 0.0000000 0.2500000

table(agretp_sig$re, agretp_sig$tp) / table(agretp$re, agretp$tp)
##                         1          2          3          4          5          6          7          8
## aRhIgG.PE.high 0.05000000 0.20000000 0.30000000 0.25000000 0.25000000 0.25000000 0.35000000 0.05000000
## aRhIgG.PE.low  0.00000000 0.09523810 0.14285714 0.09523810 0.19047619 0.09523810 0.23809524 0.28571429
## C1q            0.45000000 0.15000000 0.15000000 0.05000000 0.05000000 0.00000000 0.00000000 0.10000000
## R2A.2          0.00000000 0.00000000 0.00000000 0.33333333 0.00000000 0.00000000 0.00000000 0.00000000
## R2A.3          0.09523810 0.33333333 0.38095238 0.04761905 0.57142857 0.38095238 0.38095238 0.38095238
## R2A.4.high     0.00000000 0.04761905 0.14285714 0.09523810 0.09523810 0.04761905 0.19047619 0.14285714
## R2A.4.low      0.28571429 0.00000000 0.00000000 0.00000000 0.07142857 0.00000000 0.07142857 0.07142857
## R3A.1          0.05882353 0.29411765 0.17647059 0.05882353 0.52941176 0.35294118 0.35294118 0.41176471
## R3A.3          0.09523810 0.47619048 0.42857143 0.04761905 0.47619048 0.33333333 0.38095238 0.57142857

## 161 per tp
table(agretp_sig$tp)
##  1  2  3  4  5  6  7  8 
## 19 32 35 15 44 29 39 40 

table(agretp_sig$tp) / table(agretp$tp)
##         1         2         3         4         5         6         7         8 
## 0.1180124 0.1987578 0.2173913 0.0931677 0.2732919 0.1801242 0.2422360 0.2484472

table(agretp$ag)
##                   C1.Ak                     C1.TR                      G119                  G145.146 
##                      56                        56                        72                        64 
##                     G49                       G73 J08.V1V2.E660.084.AVI.His J08.V1V2.E660.2A5.AVI.His 
##                      64                        64                        64                        72 
## J08.V1V2.mac239.AVI.His            SIV.1A11.gp140            SIV.E543.gp140        SIVcpz.EK505.gp120 
##                      72                        56                        64                        48 
##         SIVmac239.gp120           SIVmac239.gp130           SIVmac239.gp140   SIVmac239.gp140.AVI.His 
##                      56                        56                        72                        56 
##       SIVmac251.BK.PR55            SIVsm.E660.2A5             SIVsm.E660.84           SIVsmH4.p55.Gag 
##                      56                        56                        56                        64 
##                     V1a 
##                      64 

table(agretp_sig$ag) / table(agretp$ag)
##                   C1.Ak                     C1.TR                      G119                  G145.146 
##              0.08928571                0.10714286                0.06944444                0.04687500 
##                     G49                       G73 J08.V1V2.E660.084.AVI.His J08.V1V2.E660.2A5.AVI.His 
##              0.04687500                0.06250000                0.21875000                0.23611111 
## J08.V1V2.mac239.AVI.His            SIV.1A11.gp140            SIV.E543.gp140        SIVcpz.EK505.gp120 
##              0.25000000                0.12500000                0.39062500                0.06250000 
##         SIVmac239.gp120           SIVmac239.gp130           SIVmac239.gp140   SIVmac239.gp140.AVI.His 
##              0.37500000                0.42857143                0.66666667                0.07142857 
##       SIVmac251.BK.PR55            SIVsm.E660.2A5             SIVsm.E660.84           SIVsmH4.p55.Gag 
##              0.03571429                0.35714286                0.16071429                0.09375000 
##                     V1a 
##              0.14062500 

table(agretp$re)
## aRhIgG.PE.high  aRhIgG.PE.low            C1q          R2A.2          R2A.3     R2A.4.high      R2A.4.low 
##            160            168            160             48            168            168            112 
##          R3A.1          R3A.3 
##            136            168 

table(agretp_sig$re) / table(agretp$re)
## aRhIgG.PE.high  aRhIgG.PE.low            C1q          R2A.2          R2A.3     R2A.4.high      R2A.4.low 
##     0.21250000     0.14285714     0.11875000     0.04166667     0.32142857     0.09523810     0.06250000 
##          R3A.1          R3A.3 
##     0.27941176     0.35119048 

lmfit <- lm(value ~ factor(groupId) - 1, data=longData.adj)

agretp$meanlm <- summary(lmfit)$coef[,1]


ggplot(agretp) + geom_point(aes(mean, meanlm)) +
    geom_abline(slope=1, intercept=0)



######################################################################
tsData <- longData.adj %>% select(AnimalId:value) %>%
  spread(tp, value)
names(tsData)[5:12] <- str_c("tp", 1:8)

N3 <- nrow(tsData)
D3 <- 8
K3 <- 7
y3 <- tsData[,5:12]

model_dat3 <- list(N = N3,
                   D = D3,
                   K = K3,
                   y = y3)

setwd("~/Projects/VRC332/Code/fh-vrc332/Modeling/Stan/")
fit3 <- stan("model-ts_clust.stan", data=model_dat3, iter=1000, chains=4)
## saveRDS(fit3, "model-ts_clust.rds")

fit3_summary <- as.data.frame(summary(fit3)$summary)
fit3_summary$var <- rownames(fit3_summary)

## write.table(fit3_summary, "fit3_summary.txt", quote=FALSE, row.names=FALSE)
## fit3_summary <- read.table("fit3_summary.txt", header=TRUE)

mus <- fit3_summary %>% filter(str_sub(var, 1, 3)=="mu[")

mumatrix <- matrix(nrow=7, ncol=8, byrow=TRUE, data=mus$mean)

mudf <- as.data.frame(mumatrix)
mudf$clust <- 1:7

mudf <- mudf %>% gather(clust)
names(mudf) <- c("clust", "tp", "value")
mudf <- mudf %>% mutate(tp=str_sub(tp, 2, 2))
mudf$tp <- as.numeric(mudf$tp)

ggplot(mudf) + geom_line(aes(tp, value, color=as.factor(clust)))


softzs <- fit3_summary %>% filter(str_sub(var, 1, 6) == "soft_z")
softzmatrix <- matrix(nrow=11431, ncol=7, byrow=TRUE, data=softzs$mean)

cluster_assign <- apply(softzmatrix, 1, clust_assign)




######################################################################
tsData <- longData.adj %>% select(AnimalId:value) %>%
  spread(tp, value)
names(tsData)[5:12] <- str_c("tp", 1:8)

agre <- tsData %>% select(ag, re) %>% distinct()
agre$groupId <- 1:nrow(agre)

tsData <- tsData %>%
    left_join(agre, by=c("ag" = "ag",
                         "re" = "re"))

N_subgroups4 <- length(unique(tsData$groupId))
obs_to_subgroup4 <- tsData$groupId

N4 <- nrow(tsData)
D4 <- 8
K4 <- 7
y4 <- tsData[,5:12]

model_dat4 <- list(N = N4,
                   D = D4,
                   K = K4,
                   N_subgroups = N_subgroups4,
                   obs_to_subgroup = obs_to_subgroup4,
                   y = y4)

setwd("~/Projects/VRC332/Code/fh-vrc332/Modeling/Stan/")
fit4 <- stan("model-ts_clust_group.stan", data=model_dat4, iter=1000, chains=4)

saveRDS(fit4, "model-ts_clust.rds")

fit4_summary <- as.data.frame(summary(fit4)$summary)
fit4_summary$var <- rownames(fit4_summary)

## write.table(fit4_summary, "fit4_summary.txt", quote=FALSE, row.names=FALSE)
## fit4_summary <- read.table("fit4_summary.txt", header=TRUE)

mus <- fit4_summary %>% filter(str_sub(var, 1, 3)=="mu[")

mumatrix <- matrix(nrow=7, ncol=8, byrow=TRUE, data=mus$mean)

mudf <- as.data.frame(mumatrix)
mudf$clust <- 1:7

mudf <- mudf %>% gather(clust)
names(mudf) <- c("clust", "tp", "value")
mudf <- mudf %>% mutate(tp=str_sub(tp, 2, 2))
mudf$tp <- as.numeric(mudf$tp)

ggplot(mudf %>% filter(clust %in% c(1,4,7))) + geom_line(aes(tp, value, color=as.factor(clust)))



softzs <- fit4_summary %>% filter(str_sub(var, 1, 6) == "soft_z")
softzmatrix <- matrix(nrow=161, ncol=7, byrow=TRUE, data=softzs$mean)

cluster_assignments <- apply(softzmatrix, 1, clust_assign)

agre$cluster <- cluster_assignments

table(agre$ag, agre$cluster)
##                           1 4 7
## C1.Ak                     2 0 5
## C1.TR                     2 0 5
## G119                      2 0 7
## G145.146                  2 0 6
## G49                       2 0 6
## G73                       2 0 6
## J08.V1V2.E660.084.AVI.His 0 4 4
## J08.V1V2.E660.2A5.AVI.His 0 4 5
## J08.V1V2.mac239.AVI.His   0 3 6
## SIV.1A11.gp140            0 3 4
## SIV.E543.gp140            0 6 2
## SIVcpz.EK505.gp120        0 0 6
## SIVmac239.gp120           0 5 2
## SIVmac239.gp130           0 5 2
## SIVmac239.gp140           0 7 2
## SIVmac239.gp140.AVI.His   2 0 5
## SIVmac251.BK.PR55         1 0 6
## SIVsm.E660.2A5            0 5 2
## SIVsm.E660.84             0 5 2
## SIVsmH4.p55.Gag           0 3 5
## V1a                       0 2 6

table(agre$re, agre$cluster)
##                 1  4  7
## aRhIgG.PE.high  0  9 11
## aRhIgG.PE.low   0  8 13
## C1q             8  1 11
## R2A.2           0  0  6
## R2A.3           0 12  9
## R2A.4.high      7  2 12
## R2A.4.low       0  0 14
## R3A.1           0  8  9
## R3A.3           0 12  9




######################################################################
tsData <- longData.adj %>% select(AnimalId:value) %>%
  spread(tp, value)
names(tsData)[5:12] <- str_c("tp", 1:8)

agregrp <- tsData %>% select(ag, re, GroupNm) %>% distinct()
agregrp$groupId <- 1:nrow(agregrp)

tsData <- tsData %>%
    left_join(agregrp, by=c("ag" = "ag",
                            "re" = "re",
                            "GroupNm" = "GroupNm"))

N_subgroups5 <- length(unique(tsData$groupId))
obs_to_subgroup5 <- tsData$groupId

N5 <- nrow(tsData)
D5 <- 8
K5 <- 7
y5 <- tsData[,5:12]

model_dat5 <- list(N = N5,
                   D = D5,
                   K = K5,
                   N_subgroups = N_subgroups5,
                   obs_to_subgroup = obs_to_subgroup5,
                   y = y5)

setwd("~/Projects/VRC332/Code/fh-vrc332/Modeling/Stan/")
fit5 <- stan("model-ts_clust_group.stan", data=model_dat5, iter=1000, chains=4)

saveRDS(fit5, "model-ts_clust_agregrp.rds")

fit5_summary <- as.data.frame(summary(fit5)$summary)
fit5_summary$var <- rownames(fit5_summary)

## write.table(fit5_summary, "fit5_summary.txt", quote=FALSE, row.names=FALSE)
## fit5_summary <- read.table("fit5_summary.txt", header=TRUE)

mus <- fit5_summary %>% filter(str_sub(var, 1, 3)=="mu[")

mumatrix <- matrix(nrow=7, ncol=8, byrow=TRUE, data=mus$mean)

mudf <- as.data.frame(mumatrix)
mudf$clust <- 1:7

mudf <- mudf %>% gather(clust)
names(mudf) <- c("clust", "tp", "value")
mudf <- mudf %>% mutate(tp=str_sub(tp, 2, 2))
mudf$tp <- as.numeric(mudf$tp)

ggplot(mudf %>% filter(clust %in% c(3,5,6,7))) +
    geom_line(aes(tp, value, color=as.factor(clust)))

softzs <- fit5_summary %>% filter(str_sub(var, 1, 6) == "soft_z")
softzmatrix <- matrix(nrow=805, ncol=7, byrow=TRUE, data=softzs$mean)

cluster_assignments <- apply(softzmatrix, 1, clust_assign)

agregrp$cluster <- cluster_assignments

## Cluster 7: low
## Cluster 5: low
## Cluster 3: mid
## Cluster 6: high

table(agregrp$GroupNm, agregrp$cluster)
##                 3  5  6  7
## Control        40 58  0 63
## SIV_Env        36 48 31 46
## SIV_Gag        37 58  2 64
## SIV_Mosaic_Env 54 52  6 49
## x_PARI         41 53 25 42


table(agregrp$ag, agregrp$cluster)
##                            3  5  6  7
## C1.Ak                      0 12  0 23
## C1.TR                      0  9  0 26
## G119                       1 14  0 30
## G145.146                   0 12  0 28
## G49                        0 16  0 24
## G73                        0 14  0 26
## J08.V1V2.E660.084.AVI.His 15 15  4  6
## J08.V1V2.E660.2A5.AVI.His 16 16  4  9
## J08.V1V2.mac239.AVI.His   12 15  6 12
## SIV.1A11.gp140            17 17  0  1
## SIV.E543.gp140            26  7  6  1
## SIVcpz.EK505.gp120         0 16  0 14
## SIVmac239.gp120           17  9  7  2
## SIVmac239.gp130           16  9  8  2
## SIVmac239.gp140           19  4 17  5
## SIVmac239.gp140.AVI.His    2 18  0 15
## SIVmac251.BK.PR55          0 15  0 20
## SIVsm.E660.2A5            22  7  6  0
## SIVsm.E660.84             21  9  0  5
## SIVsmH4.p55.Gag           17 17  2  4
## V1a                        7 18  4 11

table(agregrp$re, agregrp$cluster)
##                 3  5  6  7
## aRhIgG.PE.high 29 31 12 28
## aRhIgG.PE.low  35 51  3 16
## C1q             6 40  0 54
## R2A.2           0  0  0 30
## R2A.3          42 43 18  2
## R2A.4.high     25 30  2 48
## R2A.4.low       3 12  0 55
## R3A.1          30 23 10 22
## R3A.3          38 39 19  9

require(RColorBrewer)
cutColors <- brewer.pal(4, "Set1")

clustTrans <- data_frame(cluster=c(3,5:7), clusterName=c("Mid", "Low2", "High", "Low"))

agregrp <- agregrp %>% left_join(clustTrans)

multiPie.grp.ag <- ggplot(agregrp) + geom_bar(aes(x=factor(1), fill=factor(clusterName)), width=1) +
    coord_polar(theta="y") + theme_void(base_size=8) + facet_wrap(~GroupNm+ag, nrow=5) +
    scale_fill_manual(values=cutColors,
                      limits=levels(factor(agregrp$clusterName)))





######################################################################
tsData <- longData.adj %>% select(AnimalId:value) %>%
  spread(tp, value)
names(tsData)[5:12] <- str_c("tp", 1:8)
tsData <- tsData %>% select(AnimalId:tp6)

agregrp <- tsData %>% select(ag, re, GroupNm) %>% distinct()
agregrp$groupId <- 1:nrow(agregrp)

tsData <- tsData %>%
    left_join(agregrp, by=c("ag" = "ag",
                            "re" = "re",
                            "GroupNm" = "GroupNm"))

N_subgroups6 <- length(unique(tsData$groupId))
obs_to_subgroup6 <- tsData$groupId

N6 <- nrow(tsData)
D6 <- 6
K6 <- 10
y6 <- tsData[,5:10]

model_dat6 <- list(N = N6,
                   D = D6,
                   K = K6,
                   N_subgroups = N_subgroups6,
                   obs_to_subgroup = obs_to_subgroup6,
                   y = y6)

setwd("~/Projects/VRC332/Code/fh-vrc332/Modeling/Stan/")
fit6 <- stan("model-ts_clust_group.stan", data=model_dat6, iter=1000, chains=4)

## saveRDS(fit6, "model-ts1-6_clust_agregrp.rds")
## fit6 <- readRDS("model-ts1-6_clust_agregrp.rds")

fit6_summary <- as.data.frame(summary(fit6)$summary)
fit6_summary$var <- rownames(fit6_summary)

## write.table(fit6_summary, "fit6_summary.txt", quote=FALSE, row.names=FALSE)
## fit6_summary <- read.table("fit6_summary.txt", header=TRUE)

mus <- fit6_summary %>% filter(str_sub(var, 1, 3)=="mu[")

mumatrix <- matrix(nrow=K6, ncol=D6, byrow=TRUE, data=mus$mean)

mudf <- as.data.frame(mumatrix)
mudf$clust <- 1:K6

mudf <- mudf %>% gather(clust)
names(mudf) <- c("clust", "tp", "value")
mudf <- mudf %>% mutate(tp=str_sub(tp, 2, 2))
mudf$tp <- as.numeric(mudf$tp)


softzs <- fit6_summary %>% filter(str_sub(var, 1, 6) == "soft_z")
softzmatrix <- matrix(nrow=N_subgroups6, ncol=K6, byrow=TRUE, data=softzs$mean)

cluster_assignments <- apply(softzmatrix, 1, clust_assign)

agregrp$cluster <- cluster_assignments

clustColors <- brewer.pal(5, "Set1")

clustTrans <- data_frame(cluster=c(1,5,8:10),
                         clusterName=c("LowMid2", "High", "LowMid1", "Mid", "Low"))

agregrp <- agregrp %>% left_join(clustTrans)

multiPie.grp.ag <- ggplot(agregrp) + geom_bar(aes(x=factor(1), fill=factor(clusterName)), width=1) +
    coord_polar(theta="y") + theme_void(base_size=8) + facet_wrap(~GroupNm+ag, nrow=5) +
    scale_fill_manual(values=clustColors,
                      limits=levels(factor(agregrp$clusterName)))

mudf <- mudf %>% left_join(clustTrans, by=c("clust" = "cluster"))

ggplot(mudf %>% filter(clust %in% c(1,5,8:10))) +
    geom_line(aes(tp, value, color=clusterName)) +
    scale_color_manual(values=clustColors,
                      limits=levels(factor(agregrp$clusterName)))




######################################################################
tsData <- longData.adj %>% select(AnimalId:value) %>%
  spread(tp, value)
names(tsData)[5:12] <- str_c("tp", 1:8)
tsData <- tsData %>% select(AnimalId:tp6)

agregrp <- tsData %>% select(ag, re, GroupNm) %>% distinct()
agregrp$groupId <- 1:nrow(agregrp)

tsData <- tsData %>%
    left_join(agregrp, by=c("ag" = "ag",
                            "re" = "re",
                            "GroupNm" = "GroupNm"))

N_subgroups7 <- length(unique(tsData$groupId))
obs_to_subgroup7 <- tsData$groupId

N7 <- nrow(tsData)
D7 <- 6
K7 <- 10
y7 <- tsData[,5:10]

model_dat7 <- list(N = N7,
                   D = D7,
                   K = K7,
                   N_subgroups = N_subgroups7,
                   obs_to_subgroup = obs_to_subgroup7,
                   y = y7)

setwd("~/Projects/VRC332/Code/fh-vrc332/Modeling/Stan/")
fit7 <- stan("model-ts_clust_group.stan", data=model_dat7, iter=1000, chains=4)

## saveRDS(fit7, "model-ts1-6_clust_agregrp_b.rds")
## fit7 <- readRDS("model-ts1-6_clust_agregrp_b.rds")

fit7_summary <- as.data.frame(summary(fit7)$summary)
fit7_summary$var <- rownames(fit7_summary)

## write.table(fit7_summary, "fit7_summary.txt", quote=FALSE, row.names=FALSE)
## fit7_summary <- read.table("fit7_summary.txt", header=TRUE)

mus <- fit7_summary %>% filter(str_sub(var, 1, 3)=="mu[")

mumatrix <- matrix(nrow=K7, ncol=D7, byrow=TRUE, data=mus$mean)

mudf7 <- as.data.frame(mumatrix)
mudf7$clust <- 1:K7

mudf7 <- mudf7 %>% gather(clust)
names(mudf7) <- c("clust", "tp", "value")
mudf7 <- mudf7 %>% mutate(tp=str_sub(tp, 2, 2))
mudf7$tp <- as.numeric(mudf7$tp)


softzs <- fit7_summary %>% filter(str_sub(var, 1, 6) == "soft_z")
softzmatrix <- matrix(nrow=N_subgroups7, ncol=K7, byrow=TRUE, data=softzs$mean)

cluster_assignments7 <- apply(softzmatrix, 1, clust_assign)

agregrp$cluster <- cluster_assignments7

clustColors <- brewer.pal(5, "Set1")

clustTrans7 <- data_frame(cluster=c(3, 9, 6, 1, 8),
                         clusterName=c("Low1", "Low2", "LowMid", "Mid", "High"))

agregrp <- agregrp %>% left_join(clustTrans7)

multiPie.grp.ag <- ggplot(agregrp) + geom_bar(aes(x=factor(1), fill=factor(clusterName)), width=1) +
    coord_polar(theta="y") + theme_void(base_size=8) + facet_wrap(~GroupNm+ag, nrow=5) +
    scale_fill_manual(values=clustColors,
                      limits=levels(factor(agregrp$clusterName)))

mudf7 <- mudf7 %>% left_join(clustTrans7, by=c("clust" = "cluster"))

ggplot(mudf7 %>% filter(clust %in% c(1,3,6,8,9))) +
    geom_line(aes(tp, value, color=factor(clusterName))) +
    scale_color_manual(values=clustColors,
                       limits=levels(factor(agregrp$clusterName)))
