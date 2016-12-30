require(rstan)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "Modeling/LinearModels-Support.r"))

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

fit2_summary <- as.data.frame(summary(fit2)$summary)
fit2_summary$var <- rownames(fit2_summary)
agretp$mean2 <- (fit2_summary %>% filter(str_sub(var, 1, 5) == "beta["))$mean
agretp$omega <- (fit2_summary %>% filter(str_sub(var, 1, 5) == "omega"))$mean

ggplot(agretp) + geom_point(aes(mean2, omega, color=factor(tp)))

ggplot(agretp) + geom_point(aes(mean, mean2, color=factor(tp), alpha=omega))

ggplot(agretp) + geom_density(aes(omega, color=ag))
ggplot(agretp) + geom_density(aes(omega, color=re))


sigBetas <- (fit2_summary %>%
             filter(sign(`2.5%`) == sign(`97.5%`)))$var %>%
                                                  str_sub(6, -2)
sigBetas <- as.numeric(sigBetas)
sigBetas <- sigBetas[!is.na(sigBetas)]
agretp$isSig2 <- FALSE
agretp$isSig2[sigBetas] <- TRUE

ggplot(agretp) + geom_point(aes(mean, mean2, color=factor(tp), shape=isSig2))

ggplot(agretp) + geom_density(aes(mean2, color=isSig2))


agHighlight <- "SIVmac239.gp140"

ggplot(agretp %>% group_by(ag, re)) +
    geom_line(aes(tp, mean2, group=paste(ag, re)), alpha = 0.1) +
    geom_line(aes(tp, mean2, group=paste(ag, re), color=re),
              agretp %>% filter(ag==agHighlight) %>% group_by(ag, re),
              size=1.2) +
    geom_line(aes(tp, mean, group=paste(ag, re), color=re),
              agretp %>% filter(ag==agHighlight) %>% group_by(ag, re),
              size=1, lty=2)




lmfit <- lm(value ~ factor(groupId) - 1, data=longData.adj)

agretp$meanlm <- summary(lmfit)$coef[,1]


ggplot(agretp) + geom_point(aes(mean, meanlm)) +
    geom_abline(slope=1, intercept=0)
