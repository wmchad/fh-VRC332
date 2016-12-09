require(caret)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "GetRmse.r"))
source(file.path(fnFolder, "Prediction/CoxPhHelper.r"))


RunRandomPartitionPredictions_CoxPh <- function(data,
                                                pTrain        = 0.75,
                                                nRand         = 50,
                                                coxIter       = 100,
                                                outdir        = "",
                                                rmseFile      = "rmseDefault.txt",
                                                fitsFile      = "fitsDefault.rdata",
                                                verbose       = TRUE,
                                                progressEvery = 5) {
  vtimestamp <- function() {}
  vprint <- function(...) {}
  if ( verbose ) {
    vtimestamp <- function(){timestamp()}
    vprint <- function(..., sep=""){print(paste(..., sep=sep))}
  }

  curdir <- getwd()
  compRmse <- matrix(nrow=nRand, ncol=4)
  colnames(compRmse) <- c("baseline", "group", "predicted", "predicted_shrunk")
  fits <- NULL

  vtimestamp()
  vprint("Starting...")
  for ( j in 1:nRand ) {

      inTrain = createDataPartition(data$groups, p=0.75, list=FALSE)
      train.x <- data$x[inTrain,]
      train.y <- data$y[inTrain]
      train.surv <- data$surv[inTrain]
      test.x <- data$x[-inTrain,]
      test.y <- data$y[-inTrain]
      test.surv <- data$surv[-inTrain]

      res <- cv.glmnet(x=as.matrix(train.x), y=train.surv, family="cox")
      coef <- as.data.frame(summary(predict(res, type="coef", s=res$lambda.min)))
      vars <- (predSummary %>% filter(shortVarName %in% names(data$x)[coef$i]))$shortVarName
      preds <- rep(mean(train.y), length(test.y))
      preds2 <- rep(mean(train.y), length(test.y))
      if ( nrow(coef) > 0 ) {
          preds <- CoxPhPredictions(train.x[,coef$i], train.surv,
                                    test.x[,coef$i], coef$x, coxIter)
          preds2 <- CoxPhPredictions(train.x[,coef$i], train.surv,
                                     test.x[,coef$i], coef$x, 0)
      }

      fits[[j]] <- list(fit=res, vars=vars, predicted=preds,
                        predicted.shrunk=preds2, actual=test.y,
                        testAnimals=data$animalIds[-inTrain])
      compRmse[j,] <- c(rmse(mean(train.y), test.y),
                        GroupRmse(train.x, train.y, test.x, test.y),
                        rmse(preds, test.y),
                        rmse(preds2, test.y))

      if ( j %% progressEvery == 0 ) {
          vtimestamp()
          vprint(paste(j, "iterations finished"))
          setwd(outdir)
          write.table(compRmse, rmseFile,
                      quote=FALSE, row.names=FALSE, col.names=FALSE)
          save(fits, file=fitsFile)
      }
  }
  setwd(outdir)
  write.table(compRmse, rmseFile,
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  save(fits, file=fitsFile)
  setwd(curdir)
  list(compRmse=compRmse,
       fits=fits)
}
