require(caret)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "GetRmse.r"))

RunRandomPartitionPredictions <- function(data, pTrain=0.75,
                                          nRand=50, method="glmnet",
                                          ctrl=trainControl(method="repeatedcv", repeats=5),
                                          outdir="",
                                          rmseFile="rmseDefault.txt",
                                          fitsFile="fitsDefault.rdata",
                                          verbose=TRUE, progressEvery=5) {
  vtimestamp <- function() {}
  vprint <- function(...) {}
  if ( verbose ) {
    vtimestamp <- function(){timestamp()}
    vprint <- function(..., sep=""){print(paste(..., sep=sep))}
  }

  curdir <- getwd()
  compRmse <- matrix(nrow=nRand, ncol=3)
  colnames(compRmse) <- c("baseline", "group", "predicted")
  fits <- NULL

  vtimestamp()
  vprint("Starting...")
  for ( i in 1:nRand ) {
    inTrain = createDataPartition(data$groups, p=0.75, list=FALSE)
    train.x <- data$x[inTrain,]
    train.y <- data$y[inTrain]
    test.x <- data$x[-inTrain,]
    test.y <- data$y[-inTrain]

    fit <- GetRmse(train.x, train.y, test.x, test.y, method, ctrl)
    fits[[i]] <- list(fit=fit$fit, predicted=fit$predicted, actual=test.y,
                      testAnimals=data$animalIds[-inTrain])
    compRmse[i,] <- c(rmse(mean(train.y), test.y),
                      GroupRmse(train.x, train.y, test.x, test.y),
                      fit$rmse)
    if ( i %% progressEvery == 0 ) {
      vtimestamp()
      vprint(paste(i, "iterations finished"))
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

GetRandomPartitionPredictions <- function(outdir="",
                                          rmseFile="rmseDefault.txt",
                                          fitsFile="fitsDefault.rdata") {
  curdir <- getwd()
  setwd(outdir)
  compRmse <- read.table(rmseFile)
  names(compRmse) <- c("baseline", "group", "predicted")
  load(fitsFile)
  setwd(curdir)
  list(compRmse=compRmse,
       fits=fits)
}
