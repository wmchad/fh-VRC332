require(caret)

fnFolder <- "~/Projects/VRC332/Code/fh-vrc332/Functions"
source(file.path(fnFolder, "GetRmse.r"))

RunRandomPartitionPredictions_Classification<- function(data, pTrain=0.75,
                                          nRand=50, method="rf",
                                          ctrl=trainControl(method="repeatedcv", repeats=5),
                                          outdir="",
                                          resultsFile="resultsDefault.txt",
                                          verbose=TRUE, progressEvery=5, ...) {
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

  results <- data.frame(iteration=numeric(0),
                        animalId=character(0),
                        pred=character(0),
                        actual=character(0),
                        stringsAsFactors=FALSE)
  accuracy <- data.frame(train=numeric(0), test=numeric(0))
  missedTrain <- NULL
  missedTest <- NULL
  fits <- NULL

  vtimestamp()
  vprint("Starting...")
  for ( i in 1:nRand ) {
      inTrain = createDataPartition(data$groups, p=0.75, list=FALSE)
      train.x <- data$x[inTrain,]
      train.y <- data$y[inTrain]
      test.x <- data$x[-inTrain,]
      test.y <- data$y[-inTrain]
      train.ids <- data$animalIds[inTrain]
      test.ids <- data$animalIds[-inTrain]

      fit <- train(train.x, as.factor(train.y), method=method, trControl=ctrl, ...)
      predTrain <- predict(fit, train.x)
      predTest <- predict(fit, test.x)
      results <- rbind(results, data.frame(iteration=i,
                                           animalId=as.character(test.ids),
                                           pred=as.character(predTest),
                                           actual=as.character(test.y)))
      accuracy[i,] <- c(sum(predTrain==train.y)/length(train.y),
                        sum(predTest==test.y)/length(test.y))
      missedTrain <- c(missedTrain, train.ids[predTrain!=train.y])
      missedTest <- c(missedTest, test.ids[predTest!=test.y])
      fits[[i]] <- fit

      if ( i %% progressEvery == 0 ) {
          vtimestamp()
          vprint(paste(i, "iterations finished"))
          setwd(outdir)
          multRes <- list(testResults=results,
                          accuracy=accuracy,
                          missedTrain=missedTrain,
                          missedTest=missedTest,
                          fits=fits)
          save(multRes, file=resultsFile)
      }
  }
  setwd(outdir)
  multRes <- list(testResults=results,
                  accuracy=accuracy,
                  missedTrain=missedTrain,
                  missedTest=missedTest,
                  fits=fits)
  save(multRes, file=resultsFile)
  setwd(curdir)
  multRes
}

GetRandomPartitionPredictions_Classification<- function(outdir="",
                                                        resultsFile="resultsDefault.txt") {
  load(file.path(outdir, resultsFile))
  multRes
}
