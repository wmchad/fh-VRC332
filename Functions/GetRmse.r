GetRmse <- function(train.x, train.y, test.x, test.y,
                    method="glmnet",
                    ctrl=trainControl(method="repeatedcv", repeats=5)) {
  res <- train(x=train.x, y=train.y, method=method, trControl=ctrl)
  pred <- predict(res, test.x)
  list(fit=res, predicted=pred, rmse=RMSE(pred, test.y))
}

GroupRmse <- function(train.x, train.y, test.x, test.y) {
  preds <- rep(0, 5)
  for ( i in 1:4 ) {
    preds[i] <- mean(train.y[train.x[,i]==1])
  }
  preds[5] <- mean(train.y[apply(train.x[,1:4], 1, sum)==0])
  testPreds <- apply(test.x[,1:4], 1, function(x) { sum(x*preds[1:4]) })
  testPreds[apply(test.x[,1:4], 1, sum)==0] <- preds[5]
  RMSE(testPreds, test.y)
}

