require(glmnet)
require(glmnetUtils)

CoxPhPredictions <- function(train.x, y, test.x, initx=NULL, coxIter=0) {
    train.x <- as.matrix(train.x)
    test.x <- as.matrix(test.x)
    coxModel <- coxph(y~train.x, init=initx, iter=coxIter)
    coxFit <- survfit(coxModel, newData=train.x)
    s0 <- exp(-coxFit$cumhaz)
    coefs <- coxModel$coefficients
    const <- attr(predict(coxModel, type="terms"), "constant")
    sapply(1:nrow(test.x), function(i) {
        sum(s0^exp(sum(test.x[i,] * coefs) - const))
    })
}
