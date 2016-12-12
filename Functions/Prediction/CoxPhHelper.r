require(glmnet)
require(glmnetUtils)

CoxPhPredictions <- function(train.x, y, test.x, initx=NULL, coxIter=0) {
    train.x <- as.matrix(train.x)
    test.x <- as.matrix(test.x)
    ## Get the cox model
    ##   coxIter = 0 means use the initial values provided (e.g. shrunk from glmnet)
    ##   coxIter > 0 allows the coxph function to fit the model without shrinkage
    coxModel <- coxph(y~train.x, init=initx, iter=coxIter)
    coxFit <- survfit(coxModel, newData=train.x)
    ## Build the baseline hazard
    s0 <- exp(-coxFit$cumhaz)
    ## Get the model coefficients
    coefs <- coxModel$coefficients
    ## Need to get the constant term separately
    const <- attr(predict(coxModel, type="terms"), "constant")
    ## Compute the expected survival time for each item in the test set
    sapply(1:nrow(test.x), function(i) {
        sum(s0^exp(sum(test.x[i,] * coefs) - const))
    })
}
