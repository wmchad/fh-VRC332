ReplaceMissingWithMean <- function(x) {
  for ( i in 1:ncol(x) ) {
    missing <- is.na(x[,i]) | x[,i]<0
    x[missing,i] <- mean(x[!missing,i])
  }
  x
}
