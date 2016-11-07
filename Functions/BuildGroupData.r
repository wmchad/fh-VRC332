require(dplyr)

BuildGroupData <- function(vlData) {
  x <- model.matrix(~ GroupNm -1, data=vlData)
  colnames(x) <- substr(colnames(x), 8, 50)
  data.frame(x) %>% select(SIV_Mosaic_Env, SIV_Env, SIV_Gag, x_PARI)
}
