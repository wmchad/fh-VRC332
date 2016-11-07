library(VRC332)

## Load Data
data("FcCoreAssay")

## Remove NAs
FcCoreAssay <- FcCoreAssay[!is.na(FcCoreAssay$readout),]
FcCoreAssay$logreadout <- log(FcCoreAssay$readout)

FcLmData <- FcCoreAssay[,c("AnimalID", "Timepoint", "readout",
                           "DetectionReagent", "Antigen",
                           "GroupNm")]
FcLmData$logreadout <- log(FcLmData$readout)

summary(lm(logreadout ~ GroupNm, data=FcLmData))
##                       Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            6.18161    0.01381 447.727  < 2e-16 ***
## GroupNmSIV_Env         0.92271    0.02116  43.603  < 2e-16 ***
## GroupNmSIV_Gag         0.13822    0.01952   7.082 1.43e-12 ***
## GroupNmSIV_Mosaic_Env  0.58660    0.01945  30.154  < 2e-16 ***
## GroupNmx_PARI          0.68150    0.02085  32.680  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 2.216 on 116892 degrees of freedom
## Multiple R-squared:  0.02299,	Adjusted R-squared:  0.02296 
## F-statistic: 687.7 on 4 and 116892 DF,  p-value: < 2.2e-16

summary(lm(logreadout ~ Antigen, data=FcLmData))
##                                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                       6.87565    0.02945 233.439  < 2e-16 ***
## AntigenC1.TR                      0.52594    0.04165  12.626  < 2e-16 ***
## AntigenG119                      -0.44326    0.03952 -11.217  < 2e-16 ***
## AntigenG145.146                  -0.68539    0.04048 -16.932  < 2e-16 ***
## AntigenG49                        0.02387    0.04048   0.590  0.55534    
## AntigenG73                       -0.46796    0.04048 -11.560  < 2e-16 ***
## AntigenJ08.V1V2.E660.084.AVI.His -0.13103    0.04048  -3.237  0.00121 ** 
## AntigenJ08.V1V2.E660.2A5.AVI.His -0.29703    0.03952  -7.517 5.65e-14 ***
## AntigenJ08.V1V2.mac239.AVI.His   -0.32774    0.03952  -8.294  < 2e-16 ***
## AntigenSIV.1A11.gp140            -0.16674    0.04312  -3.867  0.00011 ***
## AntigenSIV.E543.gp140            -0.45924    0.04165 -11.025  < 2e-16 ***
## AntigenSIVcpz.EK505.gp120        -1.71994    0.04499 -38.226  < 2e-16 ***
## AntigenSIVmac239.gp120           -0.09673    0.04312  -2.243  0.02487 *  
## AntigenSIVmac239.gp130           -0.40292    0.04312  -9.345  < 2e-16 ***
## AntigenSIVmac239.gp140            0.54361    0.04048  13.429  < 2e-16 ***
## AntigenSIVmac239.gp140.AVI.His   -0.11327    0.04165  -2.719  0.00654 ** 
## AntigenSIVmac251.BK.PR55         -1.47421    0.04312 -34.191  < 2e-16 ***
## AntigenSIVsm.E660.2A5            -0.51562    0.04312 -11.959  < 2e-16 ***
## AntigenSIVsm.E660.84             -0.95560    0.04312 -22.163  < 2e-16 ***
## AntigenSIVsmH4.p55.Gag            0.05765    0.04165   1.384  0.16633    
## AntigenV1a                        0.79543    0.04048  19.650  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 2.172 on 116876 degrees of freedom
## Multiple R-squared:  0.0616,	Adjusted R-squared:  0.06144 
## F-statistic: 383.6 on 20 and 116876 DF,  p-value: < 2.2e-16

summary(lm(logreadout ~ DetectionReagent, data=FcLmData))
##                               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                    7.38599    0.01710 431.832  < 2e-16 ***
## DetectionReagentaRhIgG.PE.low -1.38966    0.02390 -58.148  < 2e-16 ***
## DetectionReagentC1q           -1.68865    0.02419 -69.812  < 2e-16 ***
## DetectionReagentMBL           -2.07463    0.02871 -72.254  < 2e-16 ***
## DetectionReagentR2A.2         -2.40144    0.03560 -67.448  < 2e-16 ***
## DetectionReagentR2A.3         -0.12713    0.02390  -5.319 1.04e-07 ***
## DetectionReagentR2A.4.high    -0.26396    0.02392 -11.033  < 2e-16 ***
## DetectionReagentR2A.4.low     -2.36929    0.02665 -88.889  < 2e-16 ***
## DetectionReagentR3A.1          1.05509    0.02523  41.814  < 2e-16 ***
## DetectionReagentR3A.3         -0.45636    0.02390 -19.096  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 1.995 on 116887 degrees of freedom
## Multiple R-squared:  0.2084,	Adjusted R-squared:  0.2083 
## F-statistic:  3419 on 9 and 116887 DF,  p-value: < 2.2e-16

summary(lm(logreadout ~ GroupNm + Antigen + DetectionReagent, data=FcLmData))
##                                  Estimate Std. Error  t value Pr(>|t|)    
## (Intercept)                       7.62477    0.03050  249.985  < 2e-16 ***
## GroupNmSIV_Env                    0.92188    0.01739   53.018  < 2e-16 ***
## GroupNmSIV_Gag                    0.13739    0.01604    8.567  < 2e-16 ***
## GroupNmSIV_Mosaic_Env             0.58577    0.01598   36.646  < 2e-16 ***
## GroupNmx_PARI                     0.68120    0.01713   39.755  < 2e-16 ***
## AntigenC1.TR                      0.52594    0.03492   15.061  < 2e-16 ***
## AntigenG119                      -0.52648    0.03333  -15.798  < 2e-16 ***
## AntigenG145.146                  -0.64690    0.03421  -18.908  < 2e-16 ***
## AntigenG49                       -0.24248    0.03400   -7.132 9.92e-13 ***
## AntigenG73                       -0.29408    0.03410   -8.623  < 2e-16 ***
## AntigenJ08.V1V2.E660.084.AVI.His -0.39738    0.03400  -11.689  < 2e-16 ***
## AntigenJ08.V1V2.E660.2A5.AVI.His -0.38026    0.03333  -11.410  < 2e-16 ***
## AntigenJ08.V1V2.mac239.AVI.His   -0.41096    0.03333  -12.332  < 2e-16 ***
## AntigenSIV.1A11.gp140            -0.94245    0.03643  -25.870  < 2e-16 ***
## AntigenSIV.E543.gp140            -0.93230    0.03510  -26.565  < 2e-16 ***
## AntigenSIVcpz.EK505.gp120        -2.72754    0.03809  -71.599  < 2e-16 ***
## AntigenSIVmac239.gp120           -0.87244    0.03643  -23.948  < 2e-16 ***
## AntigenSIVmac239.gp130           -1.17862    0.03643  -32.353  < 2e-16 ***
## AntigenSIVmac239.gp140            0.29702    0.03426    8.670  < 2e-16 ***
## AntigenSIVmac239.gp140.AVI.His   -0.11327    0.03492   -3.244  0.00118 ** 
## AntigenSIVmac251.BK.PR55         -2.24992    0.03643  -61.760  < 2e-16 ***
## AntigenSIVsm.E660.2A5            -1.29133    0.03643  -35.447  < 2e-16 ***
## AntigenSIVsm.E660.84             -1.73131    0.03643  -47.525  < 2e-16 ***
## AntigenSIVsmH4.p55.Gag           -0.41541    0.03510  -11.837  < 2e-16 ***
## AntigenV1a                        0.52908    0.03400   15.563  < 2e-16 ***
## DetectionReagentaRhIgG.PE.low    -1.39065    0.02185  -63.658  < 2e-16 ***
## DetectionReagentC1q              -1.79268    0.02216  -80.904  < 2e-16 ***
## DetectionReagentMBL              -2.56445    0.02685  -95.526  < 2e-16 ***
## DetectionReagentR2A.2            -2.74221    0.03373  -81.302  < 2e-16 ***
## DetectionReagentR2A.3            -0.12812    0.02185   -5.865 4.51e-09 ***
## DetectionReagentR2A.4.high       -0.26584    0.02187  -12.156  < 2e-16 ***
## DetectionReagentR2A.4.low        -2.82218    0.02475 -114.042  < 2e-16 ***
## DetectionReagentR3A.1             1.21795    0.02323   52.441  < 2e-16 ***
## DetectionReagentR3A.3            -0.45735    0.02185  -20.936  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 1.821 on 116863 degrees of freedom
## Multiple R-squared:  0.3405,	Adjusted R-squared:  0.3403 
## F-statistic:  1829 on 33 and 116863 DF,  p-value: < 2.2e-16

summary(lm(logreadout ~ GroupNm + Antigen * DetectionReagent, data=FcLmData))
##                                                                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                                                     6.580295   0.063555 103.537  < 2e-16 ***
## GroupNmSIV_Env                                                  0.921876   0.015667  58.842  < 2e-16 ***
## GroupNmSIV_Gag                                                  0.137389   0.014450   9.508  < 2e-16 ***
## GroupNmSIV_Mosaic_Env                                           0.585768   0.014402  40.672  < 2e-16 ***
## GroupNmx_PARI                                                   0.681199   0.015439  44.122  < 2e-16 ***
## AntigenC1.TR                                                    0.023702   0.088970   0.266 0.789929    
## AntigenG119                                                    -0.792286   0.088970  -8.905  < 2e-16 ***
## AntigenG145.146                                                -1.443867   0.088970 -16.229  < 2e-16 ***
## AntigenG49                                                     -0.610074   0.088970  -6.857 7.06e-12 ***
## AntigenG73                                                     -0.201680   0.088970  -2.267 0.023402 *  
## AntigenJ08.V1V2.E660.084.AVI.His                                0.481036   0.088970   5.407 6.43e-08 ***
## AntigenJ08.V1V2.E660.2A5.AVI.His                                0.783258   0.088970   8.804  < 2e-16 ***
## AntigenJ08.V1V2.mac239.AVI.His                                  0.922578   0.088970  10.370  < 2e-16 ***
## AntigenSIV.1A11.gp140                                           0.840071   0.088970   9.442  < 2e-16 ***
## AntigenSIV.E543.gp140                                           0.546913   0.088970   6.147 7.91e-10 ***
## AntigenSIVcpz.EK505.gp120                                      -1.268900   0.088970 -14.262  < 2e-16 ***
## AntigenSIVmac239.gp120                                          1.159048   0.088970  13.027  < 2e-16 ***
## AntigenSIVmac239.gp130                                          0.497558   0.088970   5.592 2.24e-08 ***
## AntigenSIVmac239.gp140                                          1.802549   0.088970  20.260  < 2e-16 ***
## AntigenSIVmac239.gp140.AVI.His                                  0.047490   0.088970   0.534 0.593497    
## AntigenSIVmac251.BK.PR55                                        0.178174   0.088970   2.003 0.045219 *  
## AntigenSIVsm.E660.2A5                                           0.591644   0.088970   6.650 2.95e-11 ***
## AntigenSIVsm.E660.84                                            0.037630   0.088970   0.423 0.672332    
## AntigenSIVsmH4.p55.Gag                                          1.060096   0.088970  11.915  < 2e-16 ***
## AntigenV1a                                                      1.436527   0.088970  16.146  < 2e-16 ***
## DetectionReagentaRhIgG.PE.low                                  -1.607196   0.088970 -18.064  < 2e-16 ***
## DetectionReagentC1q                                            -0.967905   0.088970 -10.879  < 2e-16 ***
## DetectionReagentMBL                                            -0.944692   0.088970 -10.618  < 2e-16 ***
## DetectionReagentR2A.2                                          -4.476034   0.088970 -50.310  < 2e-16 ***
## DetectionReagentR2A.3                                           1.021273   0.088970  11.479  < 2e-16 ***
## DetectionReagentR2A.4.high                                      2.165867   0.089068  24.317  < 2e-16 ***
## DetectionReagentR2A.4.low                                      -1.344543   0.088970 -15.112  < 2e-16 ***
## DetectionReagentR3A.1                                           2.764509   0.088970  31.072  < 2e-16 ***
## DetectionReagentR3A.3                                           0.617849   0.088970   6.944 3.82e-12 ***
## AntigenC1.TR:DetectionReagentaRhIgG.PE.low                     -0.083283   0.125822  -0.662 0.508031    
## AntigenG119:DetectionReagentaRhIgG.PE.low                       0.446053   0.125822   3.545 0.000393 ***
## AntigenG145.146:DetectionReagentaRhIgG.PE.low                   0.918066   0.125822   7.297 2.97e-13 ***
## AntigenG49:DetectionReagentaRhIgG.PE.low                        0.304829   0.125822   2.423 0.015408 *  
## AntigenG73:DetectionReagentaRhIgG.PE.low                        0.177073   0.125822   1.407 0.159334    
## AntigenJ08.V1V2.E660.084.AVI.His:DetectionReagentaRhIgG.PE.low  0.027710   0.125822   0.220 0.825694    
## AntigenJ08.V1V2.E660.2A5.AVI.His:DetectionReagentaRhIgG.PE.low  0.089572   0.125822   0.712 0.476531    
## AntigenJ08.V1V2.mac239.AVI.His:DetectionReagentaRhIgG.PE.low    0.162296   0.125822   1.290 0.197095    
## AntigenSIV.1A11.gp140:DetectionReagentaRhIgG.PE.low             0.800763   0.125822   6.364 1.97e-10 ***
## AntigenSIV.E543.gp140:DetectionReagentaRhIgG.PE.low             0.507275   0.125822   4.032 5.54e-05 ***
## AntigenSIVcpz.EK505.gp120:DetectionReagentaRhIgG.PE.low         0.891515   0.125822   7.085 1.39e-12 ***
## AntigenSIVmac239.gp120:DetectionReagentaRhIgG.PE.low            0.198615   0.125822   1.579 0.114447    
## AntigenSIVmac239.gp130:DetectionReagentaRhIgG.PE.low            0.651381   0.125822   5.177 2.26e-07 ***
## AntigenSIVmac239.gp140:DetectionReagentaRhIgG.PE.low            0.585159   0.125822   4.651 3.31e-06 ***
## AntigenSIVmac239.gp140.AVI.His:DetectionReagentaRhIgG.PE.low   -0.046615   0.125822  -0.370 0.711025    
## AntigenSIVmac251.BK.PR55:DetectionReagentaRhIgG.PE.low          0.131000   0.125822   1.041 0.297809    
## AntigenSIVsm.E660.2A5:DetectionReagentaRhIgG.PE.low             0.383031   0.125822   3.044 0.002333 ** 
## AntigenSIVsm.E660.84:DetectionReagentaRhIgG.PE.low              0.458373   0.125822   3.643 0.000270 ***
## AntigenSIVsmH4.p55.Gag:DetectionReagentaRhIgG.PE.low            0.076887   0.125822   0.611 0.541150    
## AntigenV1a:DetectionReagentaRhIgG.PE.low                       -0.290837   0.125822  -2.311 0.020808 *  
## AntigenC1.TR:DetectionReagentC1q                                0.456942   0.125822   3.632 0.000282 ***
## AntigenG119:DetectionReagentC1q                                 0.166004   0.125822   1.319 0.187055    
## AntigenG145.146:DetectionReagentC1q                             0.610301   0.125822   4.850 1.23e-06 ***
## AntigenG49:DetectionReagentC1q                                  0.263364   0.125822   2.093 0.036339 *  
## AntigenG73:DetectionReagentC1q                                  0.596865   0.125822   4.744 2.10e-06 ***
## AntigenJ08.V1V2.E660.084.AVI.His:DetectionReagentC1q           -0.735143   0.125822  -5.843 5.15e-09 ***
## AntigenJ08.V1V2.E660.2A5.AVI.His:DetectionReagentC1q           -1.143785   0.125822  -9.090  < 2e-16 ***
## AntigenJ08.V1V2.mac239.AVI.His:DetectionReagentC1q             -1.506705   0.125822 -11.975  < 2e-16 ***
## AntigenSIV.1A11.gp140:DetectionReagentC1q                      -1.772148   0.125822 -14.085  < 2e-16 ***
## AntigenSIV.E543.gp140:DetectionReagentC1q                      -1.017743   0.125822  -8.089 6.09e-16 ***
## AntigenSIVcpz.EK505.gp120:DetectionReagentC1q                         NA         NA      NA       NA    
## AntigenSIVmac239.gp120:DetectionReagentC1q                     -2.146463   0.125822 -17.059  < 2e-16 ***
## AntigenSIVmac239.gp130:DetectionReagentC1q                     -1.631910   0.125822 -12.970  < 2e-16 ***
## AntigenSIVmac239.gp140:DetectionReagentC1q                     -1.167896   0.125822  -9.282  < 2e-16 ***
## AntigenSIVmac239.gp140.AVI.His:DetectionReagentC1q             -0.188054   0.125822  -1.495 0.135021    
## AntigenSIVmac251.BK.PR55:DetectionReagentC1q                   -1.194844   0.125822  -9.496  < 2e-16 ***
## AntigenSIVsm.E660.2A5:DetectionReagentC1q                      -1.873102   0.125822 -14.887  < 2e-16 ***
## AntigenSIVsm.E660.84:DetectionReagentC1q                       -1.496255   0.125822 -11.892  < 2e-16 ***
## AntigenSIVsmH4.p55.Gag:DetectionReagentC1q                      0.486966   0.125822   3.870 0.000109 ***
## AntigenV1a:DetectionReagentC1q                                 -0.946348   0.125822  -7.521 5.46e-14 ***
## AntigenC1.TR:DetectionReagentMBL                                0.860036   0.125822   6.835 8.22e-12 ***
## AntigenG119:DetectionReagentMBL                                 0.154093   0.125822   1.225 0.220695    
## AntigenG145.146:DetectionReagentMBL                             0.099251   0.125822   0.789 0.430222    
## AntigenG49:DetectionReagentMBL                                 -0.314767   0.125822  -2.502 0.012362 *  
## AntigenG73:DetectionReagentMBL                                 -0.445632   0.125822  -3.542 0.000398 ***
## AntigenJ08.V1V2.E660.084.AVI.His:DetectionReagentMBL           -2.028237   0.125822 -16.120  < 2e-16 ***
## AntigenJ08.V1V2.E660.2A5.AVI.His:DetectionReagentMBL           -2.289136   0.125822 -18.193  < 2e-16 ***
## AntigenJ08.V1V2.mac239.AVI.His:DetectionReagentMBL             -2.584915   0.125822 -20.544  < 2e-16 ***
## AntigenSIV.1A11.gp140:DetectionReagentMBL                             NA         NA      NA       NA    
## AntigenSIV.E543.gp140:DetectionReagentMBL                             NA         NA      NA       NA    
## AntigenSIVcpz.EK505.gp120:DetectionReagentMBL                         NA         NA      NA       NA    
## AntigenSIVmac239.gp120:DetectionReagentMBL                            NA         NA      NA       NA    
## AntigenSIVmac239.gp130:DetectionReagentMBL                            NA         NA      NA       NA    
## AntigenSIVmac239.gp140:DetectionReagentMBL                            NA         NA      NA       NA    
## AntigenSIVmac239.gp140.AVI.His:DetectionReagentMBL             -0.238321   0.125822  -1.894 0.058213 .  
## AntigenSIVmac251.BK.PR55:DetectionReagentMBL                          NA         NA      NA       NA    
## AntigenSIVsm.E660.2A5:DetectionReagentMBL                             NA         NA      NA       NA    
## AntigenSIVsm.E660.84:DetectionReagentMBL                              NA         NA      NA       NA    
## AntigenSIVsmH4.p55.Gag:DetectionReagentMBL                            NA         NA      NA       NA    
## AntigenV1a:DetectionReagentMBL                                 -2.143908   0.125822 -17.039  < 2e-16 ***
## AntigenC1.TR:DetectionReagentR2A.2                                    NA         NA      NA       NA    
## AntigenG119:DetectionReagentR2A.2                               3.813666   0.125822  30.310  < 2e-16 ***
## AntigenG145.146:DetectionReagentR2A.2                           4.047058   0.154100  26.262  < 2e-16 ***
## AntigenG49:DetectionReagentR2A.2                                      NA         NA      NA       NA    
## AntigenG73:DetectionReagentR2A.2                                3.379341   0.125822  26.858  < 2e-16 ***
## AntigenJ08.V1V2.E660.084.AVI.His:DetectionReagentR2A.2                NA         NA      NA       NA    
## AntigenJ08.V1V2.E660.2A5.AVI.His:DetectionReagentR2A.2          1.356078   0.125822  10.778  < 2e-16 ***
## AntigenJ08.V1V2.mac239.AVI.His:DetectionReagentR2A.2            1.041460   0.125822   8.277  < 2e-16 ***
## AntigenSIV.1A11.gp140:DetectionReagentR2A.2                           NA         NA      NA       NA    
## AntigenSIV.E543.gp140:DetectionReagentR2A.2                           NA         NA      NA       NA    
## AntigenSIVcpz.EK505.gp120:DetectionReagentR2A.2                       NA         NA      NA       NA    
## AntigenSIVmac239.gp120:DetectionReagentR2A.2                          NA         NA      NA       NA    
## AntigenSIVmac239.gp130:DetectionReagentR2A.2                          NA         NA      NA       NA    
## AntigenSIVmac239.gp140:DetectionReagentR2A.2                          NA         NA      NA       NA    
## AntigenSIVmac239.gp140.AVI.His:DetectionReagentR2A.2                  NA         NA      NA       NA    
## AntigenSIVmac251.BK.PR55:DetectionReagentR2A.2                        NA         NA      NA       NA    
## AntigenSIVsm.E660.2A5:DetectionReagentR2A.2                           NA         NA      NA       NA    
## AntigenSIVsm.E660.84:DetectionReagentR2A.2                            NA         NA      NA       NA    
## AntigenSIVsmH4.p55.Gag:DetectionReagentR2A.2                          NA         NA      NA       NA    
## AntigenV1a:DetectionReagentR2A.2                                      NA         NA      NA       NA    
## AntigenC1.TR:DetectionReagentR2A.3                              0.649610   0.125822   5.163 2.44e-07 ***
## AntigenG119:DetectionReagentR2A.3                              -0.311143   0.125822  -2.473 0.013405 *  
## AntigenG145.146:DetectionReagentR2A.3                           0.075827   0.125822   0.603 0.546744    
## AntigenG49:DetectionReagentR2A.3                                0.083107   0.125822   0.661 0.508930    
## AntigenG73:DetectionReagentR2A.3                               -0.376912   0.125822  -2.996 0.002740 ** 
## AntigenJ08.V1V2.E660.084.AVI.His:DetectionReagentR2A.3         -1.044096   0.125822  -8.298  < 2e-16 ***
## AntigenJ08.V1V2.E660.2A5.AVI.His:DetectionReagentR2A.3         -1.335725   0.125822 -10.616  < 2e-16 ***
## AntigenJ08.V1V2.mac239.AVI.His:DetectionReagentR2A.3           -1.444237   0.125822 -11.478  < 2e-16 ***
## AntigenSIV.1A11.gp140:DetectionReagentR2A.3                    -0.851760   0.125822  -6.770 1.30e-11 ***
## AntigenSIV.E543.gp140:DetectionReagentR2A.3                    -1.404896   0.125822 -11.166  < 2e-16 ***
## AntigenSIVcpz.EK505.gp120:DetectionReagentR2A.3                -1.750521   0.125822 -13.913  < 2e-16 ***
## AntigenSIVmac239.gp120:DetectionReagentR2A.3                   -2.028412   0.125822 -16.121  < 2e-16 ***
## AntigenSIVmac239.gp130:DetectionReagentR2A.3                   -1.591500   0.125822 -12.649  < 2e-16 ***
## AntigenSIVmac239.gp140:DetectionReagentR2A.3                   -1.247193   0.125822  -9.912  < 2e-16 ***
## AntigenSIVmac239.gp140.AVI.His:DetectionReagentR2A.3           -0.059574   0.125822  -0.473 0.635874    
## AntigenSIVmac251.BK.PR55:DetectionReagentR2A.3                 -3.082091   0.125822 -24.496  < 2e-16 ***
## AntigenSIVsm.E660.2A5:DetectionReagentR2A.3                    -1.872615   0.125822 -14.883  < 2e-16 ***
## AntigenSIVsm.E660.84:DetectionReagentR2A.3                     -1.915880   0.125822 -15.227  < 2e-16 ***
## AntigenSIVsmH4.p55.Gag:DetectionReagentR2A.3                   -1.801965   0.125822 -14.321  < 2e-16 ***
## AntigenV1a:DetectionReagentR2A.3                               -0.985793   0.125822  -7.835 4.74e-15 ***
## AntigenC1.TR:DetectionReagentR2A.4.high                         0.554678   0.125962   4.404 1.07e-05 ***
## AntigenG119:DetectionReagentR2A.4.high                         -0.588959   0.125962  -4.676 2.93e-06 ***
## AntigenG145.146:DetectionReagentR2A.4.high                      0.198208   0.125962   1.574 0.115592    
## AntigenG49:DetectionReagentR2A.4.high                           0.620850   0.125962   4.929 8.28e-07 ***
## AntigenG73:DetectionReagentR2A.4.high                          -0.459769   0.125962  -3.650 0.000262 ***
## AntigenJ08.V1V2.E660.084.AVI.His:DetectionReagentR2A.4.high    -2.130868   0.125962 -16.917  < 2e-16 ***
## AntigenJ08.V1V2.E660.2A5.AVI.His:DetectionReagentR2A.4.high    -2.723410   0.125962 -21.621  < 2e-16 ***
## AntigenJ08.V1V2.mac239.AVI.His:DetectionReagentR2A.4.high      -2.901388   0.125962 -23.034  < 2e-16 ***
## AntigenSIV.1A11.gp140:DetectionReagentR2A.4.high               -4.722954   0.125962 -37.495  < 2e-16 ***
## AntigenSIV.E543.gp140:DetectionReagentR2A.4.high               -3.543717   0.125962 -28.133  < 2e-16 ***
## AntigenSIVcpz.EK505.gp120:DetectionReagentR2A.4.high           -2.999588   0.125962 -23.813  < 2e-16 ***
## AntigenSIVmac239.gp120:DetectionReagentR2A.4.high              -4.466397   0.125962 -35.458  < 2e-16 ***
## AntigenSIVmac239.gp130:DetectionReagentR2A.4.high              -4.253429   0.125962 -33.768  < 2e-16 ***
## AntigenSIVmac239.gp140:DetectionReagentR2A.4.high              -2.916528   0.125962 -23.154  < 2e-16 ***
## AntigenSIVmac239.gp140.AVI.His:DetectionReagentR2A.4.high      -0.466241   0.125962  -3.701 0.000214 ***
## AntigenSIVmac251.BK.PR55:DetectionReagentR2A.4.high            -4.881030   0.125962 -38.750  < 2e-16 ***
## AntigenSIVsm.E660.2A5:DetectionReagentR2A.4.high               -4.221636   0.125962 -33.515  < 2e-16 ***
## AntigenSIVsm.E660.84:DetectionReagentR2A.4.high                -3.983736   0.125962 -31.627  < 2e-16 ***
## AntigenSIVsmH4.p55.Gag:DetectionReagentR2A.4.high              -3.783988   0.125962 -30.041  < 2e-16 ***
## AntigenV1a:DetectionReagentR2A.4.high                          -1.554552   0.125962 -12.341  < 2e-16 ***
## AntigenC1.TR:DetectionReagentR2A.4.low                          0.918651   0.125822   7.301 2.87e-13 ***
## AntigenG119:DetectionReagentR2A.4.low                           0.103950   0.125822   0.826 0.408715    
## AntigenG145.146:DetectionReagentR2A.4.low                       0.273019   0.125822   2.170 0.030018 *  
## AntigenG49:DetectionReagentR2A.4.low                            0.005653   0.125822   0.045 0.964164    
## AntigenG73:DetectionReagentR2A.4.low                           -0.357306   0.125822  -2.840 0.004516 ** 
## AntigenJ08.V1V2.E660.084.AVI.His:DetectionReagentR2A.4.low     -1.756083   0.125822 -13.957  < 2e-16 ***
## AntigenJ08.V1V2.E660.2A5.AVI.His:DetectionReagentR2A.4.low     -1.847211   0.125822 -14.681  < 2e-16 ***
## AntigenJ08.V1V2.mac239.AVI.His:DetectionReagentR2A.4.low       -2.164856   0.125822 -17.206  < 2e-16 ***
## AntigenSIV.1A11.gp140:DetectionReagentR2A.4.low                       NA         NA      NA       NA    
## AntigenSIV.E543.gp140:DetectionReagentR2A.4.low                -2.174434   0.125822 -17.282  < 2e-16 ***
## AntigenSIVcpz.EK505.gp120:DetectionReagentR2A.4.low                   NA         NA      NA       NA    
## AntigenSIVmac239.gp120:DetectionReagentR2A.4.low                      NA         NA      NA       NA    
## AntigenSIVmac239.gp130:DetectionReagentR2A.4.low                      NA         NA      NA       NA    
## AntigenSIVmac239.gp140:DetectionReagentR2A.4.low               -2.063512   0.125822 -16.400  < 2e-16 ***
## AntigenSIVmac239.gp140.AVI.His:DetectionReagentR2A.4.low       -0.201974   0.125822  -1.605 0.108445    
## AntigenSIVmac251.BK.PR55:DetectionReagentR2A.4.low                    NA         NA      NA       NA    
## AntigenSIVsm.E660.2A5:DetectionReagentR2A.4.low                       NA         NA      NA       NA    
## AntigenSIVsm.E660.84:DetectionReagentR2A.4.low                        NA         NA      NA       NA    
## AntigenSIVsmH4.p55.Gag:DetectionReagentR2A.4.low               -2.075917   0.125822 -16.499  < 2e-16 ***
## AntigenV1a:DetectionReagentR2A.4.low                           -1.787953   0.125822 -14.210  < 2e-16 ***
## AntigenC1.TR:DetectionReagentR3A.1                                    NA         NA      NA       NA    
## AntigenG119:DetectionReagentR3A.1                               1.418278   0.125822  11.272  < 2e-16 ***
## AntigenG145.146:DetectionReagentR3A.1                           2.180060   0.154100  14.147  < 2e-16 ***
## AntigenG49:DetectionReagentR3A.1                                2.016470   0.125822  16.026  < 2e-16 ***
## AntigenG73:DetectionReagentR3A.1                                      NA         NA      NA       NA    
## AntigenJ08.V1V2.E660.084.AVI.His:DetectionReagentR3A.1          0.145082   0.125822   1.153 0.248885    
## AntigenJ08.V1V2.E660.2A5.AVI.His:DetectionReagentR3A.1         -0.307374   0.125822  -2.443 0.014570 *  
## AntigenJ08.V1V2.mac239.AVI.His:DetectionReagentR3A.1           -0.373632   0.125822  -2.970 0.002983 ** 
## AntigenSIV.1A11.gp140:DetectionReagentR3A.1                    -3.689456   0.125822 -29.323  < 2e-16 ***
## AntigenSIV.E543.gp140:DetectionReagentR3A.1                    -2.983289   0.125822 -23.710  < 2e-16 ***
## AntigenSIVcpz.EK505.gp120:DetectionReagentR3A.1                -3.228131   0.125822 -25.656  < 2e-16 ***
## AntigenSIVmac239.gp120:DetectionReagentR3A.1                   -3.548092   0.125822 -28.199  < 2e-16 ***
## AntigenSIVmac239.gp130:DetectionReagentR3A.1                   -3.132658   0.125822 -24.897  < 2e-16 ***
## AntigenSIVmac239.gp140:DetectionReagentR3A.1                   -2.859214   0.125822 -22.724  < 2e-16 ***
## AntigenSIVmac239.gp140.AVI.His:DetectionReagentR3A.1                  NA         NA      NA       NA    
## AntigenSIVmac251.BK.PR55:DetectionReagentR3A.1                 -4.607355   0.125822 -36.618  < 2e-16 ***
## AntigenSIVsm.E660.2A5:DetectionReagentR3A.1                    -3.513724   0.125822 -27.926  < 2e-16 ***
## AntigenSIVsm.E660.84:DetectionReagentR3A.1                     -3.452020   0.125822 -27.436  < 2e-16 ***
## AntigenSIVsmH4.p55.Gag:DetectionReagentR3A.1                   -2.941961   0.125822 -23.382  < 2e-16 ***
## AntigenV1a:DetectionReagentR3A.1                                      NA         NA      NA       NA    
## AntigenC1.TR:DetectionReagentR3A.3                              0.661502   0.125822   5.257 1.46e-07 ***
## AntigenG119:DetectionReagentR3A.3                              -0.271430   0.125822  -2.157 0.030989 *  
## AntigenG145.146:DetectionReagentR3A.3                                 NA         NA      NA       NA    
## AntigenG49:DetectionReagentR3A.3                               -0.172113   0.125822  -1.368 0.171345    
## AntigenG73:DetectionReagentR3A.3                               -0.568539   0.125822  -4.519 6.23e-06 ***
## AntigenJ08.V1V2.E660.084.AVI.His:DetectionReagentR3A.3         -0.891727   0.125822  -7.087 1.38e-12 ***
## AntigenJ08.V1V2.E660.2A5.AVI.His:DetectionReagentR3A.3         -1.164807   0.125822  -9.258  < 2e-16 ***
## AntigenJ08.V1V2.mac239.AVI.His:DetectionReagentR3A.3           -1.294087   0.125822 -10.285  < 2e-16 ***
## AntigenSIV.1A11.gp140:DetectionReagentR3A.3                    -1.748686   0.125822 -13.898  < 2e-16 ***
## AntigenSIV.E543.gp140:DetectionReagentR3A.3                    -1.152809   0.125822  -9.162  < 2e-16 ***
## AntigenSIVcpz.EK505.gp120:DetectionReagentR3A.3                -1.385233   0.125822 -11.009  < 2e-16 ***
## AntigenSIVmac239.gp120:DetectionReagentR3A.3                   -1.734044   0.125822 -13.782  < 2e-16 ***
## AntigenSIVmac239.gp130:DetectionReagentR3A.3                   -1.280172   0.125822 -10.174  < 2e-16 ***
## AntigenSIVmac239.gp140:DetectionReagentR3A.3                   -1.035319   0.125822  -8.228  < 2e-16 ***
## AntigenSIVmac239.gp140.AVI.His:DetectionReagentR3A.3           -0.086665   0.125822  -0.689 0.490956    
## AntigenSIVmac251.BK.PR55:DetectionReagentR3A.3                 -2.866817   0.125822 -22.785  < 2e-16 ***
## AntigenSIVsm.E660.2A5:DetectionReagentR3A.3                    -1.586729   0.125822 -12.611  < 2e-16 ***
## AntigenSIVsm.E660.84:DetectionReagentR3A.3                     -1.496486   0.125822 -11.894  < 2e-16 ***
## AntigenSIVsmH4.p55.Gag:DetectionReagentR3A.3                   -1.701034   0.125822 -13.519  < 2e-16 ***
## AntigenV1a:DetectionReagentR3A.3                               -0.962584   0.125822  -7.650 2.02e-14 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Residual standard error: 1.641 on 116721 degrees of freedom
## Multiple R-squared:  0.4653,	Adjusted R-squared:  0.4645 
## F-statistic: 580.3 on 175 and 116721 DF,  p-value: < 2.2e-16

additiveFit <- lm(logreadout ~ GroupNm + Antigen + DetectionReagent, data=FcLmData)
multFit <- lm(logreadout ~ GroupNm + Antigen * DetectionReagent, data=FcLmData)
