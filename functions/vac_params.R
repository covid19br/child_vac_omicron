source("functions/epi_params.R")

Yv <- rep(0,length(Y))
Y <- c(Y,Yv,Yv,Yv,Yv,Yv,Yv)
parameters$beta_vA <- 0.7*parameters$beta
parameters$ihr_vA <- 0.1*parameters$ihr
parameters$ihfr_vA <- 0.1*parameters$ihfr
parameters$asymp_vA <- parameters$asymp + (1-parameters$asymp)*0.5

parameters$beta_wA <- 0.35*parameters$beta
parameters$ihr_wA <- 0.05*parameters$ihr
parameters$ihfr_wA <- 0.05*parameters$ihfr
parameters$asymp_wA <- parameters$asymp + (1-parameters$asymp)*0.7

parameters$beta_vP <- 0.7*parameters$beta
parameters$ihr_vP <- 0.1*parameters$ihr
parameters$ihfr_vP <- 0.1*parameters$ihfr
parameters$asymp_vP <- parameters$asymp + (1-parameters$asymp)*0.5

parameters$beta_wP <- 0.35*parameters$beta
parameters$ihr_wP <- 0.05*parameters$ihr
parameters$ihfr_wP <- 0.05*parameters$ihfr
parameters$asymp_wP <- parameters$asymp + (1-parameters$asymp)*0.7

parameters$beta_vC <- 0.7*parameters$beta
parameters$ihr_vC <- 0.1*parameters$ihr
parameters$ihfr_vC <- 0.1*parameters$ihfr
parameters$asymp_vC <- parameters$asymp + (1-parameters$asymp)*0.5

parameters$beta_wC <- 0.35*parameters$beta
parameters$ihr_wC <- 0.05*parameters$ihr
parameters$ihfr_wC <- 0.05*parameters$ihfr
parameters$asymp_wC <- parameters$asymp + (1-parameters$asymp)*0.7
