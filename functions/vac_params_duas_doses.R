source("functions/epi_params.R")

Yv <- rep(0,length(Y))
Y <- c(Y,Yv,Yv)
parameters$beta_v <- 0.7*parameters$beta
parameters$ihr_v <- 0.1*parameters$ihr
parameters$ihfr_v <- 0.1*parameters$ihfr
parameters$asymp_v <- parameters$asymp + (1-parameters$asymp)*0.5

parameters$beta_w <- 0.35*parameters$beta
parameters$ihr_w<- 0.05*parameters$ihr
parameters$ihfr_w <- 0.05*parameters$ihfr
parameters$asymp_w <- parameters$asymp + (1-parameters$asymp)*0.7
