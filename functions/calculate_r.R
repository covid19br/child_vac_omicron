library("dplyr")
library("magrittr")
library(rARPACK)
library(Matrix)
library(matlib)
library(tidyverse)
library(zeallot)
library(stringr)
source("functions/calculate_rdata.R")
#r_data <- estimativa_r("São Paulo")
r_data <- estimativa_r_obs()

source("functions/vac_params_omicron.R")
source("functions/calculate_init_conditions.R")

calculate_r <- function(PREV, params){
    Y3 <- params$init.condition
    Y3 <- with(params,{
      Y3[Ru.Dindex] <- PREV*Y3[Suindex]
      Y3[Suindex] <- (1-PREV)*Y3[Suindex]
      Y3[RwP.Dindex] <- PREV*Y3[SwPindex]
      Y3[SwPindex] <- (1-PREV)*Y3[SwPindex]
      Y3[RvP.Dindex] <- PREV*Y3[SvPindex]
      Y3[SvPindex] <- (1-PREV)*Y3[SvPindex]
      Y3[RwA.Dindex] <- PREV*Y3[SwAindex]
      Y3[SwAindex] <- (1-PREV)*Y3[SwAindex]
      Y3[RvA.Dindex] <- PREV*Y3[SvAindex]
      Y3[SvAindex] <- (1-PREV)*Y3[SvAindex]
      Y3[RwC.Dindex] <- PREV*Y3[SwCindex]
      Y3[SwCindex] <- (1-PREV)*Y3[SwCindex]
      Y3[RvC.Dindex] <- PREV*Y3[SvCindex]
      Y3[SvCindex] <- (1-PREV)*Y3[SvCindex]
      Y3
    })
    params$init.condition <- Y3
    result <- calculate_init_condition(params)
    r <- result$r.O
    return(r)
}


bissec <- function(f, a, b, n = 1000, tol = 0.000000001){
  i = 0
  while(i<n){
    if(f(a)==0.0){
      result <- a
      break
    }
    if(f(b)==0.0){
      result <- b
      break
    }
    c <- (a + b)/2.0
    if(f(c)==0|(b-c)<tol){
      result <- c
      break
    }
    i <- i + 1
    if(sign(f(c))==sign(f(a))){return(bissec(f, c, b))}
    else{return(bissec(f, a, c))}
  }
  if (!result){print('Root finding failed')}
  else{return(result)}
}

Prev.estimada <- bissec(function(PREV){calculate_r(PREV,parameters) - r_data}, 2.4, 2.5) 