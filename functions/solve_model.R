library(tidyverse)
library(zeallot)

solve_system <- function(Y,t,parameters){
  classes <- c("Su","Eu","Iu","Au","Hu","Ru","Du","Sv","Ev","Iv","Av","Hv","Rv","Dv")
  names(Y) <- rep(classes,each=16)
  c(Su,Eu,Iu,Au,Hu,Ru,Du,Sv,Ev,Iv,Av,Hv,Rv,Dv) %<-% split(Y,factor(names(Y),levels= classes))
  data <- matrix(0,length(t),length(Y))
  data <- within(c(parameters,data),{
    vac <- 1-exp(-1/600)
    # vac <- 0
    for(i in t){
      N <- Su + Eu + Iu + Au + Hu + Ru + Sv + Ev + Iv + Av + Hv + Rv
      #######NAO VACINADOS#######
            ###sobreviventes*(1-vacinados)####
      Su2 <- exp(-beta*(c%*%((Iu+Iv)+omega*(Eu+Ev) + omega_a*(Au+Av) + omega_s*(Hu+Hv))/N))*(1-vac)*Su
            ###novas infecções entre não vacinados
      Eu2 <- (1-exp(-beta*(c%*%((Iu+Iv)+omega*(Eu+Ev) + omega_a*(Au+Av) + omega_s*(Hu+Hv))/N)))*(1-vac)*Su + (1-gamma)*Eu
      Iu2 <- gamma*(1-asymp)*(1-ihr)*Eu+ (1-nu)*Iu
      Au2 <- gamma*asymp*(1-ihr)*Eu + (1-nu)*Au
      Hu2 <- gamma*ihr*Eu + (1-nus)*Hu
      Ru2 <- nu*Iu + nu*Au +nus*(1-ihfr)*Hu + Ru
      Du2 <- nus*ihfr*Hu + Du
      #######VACINADOS 1ª DOSE######
             ###sobreviventes que são novos vacinados###
      Sv2 <- exp(-beta*(c%*%((Iu+Iv)+omega*(Eu+Ev) + omega_a*(Au+Av) + omega_s*(Hu+Hv))/N))*vac*Su + 
             ###sobreviventes entre vacinados###
             exp(-beta_v*(c%*%((Iu+Iv)+omega*(Eu+Ev) + omega_a*(Au+Av) + omega_s*(Hu+Hv))/N))*Sv
             ###novas infecções entre novos vacinados
      Ev2 <- (1-exp(-beta*(c%*%((Iu+Iv)+omega*(Eu+Ev) + omega_a*(Au+Av) + omega_s*(Hu+Hv))/N)))*vac*Su + 
             ###novas infecções entre vacinados antigos
             (1-exp(-beta_v*(c%*%((Iu+Iv)+omega*(Eu+Ev) + omega_a*(Au+Av) + omega_s*(Hu+Hv))/N)))*Sv + (1-gamma)*Ev
      Iv2 <- gamma*(1-asymp_v)*(1-ihr_v)*Ev+ (1-nu)*Iv
      Av2 <- gamma*asymp_v*(1-ihr_v)*Ev + (1-nu)*Av
      Hv2 <- gamma*ihr_v*Ev + (1-nus)*Hv
      Rv2 <- nu*Iv + nu*Av +nus*(1-ihfr_v)*Hv + Rv
      Dv2 <- nus*ihfr_v*Hv + Dv
      Su <- Su2
      Eu <- Eu2
      Iu <- Iu2
      Au <- Au2
      Hu <- Hu2
      Ru <- Ru2
      Du <- Du2
      Sv <- Sv2
      Ev <- Ev2
      Iv <- Iv2 
      Av <- Av2 
      Hv <- Hv2 
      Rv <- Rv2 
      Dv <- Dv2 
      data[i,] <- c(Su,Eu,Iu,Au,Hu,Ru,Du,Sv,Ev,Iv,Av,Hv,Rv,Dv)
      print(sum(data[i,]))
    }
  })$data
  data <- as.data.frame(data)
  colnames(data) <- paste0(rep(classes,each=16),rep(1:16,length(classes)))
  data$t <- t
  data <- pivot_longer(data,-c(t),"classe")
  data$classe <- factor(data$classe, levels = paste0(rep(classes,each=16),rep(1:16,length(classes))))
  return(data)
}

