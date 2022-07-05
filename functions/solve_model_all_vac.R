library(tidyverse)
library(zeallot)

solve_system <- function(parameters,t){
  Y <- parameters$init.condition
  classes <- c("Su","Eu","Iu","Au","Hu","Ru","Du",
               paste0(c("Sv","Ev","Iv","Av","Hv","Rv","Dv"),rep(c("A","P","C"),each =7)),
               paste0(c("Sw","Ew","Iw","Aw","Hw","Rw","Dw"),rep(c("A","P","C"),each =7)))
  names(Y) <- rep(classes,each=parameters$age.bins)
  c(Su,Eu,Iu,Au,Hu,Ru,Du,
    SvA,EvA,IvA,AvA,HvA,RvA,DvA,
    SvP,EvP,IvP,AvP,HvP,RvP,DvP,
    SvC,EvC,IvC,AvC,HvC,RvC,DvC,
    SwA,EwA,IwA,AwA,HwA,RwA,DwA,
    SwP,EwP,IwP,AwP,HwP,RwP,DwP,
    SwC,EwC,IwC,AwC,HwC,RwC,DwC) %<-% split(Y,factor(names(Y),levels= classes))
  
  open.cov <- NULL
  vac_in_time.P <- matrix(0,nrow = max(t)+1,ncol = 2*parameters$age.bins)
  vac_in_time.A <- matrix(0,nrow = max(t)+1,ncol = 2*parameters$age.bins)
  vac_in_time.C <- matrix(0,nrow = max(t)+1,ncol = 2*parameters$age.bins)
  Yn <- Y
  data <- with(c(parameters),{
    data <- matrix(0,length(t),length(Y))
    
    
    ######VACINAÇÂO - WORK IN PROGRESS
    VacvP <- 1-exp(-1/600)
    VacvA <- 1-exp(-1/400)
    VacvC <- 1-exp(-1/350)
    VacwP <- 1-exp(-1/300)
    VacwA <- 1-exp(-1/200)
    VacwC <- 1-exp(-1/100)
    # VacvP <- 0
    # VacvA <- 0
    # VacvC <- 0
    # VacwP <- 0
    # VacwA <- 0
    # VacwC <- 0
    Vacv <- VacvA + VacvC + VacvP
    Vacw <- VacwA + VacwC + VacwP
    ######
    for(i in t){
      #assuming no difference in infectivity in vaccinated individuals, to keep the code cleaner (and less operations)
      S <- Su + SvA + SvP + SvC + SwA + SwP + SwC 
      E <- Eu + EvA + EvP + EvC + EwA + EwP + EwC
      I <- Iu + IvA + IvP + IvC + IwA + IwP + IwC
      A <- Au + AvA + AvP + AvC + AwA + AwP + AwC
      H <- Hu + HvA + HvP + HvC + HwA + HwP + HwC
      R <- Ru + RvA + RvP + RvC + RwA + RwP + RwC 
      D <- Du + DvA + DvP + DvC + DwA + DwP + DwC #not used
      
      N <- S + E + I + A + H + R
      ######################Vacinação###########################################
      if(vax.cov$vax.cov == TRUE){
        for(age in age.bins:2){###### checa toda iteração, mais simples, faz for sobre as idades
          index <- c(SvPindex[age],EvPindex[age],IvPindex[age],AvPindex[age],HvPindex[age],RvPindex[age],DvPindex[age],
                     SwPindex[age],EwPindex[age],IwPindex[age],AwPindex[age],HwPindex[age],RwPindex[age],DwPindex[age],
                     SvAindex[age],EvAindex[age],IvAindex[age],AvAindex[age],HvAindex[age],RvAindex[age],DvAindex[age],
                     SwAindex[age],EwAindex[age],IwAindex[age],AwAindex[age],HwAindex[age],RwAindex[age],DwAindex[age],
                     SvCindex[age],EvCindex[age],IvCindex[age],AvCindex[age],HvCindex[age],RvCindex[age],DvCindex[age],
                     SwCindex[age],EwCindex[age],IwCindex[age],AwCindex[age],HwCindex[age],RwCindex[age],DwCindex[age])
          index.total <- c(index,Sindex[age],Eindex[age],Iindex[age],Aindex[age],Hindex[age],Rindex[age],Dindex[age])
          age.coverage<- sum(Yn[index])/sum(init.condition[index.total])
          if(age.coverage >= vax.cov$cov){
            open.cov <= age - 1
          }
          else{
            open.cov <- age
            break
          }
        }
      }
      ####### taxa de vacinacao primeira dose
      print(Yn[Rindex])
      
      v1<- VAX.DISTR.RATE(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i],Yn[Sindex]+Yn[Rindex],open = open.cov)
      v1.P <- v1*vac.rate.P[i]/(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i]+1)
      v1.A <- v1*vac.rate.A[i]/(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i]+1)
      v1.C <- v1*vac.rate.C[i]/(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i]+1)
      v1.S.P <- v1.P*Yn[Sindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.R.P <- v1.P*Yn[Rindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.S.A <- v1.A*Yn[Sindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.R.A <- v1.A*Yn[Rindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.S.C <- v1.C*Yn[Sindex]/(Yn[Rindex]+Yn[Sindex]+1)
      v1.R.C <- v1.C*Yn[Rindex]/(Yn[Rindex]+Yn[Sindex]+1)
      vac_in_time.P[i,] <- c(v1.S.P,v1.R.P)
      vac_in_time.A[i,] <- c(v1.S.A,v1.R.A)
      vac_in_time.C[i,] <- c(v1.S.C,v1.R.C)
      print(v1.S.P)
      ####### taxa de vacinacao segunda dose para Pfizer
      if(i > vax.window.days.P){
        v2.S.P <- (1-desist)*vac_in_time.P[i-vax.window.days.P,1:age.bins]
        v2.R.P <- (1-desist)*vac_in_time.P[i-vax.window.days.P,(age.bins+1):(2*age.bins)]
        # ###
        # v2.S.hist.P <- (1-desist)*history.P[i,1:age.bins]
        # v2.R.hist.P <- (1-desist)*history.P[i,(age.bins+1):(2*age.bins)]
        # ###
        # v2.S2.P <- ((1-vax1.beta.P*p)**vax.window.days.P)*v2.S.P + ((1-vax1.beta.P*p)**(i-1))*v2.S.hist.P
        # v2.R2.P <- v2.R.P+v2.R.hist.P+
        #   (1-vax1.ihr.P[,2]*vax1.ihfr.P[,2])*(1-(1-vax1.beta.P*p)**vax.window.days.P)*v2.S.P +
        #   (1-vax1.ihr.P[,2]*vax1.ihfr.P[,2])*(1-(1-vax1.beta.P*p)**(i-1))*v2.S.hist.P
        # ###
        v2.S.P <- pmax(pmin(v2.S.P,Yn[SvPindex]),0)
        v2.R.P <- pmax(pmin(v2.R.P,Yn[RvPindex]),0)
      }
      else{###apenas historico pode vacinar antes de "a" dias
        v2.S.P <- (1-desist)*history.P[i,1:age.bins]
        v2.R.P <- (1-desist)*history.P[i,(age.bins+1):(2*age.bins)]
        # v2.S2.P <- pmax(pmin(((1-vax1.beta.P*p)**(i-1))*v2.S.P,Yn[SvPindex]),0)
        # v2.R.P <- pmax(pmin(v2.R.P+(1-vax1.ihr.P[,2]*vax1.ihfr.P[,2])*(1-(1-vax1.beta.P*p)**(i-1))*v2.S.P,Yn[RvPindex]  ),0)
        # v2.S.P <- v2.S2.P
        v2.S.P <- pmax(pmin(v2.S.P,Yn[SvPindex]),0)
        v2.R.P <- pmax(pmin(v2.R.P,Yn[RvPindex]),0)
      }
      if(i > vax.window.days.A){ #vacinacao segunda dose Astrazeneca
        v2.S.A <- (1-desist)*vac_in_time.A[i-vax.window.days.A,1:age.bins]
        v2.R.A <- (1-desist)*vac_in_time.A[i-vax.window.days.A,(age.bins+1):(2*age.bins)]
        ###
        # v2.S.hist.A <- (1-desist)*history.A[i,1:age.bins]
        # v2.R.hist.A <- (1-desist)*history.A[i,(age.bins+1):(2*age.bins)]
        # ###
        # v2.S2.A <- ((1-vax1.beta.A*p)**vax.window.days.A)*v2.S.A + ((1-vax1.beta.A*p)**(i-1))*v2.S.hist.A
        # v2.R2.A <- v2.R.A+v2.R.hist.A+
        #   (1-vax1.ihr.A[,2]*vax1.ihfr.A[,2])*(1-(1-vax1.beta.A*p)**vax.window.days.A)*v2.S.A +
        #   (1-vax1.ihr.A[,2]*vax1.ihfr.A[,2])*(1-(1-vax1.beta.A*p)**(i-1))*v2.S.hist.A
        # ###
        v2.S.A <- pmax(pmin(v2.S.A,Yn[SvAindex]),0)
        v2.R.A <- pmax(pmin(v2.R.A,Yn[RvAindex]),0)
      }
      else{###apenas historico pode vacinar antes de "a" dias
        v2.S.A <- (1-desist)*history.A[i,1:age.bins]
        v2.R.A <- (1-desist)*history.A[i,(age.bins+1):(2*age.bins)]
        # v2.S2.A <- pmax(pmin(((1-vax1.beta.A*p)**(i-1))*v2.S.A,Yn[SvAindex]),0)
        # v2.R.A <- pmax(pmin(v2.R.A+(1-vax1.ihr.A[,2]*vax1.ihfr.A[,2])*(1-(1-vax1.beta.A*p)**(i-1))*v2.S.A,Yn[RvAindex]  ),0)
        # v2.S.A <- v2.S2.A
        v2.S.A <- pmax(pmin(v2.S.A,Yn[SvAindex]),0)
        v2.R.A <- pmax(pmin(v2.R.A,Yn[RvAindex]),0)
      }
      if(i > vax.window.days.C){
        v2.S.C <- (1-desist)*vac_in_time.C[i-vax.window.days.C,1:age.bins]
        v2.R.C <- (1-desist)*vac_in_time.C[i-vax.window.days.C,(age.bins+1):(2*age.bins)]
        ###
        # v2.S.hist.C <- (1-desist)*history.C[i,1:age.bins]
        # v2.R.hist.C <- (1-desist)*history.C[i,(age.bins+1):(2*age.bins)]
        # ###
        # v2.S2.C <- ((1-vax1.beta.C*p)**vax.window.days.C)*v2.S.C + ((1-vax1.beta.C*p)**(i-1))*v2.S.hist.C
        # v2.R2.C <- v2.R.C+v2.R.hist.C+
        #   (1-vax1.ihr.C[,2]*vax1.ihfr.C[,2])*(1-(1-vax1.beta.C*p)**vax.window.days.C)*v2.S.C +
        #   (1-vax1.ihr.C[,2]*vax1.ihfr.C[,2])*(1-(1-vax1.beta.C*p)**(i-1))*v2.S.hist.C
        # ###
        v2.S.C <- pmax(pmin(v2.S.C,Yn[SvCindex]),0)
        v2.R.C <- pmax(pmin(v2.R.C,Yn[RvCindex]),0)
      }
      else{###apenas historico pode vacinar antes de "a" dias
        v2.S.C <- (1-desist)*history.C[i,1:age.bins]
        v2.R.C <- (1-desist)*history.C[i,(age.bins+1):(2*age.bins)]
        # v2.S2.C <- pmax(pmin(((1-vax1.beta.C*p)**(i-1))*v2.S.C,Yn[SvCindex]),0)
        # v2.R.C <- pmax(pmin(v2.R.C+(1-vax1.ihr.C[,2]*vax1.ihfr.C[,2])*(1-(1-vax1.beta.C*p)**(i-1))*v2.S.C,Yn[RvCindex]  ),0)
        # v2.S.C <- v2.S2.C
        v2.S.C <- pmax(pmin(v2.S.C,Yn[SvCindex]),0)
        v2.R.C <- pmax(pmin(v2.R.C,Yn[RvCindex]),0)
      }
      # VacvP.S <- 1-exp(-1/600)
      # VacvA.S <- 1-exp(-1/600)
      # VacvC.S <- 1-exp(-1/600)
      # VacwP.S <- 1-exp(-1/600)
      # VacwA.S <- 1-exp(-1/600)
      # VacwC.S <- 1-exp(-1/600)
      # VacvP.R <- 1-exp(-1/600)
      # VacvA.R <- 1-exp(-1/600)
      # VacvC.R <- 1-exp(-1/600)
      # VacwP.R <- 1-exp(-1/600)
      # VacwA.R <- 1-exp(-1/600)
      # VacwC.R <- 1-exp(-1/600)
      
      VacvP.S <- v1.S.P/(Su+1)
      VacvA.S <- v1.S.A/(Su+1)
      VacvC.S <- v1.S.C/(Su+1)
      VacwP.S <- v2.S.P/(SvP+1)
      VacwA.S <- v2.S.A/(SvA+1)
      VacwC.S <- v2.S.C/(SvC+1)
      
      VacvP.R <- v1.R.P/(Ru+1)
      VacvA.R <- v1.R.A/(Ru+1)
      VacvC.R <- v1.R.C/(Ru+1)
      VacwP.R <- v2.R.P/(RvP+1)
      VacwA.R <- v2.R.A/(RvA+1)
      VacwC.R <- v2.R.C/(RvC+1)
      
      # VacwP.S <- 0
      # VacwA.S <- 0
      # VacwC.S <- 0
      # VacwP.S <- 0
      # VacwA.S <- 0
      # VacwC.S <- 0
    
      Vacv.S <- VacvA.S + VacvC.S + VacvP.S
      Vacw.S <- VacwA.S + VacwC.S + VacwP.S
      
      Vacv.R <- VacvA.R + VacvC.R + VacvP.R
      Vacw.R <- VacwA.R + VacwC.R + VacwP.R
      ##########################################################################
      #######NAO VACINADOS#######
      ###sobreviventes*(1-vacinados)####
      Su2 <- exp(-beta*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N))*(1-Vacv.S)*Su
      ###novas infecções entre não vacinados
      Eu2 <- (1-exp(-beta*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*(1-Vacv.S)*Su + (1-gamma)*Eu
      Iu2 <- gamma*(1-asymp)*(1-ihr)*Eu+ (1-nu)*Iu
      Au2 <- gamma*asymp*(1-ihr)*Eu + (1-nu)*Au
      Hu2 <- gamma*ihr*Eu + (1-nus)*Hu
      Ru2 <- nu*Iu + nu*Au +nus*(1-ihfr)*Hu + (1-Vacv.R)*Ru
      Du2 <- nus*ihfr*Hu + Du
      #######VACINADOS 1ª DOSE######
      #######ASTRAZENECA
      ###sobreviventes que são novos vacinados###
      SvA2 <- exp(-beta*(c%*%(I + omega*E + omega_a*A + omega_s*H)/N))*VacvA.S*Su + 
        ###sobreviventes entre vacinados D1 que não receberam D2###
        exp(-beta_vA*(c%*%(I + omega*E + omega_a*A + omega_s*H)/N))*(1-VacwA.S)*SvA
      
      ###novas infecções entre novos vacinados
      EvA2 <- (1-exp(-beta*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*VacvA.S*Su + 
        ###novas infecções entre vacinados antigos (note vac D2)
        (1-exp(-beta_vA*(c%*%(I + omega*E + omega_a*A + omega_s*H)/N)))*(1-VacwA.S)*SvA + (1-gamma)*EvA
      IvA2 <- gamma*(1-asymp_vA)*(1-ihr_vA)*EvA+ (1-nu)*IvA
      AvA2 <- gamma*asymp_vA*(1-ihr_vA)*EvA + (1-nu)*AvA
      HvA2 <- gamma*ihr_vA*EvA + (1-nus)*HvA
      RvA2 <- nu*IvA + nu*AvA +nus*(1-ihfr_vA)*HvA + (1-VacwA.R)*RvA + VacvA.R*Ru
      DvA2 <- nus*ihfr_vA*HvA + DvA
      #######PFIZER
      ###sobreviventes que são novos vacinados###
      SvP2 <- exp(-beta*(c%*%(I + omega*E + omega_a*A + omega_s*H)/N))*VacvP.S*Su + 
        ###sobreviventes entre vacinados D1 que não receberam D2###
        exp(-beta_vP*(c%*%(I + omega*E + omega_a*A + omega_s*H)/N))*(1-VacwP.S)*SvP
      
      ###novas infecções entre novos vacinados
      EvP2 <- (1-exp(-beta*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*VacvP.S*Su + 
        ###novas infecções entre vacinados antigos (note vac D2)
        (1-exp(-beta_vP*(c%*%(I + omega*E + omega_a*A + omega_s*H)/N)))*(1-VacwP.S)*SvP + (1-gamma)*EvP
      IvP2 <- gamma*(1-asymp_vP)*(1-ihr_vP)*EvP+ (1-nu)*IvP
      AvP2 <- gamma*asymp_vP*(1-ihr_vP)*EvP + (1-nu)*AvP
      HvP2 <- gamma*ihr_vP*EvP + (1-nus)*HvP
      RvP2 <- nu*IvP + nu*AvP +nus*(1-ihfr_vP)*HvP + (1-VacwP.R)*RvP + VacvP.R*Ru
      DvP2 <- nus*ihfr_vP*HvP + DvP
      #######CORONAVAC
      ###sobreviventes que são novos vacinados###
      SvC2 <- exp(-beta*(c%*%(I + omega*E + omega_a*A + omega_s*H)/N))*VacvC.S*Su + 
        ###sobreviventes entre vacinados D1 que não receberam D2###
        exp(-beta_vC*(c%*%(I + omega*E + omega_a*A + omega_s*H)/N))*(1-VacwC.S)*SvC
      
      ###novas infecções entre novos vacinados
      EvC2 <- (1-exp(-beta*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*VacvC.S*Su + 
        ###novas infecções entre vacinados antigos (note vac D2)
        (1-exp(-beta_vC*(c%*%(I + omega*E + omega_a*A + omega_s*H)/N)))*(1-VacwC.S)*SvC + (1-gamma)*EvC
      IvC2 <- gamma*(1-asymp_vC)*(1-ihr_vC)*EvC+ (1-nu)*IvC
      AvC2 <- gamma*asymp_vC*(1-ihr_vC)*EvC + (1-nu)*AvC
      HvC2 <- gamma*ihr_vC*EvC + (1-nus)*HvC
      RvC2 <- nu*IvC + nu*AvC + nus*(1-ihfr_vC)*HvC + (1-VacwC.R)*RvC + VacvC.R*Ru
      DvC2 <- nus*ihfr_vC*HvC + DvC
      ########VACINADOS 2ª DOSE#######
      ########ASTRAZENECA
      ###sobreviventes D1 que são novos vacinados D2###
      SwA2 <- exp(-beta_vA*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N))*VacwA.S*SvA + 
        ###sobreviventes entre vacinados###
        exp(-beta_wA*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N))*SwA
      
      ###novas infecções entre novos vacinados
      EwA2 <- (1-exp(-beta_vA*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*VacwA.S*SvA + 
        ###novas infecções entre vacinados antigos
        (1-exp(-beta_wA*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*SwA + (1-gamma)*EwA
      IwA2<- gamma*(1-asymp_wA)*(1-ihr_wA)*EwA+ (1-nu)*IwA
      AwA2 <- gamma*asymp_wA*(1-ihr_wA)*EwA + (1-nu)*AwA
      HwA2 <- gamma*ihr_wA*EwA + (1-nus)*HwA
      RwA2 <- nu*IwA + nu*AwA + nus*(1-ihfr_wA)*HwA + RwA + VacwA.R*RvA
      DwA2 <- nus*ihfr_wA*HwA + DwA
      ##########PFIZER
      ###sobreviventes D1 que são novos vacinados D2###
      SwP2 <- exp(-beta_vP*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N))*VacwP.S*SvP + 
        ###sobreviventes entre vacinados###
        exp(-beta_wP*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N))*SwP
      
      ###novas infecções entre novos vacinados
      EwP2 <- (1-exp(-beta_vP*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*VacwP.S*SvP + 
        ###novas infecções entre vacinados antigos
        (1-exp(-beta_wP*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*SwP + (1-gamma)*EwP
      IwP2<- gamma*(1-asymp_wP)*(1-ihr_wP)*EwP + (1-nu)*IwP
      AwP2 <- gamma*asymp_wP*(1-ihr_wP)*EwP + (1-nu)*AwP
      HwP2 <- gamma*ihr_wP*EwP + (1-nus)*HwP
      RwP2 <- nu*IwP + nu*AwP + nus*(1-ihfr_wP)*HwP + RwP + VacwP.R*RvP
      DwP2 <- nus*ihfr_wP*HwP + DwP
      
      ########CORONAVAC
      ###sobreviventes D1 que são novos vacinados D2###
      SwC2 <- exp(-beta_vC*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N))*VacwC.S*SvC + 
        ###sobreviventes entre vacinados###
        exp(-beta_wC*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N))*SwC
      
      ###novas infecções entre novos vacinados
      EwC2 <- (1-exp(-beta_vC*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*VacwC.S*SvC + 
        ###novas infecções entre vacinados antigos
        (1-exp(-beta_wC*(c%*%(I+omega*E + omega_a*A + omega_s*H)/N)))*SwC + (1-gamma)*EwC
      IwC2<- gamma*(1-asymp_wC)*(1-ihr_wC)*EwC+ (1-nu)*IwC
      AwC2 <- gamma*asymp_wC*(1-ihr_wC)*EwC + (1-nu)*AwC
      HwC2 <- gamma*ihr_wC*EwC + (1-nus)*HwC
      RwC2 <- nu*IwC + nu*AwC +nus*(1-ihfr_wC)*HwC + RwC + VacwC.R*RvC
      DwC2 <- nus*ihfr_wC*HwC + DwC
      ###################################
      Su <- Su2
      Eu <- Eu2
      Iu <- Iu2
      Au <- Au2
      Hu <- Hu2
      Ru <- Ru2
      Du <- Du2
      ##
      SvA <- SvA2
      EvA <- EvA2
      IvA <- IvA2 
      AvA <- AvA2 
      HvA <- HvA2 
      RvA <- RvA2 
      DvA <- DvA2
      ##
      SwA <- SwA2
      EwA <- EwA2
      IwA <- IwA2 
      AwA <- AwA2 
      HwA <- HwA2 
      RwA <- RwA2 
      DwA <- DwA2 
      ##
      SvP <- SvP2
      EvP <- EvP2
      IvP <- IvP2 
      AvP <- AvP2 
      HvP <- HvP2 
      RvP <- RvP2 
      DvP <- DvP2
      ##
      SwP <- SwP2
      EwP <- EwP2
      IwP <- IwP2 
      AwP <- AwP2 
      HwP <- HwP2 
      RwP <- RwP2 
      DwP <- DwP2 
      ##
      SvC <- SvC2
      EvC <- EvC2
      IvC <- IvC2 
      AvC <- AvC2 
      HvC <- HvC2 
      RvC <- RvC2 
      DvC <- DvC2
      ##
      SwC <- SwC2
      EwC <- EwC2
      IwC <- IwC2 
      AwC <- AwC2 
      HwC <- HwC2 
      RwC <- RwC2 
      DwC <- DwC2 
      data[i,] <- Yn <-   c(Su,Eu,Iu,Au,Hu,Ru,Du,
                      SvA,EvA,IvA,AvA,HvA,RvA,DvA,
                      SvP,EvP,IvP,AvP,HvP,RvP,DvP,
                      SvC,EvC,IvC,AvC,HvC,RvC,DvC,
                      SwA,EwA,IwA,AwA,HwA,RwA,DwA,
                      SwP,EwP,IwP,AwP,HwP,RwP,DwP,
                      SwC,EwC,IwC,AwC,HwC,RwC,DwC)
      # print(sum(data[i,]))
    }
    data
  })
  data <- as.data.frame(data)
  colnames(data) <- paste0(rep(classes,each=parameters$age.bins),rep(1:parameters$age.bins,length(classes)))
  data$t <- t
  data <- pivot_longer(data,-c(t),"classe")
  data$classe <- factor(data$classe, levels = paste0(rep(classes,each=parameters$age.bins),rep(1:parameters$age.bins,length(classes))))
  return(data)
}

