library(tidyverse)
library(zeallot)
library(stringr)

solve_system <- function(parameters,t){
  Y <- parameters$init.condition
  ## Definir estados
  vac.states <- c("u","vA", "vP", "vC", "wA", "wP", "wC")
  inf.states <- c("E", "A", "I", "H", "R", "D")
  # inf.states.2 <- c("E", "A", "I", "H", "R")
  var.states <- c(".D",".O")
  classes = c(paste0("S", vac.states),###susceptiveis vacinados ou não
             paste0(inf.states,rep(vac.states,each = length(inf.states)),
                    rep(var.states,each = length(inf.states)*length(vac.states))) #infectados vacinados ou não
  )
  names(Y) <- rep(classes,each=parameters$age.bins)
  c(Su,SvA,SvP,SvC,SwA,SwP,SwC,
    Eu.D,Au.D,Iu.D,Hu.D,Ru.D,Du.D,
    EvA.D,AvA.D,IvA.D,HvA.D,RvA.D,DvA.D,
    EvP.D,AvP.D,IvP.D,HvP.D,RvP.D,DvP.D,
    EvC.D,AvC.D,IvC.D,HvC.D,RvC.D,DvC.D,
    EwA.D,AwA.D,IwA.D,HwA.D,RwA.D,DwA.D,
    EwP.D,AwP.D,IwP.D,HwP.D,RwP.D,DwP.D,
    EwC.D,AwC.D,IwC.D,HwC.D,RwC.D,DwC.D,
    Eu.O,Au.O,Iu.O,Hu.O,Ru.O,Du.O,
    EvA.O,AvA.O,IvA.O,HvA.O,RvA.O,DvA.O,
    EvP.O,AvP.O,IvP.O,HvP.O,RvP.O,DvP.O,
    EvC.O,AvC.O,IvC.O,HvC.O,RvC.O,DvC.O,
    EwA.O,AwA.O,IwA.O,HwA.O,RwA.O,DwA.O,
    EwP.O,AwP.O,IwP.O,HwP.O,RwP.O,DwP.O,
    EwC.O,AwC.O,IwC.O,HwC.O,RwC.O,DwC.O)%<-% split(Y,factor(names(Y),levels= classes))
  
  open.cov <- NULL
  vac_in_time.P <- matrix(0,nrow = max(t)+1,ncol = 3*parameters$age.bins)
  vac_in_time.A <- matrix(0,nrow = max(t)+1,ncol = 3*parameters$age.bins)
  vac_in_time.C <- matrix(0,nrow = max(t)+1,ncol = 3*parameters$age.bins)
  Yn <- Y
  data <- with(c(parameters),{
    data <- matrix(0,length(t),length(Y))
    age.index <- matrix(0,nrow= age.bins, ncol = length(classes))
    age.index <- as.data.frame(age.index)
    colnames(age.index) <- classes
    for(age in age.bins:1){for(state in classes){
      eval(parse(text = paste0("age.index[",age,",","'",state,"'","]<-",state,"index[",age,"]")))
    }}
    age.index.vac <- age.index[,!str_detect(classes,"u")]
  #   
  #   ######VACINAÇÂO - WORK IN PROGRESS
  #   VacvP <- 1-exp(-1/600)
  #   VacvA <- 1-exp(-1/400)
  #   VacvC <- 1-exp(-1/350)
  #   VacwP <- 1-exp(-1/300)
  #   VacwA <- 1-exp(-1/200)
  #   VacwC <- 1-exp(-1/100)
  #   # VacvP <- 0
  #   # VacvA <- 0
  #   # VacvC <- 0
  #   # VacwP <- 0
  #   # VacwA <- 0
  #   # VacwC <- 0
  #   Vacv <- VacvA + VacvC + VacvP
  #   Vacw <- VacwA + VacwC + VacwP
  #   ######
    for(i in t){
      #assuming no difference in infectivity in vaccinated individuals, to keep the code cleaner (and less operations)
      ##### susceptibles #####
      S <- Su + SvA + SvP + SvC + SwA + SwP + SwC
      ######## delta #########
      E.D <- Eu.D + EvA.D + EvP.D + EvC.D + EwA.D + EwP.D + EwC.D 
      I.D <- Iu.D + IvA.D + IvP.D + IvC.D + IwA.D + IwP.D + IwC.D
      A.D <- Au.D + AvA.D + AvP.D + AvC.D + AwA.D + AwP.D + AwC.D
      H.D <- Hu.D + HvA.D + HvP.D + HvC.D + HwA.D + HwP.D + HwC.D
      R.D <- Ru.D + RvA.D + RvP.D + RvC.D + RwA.D + RwP.D + RwC.D
      D.D <- Du.D + DvA.D + DvP.D + DvC.D + DwA.D + DwP.D + DwC.D #not used
      ####### omicron ########
      E.O <- Eu.O + EvA.O + EvP.O + EvC.O + EwA.O + EwP.O + EwC.O
      I.O <- Iu.O + IvA.O + IvP.O + IvC.O + IwA.O + IwP.O + IwC.O
      A.O <- Au.O + AvA.O + AvP.O + AvC.O + AwA.O + AwP.O + AwC.O
      H.O <- Hu.O + HvA.O + HvP.O + HvC.O + HwA.O + HwP.O + HwC.O
      R.O <- Ru.O + RvA.O + RvP.O + RvC.O + RwA.O + RwP.O + RwC.O
      D.O <- Du.O + DvA.O + DvP.O + DvC.O + DwA.O + DwP.O + DwC.O #not used
      # print(sum(E.O + I.O + A.O + H.O + R.O + D.O))
      ###
      N <- S + E.D + I.D + A.D + H.D + R.D + E.O + I.O + A.O + H.O + R.O
      ######################Vacinação###########################################
      if(vax.cov$vax.cov == TRUE){
        for(age in age.bins:2){###### checa toda iteração, mais simples, faz for sobre as idades
          age.coverage<- sum(Yn[as.numeric(age.index.vac[age,])])/sum(init.condition[as.numeric(age.index[age,])])
          if(age.coverage >= vax.cov$cov){
            open.cov <= age - 1
          }
          else{
            open.cov <- age
            break
          }
        }
      }
  #     ####### taxa de vacinacao primeira dose  #     
      v1<- VAX.DISTR.RATE(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i],
                          Yn[Suindex]+Yn[Ru.Dindex]+Yn[Ru.Oindex],
                          open = open.cov)
      
      v1.P <- v1*vac.rate.P[i]/(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i]+1)
      v1.A <- v1*vac.rate.A[i]/(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i]+1)
      v1.C <- v1*vac.rate.C[i]/(vac.rate.P[i]+vac.rate.A[i]+vac.rate.C[i]+1)
      
      v1.S.P <- v1.P*Yn[Suindex]/(Yn[Ru.Dindex]+Yn[Ru.Oindex]+Yn[Suindex]+1)
      v1.R.P.D <- v1.P*Yn[Ru.Dindex]/(Yn[Ru.Dindex]+Yn[Ru.Oindex]+Yn[Suindex]+1)
      v1.R.P.O <- v1.P*Yn[Ru.Oindex]/(Yn[Ru.Dindex]+Yn[Ru.Oindex]+Yn[Suindex]+1)
      
      v1.S.A <- v1.A*Yn[Suindex]/(Yn[Ru.Dindex]+Yn[Ru.Oindex]+Yn[Suindex]+1)
      v1.R.A.D <- v1.A*Yn[Ru.Dindex]/(Yn[Ru.Dindex]+Yn[Ru.Oindex]+Yn[Suindex]+1)
      v1.R.A.O <- v1.A*Yn[Ru.Oindex]/(Yn[Ru.Dindex]+Yn[Ru.Oindex]+Yn[Suindex]+1)
      
      v1.S.C <- v1.C*Yn[Suindex]/(Yn[Ru.Dindex]+Yn[Ru.Oindex]+Yn[Suindex]+1)
      v1.R.C.D <- v1.C*Yn[Ru.Dindex]/(Yn[Ru.Dindex]+Yn[Ru.Oindex]+Yn[Suindex]+1)
      v1.R.C.O <- v1.C*Yn[Ru.Oindex]/(Yn[Ru.Dindex]+Yn[Ru.Oindex]+Yn[Suindex]+1)
      
      vac_in_time.P[i,] <- c(v1.S.P,v1.R.P.D,v1.R.P.O)
      vac_in_time.A[i,] <- c(v1.S.A,v1.R.A.D,v1.R.A.O)
      vac_in_time.C[i,] <- c(v1.S.C,v1.R.C.D,v1.R.C.O)

      ####### taxa de vacinacao segunda dose para Pfizer
      if(i > vax.window.days.P){
        v2.S.P <- (1-desist)*vac_in_time.P[i-vax.window.days.P,1:age.bins]
        v2.R.P.D <- (1-desist)*vac_in_time.P[i-vax.window.days.P,(age.bins+1):(2*age.bins)]
        v2.R.P.O <- (1-desist)*vac_in_time.P[i-vax.window.days.P,(2*age.bins+1):(3*age.bins)]
        ###
        v2.S.hist.P <- (1-desist)*history.P[i,1:age.bins]
        v2.R.hist.P.D <- (1-desist)*history.P[i,(age.bins+1):(2*age.bins)]
        v2.R.hist.P.O <- (1-desist)*history.P[i,(2*age.bins+1):(3*age.bins)]
        # ###
        # v2.S2.P <- ((1-vax1.beta.P*p)**vax.window.days.P)*v2.S.P + ((1-vax1.beta.P*p)**(i-1))*v2.S.hist.P
        # v2.R2.P <- v2.R.P+v2.R.hist.P+
        #   (1-vax1.ihr.P[,2]*vax1.ihfr.P[,2])*(1-(1-vax1.beta.P*p)**vax.window.days.P)*v2.S.P +
        #   (1-vax1.ihr.P[,2]*vax1.ihfr.P[,2])*(1-(1-vax1.beta.P*p)**(i-1))*v2.S.hist.P
        # ###
        v2.S.P <- pmax(pmin(v2.S.P+v2.S.hist.P,Yn[SvPindex]),0)
        v2.R.P.D <- pmax(pmin(v2.R.P.D+v2.R.hist.P.D,Yn[RvP.Dindex]),0)
        v2.R.P.O <- pmax(pmin(v2.R.P.O+v2.R.hist.P.O,Yn[RvP.Oindex]),0)
      }
      else{###apenas historico pode vacinar antes de "a" dias
        v2.S.P <- (1-desist)*history.P[i,1:age.bins]
        v2.R.P.D <- (1-desist)*history.P[i,(age.bins+1):(2*age.bins)]
        v2.R.P.O <- (1-desist)*history.P[i,(2*age.bins+1):(3*age.bins)]
        
        # v2.S2.P <- pmax(pmin(((1-vax1.beta.P*p)**(i-1))*v2.S.P,Yn[SvPindex]),0)
        # v2.R.P <- pmax(pmin(v2.R.P+(1-vax1.ihr.P[,2]*vax1.ihfr.P[,2])*(1-(1-vax1.beta.P*p)**(i-1))*v2.S.P,Yn[RvPindex]  ),0)
        # v2.S.P <- v2.S2.P
        v2.S.P <- pmax(pmin(v2.S.P,Yn[SvPindex]),0)
        v2.R.P.D <- pmax(pmin(v2.R.P.D,Yn[RvP.Dindex]),0)
        v2.R.P.O <- pmax(pmin(v2.R.P.O,Yn[RvP.Oindex]),0)
        
      }
      if(i > vax.window.days.A){ #vacinacao segunda dose Astrazeneca
        v2.S.A <- (1-desist)*vac_in_time.A[i-vax.window.days.A,1:age.bins]
        v2.R.A.D <- (1-desist)*vac_in_time.A[i-vax.window.days.A,(age.bins+1):(2*age.bins)]
        v2.R.A.O <- (1-desist)*vac_in_time.A[i-vax.window.days.A,(2*age.bins+1):(3*age.bins)]
        
        ###
        v2.S.hist.A <- (1-desist)*history.A[i,1:age.bins]
        v2.R.hist.A.D <- (1-desist)*history.A[i,(age.bins+1):(2*age.bins)]
        v2.R.hist.A.O <- (1-desist)*history.A[i,(2*age.bins+1):(3*age.bins)]
        # ###
        # v2.S2.A <- ((1-vax1.beta.A*p)**vax.window.days.A)*v2.S.A + ((1-vax1.beta.A*p)**(i-1))*v2.S.hist.A
        # v2.R2.A <- v2.R.A+v2.R.hist.A+
        #   (1-vax1.ihr.A[,2]*vax1.ihfr.A[,2])*(1-(1-vax1.beta.A*p)**vax.window.days.A)*v2.S.A +
        #   (1-vax1.ihr.A[,2]*vax1.ihfr.A[,2])*(1-(1-vax1.beta.A*p)**(i-1))*v2.S.hist.A
        # ###
        v2.S.A <- pmax(pmin(v2.S.A+v2.S.hist.A,Yn[SvAindex]),0)
        v2.R.A.D <- pmax(pmin(v2.R.A.D+v2.R.hist.A.D,Yn[RvA.Dindex]),0)
        v2.R.A.O <- pmax(pmin(v2.R.A.O+v2.R.hist.A.O,Yn[RvA.Oindex]),0)
      }
      else{###apenas historico pode vacinar antes de "a" dias
        v2.S.A <- (1-desist)*history.A[i,1:age.bins]
        v2.R.A.D <- (1-desist)*history.A[i,(age.bins+1):(2*age.bins)]
        v2.R.A.O <- (1-desist)*history.A[i,(2*age.bins+1):(3*age.bins)]
        
        # v2.S2.A <- pmax(pmin(((1-vax1.beta.A*p)**(i-1))*v2.S.A,Yn[SvAindex]),0)
        # v2.R.A <- pmax(pmin(v2.R.A+(1-vax1.ihr.A[,2]*vax1.ihfr.A[,2])*(1-(1-vax1.beta.A*p)**(i-1))*v2.S.A,Yn[RvAindex]  ),0)
        # v2.S.A <- v2.S2.A
        v2.S.A <- pmax(pmin(v2.S.A,Yn[SvAindex]),0)
        v2.R.A.D <- pmax(pmin(v2.R.A.D,Yn[RvA.Dindex]),0)
        v2.R.A.O <- pmax(pmin(v2.R.A.O,Yn[RvA.Oindex]),0)
      }
      if(i > vax.window.days.C){
        v2.S.C <- (1-desist)*vac_in_time.C[i-vax.window.days.C,1:age.bins]
        v2.R.C.D <- (1-desist)*vac_in_time.C[i-vax.window.days.C,(age.bins+1):(2*age.bins)]
        v2.R.C.O <- (1-desist)*vac_in_time.C[i-vax.window.days.C,(2*age.bins+1):(3*age.bins)]
        ###
        v2.S.hist.C <- (1-desist)*history.C[i,1:age.bins]
        v2.R.hist.C.D <- (1-desist)*history.C[i,(age.bins+1):(2*age.bins)]
        v2.R.hist.C.O <- (1-desist)*history.C[i,(2*age.bins+1):(3*age.bins)]
        # ###
        # v2.S2.C <- ((1-vax1.beta.C*p)**vax.window.days.C)*v2.S.C + ((1-vax1.beta.C*p)**(i-1))*v2.S.hist.C
        # v2.R2.C <- v2.R.C+v2.R.hist.C+
        #   (1-vax1.ihr.C[,2]*vax1.ihfr.C[,2])*(1-(1-vax1.beta.C*p)**vax.window.days.C)*v2.S.C +
        #   (1-vax1.ihr.C[,2]*vax1.ihfr.C[,2])*(1-(1-vax1.beta.C*p)**(i-1))*v2.S.hist.C
        # ###
        v2.S.C <- pmax(pmin(v2.S.C + v2.S.hist.C,Yn[SvCindex]),0)
        v2.R.C.D <- pmax(pmin(v2.R.C.D + v2.R.hist.C.D,Yn[RvC.Dindex]),0)
        v2.R.C.O <- pmax(pmin(v2.R.C.O + v2.R.hist.C.O,Yn[RvC.Oindex]),0)
        
      }
      else{###apenas historico pode vacinar antes de "a" dias
        v2.S.C <- (1-desist)*history.C[i,1:age.bins]
        v2.R.C.D <- (1-desist)*history.C[i,(age.bins+1):(2*age.bins)]
        v2.R.C.O <- (1-desist)*history.C[i,(2*age.bins+1):(3*age.bins)]
        # v2.S2.C <- pmax(pmin(((1-vax1.beta.C*p)**(i-1))*v2.S.C,Yn[SvCindex]),0)
        # v2.R.C <- pmax(pmin(v2.R.C+(1-vax1.ihr.C[,2]*vax1.ihfr.C[,2])*(1-(1-vax1.beta.C*p)**(i-1))*v2.S.C,Yn[RvCindex]  ),0)
        # v2.S.C <- v2.S2.C
        v2.S.C <- pmax(pmin(v2.S.C,Yn[SvCindex]),0)
        v2.R.C.D <- pmax(pmin(v2.R.C.D,Yn[RvC.Dindex]),0)
        v2.R.C.O <- pmax(pmin(v2.R.C.O,Yn[RvC.Oindex]),0)
      }
      VacvP.S <- v1.S.P/(Su+1)
      VacvA.S <- v1.S.A/(Su+1)
      VacvC.S <- v1.S.C/(Su+1)
      VacwP.S <- v2.S.P/(SvP+1)
      VacwA.S <- v2.S.A/(SvA+1)
      VacwC.S <- v2.S.C/(SvC+1)
      # print(sum(VacvP.S))
      VacvP.R.D <- v1.R.P.D/(Ru.D+1)
      VacvA.R.D <- v1.R.A.D/(Ru.D+1)
      VacvC.R.D <- v1.R.C.D/(Ru.D+1)
      VacwP.R.D <- v2.R.P.D/(RvP.D+1)
      VacwA.R.D <- v2.R.A.D/(RvA.D+1)
      VacwC.R.D <- v2.R.C.D/(RvC.D+1)

      VacvP.R.O <- v1.R.P.O/(Ru.O+1)
      VacvA.R.O <- v1.R.A.O/(Ru.O+1)
      VacvC.R.O <- v1.R.C.O/(Ru.O+1)
      VacwP.R.O <- v2.R.P.O/(RvP.O+1)
      VacwA.R.O <- v2.R.A.O/(RvA.O+1)
      VacwC.R.O <- v2.R.C.O/(RvC.O+1)

      Vacv.S <- VacvA.S + VacvC.S + VacvP.S
      Vacw.S <- VacwA.S + VacwC.S + VacwP.S

      Vacv.R.D <- VacvA.R.D + VacvC.R.D + VacvP.R.D
      Vacw.R.D <- VacwA.R.D + VacwC.R.D + VacwP.R.D
      
      Vacv.R.O <- VacvA.R.O + VacvC.R.O + VacvP.R.O
      Vacw.R.O <- VacwA.R.O + VacwC.R.O + VacwP.R.O
      
      ##########################################################################
      ####
      ###potential contacts with infected individuals
      C.Delta <- (c%*%(I.D+omega*E.D + omega_a*A.D + omega_s*H.D)/N)
      C.Omicron <- (c%*%(I.O+omega*E.O + omega_a*A.O + omega_s*H.O)/N)
      #######NAO VACINADOS#######
      ###sobreviventes*(1-vacinados)####
      Su2 <- exp(-beta.D*C.Delta - beta.O*C.Omicron)*(1-Vacv.S)*Su
      #######VACINADOS 1ª DOSE######
      #######ASTRAZENECA
      ###sobreviventes que são novos vacinados###
      SvA2 <- exp(-beta.D*C.Delta - beta.O*C.Omicron)*VacvA.S*Su +
        ###sobreviventes entre vacinados D1 que não receberam D2###
              exp(-beta_vA.D*C.Delta - beta_vA.O*C.Omicron)*(1-VacwA.S)*SvA
      #######PFIZER
      ###sobreviventes que são novos vacinados###
      SvP2 <- exp(-beta.D*C.Delta - beta.O*C.Omicron)*VacvP.S*Su +
      ###sobreviventes entre vacinados D1 que não receberam D2###
              exp(-beta_vP.D*C.Delta - beta_vP.O*C.Omicron)*(1-VacwP.S)*SvP
      #######CORONAVAC
      ###sobreviventes que são novos vacinados###
      SvC2 <- exp(-beta.D*C.Delta - beta.O*C.Omicron)*VacvC.S*Su +
      ###sobreviventes entre vacinados D1 que não receberam D2###
              exp(-beta_vC.D*C.Delta - beta_vC.O*C.Omicron)*(1-VacwC.S)*SvC
      ########VACINADOS 2ª DOSE#######
      ########ASTRAZENECA
      ###sobreviventes D1 que são novos vacinados D2###
      SwA2 <- exp(-beta_vA.D*C.Delta - beta_vA.O*C.Omicron)*VacwA.S*SvA +
      ###sobreviventes entre vacinados###
              exp(-beta_wA.D*C.Delta - beta_wA.O*C.Omicron)*SwA
    
      ##########PFIZER
      ###sobreviventes D1 que são novos vacinados D2###
      SwP2 <- exp(-beta_vP.D*C.Delta - beta_vP.O*C.Omicron)*VacwP.S*SvP +
      ###sobreviventes entre vacinados###
        exp(-beta_wP.D*C.Delta - beta_wP.O*C.Omicron)*SwP

      ########CORONAVAC
      ###sobreviventes D1 que são novos vacinados D2###
      SwC2 <- exp(-beta_vC.D*C.Delta - beta_vC.O*C.Omicron)*VacwC.S*SvC +
      ###sobreviventes entre vacinados###
            exp(-beta_wC.D*C.Delta - beta_wC.O*C.Omicron)*SwC
      ########################  DELTA  #########################################
      ###novas infecções entre não vacinados###
      Eu.D2 <- (beta.D*C.Delta)/(beta.D*C.Delta + beta.O*C.Omicron)*(1-exp(-beta.D*C.Delta -beta.O*C.Omicron))*(1-Vacv.S)*Su + (1-gamma)*Eu.D
      Iu.D2 <- gamma*(1-asymp)*(1-ihr)*Eu.D+ (1-nu)*Iu.D
      Au.D2 <- gamma*asymp*(1-ihr)*Eu.D + (1-nu)*Au.D
      Hu.D2 <- gamma*ihr*Eu.D + (1-nus)*Hu.D
      Ru.D2 <- nu*Iu.D + nu*Au.D +nus*(1-ihfr)*Hu.D + (1-Vacv.R.D)*Ru.D
      Du.D2 <- nus*ihfr*Hu.D + Du.D
      ################FIRST DOSE#################
      #######ASTRAZENECA######
      ###novas infecções entre novos vacinados
      EvA.D2 <- (beta.D*C.Delta)/(beta.D*C.Delta + beta.O*C.Omicron)*(1-exp(-beta.D*C.Delta - beta.O*C.Omicron))*VacvA.S*Su +
        ###novas infecções entre vacinados antigos (note vac D2)
        (beta_vA.D*C.Delta)/(beta_vA.D*C.Delta + beta_vA.O*C.Omicron)*(1-exp(-beta_vA.D*C.Delta - beta_vA.O*C.Omicron))*(1-VacwA.S)*SvA + (1-gamma)*EvA.D
      IvA.D2 <- gamma*(1-asymp_vA.D)*(1-ihr_vA.D)*EvA.D+ (1-nu)*IvA.D
      AvA.D2 <- gamma*asymp_vA.D*(1-ihr_vA.D)*EvA.D + (1-nu)*AvA.D
      HvA.D2 <- gamma*ihr_vA.D*EvA.D + (1-nus)*HvA.D
      RvA.D2 <- nu*IvA.D + nu*AvA.D +nus*(1-ihfr_vA.D)*HvA.D + (1-VacwA.R.D)*RvA.D + VacvA.R.D*Ru.D
      DvA.D2 <- nus*ihfr_vA.D*HvA.D + DvA.D

      ########PFIZER#########
      ###novas infecções entre novos vacinados
      EvP.D2 <-  (beta.D*C.Delta)/(beta.D*C.Delta + beta.O*C.Omicron)*(1-exp(-beta.D*C.Delta - beta.O*C.Omicron))*VacvP.S*Su +
        ###novas infecções entre vacinados antigos (note vac D2)
        (beta_vP.D*C.Delta)/(beta_vP.D*C.Delta + beta_vP.O*C.Omicron)*(1-exp(-beta_vP.D*C.Delta - beta_vP.O*C.Omicron))*(1-VacwP.S)*SvP + (1-gamma)*EvP.D
      IvP.D2 <- gamma*(1-asymp_vP.D)*(1-ihr_vP.D)*EvP.D+ (1-nu)*IvP.D
      AvP.D2 <- gamma*asymp_vP.D*(1-ihr_vP.D)*EvP.D + (1-nu)*AvP.D
      HvP.D2 <- gamma*ihr_vP.D*EvP.D + (1-nus)*HvP.D
      RvP.D2 <- nu*IvP.D + nu*AvP.D +nus*(1-ihfr_vP.D)*HvP.D + (1-VacwP.R.D)*RvP.D + VacvP.R.D*Ru.D
      DvP.D2 <- nus*ihfr_vP.D*HvP.D + DvP.D

      #######CORONAVAC#######
      ###novas infecções entre novos vacinados
      EvC.D2 <- (beta.D*C.Delta)/(beta.D*C.Delta + beta.O*C.Omicron)*(1-exp(-beta.D*C.Delta - beta.O*C.Omicron))*VacvC.S*Su +
        ###novas infecções entre vacinados antigos (note vac D2)
        (beta_vC.D*C.Delta)/(beta_vC.D*C.Delta + beta_vC.O*C.Omicron)*(1-exp(-beta_vC.D*C.Delta - beta_vC.O*C.Omicron))*(1-VacwC.S)*SvC + (1-gamma)*EvC.D
      IvC.D2 <- gamma*(1-asymp_vC.D)*(1-ihr_vC.D)*EvC.D+ (1-nu)*IvC.D
      AvC.D2 <- gamma*asymp_vC.D*(1-ihr_vC.D)*EvC.D + (1-nu)*AvC.D
      HvC.D2 <- gamma*ihr_vC.D*EvC.D + (1-nus)*HvC.D
      RvC.D2 <- nu*IvC.D + nu*AvC.D + nus*(1-ihfr_vC.D)*HvC.D + (1-VacwC.R.D)*RvC.D + VacvC.R.D*Ru.D
      DvC.D2 <- nus*ihfr_vC.D*HvC.D + DvC.D
      
      ############SECOND DOSE################
      #########ASTRAZENECA
      ###novas infecções entre novos vacinados
      EwA.D2 <- (beta_vA.D*C.Delta)/(beta_vA.D*C.Delta + beta_vA.O*C.Omicron)*(1-exp(-beta_vA.D*C.Delta - beta_vA.O*C.Omicron))*VacwA.S*SvA +
        ###novas infecções entre vacinados antigos
        (beta_wA.D*C.Delta)/(beta_wA.D*C.Delta + beta_wA.O*C.Omicron)*(1-exp(-beta_wA.D*C.Delta - beta_wA.O*C.Omicron))*SwA + (1-gamma)*EwA.D
      IwA.D2<- gamma*(1-asymp_wA.D)*(1-ihr_wA.D)*EwA.D+ (1-nu)*IwA.D
      AwA.D2 <- gamma*asymp_wA.D*(1-ihr_wA.D)*EwA.D + (1-nu)*AwA.D
      HwA.D2 <- gamma*ihr_wA.D*EwA.D + (1-nus)*HwA.D
      RwA.D2 <- nu*IwA.D + nu*AwA.D + nus*(1-ihfr_wA.D)*HwA.D + RwA.D + VacwA.R.D*RvA.D
      DwA.D2 <- nus*ihfr_wA.D*HwA.D + DwA.D

      #####PFIZER     
      ###novas infecções entre novos vacinados
      EwP.D2 <- (beta_vP.D*C.Delta)/(beta_vP.D*C.Delta + beta_vP.O*C.Omicron)*(1-exp(-beta_vP.D*C.Delta - beta_vP.O*C.Omicron))*VacwP.S*SvP +
        ###novas infecções entre vacinados antigos
        (beta_wP.D*C.Delta)/(beta_wP.D*C.Delta + beta_wP.O*C.Omicron)*(1-exp(-beta_wP.D*C.Delta - beta_wP.O*C.Omicron))*SwP + (1-gamma)*EwP.D
      IwP.D2<- gamma*(1-asymp_wP.D)*(1-ihr_wP.D)*EwP.D + (1-nu)*IwP.D
      AwP.D2 <- gamma*asymp_wP.D*(1-ihr_wP.D)*EwP.D + (1-nu)*AwP.D
      HwP.D2 <- gamma*ihr_wP.D*EwP.D + (1-nus)*HwP.D
      RwP.D2 <- nu*IwP.D + nu*AwP.D + nus*(1-ihfr_wP.D)*HwP.D + RwP.D + VacwP.R.D*RvP.D
      DwP.D2 <- nus*ihfr_wP.D*HwP.D + DwP.D

      #######CORONAVAC
      ###novas infecções entre novos vacinados
      EwC.D2 <- (beta_vC.D*C.Delta)/(beta_vC.D*C.Delta + beta_vC.O*C.Omicron)*(1-exp(-beta_vC.D*C.Delta - beta_vC.O*C.Omicron))*VacwC.S*SvC +
        ###novas infecções entre vacinados antigos
        (beta_wC.D*C.Delta)/(beta_wC.D*C.Delta + beta_wC.O*C.Omicron)*(1-exp(-beta_wC.D*C.Delta - beta_wC.O*C.Omicron))*SwC + (1-gamma)*EwC.D
      IwC.D2<- gamma*(1-asymp_wC.D)*(1-ihr_wC.D)*EwC.D+ (1-nu)*IwC.D
      AwC.D2 <- gamma*asymp_wC.D*(1-ihr_wC.D)*EwC.D + (1-nu)*AwC.D
      HwC.D2 <- gamma*ihr_wC.D*EwC.D + (1-nus)*HwC.D
      RwC.D2 <- nu*IwC.D + nu*AwC.D +nus*(1-ihfr_wC.D)*HwC.D + RwC.D + VacwC.R.D*RvC.D
      DwC.D2 <- nus*ihfr_wC.D*HwC.D + DwC.D
      
      ########################  OMICRON  #######################################
      ###novas infecções entre não vacinados###
      Eu.O2 <- (beta.O*C.Omicron)/(beta.D*C.Delta + beta.O*C.Omicron)*(1-exp(-beta.D*C.Delta -beta.O*C.Omicron))*(1-Vacv.S)*Su + (1-gamma)*Eu.O
      Iu.O2 <- gamma*(1-asymp)*(1-ihr)*Eu.O+ (1-nu)*Iu.O
      Au.O2 <- gamma*asymp*(1-ihr)*Eu.O + (1-nu)*Au.O
      Hu.O2 <- gamma*ihr*Eu.O + (1-nus)*Hu.O
      Ru.O2 <- nu*Iu.O + nu*Au.O +nus*(1-ihfr)*Hu.O + (1-Vacv.R.O)*Ru.O
      Du.O2 <- nus*ihfr*Hu.O + Du.O
      ################FIRST DOSE#################
      #######ASTRAZENECA######
      ###novas infecções entre novos vacinados
      EvA.O2 <- (beta.O*C.Omicron)/(beta.D*C.Delta + beta.O*C.Omicron)*(1-exp(-beta.D*C.Delta - beta.O*C.Omicron))*VacvA.S*Su +
        ###novas infecções entre vacinados antigos (note vac D2)
        (beta_vA.O*C.Omicron)/(beta_vA.D*C.Delta + beta_vA.O*C.Omicron)*(1-exp(-beta_vA.D*C.Delta - beta_vA.O*C.Omicron))*(1-VacwA.S)*SvA + (1-gamma)*EvA.O
      IvA.O2 <- gamma*(1-asymp_vA.O)*(1-ihr_vA.O)*EvA.O+ (1-nu)*IvA.O
      AvA.O2 <- gamma*asymp_vA.O*(1-ihr_vA.O)*EvA.O + (1-nu)*AvA.O
      HvA.O2 <- gamma*ihr_vA.O*EvA.O + (1-nus)*HvA.O
      RvA.O2 <- nu*IvA.O + nu*AvA.O +nus*(1-ihfr_vA.O)*HvA.O + (1-VacwA.R.O)*RvA.O + VacvA.R.O*Ru.O
      DvA.O2 <- nus*ihfr_vA.O*HvA.O + DvA.O
      
      ########PFIZER#########
      ###novas infecções entre novos vacinados
      EvP.O2 <-  (beta.O*C.Omicron)/(beta.D*C.Delta + beta.O*C.Omicron)*(1-exp(-beta.D*C.Delta - beta.O*C.Omicron))*VacvP.S*Su +
        ###novas infecções entre vacinados antigos (note vac D2)
        (beta_vP.O*C.Omicron)/(beta_vP.D*C.Delta + beta_vP.O*C.Omicron)*(1-exp(-beta_vP.D*C.Delta - beta_vP.O*C.Omicron))*(1-VacwP.S)*SvP + (1-gamma)*EvP.O
      IvP.O2 <- gamma*(1-asymp_vP.O)*(1-ihr_vP.O)*EvP.O+ (1-nu)*IvP.O
      AvP.O2 <- gamma*asymp_vP.O*(1-ihr_vP.O)*EvP.O + (1-nu)*AvP.O
      HvP.O2 <- gamma*ihr_vP.O*EvP.O + (1-nus)*HvP.O
      RvP.O2 <- nu*IvP.O + nu*AvP.O +nus*(1-ihfr_vP.O)*HvP.O + (1-VacwP.R.O)*RvP.O + VacvP.R.O*Ru.O
      DvP.O2 <- nus*ihfr_vP.O*HvP.O + DvP.O
      
      #######CORONAVAC#######
      ###novas infecções entre novos vacinados
      EvC.O2 <- (beta.O*C.Omicron)/(beta.D*C.Delta + beta.O*C.Omicron)*(1-exp(-beta.D*C.Delta - beta.O*C.Omicron))*VacvC.S*Su +
        ###novas infecções entre vacinados antigos (note vac D2)
        (beta_vC.O*C.Omicron)/(beta_vC.D*C.Delta + beta_vC.O*C.Omicron)*(1-exp(-beta_vC.D*C.Delta - beta_vC.O*C.Omicron))*(1-VacwC.S)*SvC + (1-gamma)*EvC.O
      IvC.O2 <- gamma*(1-asymp_vC.O)*(1-ihr_vC.O)*EvC.O + (1-nu)*IvC.O
      AvC.O2 <- gamma*asymp_vC.O*(1-ihr_vC.O)*EvC.O + (1-nu)*AvC.O
      HvC.O2 <- gamma*ihr_vC.O*EvC.O + (1-nus)*HvC.O
      RvC.O2 <- nu*IvC.O + nu*AvC.O + nus*(1-ihfr_vC.O)*HvC.O + (1-VacwC.R.O)*RvC.O + VacvC.R.O*Ru.O
      DvC.O2 <- nus*ihfr_vC.O*HvC.O + DvC.O
      
      ############SECOND DOSE################
      #########ASTRAZENECA
      ###novas infecções entre novos vacinados
      EwA.O2 <- (beta_vA.O*C.Omicron)/(beta_vA.D*C.Delta + beta_vA.O*C.Omicron)*(1-exp(-beta_vA.D*C.Delta - beta_vA.O*C.Omicron))*VacwA.S*SvA +
        ###novas infecções entre vacinados antigos
        (beta_wA.O*C.Omicron)/(beta_wA.D*C.Delta + beta_wA.O*C.Omicron)*(1-exp(-beta_wA.D*C.Delta - beta_wA.O*C.Omicron))*SwA + (1-gamma)*EwA.O
      IwA.O2<- gamma*(1-asymp_wA.O)*(1-ihr_wA.O)*EwA.O+ (1-nu)*IwA.O
      AwA.O2 <- gamma*asymp_wA.O*(1-ihr_wA.O)*EwA.O + (1-nu)*AwA.O
      HwA.O2 <- gamma*ihr_wA.O*EwA.O + (1-nus)*HwA.O
      RwA.O2 <- nu*IwA.O + nu*AwA.O + nus*(1-ihfr_wA.O)*HwA.O + RwA.O + VacwA.R.O*RvA.O
      DwA.O2 <- nus*ihfr_wA.O*HwA.O + DwA.O
      
      #####PFIZER     
      ###novas infecções entre novos vacinados
      EwP.O2 <- (beta_vP.O*C.Omicron)/(beta_vP.D*C.Delta + beta_vP.O*C.Omicron)*(1-exp(-beta_vP.D*C.Delta - beta_vP.O*C.Omicron))*VacwP.S*SvP +
        ###novas infecções entre vacinados antigos
        (beta_wP.O*C.Omicron)/(beta_wP.D*C.Delta + beta_wP.O*C.Omicron)*(1-exp(-beta_wP.D*C.Delta - beta_wP.O*C.Omicron))*SwP + (1-gamma)*EwP.O
      IwP.O2<- gamma*(1-asymp_wP.O)*(1-ihr_wP.O)*EwP.O + (1-nu)*IwP.O
      AwP.O2 <- gamma*asymp_wP.O*(1-ihr_wP.O)*EwP.O + (1-nu)*AwP.O
      HwP.O2 <- gamma*ihr_wP.O*EwP.O + (1-nus)*HwP.O
      RwP.O2 <- nu*IwP.O + nu*AwP.O + nus*(1-ihfr_wP.O)*HwP.O + RwP.O + VacwP.R.O*RvP.O
      DwP.O2 <- nus*ihfr_wP.O*HwP.O + DwP.O
      
      #######CORONAVAC
      ###novas infecções entre novos vacinados
      EwC.O2 <- (beta_vC.O*C.Omicron)/(beta_vC.D*C.Delta + beta_vC.O*C.Omicron)*(1-exp(-beta_vC.D*C.Delta - beta_vC.O*C.Omicron))*VacwC.S*SvC +
        ###novas infecções entre vacinados antigos
        (beta_wC.O*C.Omicron)/(beta_wC.D*C.Delta + beta_wC.O*C.Omicron)*(1-exp(-beta_wC.D*C.Delta - beta_wC.O*C.Omicron))*SwC + (1-gamma)*EwC.O
      IwC.O2<- gamma*(1-asymp_wC.O)*(1-ihr_wC.O)*EwC.O+ (1-nu)*IwC.O
      AwC.O2 <- gamma*asymp_wC.O*(1-ihr_wC.O)*EwC.O + (1-nu)*AwC.O
      HwC.O2 <- gamma*ihr_wC.O*EwC.O + (1-nus)*HwC.O
      RwC.O2 <- nu*IwC.O + nu*AwC.O + nus*(1-ihfr_wC.O)*HwC.O + RwC.O + VacwC.R.O*RvC.O
      DwC.O2 <- nus*ihfr_wC.O*HwC.O + DwC.O
      ###################################
      ##susceptibles
      Su <- Su2
      SvA <- SvA2
      SvP <- SvP2
      SvC <- SvC2
      SwA <- SwA2
      SwP <- SwP2
      SwC <- SwC2
      ##delta unvaccinated
      Eu.D <- Eu.D2
      Au.D <- Au.D2
      Iu.D <- Iu.D2
      Hu.D <- Hu.D2
      Ru.D <- Ru.D2
      Du.D <- Du.D2
      ##astrazeneca 1 dose
      EvA.D <- EvA.D2
      AvA.D <- AvA.D2
      IvA.D <- IvA.D2
      HvA.D <- HvA.D2
      RvA.D <- RvA.D2
      DvA.D <- DvA.D2
      ##pfizer 1 dose
      EvP.D <- EvP.D2
      AvP.D <- AvP.D2
      IvP.D <- IvP.D2
      HvP.D <- HvP.D2
      RvP.D <- RvP.D2
      DvP.D <- DvP.D2
      ##coronavac 1 dose
      EvC.D <- EvC.D2
      AvC.D <- AvC.D2
      IvC.D <- IvC.D2
      HvC.D <- HvC.D2
      RvC.D <- RvC.D2
      DvC.D <- DvC.D2
      ##astrazeneca 2 dose
      EwA.D <- EwA.D2
      AwA.D <- AwA.D2
      IwA.D <- IwA.D2
      HwA.D <- HwA.D2
      RwA.D <- RwA.D2
      DwA.D <- DwA.D2
      ##pfizer 2 dose
      EwP.D <- EwP.D2
      AwP.D <- AwP.D2
      IwP.D <- IwP.D2
      HwP.D <- HwP.D2
      RwP.D <- RwP.D2
      DwP.D <- DwP.D2
      ##coronavac 2 dose
      EwC.D <- EwC.D2
      AwC.D <- AwC.D2
      IwC.D <- IwC.D2
      HwC.D <- HwC.D2
      RwC.D <- RwC.D2
      DwC.D <- DwC.D2
      ##omicron unvaccinated
      Eu.O <- Eu.O2
      Au.O <- Au.O2
      Iu.O <- Iu.O2
      Hu.O <- Hu.O2
      Ru.O <- Ru.O2
      Du.O <- Du.O2
      ##astrazeneca 1 dose
      EvA.O <- EvA.O2
      AvA.O <- AvA.O2
      IvA.O <- IvA.O2
      HvA.O <- HvA.O2
      RvA.O <- RvA.O2
      DvA.O <- DvA.O2
      ##pfizer 1 dose
      EvP.O <- EvP.O2
      AvP.O <- AvP.O2
      IvP.O <- IvP.O2
      HvP.O <- HvP.O2
      RvP.O <- RvP.O2
      DvP.O <- DvP.O2
      ##coronavac 1 dose
      EvC.O <- EvC.O2
      AvC.O <- AvC.O2
      IvC.O <- IvC.O2
      HvC.O <- HvC.O2
      RvC.O <- RvC.O2
      DvC.O <- DvC.O2
      ##astrazeneca 2 dose
      EwA.O <- EwA.O2
      AwA.O <- AwA.O2
      IwA.O <- IwA.O2
      HwA.O <- HwA.O2
      RwA.O <- RwA.O2
      DwA.O <- DwA.O2
      ##pfizer 2 dose
      EwP.O <- EwP.O2
      AwP.O <- AwP.O2
      IwP.O <- IwP.O2
      HwP.O <- HwP.O2
      RwP.O <- RwP.O2
      DwP.O <- DwP.O2
      ##coronavac 2 dose
      EwC.O <- EwC.O2
      AwC.O <- AwC.O2
      IwC.O <- IwC.O2
      HwC.O <- HwC.O2
      RwC.O <- RwC.O2
      DwC.O <- DwC.O2
      # Eu.O <- rep(0,9)
      # Au.O <- rep(0,9)
      # Iu.O <- rep(0,9)
      # Hu.O <- rep(0,9)
      # Ru.O <- rep(0,9)
      # Du.O <- rep(0,9)
      # EvA.O <- rep(0,9)
      # AvA.O <- rep(0,9)
      # IvA.O <- rep(0,9)
      # HvA.O <- rep(0,9)
      # RvA.O <- rep(0,9)
      # DvA.O <- rep(0,9)
      # EvP.O <- rep(0,9)
      # AvP.O <- rep(0,9)
      # IvP.O <- rep(0,9)
      # HvP.O <- rep(0,9)
      # RvP.O <- rep(0,9)
      # DvP.O <- rep(0,9)
      # EvC.O <- rep(0,9)
      # AvC.O <- rep(0,9)
      # IvC.O <- rep(0,9)
      # HvC.O <- rep(0,9)
      # RvC.O <- rep(0,9)
      # DvC.O <- rep(0,9)
      # EwA.O <- rep(0,9)
      # AwA.O <- rep(0,9)
      # IwA.O <- rep(0,9)
      # HwA.O <- rep(0,9)
      # RwA.O <- rep(0,9)
      # DwA.O <- rep(0,9)
      # EwP.O <- rep(0,9)
      # AwP.O <- rep(0,9)
      # IwP.O <- rep(0,9)
      # HwP.O <- rep(0,9)
      # RwP.O <- rep(0,9)
      # DwP.O <- rep(0,9)
      # EwC.O <- rep(0,9)
      # AwC.O <- rep(0,9)
      # IwC.O <- rep(0,9)
      # HwC.O <- rep(0,9)
      # RwC.O <- rep(0,9)
      # DwC.O <- rep(0,9)
      data[i,] <- Yn <- c(Su,SvA,SvP,SvC,SwA,SwP,SwC,
                              Eu.D,Au.D,Iu.D,Hu.D,Ru.D,Du.D,
                              EvA.D,AvA.D,IvA.D,HvA.D,RvA.D,DvA.D,
                              EvP.D,AvP.D,IvP.D,HvP.D,RvP.D,DvP.D,
                              EvC.D,AvC.D,IvC.D,HvC.D,RvC.D,DvC.D,
                              EwA.D,AwA.D,IwA.D,HwA.D,RwA.D,DwA.D,
                              EwP.D,AwP.D,IwP.D,HwP.D,RwP.D,DwP.D,
                              EwC.D,AwC.D,IwC.D,HwC.D,RwC.D,DwC.D,
                              Eu.O,Au.O,Iu.O,Hu.O,Ru.O,Du.O,
                              EvA.O,AvA.O,IvA.O,HvA.O,RvA.O,DvA.O,
                              EvP.O,AvP.O,IvP.O,HvP.O,RvP.O,DvP.O,
                              EvC.O,AvC.O,IvC.O,HvC.O,RvC.O,DvC.O,
                              EwA.O,AwA.O,IwA.O,HwA.O,RwA.O,DwA.O,
                              EwP.O,AwP.O,IwP.O,HwP.O,RwP.O,DwP.O,
                              EwC.O,AwC.O,IwC.O,HwC.O,RwC.O,DwC.O)
      print(sum(data[i,]))
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

