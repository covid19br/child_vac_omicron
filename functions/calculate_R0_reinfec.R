library(rARPACK)
library(Matrix)
library(matlib)
library(zeallot)
calculate_R0 <- function(parameters){
  ## Definir estados
  vac.states <- c("u","vA", "vP", "vC", "wA", "wP", "wC","bP")
  inf.states <- c("E", "A", "I", "H", "R", "D","C","X")
  # inf.states.2 <- c("E", "A", "I", "H", "R")
  var.states <- c(".D",".O")
  classes = c(paste0("S", vac.states),###susceptiveis vacinados ou não
              paste0(inf.states,rep(vac.states,each = length(inf.states)),
                     rep(var.states,each = length(inf.states)*length(vac.states))) #infectados vacinados ou não
  )
  Y <- parameters$init.condition
  names(Y) <- rep(classes,each=parameters$age.bins) 
  c(Su,SvA,SvP,SvC,SwA,SwP,SwC,SbP,
    Eu.D,Au.D,Iu.D,Hu.D,Ru.D,Du.D,Cu.D, Xu.D, 
    EvA.D,AvA.D,IvA.D,HvA.D,RvA.D,DvA.D,CvA.D, XvA.D,
    EvP.D,AvP.D,IvP.D,HvP.D,RvP.D,DvP.D,CvP.D, XvP.D,
    EvC.D,AvC.D,IvC.D,HvC.D,RvC.D,DvC.D,CvC.D, XvC.D,
    EwA.D,AwA.D,IwA.D,HwA.D,RwA.D,DwA.D,CwA.D, XwA.D,
    EwP.D,AwP.D,IwP.D,HwP.D,RwP.D,DwP.D,CwP.D, XwP.D,
    EwC.D,AwC.D,IwC.D,HwC.D,RwC.D,DwC.D,CwC.D, XwC.D,
    EbP.D,AbP.D,IbP.D,HbP.D,RbP.D,DbP.D,CbP.D, XbP.D,
    Eu.O,Au.O,Iu.O,Hu.O,Ru.O,Du.O,Cu.O, Xu.O,
    EvA.O,AvA.O,IvA.O,HvA.O,RvA.O,DvA.O,CvA.O, XvA.O,
    EvP.O,AvP.O,IvP.O,HvP.O,RvP.O,DvP.O,CvP.O, XvP.O,
    EvC.O,AvC.O,IvC.O,HvC.O,RvC.O,DvC.O,CvC.O, XvC.O,
    EwA.O,AwA.O,IwA.O,HwA.O,RwA.O,DwA.O,CwA.O, XwA.O,
    EwP.O,AwP.O,IwP.O,HwP.O,RwP.O,DwP.O,CwP.O, XwP.O,
    EwC.O,AwC.O,IwC.O,HwC.O,RwC.O,DwC.O,CwC.O, XwC.O,
    EbP.O,AbP.O,IbP.O,HbP.O,RbP.O,DbP.O,CbP.O, XbP.O) %<-% split(Y,factor(names(Y),levels= classes))
  N <- Su + SvA + SvP + SvC + SwA + SwP + SwC + SbP +
    Eu.D + Au.D + Iu.D + Hu.D + Ru.D + 
    EvA.D + AvA.D + IvA.D + HvA.D + RvA.D + 
    EvP.D + AvP.D + IvP.D + HvP.D + RvP.D + 
    EvC.D + AvC.D + IvC.D + HvC.D + RvC.D + 
    EwA.D + AwA.D + IwA.D + HwA.D + RwA.D + 
    EwP.D + AwP.D + IwP.D + HwP.D + RwP.D + 
    EwC.D + AwC.D + IwC.D + HwC.D + RwC.D + 
    EbP.D + AbP.D + IbP.D + HbP.D + RbP.D + 
    Eu.O + Au.O + Iu.O + Hu.O + Ru.O + 
    EvA.O + AvA.O + IvA.O + HvA.O + RvA.O + 
    EvP.O + AvP.O + IvP.O + HvP.O + RvP.O + 
    EvC.O + AvC.O + IvC.O + HvC.O + RvC.O + 
    EwA.O + AwA.O + IwA.O + HwA.O + RwA.O + 
    EwP.O + AwP.O + IwP.O + HwP.O + RwP.O + 
    EwC.O + AwC.O + IwC.O + HwC.O + RwC.O +
    EbP.O + AbP.O + IbP.O + HbP.O + RbP.O 
    
  R0 <- with(parameters,{
    #susceptibles
    contact.Su <- diag(Su)%*%(c%*%diag(1/N))
    contact.SvA <- diag(SvA)%*%(c%*%diag(1/N))
    contact.SvP <- diag(SvP)%*%(c%*%diag(1/N))
    contact.SvC <- diag(SvC)%*%(c%*%diag(1/N))
    contact.SwA <- diag(SwA)%*%(c%*%diag(1/N))
    contact.SwP <- diag(SwP)%*%(c%*%diag(1/N))
    contact.SwC <- diag(SwC)%*%(c%*%diag(1/N))
    contact.SbP <- diag(SbP)%*%(c%*%diag(1/N))
    #recovered delta
    contact.Ru.D <- diag(Ru.D)%*%(c%*%diag(1/N))
    contact.RvA.D <- diag(RvA.D)%*%(c%*%diag(1/N))
    contact.RvP.D <- diag(RvP.D)%*%(c%*%diag(1/N))
    contact.RvC.D <- diag(RvC.D)%*%(c%*%diag(1/N))
    contact.RwA.D <- diag(RwA.D)%*%(c%*%diag(1/N))
    contact.RwP.D <- diag(RwP.D)%*%(c%*%diag(1/N))
    contact.RwC.D <- diag(RwC.D)%*%(c%*%diag(1/N))
    contact.RbP.D <- diag(RbP.D)%*%(c%*%diag(1/N))
    #recovered omicron
    contact.Ru.O <- diag(Ru.O)%*%(c%*%diag(1/N))
    contact.RvA.O <- diag(RvA.O)%*%(c%*%diag(1/N))
    contact.RvP.O <- diag(RvP.O)%*%(c%*%diag(1/N))
    contact.RvC.O <- diag(RvC.O)%*%(c%*%diag(1/N))
    contact.RwA.O <- diag(RwA.O)%*%(c%*%diag(1/N))
    contact.RwP.O <- diag(RwP.O)%*%(c%*%diag(1/N))
    contact.RwC.O <- diag(RwC.O)%*%(c%*%diag(1/N))
    contact.RbP.O <- diag(RbP.O)%*%(c%*%diag(1/N))
    
    ############################################################################
    ############################### DELTA ######################################
    ############################################################################    
    
    # X_0 = (Eu.D, EvA.D, EvP.D, EvC.D, EwA.D, EwP.D, EwC.D, EbP.D
    #        Au.D, AvA.D, AvP.D, AvC.D, AwA.D, AwP.D, AwC.D, AbP.D
    #        Iu.D, IvA.D, IvP.D, IvC.D, IwA.D, IwP.D, IwC.D, IbP.D
    #        Hu.D, HvA.D, HvP.D, HvC.D, HwA.D, HwP.D, HwC.D, HbP.D)
    # X_1 = (Su, SvA, SvP, SvP, SvC, SwA, SwP, SwC, SbP
    #        Ru.D, RvA.D, RvP.D, RvC.D, RwA.D, RwP.D, RwC.D, RbP.D
    #        Du.D, DvA.D, DvP.D, DvC.D, DwA.D, DwP.D, DwC.D, DbP.D)
    
    ############################### F matrix ###################################
    # EuDxE cols
    F_Eu.DEu.D <- F_Eu.DEvA.D <- F_Eu.DEvP.D <- F_Eu.DEvC.D <- 
      F_Eu.DEwA.D <- F_Eu.DEwP.D <- F_Eu.DEwC.D <- F_Eu.DEbP.D <- beta.D*omega*(contact.Su + 
                                                                                 reinfec.DD*contact.Ru.D+
                                                                                 reinfec.OD*contact.Ru.O)
    # EuDxA cols
    F_Eu.DAu.D <- F_Eu.DAvA.D <- F_Eu.DAvP.D <- F_Eu.DAvC.D <- 
      F_Eu.DAwA.D <- F_Eu.DAwP.D <- F_Eu.DAwC.D <- F_Eu.DAbP.D <- beta.D*omega_a*(contact.Su + 
                                                                                  reinfec.DD*contact.Ru.D+
                                                                                  reinfec.OD*contact.Ru.O) 
    # EuDxI cols
    F_Eu.DIu.D <- F_Eu.DIvA.D <- F_Eu.DIvP.D <- F_Eu.DIvC.D <- 
      F_Eu.DIwA.D <- F_Eu.DIwP.D <- F_Eu.DIwC.D <- F_Eu.DIbP.D <- beta.D*(contact.Su + 
                                                                          reinfec.DD*contact.Ru.D+
                                                                          reinfec.OD*contact.Ru.O)
    # EuDxH cols
    F_Eu.DHu.D <- F_Eu.DHvA.D <- F_Eu.DHvP.D <- F_Eu.DHvC.D <- 
      F_Eu.DHwA.D <- F_Eu.DHwP.D <- F_Eu.DHwC.D <- F_Eu.DHbP.D <- beta.D*(contact.Su + 
                                                                          reinfec.DD*contact.Ru.D+
                                                                          reinfec.OD*contact.Ru.O)
    # EuD lines
    F_Eu.D <- cbind(F_Eu.DEu.D, F_Eu.DEvA.D, F_Eu.DEvP.D, F_Eu.DEvC.D, F_Eu.DEwA.D, F_Eu.DEwP.D, F_Eu.DEwC.D, F_Eu.DEbP.D,
                    F_Eu.DAu.D, F_Eu.DAvA.D, F_Eu.DAvP.D, F_Eu.DAvC.D, F_Eu.DAwA.D, F_Eu.DAwP.D, F_Eu.DAwC.D, F_Eu.DAbP.D,
                    F_Eu.DIu.D, F_Eu.DIvA.D, F_Eu.DIvP.D, F_Eu.DIvC.D, F_Eu.DIwA.D, F_Eu.DIwP.D, F_Eu.DIwC.D, F_Eu.DIbP.D,
                    F_Eu.DHu.D, F_Eu.DHvA.D, F_Eu.DHvP.D, F_Eu.DHvC.D, F_Eu.DHwA.D, F_Eu.DHwP.D, F_Eu.DHwC.D, F_Eu.DHbP.D)
    
    # EvA.DxE cols
    F_EvA.DEu.D <- F_EvA.DEvA.D <- F_EvA.DEvP.D <- F_EvA.DEvC.D <- 
      F_EvA.DEwA.D <- F_EvA.DEwP.D <- F_EvA.DEwC.D <- F_EvA.DEbP.D <- beta_vA.D*omega*(contact.SvA + 
                                                                                       reinfec.DD*contact.RvA.D+
                                                                                       reinfec.OD*contact.RvA.O)
    #EvA.DxA cols    
    F_EvA.DAu.D <- F_EvA.DAvA.D <- F_EvA.DAvP.D <- F_EvA.DAvC.D <- 
      F_EvA.DAwA.D <- F_EvA.DAwP.D <- F_EvA.DAwC.D <- F_EvA.DAbP.D <- beta_vA.D*omega_a*(contact.SvA + 
                                                                                         reinfec.DD*contact.RvA.D+
                                                                                         reinfec.OD*contact.RvA.O)
    # EvA.DxI cols
    F_EvA.DIu.D <- F_EvA.DIvA.D <- F_EvA.DIvP.D <- F_EvA.DIvC.D <- 
      F_EvA.DIwA.D <- F_EvA.DIwP.D <- F_EvA.DIwC.D <- F_EvA.DIbP.D <- beta_vA.D*(contact.SvA + 
                                                                                 reinfec.DD*contact.RvA.D+
                                                                                 reinfec.OD*contact.RvA.O)
    #EvA.DxH cols
    F_EvA.DHu.D <- F_EvA.DHvA.D <- F_EvA.DHvP.D <- F_EvA.DHvC.D <- 
      F_EvA.DHwA.D <- F_EvA.DHwP.D <- F_EvA.DHwC.D <- F_EvA.DHbP.D <- beta_vA.D*omega_s*(contact.SvA + 
                                                                                           reinfec.DD*contact.RvA.D+
                                                                                           reinfec.OD*contact.RvA.O)
    #F_EvA.D lines
    F_EvA.D <- cbind(F_EvA.DEu.D, F_EvA.DEvA.D, F_EvA.DEvP.D, F_EvA.DEvC.D, F_EvA.DEwA.D, F_EvA.DEwP.D, F_EvA.DEwC.D, F_EvA.DEbP.D,
                     F_EvA.DAu.D, F_EvA.DAvA.D, F_EvA.DAvP.D, F_EvA.DAvC.D, F_EvA.DAwA.D, F_EvA.DAwP.D, F_EvA.DAwC.D, F_EvA.DAbP.D,
                     F_EvA.DIu.D, F_EvA.DIvA.D, F_EvA.DIvP.D, F_EvA.DIvC.D, F_EvA.DIwA.D, F_EvA.DIwP.D, F_EvA.DIwC.D, F_EvA.DIbP.D,
                     F_EvA.DHu.D, F_EvA.DHvA.D, F_EvA.DHvP.D, F_EvA.DHvC.D, F_EvA.DHwA.D, F_EvA.DHwP.D, F_EvA.DHwC.D, F_EvA.DHbP.D)
    
    #EvP.DxE cols
    F_EvP.DEu.D <- F_EvP.DEvA.D <- F_EvP.DEvP.D <- F_EvP.DEvC.D <- 
      F_EvP.DEwA.D <- F_EvP.DEwP.D <- F_EvP.DEwC.D <- F_EvP.DEbP.D <- beta_vP.D*omega*(contact.SvP + 
                                                                                         reinfec.DD*contact.RvP.D+
                                                                                         reinfec.OD*contact.RvP.O)
    #EvP.DxA cols
    F_EvP.DAu.D <- F_EvP.DAvA.D <- F_EvP.DAvP.D <- F_EvP.DAvC.D <- 
      F_EvP.DAwA.D <- F_EvP.DAwP.D <- F_EvP.DAwC.D <- F_EvP.DAbP.D <- beta_vP.D*omega_a*(contact.SvP + 
                                                                                           reinfec.DD*contact.RvP.D+
                                                                                           reinfec.OD*contact.RvP.O)
    #EvP.DxI cols
    F_EvP.DIu.D <- F_EvP.DIvA.D <- F_EvP.DIvP.D <- F_EvP.DIvC.D <- 
      F_EvP.DIwA.D <- F_EvP.DIwP.D <- F_EvP.DIwC.D <- F_EvP.DIbP.D <- beta_vP.D*(contact.SvP + 
                                                                                   reinfec.DD*contact.RvP.D+
                                                                                   reinfec.OD*contact.RvP.O)
    #EvP.DxH cols
    F_EvP.DHu.D <- F_EvP.DHvA.D <- F_EvP.DHvP.D <- F_EvP.DHvC.D <- 
      F_EvP.DHwA.D <- F_EvP.DHwP.D <- F_EvP.DHwC.D <- F_EvP.DHbP.D <- beta_vP.D*omega_s*(contact.SvP + 
                                                                                           reinfec.DD*contact.RvP.D+
                                                                                           reinfec.OD*contact.RvP.O)
    #EvP.D lines
    F_EvP.D <- cbind(F_EvP.DEu.D, F_EvP.DEvA.D, F_EvP.DEvP.D, F_EvP.DEvC.D, F_EvP.DEwA.D, F_EvP.DEwP.D, F_EvP.DEwC.D, F_EvP.DEbP.D,
                     F_EvP.DAu.D, F_EvP.DAvA.D, F_EvP.DAvP.D, F_EvP.DAvC.D, F_EvP.DAwA.D, F_EvP.DAwP.D, F_EvP.DAwC.D, F_EvP.DAbP.D,
                     F_EvP.DIu.D, F_EvP.DIvA.D, F_EvP.DIvP.D, F_EvP.DIvC.D, F_EvP.DIwA.D, F_EvP.DIwP.D, F_EvP.DIwC.D, F_EvP.DIbP.D,
                     F_EvP.DHu.D, F_EvP.DHvA.D, F_EvP.DHvP.D, F_EvP.DHvC.D, F_EvP.DHwA.D, F_EvP.DHwP.D, F_EvP.DHwC.D, F_EvP.DHbP.D)
    #EvC.DxE cols
    F_EvC.DEu.D <- F_EvC.DEvA.D <- F_EvC.DEvP.D <- F_EvC.DEvC.D <- 
      F_EvC.DEwA.D <- F_EvC.DEwP.D <- F_EvC.DEwC.D <- F_EvC.DEbP.D <- beta_vC.D*omega*(contact.SvC + 
                                                                                         reinfec.DD*contact.RvC.D+
                                                                                         reinfec.OD*contact.RvC.O)
    #EvC.DxA cols
    F_EvC.DAu.D <- F_EvC.DAvA.D <- F_EvC.DAvP.D <- F_EvC.DAvC.D <- 
      F_EvC.DAwA.D <- F_EvC.DAwP.D <- F_EvC.DAwC.D <- F_EvC.DAbP.D <- beta_vC.D*omega_a*(contact.SvC + 
                                                                                           reinfec.DD*contact.RvC.D+
                                                                                           reinfec.OD*contact.RvC.O)
    #EvC.DxI cols
    F_EvC.DIu.D <- F_EvC.DIvA.D <- F_EvC.DIvP.D <- F_EvC.DIvC.D <- 
      F_EvC.DIwA.D <- F_EvC.DIwP.D <- F_EvC.DIwC.D <- F_EvC.DIbP.D <- beta_vC.D*(contact.SvC + 
                                                                                   reinfec.DD*contact.RvC.D+
                                                                                   reinfec.OD*contact.RvC.O)
    #EvC.DxH cols
    F_EvC.DHu.D <- F_EvC.DHvA.D <- F_EvC.DHvP.D <- F_EvC.DHvC.D <- 
      F_EvC.DHwA.D <- F_EvC.DHwP.D <- F_EvC.DHwC.D <- F_EvC.DHbP.D <- beta_vC.D*omega_s*(contact.SvC + 
                                                                                           reinfec.DD*contact.RvC.D+
                                                                                           reinfec.OD*contact.RvC.O)
    #EvC.D lines
    F_EvC.D <- cbind(F_EvC.DEu.D, F_EvC.DEvA.D, F_EvC.DEvP.D, F_EvC.DEvC.D, F_EvC.DEwA.D, F_EvC.DEwP.D, F_EvC.DEwC.D, F_EvC.DEbP.D,
                     F_EvC.DAu.D, F_EvC.DAvA.D, F_EvC.DAvP.D, F_EvC.DAvC.D, F_EvC.DAwA.D, F_EvC.DAwP.D, F_EvC.DAwC.D, F_EvC.DAbP.D,
                     F_EvC.DIu.D, F_EvC.DIvA.D, F_EvC.DIvP.D, F_EvC.DIvC.D, F_EvC.DIwA.D, F_EvC.DIwP.D, F_EvC.DIwC.D, F_EvC.DIbP.D,
                     F_EvC.DHu.D, F_EvC.DHvA.D, F_EvC.DHvP.D, F_EvC.DHvC.D, F_EvC.DHwA.D, F_EvC.DHwP.D, F_EvC.DHwC.D, F_EvC.DHbP.D)
    
    #F_EwA.DxE cols
    F_EwA.DEu.D <- F_EwA.DEvA.D <- F_EwA.DEvP.D <- F_EwA.DEvC.D <- 
      F_EwA.DEwA.D <- F_EwA.DEwP.D <- F_EwA.DEwC.D <- F_EwA.DEbP.D <- beta_wA.D*omega*(contact.SwA + 
                                                                                         reinfec.DD*contact.RwA.D+
                                                                                         reinfec.OD*contact.RwA.O)
    #F_EwA.DxA cols
    F_EwA.DAu.D <- F_EwA.DAvA.D <- F_EwA.DAvP.D <- F_EwA.DAvC.D <- 
      F_EwA.DAwA.D <- F_EwA.DAwP.D <- F_EwA.DAwC.D <- F_EwA.DAbP.D <- beta_wA.D*omega_a*(contact.SwA + 
                                                                                           reinfec.DD*contact.RwA.D+
                                                                                           reinfec.OD*contact.RwA.O)
    #F_EwA.DxI cols
    F_EwA.DIu.D <- F_EwA.DIvA.D <- F_EwA.DIvP.D <- F_EwA.DIvC.D <- 
      F_EwA.DIwA.D <- F_EwA.DIwP.D <- F_EwA.DIwC.D <- F_EwA.DIbP.D <- beta_wA.D*(contact.SwA + 
                                                                                   reinfec.DD*contact.RwA.D+
                                                                                   reinfec.OD*contact.RwA.O)
    #F_EwA.DxH cols
    F_EwA.DHu.D <- F_EwA.DHvA.D <- F_EwA.DHvP.D <- F_EwA.DHvC.D <- 
      F_EwA.DHwA.D <- F_EwA.DHwP.D <- F_EwA.DHwC.D <- F_EwA.DHbP.D <- beta_wA.D*omega_s*(contact.SwA + 
                                                                                           reinfec.DD*contact.RwA.D+
                                                                                           reinfec.OD*contact.RwA.O)
    #F_EwA.D lines
    F_EwA.D <- cbind(F_EwA.DEu.D, F_EwA.DEvA.D, F_EwA.DEvP.D, F_EwA.DEvC.D, F_EwA.DEwA.D, F_EwA.DEwP.D, F_EwA.DEwC.D, F_EwA.DEbP.D,
                     F_EwA.DAu.D, F_EwA.DAvA.D, F_EwA.DAvP.D, F_EwA.DAvC.D, F_EwA.DAwA.D, F_EwA.DAwP.D, F_EwA.DAwC.D, F_EwA.DAbP.D,
                     F_EwA.DIu.D, F_EwA.DIvA.D, F_EwA.DIvP.D, F_EwA.DIvC.D, F_EwA.DIwA.D, F_EwA.DIwP.D, F_EwA.DIwC.D, F_EwA.DIbP.D,
                     F_EwA.DHu.D, F_EwA.DHvA.D, F_EwA.DHvP.D, F_EwA.DHvC.D, F_EwA.DHwA.D, F_EwA.DHwP.D, F_EwA.DHwC.D, F_EwA.DHbP.D)
    
    #F_EwP.DxE cols
    F_EwP.DEu.D <- F_EwP.DEvA.D <- F_EwP.DEvP.D <- F_EwP.DEvC.D <- 
      F_EwP.DEwA.D <- F_EwP.DEwP.D <- F_EwP.DEwC.D <- F_EwP.DEbP.D <- beta_wP.D*omega*(contact.SwP + 
                                                                                         reinfec.DD*contact.RwP.D+
                                                                                         reinfec.OD*contact.RwP.O)
    #F_EwP.DxA cols
    F_EwP.DAu.D <- F_EwP.DAvA.D <- F_EwP.DAvP.D <- F_EwP.DAvC.D <- 
      F_EwP.DAwA.D <- F_EwP.DAwP.D <- F_EwP.DAwC.D <- F_EwP.DAbP.D <- beta_wP.D*omega_a*(contact.SwP + 
                                                                                           reinfec.DD*contact.RwP.D+
                                                                                           reinfec.OD*contact.RwP.O)
    #F_EwP.DxI cols
    F_EwP.DIu.D <- F_EwP.DIvA.D <- F_EwP.DIvP.D <- F_EwP.DIvC.D <- 
      F_EwP.DIwA.D <- F_EwP.DIwP.D <- F_EwP.DIwC.D <- F_EwP.DIbP.D <- beta_wP.D*(contact.SwP + 
                                                                                   reinfec.DD*contact.RwP.D+
                                                                                   reinfec.OD*contact.RwP.O)
    #F_EwP.DxH cols
    F_EwP.DHu.D <- F_EwP.DHvA.D <- F_EwP.DHvP.D <- F_EwP.DHvC.D <- 
      F_EwP.DHwA.D <- F_EwP.DHwP.D <- F_EwP.DHwC.D <- F_EwP.DHbP.D <- beta_wP.D*omega_s*(contact.SwP + 
                                                                                           reinfec.DD*contact.RwP.D+
                                                                                           reinfec.OD*contact.RwP.O)
    #F_EwP.D lines
    F_EwP.D <- cbind(F_EwP.DEu.D, F_EwP.DEvA.D, F_EwP.DEvP.D, F_EwP.DEvC.D, F_EwP.DEwA.D, F_EwP.DEwP.D, F_EwP.DEwC.D, F_EwP.DEbP.D,
                     F_EwP.DAu.D, F_EwP.DAvA.D, F_EwP.DAvP.D, F_EwP.DAvC.D, F_EwP.DAwA.D, F_EwP.DAwP.D, F_EwP.DAwC.D, F_EwP.DAbP.D,
                     F_EwP.DIu.D, F_EwP.DIvA.D, F_EwP.DIvP.D, F_EwP.DIvC.D, F_EwP.DIwA.D, F_EwP.DIwP.D, F_EwP.DIwC.D, F_EwP.DIbP.D,
                     F_EwP.DHu.D, F_EwP.DHvA.D, F_EwP.DHvP.D, F_EwP.DHvC.D, F_EwP.DHwA.D, F_EwP.DHwP.D, F_EwP.DHwC.D, F_EwP.DHbP.D)
    
    #F_EwC.DxE cols
    F_EwC.DEu.D <- F_EwC.DEvA.D <- F_EwC.DEvP.D <- F_EwC.DEvC.D <- 
      F_EwC.DEwA.D <- F_EwC.DEwP.D <- F_EwC.DEwC.D <- F_EwC.DEbP.D <- beta_wC.D*omega*(contact.SwC + 
                                                                                         reinfec.DD*contact.RwC.D+
                                                                                         reinfec.OD*contact.RwC.O)
    #F_EwC.DxA cols
    F_EwC.DAu.D <- F_EwC.DAvA.D <- F_EwC.DAvP.D <- F_EwC.DAvC.D <- 
      F_EwC.DAwA.D <- F_EwC.DAwP.D <- F_EwC.DAwC.D <- F_EwC.DAbP.D <- beta_wC.D*omega_a*(contact.SwC + 
                                                                                           reinfec.DD*contact.RwC.D+
                                                                                           reinfec.OD*contact.RwC.O)
    #F_EwC.DxI cols
    F_EwC.DIu.D <- F_EwC.DIvA.D <- F_EwC.DIvP.D <- F_EwC.DIvC.D <- 
      F_EwC.DIwA.D <- F_EwC.DIwP.D <- F_EwC.DIwC.D <- F_EwC.DIbP.D <- beta_wC.D*(contact.SwC + 
                                                                                   reinfec.DD*contact.RwC.D+
                                                                                   reinfec.OD*contact.RwC.O)
    #F_EwC.DxH cols
    F_EwC.DHu.D <- F_EwC.DHvA.D <- F_EwC.DHvP.D <- F_EwC.DHvC.D <- 
      F_EwC.DHwA.D <- F_EwC.DHwP.D <- F_EwC.DHwC.D <- F_EwC.DHbP.D <- beta_wC.D*omega_s*(contact.SwC + 
                                                                                           reinfec.DD*contact.RwC.D+
                                                                                           reinfec.OD*contact.RwC.O)
    #F_EwC.D lines
    F_EwC.D <- cbind(F_EwC.DEu.D, F_EwC.DEvA.D, F_EwC.DEvP.D, F_EwC.DEvC.D, F_EwC.DEwA.D, F_EwC.DEwP.D, F_EwC.DEwC.D, F_EwC.DEbP.D,
                     F_EwC.DAu.D, F_EwC.DAvA.D, F_EwC.DAvP.D, F_EwC.DAvC.D, F_EwC.DAwA.D, F_EwC.DAwP.D, F_EwC.DAwC.D, F_EwC.DAbP.D,
                     F_EwC.DIu.D, F_EwC.DIvA.D, F_EwC.DIvP.D, F_EwC.DIvC.D, F_EwC.DIwA.D, F_EwC.DIwP.D, F_EwC.DIwC.D, F_EwC.DIbP.D,
                     F_EwC.DHu.D, F_EwC.DHvA.D, F_EwC.DHvP.D, F_EwC.DHvC.D, F_EwC.DHwA.D, F_EwC.DHwP.D, F_EwC.DHwC.D, F_EwC.DHbP.D)
    
    #F_EwC.DxE cols
    F_EbP.DEu.D <- F_EbP.DEvA.D <- F_EbP.DEvP.D <- F_EbP.DEvC.D <- 
      F_EbP.DEwA.D <- F_EbP.DEwP.D <- F_EbP.DEwC.D <- F_EbP.DEbP.D <- beta_bP.D*omega*(contact.SbP + 
                                                                                         reinfec.DD*contact.RbP.D+
                                                                                         reinfec.OD*contact.RbP.O)
    #F_EbP.DxA cols
    F_EbP.DAu.D <- F_EbP.DAvA.D <- F_EbP.DAvP.D <- F_EbP.DAvC.D <- 
      F_EbP.DAwA.D <- F_EbP.DAwP.D <- F_EbP.DAwC.D <- F_EbP.DAbP.D <- beta_bP.D*omega_a*(contact.SbP + 
                                                                                           reinfec.DD*contact.RbP.D+
                                                                                           reinfec.OD*contact.RbP.O)
    #F_EbP.DxI cols
    F_EbP.DIu.D <- F_EbP.DIvA.D <- F_EbP.DIvP.D <- F_EbP.DIvC.D <- 
      F_EbP.DIwA.D <- F_EbP.DIwP.D <- F_EbP.DIwC.D <- F_EbP.DIbP.D <- beta_bP.D*(contact.SbP + 
                                                                                   reinfec.DD*contact.RbP.D+
                                                                                   reinfec.OD*contact.RbP.O)
    #F_EbP.DxH cols
    F_EbP.DHu.D <- F_EbP.DHvA.D <- F_EbP.DHvP.D <- F_EbP.DHvC.D <- 
      F_EbP.DHwA.D <- F_EbP.DHwP.D <- F_EbP.DHwC.D <- F_EbP.DHbP.D <- beta_bP.D*omega_s*(contact.SbP + 
                                                                                           reinfec.DD*contact.RbP.D+
                                                                                           reinfec.OD*contact.RbP.O)
    #F_EbP.D lines
    F_EbP.D <- cbind(F_EbP.DEu.D, F_EbP.DEvA.D, F_EbP.DEvP.D, F_EbP.DEvC.D, F_EbP.DEwA.D, F_EbP.DEwP.D, F_EbP.DEwC.D, F_EbP.DEbP.D,
                     F_EbP.DAu.D, F_EbP.DAvA.D, F_EbP.DAvP.D, F_EbP.DAvC.D, F_EbP.DAwA.D, F_EbP.DAwP.D, F_EbP.DAwC.D, F_EbP.DAbP.D,
                     F_EbP.DIu.D, F_EbP.DIvA.D, F_EbP.DIvP.D, F_EbP.DIvC.D, F_EbP.DIwA.D, F_EbP.DIwP.D, F_EbP.DIwC.D, F_EbP.DIbP.D,
                     F_EbP.DHu.D, F_EbP.DHvA.D, F_EbP.DHvP.D, F_EbP.DHvC.D, F_EbP.DHwA.D, F_EbP.DHwP.D, F_EbP.DHwC.D, F_EbP.DHbP.D)
    #all E's together
    F_E.D <- rbind(F_Eu.D, F_EvA.D, F_EvP.D, F_EvC.D, F_EwA.D, F_EwP.D, F_EwC.D, F_EbP.D)
    # other lines are zero
    F_mat.D <- rbind(F_E.D, matrix(0, nrow = 3*nrow(F_E.D),ncol = ncol(F_E.D)))
    
    ################################ T matrix ################################## 
    ## Eu
    # Eu.DxE cols
    T_Eu.DEu.D <- diag(1-gamma.D, nrow = age.bins, ncol = age.bins)
    T_Eu.DEvA.D <- T_Eu.DEvP.D <- T_Eu.DEvC.D <- 
      T_Eu.DEwA.D <- T_Eu.DEwP.D <- T_Eu.DEwC.D <- T_Eu.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EuDxI cols
    T_Eu.DIu.D <- T_Eu.DIvA.D <- T_Eu.DIvP.D <- T_Eu.DIvC.D <- 
      T_Eu.DIwA.D <- T_Eu.DIwP.D <- T_Eu.DIwC.D <- T_Eu.DIbP.D  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EuDxA cols
    T_Eu.DAu.D <- T_Eu.DAvA.D <- T_Eu.DAvP.D <- T_Eu.DAvC.D <- 
      T_Eu.DAwA.D <- T_Eu.DAwP.D <- T_Eu.DAwC.D <- T_Eu.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EuDxH cols
    T_Eu.DHu.D <- T_Eu.DHvA.D <- T_Eu.DHvP.D <- T_Eu.DHvC.D <- 
      T_Eu.DHwA.D <- T_Eu.DHwP.D <- T_Eu.DHwC.D <- T_Eu.DHbP.D  <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_Eu.D <- cbind(T_Eu.DEu.D, T_Eu.DEvA.D, T_Eu.DEvP.D, T_Eu.DEvC.D, T_Eu.DEwA.D, T_Eu.DEwP.D, T_Eu.DEwC.D, T_Eu.DEbP.D,
                    T_Eu.DAu.D, T_Eu.DAvA.D, T_Eu.DAvP.D, T_Eu.DAvC.D, T_Eu.DAwA.D, T_Eu.DAwP.D, T_Eu.DAwC.D, T_Eu.DAbP.D,
                    T_Eu.DIu.D, T_Eu.DIvA.D, T_Eu.DIvP.D, T_Eu.DIvC.D, T_Eu.DIwA.D, T_Eu.DIwP.D, T_Eu.DIwC.D, T_Eu.DIbP.D,
                    T_Eu.DHu.D, T_Eu.DHvA.D, T_Eu.DHvP.D, T_Eu.DHvC.D, T_Eu.DHwA.D, T_Eu.DHwP.D, T_Eu.DHwC.D, T_Eu.DHbP.D)
    ## EvA
    # EvA.DxE cols
    T_EvA.DEvA.D <- diag(1-gamma.D, nrow = age.bins, ncol = age.bins)
    T_EvA.DEu.D <- T_EvA.DEvP.D <- T_EvA.DEvC.D <- 
      T_EvA.DEwA.D <- T_EvA.DEwP.D <- T_EvA.DEwC.D <- T_EvA.DEbP.D  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvA.DxA cols
    T_EvA.DAu.D <- T_EvA.DAvA.D <- T_EvA.DAvP.D <- T_EvA.DAvC.D <- 
      T_EvA.DAwA.D <- T_EvA.DAwP.D <- T_EvA.DAwC.D <- T_EvA.DAbP.D  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvA.DxI cols
    T_EvA.DIu.D <- T_EvA.DIvA.D <- T_EvA.DIvP.D <- T_EvA.DIvC.D <- 
      T_EvA.DIwA.D <- T_EvA.DIwP.D <- T_EvA.DIwC.D <- T_EvA.DIbP.D  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvAxH cols
    T_EvA.DHu.D <- T_EvA.DHvA.D <- T_EvA.DHvP.D <- T_EvA.DHvC.D <- 
      T_EvA.DHwA.D <- T_EvA.DHwP.D <- T_EvA.DHwC.D <- T_EvA.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EvA.D <- cbind(T_EvA.DEu.D, T_EvA.DEvA.D, T_EvA.DEvP.D, T_EvA.DEvC.D, T_EvA.DEwA.D, T_EvA.DEwP.D, T_EvA.DEwC.D, T_EvA.DEbP.D,
                     T_EvA.DAu.D, T_EvA.DAvA.D, T_EvA.DAvP.D, T_EvA.DAvC.D, T_EvA.DAwA.D, T_EvA.DAwP.D, T_EvA.DAwC.D, T_EvA.DAbP.D,
                     T_EvA.DIu.D, T_EvA.DIvA.D, T_EvA.DIvP.D, T_EvA.DIvC.D, T_EvA.DIwA.D, T_EvA.DIwP.D, T_EvA.DIwC.D, T_EvA.DIbP.D,
                     T_EvA.DHu.D, T_EvA.DHvA.D, T_EvA.DHvP.D, T_EvA.DHvC.D, T_EvA.DHwA.D, T_EvA.DHwP.D, T_EvA.DHwC.D, T_EvA.DHbP.D)
    
    ## EvP
    # EvP.DxE cols
    T_EvP.DEvP.D <- diag(1-gamma.D, nrow = age.bins, ncol = age.bins)
    T_EvP.DEu.D <- T_EvP.DEvA.D <- T_EvP.DEvC.D <- 
      T_EvP.DEwA.D <- T_EvP.DEwP.D <- T_EvP.DEwC.D <- T_EvP.DEbP.D  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvP.DxA cols
    T_EvP.DAu.D <- T_EvP.DAvA.D <- T_EvP.DAvP.D <- T_EvP.DAvC.D <- 
      T_EvP.DAwA.D <- T_EvP.DAwP.D <- T_EvP.DAwC.D <- T_EvP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvP.DxI cols
    T_EvP.DIu.D <- T_EvP.DIvA.D <- T_EvP.DIvP.D <- T_EvP.DIvC.D <- 
      T_EvP.DIwA.D <- T_EvP.DIwP.D <- T_EvP.DIwC.D <- T_EvP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvPxH cols
    T_EvP.DHu.D <- T_EvP.DHvA.D <- T_EvP.DHvP.D <- T_EvP.DHvC.D <- 
      T_EvP.DHwA.D <- T_EvP.DHwP.D <- T_EvP.DHwC.D <- T_EvP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EvP.D <- cbind(T_EvP.DEu.D, T_EvP.DEvA.D, T_EvP.DEvP.D, T_EvP.DEvC.D, T_EvP.DEwA.D, T_EvP.DEwP.D, T_EvP.DEwC.D, T_EvP.DEbP.D,
                     T_EvP.DAu.D, T_EvP.DAvA.D, T_EvP.DAvP.D, T_EvP.DAvC.D, T_EvP.DAwA.D, T_EvP.DAwP.D, T_EvP.DAwC.D, T_EvP.DAbP.D,
                     T_EvP.DIu.D, T_EvP.DIvA.D, T_EvP.DIvP.D, T_EvP.DIvC.D, T_EvP.DIwA.D, T_EvP.DIwP.D, T_EvP.DIwC.D, T_EvP.DIbP.D,
                     T_EvP.DHu.D, T_EvP.DHvA.D, T_EvP.DHvP.D, T_EvP.DHvC.D, T_EvP.DHwA.D, T_EvP.DHwP.D, T_EvP.DHwC.D, T_EvP.DHbP.D)
    
    ## EvC
    # EvC.DxE cols
    T_EvC.DEvC.D <- diag(1-gamma.D, nrow = age.bins, ncol = age.bins)
    T_EvC.DEu.D <- T_EvC.DEvA.D <- T_EvC.DEvP.D <- 
      T_EvC.DEwA.D <- T_EvC.DEwP.D <- T_EvC.DEwC.D <- T_EvC.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvC.DxA cols
    T_EvC.DAu.D <- T_EvC.DAvA.D <- T_EvC.DAvP.D <- T_EvC.DAvC.D <- 
      T_EvC.DAwA.D <- T_EvC.DAwP.D <- T_EvC.DAwC.D <- T_EvC.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvC.DxI cols
    T_EvC.DIu.D <- T_EvC.DIvA.D <- T_EvC.DIvP.D <- T_EvC.DIvC.D <- 
      T_EvC.DIwA.D <- T_EvC.DIwP.D <- T_EvC.DIwC.D <- T_EvC.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvCxH cols
    T_EvC.DHu.D <- T_EvC.DHvA.D <- T_EvC.DHvP.D <- T_EvC.DHvC.D <- 
      T_EvC.DHwA.D <- T_EvC.DHwP.D <- T_EvC.DHwC.D <- T_EvC.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EvC.D <- cbind(T_EvC.DEu.D, T_EvC.DEvA.D, T_EvC.DEvP.D, T_EvC.DEvC.D, T_EvC.DEwA.D, T_EvC.DEwP.D, T_EvC.DEwC.D, T_EvC.DEbP.D,
                     T_EvC.DAu.D, T_EvC.DAvA.D, T_EvC.DAvP.D, T_EvC.DAvC.D, T_EvC.DAwA.D, T_EvC.DAwP.D, T_EvC.DAwC.D, T_EvC.DAbP.D,
                     T_EvC.DIu.D, T_EvC.DIvA.D, T_EvC.DIvP.D, T_EvC.DIvC.D, T_EvC.DIwA.D, T_EvC.DIwP.D, T_EvC.DIwC.D, T_EvC.DIbP.D,
                     T_EvC.DHu.D, T_EvC.DHvA.D, T_EvC.DHvP.D, T_EvC.DHvC.D, T_EvC.DHwA.D, T_EvC.DHwP.D, T_EvC.DHwC.D, T_EvC.DHbP.D)
    
    ## EwA
    # EwA.DxE cols
    T_EwA.DEwA.D <- diag(1-gamma.D, nrow = age.bins, ncol = age.bins)
    T_EwA.DEu.D <- T_EwA.DEvA.D <- T_EwA.DEvP.D <- 
      T_EwA.DEvC.D <- T_EwA.DEwP.D <- T_EwA.DEwC.D <- T_EwA.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwA.DxA cols
    T_EwA.DAu.D <- T_EwA.DAvA.D <- T_EwA.DAvP.D <- T_EwA.DAvC.D <- 
      T_EwA.DAwA.D <- T_EwA.DAwP.D <- T_EwA.DAwC.D <- T_EwA.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwA.DxI cols
    T_EwA.DIu.D <- T_EwA.DIvA.D <- T_EwA.DIvP.D <- T_EwA.DIvC.D <- 
      T_EwA.DIwA.D <- T_EwA.DIwP.D <- T_EwA.DIwC.D <- T_EwA.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwAxH cols
    T_EwA.DHu.D <- T_EwA.DHvA.D <- T_EwA.DHvP.D <- T_EwA.DHvC.D <- 
      T_EwA.DHwA.D <- T_EwA.DHwP.D <- T_EwA.DHwC.D <- T_EwA.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EwA.D <- cbind(T_EwA.DEu.D, T_EwA.DEvA.D, T_EwA.DEvP.D, T_EwA.DEvC.D, T_EwA.DEwA.D, T_EwA.DEwP.D, T_EwA.DEwC.D, T_EwA.DEbP.D,
                     T_EwA.DAu.D, T_EwA.DAvA.D, T_EwA.DAvP.D, T_EwA.DAvC.D, T_EwA.DAwA.D, T_EwA.DAwP.D, T_EwA.DAwC.D, T_EwA.DAbP.D,
                     T_EwA.DIu.D, T_EwA.DIvA.D, T_EwA.DIvP.D, T_EwA.DIvC.D, T_EwA.DIwA.D, T_EwA.DIwP.D, T_EwA.DIwC.D, T_EwA.DIbP.D,
                     T_EwA.DHu.D, T_EwA.DHvA.D, T_EwA.DHvP.D, T_EwA.DHvC.D, T_EwA.DHwA.D, T_EwA.DHwP.D, T_EwA.DHwC.D, T_EwA.DHbP.D)
    
    ## EwP
    # EwP.DxE cols
    T_EwP.DEwP.D <- diag(1-gamma.D, nrow = age.bins, ncol = age.bins)
    T_EwP.DEu.D <- T_EwP.DEvA.D <- T_EwP.DEvP.D <- 
      T_EwP.DEvC.D <- T_EwP.DEwA.D <- T_EwP.DEwC.D <- T_EwP.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwP.DxA cols
    T_EwP.DAu.D <- T_EwP.DAvA.D <- T_EwP.DAvP.D <- T_EwP.DAvC.D <- 
      T_EwP.DAwA.D <- T_EwP.DAwP.D <- T_EwP.DAwC.D <- T_EwP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwP.DxI cols
    T_EwP.DIu.D <- T_EwP.DIvA.D <- T_EwP.DIvP.D <- T_EwP.DIvC.D <- 
      T_EwP.DIwA.D <- T_EwP.DIwP.D <- T_EwP.DIwC.D <- T_EwP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwPxH cols
    T_EwP.DHu.D <- T_EwP.DHvA.D <- T_EwP.DHvP.D <- T_EwP.DHvC.D <- 
      T_EwP.DHwA.D <- T_EwP.DHwP.D <- T_EwP.DHwC.D <- T_EwP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EwP.D <- cbind(T_EwP.DEu.D, T_EwP.DEvA.D, T_EwP.DEvP.D, T_EwP.DEvC.D, T_EwP.DEwA.D, T_EwP.DEwP.D, T_EwP.DEwC.D, T_EwP.DEbP.D,
                     T_EwP.DAu.D, T_EwP.DAvA.D, T_EwP.DAvP.D, T_EwP.DAvC.D, T_EwP.DAwA.D, T_EwP.DAwP.D, T_EwP.DAwC.D, T_EwP.DAbP.D,
                     T_EwP.DIu.D, T_EwP.DIvA.D, T_EwP.DIvP.D, T_EwP.DIvC.D, T_EwP.DIwA.D, T_EwP.DIwP.D, T_EwP.DIwC.D, T_EwP.DIbP.D,
                     T_EwP.DHu.D, T_EwP.DHvA.D, T_EwP.DHvP.D, T_EwP.DHvC.D, T_EwP.DHwA.D, T_EwP.DHwP.D, T_EwP.DHwC.D, T_EwP.DHbP.D)
    
    ## EwC
    # EwC.DxE cols
    T_EwC.DEwC.D <- diag(1-gamma.D, nrow = age.bins, ncol = age.bins)
    T_EwC.DEu.D <- T_EwC.DEvA.D <- T_EwC.DEvP.D <- 
      T_EwC.DEvC.D <- T_EwC.DEwA.D <- T_EwC.DEwP.D <- T_EwC.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.DxA cols
    T_EwC.DAu.D <- T_EwC.DAvA.D <- T_EwC.DAvP.D <- T_EwC.DAvC.D <- 
      T_EwC.DAwA.D <- T_EwC.DAwP.D <- T_EwC.DAwC.D <- T_EwC.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.DxI cols
    T_EwC.DIu.D <- T_EwC.DIvA.D <- T_EwC.DIvP.D <- T_EwC.DIvC.D <- 
      T_EwC.DIwA.D <- T_EwC.DIwP.D <- T_EwC.DIwC.D <- T_EwC.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.DxH cols
    T_EwC.DHu.D <- T_EwC.DHvA.D <- T_EwC.DHvP.D <- T_EwC.DHvC.D <- 
      T_EwC.DHwA.D <- T_EwC.DHwP.D <- T_EwC.DHwC.D <- T_EwC.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EwC.D <- cbind(T_EwC.DEu.D, T_EwC.DEvA.D, T_EwC.DEvP.D, T_EwC.DEvC.D, T_EwC.DEwA.D, T_EwC.DEwP.D, T_EwC.DEwC.D, T_EwC.DEbP.D,
                     T_EwC.DAu.D, T_EwC.DAvA.D, T_EwC.DAvP.D, T_EwC.DAvC.D, T_EwC.DAwA.D, T_EwC.DAwP.D, T_EwC.DAwC.D, T_EwC.DAbP.D,
                     T_EwC.DIu.D, T_EwC.DIvA.D, T_EwC.DIvP.D, T_EwC.DIvC.D, T_EwC.DIwA.D, T_EwC.DIwP.D, T_EwC.DIwC.D, T_EwC.DIbP.D,
                     T_EwC.DHu.D, T_EwC.DHvA.D, T_EwC.DHvP.D, T_EwC.DHvC.D, T_EwC.DHwA.D, T_EwC.DHwP.D, T_EwC.DHwC.D, T_EwC.DHbP.D)
    ## EwC
    # EwC.DxE cols
    T_EbP.DEbP.D <- diag(1-gamma.D, nrow = age.bins, ncol = age.bins)
    T_EbP.DEu.D <- T_EbP.DEvA.D <- T_EbP.DEvP.D <- 
      T_EbP.DEvC.D <- T_EbP.DEwA.D <- T_EbP.DEwP.D <- T_EbP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.DxA cols
    T_EbP.DAu.D <- T_EbP.DAvA.D <- T_EbP.DAvP.D <- T_EbP.DAvC.D <- 
      T_EbP.DAwA.D <- T_EbP.DAwP.D <- T_EbP.DAwC.D <- T_EbP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.DxI cols
    T_EbP.DIu.D <- T_EbP.DIvA.D <- T_EbP.DIvP.D <- T_EbP.DIvC.D <- 
      T_EbP.DIwA.D <- T_EbP.DIwP.D <- T_EbP.DIwC.D <- T_EbP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.DxH cols
    T_EbP.DHu.D <- T_EbP.DHvA.D <- T_EbP.DHvP.D <- T_EbP.DHvC.D <- 
      T_EbP.DHwA.D <- T_EbP.DHwP.D <- T_EbP.DHwC.D <- T_EbP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EbP.D <- cbind(T_EbP.DEu.D, T_EbP.DEvA.D, T_EbP.DEvP.D, T_EbP.DEvC.D, T_EbP.DEwA.D, T_EbP.DEwP.D, T_EbP.DEwC.D, T_EbP.DEbP.D,
                     T_EbP.DAu.D, T_EbP.DAvA.D, T_EbP.DAvP.D, T_EbP.DAvC.D, T_EbP.DAwA.D, T_EbP.DAwP.D, T_EbP.DAwC.D, T_EbP.DAbP.D,
                     T_EbP.DIu.D, T_EbP.DIvA.D, T_EbP.DIvP.D, T_EbP.DIvC.D, T_EbP.DIwA.D, T_EbP.DIwP.D, T_EbP.DIwC.D, T_EbP.DIbP.D,
                     T_EbP.DHu.D, T_EbP.DHvA.D, T_EbP.DHvP.D, T_EbP.DHvC.D, T_EbP.DHwA.D, T_EbP.DHwP.D, T_EbP.DHwC.D, T_EbP.DHbP.D)
    #all E's together
    T_E.D <- rbind(T_Eu.D, T_EvA.D, T_EvP.D, T_EvC.D, T_EwA.D, T_EwP.D, T_EwC.D, T_EbP.D)
    
    ## Au
    # Au.DxE cols
    T_Au.DEu.D <- diag((1-ihr.D)*gamma.D*asymp, ncol = age.bins, nrow = age.bins)
    T_Au.DEvA.D <- T_Au.DEvP.D <- 
      T_Au.DEvC.D <- T_Au.DEwA.D <- T_Au.DEwP.D <- T_Au.DEwC.D <- T_Au.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Au.DxA cols
    T_Au.DAu.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_Au.DAvA.D <- T_Au.DAvP.D <- T_Au.DAvC.D <- 
      T_Au.DAwA.D <- T_Au.DAwP.D <- T_Au.DAwC.D <- T_Au.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Au.DxI cols
    T_Au.DIu.D <- T_Au.DIvA.D <- T_Au.DIvP.D <- T_Au.DIvC.D <- 
      T_Au.DIwA.D <- T_Au.DIwP.D <- T_Au.DIwC.D <- T_Au.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Au.DxH cols
    T_Au.DHu.D <- T_Au.DHvA.D <- T_Au.DHvP.D <- T_Au.DHvC.D <- 
      T_Au.DHwA.D <- T_Au.DHwP.D <- T_Au.DHwC.D <- T_Au.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_Au.D <- cbind(T_Au.DEu.D, T_Au.DEvA.D, T_Au.DEvP.D, T_Au.DEvC.D, T_Au.DEwA.D, T_Au.DEwP.D, T_Au.DEwC.D, T_Au.DEbP.D,
                    T_Au.DAu.D, T_Au.DAvA.D, T_Au.DAvP.D, T_Au.DAvC.D, T_Au.DAwA.D, T_Au.DAwP.D, T_Au.DAwC.D, T_Au.DAbP.D,
                    T_Au.DIu.D, T_Au.DIvA.D, T_Au.DIvP.D, T_Au.DIvC.D, T_Au.DIwA.D, T_Au.DIwP.D, T_Au.DIwC.D, T_Au.DIbP.D,
                    T_Au.DHu.D, T_Au.DHvA.D, T_Au.DHvP.D, T_Au.DHvC.D, T_Au.DHwA.D, T_Au.DHwP.D, T_Au.DHwC.D, T_Au.DHbP.D)
    
    ## AvA
    # AvA.DxE cols
    T_AvA.DEvA.D <- diag((1-ihr_vA.D)*gamma.D*asymp_vA.D, ncol = age.bins, nrow = age.bins)
    T_AvA.DEu.D <- T_AvA.DEvP.D <- 
      T_AvA.DEvC.D <- T_AvA.DEwA.D <- T_AvA.DEwP.D <- T_AvA.DEwC.D <- T_AvA.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvA.DxA cols
    T_AvA.DAvA.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AvA.DAu.D <- T_AvA.DAvP.D <- T_AvA.DAvC.D <- 
      T_AvA.DAwA.D <- T_AvA.DAwP.D <- T_AvA.DAwC.D <- T_AvA.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvA.DxI cols
    T_AvA.DIu.D <- T_AvA.DIvA.D <- T_AvA.DIvP.D <- T_AvA.DIvC.D <- 
      T_AvA.DIwA.D <- T_AvA.DIwP.D <- T_AvA.DIwC.D <- T_AvA.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvA.DxH cols
    T_AvA.DHu.D <- T_AvA.DHvA.D <- T_AvA.DHvP.D <- T_AvA.DHvC.D <- 
      T_AvA.DHwA.D <- T_AvA.DHwP.D <- T_AvA.DHwC.D <- T_AvA.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AvA.D <- cbind(T_AvA.DEu.D, T_AvA.DEvA.D, T_AvA.DEvP.D, T_AvA.DEvC.D, T_AvA.DEwA.D, T_AvA.DEwP.D, T_AvA.DEwC.D, T_AvA.DEbP.D,
                     T_AvA.DAu.D, T_AvA.DAvA.D, T_AvA.DAvP.D, T_AvA.DAvC.D, T_AvA.DAwA.D, T_AvA.DAwP.D, T_AvA.DAwC.D, T_AvA.DAbP.D,
                     T_AvA.DIu.D, T_AvA.DIvA.D, T_AvA.DIvP.D, T_AvA.DIvC.D, T_AvA.DIwA.D, T_AvA.DIwP.D, T_AvA.DIwC.D, T_AvA.DIbP.D,
                     T_AvA.DHu.D, T_AvA.DHvA.D, T_AvA.DHvP.D, T_AvA.DHvC.D, T_AvA.DHwA.D, T_AvA.DHwP.D, T_AvA.DHwC.D, T_AvA.DHbP.D)
    
    ## AvP
    # AvP.DxE cols
    T_AvP.DEvP.D <- diag((1-ihr_vP.D)*gamma.D*asymp_vP.D, ncol = age.bins, nrow = age.bins)
    T_AvP.DEu.D <- T_AvP.DEvA.D <- 
      T_AvP.DEvC.D <- T_AvP.DEwA.D <- T_AvP.DEwP.D <- T_AvP.DEwC.D <- T_AvP.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvP.DxA cols
    T_AvP.DAvP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AvP.DAu.D <- T_AvP.DAvA.D <- T_AvP.DAvC.D <- 
      T_AvP.DAwA.D <- T_AvP.DAwP.D <- T_AvP.DAwC.D <- T_AvP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvP.DxI cols
    T_AvP.DIu.D <- T_AvP.DIvA.D <- T_AvP.DIvP.D <- T_AvP.DIvC.D <- 
      T_AvP.DIwA.D <- T_AvP.DIwP.D <- T_AvP.DIwC.D <- T_AvP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvP.DxH cols
    T_AvP.DHu.D <- T_AvP.DHvA.D <- T_AvP.DHvP.D <- T_AvP.DHvC.D <- 
      T_AvP.DHwA.D <- T_AvP.DHwP.D <- T_AvP.DHwC.D <- T_AvP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AvP.D <- cbind(T_AvP.DEu.D, T_AvP.DEvA.D, T_AvP.DEvP.D, T_AvP.DEvC.D, T_AvP.DEwA.D, T_AvP.DEwP.D, T_AvP.DEwC.D, T_AvP.DEbP.D,
                     T_AvP.DAu.D, T_AvP.DAvA.D, T_AvP.DAvP.D, T_AvP.DAvC.D, T_AvP.DAwA.D, T_AvP.DAwP.D, T_AvP.DAwC.D, T_AvP.DAbP.D,
                     T_AvP.DIu.D, T_AvP.DIvA.D, T_AvP.DIvP.D, T_AvP.DIvC.D, T_AvP.DIwA.D, T_AvP.DIwP.D, T_AvP.DIwC.D, T_AvP.DIbP.D,
                     T_AvP.DHu.D, T_AvP.DHvA.D, T_AvP.DHvP.D, T_AvP.DHvC.D, T_AvP.DHwA.D, T_AvP.DHwP.D, T_AvP.DHwC.D, T_AvP.DHbP.D)
    
    ## AvC
    # AvC.DxE cols
    T_AvC.DEvC.D <- diag((1-ihr_vC.D)*gamma.D*asymp_vC.D, ncol = age.bins, nrow = age.bins)
    T_AvC.DEu.D <- T_AvC.DEvA.D <- 
      T_AvC.DEvP.D <- T_AvC.DEwA.D <- T_AvC.DEwP.D <- T_AvC.DEwC.D <- T_AvC.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvC.DxA cols
    T_AvC.DAvC.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AvC.DAu.D <- T_AvC.DAvA.D <- T_AvC.DAvP.D <- 
      T_AvC.DAwA.D <- T_AvC.DAwP.D <- T_AvC.DAwC.D <- T_AvC.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvC.DxI cols
    T_AvC.DIu.D <- T_AvC.DIvA.D <- T_AvC.DIvP.D <- T_AvC.DIvC.D <- 
      T_AvC.DIwA.D <- T_AvC.DIwP.D <- T_AvC.DIwC.D <- T_AvC.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvC.DxH cols
    T_AvC.DHu.D <- T_AvC.DHvA.D <- T_AvC.DHvP.D <- T_AvC.DHvC.D <- 
      T_AvC.DHwA.D <- T_AvC.DHwP.D <- T_AvC.DHwC.D <- T_AvC.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AvC.D <- cbind(T_AvC.DEu.D, T_AvC.DEvA.D, T_AvC.DEvP.D, T_AvC.DEvC.D, T_AvC.DEwA.D, T_AvC.DEwP.D, T_AvC.DEwC.D, T_AvC.DEbP.D,
                     T_AvC.DAu.D, T_AvC.DAvA.D, T_AvC.DAvP.D, T_AvC.DAvC.D, T_AvC.DAwA.D, T_AvC.DAwP.D, T_AvC.DAwC.D, T_AvC.DAbP.D,
                     T_AvC.DIu.D, T_AvC.DIvA.D, T_AvC.DIvP.D, T_AvC.DIvC.D, T_AvC.DIwA.D, T_AvC.DIwP.D, T_AvC.DIwC.D, T_AvC.DHbP.D,
                     T_AvC.DHu.D, T_AvC.DHvA.D, T_AvC.DHvP.D, T_AvC.DHvC.D, T_AvC.DHwA.D, T_AvC.DHwP.D, T_AvC.DHwC.D, T_AvC.DIbP.D)
    
    ## AwA
    # AwA.DxE cols
    T_AwA.DEwA.D <- diag((1-ihr_wA.D)*gamma.D*asymp_wA.D, ncol = age.bins, nrow = age.bins)
    T_AwA.DEu.D <- T_AwA.DEvP.D <- 
      T_AwA.DEvC.D <- T_AwA.DEvA.D <- T_AwA.DEwP.D <- T_AwA.DEwC.D <- T_AwA.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwA.DxA cols
    T_AwA.DAwA.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AwA.DAu.D <- T_AwA.DAvP.D <- T_AwA.DAvC.D <- 
      T_AwA.DAvA.D <- T_AwA.DAwP.D <- T_AwA.DAwC.D <- T_AwA.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwA.DxI cols
    T_AwA.DIu.D <- T_AwA.DIvA.D <- T_AwA.DIvP.D <- T_AwA.DIvC.D <- 
      T_AwA.DIwA.D <- T_AwA.DIwP.D <- T_AwA.DIwC.D <- T_AwA.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwA.DxH cols
    T_AwA.DHu.D <- T_AwA.DHvA.D <- T_AwA.DHvP.D <- T_AwA.DHvC.D <- 
      T_AwA.DHwA.D <- T_AwA.DHwP.D <- T_AwA.DHwC.D <- T_AwA.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AwA.D <- cbind(T_AwA.DEu.D, T_AwA.DEvA.D, T_AwA.DEvP.D, T_AwA.DEvC.D, T_AwA.DEwA.D, T_AwA.DEwP.D, T_AwA.DEwC.D, T_AwA.DEbP.D,
                     T_AwA.DAu.D, T_AwA.DAvA.D, T_AwA.DAvP.D, T_AwA.DAvC.D, T_AwA.DAwA.D, T_AwA.DAwP.D, T_AwA.DAwC.D, T_AwA.DAbP.D,
                     T_AwA.DIu.D, T_AwA.DIvA.D, T_AwA.DIvP.D, T_AwA.DIvC.D, T_AwA.DIwA.D, T_AwA.DIwP.D, T_AwA.DIwC.D, T_AwA.DIbP.D,
                     T_AwA.DHu.D, T_AwA.DHvA.D, T_AwA.DHvP.D, T_AwA.DHvC.D, T_AwA.DHwA.D, T_AwA.DHwP.D, T_AwA.DHwC.D, T_AwA.DHbP.D)
    
    ## AwP
    # AwP.DxE cols
    T_AwP.DEwP.D <- diag((1-ihr_wP.D)*gamma.D*asymp_wP.D, ncol = age.bins, nrow = age.bins)
    T_AwP.DEu.D <- T_AwP.DEvA.D <- 
      T_AwP.DEvC.D <- T_AwP.DEwA.D <- T_AwP.DEvP.D <- T_AwP.DEwC.D <- T_AwP.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwP.DxA cols
    T_AwP.DAwP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AwP.DAu.D <- T_AwP.DAvA.D <- T_AwP.DAvC.D <- 
      T_AwP.DAwA.D <- T_AwP.DAvP.D <- T_AwP.DAwC.D <- T_AwP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwP.DxI cols
    T_AwP.DIu.D <- T_AwP.DIvA.D <- T_AwP.DIvP.D <- T_AwP.DIvC.D <- 
      T_AwP.DIwA.D <- T_AwP.DIwP.D <- T_AwP.DIwC.D <- T_AwP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwP.DxH cols
    T_AwP.DHu.D <- T_AwP.DHvA.D <- T_AwP.DHvP.D <- T_AwP.DHvC.D <- 
      T_AwP.DHwA.D <- T_AwP.DHwP.D <- T_AwP.DHwC.D <- T_AwP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AwP.D <- cbind(T_AwP.DEu.D, T_AwP.DEvA.D, T_AwP.DEvP.D, T_AwP.DEvC.D, T_AwP.DEwA.D, T_AwP.DEwP.D, T_AwP.DEwC.D, T_AwP.DEbP.D,
                     T_AwP.DAu.D, T_AwP.DAvA.D, T_AwP.DAvP.D, T_AwP.DAvC.D, T_AwP.DAwA.D, T_AwP.DAwP.D, T_AwP.DAwC.D, T_AwP.DAbP.D,
                     T_AwP.DIu.D, T_AwP.DIvA.D, T_AwP.DIvP.D, T_AwP.DIvC.D, T_AwP.DIwA.D, T_AwP.DIwP.D, T_AwP.DIwC.D, T_AwP.DIbP.D,
                     T_AwP.DHu.D, T_AwP.DHvA.D, T_AwP.DHvP.D, T_AwP.DHvC.D, T_AwP.DHwA.D, T_AwP.DHwP.D, T_AwP.DHwC.D, T_AwP.DHbP.D)
    
    ## AwC
    # AwC.DxE cols
    T_AwC.DEwC.D <- diag((1-ihr_wC.D)*gamma.D*asymp_wC.D, ncol = age.bins, nrow = age.bins)
    T_AwC.DEu.D <- T_AwC.DEvA.D <- 
      T_AwC.DEvP.D <- T_AwC.DEwA.D <- T_AwC.DEwP.D <- T_AwC.DEvC.D <- T_AwC.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwC.DxA cols
    T_AwC.DAwC.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AwC.DAu.D <- T_AwC.DAvA.D <- T_AwC.DAvP.D <- 
      T_AwC.DAwA.D <- T_AwC.DAwP.D <- T_AwC.DAvC.D <- T_AwC.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwC.DxI cols
    T_AwC.DIu.D <- T_AwC.DIvA.D <- T_AwC.DIvP.D <- T_AwC.DIvC.D <- 
      T_AwC.DIwA.D <- T_AwC.DIwP.D <- T_AwC.DIwC.D <- T_AwC.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwC.DxH cols
    T_AwC.DHu.D <- T_AwC.DHvA.D <- T_AwC.DHvP.D <- T_AwC.DHvC.D <- 
      T_AwC.DHwA.D <- T_AwC.DHwP.D <- T_AwC.DHwC.D <- T_AwC.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AwC.D <- cbind(T_AwC.DEu.D, T_AwC.DEvA.D, T_AwC.DEvP.D, T_AwC.DEvC.D, T_AwC.DEwA.D, T_AwC.DEwP.D, T_AwC.DEwC.D, T_AwC.DEbP.D,
                     T_AwC.DAu.D, T_AwC.DAvA.D, T_AwC.DAvP.D, T_AwC.DAvC.D, T_AwC.DAwA.D, T_AwC.DAwP.D, T_AwC.DAwC.D, T_AwC.DAbP.D,
                     T_AwC.DIu.D, T_AwC.DIvA.D, T_AwC.DIvP.D, T_AwC.DIvC.D, T_AwC.DIwA.D, T_AwC.DIwP.D, T_AwC.DIwC.D, T_AwC.DIbP.D,
                     T_AwC.DHu.D, T_AwC.DHvA.D, T_AwC.DHvP.D, T_AwC.DHvC.D, T_AwC.DHwA.D, T_AwC.DHwP.D, T_AwC.DHwC.D, T_AwC.DHbP.D)
    
    ## AbP
    # AbP.DxE cols
    T_AbP.DEbP.D <- diag((1-ihr_bP.D)*gamma.D*asymp_bP.D, ncol = age.bins, nrow = age.bins)
    T_AbP.DEu.D <- T_AbP.DEvA.D <- 
      T_AbP.DEvP.D <- T_AbP.DEwA.D <- T_AbP.DEwP.D <- T_AbP.DEvC.D <- T_AbP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AbP.DxA cols
    T_AbP.DAbP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AbP.DAu.D <- T_AbP.DAvA.D <- T_AbP.DAvP.D <- 
      T_AbP.DAwA.D <- T_AbP.DAwP.D <- T_AbP.DAvC.D <- T_AbP.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AbP.DxI cols
    T_AbP.DIu.D <- T_AbP.DIvA.D <- T_AbP.DIvP.D <- T_AbP.DIvC.D <- 
      T_AbP.DIwA.D <- T_AbP.DIwP.D <- T_AbP.DIwC.D <- T_AbP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AbP.DxH cols
    T_AbP.DHu.D <- T_AbP.DHvA.D <- T_AbP.DHvP.D <- T_AbP.DHvC.D <- 
      T_AbP.DHwA.D <- T_AbP.DHwP.D <- T_AbP.DHwC.D <- T_AbP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AbP.D <- cbind(T_AbP.DEu.D, T_AbP.DEvA.D, T_AbP.DEvP.D, T_AbP.DEvC.D, T_AbP.DEwA.D, T_AbP.DEwP.D, T_AbP.DEwC.D, T_AbP.DEbP.D,
                     T_AbP.DAu.D, T_AbP.DAvA.D, T_AbP.DAvP.D, T_AbP.DAvC.D, T_AbP.DAwA.D, T_AbP.DAwP.D, T_AbP.DAwC.D, T_AbP.DAbP.D,
                     T_AbP.DIu.D, T_AbP.DIvA.D, T_AbP.DIvP.D, T_AbP.DIvC.D, T_AbP.DIwA.D, T_AbP.DIwP.D, T_AbP.DIwC.D, T_AbP.DIbP.D,
                     T_AbP.DHu.D, T_AbP.DHvA.D, T_AbP.DHvP.D, T_AbP.DHvC.D, T_AbP.DHwA.D, T_AbP.DHwP.D, T_AbP.DHwC.D, T_AbP.DHbP.D)
    #all A's together
    T_A.D <- rbind(T_Au.D, T_AvA.D, T_AvP.D, T_AvC.D, T_AwA.D, T_AwP.D, T_AwC.D, T_AbP.D)
    
    ## Iu
    # Iu.DxE cols
    T_Iu.DEu.D <- diag((1-ihr.D)*gamma.D*(1-asymp), ncol = age.bins, nrow = age.bins)
    T_Iu.DEvA.D <- T_Iu.DEvP.D <- 
      T_Iu.DEvC.D <- T_Iu.DEwA.D <- T_Iu.DEwP.D <- T_Iu.DEwC.D <- T_Iu.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Iu.DxI cols
    T_Iu.DIu.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_Iu.DIvA.D <- T_Iu.DIvP.D <- T_Iu.DIvC.D <- 
      T_Iu.DIwA.D <- T_Iu.DIwP.D <- T_Iu.DIwC.D <- T_Iu.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Iu.DxA cols
    T_Iu.DAu.D <- T_Iu.DAvA.D <- T_Iu.DAvP.D <- T_Iu.DAvC.D <- 
      T_Iu.DAwA.D <- T_Iu.DAwP.D <- T_Iu.DAwC.D <- T_Iu.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Iu.DxH cols
    T_Iu.DHu.D <- T_Iu.DHvA.D <- T_Iu.DHvP.D <- T_Iu.DHvC.D <- 
      T_Iu.DHwA.D <- T_Iu.DHwP.D <- T_Iu.DHwC.D <- T_Iu.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_Iu.D <- cbind(T_Iu.DEu.D, T_Iu.DEvA.D, T_Iu.DEvP.D, T_Iu.DEvC.D, T_Iu.DEwA.D, T_Iu.DEwP.D, T_Iu.DEwC.D, T_Iu.DEbP.D,
                    T_Iu.DAu.D, T_Iu.DAvA.D, T_Iu.DAvP.D, T_Iu.DAvC.D, T_Iu.DAwA.D, T_Iu.DAwP.D, T_Iu.DAwC.D, T_Iu.DAbP.D,
                    T_Iu.DIu.D, T_Iu.DIvA.D, T_Iu.DIvP.D, T_Iu.DIvC.D, T_Iu.DIwA.D, T_Iu.DIwP.D, T_Iu.DIwC.D, T_Iu.DIbP.D,
                    T_Iu.DHu.D, T_Iu.DHvA.D, T_Iu.DHvP.D, T_Iu.DHvC.D, T_Iu.DHwA.D, T_Iu.DHwP.D, T_Iu.DHwC.D, T_Iu.DHbP.D)
    
    ## IvA
    # IvA.DxE cols
    T_IvA.DEvA.D <- diag((1-ihr_vA.D)*gamma.D*(1-asymp_vA.D), ncol = age.bins, nrow = age.bins)
    T_IvA.DEu.D <- T_IvA.DEvP.D <- 
      T_IvA.DEvC.D <- T_IvA.DEwA.D <- T_IvA.DEwP.D <- T_IvA.DEwC.D <- T_IvA.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvA.DxI cols
    T_IvA.DIvA.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IvA.DIu.D <- T_IvA.DIvP.D <- T_IvA.DIvC.D <- 
      T_IvA.DIwA.D <- T_IvA.DIwP.D <- T_IvA.DIwC.D <- T_IvA.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvA.DxA cols
    T_IvA.DAu.D <- T_IvA.DAvA.D <- T_IvA.DAvP.D <- T_IvA.DAvC.D <- 
      T_IvA.DAwA.D <- T_IvA.DAwP.D <- T_IvA.DAwC.D <- T_IvA.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvA.DxH cols
    T_IvA.DHu.D <- T_IvA.DHvA.D <- T_IvA.DHvP.D <- T_IvA.DHvC.D <- 
      T_IvA.DHwA.D <- T_IvA.DHwP.D <- T_IvA.DHwC.D <- T_IvA.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IvA.D <- cbind(T_IvA.DEu.D, T_IvA.DEvA.D, T_IvA.DEvP.D, T_IvA.DEvC.D, T_IvA.DEwA.D, T_IvA.DEwP.D, T_IvA.DEwC.D, T_IvA.DEbP.D,
                     T_IvA.DAu.D, T_IvA.DAvA.D, T_IvA.DAvP.D, T_IvA.DAvC.D, T_IvA.DAwA.D, T_IvA.DAwP.D, T_IvA.DAwC.D, T_IvA.DAbP.D,
                     T_IvA.DIu.D, T_IvA.DIvA.D, T_IvA.DIvP.D, T_IvA.DIvC.D, T_IvA.DIwA.D, T_IvA.DIwP.D, T_IvA.DIwC.D, T_IvA.DIbP.D,
                     T_IvA.DHu.D, T_IvA.DHvA.D, T_IvA.DHvP.D, T_IvA.DHvC.D, T_IvA.DHwA.D, T_IvA.DHwP.D, T_IvA.DHwC.D, T_IvA.DHbP.D)
    
    ## IvP
    # IvP.DxE cols
    T_IvP.DEvP.D <- diag((1-ihr_vP.D)*gamma.D*(1-asymp_vP.D), ncol = age.bins, nrow = age.bins)
    T_IvP.DEu.D <- T_IvP.DEvA.D <- 
      T_IvP.DEvC.D <- T_IvP.DEwA.D <- T_IvP.DEwP.D <- T_IvP.DEwC.D <- T_IvP.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvP.DxI cols
    T_IvP.DIvP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IvP.DIu.D <- T_IvP.DIvA.D <- T_IvP.DIvC.D <- 
      T_IvP.DIwA.D <- T_IvP.DIwP.D <- T_IvP.DIwC.D <- T_IvP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvP.DxA cols
    T_IvP.DAu.D <- T_IvP.DAvA.D <- T_IvP.DAvP.D <- T_IvP.DAvC.D <- 
      T_IvP.DAwA.D <- T_IvP.DAwP.D <- T_IvP.DAwC.D <- T_IvP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvP.DxH cols
    T_IvP.DHu.D <- T_IvP.DHvA.D <- T_IvP.DHvP.D <- T_IvP.DHvC.D <- 
      T_IvP.DHwA.D <- T_IvP.DHwP.D <- T_IvP.DHwC.D <- T_IvP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IvP.D <- cbind(T_IvP.DEu.D, T_IvP.DEvA.D, T_IvP.DEvP.D, T_IvP.DEvC.D, T_IvP.DEwA.D, T_IvP.DEwP.D, T_IvP.DEwC.D, T_IvP.DEbP.D,
                     T_IvP.DAu.D, T_IvP.DAvA.D, T_IvP.DAvP.D, T_IvP.DAvC.D, T_IvP.DAwA.D, T_IvP.DAwP.D, T_IvP.DAwC.D, T_IvP.DAbP.D,
                     T_IvP.DIu.D, T_IvP.DIvA.D, T_IvP.DIvP.D, T_IvP.DIvC.D, T_IvP.DIwA.D, T_IvP.DIwP.D, T_IvP.DIwC.D, T_IvP.DIbP.D,
                     T_IvP.DHu.D, T_IvP.DHvA.D, T_IvP.DHvP.D, T_IvP.DHvC.D, T_IvP.DHwA.D, T_IvP.DHwP.D, T_IvP.DHwC.D, T_IvP.DHbP.D)
    
    ## IvC
    # IvC.DxE cols
    T_IvC.DEvC.D <- diag((1-ihr_vC.D)*gamma.D*(1-asymp_vC.D), ncol = age.bins, nrow = age.bins)
    T_IvC.DEu.D <- T_IvC.DEvA.D <- 
      T_IvC.DEvP.D <- T_IvC.DEwA.D <- T_IvC.DEwP.D <- T_IvC.DEwC.D <- T_IvC.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvC.DxI cols
    T_IvC.DIvC.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IvC.DIu.D <- T_IvC.DIvA.D <- T_IvC.DIvP.D <- 
      T_IvC.DIwA.D <- T_IvC.DIwP.D <- T_IvC.DIwC.D <- T_IvC.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvC.DxA cols
    T_IvC.DAu.D <- T_IvC.DAvA.D <- T_IvC.DAvP.D <- T_IvC.DAvC.D <- 
      T_IvC.DAwA.D <- T_IvC.DAwP.D <- T_IvC.DAwC.D <- T_IvC.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvC.DxH cols
    T_IvC.DHu.D <- T_IvC.DHvA.D <- T_IvC.DHvP.D <- T_IvC.DHvC.D <- 
      T_IvC.DHwA.D <- T_IvC.DHwP.D <- T_IvC.DHwC.D <- T_IvC.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IvC.D <- cbind(T_IvC.DEu.D, T_IvC.DEvA.D, T_IvC.DEvP.D, T_IvC.DEvC.D, T_IvC.DEwA.D, T_IvC.DEwP.D, T_IvC.DEwC.D, T_IvC.DEbP.D,
                     T_IvC.DAu.D, T_IvC.DAvA.D, T_IvC.DAvP.D, T_IvC.DAvC.D, T_IvC.DAwA.D, T_IvC.DAwP.D, T_IvC.DAwC.D, T_IvC.DAbP.D,
                     T_IvC.DIu.D, T_IvC.DIvA.D, T_IvC.DIvP.D, T_IvC.DIvC.D, T_IvC.DIwA.D, T_IvC.DIwP.D, T_IvC.DIwC.D, T_IvC.DIbP.D,
                     T_IvC.DHu.D, T_IvC.DHvA.D, T_IvC.DHvP.D, T_IvC.DHvC.D, T_IvC.DHwA.D, T_IvC.DHwP.D, T_IvC.DHwC.D, T_IvC.DHbP.D)
    
    ## IwA
    # IwA.DxE cols
    T_IwA.DEwA.D <- diag((1-ihr_wA.D)*gamma.D*(1-asymp_wA.D), ncol = age.bins, nrow = age.bins)
    T_IwA.DEu.D <- T_IwA.DEvP.D <- 
      T_IwA.DEvC.D <- T_IwA.DEvA.D <- T_IwA.DEwP.D <- T_IwA.DEwC.D <- T_IwA.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwA.DxI cols
    T_IwA.DIwA.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IwA.DIu.D <- T_IwA.DIvP.D <- T_IwA.DIvC.D <- 
      T_IwA.DIvA.D <- T_IwA.DIwP.D <- T_IwA.DIwC.D <- T_IwA.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwA.DxA cols
    T_IwA.DAu.D <- T_IwA.DAvA.D <- T_IwA.DAvP.D <- T_IwA.DAvC.D <- 
      T_IwA.DAwA.D <- T_IwA.DAwP.D <- T_IwA.DAwC.D <- T_IwA.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwA.DxH cols
    T_IwA.DHu.D <- T_IwA.DHvA.D <- T_IwA.DHvP.D <- T_IwA.DHvC.D <- 
      T_IwA.DHwA.D <- T_IwA.DHwP.D <- T_IwA.DHwC.D <- T_IwA.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IwA.D <- cbind(T_IwA.DEu.D, T_IwA.DEvA.D, T_IwA.DEvP.D, T_IwA.DEvC.D, T_IwA.DEwA.D, T_IwA.DEwP.D, T_IwA.DEwC.D, T_IwA.DEbP.D,
                     T_IwA.DAu.D, T_IwA.DAvA.D, T_IwA.DAvP.D, T_IwA.DAvC.D, T_IwA.DAwA.D, T_IwA.DAwP.D, T_IwA.DAwC.D, T_IwA.DAbP.D,
                     T_IwA.DIu.D, T_IwA.DIvA.D, T_IwA.DIvP.D, T_IwA.DIvC.D, T_IwA.DIwA.D, T_IwA.DIwP.D, T_IwA.DIwC.D, T_IwA.DIbP.D,
                     T_IwA.DHu.D, T_IwA.DHvA.D, T_IwA.DHvP.D, T_IwA.DHvC.D, T_IwA.DHwA.D, T_IwA.DHwP.D, T_IwA.DHwC.D, T_IwA.DHbP.D)
    
    ## IwP
    # IwP.DxE cols
    T_IwP.DEwP.D <- diag((1-ihr_wP.D)*gamma.D*(1-asymp_wP.D), ncol = age.bins, nrow = age.bins)
    T_IwP.DEu.D <- T_IwP.DEvA.D <- 
      T_IwP.DEvC.D <- T_IwP.DEwA.D <- T_IwP.DEvP.D <- T_IwP.DEwC.D <- T_IwP.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwP.DxI cols
    T_IwP.DIwP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IwP.DIu.D <- T_IwP.DIvA.D <- T_IwP.DIvC.D <- 
      T_IwP.DIwA.D <- T_IwP.DIvP.D <- T_IwP.DIwC.D <- T_IwP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwP.DxA cols
    T_IwP.DAu.D <- T_IwP.DAvA.D <- T_IwP.DAvP.D <- T_IwP.DAvC.D <- 
      T_IwP.DAwA.D <- T_IwP.DAwP.D <- T_IwP.DAwC.D <- T_IwP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwP.DxH cols
    T_IwP.DHu.D <- T_IwP.DHvA.D <- T_IwP.DHvP.D <- T_IwP.DHvC.D <- 
      T_IwP.DHwA.D <- T_IwP.DHwP.D <- T_IwP.DHwC.D <- T_IwP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IwP.D <- cbind(T_IwP.DEu.D, T_IwP.DEvA.D, T_IwP.DEvP.D, T_IwP.DEvC.D, T_IwP.DEwA.D, T_IwP.DEwP.D, T_IwP.DEwC.D, T_IwP.DEbP.D,
                     T_IwP.DAu.D, T_IwP.DAvA.D, T_IwP.DAvP.D, T_IwP.DAvC.D, T_IwP.DAwA.D, T_IwP.DAwP.D, T_IwP.DAwC.D, T_IwP.DAbP.D,
                     T_IwP.DIu.D, T_IwP.DIvA.D, T_IwP.DIvP.D, T_IwP.DIvC.D, T_IwP.DIwA.D, T_IwP.DIwP.D, T_IwP.DIwC.D, T_IwP.DIbP.D,
                     T_IwP.DHu.D, T_IwP.DHvA.D, T_IwP.DHvP.D, T_IwP.DHvC.D, T_IwP.DHwA.D, T_IwP.DHwP.D, T_IwP.DHwC.D, T_IwP.DHbP.D)
    
    ## IwC
    # IwC.DxE cols
    T_IwC.DEwC.D <- diag((1-ihr_wC.D)*gamma.D*(1-asymp_wC.D), ncol = age.bins, nrow = age.bins)
    T_IwC.DEu.D <- T_IwC.DEvA.D <- 
      T_IwC.DEvP.D <- T_IwC.DEwA.D <- T_IwC.DEwP.D <- T_IwC.DEvC.D <- T_IwC.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.DxI cols
    T_IwC.DIwC.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IwC.DIu.D <- T_IwC.DIvA.D <- T_IwC.DIvP.D <- 
      T_IwC.DIwA.D <- T_IwC.DIwP.D <- T_IwC.DIvC.D <- T_IwC.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.DxA cols
    T_IwC.DAu.D <- T_IwC.DAvA.D <- T_IwC.DAvP.D <- T_IwC.DAvC.D <- 
      T_IwC.DAwA.D <- T_IwC.DAwP.D <- T_IwC.DAwC.D <- T_IwC.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.DxH cols
    T_IwC.DHu.D <- T_IwC.DHvA.D <- T_IwC.DHvP.D <- T_IwC.DHvC.D <- 
      T_IwC.DHwA.D <- T_IwC.DHwP.D <- T_IwC.DHwC.D <- T_IwC.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IwC.D <- cbind(T_IwC.DEu.D, T_IwC.DEvA.D, T_IwC.DEvP.D, T_IwC.DEvC.D, T_IwC.DEwA.D, T_IwC.DEwP.D, T_IwC.DEwC.D, T_IwC.DEbP.D,
                     T_IwC.DAu.D, T_IwC.DAvA.D, T_IwC.DAvP.D, T_IwC.DAvC.D, T_IwC.DAwA.D, T_IwC.DAwP.D, T_IwC.DAwC.D, T_IwC.DAbP.D,
                     T_IwC.DIu.D, T_IwC.DIvA.D, T_IwC.DIvP.D, T_IwC.DIvC.D, T_IwC.DIwA.D, T_IwC.DIwP.D, T_IwC.DIwC.D, T_IwC.DIbP.D,
                     T_IwC.DHu.D, T_IwC.DHvA.D, T_IwC.DHvP.D, T_IwC.DHvC.D, T_IwC.DHwA.D, T_IwC.DHwP.D, T_IwC.DHwC.D, T_IwC.DHbP.D)
 
    ## IwC
    # IwC.DxE cols
    T_IbP.DEbP.D <- diag((1-ihr_bP.D)*gamma.D*(1-asymp_bP.D), ncol = age.bins, nrow = age.bins)
    T_IbP.DEu.D <- T_IbP.DEvA.D <- 
      T_IbP.DEvP.D <- T_IbP.DEwA.D <- T_IbP.DEwP.D <- T_IbP.DEvC.D <- T_IbP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.DxI cols
    T_IbP.DIbP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IbP.DIu.D <- T_IbP.DIvA.D <- T_IbP.DIvP.D <- 
      T_IbP.DIwA.D <- T_IbP.DIwP.D <- T_IbP.DIvC.D <- T_IbP.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.DxA cols
    T_IbP.DAu.D <- T_IbP.DAvA.D <- T_IbP.DAvP.D <- T_IbP.DAvC.D <- 
      T_IbP.DAwA.D <- T_IbP.DAwP.D <- T_IbP.DAwC.D <- T_IbP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.DxH cols
    T_IbP.DHu.D <- T_IbP.DHvA.D <- T_IbP.DHvP.D <- T_IbP.DHvC.D <- 
      T_IbP.DHwA.D <- T_IbP.DHwP.D <- T_IbP.DHwC.D <- T_IbP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IbP.D <- cbind(T_IbP.DEu.D, T_IbP.DEvA.D, T_IbP.DEvP.D, T_IbP.DEvC.D, T_IbP.DEwA.D, T_IbP.DEwP.D, T_IbP.DEwC.D, T_IbP.DEbP.D,
                     T_IbP.DAu.D, T_IbP.DAvA.D, T_IbP.DAvP.D, T_IbP.DAvC.D, T_IbP.DAwA.D, T_IbP.DAwP.D, T_IbP.DAwC.D, T_IbP.DAbP.D,
                     T_IbP.DIu.D, T_IbP.DIvA.D, T_IbP.DIvP.D, T_IbP.DIvC.D, T_IbP.DIwA.D, T_IbP.DIwP.D, T_IbP.DIwC.D, T_IbP.DIbP.D,
                     T_IbP.DHu.D, T_IbP.DHvA.D, T_IbP.DHvP.D, T_IbP.DHvC.D, T_IbP.DHwA.D, T_IbP.DHwP.D, T_IbP.DHwC.D, T_IbP.DHbP.D)   
    #all I's together
    T_I.D <- rbind(T_Iu.D, T_IvA.D, T_IvP.D, T_IvC.D, T_IwA.D, T_IwP.D, T_IwC.D, T_IbP.D)
    
    ## Hu
    # Hu.DxE cols
    T_Hu.DEu.D <- diag(ihr.D*gamma.D, ncol = age.bins, nrow = age.bins)
    T_Hu.DEvA.D <- T_Hu.DEvP.D <- 
      T_Hu.DEvC.D <- T_Hu.DEwA.D <- T_Hu.DEwP.D <- T_Hu.DEwC.D <- T_Hu.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Hu.DxH cols
    T_Hu.DHu.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_Hu.DHvA.D <- T_Hu.DHvP.D <- T_Hu.DHvC.D <- 
      T_Hu.DHwA.D <- T_Hu.DHwP.D <- T_Hu.DHwC.D <- T_Hu.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Hu.DxI cols
    T_Hu.DIu.D <- T_Hu.DIvA.D <- T_Hu.DIvP.D <- T_Hu.DIvC.D <- 
      T_Hu.DIwA.D <- T_Hu.DIwP.D <- T_Hu.DIwC.D <- T_Hu.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Hu.DxA cols
    T_Hu.DAu.D <- T_Hu.DAvA.D <- T_Hu.DAvP.D <- T_Hu.DAvC.D <- 
      T_Hu.DAwA.D <- T_Hu.DAwP.D <- T_Hu.DAwC.D <- T_Hu.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_Hu.D <- cbind(T_Hu.DEu.D, T_Hu.DEvA.D, T_Hu.DEvP.D, T_Hu.DEvC.D, T_Hu.DEwA.D, T_Hu.DEwP.D, T_Hu.DEwC.D, T_Hu.DEbP.D,
                    T_Hu.DAu.D, T_Hu.DAvA.D, T_Hu.DAvP.D, T_Hu.DAvC.D, T_Hu.DAwA.D, T_Hu.DAwP.D, T_Hu.DAwC.D, T_Hu.DAbP.D,
                    T_Hu.DIu.D, T_Hu.DIvA.D, T_Hu.DIvP.D, T_Hu.DIvC.D, T_Hu.DIwA.D, T_Hu.DIwP.D, T_Hu.DIwC.D, T_Hu.DIbP.D,
                    T_Hu.DHu.D, T_Hu.DHvA.D, T_Hu.DHvP.D, T_Hu.DHvC.D, T_Hu.DHwA.D, T_Hu.DHwP.D, T_Hu.DHwC.D, T_Hu.DHbP.D)
    
    ## HvA
    # HvA.DxE cols
    T_HvA.DEvA.D <- diag(ihr_vA.D*gamma.D, ncol = age.bins, nrow = age.bins)
    T_HvA.DEu.D <- T_HvA.DEvP.D <- 
      T_HvA.DEvC.D <- T_HvA.DEwA.D <- T_HvA.DEwP.D <- T_HvA.DEwC.D <- T_HvA.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvA.DxH cols
    T_HvA.DHvA.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HvA.DHu.D <- T_HvA.DHvP.D <- T_HvA.DHvC.D <- 
      T_HvA.DHwA.D <- T_HvA.DHwP.D <- T_HvA.DHwC.D <- T_HvA.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvA.DxI cols
    T_HvA.DIu.D <- T_HvA.DIvA.D <- T_HvA.DIvP.D <- T_HvA.DIvC.D <- 
      T_HvA.DIwA.D <- T_HvA.DIwP.D <- T_HvA.DIwC.D <- T_HvA.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvA.DxA cols
    T_HvA.DAu.D <- T_HvA.DAvA.D <- T_HvA.DAvP.D <- T_HvA.DAvC.D <- 
      T_HvA.DAwA.D <- T_HvA.DAwP.D <- T_HvA.DAwC.D <- T_HvA.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HvA.D <- cbind(T_HvA.DEu.D, T_HvA.DEvA.D, T_HvA.DEvP.D, T_HvA.DEvC.D, T_HvA.DEwA.D, T_HvA.DEwP.D, T_HvA.DEwC.D, T_HvA.DEbP.D,
                     T_HvA.DAu.D, T_HvA.DAvA.D, T_HvA.DAvP.D, T_HvA.DAvC.D, T_HvA.DAwA.D, T_HvA.DAwP.D, T_HvA.DAwC.D, T_HvA.DAbP.D,
                     T_HvA.DIu.D, T_HvA.DIvA.D, T_HvA.DIvP.D, T_HvA.DIvC.D, T_HvA.DIwA.D, T_HvA.DIwP.D, T_HvA.DIwC.D, T_HvA.DIbP.D,
                     T_HvA.DHu.D, T_HvA.DHvA.D, T_HvA.DHvP.D, T_HvA.DHvC.D, T_HvA.DHwA.D, T_HvA.DHwP.D, T_HvA.DHwC.D, T_HvA.DHbP.D)
    
    ## HvP
    # HvP.DxE cols
    T_HvP.DEvP.D <- diag(ihr_vP.D*gamma.D, ncol = age.bins, nrow = age.bins)
    T_HvP.DEu.D <- T_HvP.DEvA.D <- 
      T_HvP.DEvC.D <- T_HvP.DEwA.D <- T_HvP.DEwP.D <- T_HvP.DEwC.D <- T_HvP.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvP.DxH cols
    T_HvP.DHvP.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HvP.DHu.D <- T_HvP.DHvA.D <- T_HvP.DHvC.D <- 
      T_HvP.DHwA.D <- T_HvP.DHwP.D <- T_HvP.DHwC.D <- T_HvP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvP.DxI cols
    T_HvP.DIu.D <- T_HvP.DIvA.D <- T_HvP.DIvP.D <- T_HvP.DIvC.D <- 
      T_HvP.DIwA.D <- T_HvP.DIwP.D <- T_HvP.DIwC.D <- T_HvP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvP.DxA cols
    T_HvP.DAu.D <- T_HvP.DAvA.D <- T_HvP.DAvP.D <- T_HvP.DAvC.D <- 
      T_HvP.DAwA.D <- T_HvP.DAwP.D <- T_HvP.DAwC.D <- T_HvP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HvP.D <- cbind(T_HvP.DEu.D, T_HvP.DEvA.D, T_HvP.DEvP.D, T_HvP.DEvC.D, T_HvP.DEwA.D, T_HvP.DEwP.D, T_HvP.DEwC.D, T_HvP.DEbP.D,
                     T_HvP.DAu.D, T_HvP.DAvA.D, T_HvP.DAvP.D, T_HvP.DAvC.D, T_HvP.DAwA.D, T_HvP.DAwP.D, T_HvP.DAwC.D, T_HvP.DAbP.D,
                     T_HvP.DIu.D, T_HvP.DIvA.D, T_HvP.DIvP.D, T_HvP.DIvC.D, T_HvP.DIwA.D, T_HvP.DIwP.D, T_HvP.DIwC.D, T_HvP.DIbP.D,
                     T_HvP.DHu.D, T_HvP.DHvA.D, T_HvP.DHvP.D, T_HvP.DHvC.D, T_HvP.DHwA.D, T_HvP.DHwP.D, T_HvP.DHwC.D, T_HvP.DHbP.D)
    
    ## HvC
    # HvC.DxE cols
    T_HvC.DEvC.D <- diag(ihr_vC.D*gamma.D, ncol = age.bins, nrow = age.bins)
    T_HvC.DEu.D <- T_HvC.DEvA.D <- 
      T_HvC.DEvP.D <- T_HvC.DEwA.D <- T_HvC.DEwP.D <- T_HvC.DEwC.D <- T_HvC.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvC.DxH cols
    T_HvC.DHvC.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HvC.DHu.D <- T_HvC.DHvA.D <- T_HvC.DHvP.D <- 
      T_HvC.DHwA.D <- T_HvC.DHwP.D <- T_HvC.DHwC.D <- T_HvC.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvC.DxI cols
    T_HvC.DIu.D <- T_HvC.DIvA.D <- T_HvC.DIvP.D <- T_HvC.DIvC.D <- 
      T_HvC.DIwA.D <- T_HvC.DIwP.D <- T_HvC.DIwC.D <- T_HvC.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvC.DxA cols
    T_HvC.DAu.D <- T_HvC.DAvA.D <- T_HvC.DAvP.D <- T_HvC.DAvC.D <- 
      T_HvC.DAwA.D <- T_HvC.DAwP.D <- T_HvC.DAwC.D <- T_HvC.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HvC.D <- cbind(T_HvC.DEu.D, T_HvC.DEvA.D, T_HvC.DEvP.D, T_HvC.DEvC.D, T_HvC.DEwA.D, T_HvC.DEwP.D, T_HvC.DEwC.D, T_HvC.DEbP.D,
                     T_HvC.DAu.D, T_HvC.DAvA.D, T_HvC.DAvP.D, T_HvC.DAvC.D, T_HvC.DAwA.D, T_HvC.DAwP.D, T_HvC.DAwC.D, T_HvC.DAbP.D,
                     T_HvC.DIu.D, T_HvC.DIvA.D, T_HvC.DIvP.D, T_HvC.DIvC.D, T_HvC.DIwA.D, T_HvC.DIwP.D, T_HvC.DIwC.D, T_HvC.DIbP.D,
                     T_HvC.DHu.D, T_HvC.DHvA.D, T_HvC.DHvP.D, T_HvC.DHvC.D, T_HvC.DHwA.D, T_HvC.DHwP.D, T_HvC.DHwC.D, T_HvC.DHbP.D)
    
    ## HwA
    # HwA.DxE cols
    T_HwA.DEwA.D <- diag(ihr_wA.D*gamma.D, ncol = age.bins, nrow = age.bins)
    T_HwA.DEu.D <- T_HwA.DEvP.D <- 
      T_HwA.DEvC.D <- T_HwA.DEvA.D <- T_HwA.DEwP.D <- T_HwA.DEwC.D <- T_HwA.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwA.DxH cols
    T_HwA.DHwA.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HwA.DHu.D <- T_HwA.DHvP.D <- T_HwA.DHvC.D <- 
      T_HwA.DHvA.D <- T_HwA.DHwP.D <- T_HwA.DHwC.D <- T_HwA.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwA.DxI cols
    T_HwA.DIu.D <- T_HwA.DIvA.D <- T_HwA.DIvP.D <- T_HwA.DIvC.D <- 
      T_HwA.DIwA.D <- T_HwA.DIwP.D <- T_HwA.DIwC.D <- T_HwA.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwA.DxA cols
    T_HwA.DAu.D <- T_HwA.DAvA.D <- T_HwA.DAvP.D <- T_HwA.DAvC.D <- 
      T_HwA.DAwA.D <- T_HwA.DAwP.D <- T_HwA.DAwC.D <- T_HwA.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HwA.D <- cbind(T_HwA.DEu.D, T_HwA.DEvA.D, T_HwA.DEvP.D, T_HwA.DEvC.D, T_HwA.DEwA.D, T_HwA.DEwP.D, T_HwA.DEwC.D, T_HwA.DEbP.D,
                     T_HwA.DAu.D, T_HwA.DAvA.D, T_HwA.DAvP.D, T_HwA.DAvC.D, T_HwA.DAwA.D, T_HwA.DAwP.D, T_HwA.DAwC.D, T_HwA.DAbP.D,
                     T_HwA.DIu.D, T_HwA.DIvA.D, T_HwA.DIvP.D, T_HwA.DIvC.D, T_HwA.DIwA.D, T_HwA.DIwP.D, T_HwA.DIwC.D, T_HwA.DIbP.D,
                     T_HwA.DHu.D, T_HwA.DHvA.D, T_HwA.DHvP.D, T_HwA.DHvC.D, T_HwA.DHwA.D, T_HwA.DHwP.D, T_HwA.DHwC.D, T_HwA.DHbP.D)
    
    ## HwP
    # HwP.DxE cols
    T_HwP.DEwP.D <- diag(ihr_wP.D*gamma.D, ncol = age.bins, nrow = age.bins)
    T_HwP.DEu.D <- T_HwP.DEvA.D <- 
      T_HwP.DEvC.D <- T_HwP.DEwA.D <- T_HwP.DEvP.D <- T_HwP.DEwC.D <- T_HwP.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwP.DxH cols
    T_HwP.DHwP.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HwP.DHu.D <- T_HwP.DHvA.D <- T_HwP.DHvC.D <- 
      T_HwP.DHwA.D <- T_HwP.DHvP.D <- T_HwP.DHwC.D <- T_HwP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwP.DxI cols
    T_HwP.DIu.D <- T_HwP.DIvA.D <- T_HwP.DIvP.D <- T_HwP.DIvC.D <- 
      T_HwP.DIwA.D <- T_HwP.DIwP.D <- T_HwP.DIwC.D <- T_HwP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwP.DxA cols
    T_HwP.DAu.D <- T_HwP.DAvA.D <- T_HwP.DAvP.D <- T_HwP.DAvC.D <- 
      T_HwP.DAwA.D <- T_HwP.DAwP.D <- T_HwP.DAwC.D <- T_HwP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HwP.D <- cbind(T_HwP.DEu.D, T_HwP.DEvA.D, T_HwP.DEvP.D, T_HwP.DEvC.D, T_HwP.DEwA.D, T_HwP.DEwP.D, T_HwP.DEwC.D, T_HwP.DEbP.D,
                     T_HwP.DAu.D, T_HwP.DAvA.D, T_HwP.DAvP.D, T_HwP.DAvC.D, T_HwP.DAwA.D, T_HwP.DAwP.D, T_HwP.DAwC.D, T_HwP.DAbP.D,
                     T_HwP.DIu.D, T_HwP.DIvA.D, T_HwP.DIvP.D, T_HwP.DIvC.D, T_HwP.DIwA.D, T_HwP.DIwP.D, T_HwP.DIwC.D, T_HwP.DIbP.D,
                     T_HwP.DHu.D, T_HwP.DHvA.D, T_HwP.DHvP.D, T_HwP.DHvC.D, T_HwP.DHwA.D, T_HwP.DHwP.D, T_HwP.DHwC.D, T_HwP.DHbP.D)
    
    ## HwC
    # HwC.DxE cols
    T_HwC.DEwC.D <- diag(ihr_wC.D*gamma.D, ncol = age.bins, nrow = age.bins)
    T_HwC.DEu.D <- T_HwC.DEvA.D <- 
      T_HwC.DEvP.D <- T_HwC.DEwA.D <- T_HwC.DEwP.D <- T_HwC.DEvC.D <- T_HwC.DEbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.DxH cols
    T_HwC.DHwC.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HwC.DHu.D <- T_HwC.DHvA.D <- T_HwC.DHvP.D <- 
      T_HwC.DHwA.D <- T_HwC.DHwP.D <- T_HwC.DHvC.D <- T_HwC.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.DxI cols
    T_HwC.DIu.D <- T_HwC.DIvA.D <- T_HwC.DIvP.D <- T_HwC.DIvC.D <- 
      T_HwC.DIwA.D <- T_HwC.DIwP.D <- T_HwC.DIwC.D <- T_HwC.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.DxA cols
    T_HwC.DAu.D <- T_HwC.DAvA.D <- T_HwC.DAvP.D <- T_HwC.DAvC.D <- 
      T_HwC.DAwA.D <- T_HwC.DAwP.D <- T_HwC.DAwC.D <- T_HwC.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HwC.D <- cbind(T_HwC.DEu.D, T_HwC.DEvA.D, T_HwC.DEvP.D, T_HwC.DEvC.D, T_HwC.DEwA.D, T_HwC.DEwP.D, T_HwC.DEwC.D, T_HwC.DEbP.D,
                     T_HwC.DAu.D, T_HwC.DAvA.D, T_HwC.DAvP.D, T_HwC.DAvC.D, T_HwC.DAwA.D, T_HwC.DAwP.D, T_HwC.DAwC.D, T_HwC.DAbP.D,
                     T_HwC.DIu.D, T_HwC.DIvA.D, T_HwC.DIvP.D, T_HwC.DIvC.D, T_HwC.DIwA.D, T_HwC.DIwP.D, T_HwC.DIwC.D, T_HwC.DIbP.D,
                     T_HwC.DHu.D, T_HwC.DHvA.D, T_HwC.DHvP.D, T_HwC.DHvC.D, T_HwC.DHwA.D, T_HwC.DHwP.D, T_HwC.DHwC.D, T_HwC.DHbP.D)
    
    ## HwC
    # HwC.DxE cols
    T_HbP.DEbP.D <- diag(ihr_bP.D*gamma.D, ncol = age.bins, nrow = age.bins)
    T_HbP.DEu.D <- T_HbP.DEvA.D <- 
      T_HbP.DEvP.D <- T_HbP.DEwA.D <- T_HbP.DEwP.D <- T_HbP.DEvC.D <- T_HbP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.DxH cols
    T_HbP.DHwC.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HbP.DHu.D <- T_HbP.DHvA.D <- T_HbP.DHvP.D <- 
      T_HbP.DHwA.D <- T_HbP.DHwP.D <- T_HbP.DHvC.D <- T_HbP.DHbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.DxI cols
    T_HbP.DIu.D <- T_HbP.DIvA.D <- T_HbP.DIvP.D <- T_HbP.DIvC.D <- 
      T_HbP.DIwA.D <- T_HbP.DIwP.D <- T_HbP.DIwC.D <- T_HbP.DIbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.DxA cols
    T_HbP.DAu.D <- T_HbP.DAvA.D <- T_HbP.DAvP.D <- T_HbP.DAvC.D <- 
      T_HbP.DAwA.D <- T_HbP.DAwP.D <- T_HbP.DAwC.D <- T_HbP.DAbP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HbP.D <- cbind(T_HbP.DEu.D, T_HbP.DEvA.D, T_HbP.DEvP.D, T_HbP.DEvC.D, T_HbP.DEwA.D, T_HbP.DEwP.D, T_HbP.DEwC.D, T_HbP.DEbP.D,
                     T_HbP.DAu.D, T_HbP.DAvA.D, T_HbP.DAvP.D, T_HbP.DAvC.D, T_HbP.DAwA.D, T_HbP.DAwP.D, T_HbP.DAwC.D, T_HbP.DAbP.D,
                     T_HbP.DIu.D, T_HbP.DIvA.D, T_HbP.DIvP.D, T_HbP.DIvC.D, T_HbP.DIwA.D, T_HbP.DIwP.D, T_HbP.DIwC.D, T_HbP.DIbP.D,
                     T_HbP.DHu.D, T_HbP.DHvA.D, T_HbP.DHvP.D, T_HbP.DHvC.D, T_HbP.DHwA.D, T_HbP.DHwP.D, T_HbP.DHwC.D, T_HbP.DHbP.D)
    
    #all H's together
    T_H.D <- rbind(T_Hu.D, T_HvA.D, T_HvP.D, T_HvC.D, T_HwA.D, T_HwP.D, T_HwC.D, T_HbP.D)
    
    
    # final matrix
    T_mat.D <- rbind(T_E.D, T_A.D, T_I.D, T_H.D)
    # identity matrix
    Id_mat <- diag(1,nrow(T_mat.D),ncol(T_mat.D))
    # sparse because why not
    F_mat.D <- Matrix(F_mat.D, sparse = T)
    # automatically computes (I - T)^-1 and converts to sparse -- faster
    T_mat.D <- Matrix(inv(Id_mat - T_mat.D), sparse = T)
    # R0 is the maximum eigenvalue (that might be complex)
    R0.D <- max(Re(eigs(F_mat.D%*%T_mat.D, 1, retvec = F)$values))
    ############################################################################
    ################################ OMICRON ###################################
    ############################################################################
    # X_0 = (Eu.O, EvA.O, EvP.O, EvC.O, EwA.O, EwP.O, EwC.O, EbP.O
    #        Au.O, AvA.O, AvP.O, AvC.O, AwA.O, AwP.O, AwC.O, AbP.O
    #        Iu.O, IvA.O, IvP.O, IvC.O, IwA.O, IwP.O, IwC.O, IbP.O
    #        Hu.O, HvA.O, HvP.O, HvC.O, HwA.O, HwP.O, HwC.O, HbP.O)
    # X_1 = (Su, SvA, SvP, SvP, SvC, SwA, SwP, SwC, SbP
    #        Ru.O, RvA.O, RvP.O, RvC.O, RwA.O, RwP.O, RwC.O, RbP.O
    #        Du.O, DvA.O, DvP.O, DvC.O, DwA.O, DwP.O, DwC.O, DbP.O)
    
    ############################### F matrix ###################################
    # EuDxE cols
    F_Eu.OEu.O <- F_Eu.OEvA.O <- F_Eu.OEvP.O <- F_Eu.OEvC.O <- 
      F_Eu.OEwA.O <- F_Eu.OEwP.O <- F_Eu.OEwC.O <- F_Eu.OEbP.O <- beta.O*omega*(contact.Su + 
                                                                                  reinfec.DO*contact.Ru.D+
                                                                                  reinfec.OO*contact.Ru.O)
    # EuDxA cols
    F_Eu.OAu.O <- F_Eu.OAvA.O <- F_Eu.OAvP.O <- F_Eu.OAvC.O <- 
      F_Eu.OAwA.O <- F_Eu.OAwP.O <- F_Eu.OAwC.O <- F_Eu.OAbP.O <- beta.O*omega_a*(contact.Su + 
                                                                                    reinfec.DO*contact.Ru.D+
                                                                                    reinfec.OO*contact.Ru.O)
    # EuDxI cols
    F_Eu.OIu.O <- F_Eu.OIvA.O <- F_Eu.OIvP.O <- F_Eu.OIvC.O <- 
      F_Eu.OIwA.O <- F_Eu.OIwP.O <- F_Eu.OIwC.O <- F_Eu.OIbP.O <- beta.O*(contact.Su + 
                                                                            reinfec.DO*contact.Ru.D+
                                                                            reinfec.OO*contact.Ru.O)
    # EuDxH cols
    F_Eu.OHu.O <- F_Eu.OHvA.O <- F_Eu.OHvP.O <- F_Eu.OHvC.O <- 
      F_Eu.OHwA.O <- F_Eu.OHwP.O <- F_Eu.OHwC.O <- F_Eu.OHbP.O <- beta.O*(contact.Su + 
                                                                            reinfec.DO*contact.Ru.D+
                                                                            reinfec.OO*contact.Ru.O)
    # EuD lines
    F_Eu.O <- cbind(F_Eu.OEu.O, F_Eu.OEvA.O, F_Eu.OEvP.O, F_Eu.OEvC.O, F_Eu.OEwA.O, F_Eu.OEwP.O, F_Eu.OEwC.O, F_Eu.OEbP.O,
                    F_Eu.OAu.O, F_Eu.OAvA.O, F_Eu.OAvP.O, F_Eu.OAvC.O, F_Eu.OAwA.O, F_Eu.OAwP.O, F_Eu.OAwC.O, F_Eu.OAbP.O,
                    F_Eu.OIu.O, F_Eu.OIvA.O, F_Eu.OIvP.O, F_Eu.OIvC.O, F_Eu.OIwA.O, F_Eu.OIwP.O, F_Eu.OIwC.O, F_Eu.OIbP.O,
                    F_Eu.OHu.O, F_Eu.OHvA.O, F_Eu.OHvP.O, F_Eu.OHvC.O, F_Eu.OHwA.O, F_Eu.OHwP.O, F_Eu.OHwC.O, F_Eu.OHbP.O)
    
    # EvA.OxE cols
    F_EvA.OEu.O <- F_EvA.OEvA.O <- F_EvA.OEvP.O <- F_EvA.OEvC.O <- 
      F_EvA.OEwA.O <- F_EvA.OEwP.O <- F_EvA.OEwC.O <- F_EvA.OEbP.O <- beta_vA.O*omega*(contact.SvA + 
                                                                                         reinfec.DO*contact.RvA.D+
                                                                                         reinfec.OO*contact.RvA.O)
    #EvA.OxA cols    
    F_EvA.OAu.O <- F_EvA.OAvA.O <- F_EvA.OAvP.O <- F_EvA.OAvC.O <- 
      F_EvA.OAwA.O <- F_EvA.OAwP.O <- F_EvA.OAwC.O <- F_EvA.OAbP.O <- beta_vA.O*omega_a*(contact.SvA + 
                                                                                           reinfec.DO*contact.RvA.D+
                                                                                           reinfec.OO*contact.RvA.O)
    # EvA.OxI cols
    F_EvA.OIu.O <- F_EvA.OIvA.O <- F_EvA.OIvP.O <- F_EvA.OIvC.O <- 
      F_EvA.OIwA.O <- F_EvA.OIwP.O <- F_EvA.OIwC.O <- F_EvA.OIbP.O <- beta_vA.O*(contact.SvA + 
                                                                                   reinfec.DO*contact.RvA.D+
                                                                                   reinfec.OO*contact.RvA.O)
    #EvA.OxH cols
    F_EvA.OHu.O <- F_EvA.OHvA.O <- F_EvA.OHvP.O <- F_EvA.OHvC.O <- 
      F_EvA.OHwA.O <- F_EvA.OHwP.O <- F_EvA.OHwC.O <- F_EvA.OHbP.O <- beta_vA.O*omega_s*(contact.SvA + 
                                                                                           reinfec.DO*contact.RvA.D+
                                                                                           reinfec.OO*contact.RvA.O)
    #F_EvA.O lines
    F_EvA.O <- cbind(F_EvA.OEu.O, F_EvA.OEvA.O, F_EvA.OEvP.O, F_EvA.OEvC.O, F_EvA.OEwA.O, F_EvA.OEwP.O, F_EvA.OEwC.O, F_EvA.OEbP.O,
                     F_EvA.OAu.O, F_EvA.OAvA.O, F_EvA.OAvP.O, F_EvA.OAvC.O, F_EvA.OAwA.O, F_EvA.OAwP.O, F_EvA.OAwC.O, F_EvA.OAbP.O,
                     F_EvA.OIu.O, F_EvA.OIvA.O, F_EvA.OIvP.O, F_EvA.OIvC.O, F_EvA.OIwA.O, F_EvA.OIwP.O, F_EvA.OIwC.O, F_EvA.OIbP.O,
                     F_EvA.OHu.O, F_EvA.OHvA.O, F_EvA.OHvP.O, F_EvA.OHvC.O, F_EvA.OHwA.O, F_EvA.OHwP.O, F_EvA.OHwC.O, F_EvA.OHbP.O)
    
    #EvP.OxE cols
    F_EvP.OEu.O <- F_EvP.OEvA.O <- F_EvP.OEvP.O <- F_EvP.OEvC.O <- 
      F_EvP.OEwA.O <- F_EvP.OEwP.O <- F_EvP.OEwC.O <- F_EvP.OEbP.O <- beta_vP.O*omega*(contact.SvP + 
                                                                                         reinfec.DO*contact.RvP.D+
                                                                                         reinfec.OO*contact.RvP.O)
    #EvP.OxA cols
    F_EvP.OAu.O <- F_EvP.OAvA.O <- F_EvP.OAvP.O <- F_EvP.OAvC.O <- 
      F_EvP.OAwA.O <- F_EvP.OAwP.O <- F_EvP.OAwC.O <- F_EvP.OAbP.O <- beta_vP.O*omega_a*(contact.SvP + 
                                                                                           reinfec.DO*contact.RvP.D+
                                                                                           reinfec.OO*contact.RvP.O)
    #EvP.OxI cols
    F_EvP.OIu.O <- F_EvP.OIvA.O <- F_EvP.OIvP.O <- F_EvP.OIvC.O <- 
      F_EvP.OIwA.O <- F_EvP.OIwP.O <- F_EvP.OIwC.O <- F_EvP.OIbP.O <- beta_vP.O*(contact.SvP + 
                                                                                   reinfec.DO*contact.RvP.D+
                                                                                   reinfec.OO*contact.RvP.O)
    #EvP.OxH cols
    F_EvP.OHu.O <- F_EvP.OHvA.O <- F_EvP.OHvP.O <- F_EvP.OHvC.O <- 
      F_EvP.OHwA.O <- F_EvP.OHwP.O <- F_EvP.OHwC.O <- F_EvP.OHbP.O <- beta_vP.O*omega_s*(contact.SvP + 
                                                                                           reinfec.DO*contact.RvP.D+
                                                                                           reinfec.OO*contact.RvP.O)
    #EvP.O lines
    F_EvP.O <- cbind(F_EvP.OEu.O, F_EvP.OEvA.O, F_EvP.OEvP.O, F_EvP.OEvC.O, F_EvP.OEwA.O, F_EvP.OEwP.O, F_EvP.OEwC.O, F_EvP.OEbP.O,
                     F_EvP.OAu.O, F_EvP.OAvA.O, F_EvP.OAvP.O, F_EvP.OAvC.O, F_EvP.OAwA.O, F_EvP.OAwP.O, F_EvP.OAwC.O, F_EvP.OAbP.O,
                     F_EvP.OIu.O, F_EvP.OIvA.O, F_EvP.OIvP.O, F_EvP.OIvC.O, F_EvP.OIwA.O, F_EvP.OIwP.O, F_EvP.OIwC.O, F_EvP.OIbP.O,
                     F_EvP.OHu.O, F_EvP.OHvA.O, F_EvP.OHvP.O, F_EvP.OHvC.O, F_EvP.OHwA.O, F_EvP.OHwP.O, F_EvP.OHwC.O, F_EvP.OHbP.O)
    #EvC.OxE cols
    F_EvC.OEu.O <- F_EvC.OEvA.O <- F_EvC.OEvP.O <- F_EvC.OEvC.O <- 
      F_EvC.OEwA.O <- F_EvC.OEwP.O <- F_EvC.OEwC.O <- F_EvC.OEbP.O <- beta_vC.O*omega*(contact.SvC + 
                                                                                         reinfec.DO*contact.RvC.D+
                                                                                         reinfec.OO*contact.RvC.O)
    #EvC.OxA cols
    F_EvC.OAu.O <- F_EvC.OAvA.O <- F_EvC.OAvP.O <- F_EvC.OAvC.O <- 
      F_EvC.OAwA.O <- F_EvC.OAwP.O <- F_EvC.OAwC.O <- F_EvC.OAbP.O <- beta_vC.O*omega_a*(contact.SvC + 
                                                                                           reinfec.DO*contact.RvC.D+
                                                                                           reinfec.OO*contact.RvC.O)
    #EvC.OxI cols
    F_EvC.OIu.O <- F_EvC.OIvA.O <- F_EvC.OIvP.O <- F_EvC.OIvC.O <- 
      F_EvC.OIwA.O <- F_EvC.OIwP.O <- F_EvC.OIwC.O <- F_EvC.OIbP.O <- beta_vC.O*(contact.SvC + 
                                                                                   reinfec.DO*contact.RvC.D+
                                                                                   reinfec.OO*contact.RvC.O)
    #EvC.OxH cols
    F_EvC.OHu.O <- F_EvC.OHvA.O <- F_EvC.OHvP.O <- F_EvC.OHvC.O <- 
      F_EvC.OHwA.O <- F_EvC.OHwP.O <- F_EvC.OHwC.O <- F_EvC.OHbP.O <- beta_vC.O*omega_s*(contact.SvC + 
                                                                                           reinfec.DO*contact.RvC.D+
                                                                                           reinfec.OO*contact.RvC.O)
    #EvC.O lines
    F_EvC.O <- cbind(F_EvC.OEu.O, F_EvC.OEvA.O, F_EvC.OEvP.O, F_EvC.OEvC.O, F_EvC.OEwA.O, F_EvC.OEwP.O, F_EvC.OEwC.O, F_EvC.OEbP.O,
                     F_EvC.OAu.O, F_EvC.OAvA.O, F_EvC.OAvP.O, F_EvC.OAvC.O, F_EvC.OAwA.O, F_EvC.OAwP.O, F_EvC.OAwC.O, F_EvC.OAbP.O,
                     F_EvC.OIu.O, F_EvC.OIvA.O, F_EvC.OIvP.O, F_EvC.OIvC.O, F_EvC.OIwA.O, F_EvC.OIwP.O, F_EvC.OIwC.O, F_EvC.OIbP.O,
                     F_EvC.OHu.O, F_EvC.OHvA.O, F_EvC.OHvP.O, F_EvC.OHvC.O, F_EvC.OHwA.O, F_EvC.OHwP.O, F_EvC.OHwC.O, F_EvC.OHbP.O)
    
    #F_EwA.OxE cols
    F_EwA.OEu.O <- F_EwA.OEvA.O <- F_EwA.OEvP.O <- F_EwA.OEvC.O <- 
      F_EwA.OEwA.O <- F_EwA.OEwP.O <- F_EwA.OEwC.O <- F_EwA.OEbP.O <- beta_wA.O*omega*(contact.SwA + 
                                                                                         reinfec.DO*contact.RwA.D+
                                                                                         reinfec.OO*contact.RwA.O)
    #F_EwA.OxA cols
    F_EwA.OAu.O <- F_EwA.OAvA.O <- F_EwA.OAvP.O <- F_EwA.OAvC.O <- 
      F_EwA.OAwA.O <- F_EwA.OAwP.O <- F_EwA.OAwC.O <- F_EwA.OAbP.O <- beta_wA.O*omega_a*(contact.SwA + 
                                                                                           reinfec.DO*contact.RwA.D+
                                                                                           reinfec.OO*contact.RwA.O)
    #F_EwA.OxI cols
    F_EwA.OIu.O <- F_EwA.OIvA.O <- F_EwA.OIvP.O <- F_EwA.OIvC.O <- 
      F_EwA.OIwA.O <- F_EwA.OIwP.O <- F_EwA.OIwC.O <- F_EwA.OIbP.O <- beta_wA.O*(contact.SwA + 
                                                                                   reinfec.DO*contact.RwA.D+
                                                                                   reinfec.OO*contact.RwA.O)
    #F_EwA.OxH cols
    F_EwA.OHu.O <- F_EwA.OHvA.O <- F_EwA.OHvP.O <- F_EwA.OHvC.O <- 
      F_EwA.OHwA.O <- F_EwA.OHwP.O <- F_EwA.OHwC.O <- F_EwA.OHbP.O <- beta_wA.O*omega_s*(contact.SwA + 
                                                                                           reinfec.DO*contact.RwA.D+
                                                                                           reinfec.OO*contact.RwA.O)
    #F_EwA.O lines
    F_EwA.O <- cbind(F_EwA.OEu.O, F_EwA.OEvA.O, F_EwA.OEvP.O, F_EwA.OEvC.O, F_EwA.OEwA.O, F_EwA.OEwP.O, F_EwA.OEwC.O, F_EwA.OEbP.O,
                     F_EwA.OAu.O, F_EwA.OAvA.O, F_EwA.OAvP.O, F_EwA.OAvC.O, F_EwA.OAwA.O, F_EwA.OAwP.O, F_EwA.OAwC.O, F_EwA.OAbP.O,
                     F_EwA.OIu.O, F_EwA.OIvA.O, F_EwA.OIvP.O, F_EwA.OIvC.O, F_EwA.OIwA.O, F_EwA.OIwP.O, F_EwA.OIwC.O, F_EwA.OIbP.O,
                     F_EwA.OHu.O, F_EwA.OHvA.O, F_EwA.OHvP.O, F_EwA.OHvC.O, F_EwA.OHwA.O, F_EwA.OHwP.O, F_EwA.OHwC.O, F_EwA.OHbP.O)
    
    #F_EwP.OxE cols
    F_EwP.OEu.O <- F_EwP.OEvA.O <- F_EwP.OEvP.O <- F_EwP.OEvC.O <- 
      F_EwP.OEwA.O <- F_EwP.OEwP.O <- F_EwP.OEwC.O <- F_EwP.OEbP.O <- beta_wP.O*omega*(contact.SwP + 
                                                                                         reinfec.DO*contact.RwP.D+
                                                                                         reinfec.OO*contact.RwP.O)
    #F_EwP.OxA cols
    F_EwP.OAu.O <- F_EwP.OAvA.O <- F_EwP.OAvP.O <- F_EwP.OAvC.O <- 
      F_EwP.OAwA.O <- F_EwP.OAwP.O <- F_EwP.OAwC.O <- F_EwP.OAbP.O <- beta_wP.O*omega_a*(contact.SwP + 
                                                                                           reinfec.DO*contact.RwP.D+
                                                                                           reinfec.OO*contact.RwP.O)
    #F_EwP.OxI cols
    F_EwP.OIu.O <- F_EwP.OIvA.O <- F_EwP.OIvP.O <- F_EwP.OIvC.O <- 
      F_EwP.OIwA.O <- F_EwP.OIwP.O <- F_EwP.OIwC.O <- F_EwP.OIbP.O <- beta_wP.O*(contact.SwP + 
                                                                                   reinfec.DO*contact.RwP.D+
                                                                                   reinfec.OO*contact.RwP.O)
    #F_EwP.OxH cols
    F_EwP.OHu.O <- F_EwP.OHvA.O <- F_EwP.OHvP.O <- F_EwP.OHvC.O <- 
      F_EwP.OHwA.O <- F_EwP.OHwP.O <- F_EwP.OHwC.O <- F_EwP.OHbP.O <- beta_wP.O*omega_s*(contact.SwP + 
                                                                                           reinfec.DO*contact.RwP.D+
                                                                                           reinfec.OO*contact.RwP.O)
    #F_EwP.O lines
    F_EwP.O <- cbind(F_EwP.OEu.O, F_EwP.OEvA.O, F_EwP.OEvP.O, F_EwP.OEvC.O, F_EwP.OEwA.O, F_EwP.OEwP.O, F_EwP.OEwC.O, F_EwP.OEbP.O,
                     F_EwP.OAu.O, F_EwP.OAvA.O, F_EwP.OAvP.O, F_EwP.OAvC.O, F_EwP.OAwA.O, F_EwP.OAwP.O, F_EwP.OAwC.O, F_EwP.OAbP.O,
                     F_EwP.OIu.O, F_EwP.OIvA.O, F_EwP.OIvP.O, F_EwP.OIvC.O, F_EwP.OIwA.O, F_EwP.OIwP.O, F_EwP.OIwC.O, F_EwP.OIbP.O,
                     F_EwP.OHu.O, F_EwP.OHvA.O, F_EwP.OHvP.O, F_EwP.OHvC.O, F_EwP.OHwA.O, F_EwP.OHwP.O, F_EwP.OHwC.O, F_EwP.OHbP.O)
    
    #F_EwC.OxE cols
    F_EwC.OEu.O <- F_EwC.OEvA.O <- F_EwC.OEvP.O <- F_EwC.OEvC.O <- 
      F_EwC.OEwA.O <- F_EwC.OEwP.O <- F_EwC.OEwC.O <- F_EwC.OEbP.O <- beta_wC.O*omega*(contact.SwC + 
                                                                                         reinfec.DO*contact.RwC.D+
                                                                                         reinfec.OO*contact.RwC.O)
    #F_EwC.OxA cols
    F_EwC.OAu.O <- F_EwC.OAvA.O <- F_EwC.OAvP.O <- F_EwC.OAvC.O <- 
      F_EwC.OAwA.O <- F_EwC.OAwP.O <- F_EwC.OAwC.O <- F_EwC.OAbP.O <- beta_wC.O*omega_a*(contact.SwC + 
                                                                                           reinfec.DO*contact.RwC.D+
                                                                                           reinfec.OO*contact.RwC.O)
    #F_EwC.OxI cols
    F_EwC.OIu.O <- F_EwC.OIvA.O <- F_EwC.OIvP.O <- F_EwC.OIvC.O <- 
      F_EwC.OIwA.O <- F_EwC.OIwP.O <- F_EwC.OIwC.O <- F_EwC.OIbP.O <- beta_wC.O*(contact.SwC + 
                                                                                   reinfec.DO*contact.RwC.D+
                                                                                   reinfec.OO*contact.RwC.O)
    #F_EwC.OxH cols
    F_EwC.OHu.O <- F_EwC.OHvA.O <- F_EwC.OHvP.O <- F_EwC.OHvC.O <- 
      F_EwC.OHwA.O <- F_EwC.OHwP.O <- F_EwC.OHwC.O <- F_EwC.OHbP.O <- beta_wC.O*omega_s*(contact.SwC + 
                                                                                           reinfec.DO*contact.RwC.D+
                                                                                           reinfec.OO*contact.RwC.O)
    #F_EwC.O lines
    F_EwC.O <- cbind(F_EwC.OEu.O, F_EwC.OEvA.O, F_EwC.OEvP.O, F_EwC.OEvC.O, F_EwC.OEwA.O, F_EwC.OEwP.O, F_EwC.OEwC.O, F_EwC.OEbP.O,
                     F_EwC.OAu.O, F_EwC.OAvA.O, F_EwC.OAvP.O, F_EwC.OAvC.O, F_EwC.OAwA.O, F_EwC.OAwP.O, F_EwC.OAwC.O, F_EwC.OAbP.O,
                     F_EwC.OIu.O, F_EwC.OIvA.O, F_EwC.OIvP.O, F_EwC.OIvC.O, F_EwC.OIwA.O, F_EwC.OIwP.O, F_EwC.OIwC.O, F_EwC.OIbP.O,
                     F_EwC.OHu.O, F_EwC.OHvA.O, F_EwC.OHvP.O, F_EwC.OHvC.O, F_EwC.OHwA.O, F_EwC.OHwP.O, F_EwC.OHwC.O, F_EwC.OHbP.O)
    
    #F_EwC.OxE cols
    F_EbP.OEu.O <- F_EbP.OEvA.O <- F_EbP.OEvP.O <- F_EbP.OEvC.O <- 
      F_EbP.OEwA.O <- F_EbP.OEwP.O <- F_EbP.OEwC.O <- F_EbP.OEbP.O <- beta_bP.O*omega*(contact.SbP + 
                                                                                         reinfec.DO*contact.RbP.D+
                                                                                         reinfec.DD*contact.RbP.O)
    #F_EbP.OxA cols
    F_EbP.OAu.O <- F_EbP.OAvA.O <- F_EbP.OAvP.O <- F_EbP.OAvC.O <- 
      F_EbP.OAwA.O <- F_EbP.OAwP.O <- F_EbP.OAwC.O <- F_EbP.OAbP.O <- beta_bP.O*omega_a*(contact.SbP + 
                                                                                           reinfec.DO*contact.RbP.D+
                                                                                           reinfec.OO*contact.RbP.O)
    #F_EbP.OxI cols
    F_EbP.OIu.O <- F_EbP.OIvA.O <- F_EbP.OIvP.O <- F_EbP.OIvC.O <- 
      F_EbP.OIwA.O <- F_EbP.OIwP.O <- F_EbP.OIwC.O <- F_EbP.OIbP.O <- beta_bP.O*(contact.SbP + 
                                                                                   reinfec.DO*contact.RbP.D+
                                                                                   reinfec.OO*contact.RbP.O)
    #F_EbP.OxH cols
    F_EbP.OHu.O <- F_EbP.OHvA.O <- F_EbP.OHvP.O <- F_EbP.OHvC.O <- 
      F_EbP.OHwA.O <- F_EbP.OHwP.O <- F_EbP.OHwC.O <- F_EbP.OHbP.O <- beta_bP.O*omega_s*(contact.SbP + 
                                                                                           reinfec.DO*contact.RbP.D+
                                                                                           reinfec.OO*contact.RbP.O)
    #F_EbP.O lines
    F_EbP.O <- cbind(F_EbP.OEu.O, F_EbP.OEvA.O, F_EbP.OEvP.O, F_EbP.OEvC.O, F_EbP.OEwA.O, F_EbP.OEwP.O, F_EbP.OEwC.O, F_EbP.OEbP.O,
                     F_EbP.OAu.O, F_EbP.OAvA.O, F_EbP.OAvP.O, F_EbP.OAvC.O, F_EbP.OAwA.O, F_EbP.OAwP.O, F_EbP.OAwC.O, F_EbP.OAbP.O,
                     F_EbP.OIu.O, F_EbP.OIvA.O, F_EbP.OIvP.O, F_EbP.OIvC.O, F_EbP.OIwA.O, F_EbP.OIwP.O, F_EbP.OIwC.O, F_EbP.OIbP.O,
                     F_EbP.OHu.O, F_EbP.OHvA.O, F_EbP.OHvP.O, F_EbP.OHvC.O, F_EbP.OHwA.O, F_EbP.OHwP.O, F_EbP.OHwC.O, F_EbP.OHbP.O)
    #all E's together
    F_E.O <- rbind(F_Eu.O, F_EvA.O, F_EvP.O, F_EvC.O, F_EwA.O, F_EwP.O, F_EwC.O, F_EbP.O)
    # other lines are zero
    F_mat.O <- rbind(F_E.O, matrix(0, nrow = 3*nrow(F_E.O),ncol = ncol(F_E.O)))
    
    ################################ T matrix ################################## 
    ## Eu
    # Eu.OxE cols
    T_Eu.OEu.O <- diag(1-gamma.O, nrow = age.bins, ncol = age.bins)
    T_Eu.OEvA.O <- T_Eu.OEvP.O <- T_Eu.OEvC.O <- 
      T_Eu.OEwA.O <- T_Eu.OEwP.O <- T_Eu.OEwC.O <- T_Eu.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EuDxI cols
    T_Eu.OIu.O <- T_Eu.OIvA.O <- T_Eu.OIvP.O <- T_Eu.OIvC.O <- 
      T_Eu.OIwA.O <- T_Eu.OIwP.O <- T_Eu.OIwC.O <- T_Eu.OIbP.O  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EuDxA cols
    T_Eu.OAu.O <- T_Eu.OAvA.O <- T_Eu.OAvP.O <- T_Eu.OAvC.O <- 
      T_Eu.OAwA.O <- T_Eu.OAwP.O <- T_Eu.OAwC.O <- T_Eu.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EuDxH cols
    T_Eu.OHu.O <- T_Eu.OHvA.O <- T_Eu.OHvP.O <- T_Eu.OHvC.O <- 
      T_Eu.OHwA.O <- T_Eu.OHwP.O <- T_Eu.OHwC.O <- T_Eu.OHbP.O  <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_Eu.O <- cbind(T_Eu.OEu.O, T_Eu.OEvA.O, T_Eu.OEvP.O, T_Eu.OEvC.O, T_Eu.OEwA.O, T_Eu.OEwP.O, T_Eu.OEwC.O, T_Eu.OEbP.O,
                    T_Eu.OAu.O, T_Eu.OAvA.O, T_Eu.OAvP.O, T_Eu.OAvC.O, T_Eu.OAwA.O, T_Eu.OAwP.O, T_Eu.OAwC.O, T_Eu.OAbP.O,
                    T_Eu.OIu.O, T_Eu.OIvA.O, T_Eu.OIvP.O, T_Eu.OIvC.O, T_Eu.OIwA.O, T_Eu.OIwP.O, T_Eu.OIwC.O, T_Eu.OIbP.O,
                    T_Eu.OHu.O, T_Eu.OHvA.O, T_Eu.OHvP.O, T_Eu.OHvC.O, T_Eu.OHwA.O, T_Eu.OHwP.O, T_Eu.OHwC.O, T_Eu.OHbP.O)
    ## EvA
    # EvA.OxE cols
    T_EvA.OEvA.O <- diag(1-gamma.O, nrow = age.bins, ncol = age.bins)
    T_EvA.OEu.O <- T_EvA.OEvP.O <- T_EvA.OEvC.O <- 
      T_EvA.OEwA.O <- T_EvA.OEwP.O <- T_EvA.OEwC.O <- T_EvA.OEbP.O  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvA.OxA cols
    T_EvA.OAu.O <- T_EvA.OAvA.O <- T_EvA.OAvP.O <- T_EvA.OAvC.O <- 
      T_EvA.OAwA.O <- T_EvA.OAwP.O <- T_EvA.OAwC.O <- T_EvA.OAbP.O  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvA.OxI cols
    T_EvA.OIu.O <- T_EvA.OIvA.O <- T_EvA.OIvP.O <- T_EvA.OIvC.O <- 
      T_EvA.OIwA.O <- T_EvA.OIwP.O <- T_EvA.OIwC.O <- T_EvA.OIbP.O  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvAxH cols
    T_EvA.OHu.O <- T_EvA.OHvA.O <- T_EvA.OHvP.O <- T_EvA.OHvC.O <- 
      T_EvA.OHwA.O <- T_EvA.OHwP.O <- T_EvA.OHwC.O <- T_EvA.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EvA.O <- cbind(T_EvA.OEu.O, T_EvA.OEvA.O, T_EvA.OEvP.O, T_EvA.OEvC.O, T_EvA.OEwA.O, T_EvA.OEwP.O, T_EvA.OEwC.O, T_EvA.OEbP.O,
                     T_EvA.OAu.O, T_EvA.OAvA.O, T_EvA.OAvP.O, T_EvA.OAvC.O, T_EvA.OAwA.O, T_EvA.OAwP.O, T_EvA.OAwC.O, T_EvA.OAbP.O,
                     T_EvA.OIu.O, T_EvA.OIvA.O, T_EvA.OIvP.O, T_EvA.OIvC.O, T_EvA.OIwA.O, T_EvA.OIwP.O, T_EvA.OIwC.O, T_EvA.OIbP.O,
                     T_EvA.OHu.O, T_EvA.OHvA.O, T_EvA.OHvP.O, T_EvA.OHvC.O, T_EvA.OHwA.O, T_EvA.OHwP.O, T_EvA.OHwC.O, T_EvA.OHbP.O)
    
    ## EvP
    # EvP.OxE cols
    T_EvP.OEvP.O <- diag(1-gamma.O, nrow = age.bins, ncol = age.bins)
    T_EvP.OEu.O <- T_EvP.OEvA.O <- T_EvP.OEvC.O <- 
      T_EvP.OEwA.O <- T_EvP.OEwP.O <- T_EvP.OEwC.O <- T_EvP.OEbP.O  <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvP.OxA cols
    T_EvP.OAu.O <- T_EvP.OAvA.O <- T_EvP.OAvP.O <- T_EvP.OAvC.O <- 
      T_EvP.OAwA.O <- T_EvP.OAwP.O <- T_EvP.OAwC.O <- T_EvP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvP.OxI cols
    T_EvP.OIu.O <- T_EvP.OIvA.O <- T_EvP.OIvP.O <- T_EvP.OIvC.O <- 
      T_EvP.OIwA.O <- T_EvP.OIwP.O <- T_EvP.OIwC.O <- T_EvP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvPxH cols
    T_EvP.OHu.O <- T_EvP.OHvA.O <- T_EvP.OHvP.O <- T_EvP.OHvC.O <- 
      T_EvP.OHwA.O <- T_EvP.OHwP.O <- T_EvP.OHwC.O <- T_EvP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EvP.O <- cbind(T_EvP.OEu.O, T_EvP.OEvA.O, T_EvP.OEvP.O, T_EvP.OEvC.O, T_EvP.OEwA.O, T_EvP.OEwP.O, T_EvP.OEwC.O, T_EvP.OEbP.O,
                     T_EvP.OAu.O, T_EvP.OAvA.O, T_EvP.OAvP.O, T_EvP.OAvC.O, T_EvP.OAwA.O, T_EvP.OAwP.O, T_EvP.OAwC.O, T_EvP.OAbP.O,
                     T_EvP.OIu.O, T_EvP.OIvA.O, T_EvP.OIvP.O, T_EvP.OIvC.O, T_EvP.OIwA.O, T_EvP.OIwP.O, T_EvP.OIwC.O, T_EvP.OIbP.O,
                     T_EvP.OHu.O, T_EvP.OHvA.O, T_EvP.OHvP.O, T_EvP.OHvC.O, T_EvP.OHwA.O, T_EvP.OHwP.O, T_EvP.OHwC.O, T_EvP.OHbP.O)
    
    ## EvC
    # EvC.OxE cols
    T_EvC.OEvC.O <- diag(1-gamma.O, nrow = age.bins, ncol = age.bins)
    T_EvC.OEu.O <- T_EvC.OEvA.O <- T_EvC.OEvP.O <- 
      T_EvC.OEwA.O <- T_EvC.OEwP.O <- T_EvC.OEwC.O <- T_EvC.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvC.OxA cols
    T_EvC.OAu.O <- T_EvC.OAvA.O <- T_EvC.OAvP.O <- T_EvC.OAvC.O <- 
      T_EvC.OAwA.O <- T_EvC.OAwP.O <- T_EvC.OAwC.O <- T_EvC.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvC.OxI cols
    T_EvC.OIu.O <- T_EvC.OIvA.O <- T_EvC.OIvP.O <- T_EvC.OIvC.O <- 
      T_EvC.OIwA.O <- T_EvC.OIwP.O <- T_EvC.OIwC.O <- T_EvC.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EvCxH cols
    T_EvC.OHu.O <- T_EvC.OHvA.O <- T_EvC.OHvP.O <- T_EvC.OHvC.O <- 
      T_EvC.OHwA.O <- T_EvC.OHwP.O <- T_EvC.OHwC.O <- T_EvC.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EvC.O <- cbind(T_EvC.OEu.O, T_EvC.OEvA.O, T_EvC.OEvP.O, T_EvC.OEvC.O, T_EvC.OEwA.O, T_EvC.OEwP.O, T_EvC.OEwC.O, T_EvC.OEbP.O,
                     T_EvC.OAu.O, T_EvC.OAvA.O, T_EvC.OAvP.O, T_EvC.OAvC.O, T_EvC.OAwA.O, T_EvC.OAwP.O, T_EvC.OAwC.O, T_EvC.OAbP.O,
                     T_EvC.OIu.O, T_EvC.OIvA.O, T_EvC.OIvP.O, T_EvC.OIvC.O, T_EvC.OIwA.O, T_EvC.OIwP.O, T_EvC.OIwC.O, T_EvC.OIbP.O,
                     T_EvC.OHu.O, T_EvC.OHvA.O, T_EvC.OHvP.O, T_EvC.OHvC.O, T_EvC.OHwA.O, T_EvC.OHwP.O, T_EvC.OHwC.O, T_EvC.OHbP.O)
    
    ## EwA
    # EwA.OxE cols
    T_EwA.OEwA.O <- diag(1-gamma.O, nrow = age.bins, ncol = age.bins)
    T_EwA.OEu.O <- T_EwA.OEvA.O <- T_EwA.OEvP.O <- 
      T_EwA.OEvC.O <- T_EwA.OEwP.O <- T_EwA.OEwC.O <- T_EwA.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwA.OxA cols
    T_EwA.OAu.O <- T_EwA.OAvA.O <- T_EwA.OAvP.O <- T_EwA.OAvC.O <- 
      T_EwA.OAwA.O <- T_EwA.OAwP.O <- T_EwA.OAwC.O <- T_EwA.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwA.OxI cols
    T_EwA.OIu.O <- T_EwA.OIvA.O <- T_EwA.OIvP.O <- T_EwA.OIvC.O <- 
      T_EwA.OIwA.O <- T_EwA.OIwP.O <- T_EwA.OIwC.O <- T_EwA.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwAxH cols
    T_EwA.OHu.O <- T_EwA.OHvA.O <- T_EwA.OHvP.O <- T_EwA.OHvC.O <- 
      T_EwA.OHwA.O <- T_EwA.OHwP.O <- T_EwA.OHwC.O <- T_EwA.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EwA.O <- cbind(T_EwA.OEu.O, T_EwA.OEvA.O, T_EwA.OEvP.O, T_EwA.OEvC.O, T_EwA.OEwA.O, T_EwA.OEwP.O, T_EwA.OEwC.O, T_EwA.OEbP.O,
                     T_EwA.OAu.O, T_EwA.OAvA.O, T_EwA.OAvP.O, T_EwA.OAvC.O, T_EwA.OAwA.O, T_EwA.OAwP.O, T_EwA.OAwC.O, T_EwA.OAbP.O,
                     T_EwA.OIu.O, T_EwA.OIvA.O, T_EwA.OIvP.O, T_EwA.OIvC.O, T_EwA.OIwA.O, T_EwA.OIwP.O, T_EwA.OIwC.O, T_EwA.OIbP.O,
                     T_EwA.OHu.O, T_EwA.OHvA.O, T_EwA.OHvP.O, T_EwA.OHvC.O, T_EwA.OHwA.O, T_EwA.OHwP.O, T_EwA.OHwC.O, T_EwA.OHbP.O)
    
    ## EwP
    # EwP.OxE cols
    T_EwP.OEwP.O <- diag(1-gamma.O, nrow = age.bins, ncol = age.bins)
    T_EwP.OEu.O <- T_EwP.OEvA.O <- T_EwP.OEvP.O <- 
      T_EwP.OEvC.O <- T_EwP.OEwA.O <- T_EwP.OEwC.O <- T_EwP.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwP.OxA cols
    T_EwP.OAu.O <- T_EwP.OAvA.O <- T_EwP.OAvP.O <- T_EwP.OAvC.O <- 
      T_EwP.OAwA.O <- T_EwP.OAwP.O <- T_EwP.OAwC.O <- T_EwP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwP.OxI cols
    T_EwP.OIu.O <- T_EwP.OIvA.O <- T_EwP.OIvP.O <- T_EwP.OIvC.O <- 
      T_EwP.OIwA.O <- T_EwP.OIwP.O <- T_EwP.OIwC.O <- T_EwP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwPxH cols
    T_EwP.OHu.O <- T_EwP.OHvA.O <- T_EwP.OHvP.O <- T_EwP.OHvC.O <- 
      T_EwP.OHwA.O <- T_EwP.OHwP.O <- T_EwP.OHwC.O <- T_EwP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EwP.O <- cbind(T_EwP.OEu.O, T_EwP.OEvA.O, T_EwP.OEvP.O, T_EwP.OEvC.O, T_EwP.OEwA.O, T_EwP.OEwP.O, T_EwP.OEwC.O, T_EwP.OEbP.O,
                     T_EwP.OAu.O, T_EwP.OAvA.O, T_EwP.OAvP.O, T_EwP.OAvC.O, T_EwP.OAwA.O, T_EwP.OAwP.O, T_EwP.OAwC.O, T_EwP.OAbP.O,
                     T_EwP.OIu.O, T_EwP.OIvA.O, T_EwP.OIvP.O, T_EwP.OIvC.O, T_EwP.OIwA.O, T_EwP.OIwP.O, T_EwP.OIwC.O, T_EwP.OIbP.O,
                     T_EwP.OHu.O, T_EwP.OHvA.O, T_EwP.OHvP.O, T_EwP.OHvC.O, T_EwP.OHwA.O, T_EwP.OHwP.O, T_EwP.OHwC.O, T_EwP.OHbP.O)
    
    ## EwC
    # EwC.OxE cols
    T_EwC.OEwC.O <- diag(1-gamma.O, nrow = age.bins, ncol = age.bins)
    T_EwC.OEu.O <- T_EwC.OEvA.O <- T_EwC.OEvP.O <- 
      T_EwC.OEvC.O <- T_EwC.OEwA.O <- T_EwC.OEwP.O <- T_EwC.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.OxA cols
    T_EwC.OAu.O <- T_EwC.OAvA.O <- T_EwC.OAvP.O <- T_EwC.OAvC.O <- 
      T_EwC.OAwA.O <- T_EwC.OAwP.O <- T_EwC.OAwC.O <- T_EwC.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.OxI cols
    T_EwC.OIu.O <- T_EwC.OIvA.O <- T_EwC.OIvP.O <- T_EwC.OIvC.O <- 
      T_EwC.OIwA.O <- T_EwC.OIwP.O <- T_EwC.OIwC.O <- T_EwC.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.OxH cols
    T_EwC.OHu.O <- T_EwC.OHvA.O <- T_EwC.OHvP.O <- T_EwC.OHvC.O <- 
      T_EwC.OHwA.O <- T_EwC.OHwP.O <- T_EwC.OHwC.O <- T_EwC.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EwC.O <- cbind(T_EwC.OEu.O, T_EwC.OEvA.O, T_EwC.OEvP.O, T_EwC.OEvC.O, T_EwC.OEwA.O, T_EwC.OEwP.O, T_EwC.OEwC.O, T_EwC.OEbP.O,
                     T_EwC.OAu.O, T_EwC.OAvA.O, T_EwC.OAvP.O, T_EwC.OAvC.O, T_EwC.OAwA.O, T_EwC.OAwP.O, T_EwC.OAwC.O, T_EwC.OAbP.O,
                     T_EwC.OIu.O, T_EwC.OIvA.O, T_EwC.OIvP.O, T_EwC.OIvC.O, T_EwC.OIwA.O, T_EwC.OIwP.O, T_EwC.OIwC.O, T_EwC.OIbP.O,
                     T_EwC.OHu.O, T_EwC.OHvA.O, T_EwC.OHvP.O, T_EwC.OHvC.O, T_EwC.OHwA.O, T_EwC.OHwP.O, T_EwC.OHwC.O, T_EwC.OHbP.O)
    ## EwC
    # EwC.OxE cols
    T_EbP.OEbP.O <- diag(1-gamma.O, nrow = age.bins, ncol = age.bins)
    T_EbP.OEu.O <- T_EbP.OEvA.O <- T_EbP.OEvP.O <- 
      T_EbP.OEvC.O <- T_EbP.OEwA.O <- T_EbP.OEwP.O <- T_EbP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.OxA cols
    T_EbP.OAu.O <- T_EbP.OAvA.O <- T_EbP.OAvP.O <- T_EbP.OAvC.O <- 
      T_EbP.OAwA.O <- T_EbP.OAwP.O <- T_EbP.OAwC.O <- T_EbP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.OxI cols
    T_EbP.OIu.O <- T_EbP.OIvA.O <- T_EbP.OIvP.O <- T_EbP.OIvC.O <- 
      T_EbP.OIwA.O <- T_EbP.OIwP.O <- T_EbP.OIwC.O <- T_EbP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # EwC.OxH cols
    T_EbP.OHu.O <- T_EbP.OHvA.O <- T_EbP.OHvP.O <- T_EbP.OHvC.O <- 
      T_EbP.OHwA.O <- T_EbP.OHwP.O <- T_EbP.OHwC.O <- T_EbP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_EbP.O <- cbind(T_EbP.OEu.O, T_EbP.OEvA.O, T_EbP.OEvP.O, T_EbP.OEvC.O, T_EbP.OEwA.O, T_EbP.OEwP.O, T_EbP.OEwC.O, T_EbP.OEbP.O,
                     T_EbP.OAu.O, T_EbP.OAvA.O, T_EbP.OAvP.O, T_EbP.OAvC.O, T_EbP.OAwA.O, T_EbP.OAwP.O, T_EbP.OAwC.O, T_EbP.OAbP.O,
                     T_EbP.OIu.O, T_EbP.OIvA.O, T_EbP.OIvP.O, T_EbP.OIvC.O, T_EbP.OIwA.O, T_EbP.OIwP.O, T_EbP.OIwC.O, T_EbP.OIbP.O,
                     T_EbP.OHu.O, T_EbP.OHvA.O, T_EbP.OHvP.O, T_EbP.OHvC.O, T_EbP.OHwA.O, T_EbP.OHwP.O, T_EbP.OHwC.O, T_EbP.OHbP.O)
    #all E's together
    T_E.O <- rbind(T_Eu.O, T_EvA.O, T_EvP.O, T_EvC.O, T_EwA.O, T_EwP.O, T_EwC.O, T_EbP.O)
    
    ## Au
    # Au.OxE cols
    T_Au.OEu.O <- diag((1-ihr.O)*gamma.O*asymp, ncol = age.bins, nrow = age.bins)
    T_Au.OEvA.O <- T_Au.OEvP.O <- 
      T_Au.OEvC.O <- T_Au.OEwA.O <- T_Au.OEwP.O <- T_Au.OEwC.O <- T_Au.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Au.OxA cols
    T_Au.OAu.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_Au.OAvA.O <- T_Au.OAvP.O <- T_Au.OAvC.O <- 
      T_Au.OAwA.O <- T_Au.OAwP.O <- T_Au.OAwC.O <- T_Au.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Au.OxI cols
    T_Au.OIu.O <- T_Au.OIvA.O <- T_Au.OIvP.O <- T_Au.OIvC.O <- 
      T_Au.OIwA.O <- T_Au.OIwP.O <- T_Au.OIwC.O <- T_Au.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Au.OxH cols
    T_Au.OHu.O <- T_Au.OHvA.O <- T_Au.OHvP.O <- T_Au.OHvC.O <- 
      T_Au.OHwA.O <- T_Au.OHwP.O <- T_Au.OHwC.O <- T_Au.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_Au.O <- cbind(T_Au.OEu.O, T_Au.OEvA.O, T_Au.OEvP.O, T_Au.OEvC.O, T_Au.OEwA.O, T_Au.OEwP.O, T_Au.OEwC.O, T_Au.OEbP.O,
                    T_Au.OAu.O, T_Au.OAvA.O, T_Au.OAvP.O, T_Au.OAvC.O, T_Au.OAwA.O, T_Au.OAwP.O, T_Au.OAwC.O, T_Au.OAbP.O,
                    T_Au.OIu.O, T_Au.OIvA.O, T_Au.OIvP.O, T_Au.OIvC.O, T_Au.OIwA.O, T_Au.OIwP.O, T_Au.OIwC.O, T_Au.OIbP.O,
                    T_Au.OHu.O, T_Au.OHvA.O, T_Au.OHvP.O, T_Au.OHvC.O, T_Au.OHwA.O, T_Au.OHwP.O, T_Au.OHwC.O, T_Au.OHbP.O)
    
    ## AvA
    # AvA.OxE cols
    T_AvA.OEvA.O <- diag((1-ihr_vA.O)*gamma.O*asymp_vA.O, ncol = age.bins, nrow = age.bins)
    T_AvA.OEu.O <- T_AvA.OEvP.O <- 
      T_AvA.OEvC.O <- T_AvA.OEwA.O <- T_AvA.OEwP.O <- T_AvA.OEwC.O <- T_AvA.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvA.OxA cols
    T_AvA.OAvA.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AvA.OAu.O <- T_AvA.OAvP.O <- T_AvA.OAvC.O <- 
      T_AvA.OAwA.O <- T_AvA.OAwP.O <- T_AvA.OAwC.O <- T_AvA.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvA.OxI cols
    T_AvA.OIu.O <- T_AvA.OIvA.O <- T_AvA.OIvP.O <- T_AvA.OIvC.O <- 
      T_AvA.OIwA.O <- T_AvA.OIwP.O <- T_AvA.OIwC.O <- T_AvA.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvA.OxH cols
    T_AvA.OHu.O <- T_AvA.OHvA.O <- T_AvA.OHvP.O <- T_AvA.OHvC.O <- 
      T_AvA.OHwA.O <- T_AvA.OHwP.O <- T_AvA.OHwC.O <- T_AvA.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AvA.O <- cbind(T_AvA.OEu.O, T_AvA.OEvA.O, T_AvA.OEvP.O, T_AvA.OEvC.O, T_AvA.OEwA.O, T_AvA.OEwP.O, T_AvA.OEwC.O, T_AvA.OEbP.O,
                     T_AvA.OAu.O, T_AvA.OAvA.O, T_AvA.OAvP.O, T_AvA.OAvC.O, T_AvA.OAwA.O, T_AvA.OAwP.O, T_AvA.OAwC.O, T_AvA.OAbP.O,
                     T_AvA.OIu.O, T_AvA.OIvA.O, T_AvA.OIvP.O, T_AvA.OIvC.O, T_AvA.OIwA.O, T_AvA.OIwP.O, T_AvA.OIwC.O, T_AvA.OIbP.O,
                     T_AvA.OHu.O, T_AvA.OHvA.O, T_AvA.OHvP.O, T_AvA.OHvC.O, T_AvA.OHwA.O, T_AvA.OHwP.O, T_AvA.OHwC.O, T_AvA.OHbP.O)
    
    ## AvP
    # AvP.OxE cols
    T_AvP.OEvP.O <- diag((1-ihr_vP.O)*gamma.O*asymp_vP.O, ncol = age.bins, nrow = age.bins)
    T_AvP.OEu.O <- T_AvP.OEvA.O <- 
      T_AvP.OEvC.O <- T_AvP.OEwA.O <- T_AvP.OEwP.O <- T_AvP.OEwC.O <- T_AvP.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvP.OxA cols
    T_AvP.OAvP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AvP.OAu.O <- T_AvP.OAvA.O <- T_AvP.OAvC.O <- 
      T_AvP.OAwA.O <- T_AvP.OAwP.O <- T_AvP.OAwC.O <- T_AvP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvP.OxI cols
    T_AvP.OIu.O <- T_AvP.OIvA.O <- T_AvP.OIvP.O <- T_AvP.OIvC.O <- 
      T_AvP.OIwA.O <- T_AvP.OIwP.O <- T_AvP.OIwC.O <- T_AvP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvP.OxH cols
    T_AvP.OHu.O <- T_AvP.OHvA.O <- T_AvP.OHvP.O <- T_AvP.OHvC.O <- 
      T_AvP.OHwA.O <- T_AvP.OHwP.O <- T_AvP.OHwC.O <- T_AvP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AvP.O <- cbind(T_AvP.OEu.O, T_AvP.OEvA.O, T_AvP.OEvP.O, T_AvP.OEvC.O, T_AvP.OEwA.O, T_AvP.OEwP.O, T_AvP.OEwC.O, T_AvP.OEbP.O,
                     T_AvP.OAu.O, T_AvP.OAvA.O, T_AvP.OAvP.O, T_AvP.OAvC.O, T_AvP.OAwA.O, T_AvP.OAwP.O, T_AvP.OAwC.O, T_AvP.OAbP.O,
                     T_AvP.OIu.O, T_AvP.OIvA.O, T_AvP.OIvP.O, T_AvP.OIvC.O, T_AvP.OIwA.O, T_AvP.OIwP.O, T_AvP.OIwC.O, T_AvP.OIbP.O,
                     T_AvP.OHu.O, T_AvP.OHvA.O, T_AvP.OHvP.O, T_AvP.OHvC.O, T_AvP.OHwA.O, T_AvP.OHwP.O, T_AvP.OHwC.O, T_AvP.OHbP.O)
    
    ## AvC
    # AvC.OxE cols
    T_AvC.OEvC.O <- diag((1-ihr_vC.O)*gamma.O*asymp_vC.O, ncol = age.bins, nrow = age.bins)
    T_AvC.OEu.O <- T_AvC.OEvA.O <- 
      T_AvC.OEvP.O <- T_AvC.OEwA.O <- T_AvC.OEwP.O <- T_AvC.OEwC.O <- T_AvC.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvC.OxA cols
    T_AvC.OAvC.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AvC.OAu.O <- T_AvC.OAvA.O <- T_AvC.OAvP.O <- 
      T_AvC.OAwA.O <- T_AvC.OAwP.O <- T_AvC.OAwC.O <- T_AvC.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvC.OxI cols
    T_AvC.OIu.O <- T_AvC.OIvA.O <- T_AvC.OIvP.O <- T_AvC.OIvC.O <- 
      T_AvC.OIwA.O <- T_AvC.OIwP.O <- T_AvC.OIwC.O <- T_AvC.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AvC.OxH cols
    T_AvC.OHu.O <- T_AvC.OHvA.O <- T_AvC.OHvP.O <- T_AvC.OHvC.O <- 
      T_AvC.OHwA.O <- T_AvC.OHwP.O <- T_AvC.OHwC.O <- T_AvC.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AvC.O <- cbind(T_AvC.OEu.O, T_AvC.OEvA.O, T_AvC.OEvP.O, T_AvC.OEvC.O, T_AvC.OEwA.O, T_AvC.OEwP.O, T_AvC.OEwC.O, T_AvC.OEbP.O,
                     T_AvC.OAu.O, T_AvC.OAvA.O, T_AvC.OAvP.O, T_AvC.OAvC.O, T_AvC.OAwA.O, T_AvC.OAwP.O, T_AvC.OAwC.O, T_AvC.OAbP.O,
                     T_AvC.OIu.O, T_AvC.OIvA.O, T_AvC.OIvP.O, T_AvC.OIvC.O, T_AvC.OIwA.O, T_AvC.OIwP.O, T_AvC.OIwC.O, T_AvC.OHbP.O,
                     T_AvC.OHu.O, T_AvC.OHvA.O, T_AvC.OHvP.O, T_AvC.OHvC.O, T_AvC.OHwA.O, T_AvC.OHwP.O, T_AvC.OHwC.O, T_AvC.OIbP.O)
    
    ## AwA
    # AwA.OxE cols
    T_AwA.OEwA.O <- diag((1-ihr_wA.O)*gamma.O*asymp_wA.O, ncol = age.bins, nrow = age.bins)
    T_AwA.OEu.O <- T_AwA.OEvP.O <- 
      T_AwA.OEvC.O <- T_AwA.OEvA.O <- T_AwA.OEwP.O <- T_AwA.OEwC.O <- T_AwA.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwA.OxA cols
    T_AwA.OAwA.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AwA.OAu.O <- T_AwA.OAvP.O <- T_AwA.OAvC.O <- 
      T_AwA.OAvA.O <- T_AwA.OAwP.O <- T_AwA.OAwC.O <- T_AwA.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwA.OxI cols
    T_AwA.OIu.O <- T_AwA.OIvA.O <- T_AwA.OIvP.O <- T_AwA.OIvC.O <- 
      T_AwA.OIwA.O <- T_AwA.OIwP.O <- T_AwA.OIwC.O <- T_AwA.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwA.OxH cols
    T_AwA.OHu.O <- T_AwA.OHvA.O <- T_AwA.OHvP.O <- T_AwA.OHvC.O <- 
      T_AwA.OHwA.O <- T_AwA.OHwP.O <- T_AwA.OHwC.O <- T_AwA.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AwA.O <- cbind(T_AwA.OEu.O, T_AwA.OEvA.O, T_AwA.OEvP.O, T_AwA.OEvC.O, T_AwA.OEwA.O, T_AwA.OEwP.O, T_AwA.OEwC.O, T_AwA.OEbP.O,
                     T_AwA.OAu.O, T_AwA.OAvA.O, T_AwA.OAvP.O, T_AwA.OAvC.O, T_AwA.OAwA.O, T_AwA.OAwP.O, T_AwA.OAwC.O, T_AwA.OAbP.O,
                     T_AwA.OIu.O, T_AwA.OIvA.O, T_AwA.OIvP.O, T_AwA.OIvC.O, T_AwA.OIwA.O, T_AwA.OIwP.O, T_AwA.OIwC.O, T_AwA.OIbP.O,
                     T_AwA.OHu.O, T_AwA.OHvA.O, T_AwA.OHvP.O, T_AwA.OHvC.O, T_AwA.OHwA.O, T_AwA.OHwP.O, T_AwA.OHwC.O, T_AwA.OHbP.O)
    
    ## AwP
    # AwP.OxE cols
    T_AwP.OEwP.O <- diag((1-ihr_wP.O)*gamma.O*asymp_wP.O, ncol = age.bins, nrow = age.bins)
    T_AwP.OEu.O <- T_AwP.OEvA.O <- 
      T_AwP.OEvC.O <- T_AwP.OEwA.O <- T_AwP.OEvP.O <- T_AwP.OEwC.O <- T_AwP.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwP.OxA cols
    T_AwP.OAwP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AwP.OAu.O <- T_AwP.OAvA.O <- T_AwP.OAvC.O <- 
      T_AwP.OAwA.O <- T_AwP.OAvP.O <- T_AwP.OAwC.O <- T_AwP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwP.OxI cols
    T_AwP.OIu.O <- T_AwP.OIvA.O <- T_AwP.OIvP.O <- T_AwP.OIvC.O <- 
      T_AwP.OIwA.O <- T_AwP.OIwP.O <- T_AwP.OIwC.O <- T_AwP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwP.OxH cols
    T_AwP.OHu.O <- T_AwP.OHvA.O <- T_AwP.OHvP.O <- T_AwP.OHvC.O <- 
      T_AwP.OHwA.O <- T_AwP.OHwP.O <- T_AwP.OHwC.O <- T_AwP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AwP.O <- cbind(T_AwP.OEu.O, T_AwP.OEvA.O, T_AwP.OEvP.O, T_AwP.OEvC.O, T_AwP.OEwA.O, T_AwP.OEwP.O, T_AwP.OEwC.O, T_AwP.OEbP.O,
                     T_AwP.OAu.O, T_AwP.OAvA.O, T_AwP.OAvP.O, T_AwP.OAvC.O, T_AwP.OAwA.O, T_AwP.OAwP.O, T_AwP.OAwC.O, T_AwP.OAbP.O,
                     T_AwP.OIu.O, T_AwP.OIvA.O, T_AwP.OIvP.O, T_AwP.OIvC.O, T_AwP.OIwA.O, T_AwP.OIwP.O, T_AwP.OIwC.O, T_AwP.OIbP.O,
                     T_AwP.OHu.O, T_AwP.OHvA.O, T_AwP.OHvP.O, T_AwP.OHvC.O, T_AwP.OHwA.O, T_AwP.OHwP.O, T_AwP.OHwC.O, T_AwP.OHbP.O)
    
    ## AwC
    # AwC.OxE cols
    T_AwC.OEwC.O <- diag((1-ihr_wC.O)*gamma.O*asymp_wC.O, ncol = age.bins, nrow = age.bins)
    T_AwC.OEu.O <- T_AwC.OEvA.O <- 
      T_AwC.OEvP.O <- T_AwC.OEwA.O <- T_AwC.OEwP.O <- T_AwC.OEvC.O <- T_AwC.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwC.OxA cols
    T_AwC.OAwC.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AwC.OAu.O <- T_AwC.OAvA.O <- T_AwC.OAvP.O <- 
      T_AwC.OAwA.O <- T_AwC.OAwP.O <- T_AwC.OAvC.O <- T_AwC.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwC.OxI cols
    T_AwC.OIu.O <- T_AwC.OIvA.O <- T_AwC.OIvP.O <- T_AwC.OIvC.O <- 
      T_AwC.OIwA.O <- T_AwC.OIwP.O <- T_AwC.OIwC.O <- T_AwC.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AwC.OxH cols
    T_AwC.OHu.O <- T_AwC.OHvA.O <- T_AwC.OHvP.O <- T_AwC.OHvC.O <- 
      T_AwC.OHwA.O <- T_AwC.OHwP.O <- T_AwC.OHwC.O <- T_AwC.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AwC.O <- cbind(T_AwC.OEu.O, T_AwC.OEvA.O, T_AwC.OEvP.O, T_AwC.OEvC.O, T_AwC.OEwA.O, T_AwC.OEwP.O, T_AwC.OEwC.O, T_AwC.OEbP.O,
                     T_AwC.OAu.O, T_AwC.OAvA.O, T_AwC.OAvP.O, T_AwC.OAvC.O, T_AwC.OAwA.O, T_AwC.OAwP.O, T_AwC.OAwC.O, T_AwC.OAbP.O,
                     T_AwC.OIu.O, T_AwC.OIvA.O, T_AwC.OIvP.O, T_AwC.OIvC.O, T_AwC.OIwA.O, T_AwC.OIwP.O, T_AwC.OIwC.O, T_AwC.OIbP.O,
                     T_AwC.OHu.O, T_AwC.OHvA.O, T_AwC.OHvP.O, T_AwC.OHvC.O, T_AwC.OHwA.O, T_AwC.OHwP.O, T_AwC.OHwC.O, T_AwC.OHbP.O)
    
    ## AbP
    # AbP.OxE cols
    T_AbP.OEbP.O <- diag((1-ihr_bP.O)*gamma.O*asymp_bP.O, ncol = age.bins, nrow = age.bins)
    T_AbP.OEu.O <- T_AbP.OEvA.O <- 
      T_AbP.OEvP.O <- T_AbP.OEwA.O <- T_AbP.OEwP.O <- T_AbP.OEvC.O <- T_AbP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AbP.OxA cols
    T_AbP.OAbP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_AbP.OAu.O <- T_AbP.OAvA.O <- T_AbP.OAvP.O <- 
      T_AbP.OAwA.O <- T_AbP.OAwP.O <- T_AbP.OAvC.O <- T_AbP.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AbP.OxI cols
    T_AbP.OIu.O <- T_AbP.OIvA.O <- T_AbP.OIvP.O <- T_AbP.OIvC.O <- 
      T_AbP.OIwA.O <- T_AbP.OIwP.O <- T_AbP.OIwC.O <- T_AbP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # AbP.OxH cols
    T_AbP.OHu.O <- T_AbP.OHvA.O <- T_AbP.OHvP.O <- T_AbP.OHvC.O <- 
      T_AbP.OHwA.O <- T_AbP.OHwP.O <- T_AbP.OHwC.O <- T_AbP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_AbP.O <- cbind(T_AbP.OEu.O, T_AbP.OEvA.O, T_AbP.OEvP.O, T_AbP.OEvC.O, T_AbP.OEwA.O, T_AbP.OEwP.O, T_AbP.OEwC.O, T_AbP.OEbP.O,
                     T_AbP.OAu.O, T_AbP.OAvA.O, T_AbP.OAvP.O, T_AbP.OAvC.O, T_AbP.OAwA.O, T_AbP.OAwP.O, T_AbP.OAwC.O, T_AbP.OAbP.O,
                     T_AbP.OIu.O, T_AbP.OIvA.O, T_AbP.OIvP.O, T_AbP.OIvC.O, T_AbP.OIwA.O, T_AbP.OIwP.O, T_AbP.OIwC.O, T_AbP.OIbP.O,
                     T_AbP.OHu.O, T_AbP.OHvA.O, T_AbP.OHvP.O, T_AbP.OHvC.O, T_AbP.OHwA.O, T_AbP.OHwP.O, T_AbP.OHwC.O, T_AbP.OHbP.O)
    #all A's together
    T_A.O <- rbind(T_Au.O, T_AvA.O, T_AvP.O, T_AvC.O, T_AwA.O, T_AwP.O, T_AwC.O, T_AbP.O)
    
    ## Iu
    # Iu.OxE cols
    T_Iu.OEu.O <- diag((1-ihr.O)*gamma.O*(1-asymp), ncol = age.bins, nrow = age.bins)
    T_Iu.OEvA.O <- T_Iu.OEvP.O <- 
      T_Iu.OEvC.O <- T_Iu.OEwA.O <- T_Iu.OEwP.O <- T_Iu.OEwC.O <- T_Iu.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Iu.OxI cols
    T_Iu.OIu.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_Iu.OIvA.O <- T_Iu.OIvP.O <- T_Iu.OIvC.O <- 
      T_Iu.OIwA.O <- T_Iu.OIwP.O <- T_Iu.OIwC.O <- T_Iu.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Iu.OxA cols
    T_Iu.OAu.O <- T_Iu.OAvA.O <- T_Iu.OAvP.O <- T_Iu.OAvC.O <- 
      T_Iu.OAwA.O <- T_Iu.OAwP.O <- T_Iu.OAwC.O <- T_Iu.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Iu.OxH cols
    T_Iu.OHu.O <- T_Iu.OHvA.O <- T_Iu.OHvP.O <- T_Iu.OHvC.O <- 
      T_Iu.OHwA.O <- T_Iu.OHwP.O <- T_Iu.OHwC.O <- T_Iu.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_Iu.O <- cbind(T_Iu.OEu.O, T_Iu.OEvA.O, T_Iu.OEvP.O, T_Iu.OEvC.O, T_Iu.OEwA.O, T_Iu.OEwP.O, T_Iu.OEwC.O, T_Iu.OEbP.O,
                    T_Iu.OAu.O, T_Iu.OAvA.O, T_Iu.OAvP.O, T_Iu.OAvC.O, T_Iu.OAwA.O, T_Iu.OAwP.O, T_Iu.OAwC.O, T_Iu.OAbP.O,
                    T_Iu.OIu.O, T_Iu.OIvA.O, T_Iu.OIvP.O, T_Iu.OIvC.O, T_Iu.OIwA.O, T_Iu.OIwP.O, T_Iu.OIwC.O, T_Iu.OIbP.O,
                    T_Iu.OHu.O, T_Iu.OHvA.O, T_Iu.OHvP.O, T_Iu.OHvC.O, T_Iu.OHwA.O, T_Iu.OHwP.O, T_Iu.OHwC.O, T_Iu.OHbP.O)
    
    ## IvA
    # IvA.OxE cols
    T_IvA.OEvA.O <- diag((1-ihr_vA.O)*gamma.O*(1-asymp_vA.O), ncol = age.bins, nrow = age.bins)
    T_IvA.OEu.O <- T_IvA.OEvP.O <- 
      T_IvA.OEvC.O <- T_IvA.OEwA.O <- T_IvA.OEwP.O <- T_IvA.OEwC.O <- T_IvA.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvA.OxI cols
    T_IvA.OIvA.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IvA.OIu.O <- T_IvA.OIvP.O <- T_IvA.OIvC.O <- 
      T_IvA.OIwA.O <- T_IvA.OIwP.O <- T_IvA.OIwC.O <- T_IvA.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvA.OxA cols
    T_IvA.OAu.O <- T_IvA.OAvA.O <- T_IvA.OAvP.O <- T_IvA.OAvC.O <- 
      T_IvA.OAwA.O <- T_IvA.OAwP.O <- T_IvA.OAwC.O <- T_IvA.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvA.OxH cols
    T_IvA.OHu.O <- T_IvA.OHvA.O <- T_IvA.OHvP.O <- T_IvA.OHvC.O <- 
      T_IvA.OHwA.O <- T_IvA.OHwP.O <- T_IvA.OHwC.O <- T_IvA.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IvA.O <- cbind(T_IvA.OEu.O, T_IvA.OEvA.O, T_IvA.OEvP.O, T_IvA.OEvC.O, T_IvA.OEwA.O, T_IvA.OEwP.O, T_IvA.OEwC.O, T_IvA.OEbP.O,
                     T_IvA.OAu.O, T_IvA.OAvA.O, T_IvA.OAvP.O, T_IvA.OAvC.O, T_IvA.OAwA.O, T_IvA.OAwP.O, T_IvA.OAwC.O, T_IvA.OAbP.O,
                     T_IvA.OIu.O, T_IvA.OIvA.O, T_IvA.OIvP.O, T_IvA.OIvC.O, T_IvA.OIwA.O, T_IvA.OIwP.O, T_IvA.OIwC.O, T_IvA.OIbP.O,
                     T_IvA.OHu.O, T_IvA.OHvA.O, T_IvA.OHvP.O, T_IvA.OHvC.O, T_IvA.OHwA.O, T_IvA.OHwP.O, T_IvA.OHwC.O, T_IvA.OHbP.O)
    
    ## IvP
    # IvP.OxE cols
    T_IvP.OEvP.O <- diag((1-ihr_vP.O)*gamma.O*(1-asymp_vP.O), ncol = age.bins, nrow = age.bins)
    T_IvP.OEu.O <- T_IvP.OEvA.O <- 
      T_IvP.OEvC.O <- T_IvP.OEwA.O <- T_IvP.OEwP.O <- T_IvP.OEwC.O <- T_IvP.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvP.OxI cols
    T_IvP.OIvP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IvP.OIu.O <- T_IvP.OIvA.O <- T_IvP.OIvC.O <- 
      T_IvP.OIwA.O <- T_IvP.OIwP.O <- T_IvP.OIwC.O <- T_IvP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvP.OxA cols
    T_IvP.OAu.O <- T_IvP.OAvA.O <- T_IvP.OAvP.O <- T_IvP.OAvC.O <- 
      T_IvP.OAwA.O <- T_IvP.OAwP.O <- T_IvP.OAwC.O <- T_IvP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvP.OxH cols
    T_IvP.OHu.O <- T_IvP.OHvA.O <- T_IvP.OHvP.O <- T_IvP.OHvC.O <- 
      T_IvP.OHwA.O <- T_IvP.OHwP.O <- T_IvP.OHwC.O <- T_IvP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IvP.O <- cbind(T_IvP.OEu.O, T_IvP.OEvA.O, T_IvP.OEvP.O, T_IvP.OEvC.O, T_IvP.OEwA.O, T_IvP.OEwP.O, T_IvP.OEwC.O, T_IvP.OEbP.O,
                     T_IvP.OAu.O, T_IvP.OAvA.O, T_IvP.OAvP.O, T_IvP.OAvC.O, T_IvP.OAwA.O, T_IvP.OAwP.O, T_IvP.OAwC.O, T_IvP.OAbP.O,
                     T_IvP.OIu.O, T_IvP.OIvA.O, T_IvP.OIvP.O, T_IvP.OIvC.O, T_IvP.OIwA.O, T_IvP.OIwP.O, T_IvP.OIwC.O, T_IvP.OIbP.O,
                     T_IvP.OHu.O, T_IvP.OHvA.O, T_IvP.OHvP.O, T_IvP.OHvC.O, T_IvP.OHwA.O, T_IvP.OHwP.O, T_IvP.OHwC.O, T_IvP.OHbP.O)
    
    ## IvC
    # IvC.OxE cols
    T_IvC.OEvC.O <- diag((1-ihr_vC.O)*gamma.O*(1-asymp_vC.O), ncol = age.bins, nrow = age.bins)
    T_IvC.OEu.O <- T_IvC.OEvA.O <- 
      T_IvC.OEvP.O <- T_IvC.OEwA.O <- T_IvC.OEwP.O <- T_IvC.OEwC.O <- T_IvC.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvC.OxI cols
    T_IvC.OIvC.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IvC.OIu.O <- T_IvC.OIvA.O <- T_IvC.OIvP.O <- 
      T_IvC.OIwA.O <- T_IvC.OIwP.O <- T_IvC.OIwC.O <- T_IvC.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvC.OxA cols
    T_IvC.OAu.O <- T_IvC.OAvA.O <- T_IvC.OAvP.O <- T_IvC.OAvC.O <- 
      T_IvC.OAwA.O <- T_IvC.OAwP.O <- T_IvC.OAwC.O <- T_IvC.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IvC.OxH cols
    T_IvC.OHu.O <- T_IvC.OHvA.O <- T_IvC.OHvP.O <- T_IvC.OHvC.O <- 
      T_IvC.OHwA.O <- T_IvC.OHwP.O <- T_IvC.OHwC.O <- T_IvC.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IvC.O <- cbind(T_IvC.OEu.O, T_IvC.OEvA.O, T_IvC.OEvP.O, T_IvC.OEvC.O, T_IvC.OEwA.O, T_IvC.OEwP.O, T_IvC.OEwC.O, T_IvC.OEbP.O,
                     T_IvC.OAu.O, T_IvC.OAvA.O, T_IvC.OAvP.O, T_IvC.OAvC.O, T_IvC.OAwA.O, T_IvC.OAwP.O, T_IvC.OAwC.O, T_IvC.OAbP.O,
                     T_IvC.OIu.O, T_IvC.OIvA.O, T_IvC.OIvP.O, T_IvC.OIvC.O, T_IvC.OIwA.O, T_IvC.OIwP.O, T_IvC.OIwC.O, T_IvC.OIbP.O,
                     T_IvC.OHu.O, T_IvC.OHvA.O, T_IvC.OHvP.O, T_IvC.OHvC.O, T_IvC.OHwA.O, T_IvC.OHwP.O, T_IvC.OHwC.O, T_IvC.OHbP.O)
    
    ## IwA
    # IwA.OxE cols
    T_IwA.OEwA.O <- diag((1-ihr_wA.O)*gamma.O*(1-asymp_wA.O), ncol = age.bins, nrow = age.bins)
    T_IwA.OEu.O <- T_IwA.OEvP.O <- 
      T_IwA.OEvC.O <- T_IwA.OEvA.O <- T_IwA.OEwP.O <- T_IwA.OEwC.O <- T_IwA.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwA.OxI cols
    T_IwA.OIwA.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IwA.OIu.O <- T_IwA.OIvP.O <- T_IwA.OIvC.O <- 
      T_IwA.OIvA.O <- T_IwA.OIwP.O <- T_IwA.OIwC.O <- T_IwA.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwA.OxA cols
    T_IwA.OAu.O <- T_IwA.OAvA.O <- T_IwA.OAvP.O <- T_IwA.OAvC.O <- 
      T_IwA.OAwA.O <- T_IwA.OAwP.O <- T_IwA.OAwC.O <- T_IwA.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwA.OxH cols
    T_IwA.OHu.O <- T_IwA.OHvA.O <- T_IwA.OHvP.O <- T_IwA.OHvC.O <- 
      T_IwA.OHwA.O <- T_IwA.OHwP.O <- T_IwA.OHwC.O <- T_IwA.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IwA.O <- cbind(T_IwA.OEu.O, T_IwA.OEvA.O, T_IwA.OEvP.O, T_IwA.OEvC.O, T_IwA.OEwA.O, T_IwA.OEwP.O, T_IwA.OEwC.O, T_IwA.OEbP.O,
                     T_IwA.OAu.O, T_IwA.OAvA.O, T_IwA.OAvP.O, T_IwA.OAvC.O, T_IwA.OAwA.O, T_IwA.OAwP.O, T_IwA.OAwC.O, T_IwA.OAbP.O,
                     T_IwA.OIu.O, T_IwA.OIvA.O, T_IwA.OIvP.O, T_IwA.OIvC.O, T_IwA.OIwA.O, T_IwA.OIwP.O, T_IwA.OIwC.O, T_IwA.OIbP.O,
                     T_IwA.OHu.O, T_IwA.OHvA.O, T_IwA.OHvP.O, T_IwA.OHvC.O, T_IwA.OHwA.O, T_IwA.OHwP.O, T_IwA.OHwC.O, T_IwA.OHbP.O)
    
    ## IwP
    # IwP.OxE cols
    T_IwP.OEwP.O <- diag((1-ihr_wP.O)*gamma.O*(1-asymp_wP.O), ncol = age.bins, nrow = age.bins)
    T_IwP.OEu.O <- T_IwP.OEvA.O <- 
      T_IwP.OEvC.O <- T_IwP.OEwA.O <- T_IwP.OEvP.O <- T_IwP.OEwC.O <- T_IwP.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwP.OxI cols
    T_IwP.OIwP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IwP.OIu.O <- T_IwP.OIvA.O <- T_IwP.OIvC.O <- 
      T_IwP.OIwA.O <- T_IwP.OIvP.O <- T_IwP.OIwC.O <- T_IwP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwP.OxA cols
    T_IwP.OAu.O <- T_IwP.OAvA.O <- T_IwP.OAvP.O <- T_IwP.OAvC.O <- 
      T_IwP.OAwA.O <- T_IwP.OAwP.O <- T_IwP.OAwC.O <- T_IwP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwP.OxH cols
    T_IwP.OHu.O <- T_IwP.OHvA.O <- T_IwP.OHvP.O <- T_IwP.OHvC.O <- 
      T_IwP.OHwA.O <- T_IwP.OHwP.O <- T_IwP.OHwC.O <- T_IwP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IwP.O <- cbind(T_IwP.OEu.O, T_IwP.OEvA.O, T_IwP.OEvP.O, T_IwP.OEvC.O, T_IwP.OEwA.O, T_IwP.OEwP.O, T_IwP.OEwC.O, T_IwP.OEbP.O,
                     T_IwP.OAu.O, T_IwP.OAvA.O, T_IwP.OAvP.O, T_IwP.OAvC.O, T_IwP.OAwA.O, T_IwP.OAwP.O, T_IwP.OAwC.O, T_IwP.OAbP.O,
                     T_IwP.OIu.O, T_IwP.OIvA.O, T_IwP.OIvP.O, T_IwP.OIvC.O, T_IwP.OIwA.O, T_IwP.OIwP.O, T_IwP.OIwC.O, T_IwP.OIbP.O,
                     T_IwP.OHu.O, T_IwP.OHvA.O, T_IwP.OHvP.O, T_IwP.OHvC.O, T_IwP.OHwA.O, T_IwP.OHwP.O, T_IwP.OHwC.O, T_IwP.OHbP.O)
    
    ## IwC
    # IwC.OxE cols
    T_IwC.OEwC.O <- diag((1-ihr_wC.O)*gamma.O*(1-asymp_wC.O), ncol = age.bins, nrow = age.bins)
    T_IwC.OEu.O <- T_IwC.OEvA.O <- 
      T_IwC.OEvP.O <- T_IwC.OEwA.O <- T_IwC.OEwP.O <- T_IwC.OEvC.O <- T_IwC.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.OxI cols
    T_IwC.OIwC.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IwC.OIu.O <- T_IwC.OIvA.O <- T_IwC.OIvP.O <- 
      T_IwC.OIwA.O <- T_IwC.OIwP.O <- T_IwC.OIvC.O <- T_IwC.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.OxA cols
    T_IwC.OAu.O <- T_IwC.OAvA.O <- T_IwC.OAvP.O <- T_IwC.OAvC.O <- 
      T_IwC.OAwA.O <- T_IwC.OAwP.O <- T_IwC.OAwC.O <- T_IwC.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.OxH cols
    T_IwC.OHu.O <- T_IwC.OHvA.O <- T_IwC.OHvP.O <- T_IwC.OHvC.O <- 
      T_IwC.OHwA.O <- T_IwC.OHwP.O <- T_IwC.OHwC.O <- T_IwC.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IwC.O <- cbind(T_IwC.OEu.O, T_IwC.OEvA.O, T_IwC.OEvP.O, T_IwC.OEvC.O, T_IwC.OEwA.O, T_IwC.OEwP.O, T_IwC.OEwC.O, T_IwC.OEbP.O,
                     T_IwC.OAu.O, T_IwC.OAvA.O, T_IwC.OAvP.O, T_IwC.OAvC.O, T_IwC.OAwA.O, T_IwC.OAwP.O, T_IwC.OAwC.O, T_IwC.OAbP.O,
                     T_IwC.OIu.O, T_IwC.OIvA.O, T_IwC.OIvP.O, T_IwC.OIvC.O, T_IwC.OIwA.O, T_IwC.OIwP.O, T_IwC.OIwC.O, T_IwC.OIbP.O,
                     T_IwC.OHu.O, T_IwC.OHvA.O, T_IwC.OHvP.O, T_IwC.OHvC.O, T_IwC.OHwA.O, T_IwC.OHwP.O, T_IwC.OHwC.O, T_IwC.OHbP.O)
    
    ## IwC
    # IwC.OxE cols
    T_IbP.OEbP.O <- diag((1-ihr_bP.O)*gamma.O*(1-asymp_bP.O), ncol = age.bins, nrow = age.bins)
    T_IbP.OEu.O <- T_IbP.OEvA.O <- 
      T_IbP.OEvP.O <- T_IbP.OEwA.O <- T_IbP.OEwP.O <- T_IbP.OEvC.O <- T_IbP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.OxI cols
    T_IbP.OIbP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
    T_IbP.OIu.O <- T_IbP.OIvA.O <- T_IbP.OIvP.O <- 
      T_IbP.OIwA.O <- T_IbP.OIwP.O <- T_IbP.OIvC.O <- T_IbP.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.OxA cols
    T_IbP.OAu.O <- T_IbP.OAvA.O <- T_IbP.OAvP.O <- T_IbP.OAvC.O <- 
      T_IbP.OAwA.O <- T_IbP.OAwP.O <- T_IbP.OAwC.O <- T_IbP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # IwC.OxH cols
    T_IbP.OHu.O <- T_IbP.OHvA.O <- T_IbP.OHvP.O <- T_IbP.OHvC.O <- 
      T_IbP.OHwA.O <- T_IbP.OHwP.O <- T_IbP.OHwC.O <- T_IbP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_IbP.O <- cbind(T_IbP.OEu.O, T_IbP.OEvA.O, T_IbP.OEvP.O, T_IbP.OEvC.O, T_IbP.OEwA.O, T_IbP.OEwP.O, T_IbP.OEwC.O, T_IbP.OEbP.O,
                     T_IbP.OAu.O, T_IbP.OAvA.O, T_IbP.OAvP.O, T_IbP.OAvC.O, T_IbP.OAwA.O, T_IbP.OAwP.O, T_IbP.OAwC.O, T_IbP.OAbP.O,
                     T_IbP.OIu.O, T_IbP.OIvA.O, T_IbP.OIvP.O, T_IbP.OIvC.O, T_IbP.OIwA.O, T_IbP.OIwP.O, T_IbP.OIwC.O, T_IbP.OIbP.O,
                     T_IbP.OHu.O, T_IbP.OHvA.O, T_IbP.OHvP.O, T_IbP.OHvC.O, T_IbP.OHwA.O, T_IbP.OHwP.O, T_IbP.OHwC.O, T_IbP.OHbP.O)   
    #all I's together
    T_I.O <- rbind(T_Iu.O, T_IvA.O, T_IvP.O, T_IvC.O, T_IwA.O, T_IwP.O, T_IwC.O, T_IbP.O)
    
    ## Hu
    # Hu.OxE cols
    T_Hu.OEu.O <- diag(ihr.O*gamma.O, ncol = age.bins, nrow = age.bins)
    T_Hu.OEvA.O <- T_Hu.OEvP.O <- 
      T_Hu.OEvC.O <- T_Hu.OEwA.O <- T_Hu.OEwP.O <- T_Hu.OEwC.O <- T_Hu.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Hu.OxH cols
    T_Hu.OHu.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_Hu.OHvA.O <- T_Hu.OHvP.O <- T_Hu.OHvC.O <- 
      T_Hu.OHwA.O <- T_Hu.OHwP.O <- T_Hu.OHwC.O <- T_Hu.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Hu.OxI cols
    T_Hu.OIu.O <- T_Hu.OIvA.O <- T_Hu.OIvP.O <- T_Hu.OIvC.O <- 
      T_Hu.OIwA.O <- T_Hu.OIwP.O <- T_Hu.OIwC.O <- T_Hu.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # Hu.OxA cols
    T_Hu.OAu.O <- T_Hu.OAvA.O <- T_Hu.OAvP.O <- T_Hu.OAvC.O <- 
      T_Hu.OAwA.O <- T_Hu.OAwP.O <- T_Hu.OAwC.O <- T_Hu.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_Hu.O <- cbind(T_Hu.OEu.O, T_Hu.OEvA.O, T_Hu.OEvP.O, T_Hu.OEvC.O, T_Hu.OEwA.O, T_Hu.OEwP.O, T_Hu.OEwC.O, T_Hu.OEbP.O,
                    T_Hu.OAu.O, T_Hu.OAvA.O, T_Hu.OAvP.O, T_Hu.OAvC.O, T_Hu.OAwA.O, T_Hu.OAwP.O, T_Hu.OAwC.O, T_Hu.OAbP.O,
                    T_Hu.OIu.O, T_Hu.OIvA.O, T_Hu.OIvP.O, T_Hu.OIvC.O, T_Hu.OIwA.O, T_Hu.OIwP.O, T_Hu.OIwC.O, T_Hu.OIbP.O,
                    T_Hu.OHu.O, T_Hu.OHvA.O, T_Hu.OHvP.O, T_Hu.OHvC.O, T_Hu.OHwA.O, T_Hu.OHwP.O, T_Hu.OHwC.O, T_Hu.OHbP.O)
    
    ## HvA
    # HvA.OxE cols
    T_HvA.OEvA.O <- diag(ihr_vA.O*gamma.O, ncol = age.bins, nrow = age.bins)
    T_HvA.OEu.O <- T_HvA.OEvP.O <- 
      T_HvA.OEvC.O <- T_HvA.OEwA.O <- T_HvA.OEwP.O <- T_HvA.OEwC.O <- T_HvA.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvA.OxH cols
    T_HvA.OHvA.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HvA.OHu.O <- T_HvA.OHvP.O <- T_HvA.OHvC.O <- 
      T_HvA.OHwA.O <- T_HvA.OHwP.O <- T_HvA.OHwC.O <- T_HvA.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvA.OxI cols
    T_HvA.OIu.O <- T_HvA.OIvA.O <- T_HvA.OIvP.O <- T_HvA.OIvC.O <- 
      T_HvA.OIwA.O <- T_HvA.OIwP.O <- T_HvA.OIwC.O <- T_HvA.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvA.OxA cols
    T_HvA.OAu.O <- T_HvA.OAvA.O <- T_HvA.OAvP.O <- T_HvA.OAvC.O <- 
      T_HvA.OAwA.O <- T_HvA.OAwP.O <- T_HvA.OAwC.O <- T_HvA.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HvA.O <- cbind(T_HvA.OEu.O, T_HvA.OEvA.O, T_HvA.OEvP.O, T_HvA.OEvC.O, T_HvA.OEwA.O, T_HvA.OEwP.O, T_HvA.OEwC.O, T_HvA.OEbP.O,
                     T_HvA.OAu.O, T_HvA.OAvA.O, T_HvA.OAvP.O, T_HvA.OAvC.O, T_HvA.OAwA.O, T_HvA.OAwP.O, T_HvA.OAwC.O, T_HvA.OAbP.O,
                     T_HvA.OIu.O, T_HvA.OIvA.O, T_HvA.OIvP.O, T_HvA.OIvC.O, T_HvA.OIwA.O, T_HvA.OIwP.O, T_HvA.OIwC.O, T_HvA.OIbP.O,
                     T_HvA.OHu.O, T_HvA.OHvA.O, T_HvA.OHvP.O, T_HvA.OHvC.O, T_HvA.OHwA.O, T_HvA.OHwP.O, T_HvA.OHwC.O, T_HvA.OHbP.O)
    
    ## HvP
    # HvP.OxE cols
    T_HvP.OEvP.O <- diag(ihr_vP.O*gamma.O, ncol = age.bins, nrow = age.bins)
    T_HvP.OEu.O <- T_HvP.OEvA.O <- 
      T_HvP.OEvC.O <- T_HvP.OEwA.O <- T_HvP.OEwP.O <- T_HvP.OEwC.O <- T_HvP.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvP.OxH cols
    T_HvP.OHvP.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HvP.OHu.O <- T_HvP.OHvA.O <- T_HvP.OHvC.O <- 
      T_HvP.OHwA.O <- T_HvP.OHwP.O <- T_HvP.OHwC.O <- T_HvP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvP.OxI cols
    T_HvP.OIu.O <- T_HvP.OIvA.O <- T_HvP.OIvP.O <- T_HvP.OIvC.O <- 
      T_HvP.OIwA.O <- T_HvP.OIwP.O <- T_HvP.OIwC.O <- T_HvP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvP.OxA cols
    T_HvP.OAu.O <- T_HvP.OAvA.O <- T_HvP.OAvP.O <- T_HvP.OAvC.O <- 
      T_HvP.OAwA.O <- T_HvP.OAwP.O <- T_HvP.OAwC.O <- T_HvP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HvP.O <- cbind(T_HvP.OEu.O, T_HvP.OEvA.O, T_HvP.OEvP.O, T_HvP.OEvC.O, T_HvP.OEwA.O, T_HvP.OEwP.O, T_HvP.OEwC.O, T_HvP.OEbP.O,
                     T_HvP.OAu.O, T_HvP.OAvA.O, T_HvP.OAvP.O, T_HvP.OAvC.O, T_HvP.OAwA.O, T_HvP.OAwP.O, T_HvP.OAwC.O, T_HvP.OAbP.O,
                     T_HvP.OIu.O, T_HvP.OIvA.O, T_HvP.OIvP.O, T_HvP.OIvC.O, T_HvP.OIwA.O, T_HvP.OIwP.O, T_HvP.OIwC.O, T_HvP.OIbP.O,
                     T_HvP.OHu.O, T_HvP.OHvA.O, T_HvP.OHvP.O, T_HvP.OHvC.O, T_HvP.OHwA.O, T_HvP.OHwP.O, T_HvP.OHwC.O, T_HvP.OHbP.O)
    
    ## HvC
    # HvC.OxE cols
    T_HvC.OEvC.O <- diag(ihr_vC.O*gamma.O, ncol = age.bins, nrow = age.bins)
    T_HvC.OEu.O <- T_HvC.OEvA.O <- 
      T_HvC.OEvP.O <- T_HvC.OEwA.O <- T_HvC.OEwP.O <- T_HvC.OEwC.O <- T_HvC.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvC.OxH cols
    T_HvC.OHvC.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HvC.OHu.O <- T_HvC.OHvA.O <- T_HvC.OHvP.O <- 
      T_HvC.OHwA.O <- T_HvC.OHwP.O <- T_HvC.OHwC.O <- T_HvC.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvC.OxI cols
    T_HvC.OIu.O <- T_HvC.OIvA.O <- T_HvC.OIvP.O <- T_HvC.OIvC.O <- 
      T_HvC.OIwA.O <- T_HvC.OIwP.O <- T_HvC.OIwC.O <- T_HvC.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HvC.OxA cols
    T_HvC.OAu.O <- T_HvC.OAvA.O <- T_HvC.OAvP.O <- T_HvC.OAvC.O <- 
      T_HvC.OAwA.O <- T_HvC.OAwP.O <- T_HvC.OAwC.O <- T_HvC.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HvC.O <- cbind(T_HvC.OEu.O, T_HvC.OEvA.O, T_HvC.OEvP.O, T_HvC.OEvC.O, T_HvC.OEwA.O, T_HvC.OEwP.O, T_HvC.OEwC.O, T_HvC.OEbP.O,
                     T_HvC.OAu.O, T_HvC.OAvA.O, T_HvC.OAvP.O, T_HvC.OAvC.O, T_HvC.OAwA.O, T_HvC.OAwP.O, T_HvC.OAwC.O, T_HvC.OAbP.O,
                     T_HvC.OIu.O, T_HvC.OIvA.O, T_HvC.OIvP.O, T_HvC.OIvC.O, T_HvC.OIwA.O, T_HvC.OIwP.O, T_HvC.OIwC.O, T_HvC.OIbP.O,
                     T_HvC.OHu.O, T_HvC.OHvA.O, T_HvC.OHvP.O, T_HvC.OHvC.O, T_HvC.OHwA.O, T_HvC.OHwP.O, T_HvC.OHwC.O, T_HvC.OHbP.O)
    
    ## HwA
    # HwA.OxE cols
    T_HwA.OEwA.O <- diag(ihr_wA.O*gamma.O, ncol = age.bins, nrow = age.bins)
    T_HwA.OEu.O <- T_HwA.OEvP.O <- 
      T_HwA.OEvC.O <- T_HwA.OEvA.O <- T_HwA.OEwP.O <- T_HwA.OEwC.O <- T_HwA.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwA.OxH cols
    T_HwA.OHwA.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HwA.OHu.O <- T_HwA.OHvP.O <- T_HwA.OHvC.O <- 
      T_HwA.OHvA.O <- T_HwA.OHwP.O <- T_HwA.OHwC.O <- T_HwA.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwA.OxI cols
    T_HwA.OIu.O <- T_HwA.OIvA.O <- T_HwA.OIvP.O <- T_HwA.OIvC.O <- 
      T_HwA.OIwA.O <- T_HwA.OIwP.O <- T_HwA.OIwC.O <- T_HwA.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwA.OxA cols
    T_HwA.OAu.O <- T_HwA.OAvA.O <- T_HwA.OAvP.O <- T_HwA.OAvC.O <- 
      T_HwA.OAwA.O <- T_HwA.OAwP.O <- T_HwA.OAwC.O <- T_HwA.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HwA.O <- cbind(T_HwA.OEu.O, T_HwA.OEvA.O, T_HwA.OEvP.O, T_HwA.OEvC.O, T_HwA.OEwA.O, T_HwA.OEwP.O, T_HwA.OEwC.O, T_HwA.OEbP.O,
                     T_HwA.OAu.O, T_HwA.OAvA.O, T_HwA.OAvP.O, T_HwA.OAvC.O, T_HwA.OAwA.O, T_HwA.OAwP.O, T_HwA.OAwC.O, T_HwA.OAbP.O,
                     T_HwA.OIu.O, T_HwA.OIvA.O, T_HwA.OIvP.O, T_HwA.OIvC.O, T_HwA.OIwA.O, T_HwA.OIwP.O, T_HwA.OIwC.O, T_HwA.OIbP.O,
                     T_HwA.OHu.O, T_HwA.OHvA.O, T_HwA.OHvP.O, T_HwA.OHvC.O, T_HwA.OHwA.O, T_HwA.OHwP.O, T_HwA.OHwC.O, T_HwA.OHbP.O)
    
    ## HwP
    # HwP.OxE cols
    T_HwP.OEwP.O <- diag(ihr_wP.O*gamma.O, ncol = age.bins, nrow = age.bins)
    T_HwP.OEu.O <- T_HwP.OEvA.O <- 
      T_HwP.OEvC.O <- T_HwP.OEwA.O <- T_HwP.OEvP.O <- T_HwP.OEwC.O <- T_HwP.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwP.OxH cols
    T_HwP.OHwP.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HwP.OHu.O <- T_HwP.OHvA.O <- T_HwP.OHvC.O <- 
      T_HwP.OHwA.O <- T_HwP.OHvP.O <- T_HwP.OHwC.O <- T_HwP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwP.OxI cols
    T_HwP.OIu.O <- T_HwP.OIvA.O <- T_HwP.OIvP.O <- T_HwP.OIvC.O <- 
      T_HwP.OIwA.O <- T_HwP.OIwP.O <- T_HwP.OIwC.O <- T_HwP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwP.OxA cols
    T_HwP.OAu.O <- T_HwP.OAvA.O <- T_HwP.OAvP.O <- T_HwP.OAvC.O <- 
      T_HwP.OAwA.O <- T_HwP.OAwP.O <- T_HwP.OAwC.O <- T_HwP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HwP.O <- cbind(T_HwP.OEu.O, T_HwP.OEvA.O, T_HwP.OEvP.O, T_HwP.OEvC.O, T_HwP.OEwA.O, T_HwP.OEwP.O, T_HwP.OEwC.O, T_HwP.OEbP.O,
                     T_HwP.OAu.O, T_HwP.OAvA.O, T_HwP.OAvP.O, T_HwP.OAvC.O, T_HwP.OAwA.O, T_HwP.OAwP.O, T_HwP.OAwC.O, T_HwP.OAbP.O,
                     T_HwP.OIu.O, T_HwP.OIvA.O, T_HwP.OIvP.O, T_HwP.OIvC.O, T_HwP.OIwA.O, T_HwP.OIwP.O, T_HwP.OIwC.O, T_HwP.OIbP.O,
                     T_HwP.OHu.O, T_HwP.OHvA.O, T_HwP.OHvP.O, T_HwP.OHvC.O, T_HwP.OHwA.O, T_HwP.OHwP.O, T_HwP.OHwC.O, T_HwP.OHbP.O)
    
    ## HwC
    # HwC.OxE cols
    T_HwC.OEwC.O <- diag(ihr_wC.O*gamma.O, ncol = age.bins, nrow = age.bins)
    T_HwC.OEu.O <- T_HwC.OEvA.O <- 
      T_HwC.OEvP.O <- T_HwC.OEwA.O <- T_HwC.OEwP.O <- T_HwC.OEvC.O <- T_HwC.OEbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.OxH cols
    T_HwC.OHwC.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HwC.OHu.O <- T_HwC.OHvA.O <- T_HwC.OHvP.O <- 
      T_HwC.OHwA.O <- T_HwC.OHwP.O <- T_HwC.OHvC.O <- T_HwC.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.OxI cols
    T_HwC.OIu.O <- T_HwC.OIvA.O <- T_HwC.OIvP.O <- T_HwC.OIvC.O <- 
      T_HwC.OIwA.O <- T_HwC.OIwP.O <- T_HwC.OIwC.O <- T_HwC.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.OxA cols
    T_HwC.OAu.O <- T_HwC.OAvA.O <- T_HwC.OAvP.O <- T_HwC.OAvC.O <- 
      T_HwC.OAwA.O <- T_HwC.OAwP.O <- T_HwC.OAwC.O <- T_HwC.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HwC.O <- cbind(T_HwC.OEu.O, T_HwC.OEvA.O, T_HwC.OEvP.O, T_HwC.OEvC.O, T_HwC.OEwA.O, T_HwC.OEwP.O, T_HwC.OEwC.O, T_HwC.OEbP.O,
                     T_HwC.OAu.O, T_HwC.OAvA.O, T_HwC.OAvP.O, T_HwC.OAvC.O, T_HwC.OAwA.O, T_HwC.OAwP.O, T_HwC.OAwC.O, T_HwC.OAbP.O,
                     T_HwC.OIu.O, T_HwC.OIvA.O, T_HwC.OIvP.O, T_HwC.OIvC.O, T_HwC.OIwA.O, T_HwC.OIwP.O, T_HwC.OIwC.O, T_HwC.OIbP.O,
                     T_HwC.OHu.O, T_HwC.OHvA.O, T_HwC.OHvP.O, T_HwC.OHvC.O, T_HwC.OHwA.O, T_HwC.OHwP.O, T_HwC.OHwC.O, T_HwC.OHbP.O)
    
    ## HwC
    # HwC.OxE cols
    T_HbP.OEbP.O <- diag(ihr_bP.O*gamma.O, ncol = age.bins, nrow = age.bins)
    T_HbP.OEu.O <- T_HbP.OEvA.O <- 
      T_HbP.OEvP.O <- T_HbP.OEwA.O <- T_HbP.OEwP.O <- T_HbP.OEvC.O <- T_HbP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.OxH cols
    T_HbP.OHwC.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
    T_HbP.OHu.O <- T_HbP.OHvA.O <- T_HbP.OHvP.O <- 
      T_HbP.OHwA.O <- T_HbP.OHwP.O <- T_HbP.OHvC.O <- T_HbP.OHbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.OxI cols
    T_HbP.OIu.O <- T_HbP.OIvA.O <- T_HbP.OIvP.O <- T_HbP.OIvC.O <- 
      T_HbP.OIwA.O <- T_HbP.OIwP.O <- T_HbP.OIwC.O <- T_HbP.OIbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    # HwC.OxA cols
    T_HbP.OAu.O <- T_HbP.OAvA.O <- T_HbP.OAvP.O <- T_HbP.OAvC.O <- 
      T_HbP.OAwA.O <- T_HbP.OAwP.O <- T_HbP.OAwC.O <- T_HbP.OAbP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
    
    T_HbP.O <- cbind(T_HbP.OEu.O, T_HbP.OEvA.O, T_HbP.OEvP.O, T_HbP.OEvC.O, T_HbP.OEwA.O, T_HbP.OEwP.O, T_HbP.OEwC.O, T_HbP.OEbP.O,
                     T_HbP.OAu.O, T_HbP.OAvA.O, T_HbP.OAvP.O, T_HbP.OAvC.O, T_HbP.OAwA.O, T_HbP.OAwP.O, T_HbP.OAwC.O, T_HbP.OAbP.O,
                     T_HbP.OIu.O, T_HbP.OIvA.O, T_HbP.OIvP.O, T_HbP.OIvC.O, T_HbP.OIwA.O, T_HbP.OIwP.O, T_HbP.OIwC.O, T_HbP.OIbP.O,
                     T_HbP.OHu.O, T_HbP.OHvA.O, T_HbP.OHvP.O, T_HbP.OHvC.O, T_HbP.OHwA.O, T_HbP.OHwP.O, T_HbP.OHwC.O, T_HbP.OHbP.O)
    
    #all H's together
    T_H.O <- rbind(T_Hu.O, T_HvA.O, T_HvP.O, T_HvC.O, T_HwA.O, T_HwP.O, T_HwC.O, T_HbP.O)
    
    
    # final matrix
    T_mat.O <- rbind(T_E.O, T_A.O, T_I.O, T_H.O)
    # identity matrix
    Id_mat <- diag(1,nrow(T_mat.O),ncol(T_mat.O))
    # sparse because why not
    F_mat.O <- Matrix(F_mat.O, sparse = T)
    # automatically computes (I - T)^-1 and converts to sparse -- faster
    T_mat.O <- Matrix(inv(Id_mat - T_mat.O), sparse = T)
    # R0 is the maximum eigenvalue (that might be complex)
    R0.O <- max(Re(eigs(F_mat.O%*%T_mat.O, 1, retvec = F)$values))
    R0 <- list(R0.D = R0.D, R0.O = R0.O)
  })
  return(R0)
}
