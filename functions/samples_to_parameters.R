library(stringr)

samples_to_parameters <- function(samples,parameters){
  for(j in names(samples)){
    if(j == "gamma.O"){
      parameters$gamma.O <- 1-exp(-1/as.numeric(samples[j]))
    }
    if(j == "R0.D"){
      rel.value <- as.numeric(samples[j])/parameters$R0.D
      parameters$beta.D <- rel.value*parameters$beta.D
      parameters$beta_vA.D <- rel.value*parameters$beta_vA.D
      parameters$beta_wA.D <- rel.value*parameters$beta_wA.D
      parameters$beta_vP.D <- rel.value*parameters$beta_vP.D
      parameters$beta_wP.D <- rel.value*parameters$beta_wP.D
      parameters$beta_bP.D <- rel.value*parameters$beta_bP.D
      parameters$beta_vC.D <- rel.value*parameters$beta_vC.D
      parameters$beta_wC.D <- rel.value*parameters$beta_wC.D
      parameters$R0.D <- as.numeric(samples[j])
    } else if(j == "R0.O"){
      rel.value <- as.numeric(samples[j])/parameters$R0.O
      parameters$beta.O <- rel.value*parameters$beta.O
      parameters$beta_vA.O <- rel.value*parameters$beta_vA.O
      parameters$beta_wA.O <- rel.value*parameters$beta_wA.O
      parameters$beta_vP.O <- rel.value*parameters$beta_vP.O
      parameters$beta_wP.O <- rel.value*parameters$beta_wP.O
      parameters$beta_bP.O <- rel.value*parameters$beta_bP.O
      parameters$beta_vC.O <- rel.value*parameters$beta_vC.O
      parameters$beta_wC.O <- rel.value*parameters$beta_wC.O
      parameters$R0.O <- as.numeric(samples[j])
    }
    
    if(str_detect(j,"vax\\d\\.(symp|beta|hosp|death)\\.effic")){###se é eficacia
      ex <- gregexpr("[0-9]+",j)
      vals <- as.numeric(regmatches(j,ex)[[1]])
      if(length(vals) == 1){###aplica a todas as idades
        parameters[[as.character(j)]] <- as.numeric(samples[j])*rep(1,parameters$age.bins)
      }
      else if(length(vals) == 2){
        parameters[[str_extract(j,"vax\\d\\.(symp|beta|hosp|death)\\.effic\\.(P|A|C)\\.(D|O)")]][vals[2]] <- as.numeric(samples[j])
      }
      else if(length(vals) == 3){
        parameters[[str_extract(j,"vax\\d\\.(symp|beta|hosp|death)\\.effic\\.(P|A|C)\\.(D|O)")]][vals[2]:vals[3]] <- as.numeric(samples[j])
      }
    }
    else{
      parameters[[as.character(j)]] <- as.numeric(samples[j])
    }
  }
  if(any(str_detect(names(samples),"vax1\\.(symp|beta|hosp|death)\\.effic\\.P\\.D"))){
    ####recomputa eficácias#####
    VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                                 vax1.beta.effic.P.D,vax1.hosp.effic.P.D,
                                                                 vax1.death.effic.P.D,vax1.symp.effic.P.D))
    VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
    VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
    VAX1.EFFIC.DEATH  <- VAX1.PARS$death
    VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
    ##### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_vP.D <- (1- VAX1.EFFIC.BETA)*parameters$beta.D
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_vP.D <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
    parameters$ihr_vP.D <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.D
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_vP.D <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr
    parameters$vax.window.days.P <- VAX.WINDOW.DAYS
  }
  if(any(str_detect(names(samples),"vax2\\.(symp|beta|hosp|death)\\.effic\\.P.D"))){
    ####recomputa eficácias#####
    VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                                 vax2.beta.effic.P.D,vax2.hosp.effic.P.D,
                                                                 vax2.death.effic.P.D,vax2.symp.effic.P.D))
    VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
    VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
    VAX2.EFFIC.DEATH  <- VAX2.PARS$death
    VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
    #### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_wP.D <- (1- VAX2.EFFIC.BETA)*parameters$beta.D
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_wP.D <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (ihr.D) \sigma_v
    parameters$ihr_wP.D <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.D
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_wP.D <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax3\\.(symp|beta|hosp|death)\\.effic\\.P.D"))){
    VAX3.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                                 vax3.beta.effic.P.D,vax3.hosp.effic.P.D,
                                                                 vax3.death.effic.P.D,vax3.symp.effic.P.D))
    VAX3.EFFIC.BETA   <- VAX3.PARS$suscep
    VAX3.EFFIC.SEVERE <- VAX3.PARS$hosp
    VAX3.EFFIC.DEATH  <- VAX3.PARS$death
    VAX3.EFFIC.CLIN   <- VAX3.PARS$symp
    #### Classification of cases related to the severity of cases of once-vaccinated individuals  
    parameters$beta_bP.D <- (1- VAX3.EFFIC.BETA)*parameters$beta.D
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_bP.D <- 1 - (1 - parameters$asymp) * (1-VAX3.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (ihr.D) \sigma_v
    parameters$ihr_bP.D <- (1-VAX3.EFFIC.SEVERE)*parameters$ihr.D
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_bP.D <- (1-VAX3.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax1\\.(symp|beta|hosp|death)\\.effic\\.P\\.O"))){
    ####recomputa eficácias#####
    VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                                 vax1.beta.effic.P.O,vax1.hosp.effic.P.O,
                                                                 vax1.death.effic.P.O,vax1.symp.effic.P.O))
    VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
    VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
    VAX1.EFFIC.DEATH  <- VAX1.PARS$death
    VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
    ##### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_vP.O <- (1- VAX1.EFFIC.BETA)*parameters$beta.O
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_vP.O <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
    parameters$ihr_vP.O <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.O
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_vP.O <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax2\\.(symp|beta|hosp|death)\\.effic\\.P.O"))){
    ####recomputa eficácias#####
    VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                                 vax2.beta.effic.P.O,vax2.hosp.effic.P.O,
                                                                 vax2.death.effic.P.O,vax2.symp.effic.P.O))
    VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
    VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
    VAX2.EFFIC.DEATH  <- VAX2.PARS$death
    VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
    #### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_wP.O <- (1- VAX2.EFFIC.BETA)*parameters$beta.O
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_wP.O <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (ihr.O) \sigma_v
    parameters$ihr_wP.O <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.O
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_wP.O <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax3\\.(symp|beta|hosp|death)\\.effic\\.P.O"))){
    VAX3.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                                 vax3.beta.effic.P.O,vax3.hosp.effic.P.O,
                                                                 vax3.death.effic.P.O,vax3.symp.effic.P.O))
    VAX3.EFFIC.BETA   <- VAX3.PARS$suscep
    VAX3.EFFIC.SEVERE <- VAX3.PARS$hosp
    VAX3.EFFIC.DEATH  <- VAX3.PARS$death
    VAX3.EFFIC.CLIN   <- VAX3.PARS$symp
    #### Classification of cases related to the severity of cases of once-vaccinated individuals  
    parameters$beta_bP.O <- (1- VAX3.EFFIC.BETA)*parameters$beta.O
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_bP.O <- 1 - (1 - parameters$asymp) * (1-VAX3.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (ihr.O) \sigma_v
    parameters$ihr_bP.O <- (1-VAX3.EFFIC.SEVERE)*parameters$ihr.O
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_bP.O <- (1-VAX3.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax1\\.(symp|beta|hosp|death)\\.effic\\.A\\.D"))){
    ####recomputa eficácias#####
    VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                                 vax1.beta.effic.A.D,vax1.hosp.effic.A.D,
                                                                 vax1.death.effic.A.D,vax1.symp.effic.A.D))
    VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
    VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
    VAX1.EFFIC.DEATH  <- VAX1.PARS$death
    VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
    ##### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_vA.D <- (1- VAX1.EFFIC.BETA)*parameters$beta.D
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_vA.D <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
    parameters$ihr_vA.D <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.D
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_vA.D <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr
    parameters$vax.window.days.P <- VAX.WINDOW.DAYS
  }
  if(any(str_detect(names(samples),"vax2\\.(symp|beta|hosp|death)\\.effic\\.A.D"))){
    ####recomputa eficácias#####
    VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                                 vax2.beta.effic.A.D,vax2.hosp.effic.A.D,
                                                                 vax2.death.effic.A.D,vax2.symp.effic.A.D))
    VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
    VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
    VAX2.EFFIC.DEATH  <- VAX2.PARS$death
    VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
    #### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_wA.D <- (1- VAX2.EFFIC.BETA)*parameters$beta.D
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_wA.D <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (ihr.D) \sigma_v
    parameters$ihr_wA.D <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.D
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_wA.D <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax1\\.(symp|beta|hosp|death)\\.effic\\.A\\.O"))){
    ####recomputa eficácias#####
    VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                                 vax1.beta.effic.A.O,vax1.hosp.effic.A.O,
                                                                 vax1.death.effic.A.O,vax1.symp.effic.A.O))

    VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
    VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
    VAX1.EFFIC.DEATH  <- VAX1.PARS$death
    VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
    ##### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_vA.O <- (1- VAX1.EFFIC.BETA)*parameters$beta.O
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_vA.O <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
    parameters$ihr_vA.O <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.O
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_vA.O <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax2\\.(symp|beta|hosp|death)\\.effic\\.A.O"))){
    ####recomputa eficácias#####
    VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                                 vax2.beta.effic.A.O,vax2.hosp.effic.A.O,
                                                                 vax2.death.effic.A.O,vax2.symp.effic.A.O))
    VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
    VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
    VAX2.EFFIC.DEATH  <- VAX2.PARS$death
    VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
    #### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_wA.O <- (1- VAX2.EFFIC.BETA)*parameters$beta.O
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_wA.O <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (ihr.O) \sigma_v
    parameters$ihr_wA.O <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.O
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_wA.O <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax1\\.(symp|beta|hosp|death)\\.effic\\.C\\.D"))){
    ####recomputa eficácias#####
    VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                                 vax1.beta.effic.C.D,vax1.hosp.effic.C.D,
                                                                 vax1.death.effic.C.D,vax1.symp.effic.C.D))
    VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
    VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
    VAX1.EFFIC.DEATH  <- VAX1.PARS$death
    VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
    ##### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_vC.D <- (1- VAX1.EFFIC.BETA)*parameters$beta.D
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_vC.D <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
    parameters$ihr_vC.D <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.D
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_vC.D <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr
    
  }
  if(any(str_detect(names(samples),"vax2\\.(symp|beta|hosp|death)\\.effic\\.C.D"))){
    ####recomputa eficácias#####
    VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                                 vax2.beta.effic.C.D,vax2.hosp.effic.C.D,
                                                                 vax2.death.effic.C.D,vax2.symp.effic.C.D))
    VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
    VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
    VAX2.EFFIC.DEATH  <- VAX2.PARS$death
    VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
    #### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_wC.D <- (1- VAX2.EFFIC.BETA)*parameters$beta.D
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_wC.D <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (ihr.D) \sigma_v
    parameters$ihr_wC.D <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.D
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_wC.D <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax1\\.(symp|beta|hosp|death)\\.effic\\.C\\.O"))){
    ####recomputa eficácias#####
    VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                                 vax1.beta.effic.C.O,vax1.hosp.effic.C.O,
                                                                 vax1.death.effic.C.O,vax1.symp.effic.C.O))
    
    VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
    VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
    VAX1.EFFIC.DEATH  <- VAX1.PARS$death
    VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
    ##### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_vC.O <- (1- VAX1.EFFIC.BETA)*parameters$beta.O
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_vC.O <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
    parameters$ihr_vC.O <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.O
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_vC.O <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr
  }
  if(any(str_detect(names(samples),"vax2\\.(symp|beta|hosp|death)\\.effic\\.C.O"))){
    ####recomputa eficácias#####
    VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                                 vax2.beta.effic.C.O,vax2.hosp.effic.C.O,
                                                                 vax2.death.effic.C.O,vax2.symp.effic.C.O))
    VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
    VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
    VAX2.EFFIC.DEATH  <- VAX2.PARS$death
    VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
    #### Classification of cases related to the severity of cases of once-vaccinated individuals
    parameters$beta_wC.O <- (1- VAX2.EFFIC.BETA)*parameters$beta.O
    # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
    parameters$asymp_wC.O <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
    # Fraction of severe cases in symptomatic cases (ihr.O) \sigma_v
    parameters$ihr_wC.O <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.O
    # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
    parameters$ihfr_wC.O <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr
  }
   if(any(names(samples) == "prevalence")){
    ####reorganiza histórico
    PREV <- parameters$prevalence
    # history <- parameters$history.P[,1:10] + parameters$history.P[,11:20]
    # parameters$history.P <- cbind((1-PREV)*history,PREV*history)
    # history <- parameters$history.A[,1:10] + parameters$history.A[,11:20]
    # parameters$history.A <- cbind((1-PREV)*history,PREV*history)
    # history <- parameters$history.C[,1:10] + parameters$history.C[,11:20]
    # parameters$history.C <- cbind((1-PREV)*history,PREV*history)
    #######reorganiza condição inicial
    Y2 <- parameters$init.condition
    parameters$init.condition <- with(parameters,{
      Y <- Y2
      Y[Suindex] <- (1-PREV)*(Y2[Suindex]+Y2[Ru.Dindex])
      Y[Ru.Dindex] <- PREV*(Y2[Suindex]+Y2[Ru.Dindex])
      Y[SvPindex] <- (1-PREV)*(Y2[SvPindex]+Y2[RvP.Dindex])
      Y[RvP.Dindex] <- PREV*(Y2[SvPindex]+Y2[RvP.Dindex])
      Y[SwPindex] <- (1-PREV)*(Y2[SwPindex]+Y2[RwP.Dindex])
      Y[RwP.Dindex] <- PREV*(Y2[SwPindex]+Y2[RwP.Dindex])
      Y[SbPindex] <- (1-PREV)*(Y2[SbPindex]+Y2[RbP.Dindex])
      Y[RbP.Dindex] <- PREV*(Y2[SbPindex]+Y2[RbP.Dindex])
      Y[SvAindex] <- (1-PREV)*(Y2[SvAindex]+Y2[RvA.Dindex])
      Y[RvA.Dindex] <- PREV*(Y2[SvAindex]+Y2[RvA.Dindex])
      Y[SwAindex] <- (1-PREV)*(Y2[SwAindex]+Y2[RwA.Dindex])
      Y[RwA.Dindex] <- PREV*(Y2[SwAindex]+Y2[RwA.Dindex])
      Y[SvCindex] <- (1-PREV)*(Y2[SvCindex]+Y2[RvC.Dindex])
      Y[RvC.Dindex] <- PREV*(Y2[SvCindex]+Y2[RvC.Dindex])
      Y[SwCindex] <- (1-PREV)*(Y2[SwCindex]+Y2[RwC.Dindex])
      Y[RwC.Dindex] <- PREV*(Y2[SwCindex]+Y2[RwC.Dindex])
      Y
    })
    
  }
  
  return(parameters)
}
