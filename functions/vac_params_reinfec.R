source("functions/epi_params_omicron.R")
source("functions/opt_vax_rate_history.R")
source("functions/VAX_DISTR_RATE.R")
library(dplyr)

#### convert efficacy from studies to parameters
OBSERVED.EFFICACY.TO.PARAMETERS <- function(HOSP.PROP,MORT.PROP,ASYMP.PROP,SUSCEP.EFF,HOSP.EFF,MORT.EFF,SYMP.EFF){
  SUSCEP.EFF.PAR <- SUSCEP.EFF
  HOSP.EFF.PAR <- 1 - (1-HOSP.EFF)/(1-SUSCEP.EFF)
  MORT.EFF.PAR <- 1 - (1-MORT.EFF)/(1-HOSP.EFF)
  SYMP.EFF.PAR <- 1 - ((1-SYMP.EFF)*(HOSP.PROP + (1 - HOSP.PROP)*(1-ASYMP.PROP))  - (1 - HOSP.EFF)*HOSP.PROP)/
    ((1-SUSCEP.EFF)*(1 - (1-HOSP.EFF)/(1-SUSCEP.EFF)*HOSP.PROP)*(1-ASYMP.PROP))
  return(list(suscep = SUSCEP.EFF.PAR,hosp = HOSP.EFF.PAR,death = MORT.EFF.PAR,symp = SYMP.EFF.PAR))
}

MAX.TIME.DAYS = 90
################################################################################
############################### Pfizer #########################################
VAX.WINDOW.DAYS = 56
parameters$vax.window.days.P <- VAX.WINDOW.DAYS
parameters$vac.rate.P <- rep(0,MAX.TIME.DAYS)
parameters$vac.rate.v2.P <- rep(0,MAX.TIME.DAYS)
prev_hist <- matrix(0,nrow = MAX.TIME.DAYS, ncol = 3*age.bins)
parameters$history.P <- prev_hist
########################## DELTA VARIANT #######################################
#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
parameters$vax2.beta.effic.P.D = 0.90*rep(1,age.bins)
parameters$vax2.symp.effic.P.D = 0.94*rep(1,age.bins)
parameters$vax2.hosp.effic.P.D = 0.87*rep(1,age.bins)
parameters$vax2.death.effic.P.D = 0.98*rep(1,age.bins) ## GUESS

parameters$vax3.beta.effic.P.D = 0.90*rep(1,age.bins)
parameters$vax3.symp.effic.P.D = 0.94*rep(1,age.bins)
parameters$vax3.hosp.effic.P.D = 0.95*rep(1,age.bins)
parameters$vax3.death.effic.P.D = 0.99*rep(1,age.bins)

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.8
parameters$vax1.beta.effic.P.D   = FIRST.DOSE.REL.EFFIC * parameters$vax2.beta.effic.P.D
parameters$vax1.symp.effic.P.D   = FIRST.DOSE.REL.EFFIC * parameters$vax2.symp.effic.P.D
parameters$vax1.hosp.effic.P.D   = FIRST.DOSE.REL.EFFIC * parameters$vax2.hosp.effic.P.D
parameters$vax1.death.effic.P.D  = FIRST.DOSE.REL.EFFIC * parameters$vax2.death.effic.P.D

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                             vax1.beta.effic.P.D,vax1.hosp.effic.P.D,
                                                             vax1.death.effic.P.D,vax1.symp.effic.P.D))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                             vax2.beta.effic.P.D,vax2.hosp.effic.P.D,
                                                             vax2.death.effic.P.D,vax2.symp.effic.P.D))
VAX3.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                             vax3.beta.effic.P.D,vax3.hosp.effic.P.D,
                                                             vax3.death.effic.P.D,vax3.symp.effic.P.D))
VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
VAX1.EFFIC.DEATH  <- VAX1.PARS$death
VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
VAX2.EFFIC.DEATH  <- VAX2.PARS$death
VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
VAX3.EFFIC.BETA   <- VAX3.PARS$suscep
VAX3.EFFIC.SEVERE <- VAX3.PARS$hosp
VAX3.EFFIC.DEATH  <- VAX3.PARS$death
VAX3.EFFIC.CLIN   <- VAX3.PARS$symp

##############################
### FIRST DOSE INFORMATION ###
##############################
##### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_vP.D <- (1- VAX1.EFFIC.BETA)*parameters$beta.D
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_vP.D <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_vP.D <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.D
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_vP.D <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr
##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_wP.D <- (1- VAX2.EFFIC.BETA)*parameters$beta.D
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_wP.D <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (ihr.D) \sigma_v
parameters$ihr_wP.D <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.D
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_wP.D <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr
##############################
#### BOOSTER  INFORMATION ####
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals  
parameters$beta_bP.D <- (1- VAX3.EFFIC.BETA)*parameters$beta.D
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_bP.D <- 1 - (1 - parameters$asymp) * (1-VAX3.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (ihr.D) \sigma_v
parameters$ihr_bP.D <- (1-VAX3.EFFIC.SEVERE)*parameters$ihr.D
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_bP.D <- (1-VAX3.EFFIC.DEATH)*parameters$ihfr
################################################################################
########################## OMICRON VARIANT #####################################
#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
######ASSUMINDO UMA REDUÇÃO ARBITRARIA DE 50% NA EFETIVIDADE, MUDAR NO FUTURO
red <- 0.5
parameters$vax2.beta.effic.P.O = red*0.90*rep(1,age.bins)
parameters$vax2.symp.effic.P.O = red*0.94*rep(1,age.bins)
parameters$vax2.hosp.effic.P.O = red*0.87*rep(1,age.bins)
parameters$vax2.death.effic.P.O = red*0.98*rep(1,age.bins) ## GUESS

parameters$vax3.beta.effic.P.O = 0.90*rep(1,age.bins)
parameters$vax3.symp.effic.P.O = 0.94*rep(1,age.bins)
parameters$vax3.hosp.effic.P.O = 0.95*rep(1,age.bins)
parameters$vax3.death.effic.P.O = 0.99*rep(1,age.bins)

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.8
parameters$vax1.beta.effic.P.O   = FIRST.DOSE.REL.EFFIC * parameters$vax2.beta.effic.P.O
parameters$vax1.symp.effic.P.O   = FIRST.DOSE.REL.EFFIC * parameters$vax2.symp.effic.P.O
parameters$vax1.hosp.effic.P.O   = FIRST.DOSE.REL.EFFIC * parameters$vax2.hosp.effic.P.O
parameters$vax1.death.effic.P.O  = FIRST.DOSE.REL.EFFIC * parameters$vax2.death.effic.P.O

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                             vax1.beta.effic.P.O,vax1.hosp.effic.P.O,
                                                             vax1.death.effic.P.O,vax1.symp.effic.P.O))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                             vax2.beta.effic.P.O,vax2.hosp.effic.P.O,
                                                             vax2.death.effic.P.O,vax2.symp.effic.P.O))
VAX3.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                             vax3.beta.effic.P.O,vax3.hosp.effic.P.O,
                                                             vax3.death.effic.P.O,vax3.symp.effic.P.O))
VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
VAX1.EFFIC.DEATH  <- VAX1.PARS$death
VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
VAX2.EFFIC.DEATH  <- VAX2.PARS$death
VAX2.EFFIC.CLIN   <- VAX2.PARS$symp
VAX3.EFFIC.BETA   <- VAX3.PARS$suscep
VAX3.EFFIC.SEVERE <- VAX3.PARS$hosp
VAX3.EFFIC.DEATH  <- VAX3.PARS$death
VAX3.EFFIC.CLIN   <- VAX3.PARS$symp

##############################
### FIRST DOSE INFORMATION ###
##############################
##### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_vP.O <- (1- VAX1.EFFIC.BETA)*parameters$beta.O
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_vP.O <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_vP.O <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.O
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_vP.O <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr
##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_wP.O <- (1- VAX2.EFFIC.BETA)*parameters$beta.O
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_wP.O <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_wP.O <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.O
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_wP.O <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr
##############################
## BOOSTER DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_bP.O <- (1- VAX3.EFFIC.BETA)*parameters$beta.O
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_bP.O <- 1 - (1 - parameters$asymp) * (1-VAX3.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_bP.O <- (1-VAX3.EFFIC.SEVERE)*parameters$ihr.O
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_bP.O <- (1-VAX3.EFFIC.DEATH)*parameters$ihfr
################################################################################
################################################################################
############################ AstraZeneca #######################################
# # Time window between first and second vaccines
VAX.WINDOW.DAYS = 56
parameters$vax.window.days.A <- VAX.WINDOW.DAYS
parameters$vac.rate.A <- rep(0,MAX.TIME.DAYS)
parameters$vac.rate.v2.A <- rep(0,MAX.TIME.DAYS)
############################### DELTA VARIANT ##################################
#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
parameters$vax2.beta.effic.A.D = 0.6*rep(1,age.bins)
parameters$vax2.symp.effic.A.D = 0.81*rep(1,age.bins)
parameters$vax2.hosp.effic.A.D = 0.90*rep(1,age.bins)
parameters$vax2.death.effic.A.D = 0.95*rep(1,age.bins) ## GUESS

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.4
parameters$vax1.beta.effic.A.D   = FIRST.DOSE.REL.EFFIC * parameters$vax2.beta.effic.A.D
parameters$vax1.symp.effic.A.D   = FIRST.DOSE.REL.EFFIC * parameters$vax2.symp.effic.A.D
parameters$vax1.hosp.effic.A.D   = FIRST.DOSE.REL.EFFIC * parameters$vax2.hosp.effic.A.D
parameters$vax1.death.effic.A.D  = FIRST.DOSE.REL.EFFIC * parameters$vax2.death.effic.A.D

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                             vax1.beta.effic.A.D,vax1.hosp.effic.A.D,
                                                             vax1.death.effic.A.D,vax1.symp.effic.A.D))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                             vax2.beta.effic.A.D,vax2.hosp.effic.A.D,
                                                             vax2.death.effic.A.D,vax2.symp.effic.A.D))
VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
VAX1.EFFIC.DEATH  <- VAX1.PARS$death
VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
VAX2.EFFIC.DEATH  <- VAX2.PARS$death
VAX2.EFFIC.CLIN   <- VAX2.PARS$symp

#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_vA.D <- (1- VAX1.EFFIC.BETA)*parameters$beta.D
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_vA.D <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (ihr.D) \sigma_v
parameters$ihr_vA.D <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.D
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_vA.D <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr

##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_wA.D <- (1- VAX2.EFFIC.BETA)*parameters$beta.D
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_wA.D <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_wA.D <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.D
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_wA.D <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr

############################## OMICRON VARIANT #################################
#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
red <- 0.5
parameters$vax2.beta.effic.A.O = red*0.6*rep(1,age.bins)
parameters$vax2.symp.effic.A.O = red*0.81*rep(1,age.bins)
parameters$vax2.hosp.effic.A.O = red*0.90*rep(1,age.bins)
parameters$vax2.death.effic.A.O = red*0.95*rep(1,age.bins) ## GUESS

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.4
parameters$vax1.beta.effic.A.O   = FIRST.DOSE.REL.EFFIC * parameters$vax2.beta.effic.A.O
parameters$vax1.symp.effic.A.O   = FIRST.DOSE.REL.EFFIC * parameters$vax2.symp.effic.A.O
parameters$vax1.hosp.effic.A.O   = FIRST.DOSE.REL.EFFIC * parameters$vax2.hosp.effic.A.O
parameters$vax1.death.effic.A.O  = FIRST.DOSE.REL.EFFIC * parameters$vax2.death.effic.A.O

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                             vax1.beta.effic.A.O,vax1.hosp.effic.A.O,
                                                             vax1.death.effic.A.O,vax1.symp.effic.A.O))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                             vax2.beta.effic.A.O,vax2.hosp.effic.A.O,
                                                             vax2.death.effic.A.O,vax2.symp.effic.A.O))
VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
VAX1.EFFIC.DEATH  <- VAX1.PARS$death
VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
VAX2.EFFIC.DEATH  <- VAX2.PARS$death
VAX2.EFFIC.CLIN   <- VAX2.PARS$symp

#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_vA.O <- (1- VAX1.EFFIC.BETA)*parameters$beta.O
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_vA.O <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_vA.O <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.O
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_vA.O <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr

##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_wA.O <- (1- VAX2.EFFIC.BETA)*parameters$beta.O
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_wA.O <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_wA.O <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.O
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_wA.O <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr

# #############separa historico usando prevalencia inicial########################
prev_hist <- matrix(0,nrow = MAX.TIME.DAYS, ncol = 3*age.bins)
parameters$history.A <- prev_hist
################################################################################
############################ CoronaVac #########################################
# Time window between first and second vaccines
VAX.WINDOW.DAYS = 28
parameters$vax.window.days.C <- VAX.WINDOW.DAYS
parameters$vac.rate.C <- rep(0,MAX.TIME.DAYS)
parameters$vac.rate.v2.C <- rep(0,MAX.TIME.DAYS)
prev_hist <- matrix(0,nrow = MAX.TIME.DAYS, ncol = 3*age.bins)
parameters$history.C <- prev_hist
######################## DELTA VARIANT #########################################
#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
parameters$vax2.beta.effic.C.D = 0.0*rep(1,age.bins)
parameters$vax2.symp.effic.C.D = 0.5*rep(1,age.bins)
parameters$vax2.hosp.effic.C.D = 0.83*rep(1,age.bins)
parameters$vax2.death.effic.C.D = 0.95*rep(1,age.bins) ## GUESS

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.8
parameters$vax1.beta.effic.C.D   = FIRST.DOSE.REL.EFFIC * parameters$vax2.beta.effic.C.D
parameters$vax1.symp.effic.C.D   = FIRST.DOSE.REL.EFFIC * parameters$vax2.symp.effic.C.D
parameters$vax1.hosp.effic.C.D   = FIRST.DOSE.REL.EFFIC * parameters$vax2.hosp.effic.C.D
parameters$vax1.death.effic.C.D  = FIRST.DOSE.REL.EFFIC * parameters$vax2.death.effic.C.D

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                             vax1.beta.effic.C.D,vax1.hosp.effic.C.D,
                                                             vax1.death.effic.C.D,vax1.symp.effic.C.D))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                             vax2.beta.effic.C.D,vax2.hosp.effic.C.D,
                                                             vax2.death.effic.C.D,vax2.symp.effic.C.D))
VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
VAX1.EFFIC.DEATH  <- VAX1.PARS$death
VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
VAX2.EFFIC.DEATH  <- VAX2.PARS$death
VAX2.EFFIC.CLIN   <- VAX2.PARS$symp

#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_vC.D <- (1- VAX1.EFFIC.BETA)*parameters$beta.D
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_vC.D <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_vC.D <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.D
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_vC.D <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr

##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_wC.D <- (1- VAX2.EFFIC.BETA)*parameters$beta.D
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_wC.D <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_wC.D <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.D
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_wC.D <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr

######################## OMICRON VARIANT #########################################
#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
red <- 0.5
parameters$vax2.beta.effic.C.O = red*0.0*rep(1,age.bins)
parameters$vax2.symp.effic.C.O = red*0.5*rep(1,age.bins)
parameters$vax2.hosp.effic.C.O = red*0.83*rep(1,age.bins)
parameters$vax2.death.effic.C.O = red*0.95*rep(1,age.bins) ## GUESS

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.8
parameters$vax1.beta.effic.C.O   = FIRST.DOSE.REL.EFFIC * parameters$vax2.beta.effic.C.O
parameters$vax1.symp.effic.C.O   = FIRST.DOSE.REL.EFFIC * parameters$vax2.symp.effic.C.O
parameters$vax1.hosp.effic.C.O   = FIRST.DOSE.REL.EFFIC * parameters$vax2.hosp.effic.C.O
parameters$vax1.death.effic.C.O  = FIRST.DOSE.REL.EFFIC * parameters$vax2.death.effic.C.O

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                             vax1.beta.effic.C.O,vax1.hosp.effic.C.O,
                                                             vax1.death.effic.C.O,vax1.symp.effic.C.O))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.O,ihfr,asymp,
                                                             vax2.beta.effic.C.O,vax2.hosp.effic.C.O,
                                                             vax2.death.effic.C.O,vax2.symp.effic.C.O))
VAX1.EFFIC.BETA   <- VAX1.PARS$suscep
VAX1.EFFIC.SEVERE <- VAX1.PARS$hosp
VAX1.EFFIC.DEATH  <- VAX1.PARS$death
VAX1.EFFIC.CLIN   <- VAX1.PARS$symp
VAX2.EFFIC.BETA   <- VAX2.PARS$suscep
VAX2.EFFIC.SEVERE <- VAX2.PARS$hosp
VAX2.EFFIC.DEATH  <- VAX2.PARS$death
VAX2.EFFIC.CLIN   <- VAX2.PARS$symp

#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_vC.O <- (1- VAX1.EFFIC.BETA)*parameters$beta.O
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_vC.O <- 1 - (1 - parameters$asymp) * (1-VAX1.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_vC.O <- (1-VAX1.EFFIC.SEVERE)*parameters$ihr.O
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_vC.O <- (1-VAX1.EFFIC.DEATH)*parameters$ihfr

##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
parameters$beta_wC.O <- (1- VAX2.EFFIC.BETA)*parameters$beta.O
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
parameters$asymp_wC.O <- 1 - (1 - parameters$asymp) * (1-VAX2.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
parameters$ihr_wC.O <- (1-VAX2.EFFIC.SEVERE)*parameters$ihr.O
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
parameters$ihfr_wC.O <- (1-VAX2.EFFIC.DEATH)*parameters$ihfr
################################################################################
parameters$vax.cov <- list(vax.cov = FALSE, cov = 1)
parameters$desist <- 0.1



