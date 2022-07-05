source("functions/epi_params_omicron.R")
source("functions/opt_vax_rate_history.R")
source("functions/VAX_DISTR_RATE.R")
library(dplyr)

## functions ##

#### convert efficacy from studies to parameters
OBSERVED.EFFICACY.TO.PARAMETERS <- function(HOSP.PROP,MORT.PROP,ASYMP.PROP,SUSCEP.EFF,HOSP.EFF,MORT.EFF,SYMP.EFF){
  SUSCEP.EFF.PAR <- SUSCEP.EFF
  HOSP.EFF.PAR <- 1 - (1-HOSP.EFF)/(1-SUSCEP.EFF)
  MORT.EFF.PAR <- 1 - (1-MORT.EFF)/(1-HOSP.EFF)
  SYMP.EFF.PAR <- 1 - ((1-SYMP.EFF)*(HOSP.PROP + (1 - HOSP.PROP)*(1-ASYMP.PROP))  - (1 - HOSP.EFF)*HOSP.PROP)/
    ((1-SUSCEP.EFF)*(1 - (1-HOSP.EFF)/(1-SUSCEP.EFF)*HOSP.PROP)*(1-ASYMP.PROP))
  return(list(suscep = SUSCEP.EFF.PAR,hosp = HOSP.EFF.PAR,death = MORT.EFF.PAR,symp = SYMP.EFF.PAR))
}


## global parameters ##

POP.CITY.REL.FRAC = 0.25 #population of SP state compared to Brazil
MAX.VAC.RATE = 5e6 # Max number of vaccine applications per day
MAX.VAC.RATE2 <- MAX.VAC.RATE
#####PREVALENCIA ARTIFICIAL ----- assumindo somente prev em rel a Delta ########
PREV <- 0.5
# Fraction of people that take the first but not the second dose \theta
SECOND.VAX.LOSS.FRAC = 0.1
parameters$desist <- SECOND.VAX.LOSS.FRAC

# tabela de entrega de vacinas
rollout <- read.csv("./DATA/rollout_vacinas.csv")

data_base <- as.Date("2021-08-09")
rollout$date <- as.Date(rollout$date)
rollout <- rollout %>% filter(date >= data_base)
rollout[,1:3] <- rollout[,1:3]
date_range <- as.integer(difftime(max(rollout$date),min(rollout$date),units = "days"))

# time until end of vaccination schedule, in days
MAX.TIME.DAYS = min(300, date_range)

# tabelas de históricos
# estados <- c("RO", "AC", "AM", "RR", "PA", "AP", "TO", "MA", "PI", "CE", "RN",
#              "PB", "PE", "AL", "SE", "BA", "MG", "ES", "RJ", "SP", "PR", "SC",
#              "RS", "MS", "MT", "GO", "DF")
estados <- "SP"
hist_D1 <- list()
hist_D2 <- list()
hist_D1_acum <- list()
for (estado in estados) {
  hist_D1[[estado]] <- read.csv(paste0("DATA/historico/historico_D1_", estado, ".csv"))
  hist_D2[[estado]] <- read.csv(paste0("DATA/historico/historico_D2_", estado, ".csv"))
  hist_D1_acum[[estado]] <- hist_D1[[estado]] %>% 
    group_by(vacina, agegroup, tempo_doses) %>%
    summarize(n = sum(n)) %>%
    as.data.frame()
}
# agrupando BR todo
hist_D1 <- bind_rows(hist_D1, .id = "estado")
hist_D1 <- hist_D1 %>% mutate(agegroup = if_else(agegroup >= 9.0, 9.0, 1.0*agegroup)) %>%
  group_by(vacina, agegroup, tempo_doses) %>%
  summarize(n = sum(n)) %>%
  as.data.frame()
hist_D1_acum <- bind_rows(hist_D1_acum, .id = "estado")
hist_D1_acum <- hist_D1_acum %>% mutate(agegroup = if_else(agegroup >= 9.0, 9.0, 1.0*agegroup)) %>%
  group_by(vacina, agegroup) %>%
  summarize(dose1 = sum(n)) %>%
  as.data.frame()
hist_D2 <- bind_rows(hist_D2, .id = "estado")
hist_D2 <- hist_D2 %>% mutate(agegroup = if_else(agegroup >= 9.0, 9.0, 1.0*agegroup)) %>%
  group_by(vacina, agegroup) %>%
  summarize(dose2 = sum(n)) %>%
  as.data.frame()
hist_D1D2 <- full_join(hist_D1_acum, hist_D2, by = c("vacina", "agegroup"))
hist_D1D2 <- hist_D1D2 %>%
  filter(!is.na(agegroup))
hist_D1D2$dose2[is.na(hist_D1D2$dose2)]<-0
## vaccine-specific parameters ##
# hist_D1$n<- hist_D1*POP.CITY.REL.FRAC
# hist_D1D2[,c("dose1","dose2")] <- hist_D1D2[,c("dose1","dose2")]*POP.CITY.REL.FRAC

list.parameters <- list()

################################################################################
############################### Pfizer #########################################

VAX.INITIAL.STORAGE.NUM = 0 # 6e6 # Number of READY vaccines on day one
#VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
VAX.PRODUCTION.RATE <- rollout$pfizer

# above numbers are for whole country, correcting (roughly) by population
VAX.INITIAL.STORAGE.NUM = POP.CITY.REL.FRAC * VAX.INITIAL.STORAGE.NUM
VAX.PRODUCTION.RATE = POP.CITY.REL.FRAC * VAX.PRODUCTION.RATE
MAX.VAC.RATE = POP.CITY.REL.FRAC * MAX.VAC.RATE2

# Time window between first and second vaccines
VAX.WINDOW.DAYS = 56

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
parameters$vax.window.days.P <- VAX.WINDOW.DAYS
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
######mexendo com historico
Y2 <- Y
dados_doses <- hist_D1D2 %>% filter(vacina == "Pfizer")
dados_historico <- hist_D1 %>% filter(vacina == "Pfizer")
dados_historico_2 <- dados_historico
dados_dose_1 <- dados_historico_2 %>% group_by(agegroup) %>% summarise(n = sum(n))
# dados_dose_1 <- rbind(c(1,0),c(2,0),dados_dose_1)

Y2[Suindex] <- Y2[Suindex] - dados_doses$dose2
Y2[SwPindex] <- Y2[SwPindex] + dados_doses$dose2
Y2[Suindex] <- Y2[Suindex] - dados_dose_1$n
Y2[SvPindex] <- Y2[SvPindex] + dados_dose_1$n
##########################
hist <- dados_historico_2 %>%group_by(tempo_doses) %>% summarize(n = sum(n))
VAX.HISTORY <- as.numeric(hist$n)

optimal_solution <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                 MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                                 MAX.TIME.DAYS,VAX.HISTORY)
vac.rate <- optimal_solution$OPT.VAX.RATE
history <- optimal_solution$HISTORY
vac.rate.v2 <- (1-SECOND.VAX.LOSS.FRAC)*c(rep(0,VAX.WINDOW.DAYS),vac.rate)[1:length(vac.rate)]
parameters$vac.rate.P <- vac.rate
parameters$vac.rate.v2.P <- vac.rate.v2
df <- generate_vac_schedule(VAX.INITIAL.STORAGE.NUM,VAX.PRODUCTION.RATE,
                            MAX.VAC.RATE,VAX.WINDOW.DAYS,SECOND.VAX.LOSS.FRAC,MAX.TIME.DAYS,VAX.HISTORY)
df$Window <- "Pfizer"
plot_vac_schedule.3(df)

#############separa historico em faixa etaria###################################
dados_historico_2<-dados_historico_2[with(dados_historico_2,order(-tempo_doses,-agegroup)),]
age_hist <- matrix(0,nrow = MAX.TIME.DAYS+1,ncol = age.bins)
j <- 2
for(i in 1:nrow(dados_historico_2)){
  val <- dados_historico_2[i,]
  while(val$n > 0){
    if(sum(age_hist[j]) >= history[j]){ ###se já preencheu a linha
      if(j < MAX.TIME.DAYS+1){##e j não está no fim, incrementa
        j <- j+1
      }
      else{###se j está no fim, quebra o loop, colocar um warning aqui?
        break
      }
    }
    if(sum(age_hist[j,]) + val$n <= history[j]){###se cabe todas as pessoas
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val$n
      val$n <- 0
    }
    else{####se não cabe
      val2 <- history[j] - sum(age_hist[j,])
      # print(val2)
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val2
      val$n <- val$n-val2
      j <- j+1
    }
  }###fim while
}
#############separa historico usando prevalencia inicial########################
prev_hist <- cbind((1-PREV)*age_hist,PREV*age_hist,0*age_hist)
parameters$history.P <- prev_hist
################################################################################
############################ AstraZeneca #######################################

VAX.INITIAL.STORAGE.NUM = 0 #6e6 # Number of READY vaccines on day one
#VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
VAX.PRODUCTION.RATE <- rollout$astrazeneca

# above numbers are for whole country, correcting (roughly) by population
VAX.INITIAL.STORAGE.NUM = POP.CITY.REL.FRAC * VAX.INITIAL.STORAGE.NUM
VAX.PRODUCTION.RATE = POP.CITY.REL.FRAC * VAX.PRODUCTION.RATE
MAX.VAC.RATE = POP.CITY.REL.FRAC * MAX.VAC.RATE2

# Time window between first and second vaccines
VAX.WINDOW.DAYS = 56
parameters$vax.window.days.A <- VAX.WINDOW.DAYS

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
VAX1.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
                                                             vax1.beta.effic.A.O,vax1.hosp.effic.A.O,
                                                             vax1.death.effic.A.O,vax1.symp.effic.A.O))
VAX2.PARS <- with(parameters,OBSERVED.EFFICACY.TO.PARAMETERS(ihr.D,ihfr,asymp,
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

###########################mexendo com historico################################
dados_doses <- hist_D1D2 %>% filter(vacina == "AZ")
dados_historico <- hist_D1 %>% filter(vacina == "AZ")
dados_historico_2 <- dados_historico 
dados_dose_1 <- dados_historico_2 %>% group_by(agegroup) %>% summarise(n = sum(n))
# dados_dose_1 <- rbind(c(1,0),c(2,0),dados_dose_1)
dados_dose_1 <- dados_dose_1[!is.na(dados_dose_1$agegroup),]
Y2[Suindex] <- Y2[Suindex] - dados_doses$dose2
Y2[SwAindex] <- Y2[SwAindex] + dados_doses$dose2
Y2[Suindex] <- Y2[Suindex] - dados_dose_1$n
Y2[SvAindex] <- Y2[SvAindex] + dados_dose_1$n
######dados outras vacinas
# dados_doses <- hist_D1D2 %>% filter(vacina != "AZ") %>% group_by(agegroup) %>% summarize(n = sum(dose1+dose2))
# Y2[Sindex] <- Y2[Sindex] - dados_doses$n

hist <- dados_historico_2 %>%group_by(tempo_doses) %>% summarize(n = sum(n))
VAX.HISTORY <- as.numeric(hist$n)

optimal_solution <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                 MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                                 MAX.TIME.DAYS,VAX.HISTORY)
vac.rate <- optimal_solution$OPT.VAX.RATE
history <- optimal_solution$HISTORY
vac.rate.v2 <- (1-SECOND.VAX.LOSS.FRAC)*c(rep(0,VAX.WINDOW.DAYS),vac.rate)[1:length(vac.rate)]
parameters$vac.rate.A <- vac.rate
parameters$vac.rate.v2.A <- vac.rate.v2
df.2 <- generate_vac_schedule(VAX.INITIAL.STORAGE.NUM,VAX.PRODUCTION.RATE,
                              MAX.VAC.RATE,VAX.WINDOW.DAYS,SECOND.VAX.LOSS.FRAC,MAX.TIME.DAYS,VAX.HISTORY)
df.2$Window <- "AstraZeneca"
plot_vac_schedule.3(df.2)
#############separa historico em faixa etaria###################################
dados_historico_2<-dados_historico_2[with(dados_historico_2,order(-tempo_doses,-agegroup)),]
age_hist <- matrix(0,nrow = MAX.TIME.DAYS+1,ncol = age.bins)
j <- 2
for(i in 1:nrow(dados_historico_2)){
  val <- dados_historico_2[i,]
  while(val$n > 0){
    if(sum(age_hist[j]) >= history[j]){ ###se já preencheu a linha
      if(j < MAX.TIME.DAYS+1){##e j não está no fim, incrementa
        j <- j+1
      }
      else{###se j está no fim, quebra o loop, colocar um warning aqui?
        break
      }
    }
    if(sum(age_hist[j,]) + val$n <= history[j]){###se cabe todas as pessoas
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val$n
      val$n <- 0
    }
    else{####se não cabe
      val2 <- history[j] - sum(age_hist[j,])
      # print(val2)
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val2
      val$n <- val$n-val2
      j <- j+1
    }
  }###fim while
}
#############separa historico usando prevalencia inicial########################
prev_hist <- cbind((1-PREV)*age_hist,PREV*age_hist,0*age_hist)
parameters$history.A <- prev_hist
################################################################################
############################ CoronaVac #########################################
VAX.INITIAL.STORAGE.NUM = 0 # 6e6 # Number of READY vaccines on day one
#VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
VAX.PRODUCTION.RATE <- rollout$coronavac

# above numbers are for whole country, correcting (roughly) by population
VAX.INITIAL.STORAGE.NUM = POP.CITY.REL.FRAC * VAX.INITIAL.STORAGE.NUM
VAX.PRODUCTION.RATE = POP.CITY.REL.FRAC * VAX.PRODUCTION.RATE
MAX.VAC.RATE = POP.CITY.REL.FRAC * MAX.VAC.RATE2

# Time window between first and second vaccines
VAX.WINDOW.DAYS = 28
parameters$vax.window.days.C <- VAX.WINDOW.DAYS

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

########################mexendo com historico###################################
dados_doses <- hist_D1D2 %>% filter(vacina == "Coronavac")
dados_historico <- hist_D1 %>% filter(vacina == "Coronavac")
dados_historico_2 <- dados_historico
dados_dose_1 <- dados_historico_2 %>% group_by(agegroup) %>% summarise(n = sum(n))
# dados_dose_1 <- rbind(c(1,0),c(2,0),dados_dose_1)

Y2[Suindex] <- Y2[Suindex] - dados_doses$dose2
Y2[SwCindex] <- Y2[SwCindex] + dados_doses$dose2
Y2[Suindex] <- Y2[Suindex] - dados_dose_1$n
Y2[SvCindex] <- Y2[SvCindex] + dados_dose_1$n
# ######dados outras vacinas
# dados_doses <- hist_D1D2 %>% filter(vacina != "Coronavac") %>% group_by(agegroup) %>% summarize(n = sum(dose1+dose2))
# Y2[Sindex] <- Y2[Sindex] - dados_doses$n

hist <- dados_historico_2 %>%group_by(tempo_doses) %>% summarize(n = sum(n))
VAX.HISTORY <- as.numeric(hist$n)

optimal_solution <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                 MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                                 MAX.TIME.DAYS,VAX.HISTORY)
vac.rate <- optimal_solution$OPT.VAX.RATE
history <- optimal_solution$HISTORY
vac.rate.v2 <- (1-SECOND.VAX.LOSS.FRAC)*c(rep(0,VAX.WINDOW.DAYS),vac.rate)[1:length(vac.rate)]
parameters$vac.rate.C <- vac.rate
parameters$vac.rate.v2.C <- vac.rate.v2
df.3 <- generate_vac_schedule(VAX.INITIAL.STORAGE.NUM,VAX.PRODUCTION.RATE,
                              MAX.VAC.RATE,VAX.WINDOW.DAYS,SECOND.VAX.LOSS.FRAC,MAX.TIME.DAYS,VAX.HISTORY)
df.3$Window <- "CoronaVac"
plot_vac_schedule.3(df.3)
#############separa historico em faixa etaria###################################
dados_historico_2<-dados_historico_2[with(dados_historico_2,order(-tempo_doses,-agegroup)),]
age_hist <- matrix(0,nrow = MAX.TIME.DAYS+1,ncol = age.bins)
j <- 2
for(i in 1:nrow(dados_historico_2)){
  val <- dados_historico_2[i,]
  while(val$n > 0){
    if(sum(age_hist[j]) >= history[j]){ ###se já preencheu a linha
      if(j < MAX.TIME.DAYS+1){##e j não está no fim, incrementa
        j <- j+1
      }
      else{###se j está no fim, quebra o loop, colocar um warning aqui?
        break
      }
    }
    if(sum(age_hist[j,]) + val$n <= history[j]){###se cabe todas as pessoas
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val$n
      val$n <- 0
    }
    else{####se não cabe
      val2 <- history[j] - sum(age_hist[j,])
      # print(val2)
      age_hist[j,val$agegroup] <- age_hist[j,val$agegroup] + val2
      val$n <- val$n-val2
      j <- j+1
    }
  }###fim while
}
#############separa historico usando prevalencia inicial########################
prev_hist <- cbind((1-PREV)*age_hist,PREV*age_hist,0*age_hist)
parameters$history.C <- prev_hist

###############################################################################
Y2[Ru.Dindex] <- PREV*Y2[Suindex]
Y2[Suindex] <- (1-PREV)*Y2[Suindex]
Y2[RwP.Dindex] <- PREV*Y2[SwPindex]
Y2[SwPindex] <- (1-PREV)*Y2[SwPindex]
Y2[RvP.Dindex] <- PREV*Y2[SvPindex]
Y2[SvPindex] <- (1-PREV)*Y2[SvPindex]
Y2[RwA.Dindex] <- PREV*Y2[SwAindex]
Y2[SwAindex] <- (1-PREV)*Y2[SwAindex]
Y2[RvA.Dindex] <- PREV*Y2[SvAindex]
Y2[SvAindex] <- (1-PREV)*Y2[SvAindex]
Y2[RwC.Dindex] <- PREV*Y2[SwCindex]
Y2[SwCindex] <- (1-PREV)*Y2[SwCindex]
Y2[RvC.Dindex] <- PREV*Y2[SvCindex]
Y2[SvCindex] <- (1-PREV)*Y2[SvCindex]
###### adiciona infectado
Y2[Eu.Dindex[3]] <- 1
Y2[Eu.Oindex[3]] <- 1
# Y[Eu.Dindex[3]] <- 10000
# Y[Eu.Oindex[3]] <- 1

parameters$init.condition <-  Y2
###############################################################################
parameters$vax.cov <- list(vax.cov = FALSE, cov = 1)

parameters$prevalence <- PREV