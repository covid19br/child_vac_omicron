source("functions/expand_matrix.R")
source("functions/aggregate_age_matrix.R")
source("functions/child_vaccine_age_groups_matrix.R")
source("functions/calculate_R0_reinfec.R")

library(tidyverse)
#########prepara dist etaria
if (! exists("estado")) estado <- "SP"
DistEt <- read_csv("DATA/DistrEtaria2020.csv")
DistEt <- DistEt[-1,]
DistEt <- DistEt[[estado]]

DistEt <- data.frame(pop = DistEt, agegroup = 1:19,agefloor = seq(0,90,5))
#########matriz contato
load('DATA/prem2020_bra.Rdata')
CONTACT.M <- new_all
##matrix with 5-year age bins
CONTACT.M <- expand_matrix(CONTACT.M,DistEt)
##reorganize to work with vaccination age-groups
CONTACT.M <- child_vaccine_age_groups_matrix(CONTACT.M,dist_et = NULL)
##aggregate in 0-4, 5-11, 12-17, 18-29, 30-39,..., 80+
aggregate_indices <- list(1, 2, 3,c(4,5,6),
                          c(7,8),c(9,10),c(11,12),
                          c(13,14),c(15,16),c(17:19))
CONTACT.M<- aggregate_age_matrix(CONTACT.M,aggregate_indices,DistEt$pop)
rm(r_all, r_home, r_others,r_school,r_work,
   u_all, u_home, u_others,u_school,u_work,
   new_all, new_home, new_others,new_school,new_work)
##############PREPARA IHR
ihr <- data.frame(read_csv("DATA/14. ihr.csv"))
ihr[,2] <- as.numeric(ihr[,2])
ihr[,2]<-ihr[,2]/100   # csv data is in percentages
ihr$agegroup <- 1:19
ihr <- left_join(ihr,DistEt)
####changing severity of new agegroups by hand
#5-11
ihr$severe[2] <- (ihr$severe[2]*ihr$pop[2] + ihr$severe[3]*0.40*ihr$pop[3])/(ihr$pop[2] + 0.40*ihr$pop[3])
ihr$pop[2] <- ihr$pop[2] + 0.40*ihr$pop[3]
ihr$pop[3] <- 0.60*ihr$pop[3]
#12-17
ihr$severe[3] <- (ihr$severe[3]*ihr$pop[3] + ihr$severe[4]*0.60*ihr$pop[4])/(ihr$pop[3] + 0.60*ihr$pop[4])
ihr$pop[3] <- ihr$pop[3] + 0.60*ihr$pop[4]
ihr$pop[4] <- 0.40*ihr$pop[4]
ihr$agefloor[3] <- 12
ihr$agefloor[4] <- 18
##the others are automatic
ihr <- ihr %>% mutate(agegroup2 = cut(agefloor, breaks = c(0,5,12,18,30,40,50,60,70,80,Inf),right = F,labels = F))
ihr <- ihr %>% group_by(agegroup2) %>% mutate(severe = severe*pop/sum(pop)) %>% 
  summarize(severe = sum(severe)) %>% rename(agegroup = agegroup2)
ihr <- as.data.frame(ihr)
ihr.D <- ihr$severe
ihr.O <- 0.5*ihr.D
###############PREPARA SINTOMATICOS
symp <- data.frame(read.csv("DATA/symptomatics_per_age.csv"))
symp[,2] <- as.numeric(symp[,2])
symp$agegroup <- 1:19
symp <- left_join(symp,DistEt)
####changing symptomatics of new agegroups by hand
#5-11
symp$symptomatics[2] <- (symp$symptomatics[2]*symp$pop[2] + symp$symptomatics[3]*0.40*symp$pop[3])/(symp$pop[2] + 0.40*symp$pop[3])
symp$pop[2] <- symp$pop[2] + 0.40*symp$pop[3]
symp$pop[3] <- 0.60*symp$pop[3]
#12-17
symp$symptomatics[3] <- (symp$symptomatics[3]*symp$pop[3] + symp$symptomatics[4]*0.60*symp$pop[4])/(symp$pop[3] + 0.60*symp$pop[4])
symp$pop[3] <- symp$pop[3] + 0.60*symp$pop[4]
symp$pop[4] <- 0.40*symp$pop[4]
symp$agefloor[3] <- 12
symp$agefloor[4] <- 18
##the others are automatic
symp <- symp %>% mutate(agegroup2 = cut(agefloor, breaks = c(0,5,12,18,30,40,50,60,70,80,Inf),right = F,labels = F))
symp <- symp %>% group_by(agegroup2) %>% mutate(symptomatics = symptomatics*pop/sum(pop)) %>% 
  summarize(symptomatics = sum(symptomatics)) %>% rename(agegroup = agegroup2)
symp <- as.data.frame(symp)
symp <- symp$symptomatics
asymp <- 1-symp
##############PREPARA OBITOS
ihfr_brasil <- read.csv("DATA/ihfrBrasil.csv")
ihfr <- data.frame(agegroup = 1:19)
ihfr$ihfr <- rep(ihfr_brasil$ihfr,each=2)[1:19]
ihfr$agefloor <- seq(0,90,5)
ihfr <- left_join(ihfr,DistEt)
####changing ihfr of new agegroups by hand
#5-11
ihfr$ihfr[2] <- (ihfr$ihfr[2]*ihfr$pop[2] + ihfr$ihfr[3]*0.40*ihfr$pop[3])/(ihfr$pop[2] + 0.40*ihfr$pop[3])
ihfr$pop[2] <- ihfr$pop[2] + 0.40*ihfr$pop[3]
ihfr$pop[3] <- 0.60*ihfr$pop[3]
#12-17
ihfr$ihfr[3] <- (ihfr$ihfr[3]*ihfr$pop[3] + ihfr$ihfr[4]*0.60*ihfr$pop[4])/(ihfr$pop[3] + 0.60*ihfr$pop[4])
ihfr$pop[3] <- ihfr$pop[3] + 0.60*ihfr$pop[4]
ihfr$pop[4] <- 0.40*ihfr$pop[4]
ihfr$agefloor[3] <- 12
ihfr$agefloor[4] <- 18

ihfr <- ihfr %>% mutate(agegroup2 = cut(agefloor, breaks = c(0,5,12,18,30,40,50,60,70,80,Inf),right = F,labels = F))
ihfr <- ihfr %>% group_by(agegroup2) %>% mutate(ihfr = ihfr*pop/sum(pop)) %>% 
  summarize(ihfr = sum(ihfr)) %>% rename(agegroup = agegroup2)
ihfr <- as.data.frame(ihfr)
ihfr <- ihfr$ihfr
########tempo hospitalização
nus <- read.csv("DATA/hosp_time_br.csv")
nus <- data.frame(agefloor =  seq(0,90,5), nus = rep(nus$mean,each=2)[1:19],agegroup=1:19)
nus <- left_join(nus,DistEt)
####changing nus of new agegroups by hand
#5-11
nus$nus[2] <- (nus$nus[2]*nus$pop[2] + nus$nus[3]*0.40*nus$pop[3])/(nus$pop[2] + 0.40*nus$pop[3])
nus$pop[2] <- nus$pop[2] + 0.40*nus$pop[3]
nus$pop[3] <- 0.60*nus$pop[3]
#12-17
nus$nus[3] <- (nus$nus[3]*nus$pop[3] + nus$nus[4]*0.60*nus$pop[4])/(nus$pop[3] + 0.60*nus$pop[4])
nus$pop[3] <- nus$pop[3] + 0.60*nus$pop[4]
nus$pop[4] <- 0.40*nus$pop[4]
nus$agefloor[3] <- 12
nus$agefloor[4] <- 18
nus <- nus %>% mutate(agegroup2 = cut(agefloor, breaks = c(0,5,12,18,30,40,50,60,70,80,Inf),right = F,labels = F))
nus <- nus %>% group_by(agegroup2) %>% mutate(nus = nus*pop/sum(pop)) %>% 
  summarize(nus = sum(nus)) %>% rename(agegroup = agegroup2)
nus <- as.data.frame(nus)
nus <- 1-exp(-1/nus$nus)
####outros params
omega <- 0.1
beta.D = 0.01
beta.O = 0.01*5
####reinfec params
reinfec.DD <- 0.10
reinfec.DO <- 0.80
reinfec.OD <- 0.10
reinfec.OO <- 0.10
# beta.O <- 3*beta.D
gamma.D <- 1-exp(-1/5.8)
gamma.O <- 1-exp(-1/3.2)
nu <- 1-exp(-1/11)
parameters <- list(omega = omega,
                   omega_s = omega,
                   omega_a = 1,
                   beta.D = beta.D,
                   beta_vA.D = 0, #dummy
                   beta_wA.D = 0,
                   beta_vP.D = 0,
                   beta_wP.D = 0,
                   beta_bP.D = 0,
                   beta_vC.D = 0,
                   beta_wC.D = 0,
                   beta.O = beta.O,
                   beta_vA.O = 0,
                   beta_wA.O = 0,
                   beta_vP.O = 0,
                   beta_wP.O = 0,
                   beta_bP.O = 0,
                   beta_vC.O = 0,
                   beta_wC.O = 0,
                   gamma.D = gamma.D,
                   gamma.O = gamma.O,
                   nu = nu,
                   nus = nus,
                   ihr.D = ihr.D,
                   ihr_vA.D = ihr.D,
                   ihr_wA.D = ihr.D,
                   ihr_vP.D = ihr.D,
                   ihr_wP.D = ihr.D,
                   ihr_bP.D = ihr.D,
                   ihr_vC.D = ihr.D,
                   ihr_wC.D = ihr.D,
                   ihr.O = ihr.O,
                   ihr_vA.O = ihr.O,
                   ihr_wA.O = ihr.O,
                   ihr_vP.O = ihr.O,
                   ihr_wP.O = ihr.O,
                   ihr_bP.O = ihr.O,
                   ihr_vC.O = ihr.O,
                   ihr_wC.O = ihr.O,
                   ihfr = ihfr,
                   ihfr_vA.D = ihfr,
                   ihfr_wA.D = ihfr,
                   ihfr_vP.D = ihfr,
                   ihfr_wP.D = ihfr,
                   ihfr_bP.D = ihfr,
                   ihfr_vC.D = ihfr,
                   ihfr_wC.D = ihfr,
                   ihfr_vA.O = ihfr,
                   ihfr_wA.O = ihfr,
                   ihfr_vP.O = ihfr,
                   ihfr_wP.O = ihfr,
                   ihfr_bP.O = ihfr,
                   ihfr_vC.O = ihfr,
                   ihfr_wC.O = ihfr,
                   asymp = asymp,
                   asymp_vA.D = asymp,
                   asymp_wA.D = asymp,
                   asymp_vP.D = asymp,
                   asymp_wP.D = asymp,
                   asymp_bP.D = asymp,
                   asymp_vC.D = asymp,
                   asymp_wC.D = asymp,
                   asymp_vA.O = asymp,
                   asymp_wA.O = asymp,
                   asymp_vP.O = asymp,
                   asymp_wP.O = asymp,
                   asymp_bP.O = asymp,
                   asymp_vC.O = asymp,
                   asymp_wC.O = asymp,
                   c = CONTACT.M,
                   reinfec.DD = 0.10,
                   reinfec.DO = 0.80,
                   reinfec.OD = 0.10,
                   reinfec.OO = 0.10)
####
DistEt$pop[2] <- DistEt$pop[2] + 0.40*DistEt$pop[3]
DistEt$pop[3] <- 0.60*DistEt$pop[3]
DistEt$pop[3] <- DistEt$pop[3] + 0.60*DistEt$pop[4]
DistEt$pop[4] <- 0.40*DistEt$pop[4]
DistEt$agefloor[3] <- 12
DistEt$agefloor[4] <- 18

DistEt<- DistEt %>% mutate(agegroup2 = cut(agefloor, breaks = c(0,5,12,18,30,40,50,60,70,80,Inf),right = F,labels = F))

DistEt2 <- DistEt %>% group_by(agegroup2) %>% summarize(pop = sum(pop))
###fim params
### cond inicial
pop <- DistEt2$pop
# pop[16] <- sum(pop[16:19])
# pop <- pop[1:16]
age.bins <- parameters$age.bins <- length(pop)
######################################
######################################
##### A === AstraZeneca
##### P === Pfizer
##### C === CoronaVac
######################################
##### D === Delta
##### O === Omicron
######################################
## Definir estados
vac.states <- c("u","vA", "vP", "vC", "wA", "wP", "wC","bP")
inf.states <- c("E", "A", "I", "H", "R", "D","C","X")
var.states <- c(".D",".O")
states = c(paste0("S", vac.states),###susceptiveis vacinados ou não
           paste0(inf.states,rep(vac.states,each = length(inf.states)),
                  rep(var.states,each = length(inf.states)*length(vac.states))) #infectados vacinados ou não
           )

for (v in 1:length(states)){
  eval(parse(text=paste0(states[v], 'index <- seq(', 1+age.bins*(v-1), ',', v*age.bins, ')')))
  eval(parse(text=paste0('init', states[v], ' <- 0*pop')))
  eval(parse(text=paste0("parameters$",states[v],"index <- ",states[v],"index")))
}

initSu <- pop
Y <- eval(parse(text=paste0('c(', paste(paste0('init', states), collapse=','), ')')))
Y <- unlist(Y)
Y <- as.numeric(Y)
names(Y) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))
parameters$init.condition <-  Y
parameters$dist.etaria <- Y
######## updates R0
R.0 <- calculate_R0(parameters)

parameters$beta.D <- parameters$beta.D*(5/R.0$R0.D)
parameters$beta.O <- parameters$beta.O*(5/R.0$R0.O)
parameters$R0.D <- 5
parameters$R0.O <- 5

