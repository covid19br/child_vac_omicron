source("functions/expand_matrix.R")
source("functions/aggregate_age_matrix.R")
library(tidyverse)
#########prepara dist etaria
DistEt <- read_csv("DATA/DistrEtaria2020.csv")
DistEt <- DistEt[-1,]
DistEt <- DistEt[["SP"]]
# DistEt[16] <- sum(DistEt[16:19])
# DistEt <- DistEt[1:16]
DistEt <- data.frame(pop = DistEt, agegroup = 1:19)
#########matriz contato
load('DATA/prem2020_bra.Rdata')
CONTACT.M <- new_all
CONTACT.M <- expand_matrix(CONTACT.M,DistEt)
aggregate_indices <- list(c(1,2), c(3,4), c(5:6),
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
ihr <- ihr %>% mutate(agegroup2 = cut(agefloor, breaks = c(0,10,20,30,40,50,60,70,80,Inf),right = F,labels = F))
ihr <- ihr %>% group_by(agegroup2) %>% mutate(severe = severe*pop/sum(pop)) %>% 
  summarize(severe = sum(severe)) %>% rename(agegroup = agegroup2)
ihr <- as.data.frame(ihr)
ihr <- ihr$severe
###############PREPARA SINTOMATICOS
symp <- data.frame(read.csv("DATA/symptomatics_per_age.csv"))
symp[,2] <- as.numeric(symp[,2])
symp$agegroup <- 1:19
symp <- left_join(symp,DistEt)
symp <- symp %>% mutate(agegroup2 = cut(agefloor, breaks = c(0,10,20,30,40,50,60,70,80,Inf),right = F,labels = F))
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
ihfr <- ihfr %>% mutate(agegroup2 = cut(agefloor, breaks = c(0,10,20,30,40,50,60,70,80,Inf),right = F,labels = F))
ihfr <- ihfr %>% group_by(agegroup2) %>% mutate(ihfr = ihfr*pop/sum(pop)) %>% 
  summarize(ihfr = sum(ihfr)) %>% rename(agegroup = agegroup2)
ihfr <- as.data.frame(ihfr)
ihfr <- ihfr$ihfr
########tempo hospitalização
nus <- read.csv("DATA/hosp_time_br.csv")
nus <- data.frame(agefloor =  seq(0,90,5), mean = rep(nus$mean,each=2)[1:19],agegroup=1:19)
nus <- left_join(nus,DistEt)
nus <- nus %>% mutate(agegroup2 = cut(agefloor, breaks = c(0,10,20,30,40,50,60,70,80,Inf),right = F,labels = F))
nus <- nus %>% group_by(agegroup2) %>% mutate(nus = mean*pop/sum(pop)) %>% 
  summarize(nus = sum(nus)) %>% rename(agegroup = agegroup2)
nus <- as.data.frame(nus)
nus <- 1-exp(-1/nus$nus)
####outros params
omega <- 0.1
beta = 0.05
gamma <- 1-exp(-1/5.8)
nu <- 1-exp(-1/11)
parameters <- list(omega = omega,omega_s = omega,omega_a = 1, beta = beta, gamma = gamma, nu = nu, nus = nus,ihr = ihr, ihfr = ihfr,asymp = asymp, c = CONTACT.M)
DistEt<- DistEt %>% mutate(agegroup2 = ceiling((agegroup)/2)) %>% mutate(agegroup2 = if_else(agegroup2 >= 9,9,agegroup2))
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
#####################################
## Definir estados
states = c(c("S", "E", "A", "I", "H", "R", "D"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"vA"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"vP"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"vC"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"wA"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"wP"),
           paste0(c("S", "E", "A", "I", "H", "R", "D"),"wC"))
for (v in 1:length(states)){
  eval(parse(text=paste0(states[v], 'index <- seq(', 1+age.bins*(v-1), ',', v*age.bins, ')')))
  eval(parse(text=paste0('init', states[v], ' <- 0*pop')))
  eval(parse(text=paste0("parameters$",states[v],"index <- ",states[v],"index")))
}

initS <- pop
Y <- eval(parse(text=paste0('c(', paste(paste0('init', states), collapse=','), ')')))
Y <- unlist(Y)
Y <- as.numeric(Y)
names(Y) <- paste0(rep(states,each=age.bins),rep(seq(age.bins),age.bins))
parameters$S_0 <- Y[Sindex]
t = 1:400