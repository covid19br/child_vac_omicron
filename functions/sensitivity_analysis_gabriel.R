library(rlist)
library(LearnBayes)
library(truncnorm)

# comentários ------------------------------------------------------------------

# 1) ICs negativos igualados a zero: ver referência abaixo (na vdd igualados a valor mt próximo a zero)
#https://www.stat.berkeley.edu/~stark/SticiGui/Text/confidenceIntervals.htm

# 2) Em caso de valores uniformes já utilizados/Barnard et al - proposta
# - D1 assumido como valor minimo 0
# - D2 assumido como a VE de D1
# - D3 assumido como a VE de D3


# ------------------------------------------------------------------------------


R0.O <- list(param = 'R0.O', distr = "truncnorm", vals = c(1,Inf,as.numeric(normal.select(list(p=0.025,x=1),
                                                                                          list(p=0.975,x=8)))))
reinfec.DO <- list(param = "reinfec.DO",distr = "beta", vals = beta.select(list(p=0.025,x=0.73),
                                                                           list(p=0.975,x=0.99)))
gamma.O <- list(param = "gamma.O", distr = "truncnorm", vals = c(1, Inf, 2.22, 1.62))
prevalence <- list(param = "prevalence", distr = "truncnorm", vals = c(0.5, 1, 0.7, 0.1))
#########################PFIZER#################################
#infeccao
#Barnard et al
vax1.beta.effic.P.O <- list(param = "vax1.beta.effic.P.O", distr = "unif", vals = c(0.187,0.187))

vax2.beta.effic.P.O <- list(param = "vax2.beta.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.139),
                                                                                              list(p=0.975,x=0.478)))

vax3.beta.effic.P.O <- list(param = "vax3.beta.effic.P.O", distr = "beta",  vals = beta.select(list(p=0.025,x=0.381),
                                                                                               list(p=0.975,x=0.478)))

#sintomatica
# D1 não temos novas estimativas

vax1.symp.effic.P.O <- list(param = "vax1.symp.effic.P.O", distr = "beta", vals = c(0.187,0.187))

vax2.symp.effic.P.O <- list(param = "vax2.symp.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.187),
                                                                                              list(p=0.975,x=0.279)))

vax3.symp.effic.P.O <- list(param = "vax3.symp.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.279),
                                                                                              list(p=0.975,x=0.469)))


#hospitalização
vax1.hosp.effic.P.O <- list(param = "vax1.hosp.effic.P.O", distr = "unif", vals = c(0.596, 0.596))

vax2.hosp.effic.P.O <- list(param = "vax2.hosp.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.543),
                                                                                              list(p=0.975,x=0.651)))

vax3.hosp.effic.P.O <- list(param = "vax3.hosp.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.826),
                                                                                              list(p=0.975,x=0.879)))
#obito
vax1.death.effic.P.O <- list(param = "vax1.death.effic.P.O", distr = "unif", vals = c(0.596, 0.596))

vax2.death.effic.P.O <- list(param = "vax2.death.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.620),
                                                                                                list(p=0.975,x=0.850)))

vax3.death.effic.P.O <- list(param = "vax3.death.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.826),
                                                                                                list(p=0.975,x=0.879)))

#########################ASTRAZENECA############################
#infeccao
# D1 não temos novas estimativas

vax1.beta.effic.A.O <- list(param = "vax1.beta.effic.A.O", distr = "unif", vals = c(0.187, 0.187))
vax2.beta.effic.A.O <- list(param = "vax2.beta.effic.A.O", distr = "beta", vals = beta.select(list(p=0.025,x= 1e-10),
                                                                                              list(p=0.975,x=0.43)))

#sintomatica
vax1.symp.effic.A.O <- list(param = "vax1.symp.effic.A.O", distr = "unif", vals = c(0.187,0.187))

vax2.symp.effic.A.O <- list(param = "vax2.symp.effic.A.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.0117),
                                                                                              list(p=0.975,x=0.4279)))

#hospitalização
# Aqui não temos novas estimativas - mas coloquei como eficácia zero na D1 ou a eficácia da primeira dose na D2

vax1.hosp.effic.A.O <- list(param = "vax1.hosp.effic.A.O", distr = "unif", vals = c(0.497, 0.497))
vax2.hosp.effic.A.O <- list(param = "vax2.hosp.effic.A.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.497),
                                                                                              list(p=0.975,x=0.837)))

#obito 
# Aqui não temos novas estimativas

vax1.death.effic.A.O <- list(param = "vax1.death.effic.A.O", distr = "unif", vals = c(0.497, 0.497))

vax2.death.effic.A.O <- list(param = "vax2.death.effic.A.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.497),
                                                                                                list(p=0.975,x=0.837)))

##########################CORONAVAC############################

# Aqui não temos novas estimativas
#infeccao
vax1.beta.effic.C.O <- list(param = "vax1.beta.effic.C.O", distr = "unif", vals = c(0.0, 0.0))
vax2.beta.effic.C.O <- list(param = "vax2.beta.effic.C.O", distr = "unif", vals = c(0.0, 0.0))

#sintomatica
vax1.symp.effic.C.O <- list(param = "vax1.symp.effic.C.O", distr = "unif", vals = c(0.00, 0.00))
vax2.symp.effic.C.O <- list(param = "vax2.symp.effic.C.O", distr = "beta", vals = beta.select(list(p=0.025,x=1e-10),
                                                                                              list(p=0.975,x=0.20)))


#hospitalização
vax1.hosp.effic.C.O <- list(param = "vax1.hosp.effic.C.O", distr = "unif", vals = c(0.0, 0.0))
vax2.hosp.effic.C.O <- list(param = "vax2.hosp.effic.C.O", distr = "beta", vals = beta.select(list(p=0.025,x=1e-10),
                                                                                              list(p=0.975,x=0.4065)))

#obito
vax1.death.effic.C.O <- list(param = "vax1.death.effic.C.O", distr = "unif", vals = c(0.0, 0.0))
vax2.death.effic.C.O <- list(param = "vax2.death.effic.C.O", distr = "beta", vals = beta.select(list(p=0.025,x=1e-10),
                                                                                                list(p=0.975,x=0.4065)))

params <- list(
  # R0.O = R0.O,
  reinfec.DO = reinfec.DO,
  prevalence = prevalence,
  # gamma.O = gamma.O,
  vax1.beta.effic.P.O = vax1.beta.effic.P.O,
  vax2.beta.effic.P.O = vax2.beta.effic.P.O,
  vax3.beta.effic.P.O = vax3.beta.effic.P.O,
  vax1.symp.effic.P.O = vax1.symp.effic.P.O,
  vax2.symp.effic.P.O = vax2.symp.effic.P.O,
  vax3.symp.effic.P.O = vax3.symp.effic.P.O,
  vax1.hosp.effic.P.O = vax1.hosp.effic.P.O,
  vax2.hosp.effic.P.O = vax2.hosp.effic.P.O,
  vax3.hosp.effic.P.O = vax3.hosp.effic.P.O,
  vax1.death.effic.P.O = vax1.death.effic.P.O,
  vax2.death.effic.P.O = vax2.death.effic.P.O,
  vax3.death.effic.P.O = vax3.death.effic.P.O,
  vax1.beta.effic.A.O = vax1.beta.effic.A.O,
  vax2.beta.effic.A.O = vax2.beta.effic.A.O,
  vax1.symp.effic.A.O = vax1.symp.effic.A.O,
  vax2.symp.effic.A.O = vax2.symp.effic.A.O,
  vax1.hosp.effic.A.O = vax1.hosp.effic.A.O,
  vax2.hosp.effic.A.O = vax2.hosp.effic.A.O,
  vax1.death.effic.A.O = vax1.death.effic.A.O,
  vax2.death.effic.A.O = vax2.death.effic.A.O,
  vax1.beta.effic.C.O = vax1.beta.effic.C.O,
  vax2.beta.effic.C.O = vax2.beta.effic.C.O,
  vax1.symp.effic.C.O = vax1.symp.effic.C.O,
  vax2.symp.effic.C.O = vax2.symp.effic.C.O,
  vax1.hosp.effic.C.O = vax1.hosp.effic.C.O,
  vax2.hosp.effic.C.O = vax2.hosp.effic.C.O,
  vax1.death.effic.C.O = vax1.death.effic.C.O,
  vax2.death.effic.C.O = vax2.death.effic.C.O
)
# p = p, gamma = gamm,
n.samples <- 1000
set.seed(42)
n.params <- length(names(params))
samples <- matrix(0,nrow = n.samples,ncol = n.params)
samples <-as.data.frame(samples)
colnames(samples) <- names(params)

for (param in params) {
  distr<- param$distr
  if(distr == "unif"){
    samples[,param$param] <- runif(n.samples,min = param$vals[1],param$vals[2])
  }
  if(distr == "normal"){
    samples[,param$param] <- rnorm(n.samples,mean = param$vals[1],sd = param$vals[2])
  }
  if(distr == "gamma"){
    samples[,param$param] <- rgamma(n.samples,shape = param$vals[1],rate = param$vals[2])
  }
  if(distr == "exp"){
    samples[,param$param] <- rexp(n.samples,rate = param$vals[1])
  }
  if(distr == "beta"){
    samples[,param$param] <- rbeta(n.samples,param$vals[1],param$vals[2])
  }
  if(distr == "betabinom"){
    samples[,param$param] <- 1-rbeta(n.samples,param$vals[1],param$vals[2])/rbeta(n.samples,param$vals[3],param$vals[4])
  }
  if(distr == "betaneg"){
    samples[,param$param] <- 2*rbeta(n.samples,param$vals[1],param$vals[2])-1
  }
  if(distr == "truncnorm"){
    samples[,param$param] <- rtruncnorm(n.samples,a = param$vals[1],b = param$vals[2], mean = param$vals[3],sd = param$vals[4])
  }
}

