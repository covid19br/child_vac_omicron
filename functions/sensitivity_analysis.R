library(rlist)
library(LearnBayes)
library(truncnorm)
R0.O <- list(param = 'R0.O', distr = "truncnorm", vals = c(1,Inf,as.numeric(normal.select(list(p=0.025,x=1),
                                                                               list(p=0.975,x=8)))))
reinfec.DO <- list(param = "reinfec.DO",distr = "beta", vals = beta.select(list(p=0.025,x=0.73),
                                                                           list(p=0.975,x=0.99)))
gamma.O <- list(param = "gamma.O", distr = "truncnorm", vals = c(1, Inf, 2.22, 1.62))
prevalence <- list(param = "prevalence", distr = "truncnorm", vals = c(0.5, 1, 0.7, 0.1))
#########################PFIZER#################################
#########infeccao
#Tseng et al assuming Moderna == Pfizer
##atualização https://www.medrxiv.org/content/10.1101/2022.01.07.22268919v3
#D1
vax1.beta.effic.P.O <- list(param = "vax1.beta.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.095),
                                                                                              list(p=0.975,x=0.300)))
#D2
#New York Health dept
vax2.beta.effic.P.O1 <- list(param = "vax2.beta.effic.P.O1-2", distr = "unif", vals = c(0.87,0.87))
#New York Health dept
vax2.beta.effic.P.O2 <- list(param = "vax2.beta.effic.P.O3", distr = "unif", vals = c(0.61,0.61))
#MMWR (CDC)
vax2.beta.effic.P.O3 <- list(param = "vax2.beta.effic.P.O4-10", distr = "beta", vals = beta.select(list(p=0.025,x=0.38),
                                                                                                list(p=0.975,x=0.435)))
#D3
#MMWR (CDC) -> changed to be able to compute IC
vax3.beta.effic.P.O <- list(param = "vax3.beta.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.81),
                                                                                              list(p=0.975,x=0.85)))
###############
#########sintomatica
###D1
##18- Powell et al
vax1.symp.effic.P.O1 <- list(param = "vax1.symp.effic.P.O1-3", distr = "beta", vals = beta.select(list(p=0.025,x=0.33),
                                                                                              list(p=0.975,x=0.469)))
#18+ Chemaitelly et al
vax1.symp.effic.P.O2 <- list(param = "vax1.symp.effic.P.O4-10", distr = "beta", vals = beta.select(list(p=0.025,x=0.058),
                                                                                              list(p=0.975,x=0.420)))
###D2
###18- Powell et al
vax2.symp.effic.P.O1 <- list(param = "vax2.symp.effic.P.O1-3", distr = "beta", vals = beta.select(list(p=0.025,x=0.644),
                                                                                                   list(p=0.975,x=0.783)))

###18+ UKHSA -> changed to be able to compute IC
vax2.symp.effic.P.O2 <- list(param = "vax2.symp.effic.P.O4-10", distr = "beta", vals = beta.select(list(p=0.025,x=0.64),
                                                                                              list(p=0.975,x=0.68)))

### Andrews et al
vax3.symp.effic.P.O <- list(param = "vax3.symp.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.561),
                                                                                              list(p=0.975,x=0.863)))
#hospitalização
#barnard et al (pessimistic)
vax1.hosp.effic.P.O <- list(param = "vax1.hosp.effic.P.O", distr = "unif", vals = c(0.596, 0.596))

#New York Health dept
vax2.hosp.effic.P.O1 <- list(param = "vax2.hosp.effic.P.O1-2", distr = "unif", vals = c(0.760, 0.760))
#New York Health dept
vax2.hosp.effic.P.O2 <- list(param = "vax2.hosp.effic.P.O3", distr = "unif", vals = c(0.810, 0.810))
#Barnard et al (pessimistic)
vax2.hosp.effic.P.O3 <- list(param = "vax2.hosp.effic.P.O4-10", distr = "beta", vals = beta.select(list(p=0.025,x=0.620),
                                                                                                   list(p=0.975,x=0.659)))
#UKHSA ->  changed to 66,70 to be able to compute
vax3.hosp.effic.P.O <- list(param = "vax3.hosp.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.66),
                                                                                              list(p=0.975,x=0.70)))
#obito
#barnard et al pessimistic
vax1.death.effic.P.O <- list(param = "vax1.death.effic.P.O", distr = "unif", vals = c(0.596, 0.596))
# young Xu et al
vax2.death.effic.P.O <- list(param = "vax2.death.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.620),
                                                                                              list(p=0.975,x=0.850)))
#Young Xu et al
vax3.death.effic.P.O <- list(param = "vax3.death.effic.P.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.910),
                                                                                              list(p=0.975,x=0.980)))

#########################ASTRAZENECA############################
#infeccao
#assumed
vax1.beta.effic.A.O <- list(param = "vax1.beta.effic.A.O", distr = "unif", vals = c(0, 0))
#Willet et al, correcting lower bound
vax2.beta.effic.A.O <- list(param = "vax2.beta.effic.A.O", distr = "beta", vals = beta.select(list(p=0.025,x=1e-10),
                                                                                              list(p=0.975,x=0.346)))
#sintomatica
##Assumed ???
vax1.symp.effic.A.O <- list(param = "vax1.symp.effic.A.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.058),
                                                                                              list(p=0.975,x=0.42)))
##UKHSA
vax2.symp.effic.A.O <- list(param = "vax2.symp.effic.A.O", distr = "beta", vals = beta.select(list(p=0.025,x=0.407),
                                                                                              list(p=0.975,x=0.575)))

#hospitalização
#Barnard et al (pessimistic)
vax1.hosp.effic.A.O <- list(param = "vax1.hosp.effic.A.O", distr = "unif", vals = c(0.497, 0.497))
#Barnard et al (pessimistic)
vax2.hosp.effic.A.O <- list(param = "vax2.hosp.effic.A.O", distr = "unif", vals = c(0.600, 0.600))

#obito
#barnard et al (pessimistic)
vax1.death.effic.A.O <- list(param = "vax1.death.effic.A.O", distr = "unif", vals = c(0.497, 0.497))
#barnard et al (pessimistic)
vax2.death.effic.A.O <- list(param = "vax2.death.effic.A.O", distr = "unif", vals = c(0.600, 0.600))

##########################CORONAVAC############################
#infeccao
vax1.beta.effic.C.O <- list(param = "vax1.beta.effic.C.O", distr = "unif", vals = c(0.0, 0.0))
vax2.beta.effic.C.O <- list(param = "vax2.beta.effic.C.O", distr = "unif", vals = c(0.0, 0.0))

#sintomatica
vax1.symp.effic.C.O <- list(param = "vax1.symp.effic.C.O", distr = "unif", vals = c(0.00, 0.00))
vax2.symp.effic.C.O <- list(param = "vax2.symp.effic.C.O", distr = "unif", vals = c(0.20, 0.20))

#hospitalização
##assumed
vax1.hosp.effic.C.O <- list(param = "vax1.hosp.effic.C.O", distr = "unif", vals = c(0.0, 0.0))
##assumed
vax2.hosp.effic.C.O <- list(param = "vax2.hosp.effic.C.O", distr = "unif", vals = c(0.4065, 0.4065))

#obito
#assumed
vax1.death.effic.C.O <- list(param = "vax1.death.effic.C.O", distr = "unif", vals = c(0.0, 0.0))
#assumed
vax2.death.effic.C.O <- list(param = "vax2.death.effic.C.O", distr = "unif", vals = c(0.4065, 0.4065))

params <- list(
               reinfec.DO = reinfec.DO,
               prevalence = prevalence,
               ##infec pfizer
               #D1
               vax1.beta.effic.P.O = vax1.beta.effic.P.O,
               #D2
               "vax2.beta.effic.P.O1-2" = vax2.beta.effic.P.O1,
               "vax2.beta.effic.P.O3" = vax2.beta.effic.P.O2,
               "vax2.beta.effic.P.O4-10" = vax2.beta.effic.P.O3,
               #D3
               vax3.beta.effic.P.O = vax3.beta.effic.P.O,
               ##symp pfizer
               #D1
               "vax1.symp.effic.P.O1-3" = vax1.symp.effic.P.O1,
               "vax1.symp.effic.P.O4-10" = vax1.symp.effic.P.O2,
               #D2
               "vax2.symp.effic.P.O1-3" = vax2.symp.effic.P.O1,
               "vax2.symp.effic.P.O4-10" = vax2.symp.effic.P.O2,
               #D3
               vax3.symp.effic.P.O = vax3.symp.effic.P.O,
               ##hosp pfizer
               #D1
               vax1.hosp.effic.P.O = vax1.hosp.effic.P.O,
               #D2
               "vax2.hosp.effic.P.O1-2" = vax2.hosp.effic.P.O1,
               "vax2.hosp.effic.P.O3" = vax2.hosp.effic.P.O2,
               "vax2.hosp.effic.P.O4-10" = vax2.hosp.effic.P.O3,
               #D3
               vax3.hosp.effic.P.O = vax3.hosp.effic.P.O,
               #death pfizer
               vax1.death.effic.P.O = vax1.death.effic.P.O,
               vax2.death.effic.P.O = vax2.death.effic.P.O,
               vax3.death.effic.P.O = vax3.death.effic.P.O,
               ####astrazeneca
               #infec
               vax1.beta.effic.A.O = vax1.beta.effic.A.O,
               vax2.beta.effic.A.O = vax2.beta.effic.A.O,
               #symp
               vax1.symp.effic.A.O = vax1.symp.effic.A.O,
               vax2.symp.effic.A.O = vax2.symp.effic.A.O,
               #hosp
               vax1.hosp.effic.A.O = vax1.hosp.effic.A.O,
               vax2.hosp.effic.A.O = vax2.hosp.effic.A.O,
               #death
               vax1.death.effic.A.O = vax1.death.effic.A.O,
               vax2.death.effic.A.O = vax2.death.effic.A.O,
               ####coronavac
               #infec
               vax1.beta.effic.C.O = vax1.beta.effic.C.O,
               vax2.beta.effic.C.O = vax2.beta.effic.C.O,
               #symp
               vax1.symp.effic.C.O = vax1.symp.effic.C.O,
               vax2.symp.effic.C.O = vax2.symp.effic.C.O,
               #hosp
               vax1.hosp.effic.C.O = vax1.hosp.effic.C.O,
               vax2.hosp.effic.C.O = vax2.hosp.effic.C.O,
               #death
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
