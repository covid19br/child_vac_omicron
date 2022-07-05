source("functions/vac_params_reinfec.R")
source("functions/organiza_cobertura_vacinal.R")
cobertura <- read_csv("DATA/total_vacinas.csv")
#colnames(cobertura) <- c("UF", "vacina_D1", "dose", "agegroup", "n")
cobertura %>% filter(UF == estado, dose != "DU") -> cobertura

parameters <- organiza_cobertura_vacinal(parameters, cobertura)

dist.children <- read_csv("DATA/proporcao_criancas.csv")

prop.children <- as.numeric(dist.children[[estado]])


#####################################################################
#cenario 1
rollout <- read_csv("DATA/rollout_children.csv")
production <- prop.children*rollout$pfizer
VAX.PRODUCTION.RATE <- c(production, rep(0, 90-length(production)))
VAX.INITIAL.STORAGE.NUM <- 0
MAX.VAC.RATE <- 1e7
VAX.WINDOW.DAYS <- 56
SECOND.VAX.LOSS.FRAC <- 0.1
MAX.TIME.DAYS <- 90
VAX.HISTORY <- rep(0,90)
optimal_solution <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                 MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                                 MAX.TIME.DAYS,VAX.HISTORY)
vac.rate.scn1 <- optimal_solution$OPT.VAX.RATE
history <- optimal_solution$HISTORY
vac.rate.v2.scn1 <- (1-SECOND.VAX.LOSS.FRAC)*c(rep(0,VAX.WINDOW.DAYS),vac.rate.scn1)[1:length(vac.rate.scn1)]
#####################################################################
#cenario 0
vac.rate.scn0 <- 0*vac.rate.scn1
vac.rate.v2.scn0 <- 0*vac.rate.v2.scn1
#####################################################################
#cenario 2
VAX.PRODUCTION.RATE <- prop.children*rep(1e6,90)
VAX.INITIAL.STORAGE.NUM <- 0
MAX.VAC.RATE <- 1e7
VAX.WINDOW.DAYS <- 56
SECOND.VAX.LOSS.FRAC <- 0.1
MAX.TIME.DAYS <- 90
VAX.HISTORY <- rep(0,90)
optimal_solution <- opt_vax_rate(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                                 MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                                 MAX.TIME.DAYS,VAX.HISTORY)
vac.rate.scn2 <- optimal_solution$OPT.VAX.RATE
history <- optimal_solution$HISTORY
vac.rate.v2.scn2 <- (1-SECOND.VAX.LOSS.FRAC)*c(rep(0,VAX.WINDOW.DAYS),vac.rate.scn2)[1:length(vac.rate.scn2)]
#####################################################################
parameters$vac.rate.children <- vac.rate.scn0
parameters$vac.rate.children.v2 <- vac.rate.v2.scn0
vac.rate.all.scn <- rbind(vac.rate.scn0,vac.rate.scn1,vac.rate.scn2)
vac.rate.v2.all.scn <- rbind(vac.rate.v2.scn0,vac.rate.v2.scn1,vac.rate.v2.scn2)
