calculate_init_from_hosp <- function(parameters, hospitalizations){
  init.condit.O <- calculate_init_condition(parameters)$init.condit.O
  #hospitalizations are the last quarter of the initial condition
  len.init <- length(init.condit.O)
 
  ##organize infection-hospitalization rate by vaccine order
  ihr <- with(parameters,{
    ihr <- c(ihr.O, 
                 ihr_vA.O, ihr_vP.O, ihr_vC.O,
                 ihr_wA.O, ihr_wP.O, ihr_wC.O,
                 ihr_bP.O)
  })
  ##recover incubation period
  gamma <- -1/log(1-parameters$gamma.O)
  gamma.week <- gamma/7
  exp.vec <- init.condit.O[1:(len.init/4)]
  hosp.factor <- hospitalizations / sum(exp.vec * ihr / gamma.week)
  ###renormalize entire vector based on exposed
  renorm.vec <- init.condit.O*(rep(hosp.factor, 4))
  
  ###insert renorm.vec in initial conditions
  Y <- parameters$init.condition
  Y2 <- with(parameters,{
    S <- Y[c(Suindex,SvAindex, SvPindex, SvCindex, SwAindex, SwPindex, SwCindex, SbPindex)]
    R.D <- Y[c(Ru.Dindex,RvA.Dindex, RvP.Dindex, RvC.Dindex, RwA.Dindex, RwP.Dindex, RwC.Dindex, RbP.Dindex)]
    E.O <- renorm.vec[1:(len.init/4)]
    A.O <- renorm.vec[(len.init/4 + 1):(2*(len.init/4))]
    I.O <- renorm.vec[((2*(len.init/4)) + 1):(3*(len.init/4))]
    H.O <- renorm.vec[((3*(len.init/4)) + 1):len.init]
    S <- S - (E.O + A.O + I.O + H.O) * S/(S+R.D+0.1)
    R.D <- R.D - (E.O + A.O + I.O + H.O) * R.D/(S+R.D+0.1)
    Y[c(Suindex,SvAindex, SvPindex, SvCindex, SwAindex, SwPindex, SwCindex, SbPindex)] <- S
    Y[c(Eu.Oindex, EvA.Oindex, EvP.Oindex, EvC.Oindex, EwA.Oindex, EwP.Oindex, EwC.Oindex, EbP.Oindex)] <- E.O
    Y[c(Au.Oindex, AvA.Oindex, AvP.Oindex, AvC.Oindex, AwA.Oindex, AwP.Oindex, AwC.Oindex, AbP.Oindex)] <- A.O
    Y[c(Iu.Oindex, IvA.Oindex, IvP.Oindex, IvC.Oindex, IwA.Oindex, IwP.Oindex, IwC.Oindex, IbP.Oindex)] <- I.O
    Y[c(Hu.Oindex, HvA.Oindex, HvP.Oindex, HvC.Oindex, HwA.Oindex, HwP.Oindex, HwC.Oindex, HbP.Oindex)] <- H.O
    Y[c(Ru.Dindex,RvA.Dindex, RvP.Dindex, RvC.Dindex, RwA.Dindex, RwP.Dindex, RwC.Dindex, RbP.Dindex)] <- R.D
    Y
  })
  parameters$init.condition <- Y2
  return(parameters)
}
