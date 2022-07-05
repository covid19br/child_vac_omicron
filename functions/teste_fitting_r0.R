source("functions/cobertura_params_reinfec.R")
source("functions/calculate_init_conditions.R")
source("functions/sensitivity_analysis.R")
source("functions/update_R0.R")
library(doParallel)

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
sol  <- foreach(i = 1:nrow(samples), .combine = rbind,.packages = c("stringr","zeallot","Matrix","rARPACK")) %dopar% {
  parameters2 <- samples_to_parameters(samples[i,],parameters)
  # prev.graf[i,] <- c(i,calculate_r(0,parameters2),calculate_r(0.5,parameters2),calculate_r(1,parameters2))
  prev <- tryCatch({uniroot(function(R0.O){return(return_r(R0.O,parameters2) - 0.1)},
                            interval = c(1,8),
                            tol = 0.01,
                            extendInt = "upX",
                            trace = 2)$root},
                   error = function(cond){
                     message("deu menor que 1")
                     return(1)
                   })
  # prev <- bissec(function(PREV){calculate_r(PREV,parameters2) - 0.1},0,1,n = 100,tol = 0.01)
  print(prev)
  prev.estimada[i] <- prev
  c(i,prev)
  # print(prev.graf[i,])
}
