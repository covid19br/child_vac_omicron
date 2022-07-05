source("functions/cobertura_params_reinfec.R") ###starts parameterization
source("functions/calculate_init_conditions.R") ###get function that calculates r and eigenvec
source("functions/calculate_init_from_hosp.R") ###recompute init cond from eigenvec and hospitalizations
source("functions/update_R0.R") ###updates R0 from new parameterization
source("functions/samples_to_parameters.R") ###transform samples into parameters
source("functions/sensitivity_analysis.R") ###randomize parameters from distributions
source("functions/solve_model_reinfec.R")
source("functions/func_utils.R")
library(doParallel)
library(tidyverse)
####simple function to return r for root finding
return_r <- function(R0.O, parameters){
  
  parameters <- update_R0(parameters, R0.O = R0.O)
  result <- calculate_init_condition(parameters)
  return(result$r.O)
}


# cl <- parallel::makeCluster(4)
# doParallel::registerDoParallel(cl)
# sol  <- foreach(i = 1:50,.packages = c("stringr","zeallot","Matrix","rARPACK","tidyr","dplyr")) %dopar% {
  parameters2 <- samples_to_parameters(samples[1,],parameters)
  # prev.graf[i,] <- c(i,calculate_r(0,parameters2),calculate_r(0.5,parameters2),calculate_r(1,parameters2))
  R0.O <- tryCatch({uniroot(function(R0.O){return(return_r(R0.O,parameters2) - 0.1)},
                            interval = c(1,8),
                            tol = 0.01,
                            extendInt = "upX",
                            trace = 2)$root},
                   error = function(cond){
                     message("deu menor que 1")
                     return(1)
                   })
  
  # prev <- bissec(function(PREV){calculate_r(PREV,parameters2) - 0.1},0,1,n = 100,tol = 0.01)
  # print(prev)
  parameters2 <- update_R0(parameters2, R0.O = R0.O)
  parameters2 <- calculate_init_from_hosp(parameters2,700)
  solution <- solve_system(parameters2,t = 1:90)
  solution2 <- group_by_dose(solution)
# }

# sol2 <- bind_rows(sol, .id = "sample")

# sol3 <- sol2 %>% group_by(t, classe) %>% summarize(mean = mean(value),
                                                   # q025 = quantile(value, 0.025),
                                                   # q975 = quantile(value, 0.975))
solution2 %>%
  filter(classe %in% c("Eu", "Ev","Ew","Eb")) %>%
  # filter(classe %in% c("Su", "Sv","Sw","Sb")) %>%
  # filter(classe %in% c("Ru", "Rv","Rw","Rb")) %>%
  ggplot(aes(x = t, y = value, color = classe, fill = classe)) + geom_line() + 
  # geom_line() + geom_ribbon(aes(ymin = q025, ymax = q975,color = NULL),  alpha = 0.3) + 
  theme_bw(base_size = 16)


