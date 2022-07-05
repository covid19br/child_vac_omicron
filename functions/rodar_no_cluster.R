source("functions/calculate_init_conditions.R") ###get function that calculates r and eigenvec
source("functions/calculate_init_from_hosp.R") ###recompute init cond from eigenvec and hospitalizations
source("functions/update_R0.R") ###updates R0 from new parameterization
source("functions/samples_to_parameters.R") ###transform samples into parameters
source("functions/sensitivity_analysis_gabriel.R") ###randomize parameters from distributions
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

r_fits <- read.csv("DATA/r_obs.csv")
hosps <- read.csv("DATA/nowcasting_10-semanas.csv")
 
estados = c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
            "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
            "RS", "SC", "SE", "SP", "TO")

cl <- parallel::makeCluster(16)
doParallel::registerDoParallel(cl)

if(!file.exists("output"))
{
  dir.create("output")
}

if(!file.exists("figures"))
{
  dir.create("figures")
}

for (estado in estados) {

   n_hosp <- hosps %>% filter(UF == estado, data == "2022-01-09", dado=="covid") %>%
       select(estimate.merged)
   n_hosp <- n_hosp[[1]]
   r_fits2 <- r_fits %>% filter(UF == estado)
   r_fits2$covid %>% which.max()-> index
   mean_val <- r_fits2[index, "covid"]
   sd_val <- r_fits2[index, "sdcovid"]
   r_fit <- rnorm(nrow(samples), mean = mean_val, sd = sd_val)
   r_fit <- rep(mean_val, nrow(samples))
   source("functions/cobertura_params_reinfec.R") ###starts parameterization
    
   for (vac_children_scenario in 0:2) {

        parameters$vac.rate.children <- vac.rate.all.scn[vac_children_scenario+1,]
        parameters$vac.rate.children.v2 <- vac.rate.v2.all.scn[vac_children_scenario+1,]

        sol  <- foreach(i = 1:nrow(samples),.packages = c("stringr","zeallot","Matrix","rARPACK","tidyr","dplyr")) %dopar% {
        parameters2 <- samples_to_parameters(samples[i,],parameters)
        # prev.graf[i,] <- c(i,calculate_r(0,parameters2),calculate_r(0.5,parameters2),calculate_r(1,parameters2))
        R0.O <- tryCatch({uniroot(function(R0.O){return(return_r(R0.O,parameters2) - r_fit[i])},
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
        parameters2 <- calculate_init_from_hosp(parameters2,n_hosp)
        solution <- solve_system(parameters2,t = 1:90)
        solution %>% wide_solution %>%  filter(t==90, compartment %in% c("D","C","X")) %>% 
            mutate(dose=substr(vac, 1, 1)) %>% group_by(dose, age_group, compartment) %>% 
            summarize(value=sum(value)) %>% ungroup()
        #solution2 <- group_by_dose(solution)
        }
        
        sol2 <- bind_rows(sol, .id = "sample")
        write.csv(sol2, paste0("output/", estado, "_scn", vac_children_scenario, "_X_v2_SA.csv"),
                  row.names=F)
        
        # sol3 <- sol2 %>%
        #     group_by(t, classe) %>%
        #     summarize(mean = mean(value),
        #               q025 = quantile(value, 0.025),
        #               q975 = quantile(value, 0.975))
        
        sol3 <- sol2 %>%
            group_by(dose, age_group, compartment) %>%
            summarize(mean = mean(value),
                      q025 = quantile(value, 0.025),
                      q975 = quantile(value, 0.975))

        write.csv(sol3, paste0("output/quant_", estado, "_scn", vac_children_scenario, "_X_v2_SA.csv"),
                  row.names=F)
        ### plotting
        # sol3 %>%
        # solution2 %>%
          # filter(classe %in% c("Eu", "Ev","Ew","Eb")) %>%
          # filter(classe %in% c("Hu", "Hv","Hw","Hb","Du","Dv","Dw","Db")) %>%
          # filter(classe %in% c("Su", "Sv","Sw","Sb")) %>%
          # filter(classe %in% c("Ru", "Rv","Rw","Rb")) %>%
          # ggplot(aes(x = t, y = mean, color = classe, fill = classe)) + geom_line() +
      #     geom_line() + geom_ribbon(aes(ymin = q025, ymax = q975,color = NULL),  alpha = 0.3) +
      #     labs(x = "data", title = paste0("estado: ", estado, "  cenario: ", vac_children_scenario)) +
      #     theme_bw(base_size = 16) -> plots
      # ggsave(paste0("figures/", estado, "_scn", vac_children_scenario, ".png"), plot = plots, width = 16, height = 12)
    }

   rm(sol, sol2, sol3, parameters2, solution, solution2, plots)
   gc()
}
