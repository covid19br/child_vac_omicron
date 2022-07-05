library(tidyverse)

efeito_indireto <- function(){
  estados <- c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
                "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
                "RS", "SC", "SE", "SP", "TO")
  DistEt <- read_csv("DATA/pop_somada.csv")
  prop_reforco <- read_csv("DATA/proporcao_reforco_estados.csv")
  colnames(prop_reforco)[1] <- "estado"
  prop_reforco %>% select(estado, porcentagem_reforco) -> prop_reforco
  all_scn <- c()
  for(uf in estados)
  {
    scn0 <- read_csv(paste0("output/", uf,"_scn0.csv"))
    scn1 <- read_csv(paste0("output/", uf,"_scn1.csv"))
    scn2 <- read_csv(paste0("output/", uf,"_scn2.csv"))
    scn0$scenario <- "scn0"
    scn1$scenario <- "scn1"
    scn2$scenario <- "scn2"
    scn0$estado <- uf
    scn1$estado <- uf
    scn2$estado <- uf
    all_scn <- rbind(all_scn, scn0, scn1,scn2)
  }
  all_scn %>% group_by(sample, compartment, age_group, scenario,estado) %>% 
    summarize(value = sum(value)) -> all_scn
  all_scn %>% mutate(age_group = case_when(age_group == 1 ~ "0-4",
                                           age_group == 2 ~ "5-11",
                                           age_group == 3 ~ "12-17",
                                           age_group == 4 ~ "18-29",
                                           age_group == 5 ~ "30-39",
                                           age_group == 6 ~ "40-49",
                                           age_group == 7 ~ "50-59",
                                           age_group == 8 ~ "60-69",
                                           age_group == 9 ~ "70-79",
                                           age_group == 10 ~ "80+")) %>%
    mutate(age_group = factor(age_group, levels = c("0-4","5-11","12-17","18-29","30-39",
                                                    "40-49","50-59","60-69","70-79","80+"))) -> all_scn
  
  all_scn %>% pivot_wider(names_from = "scenario", values_from = "value") %>%
    mutate(diff_scn1 = scn0-scn1, diff_scn2 = scn0-scn2) -> all_scn
  all_scn %>% select(sample, age_group, compartment, estado, diff_scn1, diff_scn2) %>%
    pivot_longer(c(diff_scn1, diff_scn2), names_to = "scenario", values_to = "value") -> all_scn
  all_scn %>% group_by(sample, compartment, estado, scenario) %>% summarize(value = sum(value)) %>%
                                                    left_join(DistEt) %>% left_join(prop_reforco) %>% mutate(reduction = value/pop) -> all_scn
  all_scn %>% group_by(compartment, estado, scenario) %>% summarize(porcentagem_reforco = mean(porcentagem_reforco),
                                                                    mean_reduction = mean(reduction), 
                                                                    var_reduction = var(reduction),
                                                                    q025_reduction = quantile(reduction,0.025),
                                                                    q975_reduction = quantile(reduction,0.975)) -> all_scn
  all_scn %>% mutate(scenario = if_else(scenario == "diff_scn1","Ritmo atual (lento)","Ritmo ideal")) %>%
    mutate(scenario = factor(scenario, levels = c("Ritmo atual (lento)", "Ritmo ideal"))) %>%
    mutate(compartment = if_else(compartment == "D","Óbito","Hospitalização")) %>%
    mutate(compartment = factor(compartment, levels = c("Hospitalização","Óbito"))) -> all_scn
  return(all_scn)
}

casos_totais <- function(){
  estados <- c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
                "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
                "RS", "SC", "SE", "SP", "TO")
  all_scn <- c()
  for(uf in estados)
  {
    scn0 <- read_csv(paste0("output/", uf,"_scn0.csv"))
    scn1 <- read_csv(paste0("output/", uf,"_scn1.csv"))
    scn2 <- read_csv(paste0("output/", uf,"_scn2.csv"))
    scn0$scenario <- "scn0"
    scn1$scenario <- "scn1"
    scn2$scenario <- "scn2"
    scn0$estado <- uf
    scn1$estado <- uf
    scn2$estado <- uf
    all_scn <- rbind(all_scn, scn0, scn1,scn2)
  }
  all_scn %>% group_by(sample, compartment, age_group, scenario) %>% 
    summarize(value = sum(value)) -> all_scn
  all_scn %>% mutate(age_group = case_when(age_group == 1 ~ "0-4",
                                           age_group == 2 ~ "5-11",
                                           age_group == 3 ~ "12-17",
                                           age_group == 4 ~ "18-29",
                                           age_group == 5 ~ "30-39",
                                           age_group == 6 ~ "40-49",
                                           age_group == 7 ~ "50-59",
                                           age_group == 8 ~ "60-69",
                                           age_group == 9 ~ "70-79",
                                           age_group == 10 ~ "80+")) %>%
    mutate(age_group = factor(age_group, levels = c("0-4","5-11","12-17","18-29","30-39",
                                                    "40-49","50-59","60-69","70-79","80+"))) -> all_scn
  
  all_scn %>% pivot_wider(names_from = "scenario", values_from = "value") %>%
    mutate(diff_scn1 = scn1-scn1, diff_scn2 = scn1-scn2) -> all_scn
  all_scn %>% select(sample, age_group, compartment, diff_scn1, diff_scn2) %>%
    pivot_longer(c(diff_scn1, diff_scn2), names_to = "scenario", values_to = "value") -> all_scn
  all_scn %>% group_by(sample, compartment, estado, scenario) %>% summarize(value = sum(value)) %>%
    left_join(DistEt) %>% left_join(prop_reforco) %>% mutate(reduction = value/pop) -> all_scn
  all_scn %>% group_by(compartment, estado, scenario) %>% summarize(porcentagem_reforco = mean(porcentagem_reforco),
                                                                    mean_reduction = mean(reduction), 
                                                                    var_reduction = var(reduction),
                                                                    q025_reduction = quantile(reduction,0.025),
                                                                    q975_reduction = quantile(reduction,0.975)) -> all_scn
  all_scn %>% mutate(scenario = if_else(scenario == "diff_scn1","Ritmo atual (lento)","Ritmo ideal")) %>%
    mutate(scenario = factor(scenario, levels = c("Ritmo atual (lento)", "Ritmo ideal"))) %>%
    mutate(compartment = if_else(compartment == "D","Óbito","Hospitalização")) %>%
    mutate(compartment = factor(compartment, levels = c("Hospitalização","Óbito"))) -> all_scn
  return(all_scn)
}

years_of_life_lost <- function(){
  estados <- c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
                "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
                "RS", "SC", "SE", "SP", "TO")
  all_scn <- c()
  for(uf in estados)
  {
    scn0 <- read_csv(paste0("output/", uf,"_scn0.csv"))
    scn1 <- read_csv(paste0("output/", uf,"_scn1.csv"))
    scn2 <- read_csv(paste0("output/", uf,"_scn2.csv"))
    scn0$scenario <- "scn0"
    scn1$scenario <- "scn1"
    scn2$scenario <- "scn2"
    scn0$estado <- uf
    scn1$estado <- uf
    scn2$estado <- uf
    all_scn <- rbind(all_scn, scn0, scn1,scn2)
  }
  all_scn %>% group_by(sample, compartment, age_group, scenario) %>% 
          summarize(value = sum(value)) -> all_scn
  all_scn %>% mutate(age_group = case_when(age_group == 1 ~ "0-4",
                                 age_group == 2 ~ "5-11",
                                 age_group == 3 ~ "12-17",
                                 age_group == 4 ~ "18-29",
                                 age_group == 5 ~ "30-39",
                                 age_group == 6 ~ "40-49",
                                 age_group == 7 ~ "50-59",
                                 age_group == 8 ~ "60-69",
                                 age_group == 9 ~ "70-79",
                                 age_group == 10 ~ "80+")) -> all_scn
  all_scn %>% mutate(media_fx = case_when(age_group ==  "0-4" ~ 4/2, 
                                          age_group ==  "5-11" ~ (5+11)/2,
                                          age_group ==  "12-17" ~ (12+17)/2,
                                          age_group ==  "18-29" ~ (18+29)/2,
                                          age_group ==  "30-39" ~ (30+39)/2,
                                          age_group ==  "40-49" ~ (40+49)/2,
                                          age_group ==  "50-59" ~ (50+59)/2,
                                          age_group ==  "60-69" ~ (60+69)/2,
                                          age_group ==  "70-79" ~ (70+79)/2)) -> all_scn
  
  all_scn %>% filter(compartment == "D") %>% mutate(YLL = value*(76-media_fx)) -> all_scn
  all_scn %>% group_by(age_group, scenario) %>% summarize(YLL_mean = mean(YLL, na.rm = T),
                                                          YLL_q025 = quantile(YLL,0.025, na.rm = T),
                                                          YLL_q975 = quantile(YLL,0.975, na.rm = T)) -> all_scn
  return(all_scn)
}
diff_years_of_life_lost <- function(){
  estados <- c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
                "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
                "RS", "SC", "SE", "SP", "TO")
  all_scn <- c()
  for(uf in estados)
  {
    scn0 <- read_csv(paste0("output/", uf,"_scn0.csv"))
    scn1 <- read_csv(paste0("output/", uf,"_scn1.csv"))
    scn2 <- read_csv(paste0("output/", uf,"_scn2.csv"))
    scn0$scenario <- "scn0"
    scn1$scenario <- "scn1"
    scn2$scenario <- "scn2"
    scn0$estado <- uf
    scn1$estado <- uf
    scn2$estado <- uf
    all_scn <- rbind(all_scn, scn0, scn1,scn2)
  }
  all_scn %>% group_by(sample, compartment, age_group, scenario) %>% 
    summarize(value = sum(value)) -> all_scn
  all_scn %>% pivot_wider(names_from = "scenario", values_from = "value") %>%
          mutate(diff_scn1 = scn0-scn1, diff_scn2 = scn0-scn2) -> all_scn
  all_scn %>% select(sample, age_group, compartment, diff_scn1, diff_scn2) %>%
    pivot_longer(c(diff_scn1, diff_scn2), names_to = "scenario", values_to = "value") -> all_scn
  all_scn %>% mutate(age_group = case_when(age_group == 1 ~ "0-4",
                                           age_group == 2 ~ "5-11",
                                           age_group == 3 ~ "12-17",
                                           age_group == 4 ~ "18-29",
                                           age_group == 5 ~ "30-39",
                                           age_group == 6 ~ "40-49",
                                           age_group == 7 ~ "50-59",
                                           age_group == 8 ~ "60-69",
                                           age_group == 9 ~ "70-79",
                                           age_group == 10 ~ "80+")) %>%
    mutate(age_group = factor(age_group, levels = c("0-4","5-11","12-17","18-29","30-39",
                                                    "40-49","50-59","60-69","70-79","80+"))) -> all_scn
  all_scn %>% mutate(media_fx = case_when(age_group ==  "0-4" ~ 4/2, 
                                          age_group ==  "5-11" ~ (5+11)/2,
                                          age_group ==  "12-17" ~ (12+17)/2,
                                          age_group ==  "18-29" ~ (18+29)/2,
                                          age_group ==  "30-39" ~ (30+39)/2,
                                          age_group ==  "40-49" ~ (40+49)/2,
                                          age_group ==  "50-59" ~ (50+59)/2,
                                          age_group ==  "60-69" ~ (60+69)/2,
                                          age_group ==  "70-79" ~ (70+79)/2)) -> all_scn
  
  all_scn %>% filter(compartment == "D") %>% mutate(YLL = value*(76-media_fx)) -> all_scn
  all_scn %>% group_by(age_group, scenario) %>% summarize(YLL_mean = mean(YLL, na.rm = T),
                                                          YLL_q025 = quantile(YLL,0.025, na.rm = T),
                                                          YLL_q975 = quantile(YLL,0.975, na.rm = T)) -> all_scn
  return(all_scn)
}

reducao_custos <- function(){
  all_scn <- c()
  for(uf in estados)
  {
    scn0 <- read_csv(paste0("output/", uf,"_scn0.csv"))
    scn1 <- read_csv(paste0("output/", uf,"_scn1.csv"))
    scn2 <- read_csv(paste0("output/", uf,"_scn2.csv"))
    scn0$scenario <- "scn0"
    scn1$scenario <- "scn1"
    scn2$scenario <- "scn2"
    scn0$estado <- uf
    scn1$estado <- uf
    scn2$estado <- uf
    all_scn <- rbind(all_scn, scn0, scn1,scn2)
  }
  all_scn %>% group_by(sample, compartment, age_group, scenario) %>% 
    summarize(value = sum(value)) -> all_scn
  all_scn %>% mutate(age_group = case_when(age_group == 1 ~ "0-4",
                                           age_group == 2 ~ "5-11",
                                           age_group == 3 ~ "12-17",
                                           age_group == 4 ~ "18-29",
                                           age_group == 5 ~ "30-39",
                                           age_group == 6 ~ "40-49",
                                           age_group == 7 ~ "50-59",
                                           age_group == 8 ~ "60-69",
                                           age_group == 9 ~ "70-79",
                                           age_group == 10 ~ "80+")) %>%
    mutate(age_group = factor(age_group, levels = c("0-4","5-11","12-17","18-29","30-39",
                                                    "40-49","50-59","60-69","70-79","80+"))) -> all_scn
  
  all_scn %>% pivot_wider(names_from = "scenario", values_from = "value") %>%
    mutate(diff_scn1 = scn1-scn1, diff_scn2 = scn1-scn2) -> all_scn
  all_scn %>% select(sample, age_group, compartment, diff_scn1, diff_scn2) %>%
    pivot_longer(c(diff_scn1, diff_scn2), names_to = "scenario", values_to = "value") -> all_scn
  all_scn %>% filter(compartment == "C") %>% mutate(reduced_cost = 10379.51*value) -> all_scn
  return(all_scn)
}





