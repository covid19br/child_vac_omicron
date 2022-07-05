library(dplyr)
library(ggplot2 )

deaths_by_age <- function()
{ # Plota as mortes por faixa etária ao final da simulação. Provavelmente vai
  #deixar de funcionar quando mudar o formato do output.
  estados = c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
              "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
              "RS", "SC", "SE", "SP", "TO")
  for(uf in estados)
  {
    scn0 <- read.csv(paste0("output/", uf,"_scn0.csv"))
    scn1 <- read.csv(paste0("output/", uf,"_scn1.csv"))
    scn0$scenario <- 0
    scn1$scenario <- 1
    dif <- rbind(scn0, scn1)
    dif <- dif %>% mutate(vacinado = ifelse (dose=="u", "u", "v")) %>% ungroup()
    dif$dose <- NULL
    for(d in c("u", "v"))
    {
      dif %>% filter(vacinado==d) %>% group_by(scenario, age_group, sample) %>% 
      summarize(value = sum(value)) %>% ungroup() %>% group_by(age_group, sample) %>%
      summarize(value = value[scenario==0] - value[scenario==1]) %>% ungroup() %>%
      ggplot(aes(x=as.factor(age_group), y=value)) + geom_boxplot() +
      labs(x = "Age group", y = "deaths", title = paste("Mortes",uf, d)) -> plots
      ggsave(paste0( "figures/", uf,"_", d,"_Ddiff.png" ), plot=plots)
    }
  }
}

plot_brasil <- function(){
  estados = c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
               "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
               "RS", "SC", "SE", "SP", "TO")
  # estados <- "SP"
  all_scn <- c()
  for(uf in estados)
  {
    scn0 <- read_csv(paste0("output/", uf,"_scn0_X_v2_SA.csv"))
    scn1 <- read_csv(paste0("output/", uf,"_scn1_X_v2_SA.csv"))
    scn2 <- read_csv(paste0("output/", uf,"_scn2_X_v2_SA.csv"))
    scn0$scenario <- "scn0"
    scn1$scenario <- "scn1"
    scn2$scenario <- "scn2"
    scn0$estado <- uf
    scn1$estado <- uf
    scn2$estado <- uf
    all_scn <- rbind(all_scn, scn0, scn1,scn2)
  }
  all_scn %>% group_by(sample, age_group, compartment, scenario) %>%
    summarize(value = sum(value)) -> all_scn
  all_scn %>% pivot_wider(names_from = scenario, values_from = value) -> all_scn
  all_scn %>% mutate(diff_scn1 = scn0-scn1, diff_scn2 = scn0 - scn2) %>%
              mutate(age_group = case_when(age_group == 1 ~ "0-4",
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
                                                    "40-49","50-59","60-69","70-79","80+")))%>%
    mutate(compartment = case_when(compartment == "D"~ "Deaths",
                                   compartment == "C" ~ "Hospitalizations",
                                   compartment == "X" ~ "Cases")) %>%
    mutate(compartment = factor(compartment, levels = c("Cases","Hospitalizations","Deaths"))) -> all_scn
  all_scn %>% select(sample, age_group, compartment, diff_scn1, diff_scn2) %>%
    pivot_longer(c(diff_scn1, diff_scn2), names_to = "scenario", values_to = "value") -> all_scn
  all_scn %>% mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
    mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
    filter(compartment %in% c("Hospitalizations","Deaths")) %>%
    ggplot(aes(x = age_group, y = value,fill = scenario,color = scenario, group = interaction(age_group,scenario))) + 
    facet_wrap(~compartment, scales = "free") + 
    # geom_boxplot(outlier.size = 0.5, outlier.color = "black") + 
    theme_bw(base_size = 16) + theme(legend.position = "bottom", 
                                     legend.direction = "horizontal", 
                                     axis.text.x = element_text(angle = 90, vjust = 0.5,hjust= 1)) +
    geom_violin(trim = T,scale = "width", position = "dodge", draw_quantiles = c(0.025, 0.5, 0.95))+
    # geom_boxplot(width=0.1, alpha=0.2)+ 
    labs(x = "Age Group", y = "Outcome averted by vaccination", color = "Scenario", fill = "Scenario") +   scale_color_discrete(l = 50) -> p
  return(p)
}
plot_brasil_pt <- function(){
  estados = c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
               "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
               "RS", "SC", "SE", "SP", "TO")
  # estados <- "SP"
  all_scn <- c()
  for(uf in estados)
  {
    scn0 <- read_csv(paste0("output/", uf,"_scn0_X_v2_SA.csv"))
    scn1 <- read_csv(paste0("output/", uf,"_scn1_X_v2_SA.csv"))
    scn2 <- read_csv(paste0("output/", uf,"_scn2_X_v2_SA.csv"))
    scn0$scenario <- "scn0"
    scn1$scenario <- "scn1"
    scn2$scenario <- "scn2"
    scn0$estado <- uf
    scn1$estado <- uf
    scn2$estado <- uf
    all_scn <- rbind(all_scn, scn0, scn1,scn2)
  }
  all_scn %>% group_by(sample, age_group, compartment, scenario) %>%
    summarize(value = sum(value)) -> all_scn
  all_scn %>% pivot_wider(names_from = scenario, values_from = value) -> all_scn
  all_scn %>% mutate(diff_scn1 = scn0-scn1, diff_scn2 = scn0 - scn2) %>%
    mutate(age_group = case_when(age_group == 1 ~ "0-4",
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
                                                    "40-49","50-59","60-69","70-79","80+")))%>%
    mutate(compartment = case_when(compartment == "D"~ "Óbitos",
                                   compartment == "C" ~ "Hospitalizações",
                                   compartment == "X" ~ "Casos")) %>%
    mutate(compartment = factor(compartment, levels = c("Casos","Hospitalizações","Óbitos"))) -> all_scn
  all_scn %>% select(sample, age_group, compartment, diff_scn1, diff_scn2) %>%
    pivot_longer(c(diff_scn1, diff_scn2), names_to = "scenario", values_to = "value") -> all_scn
  all_scn %>% mutate(scenario = if_else(scenario == "diff_scn1","Ritmo lento (atual)","Ritmo ideal")) %>%
    mutate(scenario = factor(scenario, levels = c("Ritmo lento (atual)", "Ritmo ideal"))) %>%
    filter(compartment %in% c("Hospitalizações","Óbitos")) %>%
    ggplot(aes(x = age_group, y = value,fill = scenario,color = scenario, group = interaction(age_group,scenario))) + 
    facet_wrap(~compartment, scales = "free") + 
    geom_boxplot(outlier.size = 0.5, outlier.color = "black") +
    theme_bw(base_size = 16) + theme(legend.position = "bottom", 
                                     legend.direction = "horizontal", 
                                     axis.text.x = element_text(angle = 90, vjust = 0.5,hjust= 1)) +
    # geom_violin(trim = T,scale = "width", position = "dodge", draw_quantiles = c(0.025, 0.5, 0.95))+
    # geom_boxplot(width=0.1, alpha=0.2)+ 
    labs(x = "Faixa etária", y = "Evento evitado por vacinação", color = "Cenário", fill = "Cenário") +   scale_color_discrete(l = 50) -> p
  return(p)
}
summary_brasil <- function(){
  estados = c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
               "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
               "RS", "SC", "SE", "SP", "TO")
  all_scn <- c()
  for(uf in estados)
  {
    scn0 <- read_csv(paste0("output/", uf,"_scn0_X_v2.csv"))
    scn1 <- read_csv(paste0("output/", uf,"_scn1_X_v2.csv"))
    scn2 <- read_csv(paste0("output/", uf,"_scn2_X_v2.csv"))
    scn0$scenario <- "scn0"
    scn1$scenario <- "scn1"
    scn2$scenario <- "scn2"
    scn0$estado <- uf
    scn1$estado <- uf
    scn2$estado <- uf
    all_scn <- rbind(all_scn, scn0, scn1,scn2)
  }
  all_scn %>% group_by(sample, age_group, compartment, scenario) %>%
    summarize(value = sum(value)) -> all_scn
  # all_scn %>% pivot_wider(names_from = scenario, values_from = value) -> all_scn
  all_scn %>%
    mutate(age_group = case_when(age_group == 1 ~ "0-4",
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
                                                    "40-49","50-59","60-69","70-79","80+")))%>%
    mutate(compartment = case_when(compartment == "D"~ "Death",
                                   compartment == "C" ~ "Hospitalization",
                                   compartment == "X" ~ "Cases")) %>%
    mutate(compartment = factor(compartment, levels = c("Cases","Hospitalization","Death")))  -> all_scn
  all_scn %>% group_by(sample,age_group,  compartment, scenario) %>% summarize(value = sum(value))%>%
          group_by(age_group, compartment, scenario) %>%
    summarize(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975)) -> all_scn
  # 
  # all_scn %>% select(sample, age_group, compartment, diff_scn1, diff_scn2) %>%
  #   pivot_longer(c(diff_scn1, diff_scn2), names_to = "scenario", values_to = "value") -> all_scn
  # all_scn %>% mutate(scenario = if_else(scenario == "diff_scn1","Realized Vaccination","Ideal Vaccination")) %>%
  #   mutate(scenario = factor(scenario, levels = c("Realized Vaccination", "Ideal Vaccination"))) %>%
  #   ggplot(aes(x = age_group, y = value,fill = scenario, group = interaction(age_group,scenario))) + 
  #   facet_wrap(~compartment, scales = "free") + geom_boxplot() + 
  #   labs(x = "Age group", y = "Outcome averted by vaccination", title = "Brazil", fill = "Scenario") + 
  #   theme_bw(base_size = 16)+
  #   theme(legend.position = "bottom",legend.direction = "horizontal") -> p
  return(all_scn)
  
}


diff_brasil <- function(){

    estados = c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
                 "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
                 "RS", "SC", "SE", "SP", "TO")
    all_scn <- c()
    for(uf in estados)
    {
      scn0 <- read_csv(paste0("output/", uf,"_scn0_X_v2.csv"))
      scn1 <- read_csv(paste0("output/", uf,"_scn1_X_v2.csv"))
      scn2 <- read_csv(paste0("output/", uf,"_scn2_X_v2.csv"))
      scn0$scenario <- "scn0"
      scn1$scenario <- "scn1"
      scn2$scenario <- "scn2"
      scn0$estado <- uf
      scn1$estado <- uf
      scn2$estado <- uf
      all_scn <- rbind(all_scn, scn0, scn1,scn2)
    }
    all_scn %>% group_by(sample, compartment, scenario,age_group, estado) %>%
      summarize(value = sum(value)) -> all_scn
    all_scn %>% pivot_wider(names_from = scenario, values_from = value) -> all_scn
    all_scn %>% mutate(diff_scn1 = (scn0-scn1)/scn0, diff_scn2 = (scn0 - scn2)/scn0) %>%
      mutate(age_group = case_when(age_group == 1 ~ "0-4",
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
                                                      "40-49","50-59","60-69","70-79","80+")))%>%
      mutate(compartment = case_when(compartment == "D"~ "Death",
                                     compartment == "C" ~ "Hospitalization",
                                     compartment == "X" ~ "Cases")) %>%
      mutate(compartment = factor(compartment, levels = c("Cases","Hospitalization","Death"))) -> all_scn
    all_scn %>% select(sample, compartment, diff_scn1, diff_scn2, age_group, estado) %>%
      pivot_longer(c(diff_scn1, diff_scn2), names_to = "scenario", values_to = "value") -> all_scn
    all_scn %>% mutate(scenario = if_else(scenario == "diff_scn1","Realized Vaccination","Ideal Vaccination")) %>%
      mutate(scenario = factor(scenario, levels = c("Realized Vaccination", "Ideal Vaccination"))) -> all_scn
    all_scn %>% group_by( age_group, compartment, scenario, estado) %>% 
      summarize(mean = mean(value), q025 = quantile(value, 0.025), q975 = quantile(value, 0.975)) -> all_scn
    return(all_scn)
}
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
                                                  "40-49","50-59","60-69","70-79","80+")))%>%
ggplot(aes(x = age_group, y = value,fill = scenario, group = interaction(age_group,scenario))) + 
  facet_wrap(~compartment, scales = "free") + geom_boxplot() + 
  labs(x = "Age group", y = "Outcome averted by vaccination", title = "São Paulo") + 
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom",legend.direction = "horizontal")

