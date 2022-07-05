library(tidyverse)

organiza_cobertura_vacinal <- function(parameters, cobertura){
  Y2 <- parameters$dist.etaria
  Y2 <- with(parameters,{
    ############### cobertura pfizer ###########################
    dados_dose <- cobertura %>% filter(vacina_D1 == "Pfizer") %>%
      mutate(agegroup = if_else(agegroup >= 10,10,agegroup)) %>%
      group_by(vacina_D1,dose,agegroup) %>% 
      summarize(n = sum(n)) %>% ungroup() %>%
      pivot_wider(names_from = "dose", values_from = "n")
    Y2[Suindex] <- Y2[Suindex] - dados_dose$D1
    Y2[SvPindex] <- Y2[SvPindex] + dados_dose$D1
    Y2[Suindex] <- Y2[Suindex] - dados_dose$D2
    Y2[SwPindex] <- Y2[SwPindex] + dados_dose$D2
    Y2[Suindex] <- Y2[Suindex] - dados_dose$R
    Y2[SbPindex] <- Y2[SbPindex] + dados_dose$R
    
    ############### cobertura astrazeneca ###########################
    dados_dose <- cobertura %>% filter(vacina_D1 == "AZ") %>%
      mutate(agegroup = if_else(agegroup >= 10,10,agegroup)) %>%
      group_by(vacina_D1,dose,agegroup) %>% 
      summarize(n = sum(n)) %>% ungroup() %>%
      pivot_wider(names_from = "dose", values_from = "n")
    Y2[Suindex] <- Y2[Suindex] - dados_dose$D1
    Y2[SvAindex] <- Y2[SvAindex] + dados_dose$D1
    Y2[Suindex] <- Y2[Suindex] - dados_dose$D2
    Y2[SwAindex] <- Y2[SwAindex] + dados_dose$D2
    Y2[Suindex] <- Y2[Suindex] - dados_dose$R
    Y2[SbPindex] <- Y2[SbPindex] + dados_dose$R ##reforço entra em pfizer
    
    ############### cobertura coronavac ###########################
    dados_dose <- cobertura %>% filter(vacina_D1 == "Coronavac") %>%
      mutate(agegroup = if_else(agegroup >= 10,10,agegroup)) %>%
      group_by(vacina_D1,dose,agegroup) %>% 
      summarize(n = sum(n)) %>% ungroup() %>%
      pivot_wider(names_from = "dose", values_from = "n")
    Y2[Suindex] <- Y2[Suindex] - dados_dose$D1
    Y2[SvCindex] <- Y2[SvCindex] + dados_dose$D1
    Y2[Suindex] <- Y2[Suindex] - dados_dose$D2
    Y2[SwCindex] <- Y2[SwCindex] + dados_dose$D2
    Y2[Suindex] <- Y2[Suindex] - dados_dose$R
    Y2[SbPindex] <- Y2[SbPindex] + dados_dose$R ##reforço entra em pfizer
    
    #################################################################
    #correção populações negativas, se menor que 1% da faixa estimada, corrige para 1%
    Y2[Suindex] <- pmax(Y2[Suindex],0.01*Y[Suindex])
    Y2
  })
  parameters$init.condition <- Y2
  return(parameters)
}
