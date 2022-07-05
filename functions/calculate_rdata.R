library("dplyr")
# [1] "Brasil"              "Região Norte"        "Região Nordeste"    
# [4] "Região Sudeste"      "Região Sul"          "Região Centro-oeste"
# [7] "Rondônia"            "Acre"                "Amazonas"           
# [10] "Roraima"             "Pará"                "Amapá"              
# [13] "Tocantins"           "Maranhão"            "Piauí"              
# [16] "Ceará"               "Rio Grande do Norte" "Paraíba"            
# [19] "Pernambuco"          "Alagoas"             "Sergipe"            
# [22] "Bahia"               "Minas Gerais"        "Espírito Santo"     
# [25] "Rio de Janeiro"      "São Paulo"           "Paraná"             
# [28] "Santa Catarina"      "Rio Grande do Sul"   "Mato Grosso do Sul" 
# [31] "Mato Grosso"         "Goiás"               "Distrito Federal"   
# [34] "Regional Norte"      "Regional Leste"      "Regional Centro"    
# [37] "Regional Sul"   

fitP.exp <- function(X, Y, only.coef = FALSE){
    fit <- try(glm(Y ~ X, family = poisson))
    if(only.coef){
            ci <- try(confint(fit))
            if(any(class(fit)=="try-error")||any(class(ci)=="try-error"))
                results <- rep(NA,6)
            else
                results <- c(coef(fit),ci[1,], ci[2,])
            names(results) <- c("intercept", "coef", "int.low", "int.upp", "coef.low", "coef.upp")
    }
    if(!only.coef){
        if(class(fit)=="try-error")
            results <- NA
        else
            results  <-  fit
    }
    return(results)
}

infogripe <- read.csv("DATA/serie_temporal_com_estimativas_recentes_sem_filtro_febre.csv", sep=";", head=TRUE)
infogripe$Casos.semanais.reportados.até.a.última.atualização <- as.numeric(gsub(",", ".", infogripe$Casos.semanais.reportados.até.a.última.atualização))
infogripe$casos.estimados <- as.numeric(gsub(",", ".", infogripe$casos.estimados))

estimativa_r <- function(estado)
{
  est <- infogripe %>% filter(Unidade.da.Federação==estado, Ano.epidemiológico >= 2020,escala=="casos", dado=="srag") %>% 
      mutate(Ano.epidemiológico = replace(Ano.epidemiológico, Ano.epidemiológico == 2020, 0)) %>%
      mutate(Ano.epidemiológico = replace(Ano.epidemiológico, Ano.epidemiológico == 2021, 53)) %>%
      mutate(Ano.epidemiológico = replace(Ano.epidemiológico, Ano.epidemiológico == 2022, 105))

  est.o <- est %>% filter(Ano.epidemiológico + Semana.epidemiológica >= 100)

  est.o$semanas <- 1:length(est.o$Semana.epidemiológica)
  
  glmp.exp.est.o <- fitP.exp(est.o$semanas, est.o$Casos.semanais.reportados.até.a.última.atualização)
  return(glmp.exp.est.o[2]/7)
}


#https://raw.githubusercontent.com/covid19br/covid19br.github.io/master/dados/municipios/SP/Sao_Paulo/tabelas_nowcasting_para_grafico/nowcasting_diario_covid_2022_01_26.csv
estimativa_r_obs_covid <- function()
{
  obsraw <- read.csv("DATA/nowcasting_diario_covid_2022_01_26.csv")
  obsraw$data <- as.Date(obsraw$data)
  obs <- obsraw %>% filter(data>=as.Date("2021-12-19") & data<=as.Date("2022-01-09"))
  obs$data <- c(1,2,3,4)
  glmp.exp.est.o.obs <- fitP.exp(obs$data, obs$estimate.merged)
  est <- coef(glmp.exp.est.o.obs)[2]/7
  sdest <- summary(glmp.exp.est.o.obs)$coefficients[2,2]/7
  return(c(est, sdest))
}

#https://raw.githubusercontent.com/covid19br/covid19br.github.io/master/dados/municipios/SP/Sao_Paulo/tabelas_nowcasting_para_grafico/nowcasting_diario_srag_2022_01_26.csv
estimativa_r_obs <- function(estado, tipo, ind)
{
  repo <- "https://raw.githubusercontent.com/covid19br/covid19br.github.io/master/dados/estado/"
  arquivo <- paste0("/tabelas_nowcasting_para_grafico/nowcasting_diario_", tipo,"_2022_01_31.csv")
  obsraw <- read.csv(paste0(repo, estado, arquivo))
  obsraw$data <- as.Date(obsraw$data)
  #obs <- obsraw %>% filter(data>=as.Date("2021-12-19") & data<=as.Date("2022-01-09"))
  obsdata <- c(1,2,3,4)
  
  if(obsraw$data[length(obsraw$data)]==as.Date("2022-01-30"))
  {
    ind = ind + 1
  }
  glmp.exp.est.o.obs <- fitP.exp(obsdata, obsraw$estimate.merged[(length(obsraw$data)-ind-3):(length(obsraw$data)-ind)])
  est <- coef(glmp.exp.est.o.obs)[2]/7
  sdest <- summary(glmp.exp.est.o.obs)$coefficients[2,2]/7
  return(c(est, sdest))
}

# UFs <- c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA","MG","MS","MT","PA","PB","PE","PI","PR","RJ","RN","RO","RR","RS","SC","SE","SP","TO")
# rest2 <- data.frame()
# for(i in 0:4)
# {
#   rsrag <- c()
#   rcovid <- c()
#   sdrsrag <- c()
#   sdrcovid <- c()
#   for(uf in UFs)
#   {
#       estsrag <- estimativa_r_obs(uf, "srag", i)
#       rsrag <- c(rsrag, estsrag[1])
#       sdrsrag <- c(sdrsrag, estsrag[2])
#       estcovid <- estimativa_r_obs(uf, "covid", i)
#       rcovid <- c(rcovid, estcovid[1])
#       sdrcovid <- c(sdrcovid, estcovid[2])
#   }
# 
#   rest <- data.frame(srag = rsrag, covid = rcovid, sdsrag = sdrsrag, sdcovid = sdrcovid, UF = UFs, janela = rep(i, length(UFs)))
#   rest2 <- rbind(rest2, rest)
# }
# 
# write.csv(rest2, file="DATA/r_obs.csv", row.names=F)

# for(uf in UFs)
# {
#   repo <- "https://raw.githubusercontent.com/covid19br/covid19br.github.io/master/dados/estado/"
#   arquivo <- paste0("/tabelas_nowcasting_para_grafico/nowcasting_diario_covid_2022_01_31.csv")
#   obsraw <- read.csv(paste0(repo, uf, arquivo))
#   obsraw$data <- as.Date(obsraw$data)
#   print(uf)
#   print(obsraw$data[length(obsraw$data)])
#   print("")
# 
# }

# rsrag <- c()
# rcovid <- c()
# df <- data.frame()
# UFs <- c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA","MG","MS","MT","PA","PB","PE","PI","PR","RJ","RN","RO","RR","RS","SC","SE","SP","TO")
# for(uf in UFs)
# {
#   repo <- "https://raw.githubusercontent.com/covid19br/covid19br.github.io/master/dados/estado/"
#   arqsrag <- paste0("/tabelas_nowcasting_para_grafico/nowcasting_diario_srag_2022_01_31.csv")
#   arqcovid <- paste0("/tabelas_nowcasting_para_grafico/nowcasting_diario_covid_2022_01_31.csv")
#   obssrag <- read.csv(paste0(repo, uf, arqsrag))
#   obscovid <- read.csv(paste0(repo, uf, arqcovid))
#   # print(uf)
#   # print(paste("SRAG: ",obssrag$data[length(obssrag$data)]))
#   # print(paste("COVID: ",obscovid$data[length(obscovid$data)]))
#   # print("")
#   obssrag$data <- as.Date(obssrag$data)
#   obscovid$data <- as.Date(obscovid$data)
#   obssrag <- obssrag %>% filter(data>=as.Date("2021-11-21") & data<=as.Date("2022-01-23"))
#   obscovid <- obscovid %>% filter(data>=as.Date("2021-11-21") & data<=as.Date("2022-01-23"))
#   obssrag$n.casos <- obssrag$estimate <- obssrag$estimate.merged.smooth <- obssrag$lower.merged.pred <- obssrag$upper.merged.pred <- NULL
#   obscovid$n.casos <- obscovid$estimate <- obscovid$estimate.merged.smooth <- obscovid$lower.merged.pred <- obscovid$upper.merged.pred <- NULL
#   # if(uf=="AC")
#   # {
#   #     obscovid <- rbind(obscovid, c(NA, NA))
#   #     obscovid <- rbind(obscovid, c(NA, NA))
#   # }
#   obssrag$semana <- 1:10
#   obscovid$semana <- 1:10
#   obssrag$dado <- c("srag")
#   obscovid$dado <- c("covid")
#   obssrag$UF <- uf
#   obscovid$UF <- uf
#   dfuf <- rbind(obssrag, obscovid)
#   df <- rbind(df, dfuf)
# }
# write.csv(df, file="DATA/nowcasting_10-semanas.csv", row.names=F)