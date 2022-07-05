library(tidyverse)
library(zeallot)
library(gdata)
source("functions/vac_params_omicron.R")
source("functions/solve_model_reinfec.R")
source("functions/func_utils.R")
t <- 1:145
keep(parameters,t,list = lsf.str(),sure = T) #### limpando o environment, tomar cuidado
solution <- solve_system(parameters2,1:90)
###group by age
solution2 <- group_by_dose(solution)
solution2 %>% 
  filter(classe %in% c("Iu","Iv","Iw","Ib")) %>%
  ggplot(aes(x = t,y = value, color = classe)) + geom_line()+theme_bw(base_size = 16)+
  viridis::scale_color_viridis(discrete = T)
