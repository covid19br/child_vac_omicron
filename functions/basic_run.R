library(tidyverse)
library(zeallot)
source("functions/solve_model.R")
source("functions/vac_params.R")
solution <- solve_system(Y, t, parameters)
solution2 <- solution %>% mutate(classe = case_when(str_detect(classe,"Su")~"Su",
                                                    str_detect(classe,"Eu")~"Eu",
                                                    str_detect(classe,"Iu")~"Iu",
                                                    str_detect(classe,"Au")~"Au",
                                                    str_detect(classe,"Hu")~"Hu",
                                                    str_detect(classe,"Ru")~"Ru",
                                                    str_detect(classe,"Du")~"Du",
                                                    str_detect(classe,"Sv")~"Sv",
                                                    str_detect(classe,"Ev")~"Ev",
                                                    str_detect(classe,"Iv")~"Iv",
                                                    str_detect(classe,"Av")~"Av",
                                                    str_detect(classe,"Hv")~"Hv",
                                                    str_detect(classe,"Rv")~"Rv",
                                                    str_detect(classe,"Dv")~"Dv")
) %>% group_by(t, classe) %>% summarize(value = sum(value)) %>% ungroup()


solution2$classe <- factor(solution2$classe, levels = c("Su","Eu","Iu","Au","Hu","Ru","Du","Sv","Ev","Iv","Av","Hv","Rv","Dv"))
ggplot(solution2, aes(x = t,y = value, color = classe)) + geom_line()+theme_bw()+theme(text = element_text(size = 16))
