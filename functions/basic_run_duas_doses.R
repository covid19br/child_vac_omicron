library(tidyverse)
library(zeallot)
source("functions/solve_model_duas_doses.R")
source("functions/vac_params_duas_doses.R")
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
                                                    str_detect(classe,"Dv")~"Dv",
                                                    str_detect(classe,"Sw")~"Sw",
                                                    str_detect(classe,"Ew")~"Ew",
                                                    str_detect(classe,"Iw")~"Iw",
                                                    str_detect(classe,"Aw")~"Aw",
                                                    str_detect(classe,"Hw")~"Hw",
                                                    str_detect(classe,"Rw")~"Rw",
                                                    str_detect(classe,"Dw")~"Dw")
) %>% group_by(t, classe) %>% summarize(value = sum(value)) %>% ungroup()


solution2$classe <- factor(solution2$classe, levels = c("Su","Eu","Iu","Au","Hu","Ru","Du",
                                                        "Sv","Ev","Iv","Av","Hv","Rv","Dv",
                                                        "Sw","Ew","Iw","Aw","Hw","Rw","Dw"))
ggplot(solution2, aes(x = t,y = value, color = classe)) + geom_line()+theme_bw()+theme(text = element_text(size = 16))
