library(dplyr)

group_by_age <- function(solution){
  solution <- solution %>% mutate(classe = case_when(str_detect(classe,"Su")~"Su",
                                                     str_detect(classe,"Eu")~"Eu",
                                                     str_detect(classe,"Iu")~"Iu",
                                                     str_detect(classe,"Au")~"Au",
                                                     str_detect(classe,"Hu")~"Hu",
                                                     str_detect(classe,"Ru")~"Ru",
                                                     str_detect(classe,"Du")~"Du",
                                                     str_detect(classe,"SvA")~"SvA",
                                                     str_detect(classe,"EvA")~"EvA",
                                                     str_detect(classe,"IvA")~"IvA",
                                                     str_detect(classe,"AvA")~"AvA",
                                                     str_detect(classe,"HvA")~"HvA",
                                                     str_detect(classe,"RvA")~"RvA",
                                                     str_detect(classe,"DvA")~"DvA",
                                                     str_detect(classe,"SvP")~"SvP",
                                                     str_detect(classe,"EvP")~"EvP",
                                                     str_detect(classe,"IvP")~"IvP",
                                                     str_detect(classe,"AvP")~"AvP",
                                                     str_detect(classe,"HvP")~"HvP",
                                                     str_detect(classe,"RvP")~"RvP",
                                                     str_detect(classe,"DvP")~"DvP",
                                                     str_detect(classe,"SvC")~"SvC",
                                                     str_detect(classe,"EvC")~"EvC",
                                                     str_detect(classe,"IvC")~"IvC",
                                                     str_detect(classe,"AvC")~"AvC",
                                                     str_detect(classe,"HvC")~"HvC",
                                                     str_detect(classe,"RvC")~"RvC",
                                                     str_detect(classe,"DvC")~"DvC",
                                                     str_detect(classe,"SwA")~"SwA",
                                                     str_detect(classe,"EwA")~"EwA",
                                                     str_detect(classe,"IwA")~"IwA",
                                                     str_detect(classe,"AwA")~"AwA",
                                                     str_detect(classe,"HwA")~"HwA",
                                                     str_detect(classe,"RwA")~"RwA",
                                                     str_detect(classe,"DwA")~"DwA",
                                                     str_detect(classe,"SwP")~"SwP",
                                                     str_detect(classe,"EwP")~"EwP",
                                                     str_detect(classe,"IwP")~"IwP",
                                                     str_detect(classe,"AwP")~"AwP",
                                                     str_detect(classe,"HwP")~"HwP",
                                                     str_detect(classe,"RwP")~"RwP",
                                                     str_detect(classe,"DwP")~"DwP",
                                                     str_detect(classe,"SwC")~"SwC",
                                                     str_detect(classe,"EwC")~"EwC",
                                                     str_detect(classe,"IwC")~"IwC",
                                                     str_detect(classe,"AwC")~"AwC",
                                                     str_detect(classe,"HwC")~"HwC",
                                                     str_detect(classe,"RwC")~"RwC",
                                                     str_detect(classe,"DwC")~"DwC",
                                                     str_detect(classe,"SbP")~"SbP",
                                                     str_detect(classe,"EbP")~"EbP",
                                                     str_detect(classe,"IbP")~"IbP",
                                                     str_detect(classe,"AbP")~"AbP",
                                                     str_detect(classe,"HbP")~"HbP",
                                                     str_detect(classe,"RbP")~"RbP",
                                                     str_detect(classe,"DbP")~"DbP")
  ) %>% group_by(t, classe) %>% summarize(value = sum(value)) %>% ungroup()
  solution$classe <- factor(solution$classe, levels =  c("Su","Eu","Iu","Au","Hu","Ru","Du",
                                                         paste0(c("Sv","Ev","Iv","Av","Hv","Rv","Dv"),rep(c("A","P","C"),each =9)),
                                                         paste0(c("Sw","Ew","Iw","Aw","Hw","Rw","Dw"),rep(c("A","P","C"),each =9))))
  
  return(solution)
}

group_by_dose <- function(solution){
  ###group by dose
  solution <- solution %>% mutate(classe = case_when(str_detect(classe,"Su")~"Su",
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
                                                     str_detect(classe,"Dw")~"Dw",
                                                     str_detect(classe,"Sb")~"Sb",
                                                     str_detect(classe,"Eb")~"Eb",
                                                     str_detect(classe,"Ib")~"Ib",
                                                     str_detect(classe,"Ab")~"Ab",
                                                     str_detect(classe,"Hb")~"Hb",
                                                     str_detect(classe,"Rb")~"Rb",
                                                     str_detect(classe,"Db")~"Db")
  ) %>% group_by(t, classe) %>% summarize(value = sum(value)) %>% ungroup()
  solution$classe <- factor(solution$classe, levels = c("Su","Eu","Iu","Au","Hu","Ru","Du",
                                                        "Sv","Ev","Iv","Av","Hv","Rv","Dv",
                                                        "Sw","Ew","Iw","Aw","Hw","Rw","Dw",
                                                        "Sb","Eb","Ib","Ab","Hb","Rb","Db"))
  return(solution)
}

wide_solution <- function(solution) {
  ### disaggregate column classe by compartment, age, variant, and vaccinal status
  solution %>% extract(col=classe,
                       into=c("compartment", "vac", "variant", "age_group"),
                       regex="^([A-Z])([A-Za-z]+)\\.?([OD]?)([0-9]+)$") %>%
               as.data.frame
}
