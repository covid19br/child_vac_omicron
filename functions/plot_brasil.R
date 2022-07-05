library(tidyverse)
install.packages("geobr")
devtools::install_github("ipeaGIT/geobr", subdir = "r-package")
library(geobr)

##tem que ler o plots.R e executar o summary
##depois fazer states <- left_join(states, data, by = c("abbrev_state" = "estado"))
states <- read_state()

states %>% mutate(compartment = factor(compartment, levels = c("Cases", "Hospitalization", "Death"), labels = c("Casos", "Hospitalizações","Óbitos")),
                  scenario = factor(scenario, levels = c("Realized Vaccination", "Ideal Vaccination"), labels = c("Ritmo atual (lento)", "Ritmo ideal"))) %>%
  filter(age_group == "5-11", compartment == "Óbitos") %>% 
  ggplot(aes(label = abbrev_state, fill = mean)) + 
  facet_wrap(~scenario)+geom_sf(lwd = 0.2) + theme_minimal(base_size = 16) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + scale_fill_distiller(direction = 1, palette ="Greens") +
  theme( plot.background = element_rect(fill ="white", color = "white"), strip.text = element_text(size = 16)) + 
  labs(fill = "Redução relativa", title = "Óbitos em 5 a 11 anos")
