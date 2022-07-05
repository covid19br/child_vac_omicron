library(tidyverse)
library(patchwork)

# Total vacinas
cobertura <- read_csv("DATA/total_vacinas.csv")

# Faixa etária
distribuicao_etaria <-read.csv("DATA/DistrEtaria2020.csv")
colnames(distribuicao_etaria)[1] <- "Idade"
distribuicao_etaria <- distribuicao_etaria[-1,]
rownames(distribuicao_etaria) <- NULL
distribuicao_etaria[,1] <- as.numeric(distribuicao_etaria[,1])

distribuicao_etaria <- rbind(distribuicao_etaria[1,],
                             distribuicao_etaria[2,] + (distribuicao_etaria[3,]*2/5), # 5 - 11
                             distribuicao_etaria[3,]*3/5 + ( distribuicao_etaria[4,]*3/5),  # 12 - 17
                             distribuicao_etaria[4,]*2/5 + colSums(distribuicao_etaria[5:6,]),  # 18 - 29
                             distribuicao_etaria[7:nrow(distribuicao_etaria),])

distribuicao_etaria$Idade <- c(1,2,3,4,rep(5:10,each = 2),10)

eta = distribuicao_etaria %>%
  group_by(Idade) %>% 
  summarise(across(everything(), list(sum))) %>%
  data.frame()

colnames(eta)[-1] <- substr(colnames(eta)[-1],1,2)
eta_long <- eta %>% gather(key = "estado", value = "pop", -Idade)


### 
doses_estados <- data.frame()
files <- list.files("../dados-vacinas/doses_estados/")


for(i in 1:length(files)) {

  da <- read.csv(paste0("../dados-vacinas/doses_estados/",files[i]))

  state = substr(files[i], 17,18)

  da2 <- da %>% filter(data <= "2022-01-17") %>% 
    group_by(vacina, agegroup, doses) %>% 
    summarise(n = sum(ag_child, na.rm = T)) %>%
    mutate(estado = state)

  doses_estados <- bind_rows(da2, doses_estados)
}

reforco <- doses_estados %>% 
  ungroup() %>% 
  filter(doses == "R") %>% 
  group_by(agegroup, estado) %>%
  summarise(m = sum(n, na.rm = T)) %>%
  drop_na(agegroup)

ref = reforco %>% 
  left_join(eta_long, by = c("agegroup" = "Idade","estado")) %>% 
  mutate(p = m/pop)

# Resultados simulações

estados = c( "AC", "AL", "AM", "AP", "BA", "CE", "DF", "ES", "GO", "MA", "MG",
             "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ", "RN", "RO", "RR",
             "RS", "SC", "SE", "SP", "TO")

all_scn_0 <- c()
for(uf in estados) {
  scn0 <- read_csv(paste0("output/cluster/", uf,"_scn0.csv"))
  scn1 <- read_csv(paste0("output/cluster/", uf,"_scn1.csv"))
  scn2 <- read_csv(paste0("output/cluster/", uf,"_scn2.csv"))
  scn0$scenario <- "scn0"
  scn1$scenario <- "scn1"
  scn2$scenario <- "scn2"
  scn0$estado <- uf
  scn1$estado <- uf
  scn2$estado <- uf
  all_scn_0 <- rbind(all_scn_0, scn0, scn1,scn2)
}

all_scn <- all_scn_0 %>% 
  group_by(sample, age_group, compartment, scenario, estado) %>%
  summarize(value = sum(value)) %>%
  pivot_wider(names_from = scenario, values_from = value) %>%
  mutate(diff_scn1 = (scn0-scn1)/scn0, diff_scn2 = (scn0 - scn2)/scn0) %>%
  select(sample, age_group, compartment, estado, diff_scn1, diff_scn2) %>%
  pivot_longer(c(diff_scn1, diff_scn2), names_to = "scenario", values_to = "value")

all2 = all_scn %>% 
  left_join(ref, by = c("age_group" = "agegroup", "estado")) %>%
  mutate(inc = value/pop)

pd <- all2 %>%  
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  filter(age_group > 3) %>%
  ggplot(aes(x = p, y = log(inc), color = factor(age_group))) +
  geom_point() +
  geom_smooth(aes(group = estado), 
              method = "lm", formula = y ~ x,
              show.legend = FALSE) +
  scale_color_viridis_d("Age group", labels = c("18 a 29",
                                                   "30 a 39",
                                                   "40 a 49",
                                                   "50 a 59",
                                                   "60 a 69",
                                                   "70 a 79",
                                                   "80+")) +
  #facet_grid(~compartment+scenario, rows = vars(compartment)) +
  facet_grid(rows = vars(compartment), cols = vars(scenario)) +
    xlab("\nCoverage of population\nwith booster shot by agegroup") +
  ylab("Log outcome averted by vaccination\n") +
  theme_minimal()

pd

ggsave(pd, file = "figures/effect_booster_cov.png", width = 15, height = 12)


pd2 <- all2 %>%  
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  filter(age_group > 3) %>%
  ggplot(aes(x = p, y = inc, color = factor(age_group))) +
  geom_point() +
 # geom_smooth(aes(group = estado), 
 #             method = "lm", formula = y ~ x,
  #            show.legend = FALSE) +
  scale_color_viridis_d("Age group", labels = c("18 a 29",
                                                "30 a 39",
                                                "40 a 49",
                                                "50 a 59",
                                                "60 a 69",
                                                "70 a 79",
                                                "80+")) +
  #facet_grid(~compartment+scenario, rows = vars(compartment)) +
  facet_grid(rows = vars(compartment), cols = vars(scenario)) +
  xlab("\nCoverage of population\nwith booster shot by agegroup") +
  ylab("Log outcome averted by vaccination\n") +
  theme_minimal()

pd2

ggsave(pd2, file = "figures/effect_booster_cov_nolog.png", width = 15, height = 12)


pd3 <- all2 %>%  
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  filter(age_group > 3) %>%
  ggplot(aes(x = p, y = value, color = factor(age_group))) +
  geom_point() +
  geom_smooth(aes(group = age_group),
              method = "lm", formula = y ~ x,
             show.legend = FALSE) +
  scale_color_viridis_d("Age group", labels = c("18 a 29",
                                                "30 a 39",
                                                "40 a 49",
                                                "50 a 59",
                                                "60 a 69",
                                                "70 a 79",
                                                "80+")) +
  #facet_grid(~compartment+scenario, rows = vars(compartment)) +
  facet_grid(rows = vars(compartment), cols = vars(scenario), scales = "free_y") +
  xlab("\nCoverage of population\nwith booster shot by agegroup") +
  ylab("Log outcome averted by vaccination\n") +
  theme_minimal()

pd3

ggsave(pd3, file = "figures/effect_booster_cov_nolog2.png", width = 15, height = 12)


page <- all2 %>%
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  ggplot(aes(x = factor(age_group), y = log(inc))) +
  geom_boxplot(fill = "steelblue") +
  facet_grid(rows = vars(compartment), cols = vars(scenario)) +
  xlab("\nAge group") +
  ylab("Log outcome averted by vaccination\n") +
  theme_minimal()

ggsave(page, file = "figures/effect_agegroup.png", width = 15, height = 12)

#####

head(all2)

diff_child <- all2 %>%
  group_by(sample, compartment, estado, scenario) %>%
  mutate(pop_total = sum(pop, na.rm = T)) %>%
  ungroup() %>%
  #select(-m, pop, p, inc) %>%
  filter(age_group %in% c(2) ) %>% 
  group_by(sample, compartment, estado, scenario) %>%
  summarise(value_child = sum(inc, na.rm = T),
            pop = pop) %>%
  ungroup()


d2 = doses_estados %>% 
  filter(doses == "D2") %>%
  filter(vacina == "AZ") %>%
  group_by(estado) %>%
  summarise(total_vac = sum(n, na.rm = T)) %>%
  ungroup() %>%
  left_join(eta_long %>%
              group_by(estado) %>%
              summarise(pop = sum(pop)) %>%
              ungroup(), by = c("estado")) %>%
  mutate(p = total_vac/pop)

total_ref <- all2 %>%
  select(-m, pop, inc) %>%
  filter(age_group > 3 & agegroup < 8) %>%
  group_by(sample, compartment, estado, scenario)  %>%
  summarise(ref = mean(p, na.rm = T)) %>%
  ungroup() %>%
  left_join(diff_child,
            by = c("sample", "compartment", "estado", "scenario"))

gref <- total_ref %>%
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  ggplot(aes(x = ref, y = log(value_child))) + #, color = estado)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              show.legend = T) +
  ylab("Log da diminuição de eventos em crianças de 5a 11 anos \n") +
  xlab("\nCobertura da população com reforço") +
  facet_wrap(~scenario + compartment, ncol = 2, scale = "free_y")

gref

ggsave(gref, file = "figures/log_children_vs_booster_all_population.png", width = 15, height = 12)
###

diff_child %>% 
  left_join(d2, by = c("estado")) %>%
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  ggplot(aes(x = p, y = log(value_child))) + #, color = estado)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              show.legend = T) +
  ylab("Log da diminuição de eventos em crianças de 5a 11 anos \n") +
  xlab("\nCobertura da população com reforço") +
  facet_wrap(~scenario + compartment, ncol = 2, scale = "free_y")

############
total_ref %>% 
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  group_by(compartment, scenario) %>%
  do(res =  r2b(summary(aov(log(value_child) ~ ref,data = .)))) %>%
  unnest_wider(res)
#####################

x = total_ref %>%
  filter(scenario == "diff_scn1" & compartment == "D") %>%
  glm((value_child) ~ ref, family = "Gamma", data = .) #%>%
       summary()
summary(x)


total_ref %>%
  filter(compartment == "C") %>%
  ggplot(aes(x = estado, y = value_child, fill = scenario)) +
  geom_boxplot() #+
  facet_wrap(~ compartment, ncol = 2, scale = "free")
 
##################################################

pc1d <- all2 %>%
  filter(scenario == "diff_scn1" & compartment == "D") %>%
  filter(age_group > 3) %>%
  mutate(age_group = case_when(age_group == 4 ~ "18 - 29",
                               age_group == 5 ~ "30 - 39",
                               age_group == 6 ~ "40 - 49",
                               age_group == 7 ~ "50 - 59",
                               age_group == 8 ~ "60 - 69",
                               age_group == 9 ~ "70 - 79",
                               age_group == 10 ~ "80+")) %>%
  ggplot(aes(x = p, y = log(inc))) +
  geom_point() +
  geom_smooth(#aes(group = estado), 
              method = "lm", formula = y ~ x,
              show.legend = FALSE) +
  facet_wrap(~age_group, scales = "free_x") +
  xlab("\nBooster shot coverage") +
  ylab("Log death incidence averted by vaccination\n") +
  theme_minimal() +
  labs(title = "Slow pace (current)")

pc2d <- all2 %>%
  filter(scenario == "diff_scn2" & compartment == "D") %>%
  filter(age_group > 3) %>%
  mutate(age_group = case_when(age_group == 4 ~ "18 - 29",
                               age_group == 5 ~ "30 - 39",
                               age_group == 6 ~ "40 - 49",
                               age_group == 7 ~ "50 - 59",
                               age_group == 8 ~ "60 - 69",
                               age_group == 9 ~ "70 - 79",
                               age_group == 10 ~ "80+")) %>%
  ggplot(aes(x = p, y = log(inc))) +
  geom_point() +
  geom_smooth(#aes(group = estado), 
    method = "lm", formula = y ~ x,
    show.legend = FALSE) +
  facet_wrap(~age_group, scales = "free_x") +
  xlab("\nBooster shot coverage") +
  ylab("Log death incidence averted by vaccination\n") +
  theme_minimal() +
  labs(title = "Ideal pace")

pc1c <- all2 %>%
  filter(scenario == "diff_scn1" & compartment == "C") %>%
  filter(age_group > 3) %>%
  mutate(age_group = case_when(age_group == 4 ~ "18 - 29",
                               age_group == 5 ~ "30 - 39",
                               age_group == 6 ~ "40 - 49",
                               age_group == 7 ~ "50 - 59",
                               age_group == 8 ~ "60 - 69",
                               age_group == 9 ~ "70 - 79",
                               age_group == 10 ~ "80+")) %>%
  ggplot(aes(x = p, y = log(inc))) +
  geom_point() +
  geom_smooth(#aes(group = estado), 
    method = "lm", formula = y ~ x,
    show.legend = FALSE) +
  facet_wrap(~age_group, scales = "free_x") +
  xlab("\nBooster shot coverage") +
  ylab("Log hospitalization incidence averted by vaccination\n") +
  theme_minimal() +
  labs(title = "Slow pace (current)")

pc2c <- all2 %>%
  filter(scenario == "diff_scn2" & compartment == "C") %>%
  filter(age_group > 3) %>%
  mutate(age_group = case_when(age_group == 4 ~ "18 - 29",
                               age_group == 5 ~ "30 - 39",
                               age_group == 6 ~ "40 - 49",
                               age_group == 7 ~ "50 - 59",
                               age_group == 8 ~ "60 - 69",
                               age_group == 9 ~ "70 - 79",
                               age_group == 10 ~ "80+")) %>%
  ggplot(aes(x = p, y = inc)) +
  geom_point() +
  geom_smooth(#aes(group = estado), 
    method = "lm", formula = y ~ x,
    show.legend = FALSE) +
  facet_wrap(~age_group, scales = "free_x") +
  xlab("\nBooster shot coverage") +
  ylab("Log hospitalization incidence averted by vaccination\n") +
  theme_minimal() +
  labs(title = "Ideal pace")

pc2c

ggsave(pc1d, file = "figures/booster_slowpace_death.png", width = 12, height = 10)
ggsave(pc2d, file = "figures/booster_idealpace_death.png", width = 12, height = 10)
ggsave(pc1c, file = "figures/booster_slowpace_hosp.png", width = 12, height = 10)
ggsave(pc2c, file = "figures/booster_idealpace_hosp.png", width = 12, height = 10)

####

pcage <- all2 %>%
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  mutate(age_group = case_when(age_group == 4 ~ "18 - 29",
                               age_group == 5 ~ "30 - 39",
                               age_group == 6 ~ "40 - 49",
                               age_group == 7 ~ "50 - 59",
                               age_group == 8 ~ "60 - 69",
                               age_group == 9 ~ "70 - 79",
                               age_group == 10 ~ "80+")) %>%
  ggplot(aes(x = p, y = log(inc), color = age_group)) +
  geom_point() +
  geom_smooth(aes(group = age_group, color = age_group), 
    method = "lm", formula = y ~ x,
    show.legend = FALSE) +
  facet_grid(rows = vars(compartment), cols = vars(scenario)) +
  xlab("\nBooster shot coverage") +
  ylab("Log outcome averted by vaccination\n") +
  theme_minimal() +
  labs(title = "Ideal pace")

pcage
ggsave(pcage, file = "figures/age_effect.png", width = 12, height = 10)
#####
plot(1:99/100, log(1:99/100))

r2 <- function(a) {
 # a = summary(mod)
  x1 = a[[1]][1,2]/sum(a[[1]][,2])
  x2 = a[[1]][2,2]/sum(a[[1]][,2])
  x3 = a[[1]][3,2]/sum(a[[1]][,2])
  return(c(x1,x2,x3))
}

r2b <- function(a) {
  # a = summary(mod)
  x1 = a[[1]][1,2]/sum(a[[1]][,2])
  x2 = a[[1]][2,2]/sum(a[[1]][,2])
  x3 = a[[1]][3,2]/sum(a[[1]][,2])
  x4 = a[[1]][4,2]/sum(a[[1]][,2])
  return(c(x1,x2,x3,x4))
}

x = all2 %>% 
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  filter(age_group > 3) %>%
  group_by(compartment, scenario) %>%
  do(res =  r2(summary(aov(log(inc) ~ factor(age_group) + p,data = .)))) %>%
  unnest_wider(res)
  

y = all2 %>% 
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  filter(age_group > 3) %>%
  group_by(compartment, scenario) %>%
  do(res =  r2(summary(aov(log(inc) ~ p + factor(age_group),data = .)))) %>%
  unnest_wider(res)

colnames(x) <- c("Outcome","Scenario","R2_age_group","R2_booster_coverage","Residual")
colnames(y) <- c("Outcome","Scenario","R2_booster_coverage","R2_age_group","Residual")

z1 = all2 %>% 
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  filter(age_group > 3) %>%
  group_by(compartment, scenario) %>%
  do(res =  r2b(summary(aov(log(inc) ~ estado + p + factor(age_group),data = .)))) %>%
  unnest_wider(res)

z2 = all2 %>% 
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  filter(age_group > 3) %>%
  group_by(compartment, scenario) %>%
  do(res =  r2b(summary(aov(log(inc) ~ estado + factor(age_group) + p,data = .)))) %>%
  unnest_wider(res)
#
z3 = all2 %>% 
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  filter(age_group > 3) %>%
  group_by(compartment, scenario) %>%
  do(res =  r2b(summary(aov(log(inc) ~ p + factor(age_group) + estado ,data = .)))) %>%
  unnest_wider(res)

z4 = all2 %>% 
  mutate(scenario = if_else(scenario == "diff_scn1","Slow pace (current)","Ideal pace")) %>%
  mutate(scenario = factor(scenario, levels = c("Slow pace (current)", "Ideal pace"))) %>%
  mutate(compartment = if_else(compartment == "D","Death incidence","Hospitalization incidence")) %>%
  mutate(compartment = factor(compartment, levels = c("Hospitalization incidence","Death incidence"))) %>%
  filter(age_group > 3) %>%
  group_by(compartment, scenario) %>%
  do(res =  r2b(summary(aov(log(inc) ~ factor(age_group) + p + estado,data = .)))) %>%
  unnest_wider(res)

colnames(z1) <- c("Outcome","Scenario","state","R2_booster_coverage","R2_age_group","Residual")
colnames(z2) <- c("Outcome","Scenario","state","R2_age_group","R2_booster_coverage","Residual")
colnames(z3) <- c("Outcome","Scenario","R2_booster_coverage","R2_age_group","state","Residual")
colnames(z4) <- c("Outcome","Scenario","R2_age_group","R2_booster_coverage","state","Residual")



x = all2 %>% mutate(inc2 = log(inc), p = p*100) %>% 
  group_by(age_group, compartment, scenario, estado) %>%
  summarise(inc = mean(inc), inc2 = mean(inc2), p = mean(p)) %>%
  arrange(age_group, compartment, scenario, inc)
plot(x$inc, x$inc2)

