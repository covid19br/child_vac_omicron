if(!require(gridExtra)){install.packages("gridExtra"); library(gridExtra)}
if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(stringi)){install.packages("stringi"); library(stringi)}

# Original IHR

ihr <- read.csv("DATA/14. ihr.csv")
colnames(ihr) <- c("agegroup","ihr_original")
ihr <- ihr[seq(1,19,2),]
ihr$agegroup <- 1:10


# Vaccination efficacy
vea = read.csv("DATA/vea.csv")
vea$eff <- vea$eff/100
vea$upper <- vea$upper/100
vea$lower <- vea$lower/100

# Distribuição etária

distribuicao_etaria <- read.csv("DATA/DistrEtaria2020.csv")
colnames(distribuicao_etaria)[1] <- "Idade"
distribuicao_etaria <- distribuicao_etaria[-1,]
rownames(distribuicao_etaria) <- NULL
distribuicao_etaria$agegroup <- c(rep(1:9, each = 2), 10)

eta = distribuicao_etaria %>% select(-Idade) %>% 
  group_by(agegroup) %>% 
  summarise(across(everything(), list(sum))) %>%
  data.frame()

colnames(eta)[-1] <- substr(colnames(eta)[-1],1,2)
eta_long <- eta %>% gather(key = "sg_uf", value = "pop", -agegroup)

eta_total <- eta_long %>% group_by(sg_uf) %>% summarise(t = sum(pop))

# Vaccination

vac_files <- list.files("DATA/doses_aplicadas_10")

vac_all <- data.frame()

for(j in vac_files) {
  df <- read.csv(paste0("DATA/doses_aplicadas_10/", j))
  state = substr(j, unlist(gregexpr("\\.csv", j))-2, nchar(j)-4)
  print(state)
  
  df2 = df %>% filter(doses == "D2") %>%
  mutate(data = as.Date(data)) %>%
  mutate(data_vac = case_when(vacina == "Pfizer" ~ data + 7, 
                              vacina != "Pfizer" ~ data + 14)) %>%
  select(-data, -doses) %>%
  mutate(sg_uf = state)
  
  vac_all <- rbind(vac_all, df2)
}

vac2 = vac_all %>% 
  filter(vacina != "Janssen") %>%
  drop_na(vacina, agegroup, sg_uf) %>%
  group_by(vacina, agegroup, sg_uf) %>% 
  mutate(total = cumsum(n)) %>%
  arrange(sg_uf, vacina, agegroup, data_vac) %>%
  left_join(eta_long, by = c("sg_uf", "agegroup"),
            na_matches = "never") %>%
  mutate(prop = total/pop) %>%
  left_join(ihr, by = c("agegroup")) %>%
  left_join(vea %>%
              select(-upper, - lower), by = c("vacina","agegroup")) %>%
  mutate(prop = case_when(prop > 1 ~ 1.0,
                          prop <= 1 ~ prop)) %>%
  #mutate(w = ((prop*(1-eff))+(1-prop))) %>%
  mutate(w = (prop*(1-eff))-prop) %>%
  group_by(agegroup, data_vac, sg_uf) %>%
  mutate(w2 = sum(w)) %>%
  mutate(w3 = 1 - w2) %>%
  mutate(ihr = ihr_original * w3)
  
  vac3 <- vac2 %>%
    select(agegroup, data_vac, sg_uf, ihr, ihr_original) %>%
    distinct()

vac3 %>% ggplot(aes(x = ihr_original, y = ihr)) +
  geom_point()


# Hospitalizations

sivep <- read_csv("DATA/sum_dados_2022_01_31.csv.xz")
sivep$agegroup <- as.numeric(substr(sivep$age_class,5,nchar(sivep$age_class)))

sivep_long <- sivep %>%
  drop_na(dt_evento, agegroup, sg_uf) %>%
  group_by(dt_evento, agegroup, sg_uf) %>%
  summarise(n = sum(srag, na.rm = T))#%>% 
 # spread(key = "sg_uf", value = "n", fill = 0)


##### Estimativa do número de casos

cases <- sivep_long %>% 
  left_join(vac3 %>%
              select(data_vac, agegroup, sg_uf, ihr), 
            by = c("dt_evento" = "data_vac", "agegroup", "sg_uf"),
            fill = NA) %>%
  left_join(ihr, by = "agegroup") %>%
  mutate(ihr_final = case_when(!is.na(ihr) ~ ihr,
                               is.na(ihr)  ~ ihr_original)) %>%
  mutate(ihr_final = ihr_final/100) %>%
  mutate(cases = n/ihr_final)

### Plots
cases %>%
  ggplot(aes(x = dt_evento, y = cases, color = factor(agegroup))) +
  geom_point() +
  facet_wrap(~sg_uf)

x = cases %>%
  group_by(dt_evento, sg_uf) %>%
  summarise(m = sum(cases, na.rm = T)) %>%
  group_by(sg_uf) %>%
  mutate(total = cumsum(m)) %>%
  left_join(eta_total, by = "sg_uf") %>%
  mutate(prev = total/t*100)

g1 <- x %>%
  ggplot(aes(x = dt_evento, y = prev)) +
  geom_line() +
  geom_hline(yintercept = 100, color = "red", linetype = "dashed") +
  facet_wrap(~sg_uf) +
  xlab("") + ylab("Prevalência (%)\n") +
  scale_x_date(date_breaks = "6 months", date_labels = "%m-%y")


ggsave(g1, file = "figures/prevalence_estimates.png", width = 15, height = 10, dpi = 200)
