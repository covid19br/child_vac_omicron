library(readxl)
library(dplyr)
library(forestplot)


VE_table_R <- read_excel("VE_table_R3.xlsx") %>%
  transform(Age_Group = ifelse(is.na(Age_Group),"Assumed",Age_Group)) %>% 
  transform(Age_Group2 = paste("Age:",Age_Group,sep = " ")) %>%
  transform(Dose_brand = paste(Brand,Dose,sep = ":")) %>%
  transform(VE_CI = paste(VE,"(",LCI,"-",UCI,")")) %>%
  transform(Name_Dose_age = paste(Name,Age_Group2,sep = " ,")) 


infection_table <- VE_table_R %>% filter(Outcome=="Infection") %>%
  arrange(Dose,Brand)

symp_inf_table <- VE_table_R %>% filter(Outcome=="Symptomatic Disease") %>%  
  arrange(Dose,Brand)

Hospital_table <- VE_table_R %>% filter(Outcome=="Hospital Admission")  %>%  
  arrange(Dose,Brand)
several_table <- VE_table_R %>% filter(Outcome=="MISC"
                                       |Outcome=="Several/Critical/Fatal Disease"|
                                         Outcome=="Hospital Admission")  %>%  
  arrange(Dose,Brand)
death_table <- VE_table_R %>% filter(Outcome=="Death")  %>%  
  arrange(Dose,Brand)

# gerar especificas para plot
infection_table2 <-  infection_table %>% dplyr::select(Name,Dose_brand,VE_CI,Age_Group,Use,VE,LCI,UCI)
symp_table2 <-  symp_inf_table %>% dplyr::select(Name,Dose_brand,VE_CI,Age_Group,Use,VE,LCI,UCI)
several_table2 <-  several_table %>% dplyr::select(Name,Dose_brand,VE_CI,Age_Group,Use,VE,LCI,UCI)
hosp_table2 <-  Hospital_table %>% dplyr::select(Name,Dose_brand,VE_CI,Age_Group,Use,VE,LCI,UCI)
death_table2 <-  death_table %>% dplyr::select(Name,Dose_brand,VE_CI,Age_Group,Use,VE,LCI,UCI)

# Colocar nomes
row_names <- cbind(
  c("Study", infection_table2$Name),
  c("VE (95% CI)", infection_table2$VE_CI),
  c("Dose/Brand", infection_table2$Dose_brand),
  c("Age", infection_table2$Age_Group),
  c("Used?", infection_table2$Use)
)

row_names2 <- cbind(
  c("Study", symp_table2$Name),
  c("VE (95% CI)", symp_table2$VE_CI),
  c("Dose/Brand", symp_table2$Dose_brand),
  c("Age", symp_table2$Age_Group),
  c("Used?", symp_table2$Use)
)

row_names3 <- cbind(
  c("Study", several_table2$Name),
  c("VE (95% CI)", several_table2$VE_CI),
  c("Dose/Brand", several_table2$Dose_brand),
  c("Age", several_table2$Age_Group),
  c("Used?", several_table2$Use)
)


row_names4 <- cbind(
  c("Study", hosp_table2$Name),
  c("VE (95% CI)", hosp_table2$VE_CI),
  c("Dose/Brand", hosp_table2$Dose_brand),
  c("Age", hosp_table2$Age_Group),
  c("Used?", hosp_table2$Use)
)


row_names5 <- cbind(
  c("Study", death_table2$Name),
  c("VE (95% CI)", death_table2$VE_CI),
  c("Dose/Brand", death_table2$Dose_brand),
  c("Age", death_table2$Age_Group),
  c("Used?", death_table2$Use)
)

infection_table2 <- rbind(rep(NA, 3), infection_table2)
symp_table2 <- rbind(rep(NA, 3), symp_table2)
several_table2 <- rbind(rep(NA, 3), several_table2)
hosp_table2 <- rbind(rep(NA, 3), hosp_table2)
death_table2 <- rbind(rep(NA, 3), death_table2)


pdf(width = 10,height = 7, pointsize = 10,
    "VE_parameters.pdf")
infection_table2 %>%
  forestplot(row_names, 
             mean = VE, lower = LCI,upper= UCI,
             title="",
             zero= 0,lwd.zero = gpar(lwd=1),
             boxsize=0.25,
             xticks=c(-50,seq(-50,100,10)),
             line.margin = .2,
             lty.ci = 2,
             clip=c(-50, 100),
             plotwidth=unit(15, "cm"),
             colgap=unit(2, "mm"),
             lineheight = "auto",
             xlab="VE for Infection (95% CI)")

symp_table2 %>%
  forestplot(row_names2, 
             mean = VE, lower = LCI,upper= UCI,
             title="",
             zero= 0,lwd.zero = gpar(lwd=1),
             boxsize=0.25,
             xticks=c(-50,seq(-50,100,10)),
             line.margin = .2,
             lty.ci = 2,
             clip=c(-50, 100),
             plotwidth=unit(15, "cm"),
             colgap=unit(2, "mm"),
             lineheight = "auto",
             xlab="VE for Symptomatic Disease (95% CI)")


several_table2 %>%
  forestplot(row_names3, 
             mean = VE, lower = LCI,upper= UCI,
             title="",
             zero= 0,lwd.zero = gpar(lwd=1),
             boxsize=0.25,
             xticks=c(-50,seq(-50,100,10)),
             line.margin = .2,
             lty.ci = 2,
             clip=c(-50, 100),
             plotwidth=unit(15, "cm"),
             colgap=unit(2, "mm"),
             lineheight = "auto",
             xlab="VE for Several outcomes (95% CI)")


hosp_table2 %>%
  forestplot(row_names4, 
             mean = VE, lower = LCI,upper= UCI,
             title="",
             zero= 0,lwd.zero = gpar(lwd=1),
             boxsize=0.25,
             xticks=c(-50,seq(-50,100,10)),
             line.margin = .2,
             lty.ci = 2,
             clip=c(-50, 100),
             plotwidth=unit(15, "cm"),
             colgap=unit(2, "mm"),
             lineheight = "auto",
             xlab="VE for Hospital Admission (95% CI)")

death_table2 %>%
  forestplot(row_names5, 
             mean = VE, lower = LCI,upper= UCI,
             title="",
             zero= 0,lwd.zero = gpar(lwd=1),
             boxsize=0.25,
             xticks=c(-50,seq(-50,100,10)),
             line.margin = .2,
             lty.ci = 2,
             clip=c(-50, 100),
             plotwidth=unit(15, "cm"),
             colgap=unit(2, "mm"),
             lineheight = "auto",
             xlab="VE for Death (95% CI)")
dev.off()


