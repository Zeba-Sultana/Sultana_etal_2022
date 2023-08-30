#!/usr/bin/env Rscript

library(readxl) 
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr) # needed for compare_means
library(egg) #needed for set_panel_size
#library(knitr) # for making HTML tables using the function kable()

source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

Mekinh_PD_Doses =c("0","0.4","1.3","4","12","37","111","333","1000")

quant_folder= "../../RAW_DATA/ValidationExperiments/MEKi_TC_DR/Analyte_Quant/"
empiria_folder = "../../RAW_DATA/ValidationExperiments/MEKi_TC_DR/TPS_Quant/"

dir.create("./OUTPUT_MekiDR") # Create folder to save all output data

# Create folder OUTPUT_PAPER to save figs used in paper
if(!(file.exists("./OUTPUT_PAPER"))){ 
  dir.create("./OUTPUT_PAPER") 
}

##################################################################
####### pMek Normalized over total protein stain   ###############
##################################################################

################ *Replicate 1 - XX : Gel24* ################
Gel24_samples = c(
  "Common_nn_nn_nn",
  "XO_R1_cntrl_PDDR",
  "XO_R1_0.4_PDDR",
  "XO_R1_1.3_PDDR",
  "XO_R1_4_PDDR",
  "XO_R1_12_PDDR",
  "XO_R1_37_PDDR",
  "XO_R1_111_PDDR",
  "XO_R1_333_PDDR",
  "XO_R1_1000_PDDR",
  "XX_R1_cntrl_PDDR",
  "XX_R1_0.4_PDDR",
  "XX_R1_1.3_PDDR",
  "XX_R1_4_PDDR",
  "XX_R1_12_PDDR",
  "XX_R1_37_PDDR",
  "XX_R1_111_PDDR",
  "XX_R1_333_PDDR",
  "XX_R1_1000_PDDR")
Gel24_quant_file = c("0000126_02_JuneGel24_PDDR_Mek.xls")
Gel24_Mek <-  read_excel(file.path(quant_folder,Gel24_quant_file))

Gel24_TPS_file = "0000141_02_June2020_MekiDR_Gel24.xlsx"
Gel24_TPS <- read_excel(file.path(empiria_folder,Gel24_TPS_file))

Gel24_MekTPS_labeled <- ReadInWBData_TPS(Gel24_samples,Gel24_Mek, Gel24_TPS, Mekinh_PD_Doses, "Gel24", "Mek")
Gel24_MekTPS_labeled <- Gel24_MekTPS_labeled %>% 
  filter(Cell_line != "Common")

Gel24_MekTPS_NormOMeanRep <- Norm_o_MeanRep(Gel24_MekTPS_labeled)

################ *Replicate 2 - XX&XO : Gel25*  ################

Gel25_samples = c(
  "Common_nn_nn_nn",
  "XX_R2_cntrl_PDDR",
  "XX_R2_0.4_PDDR",
  "XX_R2_1.3_PDDR",
  "XX_R2_4_PDDR",
  "XX_R2_12_PDDR",
  "XX_R2_37_PDDR",
  "XX_R2_111_PDDR",
  "XX_R2_333_PDDR",
  "XX_R2_1000_PDDR",
  "XO_R2_cntrl_PDDR",
  "XO_R2_0.4_PDDR",
  "XO_R2_1.3_PDDR",
  "XO_R2_4_PDDR",
  "XO_R2_12_PDDR",
  "XO_R2_37_PDDR",
  "XO_R2_111_PDDR",
  "XO_R2_333_PDDR",
  "XO_R2_1000_PDDR")
Gel25_quant_file = c("0000125_02_JuneGel25_PDDR_Mek.xls")
Gel25_Mek <-  read_excel(file.path(quant_folder,Gel25_quant_file))

Gel25_TPS_file = "0000140_02_June2020_MekiDR_Gel25_straight.xlsx"
Gel25_TPS <- read_excel(file.path(empiria_folder,Gel25_TPS_file))

Gel25_MekTPS_labeled <- ReadInWBData_TPS(Gel25_samples,Gel25_Mek, Gel25_TPS, Mekinh_PD_Doses, "Gel25", "Mek")
Gel25_MekTPS_labeled <- Gel25_MekTPS_labeled %>% 
  filter(Cell_line != "Common")

Gel25_MekTPS_NormOMeanRep <- Norm_o_MeanRep(Gel25_MekTPS_labeled)

################ *Replicate 3 - XX&XO : Gel16*  ################
Gel16_samples = c(
  "Common_nn_nn_nn",
  "XO_R3_cntrl_PDDR",
  "XO_R3_0.4_PDDR",
  "XO_R3_1.3_PDDR",
  "XO_R3_4_PDDR",
  "XO_R3_12_PDDR",
  "XO_R3_37_PDDR",
  "XO_R3_111_PDDR",
  "XO_R3_333_PDDR",
  "XO_R3_1000_PDDR",
  "XX_R3_cntrl_PDDR",
  "XX_R3_0.4_PDDR",
  "XX_R3_1.3_PDDR",
  "XX_R3_4_PDDR",
  "XX_R3_12_PDDR",
  "XX_R3_37_PDDR",
  "XX_R3_111_PDDR",
  "XX_R3_333_PDDR",
  "XX_R3_1000_PDDR")
Gel16_quant_file = c("0000039_01_JuneGel16_PDDR_Mek.xls")
Gel16_Mek <-  read_excel(file.path(quant_folder,Gel16_quant_file))

Gel16_TPS_file = "0000055_02_PDDR_R3_JuneGel16.xlsx"
Gel16_TPS <- read_excel(file.path(empiria_folder,Gel16_TPS_file))

Gel16_MekTPS_labeled <- ReadInWBData_TPS(Gel16_samples,Gel16_Mek, Gel16_TPS, Mekinh_PD_Doses, "Gel16", "Mek")
Gel16_MekTPS_labeled <- Gel16_MekTPS_labeled %>% 
  filter(Cell_line != "Common")

Gel16_MekTPS_NormOMeanRep <- Norm_o_MeanRep(Gel16_MekTPS_labeled)

################ COMBINING ALL 3 REPLICATES #####################

PDDR_MekTPS <- dplyr::bind_rows(Gel24_MekTPS_NormOMeanRep, Gel25_MekTPS_NormOMeanRep, Gel16_MekTPS_NormOMeanRep)

PDDR_MekTPS_FC <- FC_Calculation_updated(PDDR_MekTPS)

PDDR_MekTPS_plot <- PDDR_MekTPS_FC %>% 
  select(-c(Exp,Gel )) 


##################################################################
####### pcRaf Normalized over total protein stain   ###############
##################################################################

################ *Replicate 1 - XX : Gel24* ################

Gel24_quant_file = c("0000139_02_JuneGel24_PDDR_pcRaf.xls")
Gel24_cRaf <-  read_excel(file.path(quant_folder,Gel24_quant_file))

Gel24_TPS_file = "0000141_02_June2020_MekiDR_Gel24.xlsx"
Gel24_TPS <- read_excel(file.path(empiria_folder,Gel24_TPS_file))

Gel24_cRafTPS_labeled <- ReadInWBData_TPS(Gel24_samples,Gel24_cRaf, Gel24_TPS, Mekinh_PD_Doses, "Gel24", "cRaf")
Gel24_cRafTPS_labeled <- Gel24_cRafTPS_labeled %>% 
  filter(Cell_line != "Common")

Gel24_cRafTPS_NormOMeanRep <- Norm_o_MeanRep(Gel24_cRafTPS_labeled)

################ *Replicate 2 - XX&XO : Gel25*  ################
Gel25_quant_file = c("0000138_02_JuneGel25_PDDR_pcRaf.xls")
Gel25_cRaf <-  read_excel(file.path(quant_folder,Gel25_quant_file))

Gel25_TPS_file = "0000140_02_June2020_MekiDR_Gel25_straight.xlsx"
Gel25_TPS <- read_excel(file.path(empiria_folder,Gel25_TPS_file))

Gel25_cRafTPS_labeled <- ReadInWBData_TPS(Gel25_samples,Gel25_cRaf, Gel25_TPS, Mekinh_PD_Doses, "Gel25", "cRaf")
Gel25_cRafTPS_labeled <- Gel25_cRafTPS_labeled %>% 
  filter(Cell_line != "Common")

Gel25_cRafTPS_NormOMeanRep <- Norm_o_MeanRep(Gel25_cRafTPS_labeled)

################ *Replicate 3 - XX&XO : Gel16*  ################
Gel16_quant_file = c("0000041_02_JuneGel16_PDDR_pcRaf.xls")
Gel16_cRaf <-  read_excel(file.path(quant_folder,Gel16_quant_file))

Gel16_TPS_file = "0000055_02_PDDR_R3_JuneGel16.xlsx"
Gel16_TPS <- read_excel(file.path(empiria_folder,Gel16_TPS_file))

Gel16_cRafTPS_labeled <- ReadInWBData_TPS(Gel16_samples,Gel16_cRaf, Gel16_TPS, Mekinh_PD_Doses, "Gel16", "cRaf")
Gel16_cRafTPS_labeled <- Gel16_cRafTPS_labeled %>% 
  filter(Cell_line != "Common")

Gel16_cRafTPS_NormOMeanRep <- Norm_o_MeanRep(Gel16_cRafTPS_labeled)

################ COMBINING ALL 3 REPLICATES #####################

PDDR_cRafTPS <- dplyr::bind_rows(Gel24_cRafTPS_NormOMeanRep, Gel25_cRafTPS_NormOMeanRep, Gel16_cRafTPS_NormOMeanRep)

PDDR_cRafTPS_FC <- FC_Calculation_updated(PDDR_cRafTPS)

PDDR_cRaf_plot <- PDDR_cRafTPS_FC %>% 
  select(-c(Exp,Gel )) 

############# FIGURES : FC over XX ###################

##################################################################
####### PLOTS : Using pMek Over Total Protein Stain   ############
##################################################################

PDDR_All_MekTPS <- dplyr::full_join(PDDR_MekTPS_plot, PDDR_cRaf_plot, by=c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr"))

colnames(PDDR_All_MekTPS) <- gsub("\\.x","_Mek",colnames(PDDR_All_MekTPS))
colnames(PDDR_All_MekTPS) <- gsub("\\.y","_cRaf",colnames(PDDR_All_MekTPS))

PDDR_All_MekTPS <- PDDR_All_MekTPS %>%
  select(-c(grep("Analyte", colnames(PDDR_All_MekTPS))))

PDDR_All_MekTPSplot_FCoXX <- PDDR_All_MekTPS %>% 
  select(c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr",grep("FCoXX_M2",colnames(PDDR_All_MekTPS)))) %>% 
  gather(key = "Analyte", value = "Signal", -c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr"))

PDDR_All_MekTPSplot_FCoXX$Analyte = factor(PDDR_All_MekTPSplot_FCoXX$Analyte, levels =c("FCoXX_M2_Mek", "FCoXX_M2_cRaf")) # This is to get the correct sequence in the Facet_wrap
Analyte_labels <- c("FCoXX_M2_Mek" = "pMEK", "FCoXX_M2_cRaf" = "pRAF1")

##############Save this file for further Manipulation ############
write.csv(PDDR_All_MekTPSplot_FCoXX, file = "./OUTPUT_MekiDR/PDDR_All_MekTPSplot_FCoXX.csv")


##### These calculation are done to find the max y-values for the 2 panels. Inspite of using "free_y", I wanted to add some ##### 
#### space on top of the plot to accomodate the significance stars ######## 
#### The max of y-values for each of the panels is calculated and plotted on top of the plot by specifiying them to ##### 
#### the corresponding analyte panel ##### 

pMek_max <- PDDR_All_MekTPSplot_FCoXX %>%
  filter(grepl("Mek", Analyte))
pMek_max <-max(pMek_max$Signal) + 2 # the number is based on how much space is needed

pcRaf_max <- PDDR_All_MekTPSplot_FCoXX %>%
  filter(grepl("Raf", Analyte))
pcRaf_max <-max(pcRaf_max$Signal) + 0.5 # the number is based on how much space is needed

dummy_Mek <- data.frame(log2Treatment = 0, Signal = pMek_max,
                        Analyte = "FCoXX_M2_Mek", Cell_line = "Common", stringsAsFactors=FALSE)
dummy_Raf <- data.frame(log2Treatment = 0, Signal = pcRaf_max,
                        Analyte = "FCoXX_M2_cRaf", Cell_line = "Common", stringsAsFactors=FALSE) 
dummy_data <- rbind(dummy_Mek,dummy_Raf)
dummy_data$Analyte = factor(dummy_data$Analyte, levels =c("FCoXX_M2_Mek", "FCoXX_M2_cRaf")) # This is to get the correct sequence in the Facet_wrap

#################### PLOT ###########

g <- Plot_TwoPanel_ValidationPlot_updated(PDDR_All_MekTPSplot_FCoXX,"log2Treatment","Signal", Analyte_labels,"free_y")+
  geom_blank(data=dummy_data) + 
  labs(x = "\n MEKi(nM)+1 [log2]", 
       y = " Rel. phosp. (norm.)  \n",
       color = "Cell line" )

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("PDDR_pMekTPS_pcRaf_FCoXX_MEANline.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_MekiDR")
ggsave("Fig7I_PDDR_pMekTPS_pcRaf_FCoXX_MEANline.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 



############# FIGURES : FC over respective Cellline control ###################

##################################################################
####### PLOTS :  pMek Over Total Protein Stain   ############
##################################################################

PDDR_All_MekTPS <- dplyr::full_join(PDDR_MekTPS_plot, PDDR_cRaf_plot, by=c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr"))

colnames(PDDR_All_MekTPS) <- gsub("\\.x","_Mek",colnames(PDDR_All_MekTPS))
colnames(PDDR_All_MekTPS) <- gsub("\\.y","_cRaf",colnames(PDDR_All_MekTPS))

PDDR_All_MekTPS <- PDDR_All_MekTPS %>%
  select(-c(grep("Analyte", colnames(PDDR_All_MekTPS))))

PDDR_All_MekTPSplot_FCoCntrl <- PDDR_All_MekTPS %>% 
  select(c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr",grep("FCoCntrl_M2",colnames(PDDR_All_MekTPS)))) %>% 
  gather(key = "Analyte", value = "Signal", -c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr"))

PDDR_All_MekTPSplot_FCoCntrl$Analyte = factor(PDDR_All_MekTPSplot_FCoCntrl$Analyte, levels =c("FCoCntrl_M2_Mek", "FCoCntrl_M2_cRaf")) # This is to get the correct sequence in the Facet_wrap
Analyte_labels <- c("FCoCntrl_M2_Mek" = "pMEK", "FCoCntrl_M2_cRaf" = "pRAF1")

##############Save this file for further Manipulation ############
write.csv(PDDR_All_MekTPSplot_FCoCntrl, file = "./OUTPUT_MekiDR/PDDR_All_MekTPSplot_FCoCntrl.csv")


# These calculation are done to find the max y-values for the 2 panels. Inspite of using "free_y", I wanted to add some ##### 
# space on top of the plot to accommodate the significance stars ######## 
# The max of y-values for each of the panels is calculated and plotted on top of the plot by specifiying them to ##### 
# the corresponding analyte panel ##### 

pMek_max <- PDDR_All_MekTPSplot_FCoCntrl %>%
  filter(grepl("Mek", Analyte))
pMek_max <-max(pMek_max$Signal) + 5 # the number is based on how much space i need

pcRaf_max <- PDDR_All_MekTPSplot_FCoCntrl %>%
  filter(grepl("Raf", Analyte))
pcRaf_max <-max(pcRaf_max$Signal) + 0.15 # the number is based on how much space i need

dummy_Mek <- data.frame(log2Treatment = 0, Signal = pMek_max,
                        Analyte = "FCoCntrl_M2_Mek", Cell_line = "Common", stringsAsFactors=FALSE)
dummy_Raf <- data.frame(log2Treatment = 0, Signal = pcRaf_max,
                        Analyte = "FCoCntrl_M2_cRaf", Cell_line = "Common", stringsAsFactors=FALSE) 
dummy_data <- rbind(dummy_Mek,dummy_Raf)
dummy_data$Analyte = factor(dummy_data$Analyte, levels =c("FCoCntrl_M2_Mek", "FCoCntrl_M2_cRaf")) # This is to get the correct sequence in the Facet_wrap

#################### PLOT ###########

g <- Plot_TwoPanel_ValidationPlot_updated(PDDR_All_MekTPSplot_FCoCntrl,"log2Treatment","Signal", Analyte_labels,"free_y")+
  geom_blank(data=dummy_data) + 
  labs(x = "\n MEKi(nM)+1 [log2]", #\u03bc is the unicode charachter fro greek mu
       y = " Phosph. rel. to \n untreated \n",
       color = "Cell line" )

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("PDDR_pMekTPS_pcRaf_FCoCntrl_MEANline.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_MekiDR") 

########### T Tests for difference between XX and XO in case of FCoCntrl ########

Meki_dose = c(0,0.4,1.3,4,12,37,111,333,1000)

PDDR_XX_XO_Mek_pvals_FCoCntrl = list()
for(i in 1:length(Meki_dose)) {
  
  
  PDDR_XX_Mek_values <- PDDR_All_MekTPSplot_FCoCntrl %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "FCoCntrl_M2_Mek") %>% 
    filter(Treatment == Meki_dose[i])
  
  PDDR_XO_Mek_values <- PDDR_All_MekTPSplot_FCoCntrl %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "FCoCntrl_M2_Mek") %>% 
    filter(Treatment == Meki_dose[i])
  
  PDDR_XX_XO_Mek_pvals_FCoCntrl$Treatment_value[i] = Meki_dose[i]
  PDDR_XX_XO_Mek_pvals_FCoCntrl$XX_mean[i] = t.test(PDDR_XX_Mek_values$Signal,PDDR_XO_Mek_values$Signal, paired = FALSE)$estimate[1]
  PDDR_XX_XO_Mek_pvals_FCoCntrl$XO_mean[i] = t.test(PDDR_XX_Mek_values$Signal,PDDR_XO_Mek_values$Signal, paired = FALSE)$estimate[2]
  PDDR_XX_XO_Mek_pvals_FCoCntrl$p_value[i] = t.test(PDDR_XX_Mek_values$Signal,PDDR_XO_Mek_values$Signal, paired = FALSE)$p.value
  
}

PDDR_Mek_FCoCntrl_Ttest <- data.frame((sapply(PDDR_XX_XO_Mek_pvals_FCoCntrl,c)))
PDDR_Mek_FCoCntrl_Ttest <- PDDR_Mek_FCoCntrl_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))
WriteXLS::WriteXLS(PDDR_Mek_FCoCntrl_Ttest, "./OUTPUT_MekiDR/PDDR_Mek_FCoCntrl_Ttest.xls")


PDDR_XX_XO_cRaf_pvals_FCoCntrl = list()
for(i in 1:length(Meki_dose) ) {
  

  PDDR_XX_cRaf_values <- PDDR_All_MekTPSplot_FCoCntrl %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "FCoCntrl_M2_cRaf") %>% 
    filter(Treatment == Meki_dose[i])
  
  PDDR_XO_cRaf_values <- PDDR_All_MekTPSplot_FCoCntrl %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "FCoCntrl_M2_cRaf") %>% 
    filter(Treatment == Meki_dose[i])
  
  PDDR_XX_XO_cRaf_pvals_FCoCntrl$Treatment_value[i] = Meki_dose[i]
  PDDR_XX_XO_cRaf_pvals_FCoCntrl$XX_mean[i] = t.test(PDDR_XX_cRaf_values$Signal,PDDR_XO_cRaf_values$Signal, paired = FALSE)$estimate[1]
  PDDR_XX_XO_cRaf_pvals_FCoCntrl$XO_mean[i] = t.test(PDDR_XX_cRaf_values$Signal,PDDR_XO_cRaf_values$Signal, paired = FALSE)$estimate[2]
  PDDR_XX_XO_cRaf_pvals_FCoCntrl$p_value[i] = t.test(PDDR_XX_cRaf_values$Signal,PDDR_XO_cRaf_values$Signal, paired = FALSE)$p.value
  
}

PDDR_cRaf_FCoCntrl_Ttest <- data.frame((sapply(PDDR_XX_XO_cRaf_pvals_FCoCntrl,c)))
PDDR_cRaf_FCoCntrl_Ttest <- PDDR_cRaf_FCoCntrl_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))
WriteXLS::WriteXLS(PDDR_cRaf_FCoCntrl_Ttest, "./OUTPUT_MekiDR/PDDR_cRaf_FCoCntrl_Ttest.xls")


########### T Tests for difference between XX and XO in case of FCoXX ########

Meki_dose = c(0,0.4,1.3,4,12,37,111,333,1000)

PDDR_XX_XO_Mek_pvals_FCoXX = list()
for(i in 1:length(Meki_dose) ) {
  
  PDDR_XX_Mek_values <- PDDR_All_MekTPSplot_FCoXX %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "FCoXX_M2_Mek") %>% 
    filter(Treatment == Meki_dose[i])
  
  PDDR_XO_Mek_values <- PDDR_All_MekTPSplot_FCoXX %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "FCoXX_M2_Mek") %>% 
    filter(Treatment == Meki_dose[i])
  
  PDDR_XX_XO_Mek_pvals_FCoXX$Treatment_value[i] = Meki_dose[i]
  PDDR_XX_XO_Mek_pvals_FCoXX$XX_mean[i] = t.test(PDDR_XX_Mek_values$Signal,PDDR_XO_Mek_values$Signal, paired = FALSE)$estimate[1]
  PDDR_XX_XO_Mek_pvals_FCoXX$XO_mean[i] = t.test(PDDR_XX_Mek_values$Signal,PDDR_XO_Mek_values$Signal, paired = FALSE)$estimate[2]
  PDDR_XX_XO_Mek_pvals_FCoXX$p_value[i] = t.test(PDDR_XX_Mek_values$Signal,PDDR_XO_Mek_values$Signal, paired = FALSE)$p.value
  
}

PDDR_Mek_FCoXX_Ttest <- data.frame((sapply(PDDR_XX_XO_Mek_pvals_FCoXX,c)))


PDDR_XX_XO_cRaf_pvals_FCoXX = list()
for(i in 1:length(Meki_dose) ) {
  
  PDDR_XX_cRaf_values <- PDDR_All_MekTPSplot_FCoXX %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "FCoXX_M2_cRaf") %>% 
    filter(Treatment == Meki_dose[i])
  
  PDDR_XO_cRaf_values <- PDDR_All_MekTPSplot_FCoXX %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "FCoXX_M2_cRaf") %>% 
    filter(Treatment == Meki_dose[i])
  
  PDDR_XX_XO_cRaf_pvals_FCoXX$Treatment_value[i] = Meki_dose[i]
  PDDR_XX_XO_cRaf_pvals_FCoXX$XX_mean[i] = t.test(PDDR_XX_cRaf_values$Signal,PDDR_XO_cRaf_values$Signal, paired = FALSE)$estimate[1]
  PDDR_XX_XO_cRaf_pvals_FCoXX$XO_mean[i] = t.test(PDDR_XX_cRaf_values$Signal,PDDR_XO_cRaf_values$Signal, paired = FALSE)$estimate[2]
  PDDR_XX_XO_cRaf_pvals_FCoXX$p_value[i] = t.test(PDDR_XX_cRaf_values$Signal,PDDR_XO_cRaf_values$Signal, paired = FALSE)$p.value
  
}

PDDR_cRaf_FCoXX_Ttest <- data.frame((sapply(PDDR_XX_XO_cRaf_pvals_FCoXX,c)))


########## Saving the data for fig 7 F #########
PDDR_Mek_FCoXX = PDDR_All_MekTPSplot_FCoXX %>% 
  filter(grepl("Mek", Analyte))

PDDR_Raf_FCoXX = PDDR_All_MekTPSplot_FCoXX %>% 
  filter(grepl("Raf", Analyte))


############ pMek values #######
PDDR_Mek_FCoXX_Ttest <- PDDR_Mek_FCoXX_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  mutate(Treatment=as.numeric(Treatment_value)) %>% 
  select(-c(Treatment_value))

PDDR_Mek_ReplicateValues  = PDDR_Mek_FCoXX %>% #colnames()
  select(Cell_line,Replicate,Treatment,log2Treatment,Analyte,Signal) %>% 
  pivot_wider(names_from = Replicate,values_from=Signal) %>% 
  rowwise() %>% 
  mutate(Mean=mean(c(R1,R2,R3), na.rm = TRUE))

Fig7F_MekiDR_Mek_data = left_join(PDDR_Mek_ReplicateValues,PDDR_Mek_FCoXX_Ttest, by="Treatment") 
Fig7F_MekiDR_Mek_data = Fig7F_MekiDR_Mek_data %>% 
  mutate(XX_mean = ifelse(Cell_line == "XO", NA, XX_mean)) %>% 
  mutate(XO_mean = ifelse(Cell_line == "XO", NA, XO_mean)) %>% 
  mutate(p_value = ifelse(Cell_line == "XO", NA, p_value)) %>% 
  mutate(sig = ifelse(Cell_line == "XO", NA, sig))
colnames(Fig7F_MekiDR_Mek_data) = gsub("p_value","p_value(XXvsXO)",colnames(Fig7F_MekiDR_Mek_data))

############ pRaf values #######
PDDR_cRaf_FCoXX_Ttest <- PDDR_cRaf_FCoXX_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  mutate(Treatment=as.numeric(Treatment_value)) %>% 
  select(-c(Treatment_value))

PDDR_Raf_ReplicateValues  = PDDR_Raf_FCoXX %>% #colnames()
  select(Cell_line,Replicate,Treatment,log2Treatment,Analyte,Signal) %>% 
  pivot_wider(names_from = Replicate,values_from=Signal) %>% 
  rowwise() %>% 
  mutate(Mean=mean(c(R1,R2,R3), na.rm = TRUE))

Fig7F_MekiDR_Raf_data = left_join(PDDR_Raf_ReplicateValues,PDDR_cRaf_FCoXX_Ttest, by="Treatment") 
Fig7F_MekiDR_Raf_data = Fig7F_MekiDR_Raf_data %>% 
  mutate(XX_mean = ifelse(Cell_line == "XO", NA, XX_mean)) %>% 
  mutate(XO_mean = ifelse(Cell_line == "XO", NA, XO_mean)) %>% 
  mutate(p_value = ifelse(Cell_line == "XO", NA, p_value)) %>% 
  mutate(sig = ifelse(Cell_line == "XO", NA, sig))
colnames(Fig7F_MekiDR_Raf_data) = gsub("p_value","p_value(XXvsXO)",colnames(Fig7F_MekiDR_Raf_data))



Fig7F_MekiDR_All_data = dplyr::bind_rows(Fig7F_MekiDR_Mek_data,Fig7F_MekiDR_Raf_data)
Fig7F_MekiDR_All_data$Analyte = gsub("FCoXX_M2_", "",Fig7F_MekiDR_All_data$Analyte)
WriteXLS::WriteXLS(Fig7F_MekiDR_All_data, "./OUTPUT_PAPER/Fig7F_PDDR_All_data.xls")


print("Script1 : MekiDR - All steps executed ")
