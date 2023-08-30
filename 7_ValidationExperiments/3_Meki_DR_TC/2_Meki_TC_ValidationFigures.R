#!/usr/bin/env Rscript


library(dplyr)
library(stringr)
library(WriteXLS)
library(readxl) #for read_excel
library(ggplot2)
library(tidyverse) 
library(egg) #needed for set_panel_size

source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

###############  Input Data  #################

Analytes_quant_folder= "../../RAW_DATA/ValidationExperiments/MEKi_TC_DR/Analyte_Quant/"
TPS_quant_folder = "../../RAW_DATA/ValidationExperiments/MEKi_TC_DR/TPS_Quant/"

dir.create("./OUTPUT_MekiTC") # Create folder to all save output data

# Create folder OUTPUT_PAPER to save figs used in paper
if(!(file.exists("./OUTPUT_PAPER"))){ 
  dir.create("./OUTPUT_PAPER") 
}


########## Reading in Data : pMek ###########

R1_G2_TPS <- read_excel(file.path(TPS_quant_folder,"PDTC_R1_G2_TPS.xlsx"))
R1_G2_pMek <- read_excel(file.path(Analytes_quant_folder,"PDTC_R1_G2_pMek_Quant.xls"))
R1_G2_pMek_returned <- prep_analyte_G12(R1_G2_TPS, R1_G2_pMek)

R2_G2_TPS <- read_excel(file.path(TPS_quant_folder,"PDTC_R2_G2_TPS.xlsx"))
R2_G2_pMek <- read_excel(file.path(Analytes_quant_folder,"PDTC_R2_G2_pMek_Quant.xls"))
R2_G2_pMek_returned <- prep_analyte_G12(R2_G2_TPS, R2_G2_pMek)

R3_G2_TPS <- read_excel(file.path(TPS_quant_folder,"PDTC_R3_G2_TPS.xlsx"))
R3_G2_pMek <- read_excel(file.path(Analytes_quant_folder,"PDTC_R3_G2_pMek_Quant.xls"))
R3_G2_pMek_returned <- prep_analyte_G12(R3_G2_TPS, R3_G2_pMek)

pMek_XXXO <- rbind(R1_G2_pMek_returned, R2_G2_pMek_returned, R3_G2_pMek_returned)

### Normalization over Gel Mean is being done in the "prep_analyte_G12" function. 
### The normalized columns need renaming to become compatible with the function that was written later for all the dose response experiments.
colnames(pMek_XXXO) <- gsub("Signal_TNorm_M1","Norm_phosP_by_total_M1",colnames(pMek_XXXO))
colnames(pMek_XXXO) <- gsub("Signal_TNorm_M2","Norm_phosP_by_total_M2",colnames(pMek_XXXO))
pMek_XXXO_FC <- FC_Calculation_updated(pMek_XXXO) 


########## Reading in Data : pcRaf ###########

R1_G2_TPS <- read_excel(file.path(TPS_quant_folder,"PDTC_R1_G2_TPS.xlsx"))
R1_G2_pcRaf <- read_excel(file.path(Analytes_quant_folder,"PDTC_R1_G2_pcRaf_Quant.xls"))
R1_G2_pcRaf_returned <- prep_analyte_G12(R1_G2_TPS, R1_G2_pcRaf)

R2_G2_TPS <- read_excel(file.path(TPS_quant_folder,"PDTC_R2_G2_TPS.xlsx"))
R2_G2_pcRaf <- read_excel(file.path(Analytes_quant_folder,"PDTC_R2_G2_pcRaf_Quant.xls"))
R2_G2_pcRaf_returned <- prep_analyte_G12(R2_G2_TPS, R2_G2_pcRaf)

R3_G2_TPS <- read_excel(file.path(TPS_quant_folder,"PDTC_R3_G2_TPS.xlsx"))
R3_G2_pcRaf <- read_excel(file.path(Analytes_quant_folder,"PDTC_R3_G2_pcRaf_Quant.xls"))
R3_G2_pcRaf_returned <- prep_analyte_G12(R3_G2_TPS, R3_G2_pcRaf)

pcRaf_XXXO <- rbind(R1_G2_pcRaf_returned, R2_G2_pcRaf_returned, R3_G2_pcRaf_returned)

### Normalization over Gel Mean is being done in the "prep_analyte_G12" function. 
### The normalized columns need renaming to become compatible with the function that was written later for all the dose response experiments.
colnames(pcRaf_XXXO) <- gsub("Signal_TNorm_M1","Norm_phosP_by_total_M1",colnames(pcRaf_XXXO))
colnames(pcRaf_XXXO) <- gsub("Signal_TNorm_M2","Norm_phosP_by_total_M2",colnames(pcRaf_XXXO))

pcRaf_XXXO_FC <- FC_Calculation_updated(pcRaf_XXXO) 


############### PAPER FIG ###############
##### Fold change over XX control ######

pMek_XXXO_plot <- pMek_XXXO_FC %>% 
  select(!(grep("M1",colnames(pMek_XXXO_FC)))) %>% 
  select(-c("Image_Name", "Channel", "Well", "Signal", "TPS", "Analyte", "Samples"))

pMek_XXXO_plot <- pMek_XXXO_plot %>% 
  dplyr::rename(pMek_Signal_by_TPS = Signal_by_TPS,
                pMek_Signal_TNorm_M2 = Norm_phosP_by_total_M2, 
                pMek_FCoXX_M2 = FCoXX_M2,
                pMek_FCoCntrl_M2 = FCoCntrl_M2) 

pcRaf_XXXO_plot <- pcRaf_XXXO_FC %>% 
  select(!(grep("M1",colnames(pcRaf_XXXO_FC)))) %>% 
  select(-c("Image_Name", "Channel", "Well", "Signal", "TPS", "Analyte", "Samples"))

pcRaf_XXXO_plot <- pcRaf_XXXO_plot %>% 
  dplyr::rename(pcRaf_Signal_by_TPS = Signal_by_TPS,
                pcRaf_Signal_TNorm_M2 = Norm_phosP_by_total_M2,
                pcRaf_FCoXX_M2 = FCoXX_M2,
                pcRaf_FCoCntrl_M2 = FCoCntrl_M2) 


PDTC_pMek_pcRaf_Join <- dplyr::full_join(pMek_XXXO_plot, pcRaf_XXXO_plot, by=c("Cell_line","Treatment","Replicate", "Treatment_numeric"))

PDTC_pMek_pcRaf_Joinplot <- PDTC_pMek_pcRaf_Join %>% 
  select(c("Cell_line","Treatment","Replicate", "Treatment_numeric",grep("FCoXX",colnames(PDTC_pMek_pcRaf_Join)))) %>% 
  gather(key = "Analyte", value = "Signal", -c("Cell_line","Treatment","Replicate", "Treatment_numeric"))

PDTC_pMek_pcRaf_Joinplot$Treatment_factor <- factor(PDTC_pMek_pcRaf_Joinplot$Treatment_numeric)

PDTC_pMek_pcRaf_Joinplot$Analyte <- factor(PDTC_pMek_pcRaf_Joinplot$Analyte, levels = c("pMek_FCoXX_M2","pcRaf_FCoXX_M2")) 
Analyte_labels <- c("pMek_FCoXX_M2" = "pMEK", "pcRaf_FCoXX_M2" = "pRAF1")

##############Save data file  ############
write.csv(PDTC_pMek_pcRaf_Joinplot, file = "OUTPUT_MekiTC/PDTC_pMek_pcRaf_Joinplot.csv")


##### These calculation are done to find the max y-values for the 2 panels. Inspite of using "free_y", I wanted to add some ##### 
#### space on top of the plot to accomodate the significance stars ######## 
#### The max of y-values for each of the panels is calculated and plotted on top of the plot by specifiying them to ##### 
#### the corresponding analyte panel ##### 


pMek_max <- PDTC_pMek_pcRaf_Joinplot %>%
  filter(grepl("Mek", Analyte))
pMek_max <-max(pMek_max$Signal) + 1 # the number is based on how much space i need

pcRaf_max <- PDTC_pMek_pcRaf_Joinplot %>%
  filter(grepl("Raf", Analyte))
pcRaf_max <-max(pcRaf_max$Signal) + 0.2 # the number is based on how much space i need

dummy_Mek <- data.frame(Treatment_factor = 0, Signal = pMek_max,
                        Analyte = "pMek_FCoXX_M2", Cell_line = "Common", stringsAsFactors=FALSE)
dummy_Raf <- data.frame(Treatment_factor = 0, Signal = pcRaf_max,
                        Analyte = "pcRaf_FCoXX_M2", Cell_line = "Common", stringsAsFactors=FALSE) 
dummy_data <- rbind(dummy_Mek,dummy_Raf)
dummy_data$Analyte = factor(dummy_data$Analyte, levels =c("pMek_FCoXX_M2", "pcRaf_FCoXX_M2")) # To get the correct sequence in the Facet_wrap

#################### PLOT ###########

g <- Plot_TwoPanel_ValidationPlot_updated(CombinedData=PDTC_pMek_pcRaf_Joinplot, x_axis_entity="Treatment_factor", y_axis_entity="Signal", Analyte_labels=Analyte_labels,scale_choice="free")+
  geom_blank(data=dummy_data) + 
  labs(x = "\n Time(h) ",
       y = "Rel. phos. (norm.)\n",
       color = "Cell line" ) + # color within labs for user defined labels to the attribute in legend
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")


gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("PDTC_pMek_pcRaf_FCoXX.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_MekiTC") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that
ggsave("Fig7FGi_PDTC_pMek_pcRaf_FCoXX.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that



##### Fold change over respective cell line control ######

PDTC_pMek_pcRaf_Joinplot_OverCntrl <- PDTC_pMek_pcRaf_Join %>% 
  select(c("Cell_line","Treatment","Replicate", "Treatment_numeric",grep("FCoCntrl",colnames(PDTC_pMek_pcRaf_Join)))) %>% 
  gather(key = "Analyte", value = "Signal", -c("Cell_line","Treatment","Replicate", "Treatment_numeric"))

PDTC_pMek_pcRaf_Joinplot_OverCntrl$Treatment_factor <- factor(PDTC_pMek_pcRaf_Joinplot_OverCntrl$Treatment_numeric)

PDTC_pMek_pcRaf_Joinplot_OverCntrl$Analyte <- factor(PDTC_pMek_pcRaf_Joinplot_OverCntrl$Analyte, levels = c("pMek_FCoCntrl_M2","pcRaf_FCoCntrl_M2")) ## This is needed in order to plot pMek before pcRaf in the facet_wrap
Analyte_labels <- c("pMek_FCoCntrl_M2" = "pMEK", "pcRaf_FCoCntrl_M2" = "pRAF1")



##############Save data file ############
write.csv(PDTC_pMek_pcRaf_Joinplot_OverCntrl, file = "./OUTPUT_MekiTC/PDTC_pMek_pcRaf_Joinplot_OverCntrl.csv")

##### These calculation are done to find the max y-values for the 2 panels. Inspite of using "free_y", I wanted to add some ##### 
#### space on top of the plot to accomodate the significance stars ######## 
#### The max of y-values for each of the panels is calculated and plotted on top of the plot by specifiying them to ##### 
#### the corresponding analyte panel ##### 


pMek_max <- PDTC_pMek_pcRaf_Joinplot_OverCntrl %>%
  filter(grepl("Mek", Analyte))
pMek_max <-max(pMek_max$Signal) + 7 # the number is based on how much space i need

pcRaf_max <- PDTC_pMek_pcRaf_Joinplot_OverCntrl %>%
  filter(grepl("Raf", Analyte))
pcRaf_max <-max(pcRaf_max$Signal) + 0.2 # the number is based on how much space i need

dummy_Mek <- data.frame(Treatment_factor = 0, Signal = pMek_max,
                        Analyte = "pMek_FCoCntrl_M2", Cell_line = "Common", stringsAsFactors=FALSE)
dummy_Raf <- data.frame(Treatment_factor = 0, Signal = pcRaf_max,
                        Analyte = "pcRaf_FCoCntrl_M2", Cell_line = "Common", stringsAsFactors=FALSE) 
dummy_data <- rbind(dummy_Mek,dummy_Raf)
dummy_data$Analyte = factor(dummy_data$Analyte, levels =c("pMek_FCoCntrl_M2", "pcRaf_FCoCntrl_M2")) # This is to get the correct sequence in the Facet_wrap

#################### PLOT ###########


g <- Plot_TwoPanel_ValidationPlot_updated(CombinedData=PDTC_pMek_pcRaf_Joinplot_OverCntrl, x_axis_entity="Treatment_factor", y_axis_entity="Signal", Analyte_labels=Analyte_labels,scale_choice="free")+
  geom_blank(data=dummy_data) + 
  labs(x = "\n Time(h) ",
       y = "Rel. phosp. (norm.)\n Fold change over untreated\n",
       color = "Cell line" ) + # color within labs for user defined labels to the attribute in legend
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("PDTC_pMek_pcRaf_FCoCntrl.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_MekiTC") 
ggsave("Fig7FGii_PDTC_pMek_pcRaf_FCoCntrl.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER") 


########### T Tests for difference between XX and XO in case of FCoXX ########

Meki_TimePoints = c(0, 0.25, 0.50, 1,  2,  4, 24, 48)

PDTC_XX_XO_Mek_pvals_FCoXX = list()
for(i in 1:length(Meki_TimePoints) ) {

  
  PDTC_XX_Mek_values <- PDTC_pMek_pcRaf_Joinplot %>% 
    ungroup() %>% 
    select(-c(Treatment_factor)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "pMek_FCoXX_M2") %>% 
    filter(Treatment_numeric == Meki_TimePoints[i])
  
  PDTC_XO_Mek_values <- PDTC_pMek_pcRaf_Joinplot %>% 
    ungroup() %>% 
    select(-c(Treatment_factor)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "pMek_FCoXX_M2") %>% 
    filter(Treatment_numeric == Meki_TimePoints[i])
  
  PDTC_XX_XO_Mek_pvals_FCoXX$Treatment_value[i] = Meki_TimePoints[i]
  PDTC_XX_XO_Mek_pvals_FCoXX$XX_mean[i] = t.test(PDTC_XX_Mek_values$Signal,PDTC_XO_Mek_values$Signal, paired = FALSE)$estimate[1]
  PDTC_XX_XO_Mek_pvals_FCoXX$XO_mean[i] = t.test(PDTC_XX_Mek_values$Signal,PDTC_XO_Mek_values$Signal, paired = FALSE)$estimate[2]
  PDTC_XX_XO_Mek_pvals_FCoXX$p_value[i] = t.test(PDTC_XX_Mek_values$Signal,PDTC_XO_Mek_values$Signal, paired = FALSE)$p.value
  
}

PDTC_Mek_FCoXX_Ttest <- data.frame((sapply(PDTC_XX_XO_Mek_pvals_FCoXX,c)))
# PDTC_Mek_FCoXX_Ttest <- PDTC_Mek_FCoXX_Ttest %>% 
#   mutate(sig = ifelse(p_value<0.05,"*",""))
# WriteXLS::WriteXLS(PDTC_Mek_FCoXX_Ttest, "./OUTPUT_MekiTC/PDTC_Mek_FCoXX_Ttest.xls")


PDTC_XX_XO_cRaf_pvals_FCoXX = list()
for(i in 1:length(Meki_TimePoints) ) {
  
  PDTC_XX_cRaf_values <- PDTC_pMek_pcRaf_Joinplot %>% 
    ungroup() %>% 
    select(-c(Treatment_factor)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "pcRaf_FCoXX_M2") %>% 
    filter(Treatment_numeric == Meki_TimePoints[i])
  
  PDTC_XO_cRaf_values <- PDTC_pMek_pcRaf_Joinplot %>% 
    ungroup() %>% 
    select(-c(Treatment_factor)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "pcRaf_FCoXX_M2") %>% 
    filter(Treatment_numeric == Meki_TimePoints[i])
  
  PDTC_XX_XO_cRaf_pvals_FCoXX$Treatment_value[i] = Meki_TimePoints[i]
  PDTC_XX_XO_cRaf_pvals_FCoXX$XX_mean[i] = t.test(PDTC_XX_cRaf_values$Signal,PDTC_XO_cRaf_values$Signal, paired = FALSE)$estimate[1]
  PDTC_XX_XO_cRaf_pvals_FCoXX$XO_mean[i] = t.test(PDTC_XX_cRaf_values$Signal,PDTC_XO_cRaf_values$Signal, paired = FALSE)$estimate[2]
  PDTC_XX_XO_cRaf_pvals_FCoXX$p_value[i] = t.test(PDTC_XX_cRaf_values$Signal,PDTC_XO_cRaf_values$Signal, paired = FALSE)$p.value
  
}

PDTC_cRaf_FCoXX_Ttest <- data.frame((sapply(PDTC_XX_XO_cRaf_pvals_FCoXX,c)))
# PDTC_cRaf_FCoXX_Ttest <- PDTC_cRaf_FCoXX_Ttest %>% 
#   mutate(sig = ifelse(p_value<0.05,"*",""))
# WriteXLS::WriteXLS(PDTC_cRaf_FCoXX_Ttest, "./OUTPUT_MekiTC/PDTC_cRaf_FCoXX_Ttest.xls")


########### T Tests for difference between XX and XO in case of FCoCntrl ########

Meki_TimePoints = c(0, 0.25, 0.50, 1,  2,  4, 24, 48)

PDTC_XX_XO_Mek_pvals_FCoCntrl = list()
for(i in 1:length(Meki_TimePoints) ) {

  
  PDTC_XX_Mek_values <- PDTC_pMek_pcRaf_Joinplot_OverCntrl %>% 
    ungroup() %>% 
    select(-c(Treatment_factor)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "pMek_FCoCntrl_M2") %>% 
    filter(Treatment_numeric == Meki_TimePoints[i])
  
  PDTC_XO_Mek_values <- PDTC_pMek_pcRaf_Joinplot_OverCntrl %>% 
    ungroup() %>% 
    select(-c(Treatment_factor)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "pMek_FCoCntrl_M2") %>% 
    filter(Treatment_numeric == Meki_TimePoints[i])
  
  PDTC_XX_XO_Mek_pvals_FCoCntrl$Treatment_value[i] = Meki_TimePoints[i]
  PDTC_XX_XO_Mek_pvals_FCoCntrl$XX_mean[i] = t.test(PDTC_XX_Mek_values$Signal,PDTC_XO_Mek_values$Signal, paired = FALSE)$estimate[1]
  PDTC_XX_XO_Mek_pvals_FCoCntrl$XO_mean[i] = t.test(PDTC_XX_Mek_values$Signal,PDTC_XO_Mek_values$Signal, paired = FALSE)$estimate[2]
  PDTC_XX_XO_Mek_pvals_FCoCntrl$p_value[i] = t.test(PDTC_XX_Mek_values$Signal,PDTC_XO_Mek_values$Signal, paired = FALSE)$p.value
  
}

PDTC_Mek_FCoCntrl_Ttest <- data.frame((sapply(PDTC_XX_XO_Mek_pvals_FCoCntrl,c)))
# PDTC_Mek_FCoCntrl_Ttest <- PDTC_Mek_FCoCntrl_Ttest %>% 
#   mutate(sig = ifelse(p_value<0.05,"*",""))
# WriteXLS::WriteXLS(PDTC_Mek_FCoCntrl_Ttest, "./OUTPUT_MekiTC/PDTC_Mek_FCoCntrl_Ttest.xls")


PDTC_XX_XO_cRaf_pvals_FCoCntrl = list()
for(i in 1:length(Meki_TimePoints) ) {
  
  PDTC_XX_cRaf_values <- PDTC_pMek_pcRaf_Joinplot_OverCntrl %>% 
    ungroup() %>% 
    select(-c(Treatment_factor)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "pcRaf_FCoCntrl_M2") %>% 
    filter(Treatment_numeric == Meki_TimePoints[i])
  
  PDTC_XO_cRaf_values <- PDTC_pMek_pcRaf_Joinplot_OverCntrl %>% 
    ungroup() %>% 
    select(-c(Treatment_factor)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "pcRaf_FCoCntrl_M2") %>% 
    filter(Treatment_numeric == Meki_TimePoints[i])
  
  PDTC_XX_XO_cRaf_pvals_FCoCntrl$Treatment_value[i] = Meki_TimePoints[i]
  PDTC_XX_XO_cRaf_pvals_FCoCntrl$XX_mean[i] = t.test(PDTC_XX_cRaf_values$Signal,PDTC_XO_cRaf_values$Signal, paired = FALSE)$estimate[1]
  PDTC_XX_XO_cRaf_pvals_FCoCntrl$XO_mean[i] = t.test(PDTC_XX_cRaf_values$Signal,PDTC_XO_cRaf_values$Signal, paired = FALSE)$estimate[2]
  PDTC_XX_XO_cRaf_pvals_FCoCntrl$p_value[i] = t.test(PDTC_XX_cRaf_values$Signal,PDTC_XO_cRaf_values$Signal, paired = FALSE)$p.value
  
}

PDTC_cRaf_FCoCntrl_Ttest <- data.frame((sapply(PDTC_XX_XO_cRaf_pvals_FCoCntrl,c)))
# PDTC_cRaf_FCoCntrl_Ttest <- PDTC_cRaf_FCoCntrl_Ttest %>% 
#   mutate(sig = ifelse(p_value<0.05,"*",""))
# WriteXLS::WriteXLS(PDTC_cRaf_FCoCntrl_Ttest, "./OUTPUT_MekiTC/PDTC_cRaf_FCoCntrl_Ttest.xls")



########## Saving the data for fig 7 F and G : fold change over untreated XX(FCoXX) #########
PDTC_Mek_FCoXX = PDTC_pMek_pcRaf_Joinplot %>% 
  filter(grepl("Mek", Analyte))

PDTC_Raf_FCoXX = PDTC_pMek_pcRaf_Joinplot %>% 
  filter(grepl("Raf", Analyte))


############ pMek values #######
PDTC_Mek_FCoXX_Ttest <- PDTC_Mek_FCoXX_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  mutate(Treatment=as.numeric(Treatment_value)) %>% 
  select(-c(Treatment_value))

PDTC_Mek_ReplicateValues  = PDTC_Mek_FCoXX %>% #colnames()
  select(Cell_line,Replicate,Treatment,Treatment_factor,Analyte,Signal) %>% 
  mutate(Treatment=as.numeric(as.character(Treatment_factor))) %>% 
  pivot_wider(names_from = Replicate,values_from=Signal) %>% 
  rowwise() %>% 
  mutate(Mean=mean(c(R1,R2,R3), na.rm = TRUE))

Fig7FG_MekiTC_Mek_data = left_join(PDTC_Mek_ReplicateValues,PDTC_Mek_FCoXX_Ttest, by="Treatment") 
Fig7FG_MekiTC_Mek_data = Fig7FG_MekiTC_Mek_data %>% 
  mutate(XX_mean = ifelse(Cell_line == "XO", NA, XX_mean)) %>% 
  mutate(XO_mean = ifelse(Cell_line == "XO", NA, XO_mean)) %>% 
  mutate(p_value = ifelse(Cell_line == "XO", NA, p_value)) %>% 
  mutate(sig = ifelse(Cell_line == "XO", NA, sig))
colnames(Fig7FG_MekiTC_Mek_data) = gsub("p_value","p_value(XXvsXO)",colnames(Fig7FG_MekiTC_Mek_data))

############ pRaf values #######
PDTC_cRaf_FCoXX_Ttest <- PDTC_cRaf_FCoXX_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  mutate(Treatment=as.numeric(Treatment_value)) %>% 
  select(-c(Treatment_value))

PDTC_Raf_ReplicateValues  = PDTC_Raf_FCoXX %>% #colnames()
  select(Cell_line,Replicate,Treatment,Treatment_factor,Analyte,Signal) %>% 
  mutate(Treatment=as.numeric(as.character(Treatment_factor))) %>%
  pivot_wider(names_from = Replicate,values_from=Signal) %>% 
  rowwise() %>% 
  mutate(Mean=mean(c(R1,R2,R3), na.rm = TRUE))

Fig7FG_MekiTC_Raf_data = left_join(PDTC_Raf_ReplicateValues,PDTC_cRaf_FCoXX_Ttest, by="Treatment") 
Fig7FG_MekiTC_Raf_data = Fig7FG_MekiTC_Raf_data %>% 
  mutate(XX_mean = ifelse(Cell_line == "XO", NA, XX_mean)) %>% 
  mutate(XO_mean = ifelse(Cell_line == "XO", NA, XO_mean)) %>% 
  mutate(p_value = ifelse(Cell_line == "XO", NA, p_value)) %>% 
  mutate(sig = ifelse(Cell_line == "XO", NA, sig))
colnames(Fig7FG_MekiTC_Raf_data) = gsub("p_value","p_value(XXvsXO)",colnames(Fig7FG_MekiTC_Raf_data))

Fig7FGi_PDTC_Mek_Raf_FCoXX_dat = dplyr::bind_rows(Fig7FG_MekiTC_Mek_data,Fig7FG_MekiTC_Raf_data)
Fig7FGi_PDTC_Mek_Raf_FCoXX_dat$Analyte = gsub("_FCoXX_M2", "",Fig7FGi_PDTC_Mek_Raf_FCoXX_dat$Analyte)
WriteXLS::WriteXLS(Fig7FGi_PDTC_Mek_Raf_FCoXX_dat, "./OUTPUT_PAPER/Fig7FGi_PDTC_Mek_Raf_FCoXX_dat.xls")

########## Saving the data for fig 7 F and G : fold change over respective untreated cellline(FCoC) #########
PDTC_Mek_FCoCntrl = PDTC_pMek_pcRaf_Joinplot_OverCntrl %>% 
  filter(grepl("Mek", Analyte))

PDTC_Raf_FCoCntrl = PDTC_pMek_pcRaf_Joinplot_OverCntrl %>% 
  filter(grepl("Raf", Analyte))


############ pMek values #######
PDTC_Mek_FCoCntrl_Ttest <- PDTC_Mek_FCoCntrl_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  mutate(Treatment=as.numeric(Treatment_value)) %>% 
  select(-c(Treatment_value))

PDTC_Mek_ReplicateValues  = PDTC_Mek_FCoCntrl %>% #colnames()
  select(Cell_line,Replicate,Treatment,Treatment_factor,Analyte,Signal) %>% 
  mutate(Treatment=as.numeric(as.character(Treatment_factor))) %>% 
  pivot_wider(names_from = Replicate,values_from=Signal) %>% 
  rowwise() %>% 
  mutate(Mean=mean(c(R1,R2,R3), na.rm = TRUE))

Fig7FG_MekiTC_Mek_data = left_join(PDTC_Mek_ReplicateValues,PDTC_Mek_FCoCntrl_Ttest, by="Treatment") 
Fig7FG_MekiTC_Mek_data = Fig7FG_MekiTC_Mek_data %>% 
  mutate(XX_mean = ifelse(Cell_line == "XO", NA, XX_mean)) %>% 
  mutate(XO_mean = ifelse(Cell_line == "XO", NA, XO_mean)) %>% 
  mutate(p_value = ifelse(Cell_line == "XO", NA, p_value)) %>% 
  mutate(sig = ifelse(Cell_line == "XO", NA, sig))
colnames(Fig7FG_MekiTC_Mek_data) = gsub("p_value","p_value(XXvsXO)",colnames(Fig7FG_MekiTC_Mek_data))

############ pRaf values #######
PDTC_cRaf_FCoCntrl_Ttest <- PDTC_cRaf_FCoCntrl_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  mutate(Treatment=as.numeric(Treatment_value)) %>% 
  select(-c(Treatment_value))

PDTC_Raf_ReplicateValues  = PDTC_Raf_FCoCntrl %>% #colnames()
  select(Cell_line,Replicate,Treatment,Treatment_factor,Analyte,Signal) %>% 
  mutate(Treatment=as.numeric(as.character(Treatment_factor))) %>%
  pivot_wider(names_from = Replicate,values_from=Signal) %>% 
  rowwise() %>% 
  mutate(Mean=mean(c(R1,R2,R3), na.rm = TRUE))

Fig7FG_MekiTC_Raf_data = left_join(PDTC_Raf_ReplicateValues,PDTC_cRaf_FCoCntrl_Ttest, by="Treatment") 
Fig7FG_MekiTC_Raf_data = Fig7FG_MekiTC_Raf_data %>% 
  mutate(XX_mean = ifelse(Cell_line == "XO", NA, XX_mean)) %>% 
  mutate(XO_mean = ifelse(Cell_line == "XO", NA, XO_mean)) %>% 
  mutate(p_value = ifelse(Cell_line == "XO", NA, p_value)) %>% 
  mutate(sig = ifelse(Cell_line == "XO", NA, sig))
colnames(Fig7FG_MekiTC_Raf_data) = gsub("p_value","p_value(XXvsXO)",colnames(Fig7FG_MekiTC_Raf_data))

Fig7FGii_PDTC_Mek_Raf_FCoC_dat = dplyr::bind_rows(Fig7FG_MekiTC_Mek_data,Fig7FG_MekiTC_Raf_data)
Fig7FGii_PDTC_Mek_Raf_FCoC_dat$Analyte = gsub("_FCoCntrl_M2", "",Fig7FGii_PDTC_Mek_Raf_FCoC_dat$Analyte)
WriteXLS::WriteXLS(Fig7FGii_PDTC_Mek_Raf_FCoC_dat, "./OUTPUT_PAPER/Fig7FGii_PDTC_Mek_Raf_FCoC_dat.xls")




print("Script2 : MekiTC - All steps executed ")







