#!/usr/bin/env Rscript


## Libraries
library(readxl) # Needed for read_excel.Part of tidyverse, but needs to be installed separately.(Working on Mac but not my Linux - Check !)
library(tidyr)
library(ggplot2)
library(ggpubr) # needed for compare_means and some other functions
library(ggrepel) 
library(egg) #needed for set_panel_size, to get predefined sizes of figure panels.
library(plyr)
library(dplyr)


source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")


Gsk3inh_CT_Doses =c("0","1.5","3","6","12","24")


quant_folder= "../../RAW_DATA/ValidationExperiments/GSK3i_DR/Analyte_Quant/"
empiria_folder = "../../RAW_DATA/ValidationExperiments/GSK3i_DR/TPS_Quant/"


dir.create("./OUTPUT_PAPER") # Create folder to save OUTPUT data

################# pAkt/ Akt #################

## Reading in the data for the 3 replicates 

#########  Gel1 : R1   ######### 

Gel1_samples = c(
  "XX_R1_cntrl_CTDR",
  "XX_R1_1.5_CTDR",
  "XX_R1_3_CTDR",
  "XX_R1_6_CTDR",
  "XX_R1_12_CTDR",
  "XX_R1_24_CTDR",
  "XO_R1_cntrl_CTDR",
  "XO_R1_1.5_CTDR",
  "XO_R1_3_CTDR",
  "XO_R1_6_CTDR",
  "XO_R1_12_CTDR",
  "XO_R1_24_CTDR"
)

Gel1_Akt_quant_file = c("0001720_02_JuneGel1_GSK3iDR_Aktp.xls")
Gel1_AktDB <-  read_excel(file.path(quant_folder,Gel1_Akt_quant_file))

Gel1_AktDB_labeled <- ReadInWBData(Gel1_samples, Gel1_AktDB, Gsk3inh_CT_Doses, "Gel1", "Akt")
#might give an error here if plyr is loaded after dplyr

Gel1_Akt_NormOMeanRep <- Norm_o_MeanRep(Gel1_AktDB_labeled)

#########  Gel22 : R2 + R1XX(dropped)   #########  

Gel22_samples = c(
  "XO_R1_cntrl_CTDR",
  "XX_R1_cntrl_CTDR",
  "XX_R1_1.5_CTDR",
  "XX_R1_3_CTDR",
  "XX_R1_6_CTDR",
  "XX_R1_12_CTDR",
  "XX_R1_24_CTDR",
  "XO_R2_cntrl_CTDR",
  "XO_R2_1.5_CTDR",
  "XO_R2_3_CTDR",
  "XO_R2_6_CTDR",
  "XO_R2_12_CTDR",
  "XO_R2_24_CTDR",
  "XX_R2_cntrl_CTDR",
  "XX_R2_1.5_CTDR",
  "XX_R2_3_CTDR",
  "XX_R2_6_CTDR",
  "XX_R2_12_CTDR",
  "XX_R2_24_CTDR")

Gel22_Akt_quant_file = c("0000062_02_JuneGel22_GSK3iDR_Aktp.xls")
Gel22_AktDB <-  read_excel(file.path(quant_folder,Gel22_Akt_quant_file))

Gel22_AktDB_labeled <- ReadInWBData(Gel22_samples, Gel22_AktDB, Gsk3inh_CT_Doses, "Gel22", "Akt")
#Gel22_AktDB_NormOCom <- Norm_o_ComSample(Gel22_AktDB_labeled, ComSample_Exp = "CTDR", ComSample_Rep = "R1", ComSample_T = 0, ComSample_CellLine = c("XO") ) 

# Removing replicate 1, now that I have a new dataset fro that replicate
Gel22_AktDB_labeled <- Gel22_AktDB_labeled %>% 
  filter(Replicate != "R1")


Gel22_Akt_NormOMeanRep <- Norm_o_MeanRep(Gel22_AktDB_labeled)

#########  Gel23 : R3 + R1XO(dropped)    #########  

Gel23_samples = c(
  "XO_R1_cntrl_CTDR",
  "XO_R1_3_CTDR",
  "XO_R1_6_CTDR",
  "XO_R1_12_CTDR",
  "XO_R1_24_CTDR",
  "XX_R3_cntrl_CTDR",
  "XX_R3_1.5_CTDR",
  "XX_R3_3_CTDR",
  "XX_R3_6_CTDR",
  "XX_R3_12_CTDR",
  "XX_R3_24_CTDR",
  "XO_R3_cntrl_CTDR",
  "XO_R3_1.5_CTDR",
  "XO_R3_3_CTDR",
  "XO_R3_6_CTDR",
  "XO_R3_12_CTDR",
  "XO_R3_24_CTDR")
Gel23_Akt_quant_file = c("0000063_02_JuneGel23_GSK3iDR_Aktp.xls")
Gel23_AktDB <-  read_excel(file.path(quant_folder,Gel23_Akt_quant_file))

Gel23_AktDB_labeled <- ReadInWBData(Gel23_samples, Gel23_AktDB, Gsk3inh_CT_Doses, "Gel23", "Akt")
#Gel23_AktDB_NormOCom <- Norm_o_ComSample(Gel23_AktDB_labeled, ComSample_Exp = "CTDR", ComSample_Rep = "R1", ComSample_T = 0, ComSample_CellLine = c("XO") ) 

Gel23_AktDB_labeled <- Gel23_AktDB_labeled %>% 
  filter(Replicate != "R1")

Gel23_Akt_NormOMeanRep <- Norm_o_MeanRep(Gel23_AktDB_labeled)

#########  Collate the 3 REPLICATES #####

CTDR_Akt_NormOMeanRep <- dplyr::bind_rows(Gel1_Akt_NormOMeanRep,Gel22_Akt_NormOMeanRep, Gel23_Akt_NormOMeanRep) 

######## Fold change Over XX and Over resp Cntrl ##############

CTDR_Akt_FC <- FC_Calculation_updated(CTDR_Akt_NormOMeanRep)

CTDR_Akt_plot <- CTDR_Akt_FC %>% 
  select(-c(Exp,Gel)) # You can remove whichever columns are not needed


################# pP70S6k / TPS #################

## Reading in the data for the 3 replicates 

#########  Gel1 : R1re  ######### 

Gel3_samples = c(
  "XX_R1_cntrl_CTDR",
  "XX_R1_1.5_CTDR",
  "XX_R1_3_CTDR",
  "XX_R1_6_CTDR",
  "XX_R1_12_CTDR",
  "XX_R1_24_CTDR",
  "XO_R1_cntrl_CTDR",
  "XO_R1_1.5_CTDR",
  "XO_R1_3_CTDR",
  "XO_R1_6_CTDR",
  "XO_R1_12_CTDR",
  "XO_R1_24_CTDR"
)

Gel3_P70s6k_quant_file = c("0001725_02_JuneGel3_GSK3iDR_pP70S6k.xls")
Gel3_P70s6kDB <-  read_excel(file.path(quant_folder,Gel3_P70s6k_quant_file))

# Gel3_TPS_file = "June2020_Gsk3iDR_R1trial_Gel3.xlsx"
# Gel3_TPS <- read_excel(file.path(empiria_folder,Gel3_TPS_file))

Gel3_TPS_file = "0000584_02_June2020_Gsk3iDR_R1trial_Gel3.xlsx"
Gel3_TPS <- read_excel(file.path(empiria_folder,Gel3_TPS_file))


Gel3_P70s6kDB_labeled <- ReadInWBData_TPS(Gel3_samples, Gel3_P70s6kDB, Gel3_TPS,Gsk3inh_CT_Doses, "Gel3", "P70s6k")
#Gel3_P70s6kDB_NormOCom <- Norm_o_ComSample(Gel3_P70s6kDB_labeled, ComSample_Exp = "CTDR", ComSample_Rep = "R1", ComSample_T = 0, ComSample_CellLine = c("XO") ) 

Gel3_P70s6k_NormOMeanRep <- Norm_o_MeanRep(Gel3_P70s6kDB_labeled)

#########  Gel19 : R2 + R1XX(dropped)   #########  

Gel19_samples = c(
  "XO_R1_cntrl_CTDR",
  "XX_R1_cntrl_CTDR",
  "XX_R1_1.5_CTDR",
  "XX_R1_3_CTDR",
  "XX_R1_6_CTDR",
  "XX_R1_12_CTDR",
  "XX_R1_24_CTDR",
  "XO_R2_cntrl_CTDR",
  "XO_R2_1.5_CTDR",
  "XO_R2_3_CTDR",
  "XO_R2_6_CTDR",
  "XO_R2_12_CTDR",
  "XO_R2_24_CTDR",
  "XX_R2_cntrl_CTDR",
  "XX_R2_1.5_CTDR",
  "XX_R2_3_CTDR",
  "XX_R2_6_CTDR",
  "XX_R2_12_CTDR",
  "XX_R2_24_CTDR")

Gel19_P70s6k_quant_file = c("0000044_02_JuneGel19_GSK3bDR_p70Rsk.xls")
Gel19_P70s6kDB <-  read_excel(file.path(quant_folder,Gel19_P70s6k_quant_file))

Gel19_TPS_file = "0000048_02_Gel19_GSK3iDR2_R1XX_TPS.xlsx"
Gel19_TPS <- read_excel(file.path(empiria_folder,Gel19_TPS_file))

Gel19_P70s6kDB_labeled <- ReadInWBData_TPS(Gel19_samples, Gel19_P70s6kDB, Gel19_TPS,Gsk3inh_CT_Doses, "Gel19", "P70s6k")

# Removing replicate 1, now that I have a new dataset for that replicate
Gel19_P70s6kDB_labeled <- Gel19_P70s6kDB_labeled %>% 
  filter(Replicate != "R1")

#Gel19_P70s6kDB_NormOCom <- Norm_o_ComSample(Gel19_P70s6kDB_labeled, ComSample_Exp = "CTDR", ComSample_Rep = "R1", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.

Gel19_P70s6k_NormOMeanRep <- Norm_o_MeanRep(Gel19_P70s6kDB_labeled)

#########  Gel20 : R3 + R1XO(dropped)    #########  

Gel20_samples = c(
  "XX_R1_cntrl_CTDR",
  "XO_R1_cntrl_CTDR",
  "XO_R1_1.5_CTDR",
  "XO_R1_3_CTDR",
  "XO_R1_6_CTDR",
  "XO_R1_12_CTDR",
  "XO_R1_24_CTDR",
  "XX_R3_cntrl_CTDR",
  "XX_R3_1.5_CTDR",
  "XX_R3_3_CTDR",
  "XX_R3_6_CTDR",
  "XX_R3_12_CTDR",
  "XX_R3_24_CTDR",
  "XO_R3_cntrl_CTDR",
  "XO_R3_1.5_CTDR",
  "XO_R3_3_CTDR",
  "XO_R3_6_CTDR",
  "XO_R3_12_CTDR",
  "XO_R3_24_CTDR")

Gel20_P70s6k_quant_file = c("0000045_02_JuneGel20_GSK3bDR_p70Rsk.xls")
Gel20_P70s6kDB <-  read_excel(file.path(quant_folder,Gel20_P70s6k_quant_file))

Gel20_TPS_file = "0000049_02_Gel20_GSK3iDR3_R1XO_TPS.xlsx"
Gel20_TPS <- read_excel(file.path(empiria_folder,Gel20_TPS_file))

Gel20_P70s6kDB_labeled <- ReadInWBData_TPS(Gel20_samples, Gel20_P70s6kDB, Gel20_TPS, Gsk3inh_CT_Doses, "Gel20", "P70s6k")

Gel20_P70s6kDB_labeled <- Gel20_P70s6kDB_labeled %>% 
  filter(Replicate != "R1")

#Gel20_P70s6kDB_NormOCom <- Norm_o_ComSample(Gel20_P70s6kDB_labeled, ComSample_Exp = "CTDR", ComSample_Rep = "R1", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.

Gel20_P70s6k_NormOMeanRep <- Norm_o_MeanRep(Gel20_P70s6kDB_labeled)

#########  Collate the 3 REPLICATES #####

CTDR_P70s6k_NormOMeanRep <- dplyr::bind_rows(Gel3_P70s6k_NormOMeanRep,Gel19_P70s6k_NormOMeanRep,Gel20_P70s6k_NormOMeanRep)

######## Fold change Over XX and Over resp Cntrl ##############

CTDR_P70s6k_FC <- FC_Calculation_updated(CTDR_P70s6k_NormOMeanRep)

CTDR_P70s6k_plot <- CTDR_P70s6k_FC %>% 
  select(-c(Exp,Gel)) # You can remove whichever columns are not needed



################ PLOT FOR THE PAPER ##########


CTDR_All <- dplyr::full_join(CTDR_Akt_plot, CTDR_P70s6k_plot, by=c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr"), suffix = c("_Akt", "_P70s6k"))

CTDR_All <- CTDR_All %>%
  select(-c(grep("Analyte", colnames(CTDR_All))))

CTDR_All$Cell_line <- factor(CTDR_All$Cell_line, levels = c("XX", "XO", "nn"))

###### FC over XX #########

CTDR_Allplot_FCoXX <- CTDR_All %>% 
  select(c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr",grep("FCoXX_M2",colnames(CTDR_All)))) %>% 
  gather(key = "Analyte", value = "Signal", -c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr"))

Analyte_labels <- c("FCoXX_M2_Akt" = "pAKT", "FCoXX_M2_P70s6k" = "pP70S6K")

#### Setting up the plotting attributes ###
# MyCellLineColours <- c("XX" = "#FF0000", "XO" = "#0000CD", "nn" = "#000005")
# MyReplicateShapes <-  c("R1" = 2, "R2" = 0, "R3" = 6, "nn"= 4)


g <- Plot_TwoPanel_ValidationPlot_updated(CTDR_Allplot_FCoXX,"log2Treatment", "Signal" ,Analyte_labels,"fixed")+
  scale_y_continuous(limits = c(0, 1.6)) + ## In this case adding the space for putting stars was easier since the scale was fixed on y also
  labs(x = "\nGSK3i(\u03bcM)+1 [log2]", #\u03bc is the unicode charachter fro greek mu
       #y = "Rel. phosp. (norm.)\n",
       y = " Phosph. rel. to \n untreated XX \n",
       color = "Cell line" )  # color within labs,lets me give user defined labels to the attribute in legend


gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("CTDR_pAkt_pP70s6k_FCoXX.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that



###### FC over Cntrl #########

CTDR_Allplot_FCoCntrl <- CTDR_All %>% 
  select(c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr",grep("FCoCntrl_M2",colnames(CTDR_All)))) %>% 
  gather(key = "Analyte", value = "Signal", -c("Cell_line","Replicate","Treatment","log2Treatment","Treatment_Fctr"))

Analyte_labels <- c("FCoCntrl_M2_Akt" = "pAKT", "FCoCntrl_M2_P70s6k" = "pP70S6K")


g <- Plot_TwoPanel_ValidationPlot_updated(CTDR_Allplot_FCoCntrl,"log2Treatment", "Signal" ,Analyte_labels,"fixed")+
  scale_y_continuous(limits = c(0, 1.6)) + ## In this case adding the space for putting stars was easier since the scale was fixed on y also
  labs(x = "\nGSK3i(\u03bcM)+1 [log2]", #\u03bc is the unicode charachter fro greek mu
       y = "Rel. phosp. (norm.) \n Fold change over untreated",
      # y = " Phosph. rel. to \n untreated \n",
       color = "Cell line" )  # color within labs,lets me give user defined labels to the attribute in legend


gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig4C_F_CTDR_pAkt_pP70s6k_FCoCntrl.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that


########### T Tests for change from respective untreated cell lines(using FCoCntrl) ########

Gsk3i_dose = c(0,1.5,3,6,12,24)

Gsk3_XX_Akt_pvals = list()
for(i in 1:length(Gsk3i_dose) ) {
  
  #Unhash_2_Debug
  #i=1
  Gsk3_XX_Akt <- CTDR_Allplot_FCoCntrl %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "FCoCntrl_M2_Akt") %>% 
    filter(Treatment == Gsk3i_dose[i])
  
  Gsk3_XX_Akt_pvals$Treatment_value[i] = Gsk3i_dose[i]
  Gsk3_XX_Akt_pvals$mean[i] = t.test(Gsk3_XX_Akt$Signal, mu=1)$estimate
  Gsk3_XX_Akt_pvals$p_value[i] = t.test(Gsk3_XX_Akt$Signal, mu=1)$p.value
  
}

Gsk3_XX_Akt_Ttest <- data.frame((sapply(Gsk3_XX_Akt_pvals,c)))
Gsk3_XX_Akt_Ttest <- Gsk3_XX_Akt_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))
WriteXLS::WriteXLS(Gsk3_XX_Akt_Ttest, "./OUTPUT_PAPER/Gsk3_XX_Akt_Ttest.xls")

Gsk3_XO_Akt_pvals = list()
for(i in 1:length(Gsk3i_dose)) {
  
  Gsk3_XO_Akt <- CTDR_Allplot_FCoCntrl %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "FCoCntrl_M2_Akt") %>% 
    filter(Treatment == Gsk3i_dose[i]) 
  
  Gsk3_XO_Akt_pvals$Treatment_value[i] = Gsk3i_dose[i]
  Gsk3_XO_Akt_pvals$mean[i] = t.test(Gsk3_XO_Akt$Signal, mu=1)$estimate
  Gsk3_XO_Akt_pvals$p_value[i] = t.test(Gsk3_XO_Akt$Signal, mu=1)$p.value
  
}

Gsk3_XO_Akt_Ttest <- data.frame((sapply(Gsk3_XO_Akt_pvals,c)))
Gsk3_XO_Akt_Ttest <- Gsk3_XO_Akt_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))
WriteXLS::WriteXLS(Gsk3_XO_Akt_Ttest, "./OUTPUT_PAPER/Gsk3_XO_Akt_Ttest.xls")

Gsk3_XX_P70s6k_pvals = list()
for(i in 1:length(Gsk3i_dose) ) {
  
  #Unhash_2_Debug
  #i=1
  Gsk3_XX_P70s6k <- CTDR_Allplot_FCoCntrl %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XX") %>% 
    filter(Analyte == "FCoCntrl_M2_P70s6k") %>% 
    filter(Treatment == Gsk3i_dose[i])
  
  Gsk3_XX_P70s6k_pvals$Treatment_value[i] = Gsk3i_dose[i]
  Gsk3_XX_P70s6k_pvals$mean[i] = t.test(Gsk3_XX_P70s6k$Signal, mu=1)$estimate
  Gsk3_XX_P70s6k_pvals$p_value[i] = t.test(Gsk3_XX_P70s6k$Signal, mu=1)$p.value
  
}

Gsk3_XX_P70s6k_Ttest <- data.frame((sapply(Gsk3_XX_P70s6k_pvals,c)))
Gsk3_XX_P70s6k_Ttest <- Gsk3_XX_P70s6k_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))
WriteXLS::WriteXLS(Gsk3_XX_P70s6k_Ttest, "./OUTPUT_PAPER/Gsk3_XX_P70s6k_Ttest.xls")

Gsk3_XO_P70s6k_pvals = list()
for(i in 1:length(Gsk3i_dose)) {
  
  Gsk3_XO_P70s6k <- CTDR_Allplot_FCoCntrl %>% 
    ungroup() %>% 
    select(-c(log2Treatment,Treatment_Fctr)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Analyte == "FCoCntrl_M2_P70s6k") %>% 
    filter(Treatment == Gsk3i_dose[i]) 
  
  Gsk3_XO_P70s6k_pvals$Treatment_value[i] = Gsk3i_dose[i]
  Gsk3_XO_P70s6k_pvals$mean[i] = t.test(Gsk3_XO_P70s6k$Signal, mu=1)$estimate
  Gsk3_XO_P70s6k_pvals$p_value[i] = t.test(Gsk3_XO_P70s6k$Signal, mu=1)$p.value
  
}

Gsk3_XO_P70s6k_Ttest <- data.frame((sapply(Gsk3_XO_P70s6k_pvals,c)))
Gsk3_XO_P70s6k_Ttest <- Gsk3_XO_P70s6k_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))
WriteXLS::WriteXLS(Gsk3_XO_P70s6k_Ttest, "./OUTPUT_PAPER/Gsk3_XO_P70s6k_Ttest.xls")

TEMP_FILES <- list.files(pattern = "Rplots.pdf$")
file.remove(TEMP_FILES)

print("Temporary files deleted")
print("All steps executed ")
