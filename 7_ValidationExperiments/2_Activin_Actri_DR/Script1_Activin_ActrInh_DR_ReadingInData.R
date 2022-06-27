#!/usr/bin/env Rscript

library(readxl) #Part of tidyverse, but needs to be installed separately.(Working on Mac but not my Linux - Check !)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(egg) #needed for set_panel_size

source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

# quant_folder="../RAW_DATA/ValidationExperiments_WB/Validation_data1/My_Band_Quants/Analytes/"
# #image_folder="../RAW_DATA/ValidationExperiments_WB/Validation_data1/My_Pics/Analytes/"
# empiria_folder = "../RAW_DATA/ValidationExperiments_WB/Validation_data1/My_Band_Quants/TPS_Empiria/EmpiriaReports/"

quant_folder= "../../RAW_DATA/ValidationExperiments/ActA_ACTRi_DR/Analyte_Quant/"
empiria_folder = "../../RAW_DATA/ValidationExperiments/ActA_ACTRi_DR/TPS_Quant/"


dir.create("./OUTPUT") # Create folder to save output data

######### Activin DR ###########

Activin_Doses = c("0","1.11","3.33","10","30","90")

#########  pSmad2_Smad2 : Activin DR ###########
#########  Gel11 : R1 + R3XX   #########  

Gel11_samples = c(
  "XO_R3_cntrl_ActDR",
  "XX_R1_cntrl_ActDR",
  "XX_R1_1.11_ActDR",
  "XX_R1_3.33_ActDR",
  "XX_R1_10_ActDR",
  "XX_R1_30_ActDR",
  "XX_R1_90_ActDR",
  "XO_R1_cntrl_ActDR",
  "XO_R1_1.11_ActDR",
  "XO_R1_3.33_ActDR",
  "XO_R1_10_ActDR",
  "XO_R1_30_ActDR",
  "XO_R1_90_ActDR",
  "XX_R3_cntrl_ActDR",
  "XX_R3_1.11_ActDR",
  "XX_R3_3.33_ActDR",
  "XX_R3_10_ActDR",
  "XX_R3_30_ActDR",
  "XX_R3_90_ActDR")
Gel11_Smad2_quant_file = c("0000030_03_JuneGel11_ActivinDR_Smads.xls")
Gel11_Smad2DB <-  read_excel(file.path(quant_folder,Gel11_Smad2_quant_file))

Gel11_Smad2DB_labeled <- ReadInWBData(Gel11_samples, Gel11_Smad2DB, Activin_Doses, "Gel11", "Smad2")
Gel11_Smad2DB_NormOCom <- Norm_o_ComSample(Gel11_Smad2DB_labeled, ComSample_Exp = "ActDR", ComSample_Rep = "R3", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.


#########  Gel12 : R2 + R3XO    #########  

Gel12_samples = c(
  "XX_R3_cntrl_ActDR",
  "XX_R2_cntrl_ActDR",
  "XX_R2_1.11_ActDR",
  "XX_R2_3.33_ActDR",
  "XX_R2_10_ActDR",
  "XX_R2_30_ActDR",
  "XX_R2_90_ActDR",
  "XO_R2_cntrl_ActDR",
  "XO_R2_1.11_ActDR",
  "XO_R2_3.33_ActDR",
  "XO_R2_10_ActDR",
  "XO_R2_30_ActDR",
  "XO_R2_90_ActDR",
  "XO_R3_cntrl_ActDR",
  "XO_R3_1.11_ActDR",
  "XO_R3_3.33_ActDR",
  "XO_R3_10_ActDR",
  "XO_R3_30_ActDR",
  "XO_R3_90_ActDR")
Gel12_Smad2_quant_file = c("0000031_01_JuneGel12_ActivinDR_Smads.xls")
Gel12_Smad2DB <-  read_excel(file.path(quant_folder,Gel12_Smad2_quant_file))

Gel12_Smad2DB_labeled <- ReadInWBData(Gel12_samples, Gel12_Smad2DB, Activin_Doses, "Gel12", "Smad2")
Gel12_Smad2DB_NormOCom <- Norm_o_ComSample(Gel12_Smad2DB_labeled, ComSample_Exp = "ActDR", ComSample_Rep = "R3", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.

#########  Collate the 3 REPLICATES #####

ActDR_Smad2_NormOCom_all <- dplyr::bind_rows(Gel11_Smad2DB_NormOCom, Gel12_Smad2DB_NormOCom) 
ActDR_Smad2_NormOCom <-   ActDR_Smad2_NormOCom_all %>% 
  filter(Well_no!=1) # To filter out the common samples that were added only for normalization

######### Normalization OVER MEAN SIGNAL per replicate #######
ActDR_Smad2_NormOMeanRep <- Norm_o_MeanRep(ActDR_Smad2_NormOCom,Norm_o_Com = 1)

######## Fold change Over XX and Over resp Cntrl ##############
ActDR_Smad2_FC <- FC_Calculation_updated(ActDR_Smad2_NormOMeanRep)

#########  pErk_Erk : Activin DR ###########

#########  Gel11 : R1 + R3XX   #########  
Gel11_Erk_quant_file = c("0000032_02_JuneGel11_ActivinDR_Erk.xls")
Gel11_ErkDB <-  read_excel(file.path(quant_folder,Gel11_Erk_quant_file))

Gel11_ErkDB_labeled <- ReadInWBData(Gel11_samples, Gel11_ErkDB, Activin_Doses, "Gel11", "Erk")
Gel11_ErkDB_NormOCom <- Norm_o_ComSample(Gel11_ErkDB_labeled, ComSample_Exp = "ActDR", ComSample_Rep = "R3", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.

#########  Gel12 : R2 + R3XO   #########  
Gel12_Erk_quant_file = c("0000033_02_JuneGel12_ActivinDR_Erk.xls")
Gel12_ErkDB <-  read_excel(file.path(quant_folder,Gel12_Erk_quant_file))

Gel12_ErkDB_labeled <- ReadInWBData(Gel12_samples, Gel12_ErkDB, Activin_Doses, "Gel12", "Erk")
Gel12_ErkDB_NormOCom <- Norm_o_ComSample(Gel12_ErkDB_labeled, ComSample_Exp = "ActDR", ComSample_Rep = "R3", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.

#########  Collate the 3 REPLICATES #####

ActDR_Erk_NormOCom_all <- dplyr::bind_rows(Gel11_ErkDB_NormOCom, Gel12_ErkDB_NormOCom) 

ActDR_Erk_NormOCom <- ActDR_Erk_NormOCom_all %>% 
  filter(Well_no!=1) # To filter out the common samples that were added only for normalization

######### Normalization OVER MEAN SIGNAL per replicate #######

ActDR_Erk_NormOMeanRep <- Norm_o_MeanRep(ActDR_Erk_NormOCom,Norm_o_Com = 1)

######## Fold change Over XX and Over resp Cntrl ##############

ActDR_Erk_FC <- FC_Calculation_updated(ActDR_Erk_NormOMeanRep)


######### Activin Receptor Inhibitor DR ###########

ActRinh_SB_Doses = c("0","0.2","0.6","1.77","5.33","16")

#########  pSmad2_Smad2 : ActRi DR ###########

#########  Gel13 : R1 + R3XX   #########  

Gel13_samples = c(
  "XO_R3_cntrl_SBDR",
  "XX_R1_cntrl_SBDR",
  "XX_R1_0.2_SBDR",
  "XX_R1_0.6_SBDR",
  "XX_R1_1.77_SBDR",
  "XX_R1_5.33_SBDR",
  "XX_R1_16_SBDR",
  "XO_R1_cntrl_SBDR",
  "XO_R1_0.2_SBDR",
  "XO_R1_0.6_SBDR",
  "XO_R1_1.77_SBDR",
  "XO_R1_5.33_SBDR",
  "XO_R1_16_SBDR",
  "XX_R3_cntrl_SBDR",
  "XX_R3_0.2_SBDR",
  "XX_R3_0.6_SBDR",
  "XX_R3_1.77_SBDR",
  "XX_R3_5.33_SBDR",
  "XX_R3_16_SBDR")
Gel13_Smad2_quant_file = c("0000034_02_JuneGel13_SBDR_Smads.xls")

Gel13_Smad2DB <-  read_excel(file.path(quant_folder,Gel13_Smad2_quant_file))

Gel13_Smad2DB_labeled <- ReadInWBData(Gel13_samples, Gel13_Smad2DB, ActRinh_SB_Doses, "Gel13", "Smad2")
Gel13_Smad2DB_NormOCom <- Norm_o_ComSample(Gel13_Smad2DB_labeled, ComSample_Exp = "SBDR", ComSample_Rep = "R3", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.


#########  Gel14 : R2 + R3XO    #########  

Gel14_samples = c(
  "XX_R3_cntrl_SBDR",
  "XX_R2_cntrl_SBDR",
  "XX_R2_0.2_SBDR",
  "XX_R2_0.6_SBDR",
  "XX_R2_1.77_SBDR",
  "XX_R2_5.33_SBDR",
  "XX_R2_16_SBDR",
  "XO_R2_cntrl_SBDR",
  "XO_R2_0.2_SBDR",
  "XO_R2_0.6_SBDR",
  "XO_R2_1.77_SBDR",
  "XO_R2_5.33_SBDR",
  "XO_R2_16_SBDR",
  "XO_R3_cntrl_SBDR",
  "XO_R3_0.2_SBDR",
  "XO_R3_0.6_SBDR",
  "XO_R3_1.77_SBDR",
  "XO_R3_5.33_SBDR",
  "XO_R3_16_SBDR")

Gel14_Smad2_quant_file = c("0000035_02_JuneGel14_SBDR_Smads.xls")
Gel14_Smad2DB <-  read_excel(file.path(quant_folder,Gel14_Smad2_quant_file))

Gel14_Smad2DB_labeled <- ReadInWBData(Gel14_samples, Gel14_Smad2DB, ActRinh_SB_Doses, "Gel14", "Smad2")
Gel14_Smad2DB_NormOCom <- Norm_o_ComSample(Gel14_Smad2DB_labeled, ComSample_Exp = "SBDR", ComSample_Rep = "R3", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.

#########  Collate the 3 REPLICATES #####

SBDR_Smad2_NormOCom_all <- dplyr::bind_rows(Gel13_Smad2DB_NormOCom, Gel14_Smad2DB_NormOCom)

SBDR_Smad2_NormOCom <- SBDR_Smad2_NormOCom_all %>% 
  filter(Well_no!=1) # To filter out the common samples that were added only for normalization

######### Normalization OVER MEAN SIGNAL per replicate #######

SBDR_Smad2_NormOMeanRep <- Norm_o_MeanRep(SBDR_Smad2_NormOCom,Norm_o_Com = 1)

######## Fold change Over XX and Over resp Cntrl ##############

SBDR_Smad2_FC <- FC_Calculation_updated(SBDR_Smad2_NormOMeanRep)


#########  pErk_Erk : ActRi DR ###########
#########  Gel13 : R1 + R3XX   #########  

Gel13_Erk_quant_file = c("0000036_02_JuneGel13_SBDR_pErk.xls")
Gel13_ErkDB <-  read_excel(file.path(quant_folder,Gel13_Erk_quant_file))

Gel13_ErkDB_labeled <- ReadInWBData(Gel13_samples, Gel13_ErkDB, ActRinh_SB_Doses, "Gel13", "Erk")
Gel13_ErkDB_NormOCom <- Norm_o_ComSample(Gel13_ErkDB_labeled, ComSample_Exp = "SBDR", ComSample_Rep = "R3", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.

#########  Gel14 : R2 + R3XO    #########  

Gel14_Erk_quant_file = c("0000037_01_JuneGel14_SBDR_pErk.xls")
Gel14_ErkDB <-  read_excel(file.path(quant_folder,Gel14_Erk_quant_file))

Gel14_ErkDB_labeled <- ReadInWBData(Gel14_samples, Gel14_ErkDB, ActRinh_SB_Doses, "Gel14", "Erk")
Gel14_ErkDB_NormOCom <- Norm_o_ComSample(Gel14_ErkDB_labeled, ComSample_Exp = "SBDR", ComSample_Rep = "R3", ComSample_T = 0 ) # I have 2 common samples, untreated control of both cell lines of R3. Hence I dont specify cell line.

#########  Collate the 3 REPLICATES #####
SBDR_Erk_NormOCom_all <- dplyr::bind_rows(Gel13_ErkDB_NormOCom, Gel14_ErkDB_NormOCom)

SBDR_Erk_NormOCom <- SBDR_Erk_NormOCom_all %>% 
  filter(Well_no!=1) # To filter out the common samples that were added only for normalization

######### Normalization OVER MEAN SIGNAL per replicate #######

SBDR_Erk_NormOMeanRep <- Norm_o_MeanRep(SBDR_Erk_NormOCom,Norm_o_Com = 1)

######## Fold change Over XX and Over resp Cntrl ##############

SBDR_Erk_FC <- FC_Calculation_updated(SBDR_Erk_NormOMeanRep)

######## Writing the csv Files for next step ##############

write.csv(ActDR_Smad2_FC, file = "OUTPUT/ActDR_Smad2_FC.csv")
#write.csv(ActDR_Erk_FC, file = "OUTPUT/ActDR_Erk_FC.csv")
write.csv(SBDR_Smad2_FC, file = "OUTPUT/SBDR_Smad2_FC.csv")
#write.csv(SBDR_Erk_FC, file = "OUTPUT/SBDR_Erk_FC.csv")

print("All steps executed : Script1 : Data for ActDR and SBDR read-in")
