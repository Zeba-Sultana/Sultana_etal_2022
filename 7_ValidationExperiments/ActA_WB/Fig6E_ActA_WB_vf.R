#!/usr/bin/env Rscript

library(readxl) 
library(tidyverse)
library(ggplot2)
library(dplyr)
library(egg) #needed for set_panel_size

source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")
source("../ValidationPlot_Functions/Functions_Revision.R")

quant_folder= "../../RAW_DATA/ValidationExperiments/ActA_WB/"
dir.create("./OUTPUT_PAPER") # Create folder to save output data

######### Activin DR ###########

#########  pSmad2_Smad2 : Activin DR ###########
#########  BlotA : R1 + R2  #########  

BlotA_samples = c(
  "E14_XY_R1_0_cntrl",
  "E14_XY_R1_30_ActA",
  "E14_XY_R2_0_cntrl",
  "E14_XY_R2_30_ActA")

BlotA_pSmad2_quant_file = c("BlotA_230622_Activin_E14_pSMAD2.xls")
BlotA_pSmad2_DB <-  read_excel(file.path(quant_folder,BlotA_pSmad2_quant_file))
BlotA_pSmad2_DB$Sample <- BlotA_samples
BlotA_pSmad2 <- BlotA_pSmad2_DB %>% 
  dplyr::select(Signal,Sample)
colnames(BlotA_pSmad2) <- gsub("Signal","phosP",colnames(BlotA_pSmad2))


BlotA_TotSmad2_quant_file = c("BlotA_230622_Activin_E14_totSMAD2.xls")
BlotA_TotSmad2_DB <-  read_excel(file.path(quant_folder,BlotA_TotSmad2_quant_file))
BlotA_TotSmad2_DB$Sample <- BlotA_samples
BlotA_TotSmad2 <- BlotA_TotSmad2_DB %>% 
dplyr::select(Signal,Sample)
colnames(BlotA_TotSmad2) <- gsub("Signal","totalP",colnames(BlotA_TotSmad2))

BlotA_ActA = PrepareWBData_rev(BlotA_pSmad2,BlotA_TotSmad2, GelNum = "BlotA", Analyte = "SMAD2")


#########  Normalize for loading #####
BlotA_ActA_Norm = Norm_TotalP(BlotA_ActA)

######## Fold change Over untreated XX ##############
BlotA_ActA_FC <- FC_Calculation_rev(BlotA_ActA_Norm)

BlotA_ActA_FC_l <- BlotA_ActA_FC %>% 
  pivot_longer(cols = c(FCoXX,FCoCntrl), names_to = "Variable", values_to = "Value")

CellLineColours = c( "XX"="#fb9a99","XO"="#1f78b4", "XY"="black")

BlotA_ActA_FC_l_FCoCntrl = BlotA_ActA_FC_l %>% 
  filter(Variable == "FCoCntrl") %>% 
  mutate(sample = paste(Cell_line,x_status,TreatmentType, sep = "_"))

sample_order = c("E14_XY_cntrl","E14_XY_ActA")
sample_names = c("-", "ActA")

BlotA_ActA_FC_l_FCoCntrl$sample = factor(BlotA_ActA_FC_l_FCoCntrl$sample,levels=sample_order, labels = sample_names )
BlotA_ActA_FC_l_FCoCntrl$x_status = factor(BlotA_ActA_FC_l_FCoCntrl$x_status, levels = c("XX","XO","XY"))



g=ggplot(BlotA_ActA_FC_l_FCoCntrl, aes(x=sample, y = Value )) +
  geom_point(aes(colour=x_status, shape = Replicate)) +
  scale_fill_manual(values = CellLineColours) +
  scale_color_manual(values = CellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 0.02), aes(width=0.5))+
  labs(x = "", 
       y = " pSMAD2 \n Rel. phosp. (norm.) ",
       title = "XY" )+
  MyPaperTheme+
  ylim(0,NA)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  theme(legend.position = "none")

gt=egg::set_panel_size(g,width=unit(1.4,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT_PAPER/Fig6E_Activin_pSmad2_E14XY.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./") 



########## Saving source data ##########

Fig6E_ActA_pSMAD2 = BlotA_ActA_FC_l_FCoCntrl %>% 
  select(c(sample,Cell_line,x_status,Replicate,TreatmentType,Variable,Value))
WriteXLS::WriteXLS(Fig6E_ActA_pSMAD2, "./OUTPUT_PAPER/Fig6E_ActA_pSMAD2.xls")


##### t Test
untreated = BlotA_ActA_FC_l_FCoCntrl %>%
  filter(TreatmentType == "cntrl") 
Vals_untreated= untreated$Value
  
ActAtreated = BlotA_ActA_FC_l_FCoCntrl %>%
  filter(TreatmentType == "ActA") 
Vals_ActAtreated= ActAtreated$Value

t.test(Vals_untreated,Vals_ActAtreated)
