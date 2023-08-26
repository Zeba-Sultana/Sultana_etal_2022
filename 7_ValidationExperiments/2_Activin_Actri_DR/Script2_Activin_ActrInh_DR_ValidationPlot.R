#!/usr/bin/env Rscript

library(readxl)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(egg) #needed for set_panel_size

dir.create("OUTPUT_PAPER") 


source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")


######## ActDR : Smad2 #########

ActDR_Smad2_FC <- read.csv("./OUTPUT/ActDR_Smad2_FC.csv")
ActDR_Smad2_FC <- ActDR_Smad2_FC %>% 
  select(c(Cell_line,Replicate,Treatment,Exp,phosP,totalP,Gel,log2Treatment,phosP_Norm_oMeanRep,totalP_Norm_oMeanRep,Norm_phosP_by_total_M2, FCoXX_M2,FCoCntrl_M2))

colnames(ActDR_Smad2_FC) <- gsub("^phosP$","pSmad2",colnames(ActDR_Smad2_FC))
colnames(ActDR_Smad2_FC) <- gsub("^totalP$","Smad2",colnames(ActDR_Smad2_FC))
colnames(ActDR_Smad2_FC) <- gsub("^phosP_Norm_oMeanRep$","Norm_pSmad2",colnames(ActDR_Smad2_FC))
colnames(ActDR_Smad2_FC) <- gsub("^totalP_Norm_oMeanRep$","Norm_Smad2",colnames(ActDR_Smad2_FC))
colnames(ActDR_Smad2_FC) <- gsub("^Norm_phosP_by_total_M2$","pSmad2_by_Smad2",colnames(ActDR_Smad2_FC))
colnames(ActDR_Smad2_FC) <- gsub("^FCoXX_M2$","pSmad2_by_Smad2_FCoXX",colnames(ActDR_Smad2_FC))
colnames(ActDR_Smad2_FC) <- gsub("^FCoCntrl_M2$","pSmad2_by_Smad2_FCoCntrl",colnames(ActDR_Smad2_FC))

Activin_DR_Smad2 <- ActDR_Smad2_FC


######## SBDR : Smad2 #########

SBDR_Smad2_FC <- read.csv("./OUTPUT/SBDR_Smad2_FC.csv")

SBDR_Smad2_FC <- SBDR_Smad2_FC %>% 
  select(c(Cell_line,Replicate,Treatment,Exp,phosP,totalP,Gel,log2Treatment,phosP_Norm_oMeanRep,totalP_Norm_oMeanRep,Norm_phosP_by_total_M2, FCoXX_M2,FCoCntrl_M2))

colnames(SBDR_Smad2_FC) <- gsub("^phosP$","pSmad2",colnames(SBDR_Smad2_FC))
colnames(SBDR_Smad2_FC) <- gsub("^totalP$","Smad2",colnames(SBDR_Smad2_FC))
colnames(SBDR_Smad2_FC) <- gsub("^phosP_Norm_oMeanRep$","Norm_pSmad2",colnames(SBDR_Smad2_FC))
colnames(SBDR_Smad2_FC) <- gsub("^totalP_Norm_oMeanRep$","Norm_Smad2",colnames(SBDR_Smad2_FC))
colnames(SBDR_Smad2_FC) <- gsub("^Norm_phosP_by_total_M2$","pSmad2_by_Smad2",colnames(SBDR_Smad2_FC))
colnames(SBDR_Smad2_FC) <- gsub("^FCoXX_M2$","pSmad2_by_Smad2_FCoXX",colnames(SBDR_Smad2_FC))
colnames(SBDR_Smad2_FC) <- gsub("^FCoCntrl_M2$","pSmad2_by_Smad2_FCoCntrl",colnames(SBDR_Smad2_FC))

SB_DR_Smad2 <- SBDR_Smad2_FC


############### Paper Fig : FCoXX : Linear y-Scale  ####################

g <- Plot_ReplicatePoints_MeanLine_PaperFig_LOG2y(Activin_DR_Smad2,log2Treatment,pSmad2_by_Smad2_FCoXX)+
  ylim(0, 5)+
  labs(x = paste0("\n Activin(ng/ml)+1 [log2]"),
       y = "Rel.phos.(norm)\n") 
#g <- g + guides(fill = guide_legend(order=1),shape = guide_legend(order=2)) 
gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig6C_ActA_DR_pSmad2.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER")


g <- Plot_ReplicatePoints_MeanLine_PaperFig_Minus_LOG2y(SB_DR_Smad2,log2Treatment,pSmad2_by_Smad2_FCoXX)+
  ylim(0, 5)+
  labs(x = paste0("\n (-1)*ActRi(\u03bcM)+1 [log2]"),
       y = "Rel.phos.(norm)\n") 
#g <- g + guides(fill = guide_legend(order=1),shape = guide_legend(order=2)) 
gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig6C_ActRInhibitor_DR_pSmad2.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER")

########### T-tests ##########
########### T Tests for difference between XX and XO in case of FCoXX (Linear y-axis)########

Activin_dose = c(0,1.11,3.33,10,30,90)
ActDR_XX_XO_Smad2_pvals_FCoXX = list()
for(i in 1:length(Activin_dose)) {
  
  ActDR_XX_Smad2_values <- Activin_DR_Smad2 %>% 
    ungroup() %>% 
    select(c(Cell_line,Treatment,pSmad2_by_Smad2_FCoXX)) %>%
    filter(Cell_line == "XX") %>% 
    filter(Treatment == Activin_dose[i])
  
  ActDR_XO_Smad2_values <- Activin_DR_Smad2 %>% 
    ungroup() %>% 
    select(c(Cell_line,Treatment,pSmad2_by_Smad2_FCoXX)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Treatment == Activin_dose[i])
  
  ActDR_XX_XO_Smad2_pvals_FCoXX$Treatment_value[i] = Activin_dose[i]
  ActDR_XX_XO_Smad2_pvals_FCoXX$mean_XX[i] = t.test(ActDR_XX_Smad2_values$pSmad2_by_Smad2_FCoXX,ActDR_XO_Smad2_values$pSmad2_by_Smad2_FCoXX, paired = FALSE)$estimate[1]
  ActDR_XX_XO_Smad2_pvals_FCoXX$mean_XO[i] = t.test(ActDR_XX_Smad2_values$pSmad2_by_Smad2_FCoXX,ActDR_XO_Smad2_values$pSmad2_by_Smad2_FCoXX, paired = FALSE)$estimate[2]
  ActDR_XX_XO_Smad2_pvals_FCoXX$p_value[i] = t.test(ActDR_XX_Smad2_values$pSmad2_by_Smad2_FCoXX,ActDR_XO_Smad2_values$pSmad2_by_Smad2_FCoXX, paired = FALSE)$p.value
  
}

ActDR_XX_XO_Smad2_Ttest <- data.frame((sapply(ActDR_XX_XO_Smad2_pvals_FCoXX,c)))
ActDR_XX_XO_Smad2_Ttest <- ActDR_XX_XO_Smad2_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  mutate(Treatment=as.numeric(Treatment_value)) %>% 
  select(-c(Treatment_value))

ActDR_ReplicateValues  = Activin_DR_Smad2 %>% 
  select(Cell_line,Replicate,Exp,Treatment,log2Treatment,pSmad2_by_Smad2_FCoXX) %>% 
  pivot_wider(names_from = Replicate,values_from=pSmad2_by_Smad2_FCoXX) %>% 
  mutate(Mean=rowMeans(select(.,starts_with("R")), na.rm = TRUE))

Fig6C_ActA_DR_data = left_join(ActDR_ReplicateValues,ActDR_XX_XO_Smad2_Ttest, by="Treatment") 
Fig6C_ActA_DR_data = Fig6C_ActA_DR_data %>% 
  mutate(mean_XX = ifelse(Cell_line == "XO", NA, mean_XX)) %>% 
  mutate(mean_XO = ifelse(Cell_line == "XO", NA, mean_XO)) %>% 
  mutate(p_value = ifelse(Cell_line == "XO", NA, p_value)) %>% 
  mutate(sig = ifelse(Cell_line == "XO", NA, sig))
colnames(Fig6C_ActA_DR_data) = gsub("p_value","p_value(XXvsXO)",colnames(Fig6C_ActA_DR_data))

WriteXLS::WriteXLS(Fig6C_ActA_DR_data, "./OUTPUT_PAPER/Fig6C_ActA_DR_pSmad2_data.xls")


###### Activin Receptor Inhibitor(SB) dose response ###########

SB_dose = c(0,0.2,0.6,1.77,5.33,16)
SBDR_XX_XO_Smad2_pvals_FCoXX = list()
for(i in 1:length(SB_dose)) {
  
  SBDR_XX_Smad2_values <- SB_DR_Smad2 %>% 
    ungroup() %>% 
    select(c(Cell_line,Treatment,pSmad2_by_Smad2_FCoXX)) %>%
    filter(Cell_line == "XX") %>% 
    filter(Treatment == SB_dose[i])
  
  SBDR_XO_Smad2_values <- SB_DR_Smad2 %>% 
    ungroup() %>% 
    select(c(Cell_line,Treatment,pSmad2_by_Smad2_FCoXX)) %>% 
    filter(Cell_line == "XO") %>% 
    filter(Treatment == SB_dose[i])
  
  SBDR_XX_XO_Smad2_pvals_FCoXX$Treatment_value[i] = SB_dose[i]
  SBDR_XX_XO_Smad2_pvals_FCoXX$mean_XX[i] = t.test(SBDR_XX_Smad2_values$pSmad2_by_Smad2_FCoXX,SBDR_XO_Smad2_values$pSmad2_by_Smad2_FCoXX, paired = FALSE)$estimate[1]
  SBDR_XX_XO_Smad2_pvals_FCoXX$mean_XO[i] = t.test(SBDR_XX_Smad2_values$pSmad2_by_Smad2_FCoXX,SBDR_XO_Smad2_values$pSmad2_by_Smad2_FCoXX, paired = FALSE)$estimate[2]
  SBDR_XX_XO_Smad2_pvals_FCoXX$p_value[i] = t.test(SBDR_XX_Smad2_values$pSmad2_by_Smad2_FCoXX,SBDR_XO_Smad2_values$pSmad2_by_Smad2_FCoXX, paired = FALSE)$p.value
  
}

SBDR_XX_XO_Smad2_Ttest <- data.frame((sapply(SBDR_XX_XO_Smad2_pvals_FCoXX,c)))
SBDR_XX_XO_Smad2_Ttest <- SBDR_XX_XO_Smad2_Ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  mutate(Treatment=as.numeric(Treatment_value)) %>% 
  select(-c(Treatment_value))

SBDR_ReplicateValues  = SB_DR_Smad2 %>% 
  select(Cell_line,Replicate,Exp,Treatment,log2Treatment,pSmad2_by_Smad2_FCoXX) %>% 
  pivot_wider(names_from = Replicate,values_from=pSmad2_by_Smad2_FCoXX) %>% 
  mutate(Mean=rowMeans(select(.,starts_with("R")), na.rm = TRUE))

Fig6C_ActRInhibitor_DR_data = left_join(SBDR_ReplicateValues,SBDR_XX_XO_Smad2_Ttest, by="Treatment") 
Fig6C_ActRInhibitor_DR_data = Fig6C_ActRInhibitor_DR_data %>% 
  mutate(mean_XX = ifelse(Cell_line == "XO", NA, mean_XX)) %>% 
  mutate(mean_XO = ifelse(Cell_line == "XO", NA, mean_XO)) %>% 
  mutate(p_value = ifelse(Cell_line == "XO", NA, p_value)) %>% 
  mutate(sig = ifelse(Cell_line == "XO", NA, sig))
colnames(Fig6C_ActRInhibitor_DR_data) = gsub("p_value","p_value(XXvsXO)",colnames(Fig6C_ActRInhibitor_DR_data))

WriteXLS::WriteXLS(Fig6C_ActRInhibitor_DR_data, "./OUTPUT_PAPER/Fig6C_ActRInhibitor_DR_pSmad2_data.xls")



print("All steps executed : Script 2 : plots saved in OUTPUT_PAPER")


