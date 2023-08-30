#!/usr/bin/env Rscript

library(readxl) 
library(tidyverse)
library(ggplot2)
library(dplyr)
library(egg) #needed for set_panel_size

source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")
source("../ValidationPlot_Functions/Functions_Revision.R")

quant_folder= "../../RAW_DATA/ValidationExperiments/U0126_WB/"
dir.create("./OUTPUT_PAPER") # Create folder to save output data

######### U0126 : pMek ###########
#########  BlotC: R1 #########  
BlotCR1_samples = c(
  "1.8_XX_R1_0_dmso",
  "1.8_XX_R1_5_U0",
  "1.8_XO_R1_0_dmso",
  "1.8_XO_R1_5_U0",
  "E14_XY_R1_0_dmso",
  "E14_XY_R1_5_U0")

BlotCR1_pMek_quant_file = c("Blot_B_CR1_CR2_pMek_230622_1173_700.xls")
BlotCR1_pMek_DB <-  read_excel(file.path(quant_folder,BlotCR1_pMek_quant_file))
BlotCR1_pMek_DB <-  BlotCR1_pMek_DB[c(13:16,21,22),]

BlotCR1_pMek_DB$Sample <- BlotCR1_samples
BlotCR1_pMek <- BlotCR1_pMek_DB %>% 
  dplyr::select(Signal,Sample)
colnames(BlotCR1_pMek) <- gsub("Signal","phosP",colnames(BlotCR1_pMek))


BlotCR1_TPS_quant_file = c("Blot_CR1_CR2_TPS.xls")
BlotCR1_TPS_DB <-  read_excel(file.path(quant_folder,BlotCR1_TPS_quant_file))
BlotCR1_TPS_DB <- BlotCR1_TPS_DB %>% 
  filter(blott == "C_R1") %>% 
  filter(CellLine == "E14"|CellLine == "1.8") %>% 
  select(c(Replicate,CellLine,sex,treatment,sumBYmax,sum))

#BlotCR1_TPS_DB$CellLine = gsub("PgK","Pgk",BlotCR1_TPS_DB$CellLine)
BlotCR1_TPS_DB$treatment = gsub("DSMSO I","0_dmso",BlotCR1_TPS_DB$treatment)
BlotCR1_TPS_DB$treatment = gsub("UO126","5_U0",BlotCR1_TPS_DB$treatment)
BlotCR1_TPS_DB$Sample = paste(BlotCR1_TPS_DB$CellLine,BlotCR1_TPS_DB$sex,BlotCR1_TPS_DB$Replicate,BlotCR1_TPS_DB$treatment, sep = "_")

#############

BlotCR1_TPS <- BlotCR1_TPS_DB %>% 
  dplyr::select(Sample,sum)
colnames(BlotCR1_TPS) = c("Sample","totalP")

BlotCR1_pMek_data = PrepareWBData_rev(BlotCR1_pMek,BlotCR1_TPS, GelNum = "BlotCR1", Analyte = "pMek")


#########  Normalize for loading #####
BlotCR1_pMek_Norm = Norm_TotalP(BlotCR1_pMek_data)

######## Fold change Over untreated XX ##############
#BlotCR1_pMek_FC <- FC_Calculation_rev_eachXX(BlotCR1_pMek_Norm)


#################################################
#########  BlotC: R2 #########  
BlotCR2_samples = c(
  "1.8_XX_R2_0_dmso",
  "1.8_XX_R2_5_U0",
  "1.8_XO_R2_0_dmso",
  "1.8_XO_R2_5_U0",
  "E14_XY_R2_0_dmso",
  "E14_XY_R2_5_U0")

BlotCR2_pMek_quant_file = c("Blot_B_CR1_CR2_pMek_230622_1173_700.xls")
BlotCR2_pMek_DB <-  read_excel(file.path(quant_folder,BlotCR2_pMek_quant_file))
BlotCR2_pMek_DB <-  BlotCR2_pMek_DB[c(23:26,31,32),]

BlotCR2_pMek_DB$Sample <- BlotCR2_samples
BlotCR2_pMek <- BlotCR2_pMek_DB %>% 
  dplyr::select(Signal,Sample)
colnames(BlotCR2_pMek) <- gsub("Signal","phosP",colnames(BlotCR2_pMek))


BlotCR2_TPS_quant_file = c("Blot_CR1_CR2_TPS.xls")
BlotCR2_TPS_DB <-  read_excel(file.path(quant_folder,BlotCR2_TPS_quant_file))
BlotCR2_TPS_DB <- BlotCR2_TPS_DB %>% 
  filter(blott == "C_R2") %>% 
  filter(CellLine == "E14"|CellLine == "1.8") %>% 
  select(c(Replicate,CellLine,sex,treatment,sumBYmax,sum))

#BlotCR2_TPS_DB$CellLine = gsub("PgK","Pgk",BlotCR2_TPS_DB$CellLine)
BlotCR2_TPS_DB$treatment = gsub("DSMSO I","0_dmso",BlotCR2_TPS_DB$treatment)
BlotCR2_TPS_DB$treatment = gsub("UO126","5_U0",BlotCR2_TPS_DB$treatment)
BlotCR2_TPS_DB$Sample = paste(BlotCR2_TPS_DB$CellLine,BlotCR2_TPS_DB$sex,BlotCR2_TPS_DB$Replicate,BlotCR2_TPS_DB$treatment, sep = "_")

#############

BlotCR2_TPS <- BlotCR2_TPS_DB %>% 
  dplyr::select(Sample,sum)
colnames(BlotCR2_TPS) = c("Sample","totalP")

BlotCR2_pMek_data = PrepareWBData_rev(BlotCR2_pMek,BlotCR2_TPS, GelNum = "BlotCR2", Analyte = "pMek")


#########  Normalize for loading #####
BlotCR2_pMek_Norm = Norm_TotalP(BlotCR2_pMek_data)


######## Fold change Over untreated XX ##############
#BlotCR2_pMek_FC <- FC_Calculation_rev_eachXX(BlotCR2_pMek_Norm)

#################################################
################## Binding CR1 and CR2 #########
BlotC_pMek_Norm <- bind_rows(BlotCR1_pMek_Norm,BlotCR2_pMek_Norm)
BlotC_pMek_FC = FC_Calculation_rev_eachXX(BlotC_pMek_Norm)

BlotC_pMek_FC = BlotC_pMek_FC %>% 
  ungroup() %>% 
  mutate(Rel_pMek = phosP_by_total/max(phosP_by_total))

#### Plotting ###

BlotC_pMek_FC_l <- BlotC_pMek_FC %>% 
  pivot_longer(cols = c(Rel_pMek,FCoCntrl), names_to = "Variable", values_to = "Value")


BlotC_pMek_FC_l = BlotC_pMek_FC_l %>% 
  #filter(Variable == "FCoCntrl") %>% 
  mutate(sample = paste(Cell_line,x_status,TreatmentType, sep = "_"))


Fig7C_U0126_pMEK = BlotC_pMek_FC_l %>% 
select(c(sample,Cell_line,x_status,Replicate,TreatmentType,Variable,Value))
  WriteXLS::WriteXLS(Fig7C_U0126_pMEK, "./OUTPUT_PAPER/Fig7C_U0126_pMEK.xls")


sample_order = c("1.8_XX_dmso", "1.8_XX_U0","1.8_XO_dmso","1.8_XO_U0","E14_XY_dmso","E14_XY_U0")
sample_names = c("1.8_XX_dmso", "1.8_XX_U0","1.8_XO_dmso","1.8_XO_U0","E14_XY_dmso","E14_XY_U0")
BlotC_pMek_FC_l$sample = factor(BlotC_pMek_FC_l$sample,levels=sample_order, labels = sample_names )
BlotC_pMek_FC_l$Variable = factor(BlotC_pMek_FC_l$Variable, levels = c("Rel_pMek", "FCoCntrl"), labels = c("Rel. phos.", "FC over untreated"))

CellLineColours = c("XX"="#FF0000","XO"="#0000CD" , "XY"="black")

#############
#To add space at the top of the plot for putting astericks, add dummy data(which is the max value+additional space)

Rel_pMek_max <- BlotC_pMek_FC_l %>%
  filter(grepl("Rel. phos.", Variable))
Rel_pMek_max <-max(Rel_pMek_max$Value) + 0.25 # the number is based on how much space i need

FC_pMek_max <- BlotC_pMek_FC_l %>%
  filter(grepl("FC over untreated", Variable))
FC_pMek_max <-max(FC_pMek_max$Value) + 0.5 # the number is based on how much space i need

dummy_Rel_pMek <- data.frame(sample = "1.8_XX_dmso", Value = Rel_pMek_max,
                             Variable = "Rel. phos.",  stringsAsFactors=FALSE)
dummy_FC_pMek <- data.frame(sample = "1.8_XX_dmso", Value = FC_pMek_max,
                            Variable = "FC over untreated",  stringsAsFactors=FALSE) 
dummy_data <- rbind(dummy_Rel_pMek,dummy_FC_pMek)
dummy_data$sample = factor(dummy_data$sample,levels=sample_order, labels = sample_names )
dummy_data$Variable = factor(dummy_data$Variable, levels =c("Rel. phos.", "FC over untreated")) # This is to get the correct sequence in the Facet_wrap


g=ggplot(BlotC_pMek_FC_l, aes(x=sample, y = Value )) +
  geom_blank(data=dummy_data) + 
  geom_point(aes(colour=x_status, shape=Replicate)) +
  #facet_wrap(.~Variable, scales = "free_y", ncol = 1)+
  facet_grid(Variable~., scales = "free_y", switch = 'y')+
  scale_fill_manual(values = CellLineColours) +
  scale_color_manual(values = CellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 0.02), aes(width=0.5))+
  labs(x = "",
       y = "",
       title = "pMek" )+
  ylim(0,NA)+
  MyPaperTheme+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5, size=7),
        strip.placement = "outside",
        strip.text.y = element_text(size = 7, colour = "black", angle = 90,face = "plain" ))
gt=egg::set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT_PAPER/Fig_7C_U0126_pMEK.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./") 


######### U0126 DR : pcRaf ###########
#########  BlotC: R1 #########  
BlotCR1_samples = c(
  "1.8_XX_R1_0_dmso",
  "1.8_XX_R1_5_U0",
  "1.8_XO_R1_0_dmso",
  "1.8_XO_R1_5_U0",
  "E14_XY_R1_0_dmso",
  "E14_XY_R1_5_U0")

BlotCR1_pcRaf_quant_file = c("Blot_B_CR1_CR2_pcRaf_230621_1171_700.xls")
BlotCR1_pcRaf_DB <-  read_excel(file.path(quant_folder,BlotCR1_pcRaf_quant_file))
BlotCR1_pcRaf_DB <-  BlotCR1_pcRaf_DB[c(13:16,21,22),]

BlotCR1_pcRaf_DB$Sample <- BlotCR1_samples
BlotCR1_pcRaf <- BlotCR1_pcRaf_DB %>% 
  dplyr::select(Signal,Sample)
colnames(BlotCR1_pcRaf) <- gsub("Signal","phosP",colnames(BlotCR1_pcRaf))


BlotCR1_TPS_quant_file = c("Blot_CR1_CR2_TPS.xls")
BlotCR1_TPS_DB <-  read_excel(file.path(quant_folder,BlotCR1_TPS_quant_file))
BlotCR1_TPS_DB <- BlotCR1_TPS_DB %>% 
  filter(blott == "C_R1") %>% 
  filter(CellLine == "E14"|CellLine == "1.8") %>% 
  select(c(Replicate,CellLine,sex,treatment,sumBYmax,sum))

#BlotCR1_TPS_DB$CellLine = gsub("PgK","Pgk",BlotCR1_TPS_DB$CellLine)
BlotCR1_TPS_DB$treatment = gsub("DSMSO I","0_dmso",BlotCR1_TPS_DB$treatment)
BlotCR1_TPS_DB$treatment = gsub("UO126","5_U0",BlotCR1_TPS_DB$treatment)
BlotCR1_TPS_DB$Sample = paste(BlotCR1_TPS_DB$CellLine,BlotCR1_TPS_DB$sex,BlotCR1_TPS_DB$Replicate,BlotCR1_TPS_DB$treatment, sep = "_")

#############

BlotCR1_TPS <- BlotCR1_TPS_DB %>% 
  dplyr::select(Sample,sum)
colnames(BlotCR1_TPS) = c("Sample","totalP")

BlotCR1_pcRaf_data = PrepareWBData_rev(BlotCR1_pcRaf,BlotCR1_TPS, GelNum = "BlotCR1", Analyte = "pcRaf")


#########  Normalize for loading #####
BlotCR1_pcRaf_Norm = Norm_TotalP(BlotCR1_pcRaf_data)

######## Fold change Over untreated XX ##############
#BlotCR1_pcRaf_FC <- FC_Calculation_rev_eachXX(BlotCR1_pcRaf_Norm)


#################################################
#########  BlotC: R2 #########  
BlotCR2_samples = c(
  "1.8_XX_R2_0_dmso",
  "1.8_XX_R2_5_U0",
  "1.8_XO_R2_0_dmso",
  "1.8_XO_R2_5_U0",
  "E14_XY_R2_0_dmso",
  "E14_XY_R2_5_U0")

BlotCR2_pcRaf_quant_file = c("Blot_B_CR1_CR2_pcRaf_230621_1171_700.xls")
BlotCR2_pcRaf_DB <-  read_excel(file.path(quant_folder,BlotCR2_pcRaf_quant_file))
BlotCR2_pcRaf_DB <-  BlotCR2_pcRaf_DB[c(23:26,31,32),]

BlotCR2_pcRaf_DB$Sample <- BlotCR2_samples
BlotCR2_pcRaf <- BlotCR2_pcRaf_DB %>% 
  dplyr::select(Signal,Sample)
colnames(BlotCR2_pcRaf) <- gsub("Signal","phosP",colnames(BlotCR2_pcRaf))


BlotCR2_TPS_quant_file = c("Blot_CR1_CR2_TPS.xls")
BlotCR2_TPS_DB <-  read_excel(file.path(quant_folder,BlotCR2_TPS_quant_file))
BlotCR2_TPS_DB <- BlotCR2_TPS_DB %>% 
  filter(blott == "C_R2") %>% 
  filter(CellLine == "E14"|CellLine == "1.8") %>% 
  select(c(Replicate,CellLine,sex,treatment,sumBYmax,sum))

#BlotCR2_TPS_DB$CellLine = gsub("PgK","Pgk",BlotCR2_TPS_DB$CellLine)
BlotCR2_TPS_DB$treatment = gsub("DSMSO I","0_dmso",BlotCR2_TPS_DB$treatment)
BlotCR2_TPS_DB$treatment = gsub("UO126","5_U0",BlotCR2_TPS_DB$treatment)
BlotCR2_TPS_DB$Sample = paste(BlotCR2_TPS_DB$CellLine,BlotCR2_TPS_DB$sex,BlotCR2_TPS_DB$Replicate,BlotCR2_TPS_DB$treatment, sep = "_")

#############

BlotCR2_TPS <- BlotCR2_TPS_DB %>% 
  dplyr::select(Sample,sum)
colnames(BlotCR2_TPS) = c("Sample","totalP")

BlotCR2_pcRaf_data = PrepareWBData_rev(BlotCR2_pcRaf,BlotCR2_TPS, GelNum = "BlotCR2", Analyte = "pcRaf")


#########  Normalize for loading #####
BlotCR2_pcRaf_Norm = Norm_TotalP(BlotCR2_pcRaf_data)


######## Fold change Over untreated XX ##############
#BlotCR2_pcRaf_FC <- FC_Calculation_rev_eachXX(BlotCR2_pcRaf_Norm)

#################################################
################## Binding CR1 and CR2 #########
BlotC_pcRaf_Norm <- bind_rows(BlotCR1_pcRaf_Norm,BlotCR2_pcRaf_Norm)
BlotC_pcRaf_FC = FC_Calculation_rev_eachXX(BlotC_pcRaf_Norm)

BlotC_pcRaf_FC = BlotC_pcRaf_FC %>% 
  ungroup() %>% 
  mutate(Rel_pcRaf = phosP_by_total/max(phosP_by_total))

#### Plotting ###

BlotC_pcRaf_FC_l <- BlotC_pcRaf_FC %>% 
  pivot_longer(cols = c(Rel_pcRaf,FCoCntrl), names_to = "Variable", values_to = "Value")


BlotC_pcRaf_FC_l = BlotC_pcRaf_FC_l %>% 
  #filter(Variable == "FCoCntrl") %>% 
  mutate(sample = paste(Cell_line,x_status,TreatmentType, sep = "_"))

Fig7C_U0126_pRAF1 = BlotC_pcRaf_FC_l %>% 
  select(c(sample,Cell_line,x_status,Replicate,TreatmentType,Variable,Value))
WriteXLS::WriteXLS(Fig7C_U0126_pRAF1, "./OUTPUT_PAPER/Fig7D_U0126_pRAF1.xls")

sample_order = c("1.8_XX_dmso", "1.8_XX_U0","1.8_XO_dmso","1.8_XO_U0","E14_XY_dmso","E14_XY_U0")
sample_names = c("1.8_XX_dmso", "1.8_XX_U0","1.8_XO_dmso","1.8_XO_U0","E14_XY_dmso","E14_XY_U0")
BlotC_pcRaf_FC_l$sample = factor(BlotC_pcRaf_FC_l$sample,levels=sample_order, labels = sample_names )
BlotC_pcRaf_FC_l$Variable = factor(BlotC_pcRaf_FC_l$Variable, levels = c("Rel_pcRaf", "FCoCntrl"), labels = c("Rel. phos.", "FC over untreated"))


CellLineColours = c("XX"="#FF0000","XO"="#0000CD" , "XY"="black")

#############
#To add space at the top of the plot for putting astericks, add dummy data(which is the max value+additional space)


Rel_pcRaf_max <- BlotC_pcRaf_FC_l %>%
  filter(grepl("Rel. phos.", Variable))
Rel_pcRaf_max <-max(Rel_pcRaf_max$Value) + 0.5 # the number is based on how much space i need

FC_pcRaf_max <- BlotC_pcRaf_FC_l %>%
  filter(grepl("FC over untreated", Variable))
FC_pcRaf_max <-max(FC_pcRaf_max$Value) + 0.5 # the number is based on how much space i need

dummy_Rel_pcRaf <- data.frame(sample = "1.8_XX_dmso", Value = Rel_pcRaf_max,
                              Variable = "Rel. phos.",  stringsAsFactors=FALSE)
dummy_FC_pcRaf <- data.frame(sample = "1.8_XX_dmso", Value = FC_pcRaf_max,
                             Variable = "FC over untreated",  stringsAsFactors=FALSE) 
dummy_data <- rbind(dummy_Rel_pcRaf,dummy_FC_pcRaf)
dummy_data$sample = factor(dummy_data$sample,levels=sample_order, labels = sample_names )
dummy_data$Variable = factor(dummy_data$Variable, levels =c("Rel. phos.", "FC over untreated")) # This is to get the correct sequence in the Facet_wrap


g=ggplot(BlotC_pcRaf_FC_l, aes(x=sample, y = Value )) +
 geom_blank(data=dummy_data,aes(x=factor(sample, levels = sample_order), y = Value )) + 
  geom_point(aes(colour=x_status, shape=Replicate)) +
  #facet_wrap(.~Variable, scales = "free_y", ncol = 1)+
  facet_grid(Variable~., scales = "free_y", switch = 'y')+
  scale_fill_manual(values = CellLineColours) +
  scale_color_manual(values = CellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 0.02), aes(width=0.5))+
  labs(x = "",
       y = "",
       title = "pcRaf" )+
  ylim(0,NA)+
  MyPaperTheme+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5, size=7),
        strip.placement = "outside",
        strip.text.y = element_text(size = 7, colour = "black", angle = 90,face = "plain" ))
# theme(axis.text.x = element_text(angle = 90, hjust = 0.5))

gt=egg::set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT_PAPER/Fig_7D_U0126_pRAF1.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./") 


#### Statistical Tests ####
### Fig. 7C
### pMEK
data = read_excel('./OUTPUT_PAPER/Fig7C_U0126_pMEK.xls') %>%
  mutate (nr.x = ifelse(x_status=='XX',2,1)) %>% 
  pivot_wider(names_from = Variable, values_from = Value)

data2 = data %>% select(-FCoCntrl,-sample) %>% 
  pivot_wider(names_from=TreatmentType, values_from = Rel_pMek) %>%
  mutate(FC = U0/dmso)

# Paired t-test for upper panels
to.test = data2 %>% filter(Cell_line=='1.8',x_status=='XO')
t.test(to.test$U0,to.test$dmso, paired=T)

to.test = data2 %>% filter(Cell_line=='1.8',x_status=='XX')
t.test(to.test$U0,to.test$dmso, paired=T)

to.test = data2 %>% filter(Cell_line=='E14')
t.test(to.test$U0,to.test$dmso, paired=T)

# one-sample t-test of ratio for lower panel
to.test = data2 %>% filter(Cell_line=='1.8',x_status=='XO')
t.test(to.test$FC,mu=1)

to.test = data2 %>% filter(Cell_line=='1.8',x_status=='XX')
t.test(to.test$FC,mu=1)

to.test = data2 %>% filter(Cell_line=='E14')
t.test(to.test$FC,mu=1)


### pRAF
data = read_excel('./OUTPUT_PAPER/Fig7D_U0126_pRAF1.xls') %>%
  mutate (nr.x = ifelse(x_status=='XX',2,1)) %>% 
  pivot_wider(names_from = Variable, values_from = Value)

data2 = data %>% select(-FCoCntrl,-sample) %>% 
  pivot_wider(names_from=TreatmentType, values_from = Rel_pcRaf) %>%
  mutate(FC = U0/dmso)

# Paired t-test for upper panels
to.test = data2 %>% filter(Cell_line=='1.8',x_status=='XO')
t.test(to.test$U0,to.test$dmso, paired=T)

to.test = data2 %>% filter(Cell_line=='1.8',x_status=='XX')
t.test(to.test$U0,to.test$dmso, paired=T)

to.test = data2 %>% filter(Cell_line=='E14')
t.test(to.test$U0,to.test$dmso, paired=T)

# one-sample t-test of ratio for lower panel
to.test = data2 %>% filter(Cell_line=='1.8',x_status=='XO')
t.test(to.test$FC,mu=1)

to.test = data2 %>% filter(Cell_line=='1.8',x_status=='XX')
t.test(to.test$FC,mu=1)

to.test = data2 %>% filter(Cell_line=='E14')
t.test(to.test$FC,mu=1)


