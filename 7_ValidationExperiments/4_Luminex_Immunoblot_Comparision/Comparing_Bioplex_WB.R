library(gdata)
library(tidyverse)
library(egg)


# Create folder OUTPUT to save figs 
if(!(file.exists("./OUTPUT"))){
  dir.create("./OUTPUT") 
}

source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")


#############################################
########## pMek with Mek inhibitor ##########

Merged_Bioplex_WB <- gdata::read.xls("../../6_SimulationResults_vs_ExpData/INPUTS/Exp_Data/Merged_Log2FC.xls", as.is = TRUE) #as.is=T: To prevent strings to be converted to factors in column Treatment
Fig7CDii_ValidationExp_WB <- gdata::read.xls("../3_Meki_DR_TC/OUTPUT_PAPER/Fig7CDii_PDTC_Mek_Raf_FCoC_dat.xls")

Meki_pMek_Bioplex = Merged_Bioplex_WB %>% 
  filter(Treatment == "Meki") %>% 
  select(Replicate,x_status,Mek_p) %>% 
  mutate(pMek_FC = 2^(Mek_p)) %>% 
  mutate(pMek_LFC = Mek_p) %>% 
  select(-c(Mek_p)) %>% 
  rename(Cell_line = x_status) %>% 
  mutate(Exp=rep("Bioplex", nrow(.)))
Meki_pMek_Bioplex$Replicate = gsub("R3","R1",Meki_pMek_Bioplex$Replicate)
Meki_pMek_Bioplex$Replicate = gsub("R4","R2",Meki_pMek_Bioplex$Replicate)
Meki_pMek_Bioplex$Replicate = gsub("R5","R3",Meki_pMek_Bioplex$Replicate)

  
Meki_pMek_WB = Fig7CDii_ValidationExp_WB %>%  
  filter(Treatment == 0.5) %>% 
  filter(!grepl("pcRaf", Analyte)) %>% 
  select(-c(Treatment_factor,Mean,XX_mean,XO_mean,sig,Analyte,Treatment)) %>% 
  pivot_longer(cols = c("R1","R2","R3"), names_to="Replicate", values_to="pMek_FC") %>% 
  mutate(pMek_LFC = log2(pMek_FC))%>% 
  mutate(Exp=rep("WB", nrow(.))) %>% 
  select(!grep("p_value",colnames(.)))

merged_Meki = full_join(Meki_pMek_Bioplex,Meki_pMek_WB, by= c("Cell_line","Replicate"), suffix = c("_BP","_WB"))

WriteXLS::WriteXLS(merged_Meki, "./OUTPUT/merged_Meki_data.xls")


#######  Plotting fold change values #########
merged_Meki_FC = merged_Meki %>% 
  select(!grep("LFC|Exp", colnames(merged_Meki))) 

# SCATTER PLOT  
# merged_Meki_FC %>% ggplot(aes(x=pMek_FC_WB,y=pMek_FC_BP)) +
#   #geom_point(aes(color=Cell_line, shape=Replicate)) +
#   geom_point(aes(colour = Cell_line, shape = Replicate), size = 0.5)+
#   scale_fill_manual(values = MyCellLineColours) +
#   scale_color_manual(values = MyCellLineColours) +
#   scale_shape_manual(values = MyReplicateShapes, guide=FALSE) +
#   ylim(0,10)+
#   xlim(0,10)+
#   MyPaperTheme

merged_Meki_FC_long = merged_Meki_FC %>% 
  pivot_longer(cols = c("pMek_FC_BP","pMek_FC_WB"), names_to="Exp", values_to="pMek_FoldChange")

merged_Meki_FC_long$Exp <- factor(merged_Meki_FC_long$Exp, levels = c("pMek_FC_BP","pMek_FC_WB"),
                                                  labels = c("Luminex","Immunoblot."))

g = merged_Meki_FC_long %>% ggplot(aes(x=Exp,y=pMek_FoldChange)) +
  #geom_point(aes(color=Cell_line, shape=Replicate)) +
  geom_point(aes(colour = Cell_line), size = 1)+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
  #scale_shape_manual(values = MyReplicateShapes, guide=FALSE) +
  ylim(0,10)+
  labs(x = "Assay", 
       y = "Fold change over untreated",
       title = "pMek",
       color = "Cell line" ) + # color within labs,lets the user give defined labels to the attribute in legend
  MyPaperTheme
gt=egg::set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT/pMek_Mekinhibitor_FC.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./") 


#######  Plotting log2 fold change values #########
merged_Meki_LFC = merged_Meki %>% 
  select(!grep("_FC|Exp", colnames(merged_Meki))) 

merged_Meki_LFC_long = merged_Meki_LFC %>% 
  pivot_longer(cols = c("pMek_LFC_BP","pMek_LFC_WB"), names_to="Exp", values_to="pMek_Log2FoldChange")

merged_Meki_LFC_long$Exp <- factor(merged_Meki_LFC_long$Exp, levels = c("pMek_LFC_BP","pMek_LFC_WB"),
                                                    labels = c("Luminex","Immunoblot."))


g = merged_Meki_LFC_long %>% ggplot(aes(x=Exp,y=pMek_Log2FoldChange)) +
  #geom_point(aes(color=Cell_line, shape=Replicate)) +
  geom_point(aes(colour = Cell_line), size = 1)+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
 # scale_shape_manual(values = MyReplicateShapes, guide=FALSE) +
  ylim(0,4)+
  labs(x = "Assay", 
       y = "Fold change over untreated(log2)",
       title = "pMek",
       color = "Cell line" ) + # color within labs,lets the user give defined labels to the attribute in legend
  MyPaperTheme
gt=egg::set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT/pMek_Mekinhibitor_LFC.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./") 



#############################################
########## pAkt with Gsk3 inhibitor ##########

Gsk3i_pAkt_Bioplex = Merged_Bioplex_WB %>% #head()#colnames()
  filter(Treatment == "Gsk3i") %>% 
  select(Replicate,x_status,Akt_p) %>% 
  mutate(pAkt_FC = 2^(Akt_p)) %>% 
  mutate(pAkt_LFC = Akt_p) %>% 
  select(-c(Akt_p)) %>% 
  rename(Cell_line = x_status) %>% 
  mutate(Exp=rep("Bioplex", nrow(.)))
Gsk3i_pAkt_Bioplex$Replicate = gsub("R3","R1",Gsk3i_pAkt_Bioplex$Replicate)
Gsk3i_pAkt_Bioplex$Replicate = gsub("R4","R2",Gsk3i_pAkt_Bioplex$Replicate)
Gsk3i_pAkt_Bioplex$Replicate = gsub("R5","R3",Gsk3i_pAkt_Bioplex$Replicate)


Fig4C_data_pAKT <- gdata::read.xls("../1_Gsk3i_DR/OUTPUT_PAPER/Fig4C_data_pAKT.xls")

Gsk3i_pAkt_WB = Fig4C_data_pAKT %>% #head()
  filter(Treatment == 6) %>% 
  select(-c(log2Treatment,Mean,p.value,sig,Analyte,Treatment)) %>% 
  pivot_longer(cols = c("R1","R2","R3"), names_to="Replicate", values_to="pAkt_FC") %>% 
  mutate(pAkt_LFC = log2(pAkt_FC))%>% 
  mutate(Exp=rep("WB", nrow(.)))


merged_Gsk3i = full_join(Gsk3i_pAkt_Bioplex,Gsk3i_pAkt_WB, by= c("Cell_line","Replicate"), suffix = c("_BP","_WB"))

WriteXLS::WriteXLS(merged_Gsk3i, "./OUTPUT/merged_Gsk3i_data.xls")


#######  Plotting fold change values #########
merged_Gsk3i_FC = merged_Gsk3i %>% 
  select(!grep("LFC|Exp", colnames(merged_Gsk3i))) 
merged_Gsk3i_FC_long = merged_Gsk3i_FC %>% 
  pivot_longer(cols = c("pAkt_FC_BP","pAkt_FC_WB"), names_to="Exp", values_to="pAkt_FoldChange")
merged_Gsk3i_FC_long = merged_Gsk3i_FC %>% 
  pivot_longer(cols = c("pAkt_FC_BP","pAkt_FC_WB"), names_to="Exp", values_to="pAkt_FoldChange")

merged_Gsk3i_FC_long$Exp = factor(merged_Gsk3i_FC_long$Exp, levels = c("pAkt_FC_BP","pAkt_FC_WB"),
                                  labels = c("Luminex","Immunoblot."))

g = merged_Gsk3i_FC_long %>% ggplot(aes(x=Exp,y=pAkt_FoldChange)) +
  geom_point(aes(colour = Cell_line), size = 1)+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
  ylim(0,1.5)+
  labs(x = "Assay", 
       y = "Fold change over untreated",
       title = "pAkt",
       color = "Cell line" ) + # color within labs,lets the user give defined labels to the attribute in legend
  MyPaperTheme
gt=egg::set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT/pAkt_Gsk3inhibitor_FC.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./") 


#######  Plotting log2 fold change values #########
merged_Gsk3i_LFC = merged_Gsk3i %>% 
  select(!grep("_FC|Exp", colnames(merged_Gsk3i))) 

merged_Gsk3i_LFC_long = merged_Gsk3i_LFC %>% 
  pivot_longer(cols = c("pAkt_LFC_BP","pAkt_LFC_WB"), names_to="Exp", values_to="pAkt_Log2FoldChange")

merged_Gsk3i_LFC_long$Exp <- factor(merged_Gsk3i_LFC_long$Exp, levels = c("pAkt_LFC_BP","pAkt_LFC_WB"),
                                    labels = c("Luminex","Immunoblot."))


g = merged_Gsk3i_LFC_long %>% ggplot(aes(x=Exp,y=pAkt_Log2FoldChange)) +
  geom_point(aes(colour = Cell_line), size = 1)+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
  ylim(-3,0)+
  labs(x = "Assay", 
       y = "Fold change over untreated(log2)",
       title = "pAkt",
       color = "Cell line" ) + 
  MyPaperTheme
gt=egg::set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT/pAkt_Gsk3inhibitor_LFC.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./") 




