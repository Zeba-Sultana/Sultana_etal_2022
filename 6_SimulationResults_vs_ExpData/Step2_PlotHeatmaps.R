#!/usr/bin/Rscript


library(ggplot2)
library(gdata) #read.xls
library(plyr) #load plyr before dplyr
library(dplyr) 
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library("corrplot", lib.loc="~/Library/R/3.3/library") # for correlation matrices
library(scales)
library(cowplot) # for plot_grid
library(egg) #needed for set_panel_size, to get predefined sizes of figure panels.

# For setting theme now using the common theme defined in the validation plot folder :
source("../7_ValidationExperiments/ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")
source("./Functions_PlottingHeatmaps.R")

dir.create(file.path("OUTPUT_PAPER",'Comb_Heatmaps', 'csv_files'), recursive = TRUE)
dir.create(file.path("OUTPUT_PAPER",'Single_Heatmaps', 'csv_files'), recursive = TRUE)

# To define the rownames and colnames we use the Perturbation_Vector
Perturbation_Vector = c("None","Igfri","Pi3ki","Fgf4","Fgfri","Meki","NoLif","Jaki","Activin","Bmp4ri","Gsk3i") 
Treatment_Labels = c("DMSO","IGFRi","PI3Ki","FGF4","FGFRi","MEKi","NoLIF","JAKi","ActA","BMP4Ri","GSK3i") 

## Reading_In_Data
### Exp Data

ExpData_Folder <-  "./INPUTS/Exp_Data/"
#Merged_Bioplex_WB <- read.xls(file.path(ExpData_Folder,"Merged_Log2FC.xls"), as.is = TRUE) #as.is=T: To prevent strings to be converted to factors in column Treatment

t_results_all_df <- read.csv(file.path(ExpData_Folder,"t_results_all_df.csv"))
Exp_Data <- t_results_all_df %>% 
  select(-X)

colnames(Exp_Data) <- gsub("\\.a_mean","",colnames(Exp_Data))
colnames(Exp_Data) <- gsub("\\.id","ID",colnames(Exp_Data))

### SimData
SimulationData_Folder <- "./OUTPUT"

Init_Model_SimulationResults <-read.csv(file.path(SimulationData_Folder,"Init_Model_SimulationResults_Points.csv"))
Init_Model_SimulationResults <- Init_Model_SimulationResults %>% 
  select(-X)

FinalModel_SimulationResults <- read.csv(file.path(SimulationData_Folder,"Final_Model_SimulationResults_Points.csv"))
FinalModel_SimulationResults <- FinalModel_SimulationResults %>% 
  select(-X)


################### Aktp ###################
### Combination Data (Exp Data, Initial and Final Model)

Exp_XX_Aktp <- get_CombinationMatrix_new(Exp_Data, "XX","Akt_p")$Mean_values
Exp_XO_Aktp <- get_CombinationMatrix_new(Exp_Data, "XO","Akt_p")$Mean_values
Exp_XXXO_Aktp_Long <- prepMatrixView(Exp_XX_Aktp,Exp_XO_Aktp)
Exp_XXXO_Aktp_Long$Category <- rep("Exp", nrow(Exp_XXXO_Aktp_Long))

INIT_XX_Aktp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Akt_p")$Mean_values
INIT_XO_Aktp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Akt_p")$Mean_values
INIT_XXXO_Aktp_Long <- prepMatrixView(INIT_XX_Aktp,INIT_XO_Aktp)
INIT_XXXO_Aktp_Long$Category <- rep("Initial_Model", nrow(INIT_XXXO_Aktp_Long))

FINAL_XX_Aktp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Akt_p")$Mean_values
FINAL_XO_Aktp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Akt_p")$Mean_values
FINAL_XXXO_Aktp_Long <- prepMatrixView(FINAL_XX_Aktp,FINAL_XO_Aktp)
FINAL_XXXO_Aktp_Long$Category <- rep("Final_Model", nrow(FINAL_XXXO_Aktp_Long))

COMBO_Data_Aktp <- bind_rows(Exp_XXXO_Aktp_Long,INIT_XXXO_Aktp_Long,FINAL_XXXO_Aktp_Long)
COMBO_Data_Aktp$Category = factor(COMBO_Data_Aktp$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(COMBO_Data_Aktp,"./OUTPUT_PAPER/Comb_Heatmaps/csv_files/COMBO_Data_Aktp.csv")


### Single Data (Exp Data, Initial and Final Model)
SingleP_Exp_XX_Aktp <- get_CombinationMatrix_new(Exp_Data,"XX","Akt_p")$Mean_values_SinglePert
SingleP_Exp_XX_Aktp$x_status <- rep("XX",nrow(SingleP_Exp_XX_Aktp))
SingleP_Exp_XO_Aktp <- get_CombinationMatrix_new(Exp_Data,"XO","Akt_p")$Mean_values_SinglePert
SingleP_Exp_XO_Aktp$x_status <- rep("XO",nrow(SingleP_Exp_XO_Aktp))
SingleP_Exp_XXXO_Aktp <- bind_rows(SingleP_Exp_XX_Aktp,SingleP_Exp_XO_Aktp)
SingleP_Exp_XXXO_Aktp$Category <- rep("Exp", nrow(SingleP_Exp_XXXO_Aktp))

SingleP_INIT_XX_Aktp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Akt_p")$Mean_values_SinglePert
SingleP_INIT_XX_Aktp$x_status <- rep("XX",nrow(SingleP_INIT_XX_Aktp))
SingleP_INIT_XO_Aktp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Akt_p")$Mean_values_SinglePert
SingleP_INIT_XO_Aktp$x_status <- rep("XO",nrow(SingleP_INIT_XO_Aktp))
SingleP_INIT_XXXO_Aktp <- bind_rows(SingleP_INIT_XX_Aktp,SingleP_INIT_XO_Aktp)
SingleP_INIT_XXXO_Aktp$Category <- rep("Initial_Model", nrow(SingleP_INIT_XXXO_Aktp))


SingleP_FINAL_XX_Aktp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Akt_p")$Mean_values_SinglePert
SingleP_FINAL_XX_Aktp$x_status <- rep("XX",nrow(SingleP_FINAL_XX_Aktp))
SingleP_FINAL_XO_Aktp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Akt_p")$Mean_values_SinglePert
SingleP_FINAL_XO_Aktp$x_status <- rep("XO",nrow(SingleP_FINAL_XO_Aktp))
SingleP_FINAL_XXXO_Aktp <- bind_rows(SingleP_FINAL_XX_Aktp,SingleP_FINAL_XO_Aktp)
SingleP_FINAL_XXXO_Aktp$Category <- rep("Final_Model", nrow(SingleP_FINAL_XXXO_Aktp))


SINGLE_Data_Aktp <- bind_rows(SingleP_Exp_XXXO_Aktp,SingleP_INIT_XXXO_Aktp,SingleP_FINAL_XXXO_Aktp)
SINGLE_Data_Aktp$Category = factor(SINGLE_Data_Aktp$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(SINGLE_Data_Aktp,"./OUTPUT_PAPER/Single_Heatmaps/csv_files/SINGLE_Data_Aktp.csv")


##### Plotting Heatmaps ###

#### Combination Heatmap
G1 <- Plot_SmallHeatMap_OnlyComb(COMBO_Data_Aktp,"pAkt")
gt=set_panel_size(G1,width=unit(2,'cm'),height=unit(2,'cm'))
grid.arrange(gt)
ggsave("Fig3A_1_CombHeatmap_Aktp.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Comb_Heatmaps/") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that




#### Single Heatmap
G1 <- Plot_SmallHeatMap_Single_ExpOnly(COMBO_Data_Aktp,SINGLE_Data_Aktp,"pAkt")
gt=set_panel_size(G1,width=unit(1,'cm'),height=unit(4.2,'cm'))
grid.arrange(gt)
ggsave("Fig1C_1_SingleHeatmap_Aktp_ExpOnly.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Single_Heatmaps/")



############### Gsk3p ############### 
### Combination Data (Exp Data, Initial and Final Model)

Exp_XX_Gsk3p <- get_CombinationMatrix_new(Exp_Data, "XX","Gsk3_p")$Mean_values
Exp_XO_Gsk3p <- get_CombinationMatrix_new(Exp_Data, "XO","Gsk3_p")$Mean_values
Exp_XXXO_Gsk3p_Long <- prepMatrixView(Exp_XX_Gsk3p,Exp_XO_Gsk3p)
Exp_XXXO_Gsk3p_Long$Category <- rep("Exp", nrow(Exp_XXXO_Gsk3p_Long))

INIT_XX_Gsk3p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Gsk3_p")$Mean_values
INIT_XO_Gsk3p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Gsk3_p")$Mean_values
INIT_XXXO_Gsk3p_Long <- prepMatrixView(INIT_XX_Gsk3p,INIT_XO_Gsk3p)
INIT_XXXO_Gsk3p_Long$Category <- rep("Initial_Model", nrow(INIT_XXXO_Gsk3p_Long))

FINAL_XX_Gsk3p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Gsk3_p")$Mean_values
FINAL_XO_Gsk3p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Gsk3_p")$Mean_values
FINAL_XXXO_Gsk3p_Long <- prepMatrixView(FINAL_XX_Gsk3p,FINAL_XO_Gsk3p)
FINAL_XXXO_Gsk3p_Long$Category <- rep("Final_Model", nrow(FINAL_XXXO_Gsk3p_Long))

COMBO_Data_Gsk3p <- bind_rows(Exp_XXXO_Gsk3p_Long,INIT_XXXO_Gsk3p_Long,FINAL_XXXO_Gsk3p_Long)
COMBO_Data_Gsk3p$Category = factor(COMBO_Data_Gsk3p$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(COMBO_Data_Gsk3p,"./OUTPUT_PAPER/Comb_Heatmaps/csv_files/COMBO_Data_Gsk3p.csv")


### Single Data (Exp Data, Initial and Final Model)
SingleP_Exp_XX_Gsk3p <- get_CombinationMatrix_new(Exp_Data,"XX","Gsk3_p")$Mean_values_SinglePert
SingleP_Exp_XX_Gsk3p$x_status <- rep("XX",nrow(SingleP_Exp_XX_Gsk3p))
SingleP_Exp_XO_Gsk3p <- get_CombinationMatrix_new(Exp_Data,"XO","Gsk3_p")$Mean_values_SinglePert
SingleP_Exp_XO_Gsk3p$x_status <- rep("XO",nrow(SingleP_Exp_XO_Gsk3p))
SingleP_Exp_XXXO_Gsk3p <- bind_rows(SingleP_Exp_XX_Gsk3p,SingleP_Exp_XO_Gsk3p)
SingleP_Exp_XXXO_Gsk3p$Category <- rep("Exp", nrow(SingleP_Exp_XXXO_Gsk3p))

SingleP_INIT_XX_Gsk3p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Gsk3_p")$Mean_values_SinglePert
SingleP_INIT_XX_Gsk3p$x_status <- rep("XX",nrow(SingleP_INIT_XX_Gsk3p))
SingleP_INIT_XO_Gsk3p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Gsk3_p")$Mean_values_SinglePert
SingleP_INIT_XO_Gsk3p$x_status <- rep("XO",nrow(SingleP_INIT_XO_Gsk3p))
SingleP_INIT_XXXO_Gsk3p <- bind_rows(SingleP_INIT_XX_Gsk3p,SingleP_INIT_XO_Gsk3p)
SingleP_INIT_XXXO_Gsk3p$Category <- rep("Initial_Model", nrow(SingleP_INIT_XXXO_Gsk3p))


SingleP_FINAL_XX_Gsk3p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Gsk3_p")$Mean_values_SinglePert
SingleP_FINAL_XX_Gsk3p$x_status <- rep("XX",nrow(SingleP_FINAL_XX_Gsk3p))
SingleP_FINAL_XO_Gsk3p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Gsk3_p")$Mean_values_SinglePert
SingleP_FINAL_XO_Gsk3p$x_status <- rep("XO",nrow(SingleP_FINAL_XO_Gsk3p))
SingleP_FINAL_XXXO_Gsk3p <- bind_rows(SingleP_FINAL_XX_Gsk3p,SingleP_FINAL_XO_Gsk3p)
SingleP_FINAL_XXXO_Gsk3p$Category <- rep("Final_Model", nrow(SingleP_FINAL_XXXO_Gsk3p))


SINGLE_Data_Gsk3p <- bind_rows(SingleP_Exp_XXXO_Gsk3p,SingleP_INIT_XXXO_Gsk3p,SingleP_FINAL_XXXO_Gsk3p)
SINGLE_Data_Gsk3p$Category = factor(SINGLE_Data_Gsk3p$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(SINGLE_Data_Gsk3p,"./OUTPUT_PAPER/Single_Heatmaps/csv_files/SINGLE_Data_Gsk3p.csv")


##### Plotting Heatmaps ###

#### Combination Heatmap
G1 <- Plot_SmallHeatMap_OnlyComb(COMBO_Data_Gsk3p,"pGsk3")
gt=set_panel_size(G1,width=unit(2,'cm'),height=unit(2,'cm'))
grid.arrange(gt)
ggsave("Fig3A_2_CombHeatmap_Gsk3p.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Comb_Heatmaps/") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that




#### Single Heatmap
G1 <- Plot_SmallHeatMap_Single_ExpOnly(COMBO_Data_Gsk3p,SINGLE_Data_Gsk3p,"pGsk3")
gt=set_panel_size(G1,width=unit(1,'cm'),height=unit(4.2,'cm'))
grid.arrange(gt)
ggsave("Fig1C_2_SingleHeatmap_Gsk3p_ExpOnly.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Single_Heatmaps/")



############### mTorp ############### 
### Combination Data (Exp Data, Initial and Final Model)

Exp_XX_mTorp <- get_CombinationMatrix_new(Exp_Data, "XX","mTor_p")$Mean_values
Exp_XO_mTorp <- get_CombinationMatrix_new(Exp_Data, "XO","mTor_p")$Mean_values
Exp_XXXO_mTorp_Long <- prepMatrixView(Exp_XX_mTorp,Exp_XO_mTorp)
Exp_XXXO_mTorp_Long$Category <- rep("Exp", nrow(Exp_XXXO_mTorp_Long))

INIT_XX_mTorp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","mTor_p")$Mean_values
INIT_XO_mTorp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","mTor_p")$Mean_values
INIT_XXXO_mTorp_Long <- prepMatrixView(INIT_XX_mTorp,INIT_XO_mTorp)
INIT_XXXO_mTorp_Long$Category <- rep("Initial_Model", nrow(INIT_XXXO_mTorp_Long))

FINAL_XX_mTorp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","mTor_p")$Mean_values
FINAL_XO_mTorp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","mTor_p")$Mean_values
FINAL_XXXO_mTorp_Long <- prepMatrixView(FINAL_XX_mTorp,FINAL_XO_mTorp)
FINAL_XXXO_mTorp_Long$Category <- rep("Final_Model", nrow(FINAL_XXXO_mTorp_Long))

COMBO_Data_mTorp <- bind_rows(Exp_XXXO_mTorp_Long,INIT_XXXO_mTorp_Long,FINAL_XXXO_mTorp_Long)
COMBO_Data_mTorp$Category = factor(COMBO_Data_mTorp$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(COMBO_Data_mTorp,"./OUTPUT_PAPER/Comb_Heatmaps/csv_files/COMBO_Data_mTorp.csv")


### Single Data (Exp Data, Initial and Final Model)
SingleP_Exp_XX_mTorp <- get_CombinationMatrix_new(Exp_Data,"XX","mTor_p")$Mean_values_SinglePert
SingleP_Exp_XX_mTorp$x_status <- rep("XX",nrow(SingleP_Exp_XX_mTorp))
SingleP_Exp_XO_mTorp <- get_CombinationMatrix_new(Exp_Data,"XO","mTor_p")$Mean_values_SinglePert
SingleP_Exp_XO_mTorp$x_status <- rep("XO",nrow(SingleP_Exp_XO_mTorp))
SingleP_Exp_XXXO_mTorp <- bind_rows(SingleP_Exp_XX_mTorp,SingleP_Exp_XO_mTorp)
SingleP_Exp_XXXO_mTorp$Category <- rep("Exp", nrow(SingleP_Exp_XXXO_mTorp))

SingleP_INIT_XX_mTorp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","mTor_p")$Mean_values_SinglePert
SingleP_INIT_XX_mTorp$x_status <- rep("XX",nrow(SingleP_INIT_XX_mTorp))
SingleP_INIT_XO_mTorp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","mTor_p")$Mean_values_SinglePert
SingleP_INIT_XO_mTorp$x_status <- rep("XO",nrow(SingleP_INIT_XO_mTorp))
SingleP_INIT_XXXO_mTorp <- bind_rows(SingleP_INIT_XX_mTorp,SingleP_INIT_XO_mTorp)
SingleP_INIT_XXXO_mTorp$Category <- rep("Initial_Model", nrow(SingleP_INIT_XXXO_mTorp))


SingleP_FINAL_XX_mTorp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","mTor_p")$Mean_values_SinglePert
SingleP_FINAL_XX_mTorp$x_status <- rep("XX",nrow(SingleP_FINAL_XX_mTorp))
SingleP_FINAL_XO_mTorp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","mTor_p")$Mean_values_SinglePert
SingleP_FINAL_XO_mTorp$x_status <- rep("XO",nrow(SingleP_FINAL_XO_mTorp))
SingleP_FINAL_XXXO_mTorp <- bind_rows(SingleP_FINAL_XX_mTorp,SingleP_FINAL_XO_mTorp)
SingleP_FINAL_XXXO_mTorp$Category <- rep("Final_Model", nrow(SingleP_FINAL_XXXO_mTorp))


SINGLE_Data_mTorp <- bind_rows(SingleP_Exp_XXXO_mTorp,SingleP_INIT_XXXO_mTorp,SingleP_FINAL_XXXO_mTorp)
SINGLE_Data_mTorp$Category = factor(SINGLE_Data_mTorp$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(SINGLE_Data_mTorp,"./OUTPUT_PAPER/Single_Heatmaps/csv_files/SINGLE_Data_mTorp.csv")


##### Plotting Heatmaps ###

#### Combination Heatmap
G1 <- Plot_SmallHeatMap_OnlyComb(COMBO_Data_mTorp,"pmTor")
gt=set_panel_size(G1,width=unit(2,'cm'),height=unit(2,'cm'))
grid.arrange(gt)
ggsave("Fig3A_3_CombHeatmap_mTorp.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Comb_Heatmaps/") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that




#### Single Heatmap
G1 <- Plot_SmallHeatMap_Single_ExpOnly(COMBO_Data_mTorp,SINGLE_Data_mTorp,"pmTor")
gt=set_panel_size(G1,width=unit(1,'cm'),height=unit(4.2,'cm'))
grid.arrange(gt)
ggsave("Fig1C_3_SingleHeatmap_mTorp_ExpOnly.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Single_Heatmaps/")




################### Mekp ###################
### Combination Data (Exp Data, Initial and Final Model)

Exp_XX_Mekp <- get_CombinationMatrix_new(Exp_Data, "XX","Mek_p")$Mean_values
Exp_XO_Mekp <- get_CombinationMatrix_new(Exp_Data, "XO","Mek_p")$Mean_values
Exp_XXXO_Mekp_Long <- prepMatrixView(Exp_XX_Mekp,Exp_XO_Mekp)
Exp_XXXO_Mekp_Long$Category <- rep("Exp", nrow(Exp_XXXO_Mekp_Long))

INIT_XX_Mekp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Mek_p")$Mean_values
INIT_XO_Mekp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Mek_p")$Mean_values
INIT_XXXO_Mekp_Long <- prepMatrixView(INIT_XX_Mekp,INIT_XO_Mekp)
INIT_XXXO_Mekp_Long$Category <- rep("Initial_Model", nrow(INIT_XXXO_Mekp_Long))

FINAL_XX_Mekp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Mek_p")$Mean_values
FINAL_XO_Mekp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Mek_p")$Mean_values
FINAL_XXXO_Mekp_Long <- prepMatrixView(FINAL_XX_Mekp,FINAL_XO_Mekp)
FINAL_XXXO_Mekp_Long$Category <- rep("Final_Model", nrow(FINAL_XXXO_Mekp_Long))

COMBO_Data_Mekp <- bind_rows(Exp_XXXO_Mekp_Long,INIT_XXXO_Mekp_Long,FINAL_XXXO_Mekp_Long)
COMBO_Data_Mekp$Category = factor(COMBO_Data_Mekp$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(COMBO_Data_Mekp,"./OUTPUT_PAPER/Comb_Heatmaps/csv_files/COMBO_Data_Mekp.csv")


### Single Data (Exp Data, Initial and Final Model)
SingleP_Exp_XX_Mekp <- get_CombinationMatrix_new(Exp_Data,"XX","Mek_p")$Mean_values_SinglePert
SingleP_Exp_XX_Mekp$x_status <- rep("XX",nrow(SingleP_Exp_XX_Mekp))
SingleP_Exp_XO_Mekp <- get_CombinationMatrix_new(Exp_Data,"XO","Mek_p")$Mean_values_SinglePert
SingleP_Exp_XO_Mekp$x_status <- rep("XO",nrow(SingleP_Exp_XO_Mekp))
SingleP_Exp_XXXO_Mekp <- bind_rows(SingleP_Exp_XX_Mekp,SingleP_Exp_XO_Mekp)
SingleP_Exp_XXXO_Mekp$Category <- rep("Exp", nrow(SingleP_Exp_XXXO_Mekp))

SingleP_INIT_XX_Mekp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Mek_p")$Mean_values_SinglePert
SingleP_INIT_XX_Mekp$x_status <- rep("XX",nrow(SingleP_INIT_XX_Mekp))
SingleP_INIT_XO_Mekp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Mek_p")$Mean_values_SinglePert
SingleP_INIT_XO_Mekp$x_status <- rep("XO",nrow(SingleP_INIT_XO_Mekp))
SingleP_INIT_XXXO_Mekp <- bind_rows(SingleP_INIT_XX_Mekp,SingleP_INIT_XO_Mekp)
SingleP_INIT_XXXO_Mekp$Category <- rep("Initial_Model", nrow(SingleP_INIT_XXXO_Mekp))


SingleP_FINAL_XX_Mekp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Mek_p")$Mean_values_SinglePert
SingleP_FINAL_XX_Mekp$x_status <- rep("XX",nrow(SingleP_FINAL_XX_Mekp))
SingleP_FINAL_XO_Mekp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Mek_p")$Mean_values_SinglePert
SingleP_FINAL_XO_Mekp$x_status <- rep("XO",nrow(SingleP_FINAL_XO_Mekp))
SingleP_FINAL_XXXO_Mekp <- bind_rows(SingleP_FINAL_XX_Mekp,SingleP_FINAL_XO_Mekp)
SingleP_FINAL_XXXO_Mekp$Category <- rep("Final_Model", nrow(SingleP_FINAL_XXXO_Mekp))


SINGLE_Data_Mekp <- bind_rows(SingleP_Exp_XXXO_Mekp,SingleP_INIT_XXXO_Mekp,SingleP_FINAL_XXXO_Mekp)
SINGLE_Data_Mekp$Category = factor(SINGLE_Data_Mekp$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(SINGLE_Data_Mekp,"./OUTPUT_PAPER/Single_Heatmaps/csv_files/SINGLE_Data_Mekp.csv")


##### Plotting Heatmaps ###

#### Combination Heatmap
G1 <- Plot_SmallHeatMap_OnlyComb(COMBO_Data_Mekp,"pMek")
gt=set_panel_size(G1,width=unit(2,'cm'),height=unit(2,'cm'))
grid.arrange(gt)
ggsave("Fig3A_4_CombHeatmap_Mekp.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Comb_Heatmaps/") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that




#### Single Heatmap
G1 <- Plot_SmallHeatMap_Single_ExpOnly(COMBO_Data_Mekp,SINGLE_Data_Mekp,"pMek")
gt=set_panel_size(G1,width=unit(1,'cm'),height=unit(4.2,'cm'))
grid.arrange(gt)
ggsave("Fig1C_4_SingleHeatmap_Mekp_ExpOnly.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Single_Heatmaps/")

############# Erkp ############# 
### Combination Data (Exp Data, Initial and Final Model)

Exp_XX_Erkp <- get_CombinationMatrix_new(Exp_Data, "XX","Erk_p")$Mean_values
Exp_XO_Erkp <- get_CombinationMatrix_new(Exp_Data, "XO","Erk_p")$Mean_values
Exp_XXXO_Erkp_Long <- prepMatrixView(Exp_XX_Erkp,Exp_XO_Erkp)
Exp_XXXO_Erkp_Long$Category <- rep("Exp", nrow(Exp_XXXO_Erkp_Long))

INIT_XX_Erkp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Erk_p")$Mean_values
INIT_XO_Erkp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Erk_p")$Mean_values
INIT_XXXO_Erkp_Long <- prepMatrixView(INIT_XX_Erkp,INIT_XO_Erkp)
INIT_XXXO_Erkp_Long$Category <- rep("Initial_Model", nrow(INIT_XXXO_Erkp_Long))

FINAL_XX_Erkp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Erk_p")$Mean_values
FINAL_XO_Erkp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Erk_p")$Mean_values
FINAL_XXXO_Erkp_Long <- prepMatrixView(FINAL_XX_Erkp,FINAL_XO_Erkp)
FINAL_XXXO_Erkp_Long$Category <- rep("Final_Model", nrow(FINAL_XXXO_Erkp_Long))

COMBO_Data_Erkp <- bind_rows(Exp_XXXO_Erkp_Long,INIT_XXXO_Erkp_Long,FINAL_XXXO_Erkp_Long)
COMBO_Data_Erkp$Category = factor(COMBO_Data_Erkp$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(COMBO_Data_Erkp,"./OUTPUT_PAPER/Comb_Heatmaps/csv_files/COMBO_Data_Erkp.csv")


### Single Data (Exp Data, Initial and Final Model)
SingleP_Exp_XX_Erkp <- get_CombinationMatrix_new(Exp_Data,"XX","Erk_p")$Mean_values_SinglePert
SingleP_Exp_XX_Erkp$x_status <- rep("XX",nrow(SingleP_Exp_XX_Erkp))
SingleP_Exp_XO_Erkp <- get_CombinationMatrix_new(Exp_Data,"XO","Erk_p")$Mean_values_SinglePert
SingleP_Exp_XO_Erkp$x_status <- rep("XO",nrow(SingleP_Exp_XO_Erkp))
SingleP_Exp_XXXO_Erkp <- bind_rows(SingleP_Exp_XX_Erkp,SingleP_Exp_XO_Erkp)
SingleP_Exp_XXXO_Erkp$Category <- rep("Exp", nrow(SingleP_Exp_XXXO_Erkp))

SingleP_INIT_XX_Erkp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Erk_p")$Mean_values_SinglePert
SingleP_INIT_XX_Erkp$x_status <- rep("XX",nrow(SingleP_INIT_XX_Erkp))
SingleP_INIT_XO_Erkp <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Erk_p")$Mean_values_SinglePert
SingleP_INIT_XO_Erkp$x_status <- rep("XO",nrow(SingleP_INIT_XO_Erkp))
SingleP_INIT_XXXO_Erkp <- bind_rows(SingleP_INIT_XX_Erkp,SingleP_INIT_XO_Erkp)
SingleP_INIT_XXXO_Erkp$Category <- rep("Initial_Model", nrow(SingleP_INIT_XXXO_Erkp))


SingleP_FINAL_XX_Erkp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Erk_p")$Mean_values_SinglePert
SingleP_FINAL_XX_Erkp$x_status <- rep("XX",nrow(SingleP_FINAL_XX_Erkp))
SingleP_FINAL_XO_Erkp <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Erk_p")$Mean_values_SinglePert
SingleP_FINAL_XO_Erkp$x_status <- rep("XO",nrow(SingleP_FINAL_XO_Erkp))
SingleP_FINAL_XXXO_Erkp <- bind_rows(SingleP_FINAL_XX_Erkp,SingleP_FINAL_XO_Erkp)
SingleP_FINAL_XXXO_Erkp$Category <- rep("Final_Model", nrow(SingleP_FINAL_XXXO_Erkp))


SINGLE_Data_Erkp <- bind_rows(SingleP_Exp_XXXO_Erkp,SingleP_INIT_XXXO_Erkp,SingleP_FINAL_XXXO_Erkp)
SINGLE_Data_Erkp$Category = factor(SINGLE_Data_Erkp$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(SINGLE_Data_Erkp,"./OUTPUT_PAPER/Single_Heatmaps/csv_files/SINGLE_Data_Erkp.csv")


##### Plotting Heatmaps ###

#### Combination Heatmap
G1 <- Plot_SmallHeatMap_OnlyComb(COMBO_Data_Erkp,"pErk")
gt=set_panel_size(G1,width=unit(2,'cm'),height=unit(2,'cm'))
grid.arrange(gt)
ggsave("Fig3A_5_CombHeatmap_Erkp.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Comb_Heatmaps/") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that




#### Single Heatmap
G1 <- Plot_SmallHeatMap_Single_ExpOnly(COMBO_Data_Erkp,SINGLE_Data_Erkp,"pErk")
gt=set_panel_size(G1,width=unit(1,'cm'),height=unit(4.2,'cm'))
grid.arrange(gt)
ggsave("Fig1C_5_SingleHeatmap_Erkp_ExpOnly.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Single_Heatmaps/")

############# Stat3p ############# 
### Combination Data (Exp Data, Initial and Final Model)

Exp_XX_Stat3p <- get_CombinationMatrix_new(Exp_Data, "XX","Stat3_p")$Mean_values
Exp_XO_Stat3p <- get_CombinationMatrix_new(Exp_Data, "XO","Stat3_p")$Mean_values
Exp_XXXO_Stat3p_Long <- prepMatrixView(Exp_XX_Stat3p,Exp_XO_Stat3p)
Exp_XXXO_Stat3p_Long$Category <- rep("Exp", nrow(Exp_XXXO_Stat3p_Long))

INIT_XX_Stat3p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Stat3_p")$Mean_values
INIT_XO_Stat3p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Stat3_p")$Mean_values
INIT_XXXO_Stat3p_Long <- prepMatrixView(INIT_XX_Stat3p,INIT_XO_Stat3p)
INIT_XXXO_Stat3p_Long$Category <- rep("Initial_Model", nrow(INIT_XXXO_Stat3p_Long))

FINAL_XX_Stat3p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Stat3_p")$Mean_values
FINAL_XO_Stat3p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Stat3_p")$Mean_values
FINAL_XXXO_Stat3p_Long <- prepMatrixView(FINAL_XX_Stat3p,FINAL_XO_Stat3p)
FINAL_XXXO_Stat3p_Long$Category <- rep("Final_Model", nrow(FINAL_XXXO_Stat3p_Long))

COMBO_Data_Stat3p <- bind_rows(Exp_XXXO_Stat3p_Long,INIT_XXXO_Stat3p_Long,FINAL_XXXO_Stat3p_Long)
COMBO_Data_Stat3p$Category = factor(COMBO_Data_Stat3p$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(COMBO_Data_Stat3p,"./OUTPUT_PAPER/Comb_Heatmaps/csv_files/COMBO_Data_Stat3p.csv")


### Single Data (Exp Data, Initial and Final Model)
SingleP_Exp_XX_Stat3p <- get_CombinationMatrix_new(Exp_Data,"XX","Stat3_p")$Mean_values_SinglePert
SingleP_Exp_XX_Stat3p$x_status <- rep("XX",nrow(SingleP_Exp_XX_Stat3p))
SingleP_Exp_XO_Stat3p <- get_CombinationMatrix_new(Exp_Data,"XO","Stat3_p")$Mean_values_SinglePert
SingleP_Exp_XO_Stat3p$x_status <- rep("XO",nrow(SingleP_Exp_XO_Stat3p))
SingleP_Exp_XXXO_Stat3p <- bind_rows(SingleP_Exp_XX_Stat3p,SingleP_Exp_XO_Stat3p)
SingleP_Exp_XXXO_Stat3p$Category <- rep("Exp", nrow(SingleP_Exp_XXXO_Stat3p))

SingleP_INIT_XX_Stat3p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Stat3_p")$Mean_values_SinglePert
SingleP_INIT_XX_Stat3p$x_status <- rep("XX",nrow(SingleP_INIT_XX_Stat3p))
SingleP_INIT_XO_Stat3p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Stat3_p")$Mean_values_SinglePert
SingleP_INIT_XO_Stat3p$x_status <- rep("XO",nrow(SingleP_INIT_XO_Stat3p))
SingleP_INIT_XXXO_Stat3p <- bind_rows(SingleP_INIT_XX_Stat3p,SingleP_INIT_XO_Stat3p)
SingleP_INIT_XXXO_Stat3p$Category <- rep("Initial_Model", nrow(SingleP_INIT_XXXO_Stat3p))


SingleP_FINAL_XX_Stat3p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Stat3_p")$Mean_values_SinglePert
SingleP_FINAL_XX_Stat3p$x_status <- rep("XX",nrow(SingleP_FINAL_XX_Stat3p))
SingleP_FINAL_XO_Stat3p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Stat3_p")$Mean_values_SinglePert
SingleP_FINAL_XO_Stat3p$x_status <- rep("XO",nrow(SingleP_FINAL_XO_Stat3p))
SingleP_FINAL_XXXO_Stat3p <- bind_rows(SingleP_FINAL_XX_Stat3p,SingleP_FINAL_XO_Stat3p)
SingleP_FINAL_XXXO_Stat3p$Category <- rep("Final_Model", nrow(SingleP_FINAL_XXXO_Stat3p))


SINGLE_Data_Stat3p <- bind_rows(SingleP_Exp_XXXO_Stat3p,SingleP_INIT_XXXO_Stat3p,SingleP_FINAL_XXXO_Stat3p)
SINGLE_Data_Stat3p$Category = factor(SINGLE_Data_Stat3p$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(SINGLE_Data_Stat3p,"./OUTPUT_PAPER/Single_Heatmaps/csv_files/SINGLE_Data_Stat3p.csv")


##### Plotting Heatmaps ###

#### Combination Heatmap
G1 <- Plot_SmallHeatMap_OnlyComb(COMBO_Data_Stat3p,"pStat3")
gt=set_panel_size(G1,width=unit(2,'cm'),height=unit(2,'cm'))
grid.arrange(gt)
ggsave("Fig3A_6_CombHeatmap_Stat3p.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Comb_Heatmaps/") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that




#### Single Heatmap
G1 <- Plot_SmallHeatMap_Single_ExpOnly(COMBO_Data_Stat3p,SINGLE_Data_Stat3p,"pStat3")
gt=set_panel_size(G1,width=unit(1,'cm'),height=unit(4.2,'cm'))
grid.arrange(gt)
ggsave("Fig1C_6_SingleHeatmap_Stat3p_ExpOnly.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Single_Heatmaps/")


############ Smad2p ############ 
### Combination Data (Exp Data, Initial and Final Model)

Exp_XX_Smad2p <- get_CombinationMatrix_new(Exp_Data, "XX","Smad2_p")$Mean_values
Exp_XO_Smad2p <- get_CombinationMatrix_new(Exp_Data, "XO","Smad2_p")$Mean_values
Exp_XXXO_Smad2p_Long <- prepMatrixView(Exp_XX_Smad2p,Exp_XO_Smad2p)
Exp_XXXO_Smad2p_Long$Category <- rep("Exp", nrow(Exp_XXXO_Smad2p_Long))

INIT_XX_Smad2p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Smad2_p")$Mean_values
INIT_XO_Smad2p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Smad2_p")$Mean_values
INIT_XXXO_Smad2p_Long <- prepMatrixView(INIT_XX_Smad2p,INIT_XO_Smad2p)
INIT_XXXO_Smad2p_Long$Category <- rep("Initial_Model", nrow(INIT_XXXO_Smad2p_Long))

FINAL_XX_Smad2p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Smad2_p")$Mean_values
FINAL_XO_Smad2p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Smad2_p")$Mean_values
FINAL_XXXO_Smad2p_Long <- prepMatrixView(FINAL_XX_Smad2p,FINAL_XO_Smad2p)
FINAL_XXXO_Smad2p_Long$Category <- rep("Final_Model", nrow(FINAL_XXXO_Smad2p_Long))

COMBO_Data_Smad2p <- bind_rows(Exp_XXXO_Smad2p_Long,INIT_XXXO_Smad2p_Long,FINAL_XXXO_Smad2p_Long)
COMBO_Data_Smad2p$Category = factor(COMBO_Data_Smad2p$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(COMBO_Data_Smad2p,"./OUTPUT_PAPER/Comb_Heatmaps/csv_files/COMBO_Data_Smad2p.csv")


### Single Data (Exp Data, Initial and Final Model)
SingleP_Exp_XX_Smad2p <- get_CombinationMatrix_new(Exp_Data,"XX","Smad2_p")$Mean_values_SinglePert
SingleP_Exp_XX_Smad2p$x_status <- rep("XX",nrow(SingleP_Exp_XX_Smad2p))
SingleP_Exp_XO_Smad2p <- get_CombinationMatrix_new(Exp_Data,"XO","Smad2_p")$Mean_values_SinglePert
SingleP_Exp_XO_Smad2p$x_status <- rep("XO",nrow(SingleP_Exp_XO_Smad2p))
SingleP_Exp_XXXO_Smad2p <- bind_rows(SingleP_Exp_XX_Smad2p,SingleP_Exp_XO_Smad2p)
SingleP_Exp_XXXO_Smad2p$Category <- rep("Exp", nrow(SingleP_Exp_XXXO_Smad2p))

SingleP_INIT_XX_Smad2p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XX","Smad2_p")$Mean_values_SinglePert
SingleP_INIT_XX_Smad2p$x_status <- rep("XX",nrow(SingleP_INIT_XX_Smad2p))
SingleP_INIT_XO_Smad2p <- get_CombinationMatrix_Simulation(Init_Model_SimulationResults, "XO","Smad2_p")$Mean_values_SinglePert
SingleP_INIT_XO_Smad2p$x_status <- rep("XO",nrow(SingleP_INIT_XO_Smad2p))
SingleP_INIT_XXXO_Smad2p <- bind_rows(SingleP_INIT_XX_Smad2p,SingleP_INIT_XO_Smad2p)
SingleP_INIT_XXXO_Smad2p$Category <- rep("Initial_Model", nrow(SingleP_INIT_XXXO_Smad2p))


SingleP_FINAL_XX_Smad2p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XX","Smad2_p")$Mean_values_SinglePert
SingleP_FINAL_XX_Smad2p$x_status <- rep("XX",nrow(SingleP_FINAL_XX_Smad2p))
SingleP_FINAL_XO_Smad2p <- get_CombinationMatrix_Simulation(FinalModel_SimulationResults, "XO","Smad2_p")$Mean_values_SinglePert
SingleP_FINAL_XO_Smad2p$x_status <- rep("XO",nrow(SingleP_FINAL_XO_Smad2p))
SingleP_FINAL_XXXO_Smad2p <- bind_rows(SingleP_FINAL_XX_Smad2p,SingleP_FINAL_XO_Smad2p)
SingleP_FINAL_XXXO_Smad2p$Category <- rep("Final_Model", nrow(SingleP_FINAL_XXXO_Smad2p))


SINGLE_Data_Smad2p <- bind_rows(SingleP_Exp_XXXO_Smad2p,SingleP_INIT_XXXO_Smad2p,SingleP_FINAL_XXXO_Smad2p)
SINGLE_Data_Smad2p$Category = factor(SINGLE_Data_Smad2p$Category, levels = c("Exp", "Initial_Model", "Final_Model"))

write.csv(SINGLE_Data_Smad2p,"./OUTPUT_PAPER/Single_Heatmaps/csv_files/SINGLE_Data_Smad2p.csv")


##### Plotting Heatmaps ###

#### Combination Heatmap
G1 <- Plot_SmallHeatMap_OnlyComb(COMBO_Data_Smad2p,"pSmad2")
gt=set_panel_size(G1,width=unit(2,'cm'),height=unit(2,'cm'))
grid.arrange(gt)
ggsave("Fig3A_7_CombHeatmap_Smad2p.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Comb_Heatmaps/") # device = cairo_pdf was tried to get the greek letter in pdf, but not yet successful in that




#### Single Heatmap
G1 <- Plot_SmallHeatMap_Single_ExpOnly(COMBO_Data_Smad2p,SINGLE_Data_Smad2p,"pSmad2")
gt=set_panel_size(G1,width=unit(1,'cm'),height=unit(4.2,'cm'))
grid.arrange(gt)
ggsave("Fig1C_7_SingleHeatmap_Smad2p_ExpOnly.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER/Single_Heatmaps/")


print("All heatmaps saved")




















