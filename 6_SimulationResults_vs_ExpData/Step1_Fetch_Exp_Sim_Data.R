#!/usr/bin/Rscript

library(STASNet)
library(gdata) #read.xls
library(tidyverse) #required for function separate, because that comes from Stringr 

## Create OUTPUT folder to save results
if(!file.exists("OUTPUT")){
dir.create(file.path('OUTPUT'))}


#############################################################
############ Read Initial Model : XX and XO    ############ 
############################################################

Dir_Init_Model = "./INPUTS/Model_Initial/"
XX_Init_Model_file = "XX_R345_Network_100k.mra"
XO_Init_Model_file = "XO_R345_Network_100k.mra"

Dir_MIDAS_Files = "./INPUTS/MIDAS_files/"
XX_MIDAS_file = "XX_R345.csv"
XO_MIDAS_file = "XO_R345.csv"

XX_Init_Model = rebuildModel(file.path(Dir_Init_Model,XX_Init_Model_file),file.path(Dir_MIDAS_Files,XX_MIDAS_file))
XO_Init_Model = rebuildModel(file.path(Dir_Init_Model,XO_Init_Model_file),file.path(Dir_MIDAS_Files,XO_MIDAS_file))

#### Init XX Accuracy data ############
Init_XX_Accuracy <- plotModelAccuracy(XX_Init_Model)
Init_XX_Accuracy_Sim <- Init_XX_Accuracy$simulation
Init_XX_Accuracy_Sim <- data.frame(Init_XX_Accuracy_Sim)
colnames(Init_XX_Accuracy_Sim) <- paste(colnames(Init_XX_Accuracy_Sim), "p", sep = "_")
Init_XX_Accuracy_Sim$x_status <- rep("XX", nrow(Init_XX_Accuracy_Sim))
Init_XX_Accuracy_Sim$id <- rownames(Init_XX_Accuracy_Sim)

Init_XX_Accuracy_Sim <- Init_XX_Accuracy_Sim %>%
  separate(id, c("Perturbation1","Perturbation2"))

#### Init XO Accuracy data ############
Init_XO_Accuracy <- plotModelAccuracy(XO_Init_Model)
Init_XO_Accuracy_Sim <- Init_XO_Accuracy$simulation
Init_XO_Accuracy_Sim <- data.frame(Init_XO_Accuracy_Sim)
colnames(Init_XO_Accuracy_Sim) <- paste(colnames(Init_XO_Accuracy_Sim), "p", sep = "_")
Init_XO_Accuracy_Sim$x_status <- rep("XO", nrow(Init_XO_Accuracy_Sim))
Init_XO_Accuracy_Sim$id <- rownames(Init_XO_Accuracy_Sim)

Init_XO_Accuracy_Sim <- Init_XO_Accuracy_Sim %>%
  separate(id, c("Perturbation1","Perturbation2"))

#### Init model simulation ############
Init_Model_SimulationResults_Points <- bind_rows(Init_XX_Accuracy_Sim,Init_XO_Accuracy_Sim)
Init_Model_SimulationResults_Points$Category <- rep("Sim_Init", nrow(Init_Model_SimulationResults_Points))

write.csv(Init_Model_SimulationResults_Points, "./OUTPUT/Init_Model_SimulationResults_Points.csv")

print("Saved Simulation results from Initial Models")


############################################################
############ Read Final Model : XX and XO    ############ 
############################################################


Dir_Final_Model = "./INPUTS/Model_Completed/"

XX_Final_Model_file = "XX_R345_Network_LA10_100k.mra"
XO_Final_Model_file = "XO_R345_Network_LA10_100k.mra"

XX_Final_Model = rebuildModel(file.path(Dir_Final_Model,XX_Final_Model_file),file.path(Dir_MIDAS_Files,XX_MIDAS_file))
XO_Final_Model = rebuildModel(file.path(Dir_Final_Model,XO_Final_Model_file),file.path(Dir_MIDAS_Files,XO_MIDAS_file))

#### Final XX Accuracy data ############
Final_XX_Accuracy <- plotModelAccuracy(XX_Final_Model)
Final_XX_Accuracy_Sim <- Final_XX_Accuracy$simulation
Final_XX_Accuracy_Sim <- data.frame(Final_XX_Accuracy_Sim)
colnames(Final_XX_Accuracy_Sim) <- paste(colnames(Final_XX_Accuracy_Sim), "p", sep = "_")
Final_XX_Accuracy_Sim$x_status <- rep("XX", nrow(Final_XX_Accuracy_Sim))
Final_XX_Accuracy_Sim$id <- rownames(Final_XX_Accuracy_Sim)

Final_XX_Accuracy_Sim <- Final_XX_Accuracy_Sim %>%
  separate(id, c("Perturbation1","Perturbation2"))

#### Final XO Accuracy data ############
Final_XO_Accuracy <- plotModelAccuracy(XO_Final_Model)
Final_XO_Accuracy_Sim <- Final_XO_Accuracy$simulation
Final_XO_Accuracy_Sim <- data.frame(Final_XO_Accuracy_Sim)
colnames(Final_XO_Accuracy_Sim) <- paste(colnames(Final_XO_Accuracy_Sim), "p", sep = "_")
Final_XO_Accuracy_Sim$x_status <- rep("XO", nrow(Final_XO_Accuracy_Sim))
Final_XO_Accuracy_Sim$id <- rownames(Final_XO_Accuracy_Sim)

Final_XO_Accuracy_Sim <- Final_XO_Accuracy_Sim %>%
  separate(id, c("Perturbation1","Perturbation2"))

#### Final model simulation ############
Final_Model_SimulationResults_Points <- bind_rows(Final_XX_Accuracy_Sim,Final_XO_Accuracy_Sim)
Final_Model_SimulationResults_Points$Category <- rep("Sim_Final", nrow(Final_Model_SimulationResults_Points))

write.csv(Final_Model_SimulationResults_Points, "./OUTPUT/Final_Model_SimulationResults_Points.csv")

print("Saved Simulation results from Completed Models")


######################################################################################################
########## Experimental Data : Merged_Bioplex_WB : Log2FC values(Not in MIDAS format) ##############
######################################################################################################

## Experimental Data : as 3 replicates
ExpData_Folder <-  "./INPUTS/Exp_Data/"

Merged_Bioplex_WB_Points <- read.xls(file.path(ExpData_Folder,"Merged_Log2FC.xls"), as.is = TRUE) #as.is=T: To prevent strings to be converted to factors in column Treatment
# Perturbation2 column in case of single treatments did not have NA for some reason. So added NA in both comlumns first which overwritten by the correct values by using separate.
Merged_Bioplex_WB_Points$Perturbation1 <- NA
Merged_Bioplex_WB_Points$Perturbation2 <- NA

Merged_Bioplex_WB_Points <- Merged_Bioplex_WB_Points %>%
  separate(Treatment, c("Perturbation1", "Perturbation2"),
           sep="\\+", remove = TRUE, # remove = FALSE can be used to keep the Treatment column
           fill = "right") # fill=right To ensure that NAs are filled in the second column, ie Pertrubation2
Merged_Bioplex_WB_Points <- Merged_Bioplex_WB_Points %>%
  mutate(Category = ifelse(x_status == "XX",
                           "Exp_XX", "Exp_XO"))

############### Collect all together ########
my_drop_cols <- c("X", "bCatenin","GEL_Number", "Treatment")
Merged_Bioplex_WB_Points <- Merged_Bioplex_WB_Points[,!(colnames(Merged_Bioplex_WB_Points) %in% my_drop_cols)]
Init_Model_SimulationResults_Points <- Init_Model_SimulationResults_Points[,!(colnames(Init_Model_SimulationResults_Points) %in% my_drop_cols)]
Final_Model_SimulationResults_Points <- Final_Model_SimulationResults_Points[,!(colnames(Final_Model_SimulationResults_Points) %in% my_drop_cols)]

Init_Model_SimulationResults_Points$Replicate <- rep("Init_model",nrow(Init_Model_SimulationResults_Points))
Final_Model_SimulationResults_Points$Replicate <- rep("Final_model",nrow(Final_Model_SimulationResults_Points))

ALL_DATA <- bind_rows(Merged_Bioplex_WB_Points,Init_Model_SimulationResults_Points,Final_Model_SimulationResults_Points)
ALL_DATA$Category <- factor(ALL_DATA$Category, levels=c("Exp_XO","Exp_XX","Sim_Init", "Sim_Final"))


ALL_DATA <- ALL_DATA %>%
  mutate(Treatment_id = paste(Perturbation1,Perturbation2,sep = "+"))

write.csv(ALL_DATA, "./OUTPUT/ALL_DATA.csv")

print("Saved Experimental data + Simulation results")


#########################################################################################
######## Fetching the mismatch values (DataPoints that contribute to Residual) ############
#########################################################################################

#############  XX - Initial Model ###############
##################################################

Init_XX_Accuracy_MM <- Init_XX_Accuracy$mismatch
Init_XX_Accuracy_MM <- data.frame(Init_XX_Accuracy_MM)
colnames(Init_XX_Accuracy_MM) <- paste(colnames(Init_XX_Accuracy_MM), "p", sep = "_")
Init_XX_Accuracy_MM$x_status <- rep("XX", nrow(Init_XX_Accuracy_MM))
Init_XX_Accuracy_MM$Treatment_id <- rownames(Init_XX_Accuracy_MM)

Init_XX_Accuracy_MM <- Init_XX_Accuracy_MM %>%
  separate(Treatment_id, c("Perturbation1","Perturbation2"), remove = FALSE)

#############  XO - Initial Model ###############
##################################################

Init_XO_Accuracy_MM <- Init_XO_Accuracy$mismatch
Init_XO_Accuracy_MM <- data.frame(Init_XO_Accuracy_MM)
colnames(Init_XO_Accuracy_MM) <- paste(colnames(Init_XO_Accuracy_MM), "p", sep = "_")
Init_XO_Accuracy_MM$x_status <- rep("XO", nrow(Init_XO_Accuracy_MM))
Init_XO_Accuracy_MM$Treatment_id <- rownames(Init_XO_Accuracy_MM)

Init_XO_Accuracy_MM <- Init_XO_Accuracy_MM %>%
  separate(Treatment_id, c("Perturbation1","Perturbation2"), remove = FALSE)


#############  XX - Final Model ###############
##################################################
Final_XX_Accuracy_MM <- Final_XX_Accuracy$mismatch
Final_XX_Accuracy_MM <- data.frame(Final_XX_Accuracy_MM)
colnames(Final_XX_Accuracy_MM) <- paste(colnames(Final_XX_Accuracy_MM), "p", sep = "_")
Final_XX_Accuracy_MM$x_status <- rep("XX", nrow(Final_XX_Accuracy_MM))
Final_XX_Accuracy_MM$Treatment_id <- rownames(Final_XX_Accuracy_MM)

Final_XX_Accuracy_MM <- Final_XX_Accuracy_MM %>%
  separate(Treatment_id, c("Perturbation1","Perturbation2"), remove = FALSE)

#############  XO - Final Model ###############
##################################################
Final_XO_Accuracy_MM <- Final_XO_Accuracy$mismatch
Final_XO_Accuracy_MM <- data.frame(Final_XO_Accuracy_MM)
colnames(Final_XO_Accuracy_MM) <- paste(colnames(Final_XO_Accuracy_MM), "p", sep = "_")
Final_XO_Accuracy_MM$x_status <- rep("XO", nrow(Final_XO_Accuracy_MM))
Final_XO_Accuracy_MM$Treatment_id <- rownames(Final_XO_Accuracy_MM)

Final_XO_Accuracy_MM <- Final_XO_Accuracy_MM %>%
  separate(Treatment_id, c("Perturbation1","Perturbation2"), remove = FALSE)


############# Compiling together ###############
##################################################

#### Init model simulation ############
Init_Model_MismatchValues <- bind_rows(Init_XX_Accuracy_MM,Init_XO_Accuracy_MM)
Init_Model_MismatchValues$Category <- rep("Sim_Init", nrow(Init_Model_MismatchValues))

#### Final model simulation ############
Final_Model_MismatchValues <- bind_rows(Final_XX_Accuracy_MM,Final_XO_Accuracy_MM)
Final_Model_MismatchValues$Category <- rep("Sim_Final", nrow(Final_Model_MismatchValues))

Mismatch_DATA <- bind_rows(Init_Model_MismatchValues,Final_Model_MismatchValues)
Mismatch_DATA$Category <- factor(Mismatch_DATA$Category, levels=c("Sim_Init", "Sim_Final"))

Mismatch_DATA <- Mismatch_DATA %>%
  mutate(Treatment_id = paste(Perturbation1,Perturbation2,sep = "+"))

Mismatch_DATA$Treatment_id <- gsub("\\+NA","",Mismatch_DATA$Treatment_id)

write.csv(Mismatch_DATA, "./OUTPUT/Mismatch_DATA.csv")

print("Saved data on mismatch between Experimental results and Simulation results")
print("All steps executed")






