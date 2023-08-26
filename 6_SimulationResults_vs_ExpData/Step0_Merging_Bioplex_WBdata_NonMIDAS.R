#!/usr/bin/Rscript

library(STASNet)
library(gdata) # read.xls
library(openxlsx) # read.xlsx
library(plyr) # for ldply
library(tidyverse) # required for the function separate (package Stringr) 

dir.create("./INPUTS/Exp_Data")

#Read in the normalized FC data from Bioplex and WB
BP_Normalized_folder <- "../2A_Bioplex_Norm_FoldChange/OUTPUT"
BP_FC_STASNet <- read.xlsx(file.path(BP_Normalized_folder,"Bioplex_FC_STASNet.xlsx"))

WB_Normalized_folder <- "../2B_WB_Norm_FoldChange/OUTPUT"
WBN2_GelMean_Norm_FC <- read.xlsx(file.path(WB_Normalized_folder,"WBN2_GelMean_Norm_FC_Statsnet.xlsx"))

## Prepping the Bioplex Data 
BP_FC_STASNet$Treatment <- gsub(" ","",BP_FC_STASNet$Treatment) # Remove spaces
BP_FC_STASNet$Treatment <- gsub("c.\\d+","c",BP_FC_STASNet$Treatment) # To remove the .numbers[1:8] that had been addeded to rownames in case of "control" rows. These digist were aded previously by make.names(unique=TRUE) because rownames could not be duplicated.
BP_FC_STASNet$Treatment <- gsub("^c$","control",BP_FC_STASNet$Treatment)

# Replacing with NA the Akt readouts for XX5 and XO4  because these replicates were inconsistent. 
BP_FC_STASNet$Akt[BP_FC_STASNet$Replicate == "R5" & BP_FC_STASNet$x_status == "XX"] <- NA
BP_FC_STASNet$Akt[BP_FC_STASNet$Replicate == "R4" & BP_FC_STASNet$x_status == "XO"] <- NA

## Prepping the WB Data 
WBN2_GelMean_Norm_FC$Treatment <- gsub("Bmp4i","Bmp4ri",WBN2_GelMean_Norm_FC$Treatment) #Correcting the names of some inhibitors for consistency 
WBN2_GelMean_Norm_FC$Treatment <- gsub("Gsk3bi","Gsk3i",WBN2_GelMean_Norm_FC$Treatment) #Correcting the names of some inhibitors for consistency

WBN2_GelMean_Norm_FC <- WBN2_GelMean_Norm_FC %>% filter(Treatment !="Common") # removing the common samples that had been added for normalization
colnames(WBN2_GelMean_Norm_FC) <- gsub("_FC","",colnames(WBN2_GelMean_Norm_FC))
WBN2_GelMean_Norm_FC <- WBN2_GelMean_Norm_FC %>% 
  select(-GAPDH)

############# Merging the BP and WB data.frames #############
# 1.) Separated the Treatment column into perturbation1 and perturbation2, 
#so that these cane be used to create my id columns, so as to compare both possibilities : P1+P2 and P2+P1.
# 2.) Calculated the Log2 of the fold change values.


#Separate the Treatment into two columns to specify Perturbation1 and Perturbation2

BP_FC_STASNet$Treatment_bk <- BP_FC_STASNet$Treatment
BP_FC_STASNet <- BP_FC_STASNet %>% 
  separate("Treatment_bk", c("Perturbation1","Perturbation2"), sep = "\\+")

WBN2_GelMean_Norm_FC$Treatment_bk <- WBN2_GelMean_Norm_FC$Treatment
WBN2_GelMean_Norm_FC <- WBN2_GelMean_Norm_FC %>% 
  separate("Treatment_bk", c("Perturbation1","Perturbation2"), sep = "\\+")

#Calculating log2 of the fold change values :

BP_LFC <- BP_FC_STASNet %>% 
  mutate_at(vars(Gsk3:Akt),funs(LFC=log2(.)))

WB_LFC <- WBN2_GelMean_Norm_FC %>% 
  mutate_at(vars(ERKp:bCatenin),funs(LFC=log2(.)))

#Rearranging the columns 
BP_LFC <- BP_LFC %>%
  #select(Gsk3:Akt,Gsk3_LFC:Akt_LFC,Replicate,x_status,Treatment,Perturbation1,Perturbation2)
  select(Gsk3_LFC:Akt_LFC,Replicate,x_status,Treatment,Perturbation1,Perturbation2) %>% 
  filter(Treatment != "control")

WB_LFC <- WB_LFC %>% 
  #select(ERKp:bCatenin,ERKp_LFC:bCatenin_LFC,Replicate,x_status,GEL_Number,Treatment,Perturbation1,Perturbation2)
  select(ERKp_LFC:bCatenin_LFC,Replicate,x_status,GEL_Number,Treatment,Perturbation1,Perturbation2) %>% 
  filter(Treatment != "control")

#Create id for merging the BP and WB data :

WB_LFC$id1 <- paste0(WB_LFC$x_status, ".", WB_LFC$Replicate, ".", WB_LFC$Perturbation1, ".", WB_LFC$Perturbation2)
WB_LFC$id2 <- paste0(WB_LFC$x_status, ".", WB_LFC$Replicate, ".", WB_LFC$Perturbation2, ".", WB_LFC$Perturbation1)

BP_LFC$id1 <- paste0(BP_LFC$x_status, ".", BP_LFC$Replicate, ".", BP_LFC$Perturbation1, ".", BP_LFC$Perturbation2)
BP_LFC$id2 <- paste0(BP_LFC$x_status, ".", BP_LFC$Replicate, ".", BP_LFC$Perturbation2, ".", BP_LFC$Perturbation1)


#############

Merged1 = merge(BP_LFC , WB_LFC , by.x=c("id1"), by.y=c("id1"))
Merged1 <- Merged1[,!(grepl("\\.x",colnames(Merged1)))]
Merged1 <- Merged1[,!(grepl("id",colnames(Merged1)))]

Merged2 = merge(BP_LFC , WB_LFC, by.x=c("id1"), by.y=c("id2"))
Merged2 <- Merged2[,!(grepl("\\.x",colnames(Merged2)))]
Merged2 <- Merged2[,!(grepl("id",colnames(Merged2)))]


Merged_Log2FC <- rbind(Merged1, Merged2)
colnames(Merged_Log2FC) <- gsub("\\.y","",colnames(Merged_Log2FC))

#Merged_Log2FC <- Merged_Log2FC[,!(colnames(Merged_Log2FC) %in% c("TPS", "replicate"))]

colnames(Merged_Log2FC) <- gsub("_LFC","",colnames(Merged_Log2FC))

colnames(Merged_Log2FC) <- gsub("Gsk3","Gsk3_p",colnames(Merged_Log2FC))
colnames(Merged_Log2FC) <- gsub("Mek","Mek_p",colnames(Merged_Log2FC))
colnames(Merged_Log2FC) <- gsub("mTor","mTor_p",colnames(Merged_Log2FC))
colnames(Merged_Log2FC) <- gsub("Akt","Akt_p",colnames(Merged_Log2FC))
colnames(Merged_Log2FC) <- gsub("ERKp","Erk_p",colnames(Merged_Log2FC))
colnames(Merged_Log2FC) <- gsub("STAT3p","Stat3_p",colnames(Merged_Log2FC))
colnames(Merged_Log2FC) <- gsub("SMAD2p","Smad2_p",colnames(Merged_Log2FC))
#colnames(Merged_Log2FC) <- gsub("GAPDH","Gapdh",colnames(Merged_Log2FC))


WriteXLS::WriteXLS(Merged_Log2FC, "./INPUTS/Exp_Data/Merged_Log2FC.xls")

####################

Merged_Bioplex_WB <- read.xls("./INPUTS/Exp_Data/Merged_Log2FC.xls", as.is = TRUE) #as.is=T: To prevent strings to be converted to factors in column Treatment
#colnames(Merged_Bioplex_WB) <- gsub("Gsk3","Gsk3_p",colnames(Merged_Bioplex_WB))

perturbations <- unique(as.character(Merged_Bioplex_WB$Treatment)) # 53 perturbations in all
analytes <- colnames(Merged_Bioplex_WB)[1:8]

t_result_XX = list() #prep empty list to store t-test results
t_result_XO = list()


for(perturbation_name in perturbations){
  #perturbation_name = perturbations[2]
  each_perturbation_data <-  Merged_Bioplex_WB %>%
    filter(Treatment == perturbation_name)
  for (analyte in analytes){
    #analyte = analytes[3]
    subdata_XX <- subset(each_perturbation_data,x_status == "XX",select = analyte)
    t_result_XX[[perturbation_name]][[analyte]]$a_mean<- t.test(as.data.frame(subdata_XX),alternative = "two.sided", mu=0)$estimate
    t_result_XX[[perturbation_name]][[analyte]]$a_pvalue<- t.test(as.data.frame(subdata_XX),alternative = "two.sided", mu=0)$p.value
    
    subdata_XO <- subset(each_perturbation_data,x_status == "XO",select = analyte)
    t_result_XO[[perturbation_name]][[analyte]]$a_mean<- t.test(as.data.frame(subdata_XO),alternative = "two.sided", mu=0)$estimate
    t_result_XO[[perturbation_name]][[analyte]]$a_pvalue<- t.test(as.data.frame(subdata_XO),alternative = "two.sided", mu=0)$p.value
  }
}

t_result_XX_df <- ldply (t_result_XX, data.frame) #to convert list to dataframe ldply is from package (plyr)
t_result_XX_df$x_status = rep("XX",nrow(t_result_XX_df))

t_result_XO_df <- ldply (t_result_XO, data.frame)  
t_result_XO_df$x_status = rep("XO",nrow(t_result_XO_df))

t_results_all_df = rbind(t_result_XX_df,t_result_XO_df) # Combining t-test results of XX and XO 


t_results_all_df$Treatment <- t_results_all_df$.id
t_results_all_df <- t_results_all_df %>%
  separate(Treatment, c("Perturbation1", "Perturbation2"), "\\+") 
# spearate the col "Treatment" into two cols "Perturbation1" and "Perturbation2" which will used to crease the 10*10 matrix heatmap
#separate in from  package tidyr

write.csv(t_results_all_df, "./INPUTS/Exp_Data/t_results_all_df.csv")

print("All steps executed : Files saved in ./INPUTS/Exp_Data/")
