#!/usr/bin/env Rscript


library(tidyr)
library(stringr)
library(openxlsx)
library(ggplot2)
library(tidyverse) 

##### OUTPUT folder ####
dir.create(file.path('OUTPUT'))
dir.create(file.path('OUTPUT_PAPER'))

##### Functions ##########

#Removes columns for raw data and of the normalization over GAPDH. Just retains Normalization over TPS and the annotation columns.
WB_alldata_to_TNorm <- function(WB_alldata){
  
  #WB_alldata <- Replicate3_alldata #TODEBUG
  
  WB_alldata <- WB_alldata %>% filter(!Treatment == "Common2")
  WB_alldata$Treatment <- gsub("Common1","Common",WB_alldata$Treatment)
  
  annotation_cols <- c("GEL_Number","Replicate","x_status","TPS_TNorm","Treatment")
  rawdata_cols <- c("ERKp","STAT3p","TPS","SMAD2p","GAPDH", "bCatenin")
  
  WB_TNormdata <- WB_alldata %>% 
    select(-grep("GNorm",colnames(WB_alldata))) %>%
    select(-one_of(rawdata_cols))
  
  WB_TNormdata$Treatment <- gsub("Cntrl","control",WB_TNormdata$Treatment)
  
  return(WB_TNormdata)
  
}

##For Normalization Strategy 2 :
Calc_MeanperGel_perReplicate <- function(WB_TNorm,Replicate_number){
  
  WB_TNorm <- WB_TNorm %>% mutate(id_col = paste0(GEL_Number,"_",Treatment,"_",x_status))
  
  #Find rows that have NA values :
  WB_TNorm_NArows <-  WB_TNorm[!complete.cases(WB_TNorm),]
  
  #To drop these treatments from the other replicates as well, I create this vector from where id_col values have to be matched :
  WB_TNorm_NArows_ids <- WB_TNorm_NArows$id_col
  
  WB_TNorm_subset <- 
    WB_TNorm %>%
    filter(!(id_col %in% WB_TNorm_NArows_ids))
  
  #To make it per_replicate
  WB_TNorm_subset <- 
    WB_TNorm_subset %>% 
    filter(Replicate == Replicate_number)
  
  #Both the XX and XO samples that were on each gel have been retained for calculating the mean.
  WB_TNorm_subset_GelMean <- 
    WB_TNorm_subset %>% 
    group_by(GEL_Number) %>% 
    summarise_at(.vars = vars(ERKp:GAPDH),
                 .funs = funs(mean(.))) %>% 
    rename_at(.vars = vars(ERKp:GAPDH),
              .funs = funs(paste0(.,"_mean")))
  
  WB_TNorm_subset_GelMean$Replicate <- rep(Replicate_number,nrow(WB_TNorm_subset_GelMean))
  
  return(WB_TNorm_subset_GelMean)
}

WBN2_overGelMean <- function(WB_TNorm,MeanperGel_perReplicate){
  
  Replicate_number <- unique(MeanperGel_perReplicate$Replicate)
  Replicate_WB_TNorm <- 
    WB_TNorm %>% 
    filter(Replicate == Replicate_number)
  
  MeanperGel_perReplicate <- MeanperGel_perReplicate %>% 
    select(-Replicate) #Removing the replicate column because Replicate column is in both the data frame and hence one needs to be removed before leftjoin else gets a funny suffix. 
  
  ReplicateSpecific_WBdata_Gelmeanvalues <- left_join(Replicate_WB_TNorm,MeanperGel_perReplicate,by= "GEL_Number")
  
  # annotation_cols <- c("GEL_Number","Treatment","x_status","Replicate")
  # All_cols <- colnames(WB_data_Gelmeanvalues)
  
  
  MeanperGel_perReplicate_Normalized <- ReplicateSpecific_WBdata_Gelmeanvalues %>% 
    transmute(GEL_Number, # Transmute retains all the columns that are passed to it.  
              Treatment, #So u can list the columns that u want retained even when u dont use them to calculate new columns.
              x_status,
              Replicate,
              ERKp_Norm = ERKp/ERKp_mean,
              STAT3p_Norm = STAT3p/STAT3p_mean,
              SMAD2p_Norm = SMAD2p/SMAD2p_mean,
              bCatenin_Norm = bCatenin/bCatenin_mean,
              GAPDH_Norm = GAPDH/GAPDH_mean)
  
  return(MeanperGel_perReplicate_Normalized)
  
}

##### Input Data ##########

WB_folder <- "../1B_WB_SignalQuantification/OUTPUT"


Replicate3_alldata <- openxlsx::read.xlsx(file.path(WB_folder,"Replicate3_Alldata.xlsx"))
Replicate4_alldata <- openxlsx::read.xlsx(file.path(WB_folder,"Replicate4_Alldata.xlsx"))
Replicate5_alldata <- openxlsx::read.xlsx(file.path(WB_folder,"Replicate5_Alldata.xlsx"))

WB_alldata <- rbind(Replicate3_alldata,Replicate4_alldata,Replicate5_alldata)

WB_TNorm <- WB_alldata_to_TNorm(WB_alldata)
colnames(WB_TNorm) <- gsub("_TNorm", "", colnames(WB_TNorm))


#####  WBN2 : Normalization over mean signal per gel : ##########
R3_MeanperGel_perReplicate <- Calc_MeanperGel_perReplicate(WB_TNorm,"R3")
R4_MeanperGel_perReplicate <- Calc_MeanperGel_perReplicate(WB_TNorm,"R4")
R5_MeanperGel_perReplicate <- Calc_MeanperGel_perReplicate(WB_TNorm,"R5")


WBN2_R3_GelMean_Norm <- WBN2_overGelMean(WB_TNorm,R3_MeanperGel_perReplicate)
WBN2_R4_GelMean_Norm <- WBN2_overGelMean(WB_TNorm,R4_MeanperGel_perReplicate)
WBN2_R5_GelMean_Norm <- WBN2_overGelMean(WB_TNorm,R5_MeanperGel_perReplicate)

WBN2_GelMean_Norm <- rbind(WBN2_R3_GelMean_Norm,WBN2_R4_GelMean_Norm,WBN2_R5_GelMean_Norm)
colnames(WBN2_GelMean_Norm) <- gsub("_Norm","",colnames(WBN2_GelMean_Norm))


############ FOR_PAPER ##############

WBData_MeanNormalized <- WBN2_GelMean_Norm %>% 
  select(x_status,Replicate,Treatment,GEL_Number,ERKp,STAT3p,SMAD2p)

WBData_MeanNormalized$Treatment <- gsub("control","DMSO",WBData_MeanNormalized$Treatment )
WBData_MeanNormalized$Replicate <- gsub("R3","Rep1",WBData_MeanNormalized$Replicate )
WBData_MeanNormalized$Replicate <- gsub("R4","Rep2",WBData_MeanNormalized$Replicate )
WBData_MeanNormalized$Replicate <- gsub("R5","Rep3",WBData_MeanNormalized$Replicate )

WBData_MeanNormalized <- WBData_MeanNormalized %>% 
  filter(Treatment != "Common")

openxlsx::write.xlsx(WBData_MeanNormalized, "./OUTPUT_PAPER/Suppl_Table_S3_2_WBData_MeanNormalized.xlsx")


### FOLD CHANGE for STASNet INPUT after Normalization 2 

N2_Control_means_GELwise <- WBN2_GelMean_Norm %>% 
  filter(Treatment=="control") %>% 
  group_by(x_status, GEL_Number) %>% 
  summarise_at(.vars = vars(ERKp:GAPDH),
               .funs = funs(mean(.))) %>% 
  rename_at(.vars = vars(ERKp:GAPDH),
            .funs = funs(paste0(.,"_mean")))

WBN2_GelMean_Norm_forFC <- WBN2_GelMean_Norm

N2_Control_means_GELwise$id <- paste0(N2_Control_means_GELwise$x_status,"_",N2_Control_means_GELwise$GEL_Number)
WBN2_GelMean_Norm_forFC$id <- paste0(WBN2_GelMean_Norm_forFC$x_status,"_",WBN2_GelMean_Norm_forFC$GEL_Number)

WBN2_GelMean_Norm_forFC <- WBN2_GelMean_Norm_forFC %>% 
  filter(Treatment != "Common") %>% 
  select(-c(GEL_Number,x_status)) #Removing the columns which are in both the data frames else they get a funny suffix.

WBN2_GelMean_Norm_FC_calc <- left_join(WBN2_GelMean_Norm_forFC,N2_Control_means_GELwise,by= "id")

WBN2_GelMean_Norm_FC_Statsnet <- WBN2_GelMean_Norm_FC_calc %>% 
  transmute(GEL_Number, # Transmute retains all the columns that are passed to it.  
            Treatment, #So u can list the columns that u want retained even when u dont use them to calculate new columns.
            x_status,
            Replicate,
            ERKp_FC = ERKp/ERKp_mean,
            STAT3p_FC = STAT3p/STAT3p_mean,
            SMAD2p_FC = SMAD2p/SMAD2p_mean,
            bCatenin_FC = bCatenin/bCatenin_mean,
            GAPDH_FC = GAPDH/GAPDH_mean)

openxlsx::write.xlsx(WBN2_GelMean_Norm_FC_Statsnet, "./OUTPUT/WBN2_GelMean_Norm_FC_Statsnet.xlsx")


############ FOR_PAPER ##############

WBData_MeanNormalized_FC <- WBN2_GelMean_Norm_FC_Statsnet %>% 
  select(x_status,Replicate,Treatment,GEL_Number,ERKp_FC,STAT3p_FC,SMAD2p_FC)

WBData_MeanNormalized_FC$Treatment <- gsub("control","DMSO",WBData_MeanNormalized_FC$Treatment )
WBData_MeanNormalized_FC$Replicate <- gsub("R3","Rep1",WBData_MeanNormalized_FC$Replicate )
WBData_MeanNormalized_FC$Replicate <- gsub("R4","Rep2",WBData_MeanNormalized_FC$Replicate )
WBData_MeanNormalized_FC$Replicate <- gsub("R5","Rep3",WBData_MeanNormalized_FC$Replicate )

WBData_MeanNormalized_FC <- WBData_MeanNormalized_FC %>% 
  filter(Treatment != "Common")

openxlsx::write.xlsx(WBData_MeanNormalized_FC, "./OUTPUT_PAPER/Suppl_Table_S3_3_WBData_MeanNormalized_FC.xlsx")

print("Script 2B : All steps executed")






