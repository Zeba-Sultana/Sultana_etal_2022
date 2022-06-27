#!/usr/bin/env Rscript

library(readxl) 
library(openxlsx)
#library(dplyr)
#library(reshape2) # for melt
#library(ggplot2)
library(tidyverse) 
library(tibble)  # for `rownames_to_column` and `column_to_rownames`

dir.create("OUTPUT")
dir.create("OUTPUT_PAPER")

## Functions :

# 1.) Function to read in gel data : read_in_gel_data
# Function inputs : 
# Gel_samples <- vector specifying names of the treatment in the seq that they were on the gel
# quant_files <- vector specifying the name of the xls files that have quantification information. 

# Function output : 
# Returns Dataframe with following 18 cols for all 6 analytes(ERKp,STAT3p,SMAD2p,b-Catenin,GAPDH,TPS)-> Raw signal values(6cols)+ TPS normalized values(6cols) and GAPDH normalized values(6cols)

# What function does :
# Creates "Gel", a df with well number and treatment names.
# Goes through the quantification files in a for loop(S1.xls,S2.xls,S3.xls,TPS.xls)
# Reads in the "signal"" value from the correct channel and creates a column in the df "Gel".(Total of 6 columns created)
# Drops ladder and empty well rows.
# Create new columns - G_Norm = AnalyteSignal/GAPDHsignal in the well and T_Norm=AnalyteSignal/TPSsignal(12 additional columns)
read_in_gel_data <- function(Gel_samples,quant_files) {
  
  Wells = c("Well_1","Well_2","Well_3","Well_4","Well_5","Well_6","Well_7","Well_8","Well_9","Well_10","Well_11","Well_12","Well_13","Well_14","Well_15","Well_16","Well_17","Well_18","Well_19","Well_20")
  Gel = data.frame(Gel_samples,row.names = Wells) #Makes a df with just one col. This col has treatment names in correct sequence. Rownames of df are "Well numbers 1-20"
  
  Gel_data_list = list()
  Gel_analytes_data_list = list()
  
  for (ii in 1:length(quant_files)){
    
    my_file_name <- strsplit(quant_files[ii],"\\.")[[1]][1] #need two backslashes to escape the dot being interpreted as special charachter. Removed xls from name of file. Now, my_file_name is R3_G1_S1(for example)
    
    # read in data from excel file R3_G1_S1(for example)
    Gel_data_list[[as.character(parse(text = my_file_name))]] = read_excel(file.path(quant_folder,quant_files[ii]))
    
    channel800_analytes <- c("ERKp","SMAD2p","bCatenin")
    
    step_number <- strsplit(my_file_name,"_")[[1]][3] #This yields : S1/S2/S3/TPS
    if(step_number == "S1") {
      analytes <- c("ERKp", "STAT3p")
      #TO_IMPROVE : the keyword "analytes" is very generic. Change it to "Step1_analytes"/"Step2_analytes" etc
      for (i in 1:length(analytes)){
        my_analyte_rows <- paste(names(Gel_data_list[1]),analytes[i],sep="_")
        if(analytes[i] %in% channel800_analytes ){
          Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 800)
          Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
        }else{
          Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 700)
          Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
        }  
      }
    }
    if(step_number == "S2") {
      analytes <- c("SMAD2p", "GAPDH")
      for (i in 1:length(analytes)){
        my_analyte_rows <- paste(names(Gel_data_list[1]),analytes[i],sep="_")
        if(analytes[i] %in% channel800_analytes ){
          Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 800)
          Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
        }else{
          Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 700)
          Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
        }  
      }
    }
    if(step_number == "S3") {
      analytes <- c("bCatenin")
      for (i in 1:length(analytes)){
        my_analyte_rows <- paste(names(Gel_data_list[1]),analytes[i],sep="_")
        Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 800)
        Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
      }
    }
    else{
      analytes <- c("TPS")
      for (i in 1:length(analytes)){
        my_analyte_rows <- paste(names(Gel_data_list[1]),analytes[i],sep="_")
        Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == "W")
        Gel[,as.character(parse(text = analytes[i]))] <- abs(Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal)
      }
    }
  }
  
  #Drop the empty cells and ladder
  row.names(Gel) <- NULL
  row.names(Gel) <- Gel_samples
  drop_cols <- c("Gel_samples")
  drop_wells <- c("Empty1","Empty2", "Ladder")
  Gel <- Gel[!(rownames(Gel) %in% drop_wells),!(colnames(Gel) %in% drop_cols)]
  
  all_analytes <- c("ERKp", "STAT3p", "SMAD2p","bCatenin", "GAPDH", "TPS")
  #In the for loop above the term 'analytes' has been used for step specific analytes. Eg for step 1 : analytes <- c("ERKp", "STAT3p")
  #Also, channel800_analytes <- c("ERKp","SMAD2p","bCatenin")
  #So, all_analytes used here is for all 6 analyte signals collected.
  
  #Normalization for unequal loading across wells
  for (analyte in all_analytes){
    GNormalized_col <- paste(analyte,"GNorm",sep="_")
    Gel[,as.character(parse(text = GNormalized_col))] <- (Gel[,as.character(parse(text = analyte))]/Gel$GAPDH)
    
    TNormalized_col <- paste(analyte,"TNorm",sep="_")
    Gel[,as.character(parse(text = TNormalized_col))] <- (Gel[,as.character(parse(text = analyte))]/Gel$TPS)
    
  }
  return(Gel)
  
}

# The function used for reading in data from Replicate5 gels is slighly different. (In this case it is needed to define Wells outside the generic funtion because at that stage during 
# quantification empty lanes and ladder lanes were not quantified resulting in different number of lanes per gel. Later(for R3 and R4) I quantified 20 lanes per gel and so included the 
# empty/ladder lanes as well for consistently having 20 lanes per gel.)
read_in_gel_data_Rep5 <- function(Gel_samples,quant_files,Wells) {
  
  Gel = data.frame(Gel_samples,row.names = Wells) #Makes a df with just one col. This col has treatment names in correct sequence. Rownames of df are "Well numbers 1-20"
  
  Gel_data_list = list()
  Gel_analytes_data_list = list()
  
  for (ii in 1:length(quant_files)){
    
    my_file_name <- strsplit(quant_files[ii],"\\.")[[1]][1] #need two backslashes to escape the dot being interpreted as special charachter. Removed xls from name of file. Now, my_file_name is R5_G1_S1(for example)
    
    # read in data from excel file R5_G1_S1(for example)
    Gel_data_list[[as.character(parse(text = my_file_name))]] = read_excel(file.path(quant_folder,quant_files[ii]))
    
    channel800_analytes <- c("ERKp","SMAD2p","bCatenin")
    
    step_number <- strsplit(my_file_name,"_")[[1]][3] #This yields : S1/S2/S3/TPS
    if(step_number == "S1") {
      analytes <- c("ERKp", "STAT3p")
      #TO_IMPROVE : the keyword "analytes" is very generic. Change it to "Step1_analytes"/"Step2_analytes" etc
      for (i in 1:length(analytes)){
        my_analyte_rows <- paste(names(Gel_data_list[1]),analytes[i],sep="_")
        if(analytes[i] %in% channel800_analytes ){
          Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 800)
          Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
        }else{
          Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 700)
          Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
        }  
      }
    }
    if(step_number == "S2") {
      analytes <- c("SMAD2p", "GAPDH")
      for (i in 1:length(analytes)){
        my_analyte_rows <- paste(names(Gel_data_list[1]),analytes[i],sep="_")
        if(analytes[i] %in% channel800_analytes ){
          Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 800)
          Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
        }else{
          Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 700)
          Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
        }  
      }
    }
    if(step_number == "S3") {
      analytes <- c("bCatenin")
      for (i in 1:length(analytes)){
        my_analyte_rows <- paste(names(Gel_data_list[1]),analytes[i],sep="_")
        Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == 800)
        Gel[,as.character(parse(text = analytes[i]))] <- Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal
      }
    }
    else{
      analytes <- c("TPS")
      for (i in 1:length(analytes)){
        my_analyte_rows <- paste(names(Gel_data_list[1]),analytes[i],sep="_")
        Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]] <- as.data.frame(Gel_data_list[[ii]]) %>% filter(Channel == "W")
        Gel[,as.character(parse(text = analytes[i]))] <- abs(Gel_analytes_data_list[[as.character(parse(text = my_analyte_rows))]]$Signal)
      }
    }
  }
  
  #Drop the empty cells and ladder
  row.names(Gel) <- NULL
  row.names(Gel) <- Gel_samples
  drop_cols <- c("Gel_samples")
  drop_wells <- c("Empty1","Empty2", "Ladder")
  Gel <- Gel[!(rownames(Gel) %in% drop_wells),!(colnames(Gel) %in% drop_cols)]
  
  all_analytes <- c("ERKp", "STAT3p", "SMAD2p","bCatenin", "GAPDH", "TPS")
  #In the for loop above the term 'analytes' has been used for step specific analytes. Eg for step 1 : analytes <- c("ERKp", "STAT3p")
  #Also, channel800_analytes <- c("ERKp","SMAD2p","bCatenin")
  #So, all_analytes used here is for all 6 analyte signals collected.
  
  #Normalization for unequal loading across wells
  for (analyte in all_analytes){
    GNormalized_col <- paste(analyte,"GNorm",sep="_")
    Gel[,as.character(parse(text = GNormalized_col))] <- (Gel[,as.character(parse(text = analyte))]/Gel$GAPDH)
    
    TNormalized_col <- paste(analyte,"TNorm",sep="_")
    Gel[,as.character(parse(text = TNormalized_col))] <- (Gel[,as.character(parse(text = analyte))]/Gel$TPS)
    
  }
  return(Gel)
  
}


# 2.) Function to add annotations like treatment name and replicate & gel number : Gel_Annotation
Gel_Annotation <- function(my_gel_data){
  

  my_file_name <- as.character(substitute(my_gel_data))
  
  
  gel_number <- strsplit(my_file_name,"_")[[1]][2] #This yields : G1/.../G7
  my_gel_data$GEL_Number <- rep(gel_number,nrow(my_gel_data))
  
  my_gel_data$Treatment <- row.names(my_gel_data)
  
  my_gel_data <- my_gel_data %>%
    rownames_to_column("temporary") %>% # Using mutate lost rownames of the df.# So, this is to preserve the rownames when mutate is used. Convert rownames to a column and then back
    mutate(x_status = ifelse(grepl("XX",my_gel_data$Treatment),paste0("XX"),paste0("XO")))%>%
    column_to_rownames("temporary")
  my_gel_data$x_status[which(my_gel_data$Treatment=="Common")] <- "Common"
  
  my_gel_data$Treatment <- gsub("XX_","",my_gel_data$Treatment)
  my_gel_data$Treatment <-  gsub("XO_","",my_gel_data$Treatment)
  
  replicate_number <- strsplit(my_file_name,"_")[[1]][1] #This yields : R3/R4/R5
  my_gel_data$Replicate <- rep(replicate_number,nrow(my_gel_data))
  
  return(my_gel_data)
  
}


########### ########### ########### ########### ########### 
########### Reading in WB data : Replicate 3 ########### 
########### ########### ########### ########### ########### 

# Input Quantification Data and Images 

quant_folder="../RAW_DATA/BioplexLysates_WB/QUANTIFIED/Replicate3/R3_Quantification/"
image_folder="../RAW_DATA/BioplexLysates_WB/QUANTIFIED/Replicate3/R3_Images"

## Gel1_Fgf4
R3_G1_samples = c("Ladder","Empty1","XX_Activin+NoLif","XX_Fgf4+Bmp4i","XX_Fgf4+Pi3ki","XX_Fgf4+Igfri","XX_Fgf4+Jaki","XX_Fgf4+NoLif","XX_Fgf4","XX_Cntrl","Common","XO_Cntrl","XO_Fgf4","XO_Fgf4+NoLif","XO_Fgf4+Jaki","XO_Fgf4+Igfri","XO_Fgf4+Pi3ki","XO_Fgf4+Bmp4i","XO_Activin+NoLif","Empty2")
R3_G1_quant_files = c("R3_G1_S1.xls","R3_G1_S2.xls","R3_G1_S3.xls","R3_G1_TPS.xls")
R3_G1_data <- read_in_gel_data(R3_G1_samples,R3_G1_quant_files)
R3_G1_data["XX_Fgf4+Pi3ki", grepl("SMAD2|bCatenin",colnames(R3_G1_data))] <- NA  # if there was difficulty with signal quantification due to background problem, the respective analyte reading was replaced with NA

## Gel2_Fgfri
R3_G2_samples = c("Ladder","XX_Fgfri+Bmp4i","XX_Fgfri+Jaki","XX_Fgfri+Gsk3bi","XX_Fgfri+Igfri","XX_Fgfri+Pi3ki","XX_Activin+Fgfri","XX_Fgf4+Fgfri","XX_Fgfri","XX_Cntrl","Common","XO_Cntrl","XO_Fgfri","XO_Fgf4+Fgfri","XO_Activin+Fgfri","XO_Fgfri+Pi3ki","XO_Fgfri+Igfri","XO_Fgfri+Gsk3bi","XO_Fgfri+Jaki","XO_Fgfri+Bmp4i")
R3_G2_quant_files = c("R3_G2_S1.xls","R3_G2_S2.xls","R3_G2_S3.xls","R3_G2_TPS.xls")
R3_G2_data <- read_in_gel_data(R3_G2_samples,R3_G2_quant_files)

## Gel3_Gsk3bi
R3_G3_samples = c("Ladder","XX_Bmp4i+Gsk3bi","XX_Igfri+Gsk3bi","XX_Pi3ki+Gsk3bi","XX_Jaki+Gsk3bi","XX_Fgf4+Gsk3bi","XX_Jaki+Bmp4i","XX_Gsk3bi","XX_Cntrl","Common","XO_Cntrl","XO_Gsk3bi","XO_Jaki+Bmp4i","XO_Fgf4+Gsk3bi","XO_Jaki+Gsk3bi","XO_Pi3ki+Gsk3bi","XO_Igfri+Gsk3bi","XO_Bmp4i+Gsk3bi","Empty1","Empty2")
R3_G3_quant_files = c("R3_G3_S1.xls","R3_G3_S2.xls","R3_G3_S3.xls","R3_G3_TPS.xls")
R3_G3_data <- read_in_gel_data(R3_G3_samples,R3_G3_quant_files)

## Gel4_Meki
R3_G4_samples = c("Ladder","XX_Meki+Bmp4i","XX_Meki+Jaki","XX_Meki+Gsk3bi","XX_Activin+Meki","XX_Fgf4+Meki","XX_Meki+Igfri","XX_Meki+Pi3ki","XX_Meki","XX_Cntrl","Common","XO_Cntrl","XO_Meki","XO_Meki+Pi3ki","XO_Meki+Igfri","XO_Fgf4+Meki","XO_Activin+Meki","XO_Meki+Gsk3bi","XO_Meki+Jaki","XO_Meki+Bmp4i")
R3_G4_quant_files = c("R3_G4_S1.xls","R3_G4_S2.xls","R3_G4_S3.xls","R3_G4_TPS.xls")
R3_G4_data <- read_in_gel_data(R3_G4_samples,R3_G4_quant_files)

## Gel5_Activin
R3_G5_samples = c("Ladder","Empty1","XX_Activin+Bmp4i","XX_Activin+Jaki","XX_Activin+Gsk3bi","XX_Activin+Igfri","XX_Activin+Pi3ki","XX_Activin+Fgf4","XX_Activin","XX_Cntrl","Common","XO_Cntrl","XO_Activin","XO_Activin+Fgf4","XO_Activin+Pi3ki","XO_Activin+Igfri","XO_Activin+Gsk3bi","XO_Activin+Jaki","XO_Activin+Bmp4i","Empty2")
R3_G5_quant_files = c("R3_G5_S1.xls","R3_G5_S2.xls","R3_G5_S3.xls","R3_G5_TPS.xls")
R3_G5_data <- read_in_gel_data(R3_G5_samples,R3_G5_quant_files)

## Gel6_misc
R3_G6_samples = c("Ladder","XX_Bmp4i","XX_Pi3ki+Bmp4i","XX_Igfri+Bmp4i","XX_Igfri","XX_Pi3ki","XX_Pi3ki+Jaki","XX_Igfri+Jaki","XX_Jaki","XX_Cntrl","Common","XO_Cntrl","XO_Jaki","XO_Igfri+Jaki","XO_Pi3ki+Jaki","XO_Pi3ki","XO_Igfri","XO_Igfri+Bmp4i","XO_Pi3ki+Bmp4i","XO_Bmp4i")
R3_G6_quant_files = c("R3_G6_S1.xls","R3_G6_S2.xls","R3_G6_S3.xls","R3_G6_TPS.xls")
R3_G6_data <- read_in_gel_data(R3_G6_samples,R3_G6_quant_files)

## Gel7_NoLif
R3_G7_samples = c("Ladder","XX_Bmp4i+NoLif","XX_Jaki+NoLif","XX_Gsk3bi+NoLif","XX_Igfri+NoLif","XX_Pi3ki+NoLif","XX_Fgfri+NoLif","XX_Meki+NoLif","XX_NoLif","XX_Cntrl","Common","XO_Cntrl","XO_NoLif","XO_Meki+NoLif","XO_Fgfri+NoLif","XO_Pi3ki+NoLif","XO_Igfri+NoLif","XO_Gsk3bi+NoLif","XO_Jaki+NoLif","XO_Bmp4i+NoLif")
R3_G7_quant_files = c("R3_G7_S1.xls","R3_G7_S2.xls","R3_G7_S3.xls","R3_G7_TPS.xls")
R3_G7_data <- read_in_gel_data(R3_G7_samples,R3_G7_quant_files)

#Saving DF1 : 
#df of raw signal values and their loading normalized values(_TNorm and _GNorm) along with annotation
R3_G1_data_anno <- Gel_Annotation(R3_G1_data)
R3_G2_data_anno <- Gel_Annotation(R3_G2_data)
R3_G3_data_anno <- Gel_Annotation(R3_G3_data)
R3_G4_data_anno <- Gel_Annotation(R3_G4_data)
R3_G5_data_anno <- Gel_Annotation(R3_G5_data)
R3_G6_data_anno <- Gel_Annotation(R3_G6_data)
R3_G7_data_anno <- Gel_Annotation(R3_G7_data)

Replicate3_Alldata <- rbind(R3_G1_data_anno,
                            R3_G2_data_anno,
                            R3_G3_data_anno,
                            R3_G4_data_anno,
                            R3_G5_data_anno,
                            R3_G6_data_anno,
                            R3_G7_data_anno)


#WriteXLS::WriteXLS(Replicate3_Alldata, "./OUTPUT/Replicate3_Alldata.xls")
openxlsx::write.xlsx(Replicate3_Alldata, "./OUTPUT/Replicate3_Alldata.xlsx")

WBData_Rep1 <- Replicate3_Alldata %>% 
  select(!grep("GNorm",colnames(Replicate3_Alldata))) %>% 
  select(x_status,Replicate,Treatment,GEL_Number,ERKp,STAT3p,SMAD2p,GAPDH,TPS,
         ERKp_TNorm,STAT3p_TNorm,SMAD2p_TNorm)

 WBData_Rep1$Replicate <- gsub("R3","Rep1",WBData_Rep1$Replicate )
 WBData_Rep1$Treatment <- gsub("Cntrl","DMSO",WBData_Rep1$Treatment )
 

########### ########### ########### ########### ########### 
########### Reading in WB data : Replicate 4 ########### 
########### ########### ########### ########### ###########

quant_folder="../RAW_DATA/BioplexLysates_WB/QUANTIFIED/Replicate4/R4_Quantification/"
image_folder="../RAW_DATA/BioplexLysates_WB/QUANTIFIED/Replicate4/R4_Images/"


## Gel1_Fgf4
R4_G1_samples = c("Ladder","XX_Activin+NoLif","XX_Fgf4+Bmp4i","XX_Fgf4+Pi3ki","XX_Fgf4+Igfri","Empty1","XX_Fgf4+Jaki","XX_Fgf4+NoLif","XX_Fgf4","XX_Cntrl","Common","XO_Cntrl","XO_Fgf4","XO_Fgf4+NoLif","XO_Fgf4+Jaki","XO_Fgf4+Igfri","XO_Fgf4+Pi3ki","XO_Fgf4+Bmp4i","XO_Activin+NoLif","Empty2")
R4_G1_quant_files = c("R4_G1_S1.xls","R4_G1_S2.xls","R4_G1_S3.xls","R4_G1_TPS.xls")
R4_G1_data <- read_in_gel_data(R4_G1_samples,R4_G1_quant_files)

## Gel2_Fgfri
R4_G2_samples = c("Ladder","XX_Fgfri+Bmp4i","XX_Fgfri+Jaki","XX_Fgfri+Gsk3bi","XX_Fgfri+Igfri","XX_Fgfri+Pi3ki","XX_Activin+Fgfri","XX_Fgf4+Fgfri","XX_Fgfri","XX_Cntrl","Common","XO_Cntrl","XO_Fgfri","XO_Fgf4+Fgfri","XO_Activin+Fgfri","XO_Fgfri+Pi3ki","XO_Fgfri+Igfri","XO_Fgfri+Gsk3bi","XO_Fgfri+Jaki","XO_Fgfri+Bmp4i")
R4_G2_quant_files = c("R4_G2_S1.xls","R4_G2_S2.xls","R4_G2_S3.xls","R4_G2_TPS.xls")
R4_G2_data <- read_in_gel_data(R4_G2_samples,R4_G2_quant_files)

## Gel3_Gsk3bi
R4_G3_samples = c("Ladder","Empty1","XX_Bmp4i+Gsk3bi","XX_Igfri+Gsk3bi","XX_Pi3ki+Gsk3bi","XX_Jaki+Gsk3bi","XX_Fgf4+Gsk3bi","XX_Jaki+Bmp4i","XX_Gsk3bi","XX_Cntrl","Common","XO_Cntrl","XO_Gsk3bi","XO_Jaki+Bmp4i","XO_Fgf4+Gsk3bi","XO_Jaki+Gsk3bi","XO_Pi3ki+Gsk3bi","XO_Igfri+Gsk3bi","XO_Bmp4i+Gsk3bi","Empty2")
R4_G3_quant_files = c("R4_G3_S1.xls","R4_G3_S2.xls","R4_G3_S3.xls","R4_G3_TPS.xls")
R4_G3_data <- read_in_gel_data(R4_G3_samples,R4_G3_quant_files)

## Gel4_Meki
R4_G4_samples = c("Ladder","XX_Meki+Bmp4i","XX_Meki+Jaki","XX_Meki+Gsk3bi","XX_Activin+Meki","XX_Fgf4+Meki","XX_Meki+Igfri","XX_Meki+Pi3ki","XX_Meki","XX_Cntrl","Common","XO_Cntrl","XO_Meki","XO_Meki+Pi3ki","XO_Meki+Igfri","XO_Fgf4+Meki","XO_Activin+Meki","XO_Meki+Gsk3bi","XO_Meki+Jaki","XO_Meki+Bmp4i")
R4_G4_quant_files = c("R4_G4_S1.xls","R4_G4_S2.xls","R4_G4_S3.xls","R4_G4_TPS.xls")
R4_G4_data <- read_in_gel_data(R4_G4_samples,R4_G4_quant_files)

## Gel5_Activin
R4_G5_samples = c("Ladder","Empty1","XX_Activin+Bmp4i","XX_Activin+Jaki","XX_Activin+Gsk3bi","XX_Activin+Igfri","XX_Activin+Pi3ki","XX_Activin+Fgf4","XX_Activin","XX_Cntrl","Common","XO_Cntrl","XO_Activin","XO_Activin+Fgf4","XO_Activin+Pi3ki","XO_Activin+Igfri","XO_Activin+Gsk3bi","XO_Activin+Jaki","XO_Activin+Bmp4i","Empty2")
R4_G5_quant_files = c("R4_G5_S1.xls","R4_G5_S2.xls","R4_G5_S3.xls","R4_G5_TPS.xls")
R4_G5_data <- read_in_gel_data(R4_G5_samples,R4_G5_quant_files)

## Gel6_misc
R4_G6_samples = c("Ladder","XX_Bmp4i","XX_Pi3ki+Bmp4i","XX_Igfri+Bmp4i","XX_Igfri","XX_Pi3ki","XX_Pi3ki+Jaki","XX_Igfri+Jaki","XX_Jaki","XX_Cntrl","Common","XO_Cntrl","XO_Jaki","XO_Igfri+Jaki","XO_Pi3ki+Jaki","XO_Pi3ki","XO_Igfri","XO_Igfri+Bmp4i","XO_Pi3ki+Bmp4i","XO_Bmp4i")
R4_G6_quant_files = c("R4_G6_S1.xls","R4_G6_S2.xls","R4_G6_S3.xls","R4_G6_TPS.xls")
R4_G6_data <- read_in_gel_data(R4_G6_samples,R4_G6_quant_files)

## Gel7_NoLif
R4_G7_samples = c("Ladder","XX_Bmp4i+NoLif","XX_Jaki+NoLif","XX_Gsk3bi+NoLif","XX_Igfri+NoLif","XX_Pi3ki+NoLif","XX_Fgfri+NoLif","XX_Meki+NoLif","XX_NoLif","XX_Cntrl","Common","XO_Cntrl","XO_NoLif","XO_Meki+NoLif","XO_Fgfri+NoLif","XO_Pi3ki+NoLif","XO_Igfri+NoLif","XO_Gsk3bi+NoLif","XO_Jaki+NoLif","XO_Bmp4i+NoLif")
R4_G7_quant_files = c("R4_G7_S1.xls","R4_G7_S2.xls","R4_G7_S3.xls","R4_G7_TPS.xls")
R4_G7_data <- read_in_gel_data(R4_G7_samples,R4_G7_quant_files)

#Saving DF1 : 
#df of raw signal values and their loading normalized values(_TNorm and _GNorm) along with annotation
R4_G1_data_anno <- Gel_Annotation(R4_G1_data)
R4_G2_data_anno <- Gel_Annotation(R4_G2_data)
R4_G3_data_anno <- Gel_Annotation(R4_G3_data)
R4_G4_data_anno <- Gel_Annotation(R4_G4_data)
R4_G5_data_anno <- Gel_Annotation(R4_G5_data)
R4_G6_data_anno <- Gel_Annotation(R4_G6_data)
R4_G7_data_anno <- Gel_Annotation(R4_G7_data)

Replicate4_Alldata <- rbind(R4_G1_data_anno,
                            R4_G2_data_anno,
                            R4_G3_data_anno,
                            R4_G4_data_anno,
                            R4_G5_data_anno,
                            R4_G6_data_anno,
                            R4_G7_data_anno)

#WriteXLS::WriteXLS(Replicate4_Alldata, "./OUTPUT/Replicate4_Alldata.xls")
openxlsx::write.xlsx(Replicate4_Alldata, "./OUTPUT/Replicate4_Alldata.xlsx")

WBData_Rep2 <- Replicate4_Alldata %>% 
  select(!grep("GNorm",colnames(Replicate4_Alldata))) %>% 
  select(x_status,Replicate,Treatment,GEL_Number,ERKp,STAT3p,SMAD2p,GAPDH,TPS,
         ERKp_TNorm,STAT3p_TNorm,SMAD2p_TNorm)

WBData_Rep2$Replicate <- gsub("R4","Rep2",WBData_Rep2$Replicate )
WBData_Rep2$Treatment <- gsub("Cntrl","DMSO",WBData_Rep2$Treatment )


########### ########### ########### ########### ########### 
########### Reading in WB data : Replicate 5 ########### 
########### ########### ########### ########### ###########

quant_folder="../RAW_DATA/BioplexLysates_WB/QUANTIFIED/Replicate5/R5_Quantification/"
image_folder="../RAW_DATA/BioplexLysates_WB/QUANTIFIED/Replicate5/R5_Images/"


## Gel1_Fgf4
R5_G1_samples = c("XX_Activin+NoLif","XX_Fgf4+Bmp4i","XX_Fgf4+Pi3ki","XX_Fgf4+Igfri","XX_Fgf4+Jaki","XX_Fgf4+NoLif","XX_Fgf4","XX_Cntrl","Common1","XO_Cntrl","XO_Fgf4","XO_Fgf4+NoLif","XO_Fgf4+Jaki","XO_Fgf4+Igfri","XO_Fgf4+Pi3ki","XO_Fgf4+Bmp4i","XO_Activin+NoLif","Common2")
R5_G1_quant_files = c("R5_G1_S1.xls","R5_G1_S2.xls","R5_G1_S3.xls","R5_G1_TPS.xls")
Wells = c("Well_1","Well_2","Well_3","Well_4","Well_5","Well_6","Well_7","Well_8","Well_9","Well_10","Well_11","Well_12","Well_13","Well_14","Well_15","Well_16","Well_17","Well_18")

R5_G1_data <- read_in_gel_data_Rep5(R5_G1_samples,R5_G1_quant_files,Wells)
R5_G1_data["XX_Fgf4+Jaki",] <- rep(NA,ncol(R5_G1_data))


## Gel2_Fgfri
R5_G2_samples = c("XX_Fgfri+Bmp4i","XX_Fgfri+Jaki","XX_Fgfri+Gsk3bi","XX_Fgfri+Igfri","XX_Fgfri+Pi3ki","XX_Activin+Fgfri","XX_Fgf4+Fgfri","XX_Fgfri","XX_Cntrl","Common","XO_Cntrl","XO_Fgfri","XO_Fgf4+Fgfri","XO_Activin+Fgfri","XO_Fgfri+Pi3ki","XO_Fgfri+Igfri","XO_Fgfri+Gsk3bi","XO_Fgfri+Jaki","XO_Fgfri+Bmp4i")
R5_G2_quant_files = c("R5_G2_S1.xls","R5_G2_S2.xls","R5_G2_S3.xls","R5_G2_TPS.xls")
Wells = c("Well_1","Well_2","Well_3","Well_4","Well_5","Well_6","Well_7","Well_8","Well_9","Well_10","Well_11","Well_12","Well_13","Well_14","Well_15","Well_16","Well_17","Well_18","Well_19")

R5_G2_data <- read_in_gel_data_Rep5(R5_G2_samples,R5_G2_quant_files,Wells)

## Gel3_Gsk3bi
R5_G3_samples = c("XX_Bmp4i+Gsk3bi","XX_Igfri+Gsk3bi","XX_Pi3ki+Gsk3bi","XX_Jaki+Gsk3bi","XX_Fgf4+Gsk3bi","XX_Jaki+Bmp4i","XX_Gsk3bi","XX_Cntrl","Common","XO_Cntrl","XO_Gsk3bi","XO_Jaki+Bmp4i","XO_Fgf4+Gsk3bi","XO_Jaki+Gsk3bi","XO_Pi3ki+Gsk3bi","XO_Igfri+Gsk3bi","XO_Bmp4i+Gsk3bi")
R5_G3_quant_files = c("R5_G3_S1.xls","R5_G3_S2.xls","R5_G3_S3.xls","R5_G3_TPS.xls")
Wells = c("Well_1","Well_2","Well_3","Well_4","Well_5","Well_6","Well_7","Well_8","Well_9","Well_10","Well_11","Well_12","Well_13","Well_14","Well_15","Well_16","Well_17")

R5_G3_data <- read_in_gel_data_Rep5(R5_G3_samples,R5_G3_quant_files,Wells)
R5_G3_data["XO_Gsk3bi",] <- rep(NA,ncol(R5_G3_data))
R5_G3_data["XO_Fgf4+Gsk3bi",] <- rep(NA,ncol(R5_G3_data))
R5_G3_data["XO_Pi3ki+Gsk3bi",] <- rep(NA,ncol(R5_G3_data))

## Gel4_Meki
R5_G4_samples = c("XX_Meki+Bmp4i","XX_Meki+Jaki","XX_Meki+Gsk3bi","XX_Activin+Meki","XX_Fgf4+Meki","XX_Meki+Igfri","XX_Meki+Pi3ki","XX_Meki","XX_Cntrl","Common","XO_Cntrl","XO_Meki","XO_Meki+Pi3ki","XO_Meki+Igfri","XO_Fgf4+Meki","XO_Activin+Meki","XO_Meki+Gsk3bi","XO_Meki+Jaki","XO_Meki+Bmp4i")
R5_G4_quant_files = c("R5_G4_S1.xls","R5_G4_S2.xls","R5_G4_S3.xls","R5_G4_TPS.xls")
Wells = c("Well_1","Well_2","Well_3","Well_4","Well_5","Well_6","Well_7","Well_8","Well_9","Well_10","Well_11","Well_12","Well_13","Well_14","Well_15","Well_16","Well_17","Well_18","Well_19")

R5_G4_data <- read_in_gel_data_Rep5(R5_G4_samples,R5_G4_quant_files,Wells)
R5_G4_data["XO_Meki+Igfri",] <- rep(NA,ncol(R5_G4_data))
R5_G4_data["XO_Fgf4+Meki",] <- rep(NA,ncol(R5_G4_data))
R5_G4_data["XO_Meki+Jaki",] <- rep(NA,ncol(R5_G4_data))

## Gel5_Activin
R5_G5_samples = c("Ladder","XX_Activin+Bmp4i","XX_Activin+Jaki","XX_Activin+Gsk3bi","Empty1","XX_Activin+Igfri","XX_Activin+Pi3ki","XX_Activin+Fgf4","XX_Activin","XX_Cntrl","Common","XO_Cntrl","XO_Activin","XO_Activin+Fgf4","XO_Activin+Pi3ki","XO_Activin+Igfri","XO_Activin+Gsk3bi","XO_Activin+Jaki","XO_Activin+Bmp4i","Empty2")
R5_G5_quant_files = c("R5_G5_S1.xls","R5_G5_S2.xls","R5_G5_S3.xls","R5_G5_TPS.xls")
Wells = c("Well_1","Well_2","Well_3","Well_4","Well_5","Well_6","Well_7","Well_8","Well_9","Well_10","Well_11","Well_12","Well_13","Well_14","Well_15","Well_16","Well_17","Well_18","Well_19","Well_20")

R5_G5_data <- read_in_gel_data_Rep5(R5_G5_samples,R5_G5_quant_files,Wells)
R5_G5_data["XO_Activin+Pi3ki",] <- rep(NA,ncol(R5_G5_data))

## Gel6_misc
R5_G6_samples = c("Ladder","XX_Bmp4i","XX_Pi3ki+Bmp4i","XX_Igfri+Bmp4i","XX_Igfri","XX_Pi3ki","XX_Pi3ki+Jaki","XX_Igfri+Jaki","XX_Jaki","XX_Cntrl","Common","XO_Cntrl","XO_Jaki","XO_Igfri+Jaki","XO_Pi3ki+Jaki","XO_Pi3ki","XO_Igfri","XO_Igfri+Bmp4i","XO_Pi3ki+Bmp4i","XO_Bmp4i")
R5_G6_quant_files = c("R5_G6_S1.xls","R5_G6_S2.xls","R5_G6_S3.xls","R5_G6_TPS.xls")
Wells = c("Well_1","Well_2","Well_3","Well_4","Well_5","Well_6","Well_7","Well_8","Well_9","Well_10","Well_11","Well_12","Well_13","Well_14","Well_15","Well_16","Well_17","Well_18","Well_19","Well_20")

R5_G6_data <- read_in_gel_data_Rep5(R5_G6_samples,R5_G6_quant_files,Wells)
R5_G6_data["XO_Jaki",] <- rep(NA,ncol(R5_G6_data))
R5_G6_data["XO_Igfri+Bmp4i",] <- rep(NA,ncol(R5_G6_data))

## Gel7_NoLif
R5_G7_samples = c("Ladder","XX_Bmp4i+NoLif","XX_Jaki+NoLif","XX_Gsk3bi+NoLif","XX_Igfri+NoLif","XX_Pi3ki+NoLif","XX_Fgfri+NoLif","XX_Meki+NoLif","XX_NoLif","XX_Cntrl","Common","XO_Cntrl","XO_NoLif","XO_Meki+NoLif","XO_Fgfri+NoLif","XO_Pi3ki+NoLif","XO_Igfri+NoLif","XO_Gsk3bi+NoLif","XO_Jaki+NoLif","XO_Bmp4i+NoLif")
R5_G7_quant_files = c("R5_G7_S1.xls","R5_G7_S2.xls","R5_G7_S3.xls","R5_G7_TPS.xls")
Wells = c("Well_1","Well_2","Well_3","Well_4","Well_5","Well_6","Well_7","Well_8","Well_9","Well_10","Well_11","Well_12","Well_13","Well_14","Well_15","Well_16","Well_17","Well_18","Well_19","Well_20")

R5_G7_data <- read_in_gel_data_Rep5(R5_G7_samples,R5_G7_quant_files,Wells)
R5_G7_data["XO_NoLif",] <- rep(NA,ncol(R5_G7_data))


R5_G1_data_anno <- Gel_Annotation(R5_G1_data)
R5_G2_data_anno <- Gel_Annotation(R5_G2_data)
R5_G3_data_anno <- Gel_Annotation(R5_G3_data)
R5_G4_data_anno <- Gel_Annotation(R5_G4_data)
R5_G5_data_anno <- Gel_Annotation(R5_G5_data)
R5_G6_data_anno <- Gel_Annotation(R5_G6_data)
R5_G7_data_anno <- Gel_Annotation(R5_G7_data)

Replicate5_Alldata_values <- rbind(R5_G1_data_anno,
                                   R5_G2_data_anno,
                                   R5_G3_data_anno,
                                   R5_G4_data_anno,
                                   R5_G5_data_anno,
                                   R5_G6_data_anno,
                                   R5_G7_data_anno)

######Calculations for XO/Common ratio

cntrl_Common <- Replicate5_Alldata_values %>% 
  filter(grepl("Cntrl|Common",Treatment)) %>% 
  filter(!grepl("Common2",Treatment))


cntrl_Common$x_status[cntrl_Common$Treatment == "Common1"] <- "Common"
cntrl_Common$Treatment <- gsub("Common1","Common",cntrl_Common$Treatment)


XX_by_Common <- cntrl_Common %>% 
  filter(x_status!="XO") %>% 
  group_by(GEL_Number) %>% 
  summarise_if(is.numeric, funs(.[x_status=="XX"]/.[x_status=="Common"])) %>% 
  #select(grep("Norm",colnames(cntrl_Common)))
  select(c(GEL_Number,grep("Norm",colnames(cntrl_Common))))


XO_by_Common <- cntrl_Common %>% 
  filter(x_status!="XX") %>% 
  group_by(GEL_Number) %>% 
  summarise_if(is.numeric, funs(.[x_status=="XO"]/.[x_status=="Common"])) %>% 
  #select(grep("Norm",colnames(cntrl_Common)))
  select(c(GEL_Number,grep("Norm",colnames(cntrl_Common))))

ratio_XO_Common_df_G3removed <- XO_by_Common %>% 
  filter(!GEL_Number == "G3") %>% 
  select(-c(GEL_Number))

mean_XO_Common_df_G3removed <- colMeans(ratio_XO_Common_df_G3removed)

R5_G3_data_redone <- R5_G3_data
all_analytes_Norm <- colnames(Replicate5_Alldata_values)[grepl("Norm",colnames(Replicate5_Alldata_values))]

for(analyte_Norm in all_analytes_Norm){
  R5_G3_data_redone["XO_Cntrl",analyte_Norm] <-mean_XO_Common_df_G3removed[analyte_Norm]*R5_G3_data["Common",analyte_Norm]
}



###### Adding the Recalculated XO control fr G3 Replicate 5 ######

R5_G1_data_anno <- Gel_Annotation(R5_G1_data)
R5_G2_data_anno <- Gel_Annotation(R5_G2_data)
R5_G3_data_anno <- Gel_Annotation(R5_G3_data_redone)
R5_G4_data_anno <- Gel_Annotation(R5_G4_data)
R5_G5_data_anno <- Gel_Annotation(R5_G5_data)
R5_G6_data_anno <- Gel_Annotation(R5_G6_data)
R5_G7_data_anno <- Gel_Annotation(R5_G7_data)

Replicate5_Alldata <- rbind(R5_G1_data_anno,
                            R5_G2_data_anno,
                            R5_G3_data_anno,
                            R5_G4_data_anno,
                            R5_G5_data_anno,
                            R5_G6_data_anno,
                            R5_G7_data_anno)


#WriteXLS::WriteXLS(Replicate5_Alldata, "./OUTPUT/Replicate5_Alldata.xls")
openxlsx::write.xlsx(Replicate5_Alldata, "./OUTPUT/Replicate5_Alldata.xlsx")

WBData_Rep3 <- Replicate5_Alldata %>% 
  select(!grep("GNorm",colnames(Replicate5_Alldata))) %>% 
  select(x_status,Replicate,Treatment,GEL_Number,ERKp,STAT3p,SMAD2p,GAPDH,TPS,
         ERKp_TNorm,STAT3p_TNorm,SMAD2p_TNorm)

WBData_Rep3$Replicate <- gsub("R5","Rep3",WBData_Rep3$Replicate )
WBData_Rep3$Treatment <- gsub("Cntrl","DMSO",WBData_Rep3$Treatment )


WB_RawData <- rbind(WBData_Rep1,
                    WBData_Rep2,
                    WBData_Rep3)

WB_RawData <- WB_RawData %>% 
  filter(!grepl("Common", WB_RawData$Treatment))


openxlsx::write.xlsx(WB_RawData, "./OUTPUT_PAPER/Suppl_Table_S3_1_WB_RawData.xlsx")

print("Script 1B : All steps executed")



