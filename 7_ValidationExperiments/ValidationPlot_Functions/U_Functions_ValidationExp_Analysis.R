
ReadInWBData <- function(Gel_samples, Geldb_data, T_Doses, GelNum = "GelNum", Analyte = "analyte") {
  ######################  IMP NOTE  ###################
  ## CROSS CHECK THE FUNCTION BEFORE USING ELSEWHERE ##
  ## Has some steps that are specific for this setup ##
  
  # unhash For Testing #######
  # Gel_samples <- Gel11_samples
  # Geldb_data <- Gel11_Smad2DB
  # T_Doses <- Activin_Doses
  # GelNum = "GelNum"
  # Analyte = "analyte"
  ##########################
  Geldb_data$Sample <- rep(Gel_samples,2)
  
  Geldb_phosP <- Geldb_data %>% filter(Channel == 800) %>% select(Signal,Sample)
  colnames(Geldb_phosP) <- gsub("Signal","phosP",colnames(Geldb_phosP))
  
  Geldb_totalP <- Geldb_data %>% filter(Channel == 700) %>% select(Signal,Sample)
  colnames(Geldb_totalP) <- gsub("Signal","totalP",colnames(Geldb_totalP))
  
  WBData_updated <- merge(Geldb_phosP,Geldb_totalP, by=c("Sample"), sort = FALSE) #Sort false keeps the original seq of rows.
  
  WBData_updated <- WBData_updated %>% 
    separate(Sample, c("Cell_line","Replicate","Treatment","Exp"), "_") # Needs dplyr and tidyr
  
  # WBData_updated <- WBData_updated %>% 
  #   filter(Cell_line %in% c("XX","XO"))
  
  WBData_updated$Gel <- rep(GelNum,nrow(WBData_updated))
  WBData_updated$Analyte <- rep(Analyte,nrow(WBData_updated))
  WBData_updated <- WBData_updated %>% 
    mutate(Well_no=row_number())
  
  WBData_updated$Treatment <- gsub("cntrl", "0", WBData_updated$Treatment)
  WBData_updated$Treatment <- as.numeric(WBData_updated$Treatment)
  WBData_updated$log2Treatment <- log2(WBData_updated$Treatment+1)
  WBData_updated$Treatment_Fctr <- factor(WBData_updated$Treatment,levels = T_Doses )
  
  return(WBData_updated)
  
}

ReadInWBData_TPS <- function(Gel_samples, Geldb_data, TPS_data, T_Doses, GelNum = "GelNum", Analyte = "analyte") {
  ######################  IMP NOTE  ###################
  ## CROSS CHECK THE FUNCTION BEFORE USING ELSEWHERE ##
  ## Has some steps that are specific for this setup ##
  
  # unhash For Testing #######
  # Gel_samples <- Gel19_samples
  # Geldb_data <- Gel19_P70s6kDB
  # TPS_data <- Gel19_TPS
  # T_Doses <- Gsk3inh_CT_Doses
  # GelNum = "GelNum"
  # Analyte = "analyte"
  ##########################
  
  Geldb_phosP <- Geldb_data %>% filter(Channel == 800) %>% select(Signal,Name)
  colnames(Geldb_phosP) <- gsub("Signal","phosP",colnames(Geldb_phosP))
  
  Geldb_phosP$Sample <- Gel_samples
  Geldb_phosP <-  Geldb_phosP %>% 
    select(-c(Name))
  
  # Geldb_totalP <- Geldb_data %>% filter(Channel == 700) %>% select(Signal,Sample)
  # colnames(Geldb_totalP) <- gsub("Signal","totalP",colnames(Geldb_totalP))
  
  TPS_data <- TPS_data[-(1:3),]
  colnames(TPS_data) <- TPS_data[1,]
  TPS_data <- TPS_data[-1,]
  TPS_data$Sample <- Gel_samples
  TPS_data <- TPS_data %>% 
    select(Signal,Sample)
  colnames(TPS_data) <- gsub("Signal","TPS",colnames(TPS_data))
  
  WBData_updated <- merge(Geldb_phosP,TPS_data, by=c("Sample"), sort = FALSE) #Sort false keeps the original seq of rows.
  #WBData_updated <- merge(Geldb_phosP,Geldb_totalP, by=c("Sample"), sort = FALSE) #Sort false keeps the original seq of rows.
  
  WBData_updated <- WBData_updated %>% 
    separate(Sample, c("Cell_line","Replicate","Treatment","Exp"), "_") # Needs dplyr and tidyr
  
  # WBData_updated <- WBData_updated %>% 
  #   filter(Cell_line %in% c("XX","XO"))
  
  WBData_updated$Gel <- rep(GelNum,nrow(WBData_updated))
  WBData_updated$Analyte <- rep(Analyte,nrow(WBData_updated))
  WBData_updated <- WBData_updated %>% 
    mutate(Well_no=row_number())
  
  WBData_updated$Treatment <- gsub("cntrl", "0", WBData_updated$Treatment)
  WBData_updated$Treatment <- as.numeric(WBData_updated$Treatment)
  WBData_updated$log2Treatment <- log2(WBData_updated$Treatment+1)
  WBData_updated$Treatment_Fctr <- factor(WBData_updated$Treatment,levels = T_Doses )
  
  colnames(WBData_updated) <- gsub("TPS","totalP",colnames(WBData_updated))
  
  return(WBData_updated)
  
}

#Objective : To normalize over the common smaple across gels OR the mean of common samples(if more than one common sample are present)
#Arguments to the function : db of WB data readin, Specifications to identify the common sample, ComSample_CellLine is a vector of cell line names that the common samples can have for the current gel.
Norm_o_ComSample <- function(WBData_updated, ComSample_Exp = "Exp_name", ComSample_Rep = "Replicate_number", ComSample_T = 0, ComSample_CellLine = c("XO","XX","nn") ) {
  ######################  IMP NOTE  ###################
  ## CROSS CHECK THE FUNCTION BEFORE USING ELSEWHERE ##
  ## Has some steps that are specific for this setup ##
  
  #### unhash For Testing #######
  # WBData_updated <- Gel23_AktDB_labeled
  # ComSample_Exp <- "CTDR"
  # ComSample_Rep = "R1"
  # ComSample_T = 0
  # ComSample_CellLine = c("XO")
  
  ##########################
  
  WBData_updated$phosP = as.numeric(WBData_updated$phosP)
  WBData_updated$totalP = as.numeric(WBData_updated$totalP)
  
  WBData_updated_NormOCom <- WBData_updated %>% 
    mutate(Mean_Com_phosP = mean(phosP[Exp == ComSample_Exp & Replicate == ComSample_Rep & Cell_line %in% ComSample_CellLine & ((Treatment == ComSample_T) | is.na(Treatment))])) %>% 
    mutate(Mean_Com_totalP = mean(totalP[Exp == ComSample_Exp & Replicate == ComSample_Rep & Cell_line %in% ComSample_CellLine & ((Treatment == ComSample_T) | is.na(Treatment)) ])) %>% 
    mutate(phosP_Norm_oCom = phosP/Mean_Com_phosP) %>% 
    mutate(totalP_Norm_oCom = totalP/Mean_Com_totalP) 
  
  # mutate(phosP_by_total = phosP/totalP) 
  return(WBData_updated_NormOCom)
  
}

Norm_o_MeanRep <- function(WBData_NormOCom_3Replicates, Norm_o_Com = 0){
  # unhash For Testing #######
  # WBData_NormOCom <- PDTC_Mek_TPS
  # Norm_o_Com <- Norm_o_Com
  ##########################
  
  WBData_NormOCom_3Replicates$phosP = as.numeric(WBData_NormOCom_3Replicates$phosP)
  WBData_NormOCom_3Replicates$totalP = as.numeric(WBData_NormOCom_3Replicates$totalP)
  
  if(Norm_o_Com == 1){
    WBData_NormORepMean <- WBData_NormOCom_3Replicates %>% 
      mutate(phosP_by_total = phosP_Norm_oCom/totalP_Norm_oCom) %>% #each phosP divided by totalP [ Both phosP and totalP values for each sample is already norm wrt common sample(or their mean, if more than one common samples)] 
      group_by(Replicate) %>% # remaining calculations are per replicate
      mutate(phosP_Mean_Rep = mean(phosP_Norm_oCom))%>%    # For use below, calculate mean of the phosP norm wrt common sample(per replicate) 
      mutate(totalP_Mean_Rep = mean(totalP_Norm_oCom))%>%  # For use below, calculate mean of the phosP norm wrt common sample(per replicate) 
      mutate(phosP_Norm_oMeanRep = phosP_Norm_oCom/phosP_Mean_Rep) %>% # Divide phosP_Norm_oCom by its mean for each replicate. This is to normalize over the mean of all samples in each replicate.
      mutate(totalP_Norm_oMeanRep = totalP_Norm_oCom/totalP_Mean_Rep) %>% # Divide totalP_Norm_oCom by its mean for each replicate. This is to normalize over the mean of all samples in each replicate.
      ungroup()
    
    #Normalization across Gels : Method 1 & 2
    WBData_NormORepMean <- WBData_NormORepMean %>% 
      group_by(Replicate) %>% # remaining calculations are per replicate. Should M1 be per replicate ?
      mutate(Norm_phosP_by_total_M1= (phosP_by_total/mean(phosP_by_total)))%>% 
      mutate(Norm_phosP_by_total_M2= phosP_Norm_oMeanRep/totalP_Norm_oMeanRep)
    
    return(WBData_NormORepMean)
    
  }else{
    WBData_NormORepMean <- WBData_NormOCom_3Replicates %>% 
      mutate(phosP_by_total = phosP/totalP) %>% #each phosP divided by totalP [ Both phosP and totalP values for each sample is already norm wrt common sample(or their mean, if more than one common samples)] 
      group_by(Replicate) %>% # remaining calculations are per replicate
      mutate(phosP_Mean_Rep = mean(phosP))%>%    # For use below, calculate mean of the phosP norm wrt common sample(per replicate) 
      mutate(totalP_Mean_Rep = mean(totalP))%>%  # For use below, calculate mean of the phosP norm wrt common sample(per replicate) 
      mutate(phosP_Norm_oMeanRep = phosP/phosP_Mean_Rep) %>% # Divide phosP_Norm_oCom by its mean for each replicate. This is to normalize over the mean of all samples in each replicate.
      mutate(totalP_Norm_oMeanRep = totalP/totalP_Mean_Rep) %>% # Divide totalP_Norm_oCom by its mean for each replicate. This is to normalize over the mean of all samples in each replicate.
      ungroup()
    
    #Normalization across Gels : Method 1 & 2
    WBData_NormORepMean <- WBData_NormORepMean %>% 
      group_by(Replicate) %>% # remaining calculations are per replicate. Should M1 be per replicate ?
      mutate(Norm_phosP_by_total_M1= (phosP_by_total/mean(phosP_by_total)))%>% 
      mutate(Norm_phosP_by_total_M2= phosP_Norm_oMeanRep/totalP_Norm_oMeanRep)
    
    return(WBData_NormORepMean)
    
    
  }
}

# FC_Calculation <- function(NormOverGels_3Reps){
#   #Unhash To TEST
#   #NormOverGels_3Reps <- CTDR_Akt
#   
#   #### 
#   NormOverGels_3Reps_FC <- NormOverGels_3Reps %>% 
#     dplyr::filter(Cell_line %in% c("XX","XO")) %>%
#     dplyr::mutate( FCoXX_M1 = Norm_phosP_by_total_M1/mean(Norm_phosP_by_total_M1[Cell_line == "XX" & Treatment == 0])) %>% 
#     dplyr::mutate( FCoXX_M2 = Norm_phosP_by_total_M2/mean(Norm_phosP_by_total_M2[Cell_line == "XX" & Treatment == 0]))
#   
#   
#   NormOverGels_3Reps_FC <- NormOverGels_3Reps_FC %>% 
#     dplyr::filter(Cell_line %in% c("XX","XO")) %>%
#     dplyr::group_by(Replicate, Cell_line) %>% 
#     dplyr::mutate(FCoCntrl_M1 = Norm_phosP_by_total_M1/mean(Norm_phosP_by_total_M1[Treatment == 0])) %>% 
#     dplyr::mutate(FCoCntrl_M2 = Norm_phosP_by_total_M2/mean(Norm_phosP_by_total_M2[Treatment == 0]))
#   
# }

FC_Calculation_updated <- function(NormOverGels_3Reps){
  #Unhash To TEST
  #NormOverGels_3Reps <- CTDR_Akt
  
  #### 
  NormOverGels_3Reps_FC <- NormOverGels_3Reps %>% 
    dplyr::filter(Cell_line %in% c("XX","XO")) %>%
    ungroup() %>% 
    dplyr::mutate( FCoXX_M1 = Norm_phosP_by_total_M1/mean(Norm_phosP_by_total_M1[Cell_line == "XX" & Treatment == 0])) %>% 
    dplyr::mutate( FCoXX_M2 = Norm_phosP_by_total_M2/mean(Norm_phosP_by_total_M2[Cell_line == "XX" & Treatment == 0]))
  
  
  NormOverGels_3Reps_FC <- NormOverGels_3Reps_FC %>% 
    dplyr::filter(Cell_line %in% c("XX","XO")) %>%
    ungroup() %>% 
    #dplyr::group_by(Replicate, Cell_line) %>% 
    dplyr::group_by(Cell_line) %>% # remove Replicate form the group_by, because this now calculates the mean of the control from the 3 replicates and calculates fold change of everything from that mean value.
    dplyr::mutate(FCoCntrl_M1 = Norm_phosP_by_total_M1/mean(Norm_phosP_by_total_M1[Treatment == 0])) %>% 
    dplyr::mutate(FCoCntrl_M2 = Norm_phosP_by_total_M2/mean(Norm_phosP_by_total_M2[Treatment == 0]))
  
}


prep_analyte_G12 <- function(TPS_file, analyte_quant_file) {
  
  ###### TO TEST : Unhash these & Hash the lower block ########
  # TPS_file <- R1_G2_TPS
  # analyte_quant_file <- R1_G2_pMek
  # Replicate_label <- "R1"
  # Analyte_label <- "pMek"
  
  ###### TO TEST : Hash these ########
  analyte_quant_file_name <- deparse(substitute(analyte_quant_file))
  Replicate_label <- str_split(analyte_quant_file_name, "_")[[1]][1]
  Analyte_label <-str_split(analyte_quant_file_name, "_")[[1]][3]
  
  
  ##########
  Wells = c("Well_1","Well_2","Well_3","Well_4","Well_5","Well_6","Well_7","Well_8","Well_9","Well_10","Well_11","Well_12","Well_13","Well_14","Well_15","Well_16","Well_17","Well_18")
  Samples = c("D9_cntrl","DKO_cntrl","XX_cntrl","15 min","30 min","1 h","2 h","4 h","24 h","48 h","XO_cntrl","15 min","30 min","1 h","2 h","4 h","24 h","48 h")
  #Treatment_numeric <- c(0,0,0,15,30,60,120,240,1440,2880,0,15,30,60,120,240,1440,2880)
  Treatment_numeric <- c(0,0,0,0.25,0.5,1,2,4,24,48,0,0.25,0.5,1,2,4,24,48)
  
  
  
  TPS_file <- TPS_file
  TPS_file <- data.frame(TPS_file[4:nrow(TPS_file),1:4])
  colnames(TPS_file) <- TPS_file[1,]
  TPS_file <- TPS_file[-1,]
  colnames(TPS_file) <- gsub("Signal","TPS",colnames(TPS_file))
  TPS_file <- TPS_file %>% 
    mutate(Well = Wells)%>% 
    select(c(Well,TPS))
  
  colnames(analyte_quant_file) <- gsub("Image Name","Image_Name",colnames(analyte_quant_file))
  analyte_quant_file <- analyte_quant_file %>% 
    mutate(Well = Wells) %>% 
    select(c(Image_Name,Well,Channel,Signal))
  
  analyte_and_TPS <- merge(analyte_quant_file,TPS_file, by=c("Well"), sort = FALSE) #Sort false keeps the original seq of rows.
  
  analyte_and_TPS <- analyte_and_TPS %>% 
    mutate(Samples = Samples) %>%
    mutate(Cell_line = c("D9","DKO",rep("XX",8),rep("XO",8))) %>% 
    mutate(Treatment = sub(".*_cntrl$","cntrl",Samples)) %>%
    mutate(Treatment_numeric = Treatment_numeric) %>%
    mutate(Analyte = c(rep(Analyte_label,18)))%>% 
    mutate(Replicate = rep(Replicate_label,18)) %>% 
    mutate(Signal = as.numeric(Signal)) %>%
    mutate(TPS = as.numeric(TPS))%>%
    mutate(Signal_by_TPS = Signal/TPS)
  
  analyte_and_TPS$Treatment <- gsub("cntrl","0", analyte_and_TPS$Treatment)
  analyte_and_TPS$Treatment <- factor(analyte_and_TPS$Treatment,
                                      levels = c("0","15 min","30 min","1 h","2 h","4 h","24 h","48 h"))
  
  #Normalization across Gels : Method 1 & 2
  analyte_and_TPS <- analyte_and_TPS %>%
    mutate(Signal_TNorm_M1 = Signal_by_TPS/mean(Signal_by_TPS)) %>% 
    mutate(Signal_TNorm_M2 = (Signal/mean(Signal))/(TPS/mean(TPS)))
  
  
  
  return(analyte_and_TPS)
}
