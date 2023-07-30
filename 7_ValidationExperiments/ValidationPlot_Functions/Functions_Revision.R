PrepareWBData_rev <- function(Geldb_phosP, Geldb_totalP, GelNum = "GelNum", Analyte = "analyte") {
  
  ######################  IMP NOTE  ###################
  ## CROSS CHECK THE FUNCTION BEFORE USING ELSEWHERE ##
  ## Has some steps that are specific for this setup ##
  
  # unhash For Testing #######
  # Geldb_phosP <- BlotA_pSmad2
  # Geldb_totalP <- BlotA_TotSmad2
  # GelNum = "BlotA"
  # Analyte = "SMAD2"
  ##########################
  
  
  WBData_updated <- merge(Geldb_phosP,Geldb_totalP, by=c("Sample"), sort = FALSE) #Sort false keeps the original seq of rows.
  
  WBData_updated <- WBData_updated %>% 
    separate(Sample, c("Cell_line","x_status","Replicate","Treatment","TreatmentType"), "_") # Needs dplyr and tidyr
  
  # WBData_updated <- WBData_updated %>% 
  #   filter(Cell_line %in% c("XX","XO"))
  
  WBData_updated$Gel <- rep(GelNum,nrow(WBData_updated))
  WBData_updated$Analyte <- rep(Analyte,nrow(WBData_updated))
  WBData_updated <- WBData_updated %>% 
    dplyr::mutate(Well_no=row_number())
  
  WBData_updated$Treatment <- gsub("cntrl", "0", WBData_updated$Treatment)
  WBData_updated$Treatment <- as.numeric(WBData_updated$Treatment)
  WBData_updated$log2Treatment <- log2(WBData_updated$Treatment+1)
  WBData_updated$Treatment_Fctr <- factor(WBData_updated$Treatment)
  
  return(WBData_updated)}
Norm_TotalP <- function(WBData){
  # unhash For Testing #######
  #WBData <- BlotA_ActivinResponse
  ##########################
  
  WBData$phosP = as.numeric(WBData$phosP)
  WBData$totalP = as.numeric(WBData$totalP)
  
  WBData_Norm_TotalP <- WBData %>% 
    mutate(phosP_by_total = phosP/totalP) 
  
  return(WBData_Norm_TotalP)
  
  
}
FC_Calculation_rev<- function(WBData_Norm){
  #Unhash To TEST
  
  
  WBData_Norm_FC <- WBData_Norm %>% 
    ungroup() %>% 
    dplyr::mutate(Mean_cntrl_XX = mean(phosP_by_total[x_status == "XX" & Treatment == 0]))
  
  WBData_Norm_FC <- WBData_Norm_FC %>% 
    ungroup() %>% 
    dplyr::group_by(Cell_line, x_status) %>%
    dplyr::mutate(FCoCntrl = phosP_by_total/mean(phosP_by_total[Treatment == 0])) %>% 
    dplyr::mutate(FCoXX = phosP_by_total/Mean_cntrl_XX)
  #dplyr::select(-c(Gel,Analyte,Well_no,log2Treatment,Treatment_Fctr,phosP,totalP)) 
  
  return(WBData_Norm_FC)
  
}
FC_Calculation_rev_eachXX<- function(WBData_Norm){
  #Unhash To TEST
  
  WBData_Norm_FC <- WBData_Norm %>% 
    ungroup() %>% 
    dplyr::group_by(Cell_line) %>%
    dplyr::mutate(FCoXX = phosP_by_total/mean(phosP_by_total[x_status == "XX" & Treatment == 0]))
  
  WBData_Norm_FC <- WBData_Norm_FC %>% 
    ungroup() %>% 
    dplyr::group_by(Cell_line, x_status) %>%
    dplyr::mutate(FCoCntrl = phosP_by_total/mean(phosP_by_total[Treatment == 0])) %>% 
    
    return(WBData_Norm_FC)
  
}

