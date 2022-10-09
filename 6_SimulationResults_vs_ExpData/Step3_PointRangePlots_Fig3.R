#!/usr/bin/Rscript

library(STASNet)
library(ggplot2)
library(gdata) #read.xls
library(plyr) # always load plyr before dplyr
library(dplyr) 
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library("corrplot", lib.loc="~/Library/R/3.3/library") # for correlation matrices
library(cowplot) # for plot_grid , for get_legend
library(egg) #needed for set_panel_size, to get predefined sizes of figure panels.
#library(tidyverse)

source("../7_ValidationExperiments/ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")
source("./Functions_Plotting_PointRange.R")

########### Reading in ALL_DATA (ie 3 Exp replicates+ Simulation from Initial and Final Models), subset ExpData from here and do t-tests

ALL_DATA <- read.csv("./OUTPUT/ALL_DATA.csv")

ALL_DATA <- ALL_DATA %>%
  select(-c(X))

ALL_DATA <- ALL_DATA %>% 
  mutate(Cat3 = case_when(grepl("Exp", Category) ~ "Exp",
                          grepl("Init", Category) ~ "Initial_Model",
                          grepl("Final", Category) ~ "Final_Model"))

ALL_DATA$Cat3 = factor(ALL_DATA$Cat3, levels = c("Exp","Initial_Model","Final_Model"))

colnames(ALL_DATA) <- gsub("Gsk3_p", "pGSK3", colnames(ALL_DATA) )
colnames(ALL_DATA) <- gsub("Mek_p", "pMEK", colnames(ALL_DATA) )
colnames(ALL_DATA) <- gsub("mTor_p", "pmTOR", colnames(ALL_DATA) )
colnames(ALL_DATA) <- gsub("Akt_p", "pAKT", colnames(ALL_DATA) )
colnames(ALL_DATA) <- gsub("Erk_p", "pERK", colnames(ALL_DATA) )
colnames(ALL_DATA) <- gsub("Stat3_p", "pSTAT3", colnames(ALL_DATA) )
colnames(ALL_DATA) <- gsub("Smad2_p", "pSMAD2", colnames(ALL_DATA) )


## Separating the Exp and Simulation Data :
Exp_DATA <- ALL_DATA %>% 
  filter(grepl("Exp", Category))

Sim_DATA <- ALL_DATA %>% 
  filter(!grepl("Exp", Category))


#### Doing the T-Tests
perturbations <- unique(as.character(Exp_DATA$Treatment_id)) # 53 perturbations in all
analytes <- colnames(Exp_DATA)[1:7] # 9 : 4 From Bioplex("Gsk3b","Mek","mTor","Akt" - Bioplex Stat3 is not to be used) and 5 from WB("ERKp","STAT3p","SMAD2p","bCatenin","GAPDH")

t_result = list() #prep empty list to store t-test results
#res <- t.test(weight ~ group, data = my_data, var.equal = TRUE)
for(perturbation_name in perturbations){
  
  ##Unhash To test
  #perturbation_name = perturbations[2]
  
  each_perturbation_data <-  Exp_DATA %>%
    filter(Treatment_id == perturbation_name)
  
  for (analyte in analytes){
    
    ##Unhash To test
    #analyte = analytes[3]
    
    subdata <- each_perturbation_data %>%
      select(c(substituteDirect(analyte),Category))
    
    # t_result[[perturbation_name]][[analyte]]$a_pvalue <- t.test(get(substituteDirect(analyte)) ~ Category, data = subdata, var.equal = FALSE)$p.value   # This WORKS
    
    subdata_XX <- subdata %>% 
      filter(grepl("XX",Category)) %>% 
      select(substituteDirect(analyte))
    
    subdata_XO <- subdata %>% 
      filter(grepl("XO",Category))%>% 
      select(substituteDirect(analyte))
    
    # t_result[[perturbation_name]][[analyte]]$a_pvalue <- 
    #   t.test(as.data.frame(subdata_XX[,analyte]),as.data.frame(subdata_XO[,analyte]), alternative = "two.sided", var.equal = FALSE)$p.value   
    
    #ANS <- t.test(as.data.frame(subdata_XX),as.data.frame(subdata_XO),alternative = "two.sided", var.equal = FALSE)
    
    t_test_XX_XO <- t.test(as.data.frame(subdata_XX),as.data.frame(subdata_XO),alternative = "two.sided", var.equal = FALSE) #Paired = FALSE by default
    t_result[[perturbation_name]][[analyte]]$XX_mean<- t_test_XX_XO$estimate[1]
    t_result[[perturbation_name]][[analyte]]$XO_mean<- t_test_XX_XO$estimate[2]
    t_result[[perturbation_name]][[analyte]]$XX_XO_pvalue<- t_test_XX_XO$p.value
    
    #One-sample t-tests to check difference from 0
    t_result[[perturbation_name]][[analyte]]$XX_0_pvalue <- t.test(as.data.frame(subdata_XX), mu = 0, alternative = "two.sided")$p.value
    t_result[[perturbation_name]][[analyte]]$XO_0_pvalue <- t.test(as.data.frame(subdata_XO), mu = 0, alternative = "two.sided")$p.value
    
    
    #but I need to save atleast two values in this list so that later when I convert the t-test results to dataframe the column names are correctly retained. For some reason it does not work correctly if I save only pvalue
    
    
  }
}


t_result_df <- ldply(t_result, data.frame) #to convert list to dataframe ldply is from package (plyr)
colnames(t_result_df) <- gsub("\\.id","Treatment_id",colnames(t_result_df))
t_result_df$Treatment_id <- gsub("\\+NA","",t_result_df$Treatment_id)

write.csv(t_result_df, "./OUTPUT/t_result_df.csv")

print(" T-tests on Exp Data : results saved : t_result_df.csv in OUTPUT")

############## Plotting Fig3D : Max difference between XX and XO #########

#Subsets the data for one analyte 
# --> filters those rows where there is significant diff between XX and XO
# --> arranges the rows in descending order of XX minus XO values(absolute)

t_result_df_Akt_filt <- Filter_XXXOpval_perAnalyte(t_result_df, "pAKT", 0.1) # pAKT was very variable, none passed the TH of p<0.05
t_result_df_Gsk3_filt <- Filter_XXXOpval_perAnalyte(t_result_df, "pGSK3", 0.05)
t_result_df_mTor_filt <- Filter_XXXOpval_perAnalyte(t_result_df, "pmTOR", 0.05)
t_result_df_Mek_filt <- Filter_XXXOpval_perAnalyte(t_result_df, "pMEK", 0.05)
t_result_df_Erk_filt <- Filter_XXXOpval_perAnalyte(t_result_df, "pERK", 0.05)
t_result_df_Stat3_filt <- Filter_XXXOpval_perAnalyte(t_result_df, "pSTAT3", 0.05)
t_result_df_Smad2_filt <- Filter_XXXOpval_perAnalyte(t_result_df, "pSMAD2", 0.05)

### Collating the 1st row(corresponding to highest absoulte XXminusXO) from all the above analyte-wise subsets
Selected_XXXODiff_AnalyteWise <- bind_rows(t_result_df_Akt_filt[1,],
                                           t_result_df_Gsk3_filt[1,],
                                           t_result_df_mTor_filt[1,],
                                           t_result_df_Mek_filt[1,],
                                           t_result_df_Erk_filt[1,],
                                           t_result_df_Stat3_filt[1,],
                                           t_result_df_Smad2_filt[1,])

XXXODiff_AnalyteWise_Selected_Pert <- as.character(Selected_XXXODiff_AnalyteWise$Treatment_id) 
XXXODiff_AnalyteWise_Selected_Analyte <- as.character(Selected_XXXODiff_AnalyteWise$Analyte)

XXXODiff_AnalyteWise_Legend <- Loop_MaxDiff(ALL_DATA,c("XX","XO"),XXXODiff_AnalyteWise_Selected_Pert, XXXODiff_AnalyteWise_Selected_Analyte)

df_leg <- get_legend(XXXODiff_AnalyteWise_Legend[[2]]) 

nRows_plots = 1
nCols_plots = 7

XXXODiff_AnalyteWise <- RemoveLegend_func(XXXODiff_AnalyteWise_Legend)

XXXODiff_AnalyteWise <- lapply(
  XXXODiff_AnalyteWise ,
  set_panel_size,
  width = unit(1.3, "cm"),
  height = unit(1.8, "cm"))


final_plot <- grid.arrange(arrangeGrob(grobs=XXXODiff_AnalyteWise, nrow = nRows_plots), 
                           arrangeGrob(nullGrob(),df_leg,nullGrob(), nrow = 1),
                           #top = textGrob("Title Text Holder Here", gp = gpar(fontsize = 12, font = 2)),
                           nrow = 2,
                           heights = c(nRows_plots,0.5))

ggsave(file="Fig3_XXXODiff_AllPert_perAnalyte.pdf", final_plot, dpi=300, path='OUTPUT_PAPER', useDingbats=FALSE, width = (2.3*nCols_plots)+0.5, height = (2.8*nRows_plots)+3, units = "cm" )

print("Saved Fig3D in OUTPUT_PAPER")


############## Plotting Fig3B and C: Highest residual per analyte for initial model of XX and XO #########

### Fetching Mismatch Data(Quick) 
Mismatch_DATA <- read.csv("./OUTPUT/Mismatch_DATA.csv")

Mismatch_DATA <- Mismatch_DATA %>% 
  select(-X)

colnames(Mismatch_DATA) <- gsub("Gsk3_p", "pGSK3", colnames(Mismatch_DATA) )
colnames(Mismatch_DATA) <- gsub("Mek_p", "pMEK", colnames(Mismatch_DATA) )
colnames(Mismatch_DATA) <- gsub("mTor_p", "pmTOR", colnames(Mismatch_DATA) )
colnames(Mismatch_DATA) <- gsub("Akt_p", "pAKT", colnames(Mismatch_DATA) )
colnames(Mismatch_DATA) <- gsub("Erk_p", "pERK", colnames(Mismatch_DATA) )
colnames(Mismatch_DATA) <- gsub("Stat3_p", "pSTAT3", colnames(Mismatch_DATA) )
colnames(Mismatch_DATA) <- gsub("Smad2_p", "pSMAD2", colnames(Mismatch_DATA) )


### Subsetting the initial and completed models of the 2 cell lines 
Init_XX_Accuracy_MM <- Mismatch_DATA %>% 
  filter(x_status == "XX" & Category == "Sim_Init")

Init_XO_Accuracy_MM <- Mismatch_DATA %>% 
  filter(x_status == "XO" & Category == "Sim_Init")

Final_XX_Accuracy_MM <- Mismatch_DATA %>% 
  filter(x_status == "XX" & Category == "Sim_Final")

Final_XO_Accuracy_MM <- Mismatch_DATA %>% 
  filter(x_status == "XO" & Category == "Sim_Final")

############### Plotting Fig3B ######

Init_XX_Accuracy_MM_AllPert <- Init_XX_Accuracy_MM 

Init_XX_Accuracy_MM_long <- Init_XX_Accuracy_MM_AllPert %>% 
  tidyr::pivot_longer(
    cols = starts_with("p",ignore.case = FALSE), 
    names_to = "Analyte", 
    values_to = "mismatch_value")

Init_XX_Accuracy_MM_long <- Init_XX_Accuracy_MM_long %>% 
  arrange(desc(abs(mismatch_value)))


Highest_mismatch_AllPert_perAnalayte <- Init_XX_Accuracy_MM_long %>% 
  group_by(Analyte) %>% 
  slice(which.max(abs(mismatch_value))) ## Best way to get the rows corresponding to the max values !!


Analyte_seq <- c("pAKT","pGSK3","pmTOR","pMEK","pERK","pSTAT3","pSMAD2")

Highest_mismatch_AllPert_perAnalayte$Analyte <- factor(Highest_mismatch_AllPert_perAnalayte$Analyte, levels = Analyte_seq)

Highest_mismatch_AllPert_perAnalayte <- Highest_mismatch_AllPert_perAnalayte %>% 
  arrange(factor(Analyte, levels = Analyte_seq))


Init_XX_Selected_Pert <- as.character(Highest_mismatch_AllPert_perAnalayte$Treatment_id)
Init_XX_Selected_Analyte <- as.character(Highest_mismatch_AllPert_perAnalayte$Analyte)


Mismatch_AllPert_perAnalayte_Legend <- Loop_MaxDiff(ALL_DATA,c("XX"),Init_XX_Selected_Pert, Init_XX_Selected_Analyte)

df_leg <- get_legend(Mismatch_AllPert_perAnalayte_Legend[[2]]) 

nRows_plots = 1
nCols_plots = 7

Mismatch_AllPert_perAnalayte <- RemoveLegend_func(Mismatch_AllPert_perAnalayte_Legend)

Mismatch_AllPert_perAnalayte <- lapply(
  Mismatch_AllPert_perAnalayte ,
  set_panel_size,
  width = unit(1.3, "cm"),
  height = unit(1.8, "cm"))


final_plot <- grid.arrange(arrangeGrob(grobs=Mismatch_AllPert_perAnalayte, nrow = nRows_plots), 
                           arrangeGrob(nullGrob(),df_leg,nullGrob(), nrow = 1),
                           #top = textGrob("Title Text Holder Here", gp = gpar(fontsize = 12, font = 2)),
                           nrow = 2,
                           heights = c(nRows_plots,0.5))

ggsave(file="Fig3_Mismatch_AllPert_perAnalayte_XX.pdf", final_plot, dpi=300, path='OUTPUT_PAPER', useDingbats=FALSE, width = (2.3*nCols_plots)+0.5, height = (2.8*nRows_plots)+3, units = "cm" )


print("Fig3B saved in OUTPUT_PAPER")

############### Plotting Fig3C ######

Init_XO_Accuracy_MM_AllPert <- Init_XO_Accuracy_MM 

Init_XO_Accuracy_MM_long <- Init_XO_Accuracy_MM_AllPert %>% 
  tidyr::pivot_longer(
    cols = starts_with("p",ignore.case = FALSE), 
    names_to = "Analyte", 
    values_to = "mismatch_value")

Init_XO_Accuracy_MM_long <- Init_XO_Accuracy_MM_long %>% 
  arrange(desc(abs(mismatch_value)))
#arrange(desc(abs(phospho_value_Diff)))

Highest_mismatch_AllPert_perAnalayte <- Init_XO_Accuracy_MM_long %>% 
  group_by(Analyte) %>% 
  slice(which.max(abs(mismatch_value))) ## Best way to get the rows corresponding to the max values !!


Analyte_seq <- c("pAKT","pGSK3","pmTOR","pMEK","pERK","pSTAT3","pSMAD2")

Highest_mismatch_AllPert_perAnalayte$Analyte <- factor(Highest_mismatch_AllPert_perAnalayte$Analyte, levels = Analyte_seq)

Highest_mismatch_AllPert_perAnalayte <- Highest_mismatch_AllPert_perAnalayte %>% 
  arrange(factor(Analyte, levels = Analyte_seq))


Init_XO_Selected_Pert <- as.character(Highest_mismatch_AllPert_perAnalayte$Treatment_id)
Init_XO_Selected_Analyte <- as.character(Highest_mismatch_AllPert_perAnalayte$Analyte)


Mismatch_AllPert_perAnalayte_Legend <- Loop_MaxDiff(ALL_DATA,c("XO"),Init_XO_Selected_Pert, Init_XO_Selected_Analyte)

df_leg <- get_legend(Mismatch_AllPert_perAnalayte_Legend[[2]]) 

nRows_plots = 1
nCols_plots = 7

Mismatch_AllPert_perAnalayte <- RemoveLegend_func(Mismatch_AllPert_perAnalayte_Legend)

Mismatch_AllPert_perAnalayte <- lapply(
  Mismatch_AllPert_perAnalayte ,
  set_panel_size,
  width = unit(1.3, "cm"),
  height = unit(1.8, "cm"))


final_plot <- grid.arrange(arrangeGrob(grobs=Mismatch_AllPert_perAnalayte, nrow = nRows_plots), 
                           arrangeGrob(nullGrob(),df_leg,nullGrob(), nrow = 1),
                           #top = textGrob("Title Text Holder Here", gp = gpar(fontsize = 12, font = 2)),
                           nrow = 2,
                           heights = c(nRows_plots,0.5))

ggsave(file="Fig3_Mismatch_AllPert_perAnalayte_XO.pdf", final_plot, dpi=300, path='OUTPUT_PAPER', useDingbats=FALSE, width = (2.3*nCols_plots)+0.5, height = (2.8*nRows_plots)+3, units = "cm" )

print("Fig3C saved in OUTPUT_PAPER")




