#!/usr/bin/env Rscript

library(readxl)
library(tidyverse)
library(htmlTable)


# FUNCTIONS 
filter_ext_list <- function(input_ext_list) {
  
  ###UNHASH to test ####
  #input_ext_list <- INIT_XX_ext_InputFile
  
  filt_ext_list <- input_ext_list %>%
    filter(!(from %in% c("Stat3", "Smad1", "Smad2", "Bmp4","NoLif", "Fgf4","Activin"))) %>%
    filter(adj_pval < 0.005 & value != 0) %>%
    mutate(Res_delta = round(Res_delta,3)) %>%
    mutate(adj_pval = round(adj_pval,3)) %>%
    mutate(Link = paste0(from,"_",to)) %>% 
    mutate(Base_Res = residual+Res_delta) %>% 
    mutate(Perc_Res_delta = round((Res_delta/Base_Res)*100,2) )
  return(filt_ext_list)
}

# Selection_Criterion = 1 : common_links_byComRanks
# Selection_Criterion = 2 : common_links_byTotDelResidual
# Selection_Criterion = 3 : common_links_byTot_Percentage_DelResidual
Select_Link<- function(XX_ext_list, XO_ext_list, Selection_Criterion ) {
  
  ##########  Unhash to test######
  # XX_ext <-  INIT_XX_ext
  # XO_ext <-  INIT_XO_ext
  # Selection_Criterion <- 1
  
  
  ###  Hash for testing. Unhash for real work
  XX_ext <-  XX_ext_list
  XO_ext <-  XO_ext_list
  Selection_Criterion <- Selection_Criterion
  
  XX_ext_common = XX_ext[XX_ext$Link %in% XO_ext$Link,]
  XO_ext_common = XO_ext[XO_ext$Link %in% XX_ext$Link,]
  
  
  XX_order_delRes = order(XX_ext_common$Res_delta, decreasing = TRUE) #in the sorted seq, to which position should the input seqelements be moved to
  XX_ext_ordered = XX_ext_common[XX_order_delRes,]# moving the rows to those postions, keeping column seq as is
  XX_ext_ordered$Rank = rank(-XX_ext_common$Res_delta,ties.method = "min")# Now creating the Rank column - such that equal values get average rank
  
  XO_order_delRes = order(XO_ext_common$Res_delta, decreasing = TRUE)
  XO_ext_ordered = XO_ext_common[XO_order_delRes,]
  XO_ext_ordered$Rank = rank(-XO_ext_common$Res_delta,ties.method = "min")
  
  
  # XX_CL_Rank = dplyr::select(XX_ext_ordered, c(Link, Rank, Res_delta, adj_pval,Perc_Res_delta))
  # XO_CL_Rank = dplyr::select(XO_ext_ordered, c(Link, Rank, Res_delta, adj_pval,Perc_Res_delta))
  
  XX_CL_Rank = XX_ext_ordered %>% 
    dplyr::select(Link, Rank, Res_delta, adj_pval,Perc_Res_delta)
  XO_CL_Rank = XO_ext_ordered %>% 
    dplyr::select(Link, Rank, Res_delta, adj_pval,Perc_Res_delta)
  
  common_links_ranks = dplyr::full_join(XX_CL_Rank,XO_CL_Rank, by="Link",suffix = c("_XX", "_XO") )
  
  common_links_ranks = common_links_ranks %>%
    mutate(SumOfRanks = Rank_XX+Rank_XO) %>%
    mutate(CommonRank = pmax(Rank_XX,Rank_XO)) %>%  #pmax is for pair-wise max 
    mutate(Total_DelResidual = Res_delta_XX+Res_delta_XO) %>% 
    mutate(Tot_Percentage_DelResidual = Perc_Res_delta_XX+Perc_Res_delta_XO)
  
  
  #LA3_common_links_bySumOfRanks = common_links_ranks[order(common_links_ranks$SumOfRanks),]
  common_links_byComRanks = common_links_ranks[order(common_links_ranks$CommonRank),]
  common_links_byTotDelResidual = common_links_ranks[order(-common_links_ranks$Total_DelResidual),] # The minus sign is to order in descending
  common_links_byTot_Percentage_DelResidual = common_links_ranks[order(-common_links_ranks$Tot_Percentage_DelResidual),] # The minus sign is to order in descending
  
  if(Selection_Criterion == 1){
    return(common_links_byComRanks)
  }else{
    if(Selection_Criterion == 2){
      return(common_links_byTotDelResidual)
    }else{
      if(Selection_Criterion == 3){
        return(common_links_byTot_Percentage_DelResidual)
      }
    }
  }
  
  
  
}

PrintOut_Table <- function(my_ext_list) {
  htmlTable(my_ext_list,
            align.header = "c",
            align = "c",
            col.rgroup = c("#F7F7F7","none"), # to get the stripes in alternate rows
            css.cell = rbind(rep("font-size: 1em; padding-left: 1.2em; padding-right: 1.2em;", times=ncol(my_ext_list)),
                             matrix("", ncol=ncol(my_ext_list), nrow=nrow(my_ext_list))))
}

################# INPUT HERE ###############


Ext_files_dir ="/project/ag_schulz/Zeba/MODELING_STASNet/Perc_DelRes_Midpvalue_005_20201029"

######### INIT

INIT_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_100k_2020-10-29/extension_XX_R345_Network_100k+red.csv"))
INIT_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_100k_2020-10-29/extension_XO_R345_Network_100k+red.csv"))


INIT_XX_ext <- filter_ext_list(INIT_XX_ext_InputFile)
INIT_XO_ext <- filter_ext_list(INIT_XO_ext_InputFile)

INIT_common_links_byTot_Percentage_DelResidual <- Select_Link(INIT_XX_ext,INIT_XO_ext,3)
PrintOut_Table(head(INIT_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(INIT_common_links_byTot_Percentage_DelResidual, "INIT_common_links_byTot_Percentage_DelResidual.csv")

########## LA1

LA1_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA1_100k_2020-10-29/extension_XX_R345_Network_LA1_100k+red.csv"))
LA1_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA1_100k_2020-10-29/extension_XO_R345_Network_LA1_100k+red.csv"))

LA1_XX_ext <- filter_ext_list(LA1_XX_ext_InputFile)
LA1_XO_ext <- filter_ext_list(LA1_XO_ext_InputFile)

LA1_common_links_byTot_Percentage_DelResidual <- Select_Link(LA1_XX_ext,LA1_XO_ext,3)
PrintOut_Table(head(LA1_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA1_common_links_byTot_Percentage_DelResidual, "LA1_common_links_byTot_Percentage_DelResidual.csv")

########## LA2

LA2_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA2_100k_2020-10-29/extension_XX_R345_Network_LA2_100k+red.csv"))
LA2_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA2_100k_2020-10-29/extension_XO_R345_Network_LA2_100k+red.csv"))

LA2_XX_ext <- filter_ext_list(LA2_XX_ext_InputFile)
LA2_XO_ext <- filter_ext_list(LA2_XO_ext_InputFile)

LA2_common_links_byTot_Percentage_DelResidual <- Select_Link(LA2_XX_ext,LA2_XO_ext,3)
PrintOut_Table(head(LA2_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA2_common_links_byTot_Percentage_DelResidual, "LA2_common_links_byTot_Percentage_DelResidual.csv")


########## LA3

LA3_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA3_100k_2020-10-29/extension_XX_R345_Network_LA3_100k+red.csv"))
LA3_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA3_100k_2020-10-29/extension_XO_R345_Network_LA3_100k+red.csv"))

LA3_XX_ext <- filter_ext_list(LA3_XX_ext_InputFile)
LA3_XO_ext <- filter_ext_list(LA3_XO_ext_InputFile)

LA3_common_links_byTot_Percentage_DelResidual <- Select_Link(LA3_XX_ext,LA3_XO_ext,3)
PrintOut_Table(head(LA3_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA3_common_links_byTot_Percentage_DelResidual, "LA3_common_links_byTot_Percentage_DelResidual.csv")

########## LA4

LA4_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA4_100k_2020-10-29/extension_XX_R345_Network_LA4_100k+red.csv"))
LA4_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA4_100k_2020-10-29/extension_XO_R345_Network_LA4_100k+red.csv"))

LA4_XX_ext <- filter_ext_list(LA4_XX_ext_InputFile)
LA4_XO_ext <- filter_ext_list(LA4_XO_ext_InputFile)

LA4_common_links_byTot_Percentage_DelResidual <- Select_Link(LA4_XX_ext,LA4_XO_ext,3)
PrintOut_Table(head(LA4_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA4_common_links_byTot_Percentage_DelResidual, "LA4_common_links_byTot_Percentage_DelResidual.csv")

########## LA5

LA5_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA5_100k_2020-10-29/extension_XX_R345_Network_LA5_100k+red.csv"))
LA5_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA5_100k_2020-10-29/extension_XO_R345_Network_LA5_100k+red.csv"))

LA5_XX_ext <- filter_ext_list(LA5_XX_ext_InputFile)
LA5_XO_ext <- filter_ext_list(LA5_XO_ext_InputFile)

LA5_common_links_byTot_Percentage_DelResidual <- Select_Link(LA5_XX_ext,LA5_XO_ext,3)
PrintOut_Table(head(LA5_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA5_common_links_byTot_Percentage_DelResidual, "LA5_common_links_byTot_Percentage_DelResidual.csv")

########## LA6

LA6_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA6_100k_2020-10-29/extension_XX_R345_Network_LA6_100k+red.csv"))
LA6_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA6_100k_2020-10-29/extension_XO_R345_Network_LA6_100k+red.csv"))

LA6_XX_ext <- filter_ext_list(LA6_XX_ext_InputFile)
LA6_XO_ext <- filter_ext_list(LA6_XO_ext_InputFile)

LA6_common_links_byTot_Percentage_DelResidual <- Select_Link(LA6_XX_ext,LA6_XO_ext,3)
PrintOut_Table(head(LA6_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA6_common_links_byTot_Percentage_DelResidual, "LA6_common_links_byTot_Percentage_DelResidual.csv")

########## LA7

LA7_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA7_100k_2020-10-29/extension_XX_R345_Network_LA7_100k+red.csv"))
LA7_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA7_100k_2020-10-29/extension_XO_R345_Network_LA7_100k+red.csv"))

LA7_XX_ext <- filter_ext_list(LA7_XX_ext_InputFile)
LA7_XO_ext <- filter_ext_list(LA7_XO_ext_InputFile)

LA7_common_links_byTot_Percentage_DelResidual <- Select_Link(LA7_XX_ext,LA7_XO_ext,3)
PrintOut_Table(head(LA7_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA7_common_links_byTot_Percentage_DelResidual, "LA7_common_links_byTot_Percentage_DelResidual.csv")

########## LA8

LA8_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA8_100k_2020-10-29/extension_XX_R345_Network_LA8_100k+red.csv"))
LA8_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA8_100k_2020-10-29/extension_XO_R345_Network_LA8_100k+red.csv"))

LA8_XX_ext <- filter_ext_list(LA8_XX_ext_InputFile)
LA8_XO_ext <- filter_ext_list(LA8_XO_ext_InputFile)

LA8_common_links_byTot_Percentage_DelResidual <- Select_Link(LA8_XX_ext,LA8_XO_ext,3)
PrintOut_Table(head(LA8_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA8_common_links_byTot_Percentage_DelResidual, "LA8_common_links_byTot_Percentage_DelResidual.csv")

########## LA9

LA9_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA9_100k_2020-10-30/extension_XX_R345_Network_LA9_100k+red.csv"))
LA9_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA9_100k_2020-10-30/extension_XO_R345_Network_LA9_100k+red.csv"))

LA9_XX_ext <- filter_ext_list(LA9_XX_ext_InputFile)
LA9_XO_ext <- filter_ext_list(LA9_XO_ext_InputFile)

LA9_common_links_byTot_Percentage_DelResidual <- Select_Link(LA9_XX_ext,LA9_XO_ext,3)
PrintOut_Table(head(LA9_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA9_common_links_byTot_Percentage_DelResidual, "LA9_common_links_byTot_Percentage_DelResidual.csv")


########## LA10

LA10_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA10_100k_2020-10-30/extension_XX_R345_Network_LA10_100k+red.csv"))
LA10_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA10_100k_2020-10-30/extension_XO_R345_Network_LA10_100k+red.csv"))

LA10_XX_ext <- filter_ext_list(LA10_XX_ext_InputFile)
LA10_XO_ext <- filter_ext_list(LA10_XO_ext_InputFile)

LA10_common_links_byTot_Percentage_DelResidual <- Select_Link(LA10_XX_ext,LA10_XO_ext,3)
PrintOut_Table(head(LA10_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA10_common_links_byTot_Percentage_DelResidual, "LA10_common_links_byTot_Percentage_DelResidual.csv")

print("All extension links listed")






