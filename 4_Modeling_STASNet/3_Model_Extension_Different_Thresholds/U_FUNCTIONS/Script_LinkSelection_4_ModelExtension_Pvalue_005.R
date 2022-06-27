library(readxl)
library(tidyverse)
library(htmlTable)

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
  
  # To make the table simpler to see removing rank based columns because they were not used
  common_links_byTot_Percentage_DelResidual <- common_links_byTot_Percentage_DelResidual %>% 
    select(!grep("Rank", colnames(common_links_byTot_Percentage_DelResidual)))
  
  
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