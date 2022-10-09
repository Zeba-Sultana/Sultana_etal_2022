
# For colour palette
# https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=5


# function that reverses a string by words
reverse_words <- function(string){
  # split string by blank spaces
  #string <- "the big bang theory"
  
  string_split = strsplit(as.character(string), split = "\\+")
  # how many split terms?
  string_length = length(string_split[[1]])
  # decide what to do
  if (string_length == 1) {
    # one word (do nothing)
    reversed_string = string_split[[1]]
  } else {
    # more than one word (collapse them)
    reversed_split = string_split[[1]][string_length:1]
    reversed_string = paste(reversed_split, collapse = "+")
  }
  # output
  return(reversed_string)
}
# seqQ <-  c("Aa+Bb", "Cc+Dd", "Ee+Ff")
# reverse_words(seqQ)

#Plotting Function
Loop_MaxDiff <- function(ALL_DATA,CellLine,Selected_Perturbations,Selected_Analytes){
  
  # UNHASH to debug ####
  # ALL_DATA <- ALL_DATA
  # CellLine <- c("XX","XO")
  # Selected_Perturbations <- XX_minus_XO_Selected_Pert
  # Selected_Analytes <- XX_minus_XO_Selected_Analyte
  
  Selected_Perturbations_reversed <- unlist(lapply(Selected_Perturbations,reverse_words))
  
  ALL_DATA$Treatment_id <- gsub("\\+NA","",ALL_DATA$Treatment_id)
  SingPert <- ALL_DATA
  #%>% filter(is.na(Perturbation2))
  
  SingPert <- SingPert %>% 
    mutate(Category_Temp = paste(x_status,Category,sep = "_"))
  SingPert$Category <- SingPert$Category_Temp
  SingPert <- SingPert %>% 
    select(-c(Category_Temp))
  SingPert$Category = factor(SingPert$Category, levels = c("XX_Exp_XX","XO_Exp_XO","XX_Sim_Init","XO_Sim_Init","XX_Sim_Final","XO_Sim_Final"), 
                             labels= c("Exp XX","Exp XO","Initial Model-XX","Initial Model-XO","Final Model-XX","Final Model-XO"))
  
  
  
  SingPert_CellLine <- SingPert %>%
    filter(x_status %in% CellLine) %>% 
    filter(Treatment_id %in% c(Selected_Perturbations,Selected_Perturbations_reversed)) 
  
 
  
  SingPert_CellLine_long <- SingPert_CellLine  %>% 
    tidyr::pivot_longer(
      cols = starts_with("p",ignore.case = FALSE), 
      names_to = "Analyte", 
      values_to = "mismatch_value") 
  
  
  MyCellLineColours <- c("XX" = "#FF0000", "XO" = "#0000CD", "nn" = "#000005","Common" = "#000019" )
  MyReplicateShapes <- c("R3" = 2, "R4" = 0, "R5" = 6, "Init_model"= 1,"Final_model" =16)

  
  MyCategoryShapes <- c("Exp XX"=21,"Exp XO"=21,
                        "Initial Model-XX"= 23,"Initial Model-XO"=23,
                        "Final Model-XX"= 22, "Final Model-XO" = 22)
  
  MyCategoryColours <- c("Exp XX"="#252525","Exp XO"="#969696",
                         "Initial Model-XX"= "#dd1c77","Initial Model-XO"="#f768a1",
                         "Final Model-XX"= "#006d2c", "Final Model-XO" = "#74c476")
  
  # MyCategoryFill <- c("Exp XX"="#252525","Exp XO"="#FFFFFF",
  #                     "Initial Model-XX"= "#dd1c77","Initial Model-XO"="#FFFFFF",
  #                     "Final Model-XX"= "#006d2c", "Final Model-XO" = "#FFFFFF")
  
  MyCategoryFill <- c("Exp XX"="#252525","Exp XO"="#FFFFFF",
                      "Initial Model-XX"= "#dd1c77","Initial Model-XO"="#fde0dd",
                      "Final Model-XX"= "#006d2c", "Final Model-XO" = "#a1d99b")
  
  
  
  scaleFUN <- function(x) sprintf("%.1f", x)
  Max_Diff_plots = list()
  
  for(i in 1:length(Selected_Perturbations)){
    sub_data <- SingPert_CellLine_long %>% 
      filter((Treatment_id == Selected_Perturbations[i]|Treatment_id == Selected_Perturbations_reversed[i] )& Analyte == Selected_Analytes[i])
    
    Max_Diff_plots[[i]] <- ggplot(sub_data, aes(x=Cat3, y=mismatch_value, color = Category, fill = Category,shape = Category))+
      geom_hline(yintercept = 0, linetype = "dotted", color="darkgrey")+
      
      geom_pointrange(
        mapping = aes(),
        stat = "summary",
        fun.min = min,
        fun.max = max,
        fun = mean,
        fatten = 2.3, #to control the point size separately from the linerange
        size = 0.5,
        position = position_dodge(width = 0.6))+
      
      
      scale_shape_manual(values = MyCategoryShapes) +
      scale_fill_manual(values= MyCategoryFill)+ # fill works only for some point shapes , Eg 21 -25
      scale_color_manual(values= MyCategoryColours)+
      scale_y_continuous(labels=scaleFUN,
                         # position = "right",
                         expand = expansion(mult = c(0.25, 0.25)))+ # To expand the limits such that 0 zero does not lie at exactly the border
      #xlab(paste0(Selected_Perturbations[i]))+
      xlab(toupper(paste0(Selected_Perturbations[i])))+ ## Made the perturbations label to uppercase
      ggtitle(paste0(Selected_Analytes[i]))+
      
      #MyPaperTheme+
      theme( 
        
            # panel.background = element_rect(fill = "lightblue",
            #                                  colour = "lightblue",
            #                                  size = 0.5, linetype = "solid"),
            # 
            axis.text.x = element_blank(),
            #axis.text.x = element_text(angle = 90,size = 7, colour = "black",vjust = 0.5, hjust=1),
            axis.text.y = element_text(angle = 0, size = 6, colour = "black",hjust = 1),
            strip.text.y.left = element_text(angle = 0),
            #strip.text.y.right = element_text(angle = 0),
            
            
            #legend.position ="none",
            
            legend.position ="bottom",
            legend.text = element_text(color = "black", size = 7),
            legend.title = element_blank(),
            legend.key.size = unit(0.3,"cm"),
            
            #legend.justification = c("centre", "top"),
            
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            
            axis.ticks.x  = element_blank(),
            axis.ticks.y = element_line(colour = "black"),
            axis.ticks.length = unit(0.05, "cm"), # to control the length of the tick mark
            
            axis.title.y = element_blank(),  # To remove x axis title
            axis.title.x.bottom  = element_text(angle = 0, size = 7, colour = "black",hjust = 0.5),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            plot.title = element_text(color="black", size=7, hjust = 0.5,vjust = 0))
    
  }
  
  return(Max_Diff_plots)
}


# mRemove legends from individual plots
RemoveLegend_func <- function(plot_List) { 
  for(i in 1:length(plot_List)){
    #i=2
    plot_List[[i]] <- plot_List[[i]] +
      theme(legend.position = "none")
  }
  return(plot_List)
}


#Subsets the data for one analyte 
# --> filters those rows where there is significant diff between XX and XO
# --> arranges the rows in descending order of XX minus XO values(absolute)
Filter_XXXOpval_perAnalyte <- function(t_result_df, analyte_name, TH){
  
  #UNHASH to test
  #analyte_name = "pGsk3"
  
  t_result_df_analyte <- t_result_df %>% 
    select(c(Treatment_id,grep(substituteDirect(analyte_name),colnames(t_result_df)))) %>% 
    mutate(Analyte = rep(substituteDirect(analyte_name),nrow(t_result_df))) %>% 
    mutate(XXminusXO = get(paste0(substituteDirect(analyte_name),".XX_mean"))-get(paste0(substituteDirect(analyte_name),".XO_mean")) )
  
  colnames(t_result_df_analyte) <- gsub(paste0(substituteDirect(analyte_name),"\\."),"",colnames(t_result_df_analyte))
  
  t_result_df_analyte_filt <- t_result_df_analyte %>% 
    filter(XX_XO_pvalue < TH)
  
  t_result_df_analyte_filt <- t_result_df_analyte_filt %>% 
    arrange(desc(abs(XXminusXO)))
  
  return(t_result_df_analyte_filt)
  
}


#Subsets the data for one analyte 
# --> Based on which cell line in argument, filter those rows where pval th is met for that cell line
# --> arranges the rows in descending order of absolute mean value for that cell-line.
Filter_SingleSample_pval_perAnalyte <- function(t_result_df, analyte_name,cell_line){
  
  #UNHASH to test
  # analyte_name = "pGsk3"
  # cell_line = "XO"
  
  t_result_df_analyte <- t_result_df %>% 
    select(c(Treatment_id,grep(substituteDirect(analyte_name),colnames(t_result_df)))) %>% 
    mutate(Analyte = rep(substituteDirect(analyte_name),nrow(t_result_df))) 
  #mutate(XXminusXO = get(paste0(substituteDirect(analyte_name),".XX_mean"))-get(paste0(substituteDirect(analyte_name),".XO_mean")) )
  
  colnames(t_result_df_analyte) <- gsub(paste0(substituteDirect(analyte_name),"\\."),"",colnames(t_result_df_analyte))
  
  t_result_df_analyte_filt <- t_result_df_analyte %>% 
    #filter(XX_0_pvalue < 0.05)
    filter(case_when( cell_line == "XX" ~ XX_0_pvalue < 0.05, #When cell line == "XX", filter based on XX pvalue
                      cell_line == "XO" ~ XO_0_pvalue < 0.05)) #Otherwise, that of XO
  
  if(cell_line == "XX"){
    t_result_df_analyte_filt <- t_result_df_analyte_filt %>% 
      arrange(desc(abs(XX_mean)))
  }else{ 
    t_result_df_analyte_filt <- t_result_df_analyte_filt %>% 
      arrange(desc(abs(XO_mean)))
  }
  
  
  
  return(t_result_df_analyte_filt)
  
}


