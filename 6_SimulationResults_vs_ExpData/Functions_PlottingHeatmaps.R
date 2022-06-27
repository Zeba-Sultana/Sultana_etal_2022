


get_CombinationMatrix_new <- function(Data, cell_line,analyte) {
  
  ### TO DEBUG ###  
  # analyte <- "Mek_p" #DEBUG
  # cell_line <- "XX" # DEBUG
  # Data <- Exp_Data
  
  Data_analyte <- Data %>% 
    select(substituteDirect(analyte), paste0(substituteDirect(analyte),".a_pvalue"), x_status, Perturbation1,Perturbation2, ID)
  
  cell_type_data <- Data_analyte %>% 
    filter(x_status == substituteDirect(cell_line))
  
  ## Separately fetching the single treatment data and p-values as well
  SinglePert_Data <- cell_type_data %>% 
    filter(is.na(Perturbation2)) %>% 
    select(substituteDirect(analyte), Perturbation1)
  rownames(SinglePert_Data) <- SinglePert_Data$Perturbation1
  colnames(SinglePert_Data) <- c("phospho_value","Pert1")
  
  SinglePert_Exp_pvalues <- cell_type_data %>% 
    filter(is.na(Perturbation2)) %>% 
    select(paste0(substituteDirect(analyte),".a_pvalue"), Perturbation1)
  
  rownames(SinglePert_Exp_pvalues) <- SinglePert_Exp_pvalues$Perturbation1
  colnames(SinglePert_Exp_pvalues) <- c("p_value","Pert1")
  
  
  
  Perturbation_Vector = c("None","Igfri","Pi3ki","Fgf4","Fgfri","Meki","NoLif","Jaki","Activin","Bmp4ri","Gsk3i") 
  
  Matrix_Data <- data.frame(matrix(NA, nrow = 11,ncol = 11))
  row.names(Matrix_Data) <- Perturbation_Vector
  colnames(Matrix_Data) <- Perturbation_Vector
  
  Matrix_Exp_pvalues <- data.frame(matrix(NA, nrow = 11,ncol = 11))
  row.names(Matrix_Exp_pvalues) <- Perturbation_Vector
  colnames(Matrix_Exp_pvalues) <- Perturbation_Vector
  
  
  for (RN in row.names(Matrix_Data)){
    for(CN in colnames(Matrix_Data)){
      
      if((RN == "Fgfri" & CN == "Meki")||
         (RN == "Meki" & CN =="Fgfri")||
         (RN == "Igfri" & CN == "Pi3ki")||
         (RN == "Pi3ki" & CN =="Igfri")||
         (RN == CN)){
        Matrix_Data[RN,CN] <- NA
        Matrix_Exp_pvalues[RN,CN] <- NA
      }else{ 
        if(RN == "None"||CN == "None"){
          my_value <- cell_type_data %>%
            filter(Perturbation1 == CN |Perturbation1 == RN ) %>%
            filter(is.na(Perturbation2)) 
          Matrix_Data[RN,CN] <- my_value[,analyte] 
          Matrix_Exp_pvalues[RN,CN] <- my_value[,paste0(analyte,".a_pvalue")]
          #a stand for analyte. Just added this .a_ prefix to be sure these variations of mean 
          #and p-value are not R "keywords"
        }else{
          my_value <- cell_type_data %>%
            filter(Perturbation1 == RN | Perturbation1 == CN) %>%
            filter(Perturbation2 == RN | Perturbation2 == CN) 
          Matrix_Data[RN,CN] <- my_value[,analyte] 
          Matrix_Exp_pvalues[RN,CN] <- my_value[,paste0(analyte,".a_pvalue")]
        }
      }
    }
  } # Here ends the for loop for filling up the means and p-values matrices for either XX or XO
  
  
  Mean_pvalue_list = list("Mean_values"=Matrix_Data,"p_values"=Matrix_Exp_pvalues, "Mean_values_SinglePert" = SinglePert_Data, "p_values_SinglePert" = SinglePert_Exp_pvalues)
  
  ## Returing both the means and p-value matrices in the form of a named list
  return(Mean_pvalue_list)
  
  
}
get_CombinationMatrix_Simulation <- function(Data, cell_line,analyte) {
  
  ### TO DEBUG ###  
  # analyte <- "Mek_p" #DEBUG
  # cell_line <- "XX" # DEBUG
  # Data <- Init_Model_SimulationResults
  
  Data_analyte <- Data %>% 
    select(substituteDirect(analyte), x_status, Perturbation1,Perturbation2)
  
  cell_type_data <- Data_analyte %>% 
    filter(x_status == substituteDirect(cell_line))
  
  ## Separately fetching the single treatment data and p-values as well
  SinglePert_Data <- cell_type_data %>% 
    filter(is.na(Perturbation2)) %>% 
    select(substituteDirect(analyte), Perturbation1)
  rownames(SinglePert_Data) <- SinglePert_Data$Perturbation1
  colnames(SinglePert_Data) <- c("phospho_value","Pert1")
  
  SinglePert_Exp_pvalues <- cell_type_data %>%
    filter(is.na(Perturbation2)) %>%
    select(substituteDirect(analyte), Perturbation1)
  rownames(SinglePert_Exp_pvalues) <- SinglePert_Exp_pvalues$Perturbation1
  colnames(SinglePert_Exp_pvalues) <- c("p_value","Pert1")
  
  
 # Perturbation_Vector = c("None","Igfri","Pi3ki","Fgf4","Fgfri","Meki","NoLif","Jaki","Activin","Bmp4ri","Gsk3i") 
  
  Matrix_Data <- data.frame(matrix(NA, nrow = 11,ncol = 11))
  row.names(Matrix_Data) <- Perturbation_Vector
  colnames(Matrix_Data) <- Perturbation_Vector
  
  Matrix_Exp_pvalues <- data.frame(matrix(NA, nrow = 11,ncol = 11))
  row.names(Matrix_Exp_pvalues) <- Perturbation_Vector
  colnames(Matrix_Exp_pvalues) <- Perturbation_Vector
  
  for (RN in row.names(Matrix_Data)){
    for(CN in colnames(Matrix_Data)){
      
      if((RN == "Fgfri" & CN == "Meki")||
         (RN == "Meki" & CN =="Fgfri")||
         (RN == "Igfri" & CN == "Pi3ki")||
         (RN == "Pi3ki" & CN =="Igfri")||
         (RN == CN)){
        Matrix_Data[RN,CN] <- NA
        #Matrix_Exp_pvalues[RN,CN] <- NA
      }else{ 
        if(RN == "None"||CN == "None"){
          my_value <- cell_type_data %>%
            filter(Perturbation1 == CN | Perturbation1 == RN ) %>%
            filter(is.na(Perturbation2)) 
          Matrix_Data[RN,CN] <- my_value[,analyte] 
          #Matrix_Exp_pvalues[RN,CN] <- my_value[,paste0(analyte,".a_pvalue")]
          #a stand for analyte. Just added this .a_ prefix to be sure these variations of mean 
          #and p-value are not R "keywords"
        }else{
          my_value <- cell_type_data %>%
            filter(Perturbation1 == RN | Perturbation1 == CN) %>%
            filter(Perturbation2 == RN | Perturbation2 == CN) 
          Matrix_Data[RN,CN] <- my_value[,analyte] 
          #Matrix_Exp_pvalues[RN,CN] <- my_value[,paste0(analyte,".a_pvalue")]
        }
      }
    }
  } # Here ends the for loop for filling up the means and p-values matrices for either XX or XO
  
  
  
  Mean_pvalue_list = list("Mean_values"=Matrix_Data,"p_values"=Matrix_Exp_pvalues, "Mean_values_SinglePert" = SinglePert_Data, "p_values_SinglePert" = SinglePert_Exp_pvalues)
  
  ## Returing both the means and p-value matrices in the form of a named list
  return(Mean_pvalue_list)
  
  
}
prepMatrixView <- function(XX_AnalyteData,XO_AnalyteData){
  
  #### TO DEBUG ####
  # XX_AnalyteData <- Exp_XX_Mekp
  # XO_AnalyteData <- Exp_XO_Mekp
  
  
  XX_AnalyteData[upper.tri(XX_AnalyteData, diag = TRUE)] <- NA
  XO_AnalyteData[lower.tri(XO_AnalyteData, diag = TRUE)] <- NA
  
  XXXO_AnalyteData <- data.frame(matrix(NA, nrow = 11,ncol = 11))
  row.names(XXXO_AnalyteData) <- Perturbation_Vector
  colnames(XXXO_AnalyteData) <- Perturbation_Vector
  
  XXXO_AnalyteData[lower.tri(XXXO_AnalyteData)] <- XX_AnalyteData[lower.tri(XX_AnalyteData)]
  XXXO_AnalyteData[upper.tri(XXXO_AnalyteData)] <- XO_AnalyteData[upper.tri(XO_AnalyteData)]
  
  # XXXO_AnalyteData_Long <- XXXO_AnalyteData
  # 
  # Comb_Heatmap_pvalues <- as.data.frame(XXXO_AnalyteData) %>% 
  #   mutate(Pert1 = factor(row.names(.), levels=row.names(.))) %>% 
  #   gather(key = Pert2, value = analyte_pvalue, -Pert1, na.rm = TRUE, factor_key = TRUE)
  
  XXXO_AnalyteData <- XXXO_AnalyteData  %>% 
    mutate(Pert1 = factor(row.names(.), levels=row.names(.))) 
  
  XXXO_AnalyteData_Long <- XXXO_AnalyteData %>%
    tidyr::pivot_longer(
      cols = -c(Pert1),
      names_to = "Pert2",
      values_to = "phospho_value")
  
  XXXO_AnalyteData_Long$Pert1 = factor(XXXO_AnalyteData_Long$Pert1, levels = Perturbation_Vector )
  XXXO_AnalyteData_Long$Pert2 = factor(XXXO_AnalyteData_Long$Pert2, levels = Perturbation_Vector )
  
  return(XXXO_AnalyteData_Long)
  
}
Plot_SmallHeatMap_OnlyComb <- function(COMBO_DATA,analyte_name) {
  
  ## TO DEBUG
  # COMBO_DATA <- COMBO_Data_Mekp
  # ALL_single_DATA <- SINGLE_Data_Mekp
  # analyte_name <- "pMek"
  
  COMBO_DATA$Pert1 = factor(COMBO_DATA$Pert1, levels = Perturbation_Vector, labels= Treatment_Labels )
  COMBO_DATA$Pert2 = factor(COMBO_DATA$Pert2, levels = Perturbation_Vector, labels= Treatment_Labels )
  COMBO_DATA$Category = factor(COMBO_DATA$Category, levels = c("Exp","Initial_Model","Final_Model"), labels= c("Exp Data","Initial Model","Final Model"))
  
  
  my_limit_point <- max(abs(COMBO_DATA$phospho_value), na.rm = TRUE)
  My_palette = brewer.pal(n = 7, name = "PuOr")
  
  #  My_cellline_label <- textGrob("XO")
  # ggplot(data = df, aes(x = group, y = a)) + 
  #   geom_bar(stat = "identity") +
  #   annotation_custom(grob = Text1,  xmin = 2, xmax = 3, ymin = 20, ymax = 22)
  
  Comb_plot <- ggplot(data = COMBO_DATA, aes(Pert2, Pert1, fill = phospho_value)) + 
    #scale_y_discrete(expand=c(0,0))+
    scale_y_discrete(limits = rev(levels(COMBO_DATA$Pert1)))+ # To have the seq of treatments going from T1 - T10 from top 
    scale_x_discrete(position = "top") +
    facet_grid(rows = vars(Category), switch = "y") + # switch = "y" moves the y panel label to left side
    #facet_grid(rows = vars(Category)) +
    #annotate(geom="text", x = 11, y = 1, label="XO",color="red")+
    #annotation_custom(grob = "My_cellline_label",  xmin = 11, xmax = 11, ymin = 1, ymax = 1)+
    geom_tile(color = "black",size=0.1) + 
    scale_fill_gradientn(colours=My_palette,
                         values=rescale(c(my_limit_point, (my_limit_point/4), (my_limit_point/8),
                                          0,
                                          -my_limit_point/8, -my_limit_point/4, -my_limit_point)),
                         limits = c(-my_limit_point,my_limit_point),
                         # breaks=c(my_limit_point, (my_limit_point/2),
                         #                  0,
                         #         -my_limit_point/2, -my_limit_point), # breaks only defines which values are labeled on the colourbar
                         guide="colorbar",
                         name="fold change\n  (log2)")+
    #geom_text(aes(Pert2, Pert1, label = round(phospho_value,2)), color = "black", size = 3,nudge_y = -0.1) +
    # geom_text(aes(label=stars), color="black", size=6,nudge_y = 0.1) +
    #ggtitle(pasteDirect(analyte_name)) +
    ggtitle(analyte_name) +
    MyPaperTheme+
    theme(
      plot.title = element_text(color="black", size=7),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(color="black", size = 6, angle = 90,hjust= 0),
      axis.text.x.top = element_text(vjust = 0.5),#axis.text.x.top inherits from axis.text.x. But               vjust needs to separately specified when label on top
      #The x-axis label is rotated 90deg so the hjust and vjust do not behave as expected, they are reversed in a way.
      axis.text.y = element_text(color="black", size = 6),
      #strip.text.y = element_blank(),
      strip.background = element_blank(),
      #strip.text.y = element_blank(),
      strip.text.y = element_text(color='black', size=6),#set the size and bold face of panel headings of the facet grid.
      strip.text.y.left = element_text(angle = 90),
      strip.placement = "outside", # to place the facet label outside of the y-axis tick-text
      legend.position="bottom",
      aspect.ratio = 0.1,
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill="black", colour="black"),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.ticks = element_blank())
  
  
  return(Comb_plot)
  
  
  #print(Single_plot)
  
  #G1 <- plot_grid(Single_plot, Comb_plot, labels=c("", ""), ncol = 2, nrow = 1)
  
  
}
Plot_SmallHeatMap_Single_ExpOnly <- function(COMBO_DATA,ALL_single_DATA,analyte_name) {
  
  # This function is modified from that which was used to plot both the single and combination heatmaps together.
  # I am also not using the facet_grid functionality that was used to plot the Exp data, Init and Finla model simulations in rows
  
  ## TO DEBUG
  # COMBO_DATA <- COMBO_Data_Mekp
  # ALL_single_DATA <- SINGLE_Data_Mekp
  # analyte_name <- "pMek"
  
  
  #Combo Data is being read in only to define common limits of the fold chnage values between analytes 
  my_limit_point <- max(abs(COMBO_DATA$phospho_value), na.rm = TRUE)
  My_palette = brewer.pal(n = 7, name = "PuOr")
  
  # Single_Perturbation_Vector = c("Igfri","Pi3ki","Fgf4","Fgfri","Meki","NoLif","Jaki","Activin","Bmp4ri","Gsk3i") # Removed "None"
  # Single_Treatment_Labels = c("Igfri","Pi3ki","Fgf4","Fgfri","Meki","NoLif","Jaki","Activin","Bmp4ri","Gsk3i") 
  # 
  Single_Perturbation_Vector = c("Igfri","Pi3ki","Fgf4","Fgfri","Meki","NoLif","Jaki","Activin","Bmp4ri","Gsk3i") 
  Single_Treatment_Labels = c("IGFRi","PI3Ki","FGF4","FGFRi","MEKi","NoLIF","JAKi","ActA","BMP4Ri","GSK3i") 
  
  
  # Comb_plot <- ggplot(data = COMBO_DATA, aes(Pert2, Pert1, fill = phospho_value)) + 
  #   scale_y_discrete(limits = rev(levels(COMBO_DATA$Pert1)))+ # To have the seq of treatments going from T1 - T10 from top 
  #   scale_x_discrete(position = "top") +
  #   facet_grid(rows = vars(Category)) +
  #   #annotate(geom="text", x = 11, y = 1, label="XO",color="red")+
  #   #annotation_custom(grob = "My_cellline_label",  xmin = 11, xmax = 11, ymin = 1, ymax = 1)+
  #   geom_tile(color = "black",size=0.1) + 
  #   scale_fill_gradientn(colours=My_palette,
  #                        values=rescale(c(my_limit_point, (my_limit_point/4), (my_limit_point/8),
  #                                         0,
  #                                         -my_limit_point/8, -my_limit_point/4, -my_limit_point)),
  #                        limits = c(-my_limit_point,my_limit_point),
  #                        # breaks=c(my_limit_point, (my_limit_point/2),
  #                        #                  0,
  #                        #         -my_limit_point/2, -my_limit_point), # breaks only defines which values are labeled on the colourbar
  #                        guide="colorbar",
  #                        name="Log FC")+
  #   #geom_text(aes(Pert2, Pert1, label = round(phospho_value,2)), color = "black", size = 3,nudge_y = -0.1) +
  #   # geom_text(aes(label=stars), color="black", size=6,nudge_y = 0.1) +
  #   #ggtitle(pasteDirect(analyte_name)) +
  #   ggtitle(analyte_name) +
  #   MyPaperTheme+
  #   theme(
  #     plot.title = element_text(color="black", size=7),
  #     axis.title.x = element_blank(),
  #     axis.title.y = element_blank(),
  #     axis.text.x = element_text(color="black", size = 6),
  #     axis.text.y = element_text(color="black", size = 6),
  #     #strip.text.y = element_blank(),
  #     strip.background = element_blank(),
  #     strip.text.y = element_blank(),
  #     legend.position="bottom",
  #     aspect.ratio = 1,
  #     panel.grid.major = element_blank(),
  #     panel.border = element_blank(),
  #     panel.background = element_rect(fill="black", colour="black"),
  #     plot.margin = unit(c(0, 0, 0, 0), "cm"),
  #     axis.ticks = element_blank())
  # 
  
  #print(Comb_plot)
  
  ExpOnly_single_DATA <- ALL_single_DATA %>%
    filter(Category == "Exp") 
  
  ExpOnly_single_DATA$Pert1 = factor(ExpOnly_single_DATA$Pert1, levels = Single_Perturbation_Vector, labels= Single_Treatment_Labels )
  ExpOnly_single_DATA$Pert2 = factor(ExpOnly_single_DATA$Pert1, levels = Single_Perturbation_Vector, labels= Single_Treatment_Labels )
  ExpOnly_single_DATA$Category = factor(ExpOnly_single_DATA$Category, levels = c("Exp","Initial_Model" ,"Final_Model"), labels= c("Exp Data","Initial Model","Final Model"))
  ExpOnly_single_DATA$x_status = factor(ExpOnly_single_DATA$x_status, levels = c("XX","XO"), labels= c("XX","XO") )
  
  Single_plot <- ggplot(data = ExpOnly_single_DATA, aes(x_status, Pert1, fill = phospho_value)) + 
    scale_y_discrete(position = "left", limits = rev(levels(ExpOnly_single_DATA$Pert1)))+ # To have the seq of treatments going from T1 - T10 from top 
    scale_x_discrete(position = "top") +
    #scale_y_discrete(position = "right")
    facet_grid(rows = vars(Category), switch = "y") +
    #geom_tile(aes( width = 1),color = "white") + 
    geom_tile(color = "black",size=0.1) + 
    scale_fill_gradientn(colours=My_palette,
                         values=rescale(c(my_limit_point, (my_limit_point/4), (my_limit_point/8),
                                          0,
                                          -my_limit_point/8, -my_limit_point/4, -my_limit_point)),
                         limits = c(-my_limit_point,my_limit_point),
                         # breaks=c(my_limit_point, (my_limit_point/2),
                         #                  0,
                         #         -my_limit_point/2, -my_limit_point), # breaks only defines which values are labeled on the colourbar
                         guide="colorbar",
                         #name="fold change\n   (log2)")+ # removed the title of the colour bar since it was distorting the alignment
                         name=element_blank())+
    #geom_text(aes(x_status, Pert1, label = round(phospho_value,2)), color = "black", size = 3,nudge_y = -0.1) +
    # geom_text(aes(label=stars), color="black", size=6,nudge_y = 0.1) +
    ggtitle(analyte_name) +
    MyPaperTheme+
    theme(
      plot.title = element_text(color="black", size=7),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #axis.text.x = element_blank(),
      # axis.text.x = element_text(color="black", size = 6, angle = 90 ),
      axis.text.x = element_text(color="black", size = 7 ),
      #axis.text.y = element_blank(),
      axis.text.y  = element_text(color="black", size = 7),
      # strip.text.y = element_text(color='black', size=6),#set the size and bold face of panel headings of the facet grid.
      strip.text.y = element_blank(),# Remove the facet grid label.
      #strip.text.y.left = element_text(angle = 0), # For facet grid label.
      #strip.placement = "outside",# To place facet grid label outside of the yaxis text
      legend.position="bottom",
      legend.direction = "vertical",
      legend.justification="right",
      legend.margin=margin(0,0,0,0), # this and the box.margin control the diatnce between the legend and the plot panel
      legend.box.margin=margin(-5,-5,-5,-5),
      aspect.ratio = 5,
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill="black", colour="black"),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.ticks = element_blank())
  
  return(Single_plot)
  
  # Putting the sigle and combined treatments together in one row and 2 cols
  # G1 <- plot_grid(Single_plot, Comb_plot, labels=c("", ""), ncol = 2, nrow = 1)
  # 
  # return(G1)
  
}