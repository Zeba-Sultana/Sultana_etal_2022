#### Setting up the plotting attributes ###
MyCellLineColours <- c("XX" = "#FF0000", "XO" = "#0000CD", "nn" = "#000005","Common" = "#000019" )
MyReplicateShapes <- c("R1" = 2, "R2" = 0, "R3" = 6, "nn"= 4)

Plot_TwoPanel_ValidationPlot_updated <- function(CombinedData, x_axis_entity, y_axis_entity, Analyte_labels, scale_choice) {
  
  ##TO CHECK - UNHASH
  # CombinedData <- PDTC_pMek_pcRaf_Joinplot
  # x_axis_entity <-x_axis_entity
  # y_axis_entity <- y_axis_entity
  # Analyte_labels <- Analyte_labels
  # scale_choice <- "fixed"
  
  
  g <- ggplot(data=CombinedData,aes(x=get(substituteDirect(x_axis_entity)),y=get(substituteDirect(y_axis_entity)),group=Cell_line))+ # When the x-axis is a factor, grouping variable is needed else, a line is not drawn.
    facet_wrap( ~ Analyte,labeller = labeller(Analyte = Analyte_labels), scales = scale_choice) +
    
    geom_point(aes(colour = Cell_line, shape = Replicate), size = 0.5)+
    scale_fill_manual(values = MyCellLineColours) +
    scale_color_manual(values = MyCellLineColours) +
    scale_shape_manual(values = MyReplicateShapes, guide=FALSE) + # to remove the legend for shapes
    stat_summary(fun = mean,geom = "line",aes(color= Cell_line), size =0.5) +
    #stat_summary(fun = mean,geom = "point",aes(color= Cell_line), size= 1.5, alpha = 0.5) +
    
    ########################
  
  scale_y_continuous(limits = c(0, NA)) +
    #labs(x = expression(mu*"M"), # was complicated to get mu with expression()
    #labs(x = parse(text = my_xaxis_label),
    # color within labs,lets me give user defined labels to the attribute in legend
    guides(colour = guide_legend(order = 1, override.aes = list(size = 1.25)))+
           #shape = guide_legend(order = 2, override.aes = list(size = 1.25))) + # override.aes lets you have a different size in the legend key than that used in the plotting 
    MyPaperTheme
  
  return(g)
  
}

# Plot_TwoPanel_ValidationPlot <- function(CombinedData, Analyte_labels, scale_choice) {
#   
#   g <- ggplot(data=CombinedData,aes(x=log2Treatment,y=Signal))+
#     facet_wrap( ~ Analyte,labeller = labeller(Analyte = Analyte_labels), scales = scale_choice) +
#     #facet_wrap( ~ Analyte,labeller = labeller(Analyte = Analyte_labels), scales = "fixed") +
#     
#     geom_point(aes(colour = Cell_line, shape = Replicate), size = 0.5)+
#     scale_fill_manual(values = MyCellLineColours) +
#     scale_color_manual(values = MyCellLineColours) +
#     scale_shape_manual(values = MyReplicateShapes) +
#     stat_summary(fun = mean,geom = "line",aes(color= Cell_line), size =0.5) +
#     #stat_summary(fun = mean,geom = "point",aes(color= Cell_line), size= 1.5, alpha = 0.5) +
#     
#     
#     ########################
#   
#   scale_y_continuous(limits = c(0, NA)) +
#     #labs(x = expression(mu*"M"), # was complicated to get mu with expression()
#     #labs(x = parse(text = my_xaxis_label),
#     # color within labs,lets me give user defined labels to the attribute in legend
#     guides(colour = guide_legend(order = 1, override.aes = list(size = 1.25)),
#            shape = guide_legend(order = 2, override.aes = list(size = 1.25))) + # override.aes lets you have a different size in the legend key than that used in the plotting 
#     MyPaperTheme
#   
#   return(g)
#   
# }

Plot_ReplicatePoints_MeanLine<- function(plot_data,xValues,yValues){
  xValues=as.character(substitute(xValues))
  yValues=as.character(substitute(yValues))
  
  ######## FOR TESTING ###########
  # plot_data <- Fgf4_TC_Erk_plot_Num
  # xValues=as.character(substitute(Treatment))
  # yValues=as.character(substitute(pErk_by_Erk))
  ######## FOR TESTING ###########
  
  #Cell_Line_Colors <- c("XX" = "#D9717D", "XO" = "#4DB6D0", "8" = "#BECA55")
  
  # plot_data_MeanLine=summarySE(plot_data, measurevar=yValues, groupvars=c("Cell_line","Treatment"))
  # plot_data_MeanLine <- plot_data_MeanLine %>% 
  #   dplyr::mutate(log2Treatment = log2(Treatment+0.5))
  
  ggplot(data=plot_data,aes(plot_data[[xValues]],plot_data[[yValues]], group = Cell_line))+
    geom_point(aes(colour = Cell_line, shape = Replicate), size = 2)+
    #scale_shape(solid = FALSE)+ #This was just to remove the solid filling that came by default
    # scale_shape_manual(values=c(2, 0, 6))+ # By this I manually specify the shape numbers I want(which were hollow, so solid=FALSE not needed anymore)
    # scale_fill_manual(values=c("#0000CD","#FF0000"))+
    # scale_color_manual(values=c("#0000CD","#FF0000"))+
    scale_shape_manual(values= MyReplicateShapes)+ # By this I assign shape numbers to replicates. These are defined at the top of the plotting script
    scale_fill_manual(values= MyCellLineColours)+ # By this I assign colors to cell lines. These are defined at the topr of the plotting script
    scale_color_manual(values= MyCellLineColours)+
    
    #geom_line(aes(x= plot_data_MeanLine[[xValues]],y = plot_data_MeanLine[[yValues]], colour = Cell_line ))+
    stat_summary(fun = mean,geom = "line",aes(color= Cell_line), size =1) +
    stat_summary(fun = mean,geom = "point",aes(color= Cell_line), size= 3) +
    
    scale_y_continuous(limits = c(0, NA)) +
    # labs(x = "\nlog2(CHIR conc+1)",
    #      y = "Phospho-Akt norm. over XX cntrl\n")+
    MyHTMLTheme
  
}

Plot_Ribbons<- function(plot_data,T_Doses,xValues,yValues ){
  
  xValues=as.character(substitute(xValues))
  yValues=as.character(substitute(yValues))
  
  ######## FOR TESTING ###########
  # plot_data <- PD_DR_Mek_plot_Num
  # xValues=as.character(substitute(log2Treatment))
  # yValues=as.character(substitute(pMek_by_Mek_FCoXX))
  ######## FOR TESTING ###########
  
  
  plot_data_ribbon=summarySE(plot_data, measurevar=yValues, groupvars=c("Cell_line","Treatment"))
  plot_data_ribbon <- plot_data_ribbon %>% 
    dplyr::mutate(log2Treatment = log2(Treatment+1)) %>% 
    dplyr::mutate(Treatment_Fctr= factor(Treatment, levels = T_Doses))
  
  plot_data_ribbon %>% 
    ggplot(aes(.data[[xValues]],.data[[yValues]], colour = Cell_line, group = Cell_line)) +
    # scale_fill_manual(values=c("#0000CD","#FF0000"))+
    # scale_color_manual(values=c("#0000CD","#FF0000"))+
    scale_fill_manual(values= MyCellLineColours)+ # By this I assign colors to cell lines. These are defined at the topr of the plotting script
    scale_color_manual(values= MyCellLineColours)+
    geom_ribbon(aes(ymax = .data[[yValues]] + sd, ymin = .data[[yValues]] - sd, fill=Cell_line ),
                alpha = 0.3,
                colour=NA) +
    geom_line(aes(y = .data[[yValues]], colour = Cell_line ))+
    geom_point(aes(colour = Cell_line), size = 2) +
    scale_y_continuous(limits = c(0, NA)) +
    # labs(x = "\nTreatment",
    #      y = "Phosphorylated/Total\n")+
    MyHTMLTheme
  
}

Plot_RawValues <- function(plot_data,xValues,yValues,HLINE_value) {
  
  xValues=as.character(substitute(xValues))
  yValues=as.character(substitute(yValues))
  HLINE_value = as.character(substitute(HLINE_value))
  
  ######## FOR TESTING ###########
  # plot_data <- Fgf4TC_Erk_NormOCom_all
  # xValues=as.character(substitute(Well_no))
  # yValues=as.character(substitute(phosP))
  # HLINE_value= as.character(substitute(Mean_Com_phosP))
  ######## FOR TESTING ###########
  
  
  ###########
  
  ggplot(plot_data, aes(x= factor(.data[[xValues]]), y = .data[[yValues]]))+
    geom_point(aes(colour = Cell_line, shape=Replicate))+
    #scale_shape_manual(values=c(4, 2, 0, 6)) + # By this I manually specify the shape numbers I want(which were hollow, so solid=FALSE not needed anymore)
    scale_shape_manual(values = MyReplicateShapes) + 
    scale_fill_manual(values = MyCellLineColours) +
    scale_color_manual(values = MyCellLineColours) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    facet_wrap(.~Gel) +
    geom_hline(aes(yintercept = .data[[HLINE_value]])) +
    MyHTMLTheme
}

Plot_NormValues_Check <- function(plot_data,xValues,yValues) {
  
  xValues=as.character(substitute(xValues))
  yValues=as.character(substitute(yValues))
  
  ######## FOR TESTING ###########
  # plot_data <- ActDR_Smad2_NormOCom_all
  # xValues=as.character(substitute(Well_no))
  # yValues=as.character(substitute(phosP))
  # HLINE_value= as.character(substitute(Mean_Com_phosP))
  ######## FOR TESTING ###########
  
  ggplot(plot_data, aes(x= factor(.data[[xValues]]), y = .data[[yValues]]))+
    geom_point(aes(colour = Cell_line, shape=Replicate))+
    scale_shape_manual(values= MyReplicateShapes)+ # By this I assign shape numbers to replicates. These are defined at the top of the plotting script
    scale_fill_manual(values= MyCellLineColours)+ # By this I assign colors to cell lines. These are defined at the topr of the plotting script
    scale_color_manual(values= MyCellLineColours)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    facet_wrap(.~Gel)+
    ylim(0,NA)+
    MyHTMLTheme
}

####### FUNCTIONS for Activin and Activin receptor inhibitor ##########
Plot_ReplicatePoints_MeanLine_PaperFig_LOG2y<- function(plot_data,xValues,yValues){
  xValues=as.character(substitute(xValues))
  yValues=as.character(substitute(yValues))
  
  ######## FOR TESTING ###########
  # plot_data <- SB_DR_Smad2
  # xValues=as.character(substitute(log2Treatment))
  # yValues=as.character(substitute(pSmad2_by_Smad2_FCoXX))
  ######## FOR TESTING ###########
  
  
  # plot_data_MeanLine=summarySE(plot_data, measurevar=yValues, groupvars=c("Cell_line","Treatment","Activin_Signaling","log2Activin_Signaling" ))
  # plot_data_MeanLine <- plot_data_MeanLine %>%
  #   dplyr::mutate(log2Treatment = log2(Treatment+1))%>%
  #   dplyr::mutate(log2Activin_Signaling = log2(Activin_Signaling+1))
  
  ggplot(data=plot_data,aes(x=.data[[xValues]],y=.data[[yValues]]))+
    geom_point(aes(colour = Cell_line, shape = Replicate), size = 0.5)+
    #scale_shape(solid = FALSE)+ #This was just to remove the solid filling that came by default
    #scale_shape_manual(values=c(2, 0, 6))+ # By this I manually specify the shape numbers I want(which were hollow, so solid=FALSE not needed anymore)
    scale_shape_manual(values=MyReplicateShapes, guide=FALSE)+
    # scale_fill_manual(values=c("#0000CD","#FF0000"))+
    # scale_color_manual(values=c("#0000CD","#FF0000"))+
    scale_fill_manual(values=MyCellLineColours)+
    scale_color_manual(values=MyCellLineColours)+
    #scale_x_continuous(labels=My_Activin_dose_labels)+
    scale_x_continuous()+
    #scale_x_continuous(breaks=unique(plot_data$log2Treatment),labels = as.character(unique(plot_data$Treatment)))+ # To get used defined labels on the x axis
    
    #geom_line(aes(x= plot_data_MeanLine[[xValues]],y = plot_data_MeanLine[[yValues]], colour = Cell_line ))+
    stat_summary(fun = mean,geom = "line",aes(color= Cell_line), size =0.5) +
    #stat_summary(fun = mean,geom = "point",aes(color= Cell_line), size= 1.5, alpha = 0.5) +
    
    #scale_y_continuous(limits = c(0, NA)) +
    ylim(-1.5, 2)+
    labs(x = paste0("\n Activin(ng/ml)+1 [log2]"),
         y = "Rel.phos.(norm)[log2]\n")+
    guides(colour = guide_legend(order = 1, override.aes = list(size = 1.25)))+
           #shape = guide_legend(order = 2, override.aes = list(size = 1.5))) + # override.aes lets you have a different size in the legend key than that used in the plotting
    MyPaperTheme
  
}
Plot_ReplicatePoints_MeanLine_PaperFig_Minus_LOG2y<- function(plot_data,xValues,yValues){
  xValues=as.character(substitute(xValues))
  yValues=as.character(substitute(yValues))
  
  ######## FOR TESTING ###########
  # plot_data <- SB_DR_Smad2
  # xValues=as.character(substitute(log2Treatment))
  # yValues=as.character(substitute(pSmad2_by_Smad2_FCoXX))
  ######## FOR TESTING ###########
  
  
  # plot_data_MeanLine=summarySE(plot_data, measurevar=yValues, groupvars=c("Cell_line","Treatment","Activin_Signaling","log2Activin_Signaling" ))
  # plot_data_MeanLine <- plot_data_MeanLine %>%
  #   dplyr::mutate(log2Treatment = log2(Treatment+1))%>%
  #   dplyr::mutate(log2Activin_Signaling = log2(Activin_Signaling+1))
  
  ggplot(data=plot_data,aes(x=-1*(.data[[xValues]]),y=.data[[yValues]]))+
    geom_point(aes(colour = Cell_line, shape = Replicate), size = 0.5)+
    #scale_shape(solid = FALSE)+ #This was just to remove the solid filling that came by default
    #scale_shape_manual(values=c(2, 0, 6))+ # By this I manually specify the shape numbers I want(which were hollow, so solid=FALSE not needed anymore)
    scale_shape_manual(values=MyReplicateShapes, guide=FALSE)+
    # scale_fill_manual(values=c("#0000CD","#FF0000"))+
    # scale_color_manual(values=c("#0000CD","#FF0000"))+
    scale_fill_manual(values=MyCellLineColours)+
    scale_color_manual(values=MyCellLineColours)+
    #scale_x_continuous(labels=My_Activin_dose_labels)+
    scale_x_continuous()+
    #scale_x_continuous(breaks=-1*(unique(plot_data$log2Treatment)),labels = as.character(unique(plot_data$Treatment)))+ # To get used defined labels on the x axis
    
    #geom_line(aes(x= plot_data_MeanLine[[xValues]],y = plot_data_MeanLine[[yValues]], colour = Cell_line ))+
    stat_summary(fun = mean,geom = "line",aes(color= Cell_line), size =0.5) +
    #stat_summary(fun = mean,geom = "point",aes(color= Cell_line), size= 1.5, alpha = 0.5) +
    
    #scale_y_continuous(limits = c(0, NA)) +
    ylim(-1.5, 2)+
    labs(x = paste0("\n ActRi(\u03bcM)+1 [log2]"),
         y = "Rel.phos.(norm)[log2]\n")+
    guides(colour = guide_legend(order = 1, override.aes = list(size = 1.25)))+
           #shape = guide_legend(order = 2, override.aes = list(size = 1.5))) + # override.aes lets you have a different size in the legend key than that used in the plotting
    MyPaperTheme
  
}

MyPaperTheme <-  theme_set(theme_light(base_size=8))+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        #strip.background = element_rect(fill="grey90", colour="grey90",size=1),
        panel.background = element_rect(colour="black",size=1, fill = NA),
        # strip.text.x = element_text(color='black',margin = margin(.05, 0, .05, 0, "cm")), 
        # strip.text.y = element_text(color='black',margin = margin(.05, 0.05, .05, 0.05, "cm")),
        #axis.line = element_line(colour = 'black'), # This is used to have only the x and y axes of the plot in black.
        legend.text = element_text(size=6, margin = margin(l = 1, unit = "pt")),
        legend.title=element_text(size=8),
        legend.spacing.y = unit(0.01, 'cm'),
        #legend.key.size = unit(0.5, 'cm'),
        legend.key.size = unit(0.5, "lines"),    # key size (unit) # This controls the distance between legend icons
        legend.key.height = NULL,                # key height (unit)
        legend.key.width = NULL,                 # key width (unit)
        #legend.position = c(.95, .95),         # to place legend within the plot
        legend.position="bottom",
        legend.justification = c("centre", "top"),
        legend.box.background = element_rect(fill=NA, colour = NA),
        legend.background=element_blank(), ### ADDED to get rid of grey bg in legend key
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        axis.text=element_text(colour = "black"),
        axis.text.x = element_text(size = 7, face = "plain"), # previous option was face= "bold"
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 7, face = "plain"),
        axis.title.y = element_text(size = 8),
        axis.ticks = element_line(color='black'),
        axis.ticks.length = unit(0.1, "cm"), # to control the length of the tick mark
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text.x = element_text(color='black', size=8, face="bold"),#set the size and bold face of panel headings of the facet grid.
        strip.text.y = element_text(color='black', size=8, face="bold"),# Not needed, since we dont have a facet grid panel heading on y axis
        strip.background = element_rect(colour=NA, fill=NA))

MyHTMLTheme <-  theme_set(theme_light(base_size=15))+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #strip.background = element_rect(fill="grey90", colour="grey90",size=1),
        panel.background = element_rect(colour="black",size=1, fill = NA),
        # strip.text.x = element_text(color='black',margin = margin(.05, 0, .05, 0, "cm")), 
        # strip.text.y = element_text(color='black',margin = margin(.05, 0.05, .05, 0.05, "cm")),
        #axis.line = element_line(colour = 'black'), # This is used to have only the x and y axes of the plot in black.
        legend.text = element_text(size=12, margin = margin(l = 1, unit = "pt")),
        legend.title=element_text(size=15),
        legend.spacing.y = unit(0.01, 'cm'),
        #legend.key.size = unit(0.5, 'cm'),
        legend.key.size = unit(0.5, "lines"),    # key size (unit) # This controls the distance between legend icons
        legend.key.height = NULL,                # key height (unit)
        legend.key.width = NULL,                 # key width (unit)
        #legend.position = c(.95, .95),         # to place legend within the plot
        legend.position="bottom",
        legend.justification = c("centre", "top"),
        legend.box.background = element_rect(fill=NA, colour = NA),
        legend.background=element_blank(), ### ADDED to get rid of grey bg in legend key
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(size=12, face = "bold"),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=12, face = "bold"),
        axis.title.y = element_text(size=15),
        axis.ticks = element_line(color='black'),
        strip.text.x = element_text(color='black', size=15, face="bold"),#set the size and bold face of panel headings of the facet grid.
        strip.text.y = element_text(color='black', size=15, face="bold"),# Not needed, since we dont have a facet grid panel heading on y axis
        strip.background = element_rect(colour=NA, fill=NA))


