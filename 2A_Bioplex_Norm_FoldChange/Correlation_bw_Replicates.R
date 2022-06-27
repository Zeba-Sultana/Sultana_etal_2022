#!/usr/bin/env Rscript

library(ggplot2)
library(tidyverse) 
library(gridExtra)
library(RColorBrewer)
library(egg)

if(!(dir.exists("OUTPUT_PAPER"))){
  dir.create("OUTPUT_PAPER")
}

######### CORRELATION PLOTS ##########

#############
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


######### CORRELATION PLOTS ##########

#Bioplex_MFI_MeanNormalized <- read_excel("./OUTPUT/Bioplex_MFI_MEAN_Normalized.xls")

#BP_Normalized_folder <- "/project/ag_schulz/Zeba/SCRIPTS_Sultana_etal_2021/PS2A_BP_Norm_FC_updated/OUTPUT"
BP_Normalized_folder <- "./OUTPUT"

Bioplex_MFI_MeanNormalized <- openxlsx::read.xlsx(file.path(BP_Normalized_folder,"Bioplex_MFI_MEAN_Normalized.xlsx"))

Bioplex_MFI_MeanNormalized$Replicate <- str_replace(Bioplex_MFI_MeanNormalized$Replicate, "R3", "Rep1") #str_replace is a function from package "stringr"
Bioplex_MFI_MeanNormalized$Replicate <- str_replace(Bioplex_MFI_MeanNormalized$Replicate, "R4", "Rep2")
Bioplex_MFI_MeanNormalized$Replicate <- str_replace(Bioplex_MFI_MeanNormalized$Replicate, "R5", "Rep3")

#colnames(Bioplex_MFI_MeanNormalized)
annotation_cols <- c("Treatment","x_status","Replicate")
analyte_cols <- colnames(Bioplex_MFI_MeanNormalized)[!colnames(Bioplex_MFI_MeanNormalized) %in% annotation_cols]
Cell_Lines_tested <- unique(Bioplex_MFI_MeanNormalized$x_status)

for (cell_line in Cell_Lines_tested ) {
  
  #cell_line = "XX" #TO_DEBUG
  
  for(analyte_name in analyte_cols){
    
    #analyte_name = "Akt" #TO_DEBUG
    
    N_XX_analyte <-  Bioplex_MFI_MeanNormalized %>%
      filter(x_status == cell_line) %>%
      select(analyte_name,annotation_cols)
    
    N_XX_analyte_w <- N_XX_analyte %>%
      spread(Replicate, analyte_name)
    
    N_XX_analyte_w_values <- N_XX_analyte_w %>%
      ungroup() %>%
      select(-c(Treatment,x_status))
    
    N_analyte_cor <- cor(N_XX_analyte_w_values, use = "complete.obs")
    
    N_analyte_cor[lower.tri(N_analyte_cor)] <- NA
    N_analyte_cor_l <- reshape2::melt(N_analyte_cor)
    
    g <- ggplot(N_analyte_cor_l, aes(y = factor(Var1),
                                     x = factor(Var2))) +        ## global aes
      geom_tile(alpha=0) +
      geom_point(aes(colour = value, size = value))  +    ## geom_point for circle illusion
      geom_text(aes(label=round(value,2)), color="black", size=3) +
      scale_color_gradientn(
        name= "Correlation\n",
        colors=c("#D7191C","white","#2C7BB6"),
        breaks=c(-1,0,1),
        labels=c(-1,0,1),
        limits=c(-1, 1)) +
      scale_size(range = c(1, 10), limits=c(0,1))+ ## range is to tune the size of circles, by setting limits between 0 and 1 I could control the size across plots. When setting limits=c(-1,1) the change isn sizes was very little so got rid of the -1 upto 0.
      guides(size=FALSE)+
      xlab("")+
      ylab("")+
      ggtitle(paste0(analyte_name ," : " , cell_line , "  Replicates"))+
      MyPaperTheme+
      theme(legend.position="right",legend.title.align=0.5)
    
    gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
    grid.arrange(gt)
    #ggsave(paste0("Correlation_",analyte_name,"_",cell_line,".png"), gt, dpi=300, path = "./")
    ggsave(paste0("Correlation_",analyte_name,"_",cell_line,".pdf"), gt, dpi=300, useDingbats=FALSE ,path = "./OUTPUT_PAPER")
    
    
  }
}





print("Correlation plots saved in OUTPUT_PAPER")

