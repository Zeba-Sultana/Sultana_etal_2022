
source("/Users/sultana/PHD_PROJECT_Zeba/Zeba_PhD_PAPER/SCRIPTS_Sultana_etal_2021/PS10_ValidationExp/ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

#### FUNCTION 
Plot_OnePanel_ValidationPlot_updated <- function(CombinedData, x_axis_entity, y_axis_entity, Analyte_labels, scale_choice) {
  
  ##TO CHECK - UNHASH
  # CombinedData <- PDDR_All_MekTPSplot_FCoCntrl_pMek
  # x_axis_entity <- "log2Treatment"
  # y_axis_entity <- "Signal"
  # Analyte_labels <- Analyte_labels
  # scale_choice <- "free_y"
  
  
  g <- ggplot(data=CombinedData,aes(x=get(substituteDirect(x_axis_entity)),y=get(substituteDirect(y_axis_entity)),group=Cell_line))+ # When the x-axis is a factor, grouping variable is needed else, a line is not drawn.
    #facet_wrap( ~ Analyte,labeller = labeller(Analyte = Analyte_labels), scales = scale_choice) +
    
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



#### READING IN ALL DATA 
#Read in PDDR - FCoverControl
PDDR_All_FCoCntrl <- read.csv("./OutputFiles/PDDR_All_MekTPSplot_FCoCntrl.csv")

#Read in PDDR - FCoverXX
PDDR_All_FCoXX <- read.csv("./OutputFiles/PDDR_All_MekTPSplot_FCoXX.csv")


#Read in PDTC - FCoverControl
PDTC_All_FCoCntrl <- read.csv("/Users/sultana/PHD_PROJECT_Zeba/Zeba_PhD_PAPER/SCRIPTS_Sultana_etal_2021/PS10_ValidationExp/Meki_TC/OutputFiles/PDTC_pMek_pcRaf_Joinplot_OverCntrl.csv")
PDTC_All_FCoCntrl$Treatment_factor <- as.factor(PDTC_All_FCoCntrl$Treatment_factor)

#Read in PDTC - FCoverXX
PDTC_All_FCoXX <- read.csv("/Users/sultana/PHD_PROJECT_Zeba/Zeba_PhD_PAPER/SCRIPTS_Sultana_etal_2021/PS10_ValidationExp/Meki_TC/OutputFiles/PDTC_pMek_pcRaf_Joinplot.csv")
PDTC_All_FCoXX$Treatment_factor <- as.factor(PDTC_All_FCoXX$Treatment_factor)

############## PDDR+TC : FC over Control : Mek ##########

PDDR_Mek_FCoCntrl <- PDDR_All_FCoCntrl %>% 
  filter(grepl("Mek",Analyte))
PDTC_Mek_FCoCntrl <- PDTC_All_FCoCntrl %>% 
  filter(grepl("Mek",Analyte))

############## PDDR+TC : FC over Control : Raf ##########

PDDR_Raf_FCoCntrl <- PDDR_All_FCoCntrl %>% 
  filter(grepl("Raf",Analyte))
PDTC_Raf_FCoCntrl <- PDTC_All_FCoCntrl %>% 
  filter(grepl("Raf",Analyte))

 
 
############# CREATING PLOTS ############
##### pMEK with PDDR+TC (FC over control)

#pMEK : PDDR
Analyte_labels <- c("FCoCntrl_M2_Mek" = "pMEK", "FCoCntrl_M2_cRaf" = "pcRAF")
PDDR_Mek_max <-max(PDDR_Mek_FCoCntrl$Signal) + 5 # the number is based on how much space i need
dummy_PDDR_Mek <- data.frame(log2Treatment = 0, Signal = PDDR_Mek_max,
                        Analyte = "FCoCntrl_M2_Mek", Cell_line = "Common", stringsAsFactors=FALSE)

p1 <- Plot_OnePanel_ValidationPlot_updated(PDDR_Mek_FCoCntrl,"log2Treatment","Signal", Analyte_labels,"free_y")+
  geom_blank(data=dummy_PDDR_Mek) + 
  labs(x = "\n MEKi(nM)+1 [log2]", #\u03bc is the unicode charachter fro greek mu
       y = " Rel. phosp. (norm.) \n Fold Change over DMSO \n",
       title = Analyte_labels,
       color = "Cell line" )+
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")

g1=set_panel_size(p1,width=unit(2.8,'cm'),height=unit(2.8,'cm'))

#pMEK : PDTC
Analyte_labels <- c("pMek_FCoCntrl_M2" = "pMEK", "pcRaf_FCoCntrl_M2" = "pcRAF")
PDTC_Mek_max <-max(PDTC_Mek_FCoCntrl$Signal) + 5 # the number is based on how much space i need
dummy_PDTC_Mek <- data.frame(Treatment_factor = 0, Signal = PDTC_Mek_max,
                             Analyte = "pMek_FCoCntrl_M2", Cell_line = "Common", stringsAsFactors=FALSE)
p2 <- Plot_OnePanel_ValidationPlot_updated(PDTC_Mek_FCoCntrl,"Treatment_factor","Signal", Analyte_labels,"free_y")+
 geom_blank(data=dummy_PDTC_Mek) + 
  labs(x = "\n Time(h)", #\u03bc is the unicode charachter fro greek mu
       y = "",
       title = Analyte_labels,
       color = "Cell line" )+
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")

g2=set_panel_size(p2,width=unit(2.8,'cm'),height=unit(2.8,'cm'))

plot_list <- list(g1,g2)
nRows_plots = 1
nCols_plots = 2

final_plot <- grid.arrange(arrangeGrob(grobs=plot_list, nrow = nRows_plots), 
                           heights = c(nRows_plots))

ggsave(file="./OutputFiles/pMEK_PDDR_TC_FCoCntrl.pdf", final_plot, dpi=300, path='./', useDingbats=FALSE, width = (2.8*nCols_plots)+3.0, height = (2.8*nRows_plots)+3.5, units = "cm" )

##### pcRaf with PDDR+TC (FC over control)
#pRaf : PDDR
Analyte_labels <- c("FCoCntrl_M2_Mek" = "pMEK", "FCoCntrl_M2_cRaf" = "pcRAF")
PDDR_Raf_max <-max(PDDR_Raf_FCoCntrl$Signal) + 0.5 # the number is based on how much space i need
dummy_PDDR_Raf <- data.frame(log2Treatment = 0, Signal = PDDR_Raf_max,
                             Analyte = "FCoCntrl_M2_Raf", Cell_line = "Common", stringsAsFactors=FALSE)

p1 <- Plot_OnePanel_ValidationPlot_updated(PDDR_Raf_FCoCntrl,"log2Treatment","Signal", Analyte_labels,"free_y")+
  geom_blank(data=dummy_PDDR_Raf) + 
  labs(x = "\n Meki(nM)+1 [log2]", #\u03bc is the unicode charachter fro greek mu
       y = " Rel. phosp. (norm.) \n Fold Change over DMSO \n",
       title = "pcRAF",
       color = "Cell line" )+
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")

g1=set_panel_size(p1,width=unit(2.8,'cm'),height=unit(2.8,'cm'))

#pRaf : PDTC
Analyte_labels <- c("pMek_FCoCntrl_M2" = "pMEK", "pcRaf_FCoCntrl_M2" = "pcRAF")
PDTC_Raf_max <-max(PDTC_Raf_FCoCntrl$Signal) + 0.5 # the number is based on how much space i need
dummy_PDTC_Raf <- data.frame(Treatment_factor = 0, Signal = PDTC_Raf_max,
                             Analyte = "pRaf_FCoCntrl_M2", Cell_line = "Common", stringsAsFactors=FALSE)
p2 <- Plot_OnePanel_ValidationPlot_updated(PDTC_Raf_FCoCntrl,"Treatment_factor","Signal", Analyte_labels,"free_y")+
  geom_blank(data=dummy_PDTC_Raf) + 
  labs(x = "\n Time(h)", #\u03bc is the unicode charachter fro greek mu
       y = "",
       title = "pcRAF",
       color = "Cell line" )+
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")

g2=set_panel_size(p2,width=unit(2.8,'cm'),height=unit(2.8,'cm'))

plot_list <- list(g1,g2)
nRows_plots = 1
nCols_plots = 2

final_plot <- grid.arrange(arrangeGrob(grobs=plot_list, nrow = nRows_plots), 
                           heights = c(nRows_plots))

ggsave(file="./OutputFiles/pRaf_PDDR_TC_FCoCntrl.pdf", final_plot, dpi=300, path='./', useDingbats=FALSE, width = (2.8*nCols_plots)+3.0, height = (2.8*nRows_plots)+3.5, units = "cm" )


################ PLOTS : Fold Change over XX
############## PDDR+TC : FC over XX : Mek ##########

PDDR_Mek_FCoXX <- PDDR_All_FCoXX %>% 
  filter(grepl("Mek",Analyte))
PDTC_Mek_FCoXX <- PDTC_All_FCoXX %>% 
  filter(grepl("Mek",Analyte))

############## PDDR+TC : FC over XX : Raf ##########

PDDR_Raf_FCoXX <- PDDR_All_FCoXX %>% 
  filter(grepl("Raf",Analyte))
PDTC_Raf_FCoXX <- PDTC_All_FCoXX %>% 
  filter(grepl("Raf",Analyte))



##### pMEK with PDDR+TC (FC over XX)

#pMEK : PDDR
Analyte_labels <- c("FCoXX_M2_Mek" = "pMEK", "FCoXX_M2_cRaf" = "pcRAF")
PDDR_Mek_max <-max(PDDR_Mek_FCoXX$Signal) + 2 # the number is based on how much space i need
dummy_PDDR_Mek <- data.frame(log2Treatment = 0, Signal = PDDR_Mek_max,
                             Analyte = "FCoXX_M2_Mek", Cell_line = "Common", stringsAsFactors=FALSE)

p1 <- Plot_OnePanel_ValidationPlot_updated(PDDR_Mek_FCoXX,"log2Treatment","Signal", Analyte_labels,"free_y")+
  geom_blank(data=dummy_PDDR_Mek) + 
  labs(x = "\n MEKi(nM)+1 [log2]", #\u03bc is the unicode charachter fro greek mu
       y = " Rel. phosp. (norm.) \n",
       title = "pMEK",
       color = "Cell line" )+
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")

g1=set_panel_size(p1,width=unit(2.8,'cm'),height=unit(2.8,'cm'))

#pMEK : PDTC
Analyte_labels <- c("pMek_FCoXX_M2" = "pMEK", "pcRaf_FCoXX_M2" = "pcRAF")
PDTC_Mek_max <-max(PDTC_Mek_FCoXX$Signal) + 2 # the number is based on how much space i need
dummy_PDTC_Mek <- data.frame(Treatment_factor = 0, Signal = PDTC_Mek_max,
                             Analyte = "pMek_FCoXX_M2", Cell_line = "Common", stringsAsFactors=FALSE)
p2 <- Plot_OnePanel_ValidationPlot_updated(PDTC_Mek_FCoXX,"Treatment_factor","Signal", Analyte_labels,"free_y")+
  geom_blank(data=dummy_PDTC_Mek) + 
  labs(x = "\n Time(h)", #\u03bc is the unicode charachter fro greek mu
       y = "",
       title = "pMEK",
       color = "Cell line" )+
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")

g2=set_panel_size(p2,width=unit(2.8,'cm'),height=unit(2.8,'cm'))

plot_list <- list(g1,g2)
nRows_plots = 1
nCols_plots = 2

final_plot <- grid.arrange(arrangeGrob(grobs=plot_list, nrow = nRows_plots), 
                           heights = c(nRows_plots))

ggsave(file="./OutputFiles/pMEK_PDDR_TC_FCoXX.pdf", final_plot, dpi=300, path='./', useDingbats=FALSE, width = (2.8*nCols_plots)+3.0, height = (2.8*nRows_plots)+3.5, units = "cm" )

##### pcRaf with PDDR+TC (FC over XX)
#pRaf : PDDR
Analyte_labels <- c("FCoXX_M2_Mek" = "pMEK", "FCoXX_M2_cRaf" = "pcRAF")
PDDR_Raf_max <-max(PDDR_Raf_FCoXX$Signal) + 0.5 # the number is based on how much space i need
dummy_PDDR_Raf <- data.frame(log2Treatment = 0, Signal = PDDR_Raf_max,
                             Analyte = "FCoXX_M2_Raf", Cell_line = "Common", stringsAsFactors=FALSE)

p1 <- Plot_OnePanel_ValidationPlot_updated(PDDR_Raf_FCoXX,"log2Treatment","Signal", Analyte_labels,"free_y")+
  geom_blank(data=dummy_PDDR_Raf) + 
  labs(x = "\n Meki(nM)+1 [log2]", #\u03bc is the unicode charachter fro greek mu
       y = " Rel. phosp. (norm.) \n",
       title = "pcRAF",
       color = "Cell line" )+
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")

g1=set_panel_size(p1,width=unit(2.8,'cm'),height=unit(2.8,'cm'))

#pRaf : PDTC
Analyte_labels <- c("pMek_FCoXX_M2" = "pMEK", "pcRaf_FCoXX_M2" = "pcRAF")
PDTC_Raf_max <-max(PDTC_Raf_FCoXX$Signal) + 0.5 # the number is based on how much space i need
dummy_PDTC_Raf <- data.frame(Treatment_factor = 0, Signal = PDTC_Raf_max,
                             Analyte = "pRaf_FCoXX_M2", Cell_line = "Common", stringsAsFactors=FALSE)
p2 <- Plot_OnePanel_ValidationPlot_updated(PDTC_Raf_FCoXX,"Treatment_factor","Signal", Analyte_labels,"free_y")+
  geom_blank(data=dummy_PDTC_Raf) + 
  labs(x = "\n Time(h)", #\u03bc is the unicode charachter fro greek mu
       y = "",
       title = "pcRAF",
       color = "Cell line" )+
  theme(axis.text.x = element_text(size = 7, face = "plain", angle = 45, vjust = 1.2, hjust = 1)) # previous option was face= "bold")

g2=set_panel_size(p2,width=unit(2.8,'cm'),height=unit(2.8,'cm'))

plot_list <- list(g1,g2)
nRows_plots = 1
nCols_plots = 2

final_plot <- grid.arrange(arrangeGrob(grobs=plot_list, nrow = nRows_plots), 
                           heights = c(nRows_plots))

ggsave(file="./OutputFiles/pRaf_PDDR_TC_FCoXX.pdf", final_plot, dpi=300, path='./', useDingbats=FALSE, width = (2.8*nCols_plots)+3.0, height = (2.8*nRows_plots)+3.5, units = "cm" )




