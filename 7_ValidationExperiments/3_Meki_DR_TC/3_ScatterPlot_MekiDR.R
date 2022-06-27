library(ggplot2)
library(tidyverse) 
library(egg) #needed for set_panel_size


# Create folder OUTPUT_PAPER to save figs used in paper
if(!(file.exists("./OUTPUT_PAPER"))){ # create only if it does not already exist
  dir.create("./OUTPUT_PAPER") 
}


source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

PDDR_All_MekTPSplot_FCoXX <- read.csv("OUTPUT_MekiDR/PDDR_All_MekTPSplot_FCoXX.csv")

PDDR_All_MekTPSplot_FCoXX$Analyte = factor(PDDR_All_MekTPSplot_FCoXX$Analyte, levels =c("FCoXX_M2_Mek", "FCoXX_M2_cRaf")) # This is to get the correct sequence in the Facet_wrap
Analyte_labels <- c("FCoXX_M2_Mek" = "pMEK", "FCoXX_M2_cRaf" = "pRAF1")
PDDR_All_MekTPSplot_FCoXX$Cell_line = factor(PDDR_All_MekTPSplot_FCoXX$Cell_line, levels =c("XX", "XO")) # This is to get the correct sequence in the Facet_wrap



# g <- Plot_TwoPanel_ValidationPlot_updated(PDDR_All_MekTPSplot_FCoXX,"log2Treatment","Signal", Analyte_labels,"free_y")+
#   geom_blank(data=dummy_data) + 
#   labs(x = "\n MEKi(nM)+1 [log2]", #\u03bc is the unicode charachter fro greek mu
#        y = " Phosph. rel. to \n untreated XX \n",
#        color = "Cell line" )
# 
# gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
# grid.arrange(gt)
# ggsave("PDDR_pMekTPS_pcRaf_FCoXX_MEANline.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 

PDDR_All_MekTPSplot_FCoXX_subset_Mek <- PDDR_All_MekTPSplot_FCoXX %>% 
  filter((Cell_line == "XX" & (Treatment_Fctr == 0|Treatment_Fctr == 1000))|(Cell_line == "XO" & (Treatment_Fctr == 0|Treatment_Fctr == 4|Treatment_Fctr == 12|Treatment_Fctr == 1000))) %>%
 # filter(Cell_line == "XO" & (Treatment == 0|Treatment == 4|Treatment == 12)) %>% 
  filter(grepl("Mek", Analyte))

PDDR_All_MekTPSplot_ComparativeFC_subset_Mek <- PDDR_All_MekTPSplot_FCoXX_subset_Mek %>% 
  mutate(Baseline_XX = mean(Signal[Treatment == 0 & Cell_line == "XX"])) %>% 
  mutate(Baseline_XO = mean(Signal[Treatment == 4 & Cell_line == "XO"])) %>% 
  mutate(Comparative_FC = ifelse(Cell_line == "XX", Signal/Baseline_XX, Signal/Baseline_XO) )

  
#PDDR_All_MekTPSplot_FCoXX_subset_Mek$Treatment_Fctr = factor(PDDR_All_MekTPSplot_FCoXX_subset_Mek$Treatment_Fctr,levels = c("0","4","12"))
PDDR_All_MekTPSplot_ComparativeFC_subset_Mek$Treatment_Fctr = factor(PDDR_All_MekTPSplot_ComparativeFC_subset_Mek$Treatment_Fctr,levels = c(0,4,12,1000))
PDDR_All_MekTPSplot_ComparativeFC_subset_Mek$x_axis_factor <- paste0(PDDR_All_MekTPSplot_ComparativeFC_subset_Mek$Cell_line,":",PDDR_All_MekTPSplot_ComparativeFC_subset_Mek$Treatment_Fctr)
PDDR_All_MekTPSplot_ComparativeFC_subset_Mek$x_axis_factor <- factor(PDDR_All_MekTPSplot_ComparativeFC_subset_Mek$x_axis_factor, 
                                                             levels = c("XX:0",
                                                                        "XX:1000",
                                                                        "XO:0",
                                                                        "XO:4",
                                                                        "XO:12",
                                                                        "XO:1000"),
                                                             labels = c("XX:0",
                                                                        "XX:1000",
                                                                        "0",
                                                                        "4",
                                                                        "12",
                                                                        "1000"))




####################

PDDR_All_MekTPSplot_FCoXX_subset_no1000_Mek <- PDDR_All_MekTPSplot_ComparativeFC_subset_Mek %>% 
  filter(Treatment != 1000)

g <- ggplot(PDDR_All_MekTPSplot_FCoXX_subset_no1000_Mek, aes(x=x_axis_factor,y=Signal)) +
  #geom_point(alpha = 0.8,aes(x=x_axis_factor,y=Signal,color=Cell_line))+
  geom_point(aes(x=x_axis_factor,y=Signal,color=Cell_line))+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
  #facet_grid(~Cell_line, space = "free",scales = "free", switch = "x")+
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 1))+
  ylim(0,max(PDDR_All_MekTPSplot_FCoXX_subset_no1000_Mek$Signal)+0.5)+
  labs(x = "\n MEKi(nM)", #\u03bc is the unicode charachter fro greek
       y = " Rel. phosp. (norm.) ",
       title = "pMEK",
       color = "Cell line" )+
  MyPaperTheme

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig7Gi_PDDR_XX_XO_compared_pMek.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 


####### FC compared : NEW #########



PDDR_All_MekTPSplot_FCoXX_subset_newFCPlot_Mek <- PDDR_All_MekTPSplot_ComparativeFC_subset_Mek %>% 
  filter(Treatment != 12) %>% 
  filter(!(Cell_line =="XO" & Treatment == 0))

g <- ggplot(PDDR_All_MekTPSplot_FCoXX_subset_newFCPlot_Mek, aes(x=x_axis_factor,y=Signal)) +
  #geom_point(alpha = 0.8,aes(x=x_axis_factor,y=Signal,color=Cell_line))+
  geom_point(aes(x=x_axis_factor,y=Signal,color=Cell_line))+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
  #facet_grid(~Cell_line, space = "free",scales = "free", switch = "x")+
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 1))+
  ylim(0,max(PDDR_All_MekTPSplot_FCoXX_subset_newFCPlot_Mek$Signal)+0.5)+
  labs(x = "\n MEKi(nM)", #\u03bc is the unicode charachter fro greek
       y = " Rel. phos. (norm) ",
       title = "pMEK",
       color = "Cell line" )+
  MyPaperTheme

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig7H_PDDR_XX_XO_compared_pMek_FC_NEW.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 




####### FC compared : previous #########

PDDR_All_MekTPSplot_ComparativeFC_subset_noXO0_Mek <- PDDR_All_MekTPSplot_ComparativeFC_subset_Mek %>%
  filter(!(Cell_line == "XO" & Treatment == 0))

# g <- ggplot(PDDR_All_MekTPSplot_ComparativeFC_subset_noXO0_Mek, aes(x=x_axis_factor, y=Comparative_FC))+
#   #geom_point(alpha = 0.8,aes(x=x_axis_factor,y=Comparative_FC,color=Cell_line))+
#   geom_point(aes(x=x_axis_factor,y=Comparative_FC,color=Cell_line))+
#   scale_fill_manual(values = MyCellLineColours) +
#   scale_color_manual(values = MyCellLineColours) +
#   ylim(0,max(PDDR_All_MekTPSplot_ComparativeFC_subset_noXO0_Mek$Comparative_FC+2))+
#   #facet_grid(~Cell_line, space = "free",scales = "free", switch = "x")+
#   stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 1))+
#   labs(x = "\n MEKi(nM)", #\u03bc is the unicode charachter fro greek
#        y = " pMEK ",
#        title = "Fold change",
#        color = "Cell line" )+
#   MyPaperTheme
# 
# gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
# grid.arrange(gt)
# ggsave("PDDR_XX_XO_compared_pMek_FoldIncrease_updated_small.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 
# 


  




############# RAF #######

PDDR_All_MekTPSplot_FCoXX_subset_Raf <- PDDR_All_MekTPSplot_FCoXX %>% 
  filter((Cell_line == "XX" & Treatment_Fctr == 0)|(Cell_line == "XO" & (Treatment_Fctr == 0|Treatment_Fctr == 4|Treatment_Fctr == 12))) %>%
  # filter(Cell_line == "XO" & (Treatment == 0|Treatment == 4|Treatment == 12)) %>% 
  filter(grepl("Raf", Analyte))

#PDDR_All_MekTPSplot_FCoXX_subset_Raf$Treatment_Fctr = factor(PDDR_All_MekTPSplot_FCoXX_subset_Raf$Treatment_Fctr,levels = c("0","4","12"))
PDDR_All_MekTPSplot_FCoXX_subset_Raf$Treatment_Fctr = factor(PDDR_All_MekTPSplot_FCoXX_subset_Raf$Treatment_Fctr,levels = c(0,4,12))

PDDR_All_MekTPSplot_FCoXX_subset_Raf$x_axis_factor <- paste0(PDDR_All_MekTPSplot_FCoXX_subset_Raf$Cell_line,":",PDDR_All_MekTPSplot_FCoXX_subset_Raf$Treatment_Fctr)
PDDR_All_MekTPSplot_FCoXX_subset_Raf$x_axis_factor <- factor(PDDR_All_MekTPSplot_FCoXX_subset_Raf$x_axis_factor, 
                                                             levels = c("XX:0",
                                                                        "XO:0",
                                                                        "XO:4",
                                                                        "XO:12"),
                                                             labels = c("XX:0",
                                                                        "0",
                                                                        "4",
                                                                        "12"))


g <- ggplot(PDDR_All_MekTPSplot_FCoXX_subset_Raf, aes(x=x_axis_factor,y=Signal)) +
  #geom_point(alpha = 0.8,aes(x=x_axis_factor,y=Signal,color=Cell_line))+
  geom_point(aes(x=x_axis_factor,y=Signal,color=Cell_line))+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
  #facet_grid(~Cell_line, space = "free",scales = "free", switch = "x")+
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 1))+
  ylim(0,max(PDDR_All_MekTPSplot_FCoXX_subset_Raf$Signal+0.2))+
  labs(x = "\n MEKi(nM)", #\u03bc is the unicode charachter fro greek
       y = " Rel. phosp. (norm.)  ",
       title = "cRAF1",
       color = "Cell line" )+
  MyPaperTheme+
  theme(strip.placement = "outside")

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig7Gii_PDDR_XX_XO_compared_pRaf.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 

############## FOLD CHANGE over UNTREATED CONTROL #########

PDDR_All_MekTPSplot_FCoCntrl <- read.csv("OUTPUT_MekiDR/PDDR_All_MekTPSplot_FCoCntrl.csv")


PDDR_All_MekTPSplot_FCoCntrl$Analyte = factor(PDDR_All_MekTPSplot_FCoCntrl$Analyte, levels =c("FCoCntrl_M2_Mek", "FCoCntrl_M2_cRaf")) # This is to get the correct sequence in the Facet_wrap
Analyte_labels <- c("FCoCntrl_M2_Mek" = "pMEK", "FCoCntrl_M2_cRaf" = "pRAF1")
PDDR_All_MekTPSplot_FCoCntrl$Cell_line = factor(PDDR_All_MekTPSplot_FCoCntrl$Cell_line, levels =c("XX", "XO")) # This is to get the correct sequence in the Facet_wrap

PDDR_All_MekTPSplot_FCoCntrl_subset_Mek <- PDDR_All_MekTPSplot_FCoCntrl %>% 
  #filter((Cell_line == "XX" & (Treatment_Fctr == 0 |Treatment_Fctr == 1000))|(Cell_line == "XO" & (Treatment_Fctr == 0|Treatment_Fctr == 4|Treatment_Fctr == 12|Treatment_Fctr == 1000))) %>%
  filter((Cell_line == "XX" & (Treatment_Fctr == 0|Treatment_Fctr == 1000))|(Cell_line == "XO" & (Treatment_Fctr == 0|Treatment_Fctr == 4|Treatment_Fctr == 12|Treatment_Fctr == 1000))) %>%
  # filter(Cell_line == "XO" & (Treatment == 0|Treatment == 4|Treatment == 12)) %>% 
  filter(grepl("Mek", Analyte))

#PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$Treatment_Fctr = factor(PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$Treatment_Fctr,levels = c("0","4","12"))
PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$Treatment_Fctr = factor(PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$Treatment_Fctr,levels = c(0,4,12,1000))

# g <-ggplot(PDDR_All_MekTPSplot_FCoCntrl_subset_Mek, aes(x=Treatment_Fctr,y=Signal)) +
#   geom_point(aes(x=Treatment_Fctr,y=Signal,color=Cell_line))+
#   scale_fill_manual(values = MyCellLineColours) +
#   scale_color_manual(values = MyCellLineColours) +
#   facet_grid(~Cell_line, space = "free",scales = "free", switch = "x")+
#   stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 1))+
#   labs(x = "\n MEKi(nM)", #\u03bc is the unicode charachter fro greek
#        y = " pMek : Fold change ",
#        title = "MEK inhibitor dose response",
#        color = "Cell line" )+
#   MyPaperTheme+
#   theme(strip.placement = "outside")
# 
# gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
# grid.arrange(gt)
# ggsave("PDDR_XX_XO_compared_pMek_FoldIncrease.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 



############ FOLD CHANGE over respective control ###########

# PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$x_axis_factor <- paste0(PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$Cell_line,":",PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$Treatment_Fctr)
# PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$x_axis_factor <- factor(PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$x_axis_factor, 
#                                                                 levels = c("XX:0",
#                                                                            "XX:1000",
#                                                                            "XO:0",
#                                                                            "XO:4",
#                                                                            "XO:12",
#                                                                            "XO:1000"),
#                                                                 labels = c("XX:0",
#                                                                            "XX:1000",
#                                                                            "0",
#                                                                            "4",
#                                                                            "12",
#                                                                            "1000"))

# g <-ggplot(PDDR_All_MekTPSplot_FCoCntrl_subset_Mek, aes(x=x_axis_factor,y=Signal)) +
#   geom_point(aes(x=x_axis_factor,y=Signal,color=Cell_line))+
#   #geom_text(position = position_dodge(width = 1), aes(x=x_axis_factor, y=-5))+
#   scale_fill_manual(values = MyCellLineColours) +
#   scale_color_manual(values = MyCellLineColours) +
#   #facet_grid(~Cell_line, space = "free",scales = "free", switch = "x")+
#   stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 1))+
#   ylim(0,max(PDDR_All_MekTPSplot_FCoCntrl_subset_Mek$Signal)+2)+
#   labs(x = "\n MEKi(nM)", #\u03bc is the unicode charachter fro greek
#        y = " pMEK : Fold change ",
#        title = "MEK inhibitor dose response",
#        color = "Cell line" )+
#   MyPaperTheme
#   #+theme(axis.text.x=element_blank())
#   # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
#   # NULL
# 
# gt=set_panel_size(g,width=unit(4.5,'cm'),height=unit(2.8,'cm'))
# grid.arrange(gt)
# ggsave("PDDR_XX_XO_compared_pMek_FoldIncrease.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 

print("Script 3 executed : Scatter plots saved")

