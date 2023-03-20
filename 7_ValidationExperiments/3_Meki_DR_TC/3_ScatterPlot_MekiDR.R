library(ggplot2)
library(tidyverse) 
library(egg) #needed for set_panel_size


# Create folder OUTPUT_PAPER to save figs used in paper
if(!(file.exists("./OUTPUT_PAPER"))){
  dir.create("./OUTPUT_PAPER") 
}


source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

PDDR_All_MekTPSplot_FCoXX <- read.csv("OUTPUT_MekiDR/PDDR_All_MekTPSplot_FCoXX.csv")

PDDR_All_MekTPSplot_FCoXX$Analyte = factor(PDDR_All_MekTPSplot_FCoXX$Analyte, levels =c("FCoXX_M2_Mek", "FCoXX_M2_cRaf")) # To get the correct sequence in Facet_wrap
Analyte_labels <- c("FCoXX_M2_Mek" = "pMEK", "FCoXX_M2_cRaf" = "pRAF1")
PDDR_All_MekTPSplot_FCoXX$Cell_line = factor(PDDR_All_MekTPSplot_FCoXX$Cell_line, levels =c("XX", "XO")) # To get the correct sequence in Facet_wrap


PDDR_All_MekTPSplot_FCoXX_subset_Mek <- PDDR_All_MekTPSplot_FCoXX %>% 
  filter((Cell_line == "XX" & (Treatment_Fctr == 0|Treatment_Fctr == 1000))|(Cell_line == "XO" & (Treatment_Fctr == 0|Treatment_Fctr == 4|Treatment_Fctr == 12|Treatment_Fctr == 1000))) %>%
  filter(grepl("Mek", Analyte))

PDDR_All_MekTPSplot_ComparativeFC_subset_Mek <- PDDR_All_MekTPSplot_FCoXX_subset_Mek %>% 
  mutate(Baseline_XX = mean(Signal[Treatment == 0 & Cell_line == "XX"])) %>% 
  mutate(Baseline_XO = mean(Signal[Treatment == 4 & Cell_line == "XO"])) %>% 
  mutate(Comparative_FC = ifelse(Cell_line == "XX", Signal/Baseline_XX, Signal/Baseline_XO) )

  
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



#################### Fig7Gi : Mek

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
  labs(x = "\n MEKi(nM)", 
       y = " Rel. phosp. (norm.) ",
       title = "pMEK",
       color = "Cell line" )+
  MyPaperTheme

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig7Gi_PDDR_XX_XO_compared_pMek.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER")


# Test for significance of results : Fig7G_Mek

Mek_XX_0= PDDR_All_MekTPSplot_FCoXX_subset_Mek %>% 
  filter(Cell_line=="XX" & Treatment_Fctr ==0) 

Mek_XO_0 = PDDR_All_MekTPSplot_FCoXX_subset_Mek %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 0) 

Mek_XO_4 = PDDR_All_MekTPSplot_FCoXX_subset_Mek %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 4)

Mek_XO_12 = PDDR_All_MekTPSplot_FCoXX_subset_Mek %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 12) 


Fig7G_Mek_compared_groups=list("Mek_XO_0"=Mek_XO_0,"Mek_XO_4"=Mek_XO_4,"Mek_XO_12"=Mek_XO_12)
Fig7G_Mek_compared_groups_ttest=list()

for(i in 1:length(Fig7G_Mek_compared_groups)) {
  
  Fig7G_Mek_compared_groups_ttest$compared_group[i] = names(Fig7G_Mek_compared_groups[i])
  
  unpaired_Ttest = t.test(Mek_XX_0$Signal,Fig7G_Mek_compared_groups[[i]]$Signal, paired = FALSE)
  
  Fig7G_Mek_compared_groups_ttest$XX_untreated_mean[i] = unpaired_Ttest$estimate[1]
  Fig7G_Mek_compared_groups_ttest$XO_treated_mean[i] = unpaired_Ttest$estimate[2]
  Fig7G_Mek_compared_groups_ttest$p_value[i] = unpaired_Ttest$p.value
  
}

Fig7G_Mek_compared_groups_ttest <- data.frame((sapply(Fig7G_Mek_compared_groups_ttest,c)))
Fig7G_Mek_compared_groups_ttest$p_value = as.numeric(as.character(Fig7G_Mek_compared_groups_ttest$p_value))
Fig7G_Mek_compared_groups_ttest <- Fig7G_Mek_compared_groups_ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))

WriteXLS::WriteXLS(Fig7G_Mek_compared_groups_ttest, "./OUTPUT_PAPER/Fig7Gi_Mek_Ttest.xls")



####### FC compared : Fig7H #########

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

# Test for significance of results : Fig7H

Fig7H_XX = PDDR_All_MekTPSplot_FCoXX_subset_newFCPlot_Mek 

Fig7H_ttest = list()
the_cell_lines = c("XX","XO")
for(i in 1:length(the_cell_lines)){
  
  data = Fig7H_XX %>% 
    filter(Cell_line == the_cell_lines[i])
  
  treated_signal = data %>% filter(Treatment == 1000) %>% select(Signal)
  base_signal = data %>% filter(Treatment != 1000) %>% select(Signal)
  
  test_result = t.test(base_signal,treated_signal, paired = F)
  
  Fig7H_ttest$cell_line[i] = the_cell_lines[i]
  Fig7H_ttest$base_mean[i] = test_result$estimate[1]
  Fig7H_ttest$treated_mean[i] = test_result$estimate[2]
  Fig7H_ttest$p_value[i] = test_result$p.value
  
}

Fig7H_ttest <- data.frame((sapply(Fig7H_ttest,c)))
Fig7H_ttest$p_value = as.numeric(as.character(Fig7H_ttest$p_value))
Fig7H_ttest <- Fig7H_ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))

WriteXLS::WriteXLS(Fig7H_ttest, "./OUTPUT_PAPER/Fig7H_Ttest.xls")



############# Fig7Gii : RAF 

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

# Test for significance of results : Fig7G_Raf

Raf_XX_0= PDDR_All_MekTPSplot_FCoXX_subset_Raf %>% 
  filter(Cell_line=="XX" & Treatment_Fctr ==0) 

Raf_XO_0 = PDDR_All_MekTPSplot_FCoXX_subset_Raf %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 0) 

Raf_XO_4 = PDDR_All_MekTPSplot_FCoXX_subset_Raf %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 4)

Raf_XO_12 = PDDR_All_MekTPSplot_FCoXX_subset_Raf %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 12) 


Fig7G_Raf_compared_groups=list("Raf_XO_0"=Raf_XO_0,"Raf_XO_4"=Raf_XO_4,"Raf_XO_12"=Raf_XO_12)
Fig7G_Raf_compared_groups_ttest=list()

for(i in 1:length(Fig7G_Raf_compared_groups)) {
  
  Fig7G_Raf_compared_groups_ttest$compared_group[i] = names(Fig7G_Raf_compared_groups[i])
  
  unpaired_Ttest = t.test(Raf_XX_0$Signal,Fig7G_Raf_compared_groups[[i]]$Signal, paired = FALSE)
  
  Fig7G_Raf_compared_groups_ttest$XX_untreated_mean[i] = unpaired_Ttest$estimate[1]
  Fig7G_Raf_compared_groups_ttest$XO_treated_mean[i] = unpaired_Ttest$estimate[2]
  Fig7G_Raf_compared_groups_ttest$p_value[i] = unpaired_Ttest$p.value
  
}

Fig7G_Raf_compared_groups_ttest <- data.frame((sapply(Fig7G_Raf_compared_groups_ttest,c)))
Fig7G_Raf_compared_groups_ttest$p_value = as.numeric(as.character(Fig7G_Raf_compared_groups_ttest$p_value))
Fig7G_Raf_compared_groups_ttest <- Fig7G_Raf_compared_groups_ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))

WriteXLS::WriteXLS(Fig7G_Raf_compared_groups_ttest, "./OUTPUT_PAPER/Fig7Gii_Raf_Ttest.xls")



print("Script 3 executed : Scatter plots saved")
