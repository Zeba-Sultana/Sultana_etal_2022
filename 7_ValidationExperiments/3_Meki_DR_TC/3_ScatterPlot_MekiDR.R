library(ggplot2)
library(tidyverse) 
library(egg) #needed for set_panel_size


# Create folder OUTPUT_PAPER to save figs used in paper
if(!(file.exists("./OUTPUT_PAPER"))){
  dir.create("./OUTPUT_PAPER") 
}


source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

PDDR_All_FCoXX <- read.csv("OUTPUT_MekiDR/PDDR_All_MekTPSplot_FCoXX.csv")

PDDR_All_FCoXX$Analyte = factor(PDDR_All_FCoXX$Analyte, levels =c("FCoXX_M2_Mek", "FCoXX_M2_cRaf")) # To get the correct sequence in the Facet_wrap
Analyte_labels <- c("FCoXX_M2_Mek" = "pMEK", "FCoXX_M2_cRaf" = "pRAF1")
PDDR_All_FCoXX$Cell_line = factor(PDDR_All_FCoXX$Cell_line, levels =c("XX", "XO")) 


PDDR_All_FCoXX_SUBSET = PDDR_All_FCoXX %>% 
  filter(Treatment %in% c(0,4,12,1000))

PDDR_All_FCoXX_SUBSET$Treatment_Fctr = factor(PDDR_All_FCoXX_SUBSET$Treatment_Fctr,levels = c(0,4,12,1000))
PDDR_All_FCoXX_SUBSET$x_axis_factor <- paste0(PDDR_All_FCoXX_SUBSET$Cell_line,":",PDDR_All_FCoXX_SUBSET$Treatment_Fctr)
PDDR_All_FCoXX_SUBSET$x_axis_factor <- factor(PDDR_All_FCoXX_SUBSET$x_axis_factor, 
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


########### Fig 7J : Mek ###########

Fig7J_pMek = PDDR_All_FCoXX_SUBSET %>% 
  filter(Analyte == "FCoXX_M2_Mek") %>% 
  filter(Treatment != 1000) %>% 
  filter(!(Cell_line == "XX" & (Treatment == 4|Treatment == 12)))

g <- ggplot(Fig7J_pMek, aes(x=x_axis_factor,y=Signal)) +
  geom_point(aes(x=x_axis_factor,y=Signal,color=Cell_line))+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 1))+
  ylim(0,max(Fig7J_pMek$Signal)+0.5)+
  labs(x = "\n MEKi(nM)", 
       y = " Rel. phosp. (norm.) ",
       title = "pMEK",
       color = "Cell line" )+
  MyPaperTheme

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig7J_pMek.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 

# Test for significance of results : Fig7J_Mek

Mek_XX_0= Fig7J_pMek %>% 
  filter(Cell_line=="XX" & Treatment_Fctr ==0) 

Mek_XO_0 = Fig7J_pMek %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 0) 

Mek_XO_4 = Fig7J_pMek %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 4)

Mek_XO_12 = Fig7J_pMek %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 12) 


Fig7J_Mek_compared_groups=list("Mek_XO_0"=Mek_XO_0,"Mek_XO_4"=Mek_XO_4,"Mek_XO_12"=Mek_XO_12)
Fig7J_Mek_compared_groups_ttest=list()

for(i in 1:length(Fig7J_Mek_compared_groups)) {
  
  Fig7J_Mek_compared_groups_ttest$compared_group[i] = names(Fig7J_Mek_compared_groups[i])
  
  unpaired_Ttest = t.test(Mek_XX_0$Signal,Fig7J_Mek_compared_groups[[i]]$Signal, paired = FALSE)
  
  Fig7J_Mek_compared_groups_ttest$XX_untreated_mean[i] = unpaired_Ttest$estimate[1]
  Fig7J_Mek_compared_groups_ttest$XO_treated_mean[i] = unpaired_Ttest$estimate[2]
  Fig7J_Mek_compared_groups_ttest$p_value[i] = unpaired_Ttest$p.value
  
}

Fig7J_Mek_compared_groups_ttest <- data.frame((sapply(Fig7J_Mek_compared_groups_ttest,c)))
Fig7J_Mek_compared_groups_ttest$p_value = as.numeric(as.character(Fig7J_Mek_compared_groups_ttest$p_value))
Fig7J_Mek_compared_groups_ttest <- Fig7J_Mek_compared_groups_ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  separate(compared_group, into = c("Analyte1","Cell_line","Treatment"), sep = "_", fill = "right") %>% 
  select(-c(Analyte1)) %>% 
  mutate(Treatment=as.numeric(as.character(Treatment)))

Fig7J_Mek_Replicates = Fig7J_pMek %>% 
  select(Cell_line,Replicate,Treatment,log2Treatment,Analyte,Signal) %>% 
  mutate(Treatment=as.numeric(as.character(Treatment))) %>%
  pivot_wider(names_from = Replicate,values_from=Signal) %>% 
  rowwise() %>% 
  mutate(Mean=mean(c(R1,R2,R3), na.rm = TRUE))

Fig7J_Mek_data = left_join(Fig7J_Mek_Replicates,Fig7J_Mek_compared_groups_ttest, by=c("Treatment","Cell_line")) %>% 
  arrange(desc(Cell_line)) %>% 
  mutate(Analyte= gsub("FCoXX_M2_","",Analyte))

########### Fig 7G : cRaf ###########

Fig7J_pcRaf = PDDR_All_FCoXX_SUBSET %>% 
  filter(Analyte == "FCoXX_M2_cRaf") %>% 
  filter(Treatment != 1000) %>% 
  filter(!(Cell_line == "XX" & (Treatment == 4|Treatment == 12)))

g <- ggplot(Fig7J_pcRaf, aes(x=x_axis_factor,y=Signal)) +
  geom_point(aes(x=x_axis_factor,y=Signal,color=Cell_line))+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 1))+
  ylim(0.5,max(Fig7J_pcRaf$Signal)+0.5)+
  labs(x = "\n cRafi(nM)", 
       y = " Rel. phosp. (norm.) ",
       title = "pcRaf",
       color = "Cell line" )+
  MyPaperTheme

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig7J_pcRaf.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 

# Test for significance of results : Fig7J_cRaf

cRaf_XX_0= Fig7J_pcRaf %>% 
  filter(Cell_line=="XX" & Treatment_Fctr ==0) 

cRaf_XO_0 = Fig7J_pcRaf %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 0) 

cRaf_XO_4 = Fig7J_pcRaf %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 4)

cRaf_XO_12 = Fig7J_pcRaf %>% 
  filter(Cell_line=="XO" & Treatment_Fctr == 12) 


Fig7J_cRaf_compared_groups=list("cRaf_XO_0"=cRaf_XO_0,"cRaf_XO_4"=cRaf_XO_4,"cRaf_XO_12"=cRaf_XO_12)
Fig7J_cRaf_compared_groups_ttest=list()

for(i in 1:length(Fig7J_cRaf_compared_groups)) {
  
  Fig7J_cRaf_compared_groups_ttest$compared_group[i] = names(Fig7J_cRaf_compared_groups[i])
  
  unpaired_Ttest = t.test(cRaf_XX_0$Signal,Fig7J_cRaf_compared_groups[[i]]$Signal, paired = FALSE)
  
  Fig7J_cRaf_compared_groups_ttest$XX_untreated_mean[i] = unpaired_Ttest$estimate[1]
  Fig7J_cRaf_compared_groups_ttest$XO_treated_mean[i] = unpaired_Ttest$estimate[2]
  Fig7J_cRaf_compared_groups_ttest$p_value[i] = unpaired_Ttest$p.value
  
}

Fig7J_cRaf_compared_groups_ttest <- data.frame((sapply(Fig7J_cRaf_compared_groups_ttest,c)))
Fig7J_cRaf_compared_groups_ttest$p_value = as.numeric(as.character(Fig7J_cRaf_compared_groups_ttest$p_value))
Fig7J_cRaf_compared_groups_ttest <- Fig7J_cRaf_compared_groups_ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*","")) %>% 
  separate(compared_group, into = c("Analyte1","Cell_line","Treatment"), sep = "_", fill = "right") %>% 
  select(-c(Analyte1)) %>% 
  mutate(Treatment=as.numeric(as.character(Treatment)))

Fig7J_cRaf_Replicates = Fig7J_pcRaf %>% 
  select(Cell_line,Replicate,Treatment,log2Treatment,Analyte,Signal) %>% 
  mutate(Treatment=as.numeric(as.character(Treatment))) %>%
  pivot_wider(names_from = Replicate,values_from=Signal) %>% 
  rowwise() %>% 
  mutate(Mean=mean(c(R1,R2,R3), na.rm = TRUE))

Fig7J_cRaf_data = left_join(Fig7J_cRaf_Replicates,Fig7J_cRaf_compared_groups_ttest, by=c("Treatment","Cell_line")) %>% 
  arrange(desc(Cell_line)) %>% 
  mutate(Analyte= gsub("FCoXX_M2_","",Analyte))


########### Fig 7G : writing out the data as xls ###########
Fig7J_Mek_cRaf_data = dplyr::bind_rows(Fig7J_Mek_data,Fig7J_cRaf_data)
WriteXLS::WriteXLS(Fig7J_Mek_cRaf_data, "./OUTPUT_PAPER/Fig7J_Mek_cRaf_data.xls")


######################## FIG 7K ###############

Fig7K_pMek = PDDR_All_FCoXX_SUBSET %>% 
  filter(Analyte == "FCoXX_M2_Mek") %>% 
  filter(Treatment != 12) %>% 
  filter((Cell_line =="XX" & (Treatment == 0|Treatment == 1000)) | (Cell_line =="XO" & (Treatment == 4|Treatment == 1000)))



g <- ggplot(Fig7K_pMek, aes(x=x_axis_factor,y=Signal)) +
  geom_point(aes(x=x_axis_factor,y=Signal,color=Cell_line))+
  scale_fill_manual(values = MyCellLineColours) +
  scale_color_manual(values = MyCellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 1))+
  ylim(0,max(Fig7K_pMek$Signal)+0.5)+
  labs(x = "\n MEKi(nM)", 
       y = " Rel. phos. (norm) ",
       title = "pMEK",
       color = "Cell line" )+
  MyPaperTheme

gt=set_panel_size(g,width=unit(2.8,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("Fig7K_Mek.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER") 


Fig7K_ttest = list()
the_cell_lines = c("XX","XO")
for(i in 1:length(the_cell_lines)){
  #i =2
  
  data = Fig7K_pMek %>% 
    filter(Cell_line == the_cell_lines[i])
  
  treated = data %>% filter(Treatment == 1000) 
  base = data %>% filter(Treatment != 1000) 
  
  test_result = t.test(base$Signal,treated$Signal, paired = F)
  
  Fig7K_ttest$Cell_line[i] = the_cell_lines[i]
  Fig7K_ttest$base_conc[i] = unique(base$Treatment)
  Fig7K_ttest$treated_conc[i] = unique(treated$Treatment)
  Fig7K_ttest$base_mean[i] = test_result$estimate[1]
  Fig7K_ttest$treated_mean[i] = test_result$estimate[2]
  Fig7K_ttest$p_value[i] = test_result$p.value
  
}

Fig7K_ttest <- data.frame((sapply(Fig7K_ttest,c)))
Fig7K_ttest$p_value = as.numeric(as.character(Fig7K_ttest$p_value))
Fig7K_ttest <- Fig7K_ttest %>% 
  mutate(sig = ifelse(p_value<0.05,"*",""))

Fig7K_ttest <- Fig7K_ttest %>% 
  pivot_longer(cols = c("base_conc", "treated_conc"), names_to = "Treatment_conc", values_to = "Treatment") %>% 
  mutate(Treatment = as.numeric(as.character(Treatment))) %>% 
  mutate(base_mean = as.numeric(as.character(base_mean))) %>% 
  mutate(treated_mean = as.numeric(as.character(treated_mean))) %>% 
  mutate(p_value = as.numeric(as.character(p_value))) %>% 
  select(-c(Treatment_conc))


Fig7K_Mek_Replicates = Fig7K_pMek %>% 
  select(Cell_line,Replicate,Treatment,log2Treatment,Analyte,Signal) %>% 
  mutate(Treatment=as.numeric(as.character(Treatment))) %>%
  pivot_wider(names_from = Replicate,values_from=Signal) %>% 
  rowwise() %>% 
  mutate(Mean=mean(c(R1,R2,R3), na.rm = TRUE))

Fig7K_Mek_data = left_join(Fig7K_Mek_Replicates,Fig7K_ttest, by=c("Treatment","Cell_line")) %>% 
  arrange(desc(Cell_line)) %>% 
  mutate(Analyte= gsub("FCoXX_M2_","",Analyte)) %>%
  mutate(base_mean = ifelse(Treatment != 1000, NA, base_mean)) %>% 
  mutate(treated_mean = ifelse(Treatment != 1000, NA, treated_mean)) %>% 
  mutate(p_value = ifelse(Treatment != 1000, NA, p_value)) %>% 
  mutate(sig = ifelse(Treatment != 1000, NA, sig)) %>% 
  arrange(Cell_line) %>% 
  rowwise() %>% 
  mutate(FoldChange=treated_mean/base_mean)

WriteXLS::WriteXLS(Fig7K_Mek_data, "./OUTPUT_PAPER/Fig7K_Mek_data.xls")

