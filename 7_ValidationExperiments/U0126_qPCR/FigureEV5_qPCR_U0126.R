library(gdata)
library(tidyverse)
library(egg)


library(readxl)
library('gridExtra')
library(ggpubr)


source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

# Plot Function
CellLineColours = c("1.8_XX"="#FF0000","1.8_XO"="#0000CD", "Pgk_XX"="#fb9a99","Pgk_XO"="#1f78b4", "E14_XY"="black")



dir.create("./OUTPUT_PAPER")
qPCR_folder= "../../RAW_DATA/ValidationExperiments/ActA_U0126_qPCR/"
Results = read.xls(file.path(qPCR_folder,"Results_Compiled_230616_Rinput.xls"))

#colnames(Results)
#Processig the Results file
Results_filt = Results %>% 
  select(grep("Oct4|Nanog|Fgf5|Otx2|Prdm14|Dnmt3b", colnames(Results))) 
colnames(Results_filt) = gsub("Oct4_SampleName", "Samples",colnames(Results_filt) )
Results_filt$Samples = gsub("PgK", "Pgk",Results_filt$Samples) 
Results_filt = Results_filt %>% 
  select(!grep("_SampleName",colnames(Results_filt))) 
Results_filt = Results_filt %>% 
  select(c(Samples,grep("Power2",colnames(Results_filt))))


Results_filt = Results_filt %>% 
  separate(Samples, into = c("Rep", "cell_line", "x_status", "treatment_done", "time"), sep = "_") %>% 
  unite(Treatment, treatment_done, time, sep = "_") %>% 
  unite(CellLine, cell_line, x_status, sep = "_", remove = FALSE)

Results_filt = Results_filt %>% 
  unite(Sample, cell_line,x_status,Treatment, sep = "_", remove = FALSE)

colnames(Results_filt) = gsub("_Power2_del_CT", "",colnames(Results_filt))

Results_All = distinct(Results_filt)

########## Results U0126 ########## 

############ ############ ############  
############     UO126    ############ 
############ ############ ############

Results_All_UO126 = Results_All %>% 
  filter(grepl("UO126|DMSO", Treatment)) %>% 
  mutate(Sample=paste(CellLine,Treatment, sep="_"))

Results_All_UO126 = Results_All_UO126 %>% 
  separate(CellLine, into = c("cell_line", "x_status"), sep = "_", remove = FALSE)

Results_All_UO126_long =  Results_All_UO126 %>% 
  pivot_longer(cols=-c(Sample,Rep,CellLine,cell_line, x_status, Treatment), names_to = "gene", values_to = "Rawvalue")

Sample_order <- c("1.8_XX_DMSO1_24h", 
                  "1.8_XX_UO126_24h",
                  "1.8_XO_DMSO1_24h", 
                  "1.8_XO_UO126_24h", 
                  "Pgk_XX_DMSO1_24h", 
                  "Pgk_XX_UO126_24h",
                  "Pgk_XO_DMSO1_24h", 
                  "Pgk_XO_UO126_24h",
                  "E14_XY_DMSO1_24h", 
                  "E14_XY_UO126_24h" )
Sample_labels <- c("1.8_XX_0", 
                   "1.8_XX_UO126",
                   "1.8_XO_0", 
                   "1.8_XO_UO126", 
                   "Pgk_XX_0", 
                   "Pgk_XX_UO126",
                   "Pgk_XO_0", 
                   "Pgk_XO_UO126",
                   "E14_XY_0", 
                   "E14_XY_UO126" )
Results_All_UO126_long$Sample <- factor(Results_All_UO126_long$Sample, levels = Sample_order, labels = Sample_labels)


# Treatment_order <- c("DMSO1_24h", "UO126_24h")
# Treatment_labels <- c("DMSO (24h)", "UO126 (24h)")
# Results_All_UO126_long$Treatment <- factor(Results_All_UO126_long$Treatment, levels = Treatment_order, labels = Treatment_labels)
# 
# Genes_order <- c("Oct4","Nanog","Fgf5","Otx2","Prdm14","Dnmt3b")
# Genes_labels <- c("Oct4","Nanog", "Fgf5","Otx2","Prdm14","Dnmt3b")
# Results_All_UO126_long$gene <- factor(Results_All_UO126_long$gene, levels = Genes_order, labels = Genes_labels)
# 
# cell_line_order <- c("1.8", "Pgk", "E14")
# Results_All_UO126_long$cell_line <- factor(Results_All_UO126_long$cell_line, levels = cell_line_order)
# 
# x_status_order <- c("XX", "XO")
# Results_All_UO126_long$x_status <- factor(Results_All_UO126_long$x_status, levels = x_status_order)

################

FigEV5_U0126_Fgf5 = Results_All_UO126_long %>% 
  filter(gene == "Fgf5")

FigEV5_U0126_Fgf5 = FigEV5_U0126_Fgf5 %>% 
  group_by(CellLine) %>% 
  mutate(mean_cntrl= mean(Rawvalue[grepl("DMSO",Treatment)])) %>% 
  mutate(FCoCntrl = Rawvalue/(mean_cntrl+0.00001))

FigEV5_U0126_Fgf5_data = FigEV5_U0126_Fgf5 %>% 
select(-c(cell_line,x_status))
WriteXLS::WriteXLS(FigEV5_U0126_Fgf5_data, "./OUTPUT_PAPER/FigEV5_U0126_Fgf5_data.xls")

FigEV5_U0126_Fgf5 = FigEV5_U0126_Fgf5  %>% 
  select(-c(mean_cntrl))

FigEV5_U0126_Fgf5_l = FigEV5_U0126_Fgf5 %>% 
  pivot_longer(-c(Sample, Rep, CellLine, cell_line, x_status, Treatment, gene),  names_to = "Raw_FC", values_to = "yvalues")

FigEV5_U0126_Fgf5_l$Raw_FC = factor(FigEV5_U0126_Fgf5_l$Raw_FC, levels = c("Rawvalue", "FCoCntrl"), labels = c("Rel. expn.", "Fold change over untreated"))

g = ggplot(FigEV5_U0126_Fgf5_l, aes(Sample, yvalues)) +
  geom_point(aes(shape=Rep, color=CellLine, group = x_status),position = position_dodge(width = 0.75), size = 2) +
  scale_fill_manual(values = CellLineColours) +
  scale_color_manual(values = CellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 0.5))+
  facet_grid(Raw_FC~., scales = "free_y", switch = 'y')+
  labs(y = expression(paste("2^delCT")), x="")+
  labs(x = "\n U0126(uM)", 
       y = "",
       title = "Fgf5" )+
  MyPaperTheme+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5, size=8),
        strip.placement = "outside",
        strip.text.y = element_text(size = 8, colour = "black", angle = 90,face = "plain" ))
gt=egg::set_panel_size(g,width=unit(5.6,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT_PAPER/FigEV5_U0126_Fgf5.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

########## Otx2 ##########

FigEV5_U0126_Otx2 = Results_All_UO126_long %>% 
  filter(gene == "Otx2")

FigEV5_U0126_Otx2 = FigEV5_U0126_Otx2 %>% 
  group_by(CellLine) %>% 
  mutate(mean_cntrl= mean(Rawvalue[grepl("DMSO",Treatment)])) %>% 
  mutate(FCoCntrl = Rawvalue/(mean_cntrl+0.00001))

FigEV5_U0126_Otx2_data = FigEV5_U0126_Otx2 %>% 
  select(-c(cell_line,x_status))
WriteXLS::WriteXLS(FigEV5_U0126_Otx2_data, "./OUTPUT_PAPER/FigEV5_U0126_Otx2_data.xls")

FigEV5_U0126_Otx2 = FigEV5_U0126_Otx2 %>% 
  select(-c(mean_cntrl))

FigEV5_U0126_Otx2_l = FigEV5_U0126_Otx2 %>% 
  pivot_longer(-c(Sample, Rep, CellLine, cell_line, x_status, Treatment, gene),  names_to = "Raw_FC", values_to = "yvalues")

FigEV5_U0126_Otx2_l$Raw_FC = factor(FigEV5_U0126_Otx2_l$Raw_FC, levels = c("Rawvalue", "FCoCntrl"), labels = c("Rel. expn.", "Fold change over untreated"))

g = ggplot(FigEV5_U0126_Otx2_l, aes(Sample, yvalues)) +
  geom_point(aes(shape=Rep, color=CellLine, group = x_status),position = position_dodge(width = 0.75), size = 2) +
  scale_fill_manual(values = CellLineColours) +
  scale_color_manual(values = CellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 0.5))+
  facet_grid(Raw_FC~., scales = "free_y", switch = 'y')+
  labs(y = expression(paste("2^delCT")), x="")+
  labs(x = "\n U0126(uM)", 
       y = "",
       title = "Otx2" )+
  MyPaperTheme+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5, size=8),
        strip.placement = "outside",
        strip.text.y = element_text(size = 8, colour = "black", angle = 90,face = "plain" ))
gt=egg::set_panel_size(g,width=unit(5.6,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT_PAPER/FigEV5_U0126_Otx2.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 


########## Oct4 ##########

FigEV5_U0126_Oct4 = Results_All_UO126_long %>% 
  filter(gene == "Oct4")

FigEV5_U0126_Oct4 = FigEV5_U0126_Oct4 %>% 
  group_by(CellLine) %>% 
  mutate(mean_cntrl= mean(Rawvalue[grepl("DMSO",Treatment)])) %>% 
  mutate(FCoCntrl = Rawvalue/(mean_cntrl+0.00001))

FigEV5_U0126_Oct4_data = FigEV5_U0126_Oct4 %>% 
  select(-c(cell_line,x_status))
WriteXLS::WriteXLS(FigEV5_U0126_Oct4_data, "./OUTPUT_PAPER/FigEV5_U0126_Oct4_data.xls")

FigEV5_U0126_Oct4 = FigEV5_U0126_Oct4 %>% 
  select(-c(mean_cntrl))

FigEV5_U0126_Oct4_l = FigEV5_U0126_Oct4 %>% 
  pivot_longer(-c(Sample, Rep, CellLine, cell_line, x_status, Treatment, gene),  names_to = "Raw_FC", values_to = "yvalues")

FigEV5_U0126_Oct4_l$Raw_FC = factor(FigEV5_U0126_Oct4_l$Raw_FC, levels = c("Rawvalue", "FCoCntrl"), labels = c("Rel. expn.", "Fold change over untreated"))

g = ggplot(FigEV5_U0126_Oct4_l, aes(Sample, yvalues)) +
  geom_point(aes(shape=Rep, color=CellLine, group = x_status),position = position_dodge(width = 0.75), size = 2) +
  scale_fill_manual(values = CellLineColours) +
  scale_color_manual(values = CellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 0.5))+
  facet_grid(Raw_FC~., scales = "free_y", switch = 'y')+
  labs(y = expression(paste("2^delCT")), x="")+
  labs(x = "\n U0126(uM)", 
       y = "",
       title = "Oct4" )+
  MyPaperTheme+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5, size=8),
        strip.placement = "outside",
        strip.text.y = element_text(size = 8, colour = "black", angle = 90,face = "plain" ))
gt=egg::set_panel_size(g,width=unit(5.6,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT_PAPER/FigEV5_U0126_Oct4.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 


########## Nanog ##########

FigEV5_U0126_Nanog = Results_All_UO126_long %>% 
  filter(gene == "Nanog")

FigEV5_U0126_Nanog = FigEV5_U0126_Nanog %>% 
  group_by(CellLine) %>% 
  mutate(mean_cntrl= mean(Rawvalue[grepl("DMSO",Treatment)])) %>% 
  mutate(FCoCntrl = Rawvalue/(mean_cntrl+0.00001))

FigEV5_U0126_Nanog_data = FigEV5_U0126_Nanog %>% 
  select(-c(cell_line,x_status))
WriteXLS::WriteXLS(FigEV5_U0126_Nanog_data, "./OUTPUT_PAPER/FigEV5_U0126_Nanog_data.xls")

FigEV5_U0126_Nanog = FigEV5_U0126_Nanog %>% 
  select(-c(mean_cntrl))

FigEV5_U0126_Nanog_l = FigEV5_U0126_Nanog %>% 
  pivot_longer(-c(Sample, Rep, CellLine, cell_line, x_status, Treatment, gene),  names_to = "Raw_FC", values_to = "yvalues")

FigEV5_U0126_Nanog_l$Raw_FC = factor(FigEV5_U0126_Nanog_l$Raw_FC, levels = c("Rawvalue", "FCoCntrl"), labels = c("Rel. expn.", "Fold change over untreated"))

g = ggplot(FigEV5_U0126_Nanog_l, aes(Sample, yvalues)) +
  geom_point(aes(shape=Rep, color=CellLine, group = x_status),position = position_dodge(width = 0.75), size = 2) +
  scale_fill_manual(values = CellLineColours) +
  scale_color_manual(values = CellLineColours) +
  stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 0.5))+
  facet_grid(Raw_FC~., scales = "free_y", switch = 'y')+
  labs(y = expression(paste("2^delCT")), x="")+
  labs(x = "\n U0126(uM)", 
       y = "",
       title = "Nanog" )+
  MyPaperTheme+
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5, size=8),
        strip.placement = "outside",
        strip.text.y = element_text(size = 8, colour = "black", angle = 90,face = "plain" ))
gt=egg::set_panel_size(g,width=unit(5.6,'cm'),height=unit(2.8,'cm'))
grid.arrange(gt)
ggsave("./OUTPUT_PAPER/FigEV5_U0126_Nanog.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 


######## Statistical Tests : ANOVA ######

## Fgf5
data = read_excel('./OUTPUT_PAPER/FigEV5_U0126_Fgf5_data.xls') %>%
  separate(CellLine,into=c('line','geno'),sep='_') %>%
  mutate (nr.x = ifelse(geno=='XX',2,1))

two.way <- aov(Rawvalue ~ Treatment * nr.x, data = data)
summary(two.way)

## Otx2
data = read_excel('./OUTPUT_PAPER/FigEV5_U0126_Otx2_data.xls') %>%
  separate(CellLine,into=c('line','geno'),sep='_') %>%
  mutate (nr.x = ifelse(geno=='XX',2,1))

two.way <- aov(Rawvalue ~ Treatment * nr.x, data = data)
summary(two.way)


## Oct4
data = read_excel('./OUTPUT_PAPER/FigEV5_U0126_Oct4_data.xls') %>%
  separate(CellLine,into=c('line','geno'),sep='_') %>%
  mutate (nr.x = ifelse(geno=='XX',2,1))

two.way <- aov(Rawvalue ~ Treatment * nr.x, data = data)
summary(two.way)


## Nanog
data = read_excel('./OUTPUT_PAPER/FigEV5_U0126_Nanog_data.xls') %>%
  separate(CellLine,into=c('line','geno'),sep='_') %>%
  mutate (nr.x = ifelse(geno=='XX',2,1))

two.way <- aov(Rawvalue ~ Treatment * nr.x, data = data)
summary(two.way)



