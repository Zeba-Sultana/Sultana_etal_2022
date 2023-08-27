library(gdata)
library(tidyverse)
library(egg)
library(readxl)


source("../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R")
source("../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")

# Plot Function
CellLineColours = c("1.8_XX"="#FF0000","1.8_XO"="#0000CD", "Pgk_XX"="#fb9a99","Pgk_XO"="#1f78b4", "E14_XY"="black")
plot_qPCR_Activin = function(cell_line_data,cell_line,gene_name,ymax,panel_width=2.8){
  
  plotdata = cell_line_data %>% 
    filter(gene == gene_name)
  
  g = ggplot(plotdata, aes(Sample, value)) +
    geom_point(aes(shape=Rep, color=CellLine, group = x_status),position = position_dodge(width = 0.75), size = 2) +
    scale_fill_manual(values = CellLineColours) +
    scale_color_manual(values = CellLineColours) +
    stat_summary(fun = mean,geom = "crossbar", size = 0.1, position = position_dodge(width = 0.02), aes(width=0.5))+
    ylim(0,ymax)+
    labs(x = "", 
         y = paste0(gene_name," (rel. expr.) "),
         title = cell_line )+
    MyPaperTheme+
    theme(axis.text.x = element_text(angle = 90, hjust =1, vjust =0.5, size=7))
  gt=egg::set_panel_size(g,width=unit(panel_width,'cm'),height=unit(2.8,'cm'))
  grid.arrange(gt)
  return(gt)
}


dir.create("./OUTPUT_PAPER")
qPCR_folder= "../../RAW_DATA/ValidationExperiments/ActA_U0126_qPCR/"
Results = read.xls(file.path(qPCR_folder,"Results_Compiled_230616_Rinput.xls"))

#colnames(Results)
#Processig the Results file
Results_filt = Results %>% 
  select(grep("Oct4|Nanog|Fgf5|Otx2", colnames(Results))) 
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

########## Results ActA ########## 

############ ############ ############  
############     ActA    ############ 
############ ############ ############

## Activin ##
Results_All_Activin = Results_All %>% 
  filter(grepl("ActA|untreated", Treatment)) %>% 
  filter(!grepl("LIF", Treatment))

Results_All_Activin_long =  Results_All_Activin %>% 
  pivot_longer(cols=-c(Rep,Sample,CellLine,cell_line,x_status,Treatment), names_to = "gene", values_to = "value")
Results_All_Activin_long$value = round(Results_All_Activin_long$value,digits = 6)
Results_All_Activin_long = distinct(Results_All_Activin_long)

############ Saving source data ###########
Fig6F_ActA_Fgf5_qPCR = Results_All_Activin_long %>% 
  filter(gene == "Fgf5" ) %>% 
  select(-c(cell_line,x_status))
WriteXLS::WriteXLS(Fig6F_ActA_Fgf5_qPCR, "./OUTPUT_PAPER/Fig6F_ActA_Fgf5_qPCR.xls")

FigEV4A_ActA_Otx2_qPCR = Results_All_Activin_long %>% 
  filter(gene == "Otx2" ) %>% 
  select(-c(cell_line,x_status))
WriteXLS::WriteXLS(FigEV4A_ActA_Otx2_qPCR, "./OUTPUT_PAPER/FigEV4A_ActA_Otx2_qPCR.xls")

FigEV4B_ActA_Oct4_qPCR = Results_All_Activin_long %>% 
  filter(gene == "Oct4" ) %>% 
  select(-c(cell_line,x_status))
WriteXLS::WriteXLS(FigEV4B_ActA_Oct4_qPCR, "./OUTPUT_PAPER/FigEV4B_ActA_Oct4_qPCR.xls")

FigEV4C_ActA_Nanog_qPCR = Results_All_Activin_long %>% 
  filter(gene == "Nanog" ) %>% 
  select(-c(cell_line,x_status))
WriteXLS::WriteXLS(FigEV4C_ActA_Nanog_qPCR, "./OUTPUT_PAPER/FigEV4C_ActA_Nanog_qPCR.xls")

################ 

Activin_1.8_long = Results_All_Activin_long %>% 
  filter(cell_line == "1.8")

Activin_Pgk_long = Results_All_Activin_long %>% 
  filter(cell_line == "Pgk")

Activin_E14_long = Results_All_Activin_long %>% 
  filter(cell_line == "E14")


###### Plots : 1.8 XX and XO ####
Samples_order_1.8 = c("1.8_XX_untreated_48h", "1.8_XX_ActA_24h", "1.8_XX_ActA_48h",
                      "1.8_XO_untreated_48h", "1.8_XO_ActA_24h", "1.8_XO_ActA_48h")
Samples_labels_1.8 = c("-", "ActA (24h)", "ActA (48h)",
                       "- XO", "ActA (24h) XO", "ActA (48h) XO")

Activin_1.8_long$Sample <- factor(Activin_1.8_long$Sample, levels = Samples_order_1.8, labels = Samples_labels_1.8)

gt=plot_qPCR_Activin(Activin_1.8_long, "1.8","Fgf5", 0.07)
ggsave("./OUTPUT_PAPER/Fig6F_ActA_Fgf5_1.8.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

gt=plot_qPCR_Activin(Activin_1.8_long,"1.8", "Otx2", 0.25)
ggsave("./OUTPUT_PAPER/FigEV4_ActA_Otx2_1.8.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

gt=plot_qPCR_Activin(Activin_1.8_long,"1.8", "Oct4",17)
ggsave("./OUTPUT_PAPER/FigEV4_ActA_Oct4_1.8.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

gt=plot_qPCR_Activin(Activin_1.8_long, "1.8","Nanog",3)
ggsave("./OUTPUT_PAPER/FigEV4_ActA_Nanog_1.8.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

###### Plots : Pgk XX and XO ####

Samples_order_Pgk = c("Pgk_XX_untreated_48h", "Pgk_XX_ActA_24h", "Pgk_XX_ActA_48h",
                      "Pgk_XO_untreated_48h", "Pgk_XO_ActA_24h", "Pgk_XO_ActA_48h")
Samples_labels_Pgk = c("-", "ActA (24h)", "ActA (48h)",
                       "- XO", "ActA (24h) XO", "ActA (48h) XO")
Activin_Pgk_long$Sample <- factor(Activin_Pgk_long$Sample, levels = Samples_order_Pgk, labels = Samples_labels_Pgk)

gt=plot_qPCR_Activin(Activin_Pgk_long,"Pgk","Fgf5", 0.07)
ggsave("./OUTPUT_PAPER/Fig6F_ActA_Fgf5_Pgk.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

gt=plot_qPCR_Activin(Activin_Pgk_long,"Pgk", "Otx2", 0.25)
ggsave("./OUTPUT_PAPER/FigEV4_ActA_Otx2_Pgk.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

gt=plot_qPCR_Activin(Activin_Pgk_long,"Pgk", "Oct4",17)
ggsave("./OUTPUT_PAPER/FigEV4_ActA_Oct4_Pgk.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

gt=plot_qPCR_Activin(Activin_Pgk_long,"Pgk", "Nanog",3)
ggsave("./OUTPUT_PAPER/FigEV4_ActA_Nanog_Pgk.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./")

###### Plots : E14 XY ####

Samples_order_E14 = c("E14_XY_untreated_48h", "E14_XY_ActA_24h", "E14_XY_ActA_48h")
Samples_labels_E14 = c("-", "ActA (24h)", "ActA (48h)")
Activin_E14_long$Sample <- factor(Activin_E14_long$Sample, levels = Samples_order_E14, labels = Samples_labels_E14)

gt=plot_qPCR_Activin(Activin_E14_long,"E14","Fgf5", 0.07, 1.4)
ggsave("./OUTPUT_PAPER/Fig6F_ActA_Fgf5_E14.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

gt=plot_qPCR_Activin(Activin_E14_long,"E14", "Otx2", 0.25, 1.4)
ggsave("./OUTPUT_PAPER/FigEV4_ActA_Otx2_E14.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

gt=plot_qPCR_Activin(Activin_E14_long,"E14", "Oct4",17, 1.4)
ggsave("./OUTPUT_PAPER/FigEV4_ActA_Oct4_E14.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 

gt=plot_qPCR_Activin(Activin_E14_long,"E14", "Nanog",3, 1.4)
ggsave("./OUTPUT_PAPER/FigEV4_ActA_Nanog_E14.pdf", gt, dpi=300,height=6, useDingbats=FALSE ,path = "./") 


########### Statistical Test : ANOVA ###########
### Fig. 6  + EV4
## Fgf5
data = read_excel('./OUTPUT_PAPER/Fig6F_ActA_Fgf5_qPCR.xls') %>%
  separate(CellLine,into=c('line','geno'),sep='_') %>%
  mutate (nr.x = ifelse(geno=='XX',2,1))

two.way <- aov(value ~ Treatment * nr.x, data = data)
summary(two.way)

#Otx2
data = read_excel('./OUTPUT_PAPER/FigEV4A_ActA_Otx2_qPCR.xls') %>%
  separate(CellLine,into=c('line','geno'),sep='_') %>%
  mutate (nr.x = ifelse(geno=='XX',2,1))

two.way <- aov(value ~ Treatment * nr.x, data = data)
summary(two.way)

#Oct4
data = read_excel('./OUTPUT_PAPER/FigEV4B_ActA_Oct4_qPCR.xls') %>%
  separate(CellLine,into=c('line','geno'),sep='_') %>%
  mutate (nr.x = ifelse(geno=='XX',2,1))

two.way <- aov(value ~ Treatment * nr.x, data = data)
summary(two.way)

#Nanog
data = read_excel('./OUTPUT_PAPER/FigEV4C_ActA_Nanog_qPCR.xls') %>%
  separate(CellLine,into=c('line','geno'),sep='_') %>%
  mutate (nr.x = ifelse(geno=='XX',2,1))

two.way <- aov(value ~ Treatment * nr.x, data = data)
summary(two.way)

