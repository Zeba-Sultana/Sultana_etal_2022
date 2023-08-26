library(tidyverse)
library(egg) #needed for set_panel_size

source("../7_ValidationExperiments/ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")


Coeff_Determination <- read.csv("../4_Modeling_STASNet/2_Model_Extension/Coeff_Determination_df.csv")
colnames(Coeff_Determination) <-  gsub("^X$","analyte_name",colnames(Coeff_Determination))

Coeff_Determination$analyte_name <- gsub("Gsk3","pGSK3",Coeff_Determination$analyte_name)
Coeff_Determination$analyte_name <- gsub("Mek","pMEK",Coeff_Determination$analyte_name)
Coeff_Determination$analyte_name <- gsub("mTor","pmTOR",Coeff_Determination$analyte_name)
Coeff_Determination$analyte_name <- gsub("Akt","pAKT",Coeff_Determination$analyte_name)
Coeff_Determination$analyte_name <- gsub("Erk","pERK",Coeff_Determination$analyte_name)
Coeff_Determination$analyte_name <- gsub("Stat3","pSTAT3",Coeff_Determination$analyte_name)
Coeff_Determination$analyte_name <- gsub("Smad2","pSMAD2",Coeff_Determination$analyte_name)

Coeff_Determination$analyte_name <- factor(Coeff_Determination$analyte_name, levels = c("pAKT","pGSK3","pmTOR","pMEK","pERK","pSTAT3","pSMAD2"))

Coeff_Determination_long <- Coeff_Determination %>% 
  pivot_longer(cols= -c(analyte_name),names_to = "Model", values_to = "Coeff")

Coeff_Determination_long$Cell_type <- ifelse(grepl("XX",Coeff_Determination_long$Model), "XX","XO")
Coeff_Determination_long$Model = factor(Coeff_Determination_long$Model, levels = c("Init_XX","Final_XX","Init_XO","Final_XO"), labels = c("Init-XX","Comp-XX","Init-XO","Comp-XO"))

MyCategoryColours <- c("Init-XX"= "#dd1c77","Init-XO"="#f768a1",
                       "Comp-XX"= "#006d2c", "Comp-XO" = "#74c476")

# MyCategoryFill <- c("Init-XX"= "#dd1c77","Init-XO"="#FFFFFF",
#                     "Comp-XX"= "#006d2c", "Comp-XO" = "#FFFFFF")

MyCategoryFill <- c("Init-XX"= "#dd1c77","Init-XO"="#fde0dd",
                    "Comp-XX"= "#006d2c", "Comp-XO" = "#a1d99b")

g <- ggplot(Coeff_Determination_long, aes(x=Model, y=Coeff,color = Model, fill = Model))+
  facet_grid(.~analyte_name)+
  geom_bar(stat = "identity",width=0.6, position = position_dodge(width=0.2),size=0.7)+ #To make thinner outline ;size=0.7 (instead of 1)
  scale_fill_manual(values= MyCategoryFill, labels=c("Initial Model:XX","Completed Model:XX","Initial Model:XO","Completed Model:XO"))+ 
  scale_color_manual(values= MyCategoryColours, labels=c("Initial Model:XX","Completed Model:XX","Initial Model:XO","Completed Model:XO"))+
  ylim(-0.5,1)+
  MyPaperTheme+
  xlab("")+
  ylab(expression(paste("Exp. data vs. model \n coeff. of determination, ", R^2)))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, size = 6, colour = "black",hjust = 1),
        axis.ticks.x = element_blank(),
        
        legend.position ="bottom",
        legend.text = element_text(color = "black", size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(0.25,"cm"),
        
        strip.text.x = element_text(color='black', size=7, face = 'bold'),#set the size and bold face of panel headings of the facet grid.
        
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.spacing.x = unit(0.5, "lines")) #This is to control that space between panels horizontally

gt=set_panel_size(g,width=unit(1.95,'cm'),height=unit(2.0,'cm'))
grid.arrange(gt)
ggsave("Fig3B_Coeff_Determination.pdf", gt, dpi=300, useDingbats=FALSE, path = "./OUTPUT_PAPER",width = 25, height = 20, units = "cm") 


