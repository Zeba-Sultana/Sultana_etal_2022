#!/usr/bin/env Rscript

library(ggplot2)
library(egg)
library(tidyverse)

### For Setting the Paper Theme
source("/project/ag_schulz/Zeba/SCRIPTS_Sultana_etal_2021/PS10_ValidationExp_CAPS/ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")
#source("/Users/sultana/PHD_PROJECT_Zeba/Zeba_PhD_PAPER/SCRIPTS_Sultana_etal_2021/PS10_ValidationExp_CAPS/ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R")


### ONLY FINAL MODEL ###
Dir_p_005 = "../3_Model_Extension_Different_Thresholds"
Residuals_file_p_005 <- "Links_Residuals_p_005.csv"
Residuals_p_005 = read.csv(file.path(Dir_p_005,Residuals_file_p_005),header = TRUE)
Residuals_p_005_Links <- Residuals_p_005 %>% 
  select(-c(X))

# Model_Set_Residuals <- c("Model_iteration"= 11,"res_p_005_XX" =909,"res_p_005_XO" =443)
# rbind(Residuals_p_005_Links,Model_Set_Residuals) 

Model_iteration=c("Init","1","2","3","4","5","6","7","8", "9","10")
Residuals_p_005_Links$Model_iteration = factor(Residuals_p_005_Links$Model_iteration, levels =Model_iteration )

Seq_Links=c("Init Network","Jak->Fgfr","Gsk3->Igfr","Erk->Raf","Bmp4r->Fgfr","Bmp4r->Smad2","Jak->Gsk3","Actr->Ras","Lifr->Akt", "Gsk3->Lifr","Jak->Bmp4r")
Label_Links=c("Initial Network","Jak -> Fgfr","Gsk3 -> Igfr","Erk -> Raf","Bmp4r -> Fgfr","Bmp4r -> Smad2","Jak -> Gsk3","Actr -> Ras","Lifr -> Akt", "Gsk3 -> Lifr","Jak -> Bmp4r")

Residuals_p_005_Links$Links = factor(Residuals_p_005_Links$Links, levels = Seq_Links, labels = Label_Links )

Residuals_p_005_Links_long <- Residuals_p_005_Links %>% 
  pivot_longer(grep("XX|XO",colnames(Residuals_p_005_Links)),names_to = "Model", values_to = "residual")


MyModelColours <-     c("XX" = "#fc9272", 
                        "res_p_05_XX" = "#fb6a4a", 
                        "res_p_01_XX" = "#ef3b2c",
                        "res_p_005_XX" = "#FF0000", # The XX colour
                        "res_p_001_XX" = "#a50f15",
                        
                        "XO" = "#9ecae1",
                        "res_p_05_XO" = "#6baed6",
                        "res_p_01_XO" = "#4292c6",
                        "res_p_005_XO" = "#0000CD",  # The XO colour
                        "res_p_001_XO" = "#08519c")

p <- ggplot(Residuals_p_005_Links_long, aes(x=Links, y=residual, color=Model, group=Model))+
  MyPaperTheme+
  geom_point(size=2)+
  geom_line(size=1)+
  scale_color_manual(values = MyModelColours, 
                     breaks=c("XX", "res_p_05_XX", "res_p_01_XX","res_p_005_XX","res_p_001_XX",
                              "XO", "res_p_05_XO", "res_p_01_XO","res_p_005_XO","res_p_001_XO"),
                     labels=c("XX separate", "XX : p < 0.05", "XX : p < 0.01","XX : p < 0.005","XX : p < 0.001",
                              "XO separate", "XO : p < 0.05", "XO : p < 0.01","XO : p < 0.005","XO : p < 0.001")) +
  labs(x = "Added Links",
       y = "Model Residual\n")+
  ggtitle("Network Extension")+
  ylim(0,NA)+
  geom_hline(yintercept=371, linetype="dashed", color = "black")+
  theme(text=element_text(size = 8,color = "black"),
        legend.position="right",
        legend.justification = c("centre", "centre"),
        legend.title = element_blank(),
       # axis.title.x = element_text(size = 10),
       # axis.title.y = element_text(size = 10),
        panel.border = element_rect(fill=NA, colour = "black", size=0.1), #bump 
        axis.text.y = element_text(size = 7,color = "black", face = "plain"),
        axis.text.x = element_text(size = 7,color = "black", face = "plain",angle=45, vjust = 1, hjust = 1)) #vjust : po value moves label up #hjust: pos value moves label left

print(p)


gt=set_panel_size(p,height=unit(3.5,'cm'),width=unit(4.5,'cm'))
grid.arrange(gt)
ggsave("Fig2B_ModelResiduals_p005.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./")



################  p < 0.05 #########

Dir_p_05 = "../3_Model_Extension_Different_Thresholds"

Residuals_file_p_05 <- "Links_Residuals_p_05.csv"
Residuals_p_05 = read.csv(file.path(Dir_p_05,Residuals_file_p_05),header = TRUE)
Residuals_p_05_WOLinks <- Residuals_p_05 %>% 
  select(-c(X,Links))



################ p < 0.01 #########
Dir_p_01 = "../3_Model_Extension_Different_Thresholds"
#Dir_p_01 = "/Users/sultana/PHD_PROJECT_Zeba/Zeba_PhD_PAPER/SCRIPTS_Sultana_etal_2021/PS4_STASNet_Modeling"
Residuals_file_p_01 <- "Links_Residuals_p_01.csv"
Residuals_p_01 = read.csv(file.path(Dir_p_01,Residuals_file_p_01),header = TRUE)
Residuals_p_01_WOLinks <- Residuals_p_01 %>% 
  select(-c(X,Links))


################ p < 0.005 #########
Dir_p_005 = "../3_Model_Extension_Different_Thresholds"
#Dir_p_005 = "/Users/sultana/PHD_PROJECT_Zeba/Zeba_PhD_PAPER/SCRIPTS_Sultana_etal_2021/PS4_STASNet_Modeling"
Residuals_file_p_005 <- "Links_Residuals_p_005.csv"
Residuals_p_005 = read.csv(file.path(Dir_p_005,Residuals_file_p_005),header = TRUE)
Residuals_p_005_Links <- Residuals_p_005 %>% 
  select(-c(X,Links))



################ p < 0.001 #########
Dir_p_001 = "../3_Model_Extension_Different_Thresholds"
#Dir_p_001 = "/Users/sultana/PHD_PROJECT_Zeba/Zeba_PhD_PAPER/SCRIPTS_Sultana_etal_2021/PS4_STASNet_Modeling"
Residuals_file_p_001 <- "Links_Residuals_p_001.csv"
Residuals_p_001 = read.csv(file.path(Dir_p_001,Residuals_file_p_001),header = TRUE)
Residuals_p_001_WOLinks <- Residuals_p_001 %>% 
  select(-c(X,Links)) 
Residuals_p_001_WOLinks[nrow(Residuals_p_001_WOLinks)+1,] <- NA
Residuals_p_001_WOLinks[nrow(Residuals_p_001_WOLinks)+1,] <- NA
Residuals_p_001_WOLinks[nrow(Residuals_p_001_WOLinks)+1,] <- NA

############ XX and XO extended separately ########

SepXXXO_residuals <- read.csv("/Users/sultana/PHD_PROJECT_Zeba/ANALYSIS_WORKFLOW_Jan2019/S4_Modeling/All_residuals.csv")
SepXXXO_residuals_WOlinks <- SepXXXO_residuals %>% 
  select(XX,XO)


################ All residuals #########
Residuals_all <-  Residuals_p_05_WOLinks %>% 
  bind_cols(Residuals_p_01_WOLinks %>% select(-Model_iteration)) %>% 
  bind_cols(Residuals_p_005_Links %>% select(-Model_iteration)) %>% 
  bind_cols(Residuals_p_001_WOLinks %>% select(-Model_iteration)) %>% 
  bind_cols(SepXXXO_residuals_WOlinks)


########### Plottimg Residuals ##########

Model_iteration=c("Init","1","2","3","4","5","6","7","8", "9","10")
Residuals_all$Model_iteration = factor(Residuals_all$Model_iteration, levels =Model_iteration )

Residuals_all_long <- Residuals_all %>% 
  pivot_longer(grep("XX|XO",colnames(Residuals_all)),names_to = "Model", values_to = "residual")


MyModelColours <-     c("XX" = "#fc9272", 
                       "res_p_05_XX" = "#fb6a4a", 
                       "res_p_01_XX" = "#ef3b2c",
                       "res_p_005_XX" = "#FF0000", # The XX colour
                       "res_p_001_XX" = "#a50f15",
                       
                       "XO" = "#9ecae1",
                       "res_p_05_XO" = "#6baed6",
                       "res_p_01_XO" = "#4292c6",
                       "res_p_005_XO" = "#0000CD",  # The XO colour
                       "res_p_001_XO" = "#08519c")

p <- ggplot(Residuals_all_long, aes(x=Model_iteration, y=residual, color=Model, group=Model))+
  MyPaperTheme+
  geom_point(size=1)+
  geom_line(size=0.5)+
  #scale_color_manual(values = MyModelColours)+
  scale_color_manual(values = MyModelColours, 
                     breaks=c("XX", "res_p_05_XX", "res_p_01_XX","res_p_005_XX","res_p_001_XX",
                              "XO", "res_p_05_XO", "res_p_01_XO","res_p_005_XO","res_p_001_XO"),
                     labels=c("XX separate", "XX : p < 0.05", "XX : p < 0.01","XX : p < 0.005","XX : p < 0.001",
                              "XO separate", "XO : p < 0.05", "XO : p < 0.01","XO : p < 0.005","XO : p < 0.001")) +
  labs(x = "Iteration number",
       y = "Model Residual\n")+
  ggtitle("Network Extension")+
  ylim(0,NA)+
  geom_hline(yintercept=371, linetype="dashed", color = "black")+
  theme(text=element_text(size = 8,color = "black"),
        legend.position="right",
        legend.justification = c("left", "centre"),
        legend.title = element_blank(),
        # axis.title.x = element_text(size = 10),
        # axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 7,color = "black", face = "plain"),
        axis.text.x = element_text(size = 7,color = "black", face = "plain")) #vjust : po value moves label up #hjust: pos value moves label left


print(p)

gt=set_panel_size(p,height=unit(3.5,'cm'),width=unit(4.5,'cm'))
grid.arrange(gt)
ggsave("FigS2C_ModelResiduals_all.pdf", gt, dpi=300, useDingbats=FALSE ,path = "./")



