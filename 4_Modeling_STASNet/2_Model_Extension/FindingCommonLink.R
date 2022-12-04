

#!/usr/bin/env Rscript

library(readxl)
library(tidyverse)
library(htmlTable)


source("./FindingCommonLink_FuncDefinitions.R")
################# INPUT HERE ###############
#Execute from the folder that carries the outputs from XX and XO modeling
Ext_files_dir ="."
######### INIT
INIT_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_100k_2020-10-29/extension_XX_R345_Network_100k+red.csv"))
INIT_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_100k_2020-10-29/extension_XO_R345_Network_100k+red.csv"))
INIT_XX_ext <- filter_ext_list(INIT_XX_ext_InputFile)
INIT_XO_ext <- filter_ext_list(INIT_XO_ext_InputFile)
INIT_common_links_byTot_Percentage_DelResidual <- Select_Link(INIT_XX_ext,INIT_XO_ext,3)
PrintOut_Table(head(INIT_common_links_byTot_Percentage_DelResidual, n=10))
write.csv(INIT_common_links_byTot_Percentage_DelResidual, "INIT_common_links_byTot_Percentage_DelResidual.csv")
