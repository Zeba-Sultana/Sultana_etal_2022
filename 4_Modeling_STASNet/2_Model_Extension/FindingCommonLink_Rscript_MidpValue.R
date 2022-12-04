#!/usr/bin/env Rscript

library(readxl)
library(tidyverse)
library(htmlTable)


source(./FindingCommonLink_FuncDefinitions.R")

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

########## LA1

LA1_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA1_100k_2020-10-29/extension_XX_R345_Network_LA1_100k+red.csv"))
LA1_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA1_100k_2020-10-29/extension_XO_R345_Network_LA1_100k+red.csv"))

LA1_XX_ext <- filter_ext_list(LA1_XX_ext_InputFile)
LA1_XO_ext <- filter_ext_list(LA1_XO_ext_InputFile)

LA1_common_links_byTot_Percentage_DelResidual <- Select_Link(LA1_XX_ext,LA1_XO_ext,3)
PrintOut_Table(head(LA1_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA1_common_links_byTot_Percentage_DelResidual, "LA1_common_links_byTot_Percentage_DelResidual.csv")

########## LA2

LA2_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA2_100k_2020-10-29/extension_XX_R345_Network_LA2_100k+red.csv"))
LA2_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA2_100k_2020-10-29/extension_XO_R345_Network_LA2_100k+red.csv"))

LA2_XX_ext <- filter_ext_list(LA2_XX_ext_InputFile)
LA2_XO_ext <- filter_ext_list(LA2_XO_ext_InputFile)

LA2_common_links_byTot_Percentage_DelResidual <- Select_Link(LA2_XX_ext,LA2_XO_ext,3)
PrintOut_Table(head(LA2_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA2_common_links_byTot_Percentage_DelResidual, "LA2_common_links_byTot_Percentage_DelResidual.csv")


########## LA3

LA3_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA3_100k_2020-10-29/extension_XX_R345_Network_LA3_100k+red.csv"))
LA3_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA3_100k_2020-10-29/extension_XO_R345_Network_LA3_100k+red.csv"))

LA3_XX_ext <- filter_ext_list(LA3_XX_ext_InputFile)
LA3_XO_ext <- filter_ext_list(LA3_XO_ext_InputFile)

LA3_common_links_byTot_Percentage_DelResidual <- Select_Link(LA3_XX_ext,LA3_XO_ext,3)
PrintOut_Table(head(LA3_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA3_common_links_byTot_Percentage_DelResidual, "LA3_common_links_byTot_Percentage_DelResidual.csv")

########## LA4

LA4_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA4_100k_2020-10-29/extension_XX_R345_Network_LA4_100k+red.csv"))
LA4_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA4_100k_2020-10-29/extension_XO_R345_Network_LA4_100k+red.csv"))

LA4_XX_ext <- filter_ext_list(LA4_XX_ext_InputFile)
LA4_XO_ext <- filter_ext_list(LA4_XO_ext_InputFile)

LA4_common_links_byTot_Percentage_DelResidual <- Select_Link(LA4_XX_ext,LA4_XO_ext,3)
PrintOut_Table(head(LA4_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA4_common_links_byTot_Percentage_DelResidual, "LA4_common_links_byTot_Percentage_DelResidual.csv")

########## LA5

LA5_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA5_100k_2020-10-29/extension_XX_R345_Network_LA5_100k+red.csv"))
LA5_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA5_100k_2020-10-29/extension_XO_R345_Network_LA5_100k+red.csv"))

LA5_XX_ext <- filter_ext_list(LA5_XX_ext_InputFile)
LA5_XO_ext <- filter_ext_list(LA5_XO_ext_InputFile)

LA5_common_links_byTot_Percentage_DelResidual <- Select_Link(LA5_XX_ext,LA5_XO_ext,3)
PrintOut_Table(head(LA5_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA5_common_links_byTot_Percentage_DelResidual, "LA5_common_links_byTot_Percentage_DelResidual.csv")

########## LA6

LA6_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA6_100k_2020-10-29/extension_XX_R345_Network_LA6_100k+red.csv"))
LA6_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA6_100k_2020-10-29/extension_XO_R345_Network_LA6_100k+red.csv"))

LA6_XX_ext <- filter_ext_list(LA6_XX_ext_InputFile)
LA6_XO_ext <- filter_ext_list(LA6_XO_ext_InputFile)

LA6_common_links_byTot_Percentage_DelResidual <- Select_Link(LA6_XX_ext,LA6_XO_ext,3)
PrintOut_Table(head(LA6_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA6_common_links_byTot_Percentage_DelResidual, "LA6_common_links_byTot_Percentage_DelResidual.csv")

########## LA7

LA7_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA7_100k_2020-10-29/extension_XX_R345_Network_LA7_100k+red.csv"))
LA7_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA7_100k_2020-10-29/extension_XO_R345_Network_LA7_100k+red.csv"))

LA7_XX_ext <- filter_ext_list(LA7_XX_ext_InputFile)
LA7_XO_ext <- filter_ext_list(LA7_XO_ext_InputFile)

LA7_common_links_byTot_Percentage_DelResidual <- Select_Link(LA7_XX_ext,LA7_XO_ext,3)
PrintOut_Table(head(LA7_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA7_common_links_byTot_Percentage_DelResidual, "LA7_common_links_byTot_Percentage_DelResidual.csv")

########## LA8

LA8_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA8_100k_2020-10-29/extension_XX_R345_Network_LA8_100k+red.csv"))
LA8_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA8_100k_2020-10-29/extension_XO_R345_Network_LA8_100k+red.csv"))

LA8_XX_ext <- filter_ext_list(LA8_XX_ext_InputFile)
LA8_XO_ext <- filter_ext_list(LA8_XO_ext_InputFile)

LA8_common_links_byTot_Percentage_DelResidual <- Select_Link(LA8_XX_ext,LA8_XO_ext,3)
PrintOut_Table(head(LA8_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA8_common_links_byTot_Percentage_DelResidual, "LA8_common_links_byTot_Percentage_DelResidual.csv")

########## LA9

LA9_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA9_100k_2020-10-30/extension_XX_R345_Network_LA9_100k+red.csv"))
LA9_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA9_100k_2020-10-30/extension_XO_R345_Network_LA9_100k+red.csv"))

LA9_XX_ext <- filter_ext_list(LA9_XX_ext_InputFile)
LA9_XO_ext <- filter_ext_list(LA9_XO_ext_InputFile)

LA9_common_links_byTot_Percentage_DelResidual <- Select_Link(LA9_XX_ext,LA9_XO_ext,3)
PrintOut_Table(head(LA9_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA9_common_links_byTot_Percentage_DelResidual, "LA9_common_links_byTot_Percentage_DelResidual.csv")


########## LA10

LA10_XX_ext_InputFile =read.delim(file.path(Ext_files_dir,"XX/run_XX_R345_Network_LA10_100k_2020-10-30/extension_XX_R345_Network_LA10_100k+red.csv"))
LA10_XO_ext_InputFile =read.delim(file.path(Ext_files_dir,"XO/run_XO_R345_Network_LA10_100k_2020-10-30/extension_XO_R345_Network_LA10_100k+red.csv"))

LA10_XX_ext <- filter_ext_list(LA10_XX_ext_InputFile)
LA10_XO_ext <- filter_ext_list(LA10_XO_ext_InputFile)

LA10_common_links_byTot_Percentage_DelResidual <- Select_Link(LA10_XX_ext,LA10_XO_ext,3)
PrintOut_Table(head(LA10_common_links_byTot_Percentage_DelResidual, n=10))

write.csv(LA10_common_links_byTot_Percentage_DelResidual, "LA10_common_links_byTot_Percentage_DelResidual.csv")

print("All extension links listed")






