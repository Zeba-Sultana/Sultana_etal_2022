#!/usr/bin/env Rscript

setwd("./1A_Bioplex_Lxb2MIDAS")
source("Script_1A_AllReplicates_perform_lxb_extraction.R")

setwd("../")
setwd("./1B_WB_SignalQuantification")
source("Script_1B_WB_Data_ReadIn.R")

setwd("../")
setwd("./2A_Bioplex_Norm_FoldChange")
source("Script_2A_MeanNormalization_FCcalculation.R")

setwd("../")
setwd("./2B_WB_Norm_FoldChange")
source("Script_2B_WB_GelMeanNormalization_FCcalculation.R")

setwd("../")
setwd("./3_Merge_Bioplex_WB_data")
source("Script_3_Merge_BP_WB.R")

setwd("../")
setwd("./5_ProfileLikelihood")
source("2_Analyse_ProfileLieklihood_PlotResults.R")

setwd("../")
setwd("./6_SimulationResults_vs_ExpData")
source("Step0_Merging_Bioplex_WBdata_NonMIDAS.R")
source("Step1_Fetch_Exp_Sim_Data.R")
source("Step2_PlotHeatmaps.R")
source("Step3_PointRangePlots_Fig3.R")

setwd("../")
setwd("./7_ValidationExperiments/1_Gsk3i_DR")
source("Script_ValidationFigures_Gsk3iDR.R")

setwd("../../")
setwd("./7_ValidationExperiments/2_Activin_Actri_DR")
source("Script1_Activin_ActrInh_DR_ReadingInData.R")
source("Script2_Activin_ActrInh_DR_ValidationPlot.R")

setwd("../../")
setwd("./7_ValidationExperiments/3_Meki_DR_TC")
source("1_Meki_DR_ValidationFigures.R")
source("2_Meki_TC_ValidationFigures.R")
source("3_ScatterPlot_MekiDR.R")

setwd("../../")


