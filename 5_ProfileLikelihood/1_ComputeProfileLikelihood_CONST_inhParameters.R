#!/usr/bin/Rscript
# -*- coding:utf-8 -*-


library("STASNet")
library("knitHelpersSTASNet") #This is needed for checkProfileLikelihood. 


nb_steps = 1000
# Autodetection of the cores
cores = 80


dir.create("OUTPUT")


XX_LA10_Folder = "./INPUT/XX_Network_Completed/XX_MRA"
XO_LA10_Folder = "./INPUT/XO_Network_Completed/XO_MRA"

XX_Model = "XX_R345_Network_LA10_100k.mra"
XO_Model = "XO_R345_Network_LA10_100k.mra"

XX_MIDAS_file = "./INPUT/XX_Network_Completed/XX_R345.csv"
XO_MIDAS_file = "./INPUT/XO_Network_Completed/XO_R345.csv"


XX_model_completed = rebuildModel(file.path(XX_LA10_Folder,XX_Model), XX_MIDAS_file)
XO_model_completed = rebuildModel(file.path(XO_LA10_Folder,XO_Model), XO_MIDAS_file)


#mkdir _ "cache" : create a directory called cache for the profiles saved as rds files
folder = "cache/"
dir.create(folder)


#profiles_XX = profileLikelihood(mramodel_XX, nb_steps, nb_cores=min(cores, length(mramodel_XX$parameters)), const_params = c("iPi3k","iMek","iJak","iIgfr","iGsk3","iFgfr","iBmp4r")); 

# Need to specify  const_params with parameter IDs for the parameters to be kept constant in this function(checkProfileLikelihood) 
Inh_CONST_profiles_XX = checkProfileLikelihood("Inh_CONST_XX_PL",XX_model_completed, nb_steps, nb_cores=min(cores,length(XX_model_completed$parameters)), const_params = c(22,23,24,25,26,27,28));
# When using the checkProfileLikelihood function, the profiles get saved as rds files in the folder called "cache" that was created earlier ?

XX_model_completed_pl = addPLinfos(XX_model_completed, Inh_CONST_profiles_XX);
exportModel(XX_model_completed_pl,"XX_model_completed_pl.mra")
niplotPL(Inh_CONST_profiles_XX, data_name="Inh_CONST_XX_PL", folder="./OUTPUT/")

Inh_CONST_profiles_XO = checkProfileLikelihood("Inh_CONST_XO_PL",XO_model_completed, nb_steps, nb_cores=min(cores,length(XO_model_completed$parameters)), const_params = c(22,23,24,25,26,27,28));
XO_model_completed_pl = addPLinfos(XO_model_completed, Inh_CONST_profiles_XO);
exportModel(XO_model_completed_pl,"XO_model_completed_pl.mra")
niplotPL(Inh_CONST_profiles_XO, data_name="Inh_CONST_XO_PL", folder="./OUTPUT/")




