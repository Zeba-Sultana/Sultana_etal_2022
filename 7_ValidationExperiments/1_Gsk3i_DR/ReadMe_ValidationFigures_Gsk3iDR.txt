In this folder, I make the validation figure for GSK3i DR treatment time-course :

The functions needed for data preprocessing are defined in : "../ValidationPlot_Functions/U_Functions_ValidationExp_Analysis.R"
The functions needed for data plotting are defined in : "../ValidationPlot_Functions/U_Functions_ValidationExp_Plotting.R"

Summary of the steps followed in the R script are as follows :

1.) Read in the WB data quantification using the function "ReadInWBData". Samples were distributed across 3 gels for each analyte.
2.) For comparison across gels, for each analyte, I normalized the signal of each lane by the mean signal on that gel using the function "Norm_o_MeanRep"
3.) Then I calculate the fold change over XX control(FCoXX) or fold change over respective cell line control(FCoCntrl). In both cases, first the mean of the control sample across the 3 replicates is calculated, then every signal is divided by this mean to obtain the FC.
4.) Finally  plotting the figures for the paper is done by using the function "Plot_TwoPanel_ValidationPlot_updated"

