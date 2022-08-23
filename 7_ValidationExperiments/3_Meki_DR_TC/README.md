In this folder, I make the validation figure for Mek-Inhibitor DoseResponse. 


Summary of the steps followed in the R script "1_Meki_DR_ValidationFigures.R" are as follows.

######### STEP 1 : Reading in the data and normalization over the mean signal per gel/replicate ######### 
1.) read in the quantification file and assign the sample names.
2.) filter channel 800 values as phorphorylated, and 700 values as unphorphorylated signal.
Samples were distributed across 3 gels for each analyte.
3.) Use the mean signal of each gel to normalize the signal values across the gels.
4.) Finally, for each replicate, I calculate the FC over XX control and FC over respective cell line control.


######### STEP 2 : Plotting ######### 

Plotted using the function Plot_TwoPanel_ValidationPlot_updated

######### STEP 3 : Doing the t-test #########

Compare the mean of XX and XO values and save the results as a df

In this folder, I make the validation figure for Mek-Inhibitor Time Course. 

Summary of the steps followed in the R script "2_Meki_TC_ValidationFig.R" are as follows

######### STEP 1 : Reading in the data and normalization over the mean signal per gel/replicate ######### 

For the MekiTC expeiments this was done using the function "prep_analyte_G12". 
Tasks done in the function :
1.) Read in the TPS quantification file and add well numbers.
2.) Read in the analyte quantification file and add well numbers. 
3.) Merge based on common well numbers.
4.) Add annotation such as analyte name, replicate number etc and calculate Signal_by_TPS = Signal/TPS
5.) Last step : Normalization over the mean signal of the replicate.

Next, using the function "FC_Calculation", fold change over XX control and that over untreated control of the respective cell line was calculated.

######### STEP 2 : Plotting ######### 

Plotted using the function Plot_TwoPanel_ValidationPlot

######### STEP 3 : Doing the t-test #########

Compare the mean of XX and XO values and save the results as a df









