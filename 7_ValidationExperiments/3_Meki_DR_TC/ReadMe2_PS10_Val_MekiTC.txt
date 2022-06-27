In this folder, I make the validation figure for Mek-Inhibitor Time Course. 

Summary of the steps followed in the R script "2_Meki_TC_ValidationFig.R" are as follows

######### STEP 1 : Reading in the data and normalization over the mean signal per gel/replicate ######### 

For the MekiTC expeiments this was done using the function "prep_analyte_G12". 
Tasks done in the function :
1.) Read in the TPS quantification file --> remove extra rows and columns --> Add well numbers
2.) Read in the analyte quantification file --> select relevant columns --> Add well numbers 
3.) Merge based on common well numbers.
4.) Add annotation such as analyte name, replicate number etc and calculate Signal_by_TPS = Signal/TPS
5.) Last step : Normalization over the mean signal of the replicate.

Next, using the function "FC_Calculation", fold change over XX control and that over untreated control of the respective cell line was calculated.

######### STEP 2 : Plotting ######### 

To make the figure in paper(pMek and pcRaf together as facetwrap), I did the following :
1.) Selected relevant columns and renamed them to make then compatible with the plot function written for the DR experimental data analysis. Joined the data of the two analytes, using the annotation columns.
2.) Plotted using the function Plot_TwoPanel_ValidationPlot

######### STEP 3 : Doing the t-test #########

Compare the mean of XX and XO values and save the results as a df






