In this folder, I make the validation figure for Mek-Inhibitor DoseResponse. 


Summary of the steps followed in the R script "1_Meki_DR_ValidationFigures.R" are as follows.

######### STEP 1 : Reading in the data and normalization over the mean signal per gel/replicate ######### 
1.) read in the quantification file and assign the sample names.
2.) filter channel 800 values as phorphorylated, and 700 values as unphorphorylated signal.
Samples were distributed across 3 gels for each analyte :
Mek+pMek and TPS : Gel24(R1), Gel25(R2), Gel16(R3)
pcRaf+TPS : Gel24(R1), Gel25(R2), Gel16(R3)

3.) use the mean signal of each gel to normalize the signal values across the gels : I had done it by two methods to cross-check results.Both give almost the same values. I have selected M2 for plotting :
M1 : I divide Signal by total(protein/stain) for each well. Then I normalize each of these Signal_by_Total over the mean Signal_by_Total.
M2 : I divide Signal of each well by mean signal for the gel. Similarly,I divide total(protein/stain) of each well by mean total(protein/stain) for the gel. Then I divide normalized signal for each well over the normalized total(protein/stain) for that well.

4.) Finally, for each replicate, I calculate the FC over XX control and FC over respective cell line control.


######### STEP 2 : Plotting ######### 

Plotted using the function Plot_TwoPanel_ValidationPlot_updated

######### STEP 3 : Doing the t-test #########

Compare the mean of XX and XO values and save the results as a df




