**Aim : Normalization of Bioplex readings by Mean of the replicate and then calculation of Fold change over untreated control**


1.) The input data is the MIDAS files generated from lxb files

2.) The script for normalisation and fold change over DMSO : PS2A_Script_MeanNormalization_FCcalculation.R

3.) Tasks accomplished :

*S2A(i) : Bioplex : Normalization over mean signal of replicate*

The bioplex data is normalized over the the mean signal of the replicate. The calculations are done as follows :
1.) Read in the raw MFI values.
2.) Subset the data for replicate and x-status. Eg Replicate3-XX
3.) There are different number of control samples across replicates. R3 and R4 have 9 controls while R5 has 2. Therefore the controls are treated as one kind of treatment. Calculate the mean of the controls. Take these mean of control values to be 1 treatment.
4.) Take the values for the remaining 53 treatments along with mean of controls calculated above. Calculate the overall mean. This is mean_replicate_x_status.
5.) Divide each data point of the subdata(Eg Replicate3-XX) by the mean values calculated in step 4. 
This calculation also ensures that the control values are not exactly 1 but vary around 1.


*S2A(ii) : Calculation of Fold change over untreated control*

The normalized data(over mean signal in this case), is used to calculate fold change over untreated control - By diving each data point by the mean of the control values. The Fold change values are then taken as input to STASNet.


4.) Fold Change over control to be used as input to STASNet is saved at : ./OUTPUT

5.) Supplmentary Table 2 for anaytes measured by Bioplex, comprising Raw MFI values, normalisation over mean signal and fold change values are saved at : ./OUTPUT_PAPER. They have the corrected replicate numbers.
