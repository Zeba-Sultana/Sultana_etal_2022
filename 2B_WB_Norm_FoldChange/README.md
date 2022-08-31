
**Aim : Normalization of Western Blot data by Mean of the Gel and then calculation of Fold change over untreated control**

1.) The input data is at : 
../1B_WB_SignalQuantification/OUTPUT/Replicate3_Alldata.xls
../1B_WB_SignalQuantification/OUTPUT/Replicate4_Alldata.xls
../1B_WB_SignalQuantification/OUTPUT/Replicate5_Alldata.xls

2.) The script for this calculation : Script_2B_WB_GelMeanNormalization_FCcalculation.R

3.) Tasks accomplished in the script :
* In Processing Step 1B the Raw signal from immunoblots have already been read in. The data had been annotated and normalization for loading differences between lanes had been done using Total Protein Stain(_TNorm). This normalisation was found to be better than that over house keeping protein GAPDH(GNorm). 

* Since the treatment samples for each replicate had to be spread over 7 gels, we needed to normalise the signals over gels so as to make them comparable.
* Two normalization strategies were compared for this : 
(i)Normalization over a common sample that was loaded on every gel. 
(ii)Normalization over mean signal of each gel.

We decided to use the second normalization method(ie normalization over mean signal per gel), since it provided more consistent normalisation between replicates.
*Normalization over Gel mean(Steps followed)*     
*step1 : For each gel of each of the 3 replicates, the mean value for each of the analytes was calculated  (Detail : However, before calculating the mean values, the samples which had NA values (1 in XX3,1 in XX5 and 10 in XO5) were also removed from the other 2 replicates before calculating the mean, so that there is no skew.)    
*step 2 : For all analytes, the T_Norm value was divided by the mean T_Norm value of that analyte for the Gel.

## Calculation of fold change over untreated control(Steps followed)
In order to use Fold change values as input for STASNet, we needed to retain variation in control values.
If we simply calculate fold change over the untreated control present on each gel, the control values on each of those gels becomes 1 and there is no variation in the control values.(Here we needed to take fold change values since analytes have been measured with different techniques - 2 different assay kits of Bioplex and WB )
Therefore FC is calculated as follows :
step1 : Calculate the mean of control values across the 3 replicates of each gel.
step2 : Use the mean calculated above to divide all values on the 3 replicates of the gel to get fold change values


********************************

4.) Output saved at : OUTPUT and OUTPUT_PAPER

