In this folder PS3 the following steps are accomplished :  

Aim : Convert the dataframes containing the normalized Foldchange over control values of BP and WB into MIDAS files and then merge them.

1.) The input data is at : 
BP_Normalized <- "../2A_Bioplex_Norm_FoldChange/OUTPUT/Bioplex_FC_STASNet.xls"
WB_Normalized <- "../2B_WB_Norm_FoldChange/OUTPUT/WBN2_GelMean_Norm_FC_Statsnet.xls"

2.) The script for this : Script_3_Merge_BP_WB.R

3.) Tasks accomplished in the script :
************************************************
##Bioplex Data
(i) Read in Bioplex data(fold change values) 
(ii) Replace the Akt values for XXR5 and XOR4 with NA because these replicates were not consistent with the other two replicates due to very high variability.
(iii) Covert there dataframes to matrices and drop the annotation columns.
(iv) This was needed so that the matrices can be used to generate the MIDAS file using the STASNet function midasFromData.


##Western Blot Data
(i) Read in WB data(fold change values) 
(ii) Correct the nomenclature where needed to make consistent.
(iii) Covert there dataframes to matrices and drop the annotation columns.
(iv) This was needed so that the matrices can be used to generate the MIDAS file using the STASNet function midasFromData.

## Merge the MIDAS files generated from BP and WB data using the function Merge_Bioplex_WB_MIDAS
Since there are different number of controls in the Bioplex and WB data, first the treatment and control rows are separated. Treatment rows are easily merged using left join. But for control rows, I check which one has more rows, make the other one of equal size by filling in with NA rows and then merge them.
Lastly I make sure to drop the bioplex analytes that did not give correct signal for the remaining analytes (such as p38 etc). bCatenin as a Western Blot read out is also dropped because of very large variation in data.

## finally,
Merged MIDAS files for XX Replicates(XX_R345.csv)
Merged MIDAS files for XO Replicates(XO_R345.csv)


********************************
4.) Output saved at : ./OUTPUT
XX_R345.csv and XO_R345.csv will be used as input to STASNet

4.) Output for supplementary tables S4 and S5 used in paper saved at : ./OUTPUT_PAPER
Suppl_Table_S4_XX_MIDAS.csv and Suppl_Table_S5_XO_MIDAS.csv





