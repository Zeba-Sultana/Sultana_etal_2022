
Aim : Read in the data from WB 

Common Note : I had done 5 replicates of the perturbation experiments for testing various factors. The first 2 replicates have not been used for any of the final data and model building. The replicates used for the paper were internally called R3,R4 and R5.(Hence scripts and file names have these reference. Later they were renamed as follows :
R3 --> Rep1
R4 --> Rep2
R5 --> Rep3

1.) WB_DATA_ReadIn.R : The script used to read in the quantification values from xls files that I manually exported form Odyssey image quantification and TPS quantification. 
2.) I use a common function to do this for replicates 3 and 4.

3.) The function used for reading in data from Replicate5 gels is slightly different. (In this case it is needed to define Wells outside the generic function because at that stage during quantification empty lanes and ladder lanes were not quantified resulting in different number of lanes per gel. Later(for R3 and R4) I quantified 20 lanes per gel and so included the empty/ladder lanes as well for consistently having 20 lanes per gel.)   

4.) In replicate 5, some lysates had been lost during filtration. These were less in amount and not enough for Western Blotting. The data rows corresponding to these treatments have been filled with NA.

5.) Finally for Gel3 of replicate 5, the Gel3 had a problem with the XO control lane. So the normalised signal for control was estimated using the ratio of XO_cntrl/Common seen in the remaining 6 gels and using the mean of this ratio to estimate XO_cntrl on Gel3 using common lysate signals.

6.) OUTPUT : Finally The WB signals, along with their normalisation over GAPDH(GNorm) and over TPS(TNorm) and associated annotations are written out as xls files in the folder "OUTPUT". TNorm was found to be better and hence that was used going forward

These have been used further along with the Luminex Assay read outs

