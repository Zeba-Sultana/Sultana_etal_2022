## Aim : Read in the data from WB 

Common Note : I had done 5 replicates of the perturbation experiments for testing various factors. The first 2 replicates have not been used for any of the final data and model building. 
The replicates used for the paper were internally called R3,R4 and R5.(Hence scripts and file names have these reference. Later they were renamed as follows :    
R3 --> Rep1  
R4 --> Rep2   
R5 --> Rep3   

- Script_1B_WB_DATA_ReadIn.R : The script used to read in the quantification values from xls files that were exported form Odyssey image quantification and TPS quantification.    

- In replicate 5, some lysates had been lost during filtration. These were less in amount and not enough for Western Blotting. The data rows corresponding to these treatments have been filled with NA. Similarly if there was a problem with accurate quantification of a signal in Western blot image due to problem in the background, they were replaced with NA.

- OUTPUT : The WB signals, along with their normalisation over house keeping protein GAPDH(GNorm) and over Total protein stain(TNorm) and associated annotations are written out as xls files in the folder "OUTPUT". TNorm was found to be better and hence that was used going forward


