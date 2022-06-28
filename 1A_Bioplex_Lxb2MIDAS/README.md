### Aim : Read in the data from Bioplex Assay(as lxb files) and convert them into MIDAS file

Common Note : I had done 5 replicates of the experiments for testing various factors. The first 2 replicates have not been used for any of the final data and model building. The replicates used for the paper were internalising called R3,R4 and R5.(Hence scripts and file names have these reference. Later they were renamed as follows :   
R3 --> Rep1  
R4 --> Rep2  
R5 --> Rep3   

**SUMMARY**

Files in this folder and their functions :

* Custom script for converting lxb files to MIDAS files : AllReplicates_perform_lxb_extraction.R
It can be executed from the terminal using "Script AllReplicates_perform_lxb_extraction.R"
The R package "lxb" is required for the execution.

* This script calls the input files and other scripts that are in 2 folders in this directory :
(i) INPUTfiles_BioplexR3_R4 : For replicates 3 and 4
(ii) INPUTfiles_BioplexR5 : For replicate 5

The input files being used for the Replicates3/4 are different from that of replicate5 because the experimental plate design was different in R5 as compared to R3 and R4. Secondly, lesser number of analytes were assayed in R5( We retained only those analyses in the Bioplex Assay that had shown a reliable signal in the 11-plex assay done on R3 and R4)

* The MIDAS FILES to be used in the next steps are in OUTPUT_MIDASfiles


**DETAILED NOTES**

1.) AllReplicates_perform_lxb_extraction.R : Script to extract data from lxb folders in the form of corresponding MIDAS files. This needs to be done for each of the 6 lxb folders(3 replicates of XX and XO) - being done in for loops in this script.

2.) extract_lxb_data.R : A script which contains the definitions of functions required for this process.

3.) The output from the Bioplex machine(Magpix) are in the form of lxb files in “lxb” folder @ ../RAW_DATA/BioplexLysates_LUMINEX

4.) Files required to pull information on the Bioplex assay plate layout:
*********************************************************
File 1.) experimental_layout.xls : This is an excel file which summarizes the plate design of the Bioplex assay as a table, the location of different samples, control and blank in the 96 well plate(ie the plate design)

File 2.) extract_experiment_annotations.R : In this file the following information needs to be updated based on current inputs :
(i) annotations_file = The excel file described in 1 (For eg : “experimental_layout.xls”) and input the name of the sheet as the second argument.
(ii) Define “Inhibitors” and “Stimulations”. The inhibitors should be name of the species(proteins) followed by “i”. (Note :Take care to have the names of the proteins spelled exactly the same here and in the Network.tab and Basal_activity.dat files will be used along with the MIDAS files for model creation.)

Define the bead-region and corresponding analyte names from the bioplex assay used by you.
Executing the above results in the following :
-the table of experimental layout is read into a R-matrix which captures the plate design.
-the stimulators and inhibitors are defined.
-the bead-regions are mapped to the correct analyte names.
*******************************************************
EXTRACTING INFORMATION FROM LXB to get MIDAS : Excution using all the files listed above from (1 to 4) : The input files being used for the Replicates3/4 are different from that of replicate5 because the experimental plate design was different in R5 as compared to R3 and R4. Secondly, lesser number of analytes were assayed in R5.

INPUTfiles_BioplexR3_R4 : Folder with files specific for replicate 3 and 4 lxb files
(Plate_Layout_BioRad_R3_R4.xls and extract_experiment_annotations_R3_R4.R)

INPUTfiles_BioplexR5 : Folder with files specific for replicate 5 lxb files
(Plate_Layout_BioRad_R5.xls and extract_experiment_annotations_R5.R)


5.) OUTPUT : All the MIDAS files (which will be needed for model input) are copied into one folder called OUTPUT_MIDASfiles and are renamed for consistent nomenclature. 






