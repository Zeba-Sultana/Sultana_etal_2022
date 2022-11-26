### Sultana_etal_2022 : 
## Modeling the cell signaling network in mouse embryonic stem cells provides insights into sex-based differences in pluripotency maintenance and differentiation
Zeba Sultana(1), Mathurin Dorel(2), Bertram Klinger(2) , Anja Sieber(2) , Nils Blüthgen(2) , Edda G. Schulz(1)    

1 Regulatory Networks in Stem Cells, Max Planck Institute for Molecular Genetics, Berlin    
2 Computational Modelling in Medicine, Charite - Universitätsmedizin, Berlin 




The complete analysis and figures from the paper can be generated using the data and scripts in this repository.
The folder has been organised into subfolders, each of which accomplishes one step of the data preprocessing and analysis.
The following steps can be executed sequentially for the complete analysis.  

**Step 1 :** Extracting the data from perturbation experiments and its quantification
1A.) Bioplex Assay results(lxb) converted to MIDAS format files.
1B.) Western Blot Assay analyte quantification read in and annotated

**Step 2:** Normalisation of the 2 datasets and calculation of fold change over untreated control:
2A.) Normalisation of Bioplex Assay analytes and fold change over untreated control.
2B.) Normalisation of Western Blot analytes and fold change over untreated control.

**Step 3 :** Merging the fold change data from the two assays and preparation of the combined MIDAS files that are used as input to STASNet for modeling.

**Step 4 :** Building the model and network extension that resulted in addition of 10 links.

**Step 5 :** Profile Likelihood analysis to compute confidence interval for the model parameters.

**Step 6 :** Comparing the experimental data with simulation results from the initial network model and the completed network model.

**Step 7:** Analysis of results from validation experiments.

Each of the subfolders has its readme file that provides further details on that step.

The entire analysis can also be run using the script provided in this folder "Script_AllSteps", which moves into these folders and executes the required scripts sequentially. 

Note :   
1.) In the perturbation experiments, the replicates have been referred to as R3,R4 and R5 for all the preprocessing steps.(The first two replicates of these experiments were not used). These correspond to replicates Rep1 , Rep2 and Rep3 in the paper.

2.) The output from scripts is saved into folders called OUTPUT or OUTPUT_PAPER.   
OUTPUT : Has outputs that are used by later scripts.   
OUTPUT_PAPER : Has output(Figures or Supplementary Tables) that are part of the paper. These also have the correct replicate number annotation that is used in the paper.

