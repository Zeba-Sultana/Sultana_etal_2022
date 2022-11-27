
**Aim : To model the signalling networks in the two cell types, using experimentally determined perturbation data set and literature derived network. Secondly we extended the literature network to find links that can improve the fit to experimental data.**

Documented in various subfolders here are the input and major steps followed to achieve the above aims.

* Subfolder : 1_Modeling_Inputs_XX_XO    
The Script "fit_script_withExtRed.R" was used for building the model. This script was used separately with perturbation input file for XX (Suppl_Table_S4_XX_MIDAS.csv) and that for XO(Suppl_Table_S4_XO_MIDAS.csv. These are the input files carrying results from the perturbation experiments in the MIDAS format.   
The other input files that were common for both cell lines inlcude :   
(i) Network.tab : This file describes the network. It is tab separted file where each link present in the network is described on one row - first column has the name of the outgoing node, followed by a tab and name of the incoming node in the second column. Initailly, Network.tab will enlist all the links present in the literature-derived netwrok. During the network extension preocedure the selected link is added to the Network.tab.   
(ii) BasalActivity.dat : File that enlists all nodes that are assumed to have basal signaling activity, which are all nodes except the ligands(Lif, Igf1,Fgf4,ActA and Bmp4)   


* Subfolder : 2_Model_Extension   
Buiding the initial models for XX and XO using the "fit_script_withExtRed.R", results in the model(.mra) file along other files that decribe the overall fit to data achieved.  Additionally, since the fit script has "extension = TRUE", additional links that can improve the model fit are tested and a file (Extension.csv) enlisting links that can improve the fit to the model in generated.   
The script "FindingCommonLink_Rscript_MidpValue.R" was used for finding the best common link from the links suggested separately for the XX and XO models.   
This link was added to the Network.tab file and saved as Network_LA1.tab and a new model using this network was built. This procedure was repeated until no common links that cross the selection threshold were found.    
The R markdown file "Model_Extension_Threshold_p_005.Rmd" summarises the links added to the network, when filtering the links based on a significance threshold of p<0.005 in a likelihood ratio test for the improvement in fit to data and finding the best common link based on highest sum of percentage decrease in model residuals.    

* Subfolder : 3_Model_Extension_Different_Thresholds   
The results when using different p-value thresholds for filtering the links. Each of these also write out the residual values associated with the models at each link addition step with the different thrshold values :
Threshold_p_05.Rmd (writes out the residuals in : Links_Residuals_p_05.csv)
Threshold_p_01.Rmd (-> Links_Residuals_p_01.csv)
Threshold_p_005.Rmd (-> Links_Residuals_p_005.csv)
Threshold_p_001.Rmd (-> Links_Residuals_p_001.csv)

* Subfolder : 4_Plotting_ModelResiduals_Fig_2b_S2c   
Script for plotting the model residuals(Figures 2B and Suppl2C in the paper)



