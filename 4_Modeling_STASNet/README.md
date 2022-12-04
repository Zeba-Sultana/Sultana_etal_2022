
**Aim : To model the signalling networks in the two cell types, using experimentally determined perturbation data set and literature derived network. Secondly we extended the literature network to find links that can improve the fit to experimental data.**

Documented in various subfolders here are the input and major steps followed to achieve the above aims.

* **Subfolder : 1_Modeling_Inputs_XX_XO**    
 Input files used for modeling :      
(i) Network.tab : This file describes the network. It is a tab separated file where each link present in the network is described on one row - first column has the name of the outgoing node, followed by a tab and name of the incoming node in the second column. Initailly, Network.tab will enlist all the links present in the literature-derived netwrok. During the network extension procedure, the selected link is added to the Network.tab.   
(ii) BasalActivity.dat : File that enlists all nodes that are assumed to have basal signaling activity, which are all nodes except the ligands(Lif, Igf1,Fgf4,ActA and Bmp4)      
(iii) Suppl_Table_S4_XX_MIDAS.csv(for XX) and Suppl_Table_S4_XO_MIDAS.csv(for XO) : Results from systematic perturbation experimentas in the two cell lines in MIDAS format.    
(iv) fit_script_withExtRed.R : Script for building the model, plotting model accuracy, finding non-essential links and finding links that improve model fit 

The important functions used in the "fitscript" for this modeling project are  briefly listed below :
* createModel(): Creates a maximum likelihood MRAmodel using the using the perturbation data(.csv file), literature derived starting network topology(.tab file) and information about basal activity of nodes(.dat file).  
* plotModelGraph() : Creates a graphical representation of model(the network structure along with the associated parameter values.)
* plotModelAccuracy() : Produces heatmaps to help evaluate how simulation using the model compares to experimental data.
* plotModelScores() : Plots the individual goodness of fit for each analyte as coefficient of Determination (R²) 
* selectMinimalModel() : Function to identify and remove dispensible links in the model which do not contribute significantly to overall fit. In other words, those links which when removed do not result in a significant increase in model residuals.  
* suggestExtension() : Function to find links that when added to model will result in improvement of fit. The function outputs a list of suggested links along with difference in model residual and the associated p-values.

*Note :* For complete documentation of basic STASNet functions and the workflow and for modeling biological networks using perturbation data please refer the following url
https://github.com/molsysbio/STASNet/blob/master/vignettes/STASNet_vignette.html    
It describes the input files and steps for modeling a dataset using a toy example of 5 nodes.


* **Subfolder : 2_Model_Extension**   
Buiding the initial models for XX and XO using the "fit_script_withExtRed.R", results in the model(.mra) file along other files that decribe the overall fit to data achieved.  Additionally, since the fit script has "extension = TRUE", additional links that can improve the model fit are tested and a file (Extension.csv) enlisting links that can improve the fit to the model in generated.   
The script "FindingCommonLink.R" (which calls "FindingCommonLink_FuncDefinitions.R") was used to find the best common link from the links suggested separately for the XX and XO models. The extension list for the two models is provided as input to "FindingCommonLink.R", which does the follwoing steps :
1.) Read in the extension.csv for both cell lines   
2.) Remove links from TFs and Ligands   
3.) Retain only links that have adjusted p-value < selected threshold and value !=0   
4.) Calculate %improvement in residual for each cell line   
5.) Find the intersection list,(Links that appear in both lists)   
6.) Calculate the sum of %improvement in residual   
7.) Order the intersection list by descending order of sum of %improvement in residual.   
The topmost link of this oredered list was added to the Network.tab file and saved as Network_LA1.tab.  A new model using this updated network was built. This procedure was repeated until no common links that cross the selection threshold were found. 

FindingCommonLink_Rscript_MidpValue.R : Summarizes the link selection at all link addition steps and saves the orderred list of link suggestions at each step.
The R markdown file "Model_Extension_Threshold_p_005.Rmd" summarises the links added to the network, when filtering the links based on a significance threshold of p<0.005 in a likelihood ratio test for the improvement in fit to data and finding the best common link based on highest sum of percentage decrease in model residuals.    

* **Subfolder : 3_Model_Extension_Different_Thresholds**  
The results when using different p-value thresholds for filtering the links. Each of these also write out the residual values associated with the models at each link addition step with the different thrshold values :
Threshold_p_05.Rmd (writes out the residuals in : Links_Residuals_p_05.csv)
Threshold_p_01.Rmd (writes out the residuals in : Links_Residuals_p_01.csv)
Threshold_p_005.Rmd (writes out the residuals in : Links_Residuals_p_005.csv)
Threshold_p_001.Rmd (writes out the residuals in : Links_Residuals_p_001.csv)

* **Subfolder : 4_Plotting_ModelResiduals_Fig_2b_S2c**   
Script for plotting the model residuals(Figures 2B and Suppl2C in the paper)



