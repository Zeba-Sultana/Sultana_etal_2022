
Aim : To create confidence intervals for the parameters of the model using profile likelihood.

The profile likelihood analysis is done separately for the two models and later the CI for the parameters can be plotted together for comparison.

Inputs needed are in folder "INPUT" 
1.) The model(mra file) - We use the completed network model after addition of the 10 links.
2.) The perturbation data(MIDAS file)
We rebuild the model using these inputs in the STASnet function rebuildModel()

Scripts needed : 
1.) 1_ComputeProfileLikelihood_CONST_inhParameters.R : First Profile likelihood computation is done on the models created using rebuild model. This can be parallelised and therefore running on the server is better. 
2.) 2_Analyse_ProfileLieklihood_PlotResults.R : Next, the saved profile information in the cache folder is accessed, the 95% thresholds at 8 degrees of freedom is updated and then results shown in the paper are plotted.
2a_Functions_ProfileLikelihood.R is a helper script that defines some functions used for the analysis. Needs to be in the same folder.






