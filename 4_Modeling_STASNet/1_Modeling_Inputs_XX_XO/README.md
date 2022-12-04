
## Important functions from the STASNet package used in the fitscript 

https://github.com/molsysbio/STASNet/blob/master/vignettes/STASNet_vignette.html    
*Note :* the above link provdies documentation of basic STASNet functions and the workflow and for modeling biological networks using perturbation data. It describes the input files and steps for modeling a dataset using a toy example of 5 nodes.

The important functions used in the "fitscript" for this modeling project are also briefly listed below :
  
* createModel(): Creates a maximum likelihood MRAmodel using the using the perturbation data(.csv file), literature derived starting network topology(.tab file) and information about basal activity of nodes(.dat file).  
* plotModelGraph() : Creates a graphical representation of model(the network structure along with the associated parameter values.)
* plotModelAccuracy() : Produces heatmaps to help evaluate how simulation using the model compares to experimental data.
* plotModelScores() : Plots the individual goodness of fit for each analyte as coefficient of Determination (R²) 

* selectMinimalModel() : Function to identify and remove dispensible links in the model which do not contribute significantly to overall fit. In other words, those links which when removed do not result in a significant increase in model residuals.  
* suggestExtension() : Function to find links that when added to model will result in improvement of fit. The function outputs a list of suggested links along with difference in model residual and the associated p-values.





  
