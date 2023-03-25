In this folder, I make the validation figure for dose response to Activin and Activin receptor inhibitor  :  

Aims :    
(I) Read in the data from experiments measuring dose response to Activin and Activin receptor inhibitor       
(II) Plot the two toegther in such a way that it shows effect of increased signaling via the activin pathway going from left to right.

Accordingly the steps needed to accomplish these two aims are in 2 files:

**Script 1 : Reading in the data**
Using the function ReadInWBData

**Script 2 : Plotting the data**
1.) Read in the files output by script 1(Smad2p values for ActDR and SBDR) into the following dataframes : Activin_DR_Smad2 and SB_DR_Smad2
2.) Select the required columns and rename for ease of identfication
3.) Plot the results using cutom function.
Note on the x-axis :
For increasing Activin signaling, Activin dose response, it is straight forward. Plotted the log2(Activin concentration+1) on x-axis.
For decreasing signaling via Activin pathway we used Activin receptor inhibitor. So to have the plot such that moving from left to right increases the signaling via Activin pathway, multiplied log2(SB conc+1) with -1 to invert the sequence. So the values corresponding to the highest concentration of the inhibitor lie at the left most side of the plot.



