Dose response to Activin and Activin receptor inhibitor (Fig 6C) :  

Steps Followed :

**Script 1 : Reading in the data**      
Using the function ReadInWBData         

**Script 2 : Plotting the data**    
1.) Read in the files output by script 1(Smad2p values for ActDR and SBDR) into the following dataframes : Activin_DR_Smad2 and SB_DR_Smad2     
2.) Select the required columns and rename for ease of identfication        
3.) Plot the results using custom function.        
Note on the x-axis : Aim was to have the plot such that moving from left to right increases the signaling via Activin pathway.         
For increasing Activin signaling it is straight forward. Plotted the log2(Activin concentration+1).
For decreasing signaling via Activin pathway we used Activin receptor inhibitor.  In this case plotted -1*log2(SB conc+1) on the x-axis, to invert the sequence (highest concentration of the inhibitor lies at the left most side of the plot).        
        