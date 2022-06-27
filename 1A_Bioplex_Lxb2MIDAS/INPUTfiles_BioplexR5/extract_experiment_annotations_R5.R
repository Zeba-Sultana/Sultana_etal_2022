# Script specific to the dataset at hand,
#process the annotations to get an format usable by the "extract_lxb_data.R" script and 
#provide the correspondence between the bead IDs ('region') and the proteins they identify ('analyte')
library("gdata") # to read excel files

# Load the annotations of the wells
#annotations_file = "/Users/sultana/PHD_PROJECT_Zeba/Zeba_PhD_PAPER/SCRIPTS_Sultana_etal_2021/PS1A_BP_Lxb2MIDAS/INPUTfiles_BioplexR5/Plate_Layout_BioRad_R5.xls" # Name of the xls file that captures the experimental design
#annotations_file = "/project/ag_schulz/Zeba/SCRIPTS_Sultana_etal_2021/PS1A_BP_Lxb2MIDAS_updated/INPUTfiles_BioplexR5/Plate_Layout_BioRad_R5.xls" # Name of the xls file that captures the experimental design
annotations_file = "./INPUTfiles_BioplexR5/Plate_Layout_BioRad_R5.xls" # Name of the xls file that captures the experimental design


annot = read.xls(annotations_file, "Layout")    # Name of the sheet in the xls file above which has the data : "Layout" in this case
cell_lines = 
plate = matrix( as.matrix(annot[1:8,2:13]), nrow=8, ncol=12, dimnames=list(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12) )
plate = gsub("\\(?No Lif\\)?", "+NoLif", plate) # NoLif pseudostimulation : Removing brackets and spaces in the name
plate[plate=="Blank"] = "blank"                 # Blank
plate[plate=="DMSO"] = "control"                # Control
# plate = gsub("blank", "", plate)
# plate = gsub("Blank", "blank", plate)


plate = gsub(" \\+ ", "+", plate)              # No space
plate = gsub(" ", "", plate)                   # No space
plate = gsub("^\\+", "", plate)               # Removing + from beginning
plate = gsub("Fgf4i", "Fgfri", plate)         # To correct those places where Fgfri had been spelled incorrectly
plate = gsub("Bmp4i", "Bmp4ri", plate)       # Bmp4-inh being corrected to Bmp4 receptor-inh
plate = gsub("Gsk3bi", "Gsk3i", plate)               # Gsk3b-inh being corrected to Gsk3-inh(inhibits both alpha and beta isoforms)



inhibitions = c("Jaki", "Bmp4ri", "Gsk3i", "Fgfri", "Igfri", "Meki", "Pi3ki") # Define here the name of the inhibitors #Changed name to just Pi3ki and Gsk3i
receptors = c("Activin", "Fgf4", "NoLif")

region <- c( 18, 27, 46,75)
analyte <- c('Gsk3', 'Mek','mTor','Akt') 

#Check_mydf <- data.frame(analyte,region) # to check if the bead regions are being associated with the correct analyte.

#plate_stim = matrix( "blank", nrow=8, ncol=12, dimnames=list(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12) )
#plate_stim[2:4,] = matrix(rep(receptors[1:ncol(plate_stim)], 3), nrow=3, byrow=T)
#plate_stim[5:7,] = matrix(rep(receptors[(ncol(plate_stim)+1):(2*ncol(plate_stim))], 3), nrow=3, byrow=T)

#plate_inhib = matrix( "blank", nrow=8, ncol=12, dimnames=list(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12) )
#plate_inhib[2:4,] = matrix(rep(inhibitions[1:ncol(plate_stim)], 3), nrow=3, byrow=T)
#plate_inhib[5:7,] = matrix(rep(inhibitions[(ncol(plate_stim)+1):(2*ncol(plate_stim))], 3), nrow=3, byrow=T)



				



































