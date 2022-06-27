# Script specific to the dataset at hand, process the annotations to get an format usable by the "extract_lxb_data.R" script and provide the correspondence between the bead IDs ('region') and the proteins they identify ('analyte')
library("gdata")

# Load the annotations of the wells
#annotations_file = "/Users/sultana/PHD_PROJECT_Zeba/Zeba_PhD_PAPER/SCRIPTS_Sultana_etal_2021/PS1A_BP_Lxb2MIDAS/INPUTfiles_BioplexR3_R4/Plate_Layout_BioRad_R3_R4.xls"
#annotations_file = "/project/ag_schulz/Zeba/SCRIPTS_Sultana_etal_2021/PS1A_BP_Lxb2MIDAS_updated/INPUTfiles_BioplexR3_R4/Plate_Layout_BioRad_R3_R4.xls"
annotations_file = "./INPUTfiles_BioplexR3_R4/Plate_Layout_BioRad_R3_R4.xls"

annot = read.xls(annotations_file, "Layout")
#annot = read.xls(annotations_file)

cell_lines = 
plate = matrix( as.matrix(annot[1:8,2:13]), nrow=8, ncol=12, dimnames=list(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12) )
plate = gsub("\\(?No Lif\\)?", "+NoLif", plate) # NoLif pseudostimulation
#plate[plate==""] = "blank" # Blank
plate[plate=="DMSO"] = "control" # Control

plate = gsub(" \\+ ", "+", plate) # No space
plate = gsub(" ", "", plate) # No space
plate = gsub("^\\+", "", plate)
plate = gsub("Fgf4i", "Fgfri", plate)         # To correct those places where Fgfri had been spelled incorrectly
plate = gsub("Bmp4i", "Bmp4ri", plate)       # Bmp4-inh being corrected to Bmp4 receptor-inh
plate = gsub("Gsk3bi", "Gsk3i", plate)               # Gsk3b-inh being corrected to Gsk3-inh(inhibits both alpha and beta isoforms)


inhibitions = c("Jaki", "Bmp4ri", "Gsk3i", "Fgfri", "Igfri", "Meki", "Pi3ki") # Define here the name of the inhibitors #Changed name to just Pi3ki and Gsk3i
receptors = c("Activin", "Fgf4", "NoLif")

#plate_stim = matrix( "blank", nrow=8, ncol=12, dimnames=list(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12) )
#plate_stim[2:4,] = matrix(rep(receptors[1:ncol(plate_stim)], 3), nrow=3, byrow=T)
#plate_stim[5:7,] = matrix(rep(receptors[(ncol(plate_stim)+1):(2*ncol(plate_stim))], 3), nrow=3, byrow=T)

#plate_inhib = matrix( "blank", nrow=8, ncol=12, dimnames=list(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12) )
#plate_inhib[2:4,] = matrix(rep(inhibitions[1:ncol(plate_stim)], 3), nrow=3, byrow=T)
#plate_inhib[5:7,] = matrix(rep(inhibitions[(ncol(plate_stim)+1):(2*ncol(plate_stim))], 3), nrow=3, byrow=T)


region <- c(14, 18, 27, 36, 38, 42, 46, 52, 54, 56, 75)
analyte <- c('Smad2', 'Gsk3', 'Mek', 'p38', 'Erk', 'Src', 'mTor', 'Stat3', 'Pi3kp85', 'cJun', 'Akt')
				



































