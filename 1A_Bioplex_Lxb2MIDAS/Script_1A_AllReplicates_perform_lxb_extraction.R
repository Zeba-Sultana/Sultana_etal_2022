#!/usr/bin/env Rscript

# Wrapper script to perform the extraction of lxb data, simply calls functions defined in the 'extract_lxb_data.R' script in the correct order

############## For Replicates 3 and 4 ############
# Since the plate setup for Replicates and 4 was the same they are processed together. Replicate5 is done later separately

Folder_R34 <- "./INPUTfiles_BioplexR3_R4"

source(file.path(Folder_R34,"extract_experiment_annotations_R3_R4.R"))

# Functions to process the raw lxb data
source(file.path(Folder_R34,"extract_lxb_data.R"))

print("Extracting treatments")
esetup = extractExperimentalSetup(plate)

# Load the data
print("Loading the lxb files for R3 and R4")


R34_lxb_folders = c("../RAW_DATA/BioplexLysates_LUMINEX/BIOPLEX_R3_R4/170126_Zeba_XX3_lxb",
                    "../RAW_DATA/BioplexLysates_LUMINEX/BIOPLEX_R3_R4/170126_Zeba_XO3_lxb",
                    "../RAW_DATA/BioplexLysates_LUMINEX/BIOPLEX_R3_R4/170126_Zeba_XX4_lxb",
                    "../RAW_DATA/BioplexLysates_LUMINEX/BIOPLEX_R3_R4/170126_Zeba_XO4_lxb" )


for(lxb_folder in R34_lxb_folders){
    if (!grepl("/$", lxb_folder)) {
        lxb_folder = paste0(lxb_folder, "/")
    }
    
do_plots = TRUE
    
full_data = readLxb(paste0(lxb_folder, "*.lxb"), text=FALSE)
lxb_name = gsub("[0-9]{6}_", "", lxb_folder)
lxb_name = gsub("_lxb", "", gsub("/$", "", lxb_name))

print("Processing bead data")
beads = readBeads(full_data, region, analyte, esetup$controls, esetup$blanks, esetup$wptb, lxb_folder, esetup$externals, do_plots)
print("Bead data processed")
if (file.exists("post_processing.R")) {
    print("Post processing the data")
    source("post_processing.R")
}
lxb_name = basename(lxb_name)
writeMIDAS(beads, lxb_name, esetup$inhibitors, esetup$stimulators, lxb_name)
print(paste0("Data writen in MIDAS format in folder ", lxb_name))

pdf(paste0(lxb_name, "/", "beads_distribution_", lxb_name, ".pdf"), width=10)
analyseBeadDistributions(full_data, region, analyte, esetup$treatments)
dev.off()

}

print("R3 and 4 : Done")


############## For Replicate 5 ############

Folder_R5 <- "./INPUTfiles_BioplexR5"


source(file.path(Folder_R5,"extract_experiment_annotations_R5.R"))
source(file.path(Folder_R5,"extract_lxb_data.R"))

print("R5 : Extracting treatments")
esetup = extractExperimentalSetup(plate)

# Load the data
print("R5 : Loading the lxb files")

R5_lxb_folders = c("../RAW_DATA/BioplexLysates_LUMINEX/BIOPLEX_R5/171012_Zeba_XX_5_lxb",
                   "../RAW_DATA/BioplexLysates_LUMINEX/BIOPLEX_R5/171012_Zeba_XO_5_lxb" )


for(lxb_folder in R5_lxb_folders){
    if (!grepl("/$", lxb_folder)) {
        lxb_folder = paste0(lxb_folder, "/")
    }
    
    do_plots = TRUE
    
    full_data = readLxb(paste0(lxb_folder, "*.lxb"), text=FALSE)
    lxb_name = gsub("[0-9]{6}_", "", lxb_folder)
    lxb_name = gsub("_lxb", "", gsub("/$", "", lxb_name))
    
    print("Processing bead data")
    beads = readBeads(full_data, region, analyte, esetup$controls, esetup$blanks, esetup$wptb, lxb_folder, esetup$externals, do_plots)
    print("Bead data processed")
    if (file.exists("post_processing.R")) {
        print("Post processing the data")
        source("post_processing.R")
    }
    lxb_name = basename(lxb_name)
    writeMIDAS(beads, lxb_name, esetup$inhibitors, esetup$stimulators, lxb_name)
    print(paste0("Data writen in MIDAS format in folder ", lxb_name))
    
    pdf(paste0(lxb_name, "/", "beads_distribution_", lxb_name, ".pdf"), width=10)
    analyseBeadDistributions(full_data, region, analyte, esetup$treatments)
    dev.off()
    
}

# CREATING the OUTPUT Folder for all MIDAS files and copying the files there
print("CREATING the OUTPUT Folder for all MIDAS files and copying the files there")

dir.create("./OUTPUT_MIDASfiles")

# identify the folders
all_folders <- c("Zeba_XO3", "Zeba_XO4", "Zeba_XO_5", "Zeba_XX3", "Zeba_XX4", "Zeba_XX_5")
new_folder <- "./OUTPUT_MIDASfiles"


for (current_folder in all_folders) {
    # find the file that you want
    MIDAS_file <- list.files(current_folder, ".csv$")
    
    # copy the files to the new folder
    file.copy(file.path(current_folder,MIDAS_file), new_folder)
}

#Renaming all MIDAS files for consistency 
setwd("./OUTPUT_MIDASfiles")

midas_files <- list.files("./")
new_names <- sub("blunt_Zeba_", "", midas_files)
file.rename(from = midas_files, to = new_names)

R5_files <- c( "XX_5_MIDAS.csv", "XO_5_MIDAS.csv")
R5_new_names <- sub("_5", "5", R5_files)
file.rename(from = R5_files, to = R5_new_names)

setwd("../")


print("All data extracted as MIDAS files and renamed appropriately")

#Removing Temp folders not needed
Temp_Folders <- list.files(path=".",pattern = "Zeba_")
unlink(Temp_Folders, recursive = TRUE)

print("Temporary folders deleted")

print("Script_1A : All steps executed")








