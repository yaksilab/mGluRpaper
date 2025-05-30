# mGluRpaper

Collection of the code used for analysis and plotting of the figures for the Ostenrath et al 2024 titled "Inhibition mediated by group III mGluRs regulates habenula activity and defensive behaviors"
Each Figure has it's own matlab script and the implemented functions are collected in the "Subfunction" folder. 


## Requirements: 
-	All code was written in MATLAB and was tested on MATLAB version 2022a, 2023a and 2024a
-	Code runs on the standard MATLAB version (https://www.mathworks.com/products/matlab.html) please follow the instruction to install it here
-	Any extra function are part of the code delivered and is found in the the script or the "Subfunction" folder

## General instructions
-	Each figure has its own code script with instruction on how to run the code for each figure/experiment.
-	The data path (plus other required paths) needs to be replaced for your system and then the code can be run either as a whole or per section. Each code has more information in the actual script file.
-	Due to the size of the data, no demo data is included. To get access to the data please check the publication.
-	Details are provided in the method section of the manuscript. 
## Specific Instructions: 
### Fig 1:
  - This script will generate the Figure 1 and Suppl. Figure 1. 
  - Load in the corresponding data / replace the data path and replace the save path for your system
  -	When running the script each subfigure will be generated and saved.
### Fig 2E-F and Fig 5D-H
- Please refer to the READ ME file in the subfolder.
### Fig 3
  -	This script with generate Figure 3 and Suppl. Figure 3A-C.
  -	Load in the corresponding data / replace the data path and replace the save path for your system (note: there is a data file for each experiment)
### Fig 4, 5A-C and 6
  -	This script will generate Figure 4, 5 A-C and 6 as well as Suppl. Figure 4 – 7 and Suppl. Figure 9
  -	Load in the corresponding data / replace the data path and replace the save path for your system (note: for Suppl. Figure 9 it is a different data set)
  -	For the main figures the code can be run as is (using brain_region = 11 for Hb). 
  -	Suppl. Figure 4 will also be automatically generated as the non-injected fish were added as group 1 in the data. 
  -	For Suppl. Figure 5, the script “finding_midbrain” found in the “Subfunction” folder needs to be run and then the “brainnumber” variable replace with 15. 
  -	For Suppl. Figure 6, the brain region variable is replaced with a list indicating dorsal (1) and ventral (2) in Hb. Run the code snippet indicating this and then replace the “brainnumber” variable with the either 1 or 2 depending on which neuron group you   want to focus. 
  -	For Suppl. Figure 7. Replace the “std_factor” variable in the responding cells section with the number you are interested in (e.g. 1,3 or 4) and run the rest as normal (note: brainnumber here is equal to 11 indicating all of Hb)
  -	For Suppl. Figure 9. Load the dataset belonging to the mGluR mutant experiment and use the “brainnumber” 11 to focus on Hb. 
### Fig 7
  -	The scripts are divided into the different experiments 
  -	Load in the data for the corresponding experiment and replace the paths (e.g. save_path etc. as indicated in the script). 
### Suppl. Fig 10
  -  Load in the data for the corresponding experiment and replace the paths (e.g. save_path etc. as indicated in the script). 
-	The expected output for each script/ section can be seen as a figure panel in the manuscript or supplementary information. 

