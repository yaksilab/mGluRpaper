mGluR Paper
-----------  
Code for "Inhibition mediated by group III mGluRs regulates habenula activity and defensive behaviors"  

Email: nick-dot-faturos-@-gmail-dot-com

Overview  
--------  
This package contains MATLAB scripts and experimental data for analyzing and visualizing figures 2E-F and 5D-H. 
The statistical tests for each figure are performed and stored within each figure script. 

Usage Instructions  
------------------  
- *Update* the figure scripts folder name to wherever this package is currently located on your computer.
- Add the general path of '.../code' to the MATLAB path.
- Use the scripts in 'figure scripts/' to generate specific figures from the processed data.  
- 'process_cppg_acsf_mV.m' and 'process_stim_data.m' pull from the raw data not provided here, but shows how the data was organized.

Organization  
------------  
/["Base Directory"]  
  ├── README.txt                    # Documentation file (this file)  

  ├── code/                          # MATLAB scripts for data processing and analysis  
  │   ├── ephys/                     	# Scripts for organizing ephys data  
  │   │   ├── process_cppg_acsf_mV.m 	# Organizes electrophysiology data from raw files  
  │   ├── figure scripts/            # Scripts for generating figures in the paper  
  │   │   ├── figure_acsf_cdld2_lap4_mv.m         # Compares membrane potential voltage b/w acsf/cdld2/lap4
  │   │   ├── figure_acsf_cell_responses.m  	    # Shows the calcium traces of the cells around the patched neuron
  │   │   ├── figure_acsf_exc_and_inb_distance.m  # Compares the distances b/w excited and inhibited cells in aCSF
  │   │   ├── figure_acsf_exc_and_inb_percent.m   # Compares the percent of excited and inhibited cells in aCSF 
  │   │   ├── figure_cppg_acsf_percent_response.m # Compares the percent of excited and inhibited cells b/w aCSF and CPPG 
  │   │   ├── figure_cppg_and_acsf_mv.m  	  # Compares membrane potential voltage b/w acsf/cppg
  │   │   ├── figure_cppg_exc_and_inb_distance.m  # Compares the distances b/w excited and inhibited cells in CPPG
  │   │   ├── figure_cppg_exc_and_inb_percent.m   # Compares the percent of excited and inhibited cells in aCSF 
  │   ├── imaging/                   # Scripts for processing raw imaging data  
  │   │   ├── process_stim_data.m     # Main script for imaging data organization  
  │   │   ├── get_stim_responses.m    # Fn that organizes the stimulation data and analyzes a cell's response
  │   │   ├── load_stim_data.m        # Fn that loads the imaging data
  │   ├── utilities/                   # General utility functions for plotting and statistics  
  │   │   ├── add_significance_line.m   # Adds a line and asterisk(s) based on pvalue
  │   │   ├── align_yyaxis_zero.m  
  │   │   ├── compute_ranksum_stats.m   # Exports all of the associated test statistics
  │   │   ├── compute_signrank_stats.m  # Exports all of the associated test statistics
  │   │   ├── init_figure.m  			  # Opens a figure with set fonts and sizes (see fn notes)
  │   │   ├── plot_density_scatter.m  

  ├── experiments/                   # Folders containing the data used to generate figures  
  │   ├── ACSF stim imaging/         # Imaging data collected under ACSF conditions  
  │   │   ├── data.mat               # Processed imaging data  
  │   ├── ACSF+CPPG ephys/           # Electrophysiology data collected under ACSF+CPPG conditions  
  │   │   ├── data.mat               # Processed ephys data  
  │   ├── CPPG stim imaging/         # Imaging data collected under CPPG conditions  
  │   │   ├── data.mat               # Processed imaging data  
  │   ├── CdCl2+L-AP4 ephys/         # Electrophysiology data under CdCl2+L-AP4 conditions  
  │   │   ├── data.mat               # Processed ephys data 