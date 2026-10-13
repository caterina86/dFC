#!/bin/bash

# Bandpass or high-pass filtering of the functional data

# fslmaths preproc_data -bptf <hp_sigma> <lp_sigma> filtered_preproc_data

# Recommendation: <hp_sigma> = 1/(2*f*TR), <lp_sigma> = 1/(18*f*TR); where f is frequency in Hz, 
# TR is the number of seconds between volumes, and the numbers 2 and 18 are the 'magic numbers' 
# related to keeping 95% of signal in these frequencies


## Dataset 2023

baseDir=$1

inputDir=${baseDir}/dataset2023
outputDir=${baseDir}/temporal_filtering # define the output directory
value=`cat ${inputDir}/list_subjects_dataset2023.txt`

# Define the sigma values using the given formula
hp_sigma=$(echo "scale=10; 1 / (2 * 0.01 * 1)" | bc) # this is correct
lp_sigma=$(echo "scale=10; 1 / (18 * 0.1 * 1)" | bc)  # this is correct

# Print the sigma values
echo "High-pass sigma value: $hp_sigma"
echo "Low-pass sigma value: $lp_sigma"

for subj in $value; do

	filename="${outputDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func_data_clean_0.01_0.1HZ.nii.gz"
	
	if [ -e "$filename" ]; then
    		echo "File $filename exists."
	else
    		echo "$subj"
    		mkdir ${outputDir}/derivatives/ica/sub-${subj}
    		fslmaths ${inputDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func.ica/sub-${subj}_filtered_func_data_clean.nii.gz -bptf ${hp_sigma} ${lp_sigma} ${outputDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func_data_clean_0.01_0.1HZ.nii.gz
	
	fi

done



## Dataset 2017

inputDir=${baseDir}/dataset2017
value=`cat ${inputDir}/list_subjects_dataset2017.txt`

# Define the sigma values using the given formula
hp_sigma=$(echo "scale=10; 1 / (2 * 0.01 * 0.785)" | bc) # this is correct
lp_sigma=$(echo "scale=10; 1 / (18 * 0.1 * 0.785)" | bc)  # this is correct

# Print the sigma values
echo "High-pass sigma value: $hp_sigma"
echo "Low-pass sigma value: $lp_sigma"

for subj in $value; do
    
	filename="${outputDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func_data_clean_0.01_0.1HZ.nii.gz"
	
	if [ -e "$filename" ]; then
    		echo "File $filename exists."
	else
	    echo "$subj"
    		mkdir ${outputDir}/derivatives/ica/sub-${subj}
    		fslmaths ${inputDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func.ica/sub-${subj}_filtered_func_data_clean.nii.gz -bptf ${hp_sigma} ${lp_sigma} ${outputDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func_data_clean_0.01_0.1HZ.nii.gz

	fi

done





