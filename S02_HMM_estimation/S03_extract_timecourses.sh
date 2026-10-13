#!/bin/bash

# Extract the timecourse from each parcel

baseDir=$1
atlasDir=$2
projectDir=$3
filter=$4

value=`cat ${baseDir}/list_subjects_dataset2023_2017.txt`

echo "$value"

for subj in $value; do
    
    
	if [ ! -f "${atlasDir}/sub-${subj}/sub-${subj}_filtered_func_data_clean_${filter}HZ_atlas.nii.gz" ]; then
	
		mkdir ${atlasDir}/sub-${subj}/

		echo "coregister functional to MNI space for subject ${subj}"

		fslmaths ${projectDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func.ica/sub-${subj}_filtered_func_data_clean.nii.gz -Tmean ${projectDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func.ica/sub-${subj}_filtered_func_data_clean_mean.nii.gz

		flirt -in ${projectDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func.ica/sub-${subj}_filtered_func_data_clean_mean.nii.gz -ref /Users/cp3488/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -omat ${projectDir}/derivatives/atlas/sub-${subj}/atlas2func.mat

		flirt -in ${projectDir}/derivatives/ica/sub-${subj}/sub-${subj}_filtered_func_data_clean_${filter}HZ.nii.gz -ref /Users/cp3488/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz -out ${atlasDir}/sub-${subj}/sub-${subj}_filtered_func_data_clean_${filter}HZ_atlas.nii.gz -init /Users/cp3488/Documents/Projects/dFC/dataset2023/derivatives/atlas/sub-${subj}/atlas2func.mat -applyxfm

	else
        echo "The file already exists. Skipping the loop."
	fi

	echo "extract timecourse from parcels"
	mkdir ${atlasDir}/sub-${subj}/timecourses/

	for i in {1..100}; do

	if [ ! -f "${atlasDir}/sub-${subj}/timecourses/timecourse_parcel$i.txt" ]; then

		fslmeants -i ${atlasDir}/sub-${subj}/sub-${subj}_filtered_func_data_clean_${filter}HZ_atlas.nii.gz -m ${projectDir}/derivatives/parcels/parcel$i.nii.gz -o ${atlasDir}/sub-${subj}/timecourses/timecourse_parcel$i.txt --eig

	fi

	done


done



