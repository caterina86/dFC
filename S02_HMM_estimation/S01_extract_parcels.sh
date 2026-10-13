#!/bin/bash

baseDir=$1

mkdir ${baseDir}/osf/derivatives/parcels

# Extract the parcels -> do it for all the parcels (100)

echo "extract the parcels"

for i in {1..100}; do

	if [ ! -f "${baseDir}/osf/derivatives/parcels/parcel$i.nii.gz" ]; then

        fslmaths ${baseDir}/osf/derivatives/Atlas/MNI/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz -thr $i -uthr $i -bin ${baseDir}/osf/derivatives/parcels/parcel$i.nii.gz

	else

        echo "The file already exists. Skipping the loop."

	fi

done

