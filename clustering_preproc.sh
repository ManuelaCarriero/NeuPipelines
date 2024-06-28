#!/bin/bash

source project.sh

mapfile -t SCANS < allnewsubjs.txt # read subjects from text file


function mask_mp2rage() {
  #bet INV2
  bet $INV2 brain040_inv2.nii.gz -R -f 0.40 -g 0 -n -m
  
  echo "masking UNI with INV2"
  
  #uni brain masked with inv2
  fslmaths $UNI -mas ${subj}_brain040_inv2.nii.gz ${subj}_t1_mp2rage_UNIMaskedwithINV2.nii.gz
}

function register_DWIMaps() {
  echo "registering b0 to T1 weighted image"
  
  fslroi $b0_topup/my_hifi_b0.nii.gz  first_hifi_b0.nii.gz 0 1
  
  #register diffusion to T1 weighted image 
  flirt -in first_hifi_b0.nii.gz -ref t1_mp2rage_UNIMaskedwithINV2.nii.gz -out b0_to_T1.nii.gz -omat b0_to_T1.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp $interp_method
  
  
  #register T1 weighted image to MNI 
  echo "registering T1 weighted image to MNI"
  bet ${subj}_t1_mp2rage_UNIMaskedwithINV2.nii.gz ${subj}_t1_mp2rage_UNIMaskedwithINV2_brain.nii.gz -R -f 0.40 -g 0 -n -m
  
  flirt -in ${subj}_t1_mp2rage_UNIMaskedwithINV2_brain.nii.gz -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -dof 12 -out T1toMNIlin -omat T1toMNIlin.mat
  
  
  fnirt --in=${subj}_t1_mp2rage_UNIMaskedwithINV2.nii.gz --aff=T1toMNIlin.mat --config=T1_2_MNI152_2mm.cnf --iout=T1toMNInonlin --cout=T1toMNI_coef --fout=T1toMNI_warp
   
   
   
   #Concat FLIRT transform
   echo "concat matrix of first trasformation with the second matrix"
   
   #Apply FLIRT transform
   echo "register dwi maps to MNI"

} 
 
 for SUBJ in ${SCANS[@]}; do
   echo "Processing: ${SUBJ}"
   # folder definitions

   clustering_DIR="${PROC_DIR}/${SUBJ}/dwi/clustering"

   mkdir -p $clustering_DIR

   cd $clustering_DIR
   
   #set paths to filenames
   INV2="/storage/shared/PRINAntonello2022/BIDS/derivatives/${subj}/anat/${subj}_mp2rage_INV2_brain.nii.gz"
   UNI="/storage/shared/PRINAntonello2022/BIDS/derivatives/${subj}/anat/${subj}_mp2rage_UNI_cropped.nii.gz"  
  #DWI="/storage/shared/SANDI_240229+Results/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240418tr3000woPF2/rawdata/cmrr_mbep2d_diff_7shell_tr3000_20240418155310_6.nii.gz" #without  
   interp_method='nearestneighbour'  
   b0_topup="/storage/shared/PRINAntonello2022/BIDS/derivatives/${subj}/dwi/preprocessed/${subj}_my_hifi_b0.nii.gz"   
   FSLDIR="/usr/local/fsl_v6062/data/standard"
   
   mask_mp2rage
    
   register_DWImaps 
   

   
   
done