#!/bin/bash

INV2='/storage/shared/PRINAntonello2022/BIDS/sub-03/anat/sub-03_t1_mp2rage_sag_p3_iso_fast_INV2.nii.gz'
UNI='/storage/shared/PRINAntonello2022/BIDS/sub-03/anat/sub-03_t1_mp2rage_sag_p3_iso_fast_UNI_Images.nii.gz'

DWI='/storage/shared/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240523/rawdata/sub-03_cmrr_mbep2d_diff_7shell_tr3000.nii.gz'

PVE_FOLDER='/media/nas_rete/Work_manuela/PVE240523'

interp_method='nearestneighbour'

b0_topup='/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240523/derivatives'


bet $INV2 $PVE_FOLDER/brain040_inv2.nii.gz -R -f 0.40 -g 0 -n -m

#uni brain masked with inv2
fslmaths $UNI -mas $PVE_FOLDER/brain040_inv2.nii.gz $PVE_FOLDER/t1_mp2rage_UNIMaskedwithINV2.nii.gz

echo "starting image segmentation"

#image segmentation
fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o $PVE_FOLDER/t1_mp2rage_UNIMaskedwithINV2_seg $PVE_FOLDER/t1_mp2rage_UNIMaskedwithINV2

fslroi $b0_topup/my_hifi_b0.nii.gz  $PVE_FOLDER/first_hifi_b0.nii.gz 0 1
echo "registering pve to diffusion image"

#register diffusion to T1 weighted image 
flirt -in $PVE_FOLDER/first_hifi_b0.nii.gz -ref $PVE_FOLDER/t1_mp2rage_UNIMaskedwithINV2.nii.gz -out $PVE_FOLDER/b0_to_T1.nii.gz -omat $PVE_FOLDER/b0_to_T1.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp $interp_method

echo "registered pve to diffusion image. Now inverting matrix."

#invert matrix (without Partial Fourier)
convert_xfm -omat $PVE_FOLDER/inverse_b0_to_T1.mat -inverse $PVE_FOLDER/b0_to_T1.mat

#apply transformation matrix to labeled image (without Partial Fourier)
echo "apply transformation matrix to labeled image"


#pve0
flirt -in $PVE_FOLDER/t1_mp2rage_UNIMaskedwithINV2_seg_pve_0.nii.gz -ref $PVE_FOLDER/first_hifi_b0.nii.gz -out $PVE_FOLDER/pve0_to_b0.nii.gz -omat $PVE_FOLDER/pve0_to_b0.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp $interp_method

#pve1
flirt -in $PVE_FOLDER/t1_mp2rage_UNIMaskedwithINV2_seg_pve_1.nii.gz -ref $PVE_FOLDER/first_hifi_b0.nii.gz -out $PVE_FOLDER/pve1_to_b0.nii.gz -omat $PVE_FOLDER/pve1_to_b0.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp $interp_method

#pve2
flirt -in $PVE_FOLDER/t1_mp2rage_UNIMaskedwithINV2_seg_pve_2.nii.gz -ref $PVE_FOLDER/first_hifi_b0.nii.gz -out $PVE_FOLDER/pve2_to_b0.nii.gz -omat $PVE_FOLDER/pve2_to_b0.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp $interp_method





