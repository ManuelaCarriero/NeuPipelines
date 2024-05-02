%Program to estimate partial volume effect in images with Partial Fourier
%and without Partial Fourier

%Insert the filename if you are in the same folder of the data.

%Rsoma map to T1
img_path_pve1 = '/storage/shared/SANDI_240418/NIFTI/t1_mp2rage_sag_p3_iso_fast/Segmented/t1_mp2rage_UNIMaskedwithINV2_pve_1.nii.gz';
img_path_pve0 = '/storage/shared/SANDI_240418/NIFTI/t1_mp2rage_sag_p3_iso_fast/Segmented/t1_mp2rage_UNIMaskedwithINV2_pve_0.nii.gz';
img_path_pve2 = '/storage/shared/SANDI_240418/NIFTI/t1_mp2rage_sag_p3_iso_fast/Segmented/t1_mp2rage_UNIMaskedwithINV2_pve_2.nii.gz';
img_path_map = '/storage/shared/SANDI_240418/NIFTI/RsomawPF_to_T1.nii.gz';
title_plot = 'Rsoma map With Partial Fourier';

%T1 to diffusion weighted image
img_path_pve1 = '/storage/shared/SANDI_240418/NIFTI/pve1_to_DWIwoPF.nii.gz';
img_path_pve0 = '/storage/shared/SANDI_240418/NIFTI/pve0_to_DWIwoPF.nii.gz';
img_path_pve2 = '/storage/shared/SANDI_240418/NIFTI/pve2_to_DWIwoPF.nii.gz';
img_path_map = '/storage/shared/SANDI_240229/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240418tr3000woPF/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_Rsoma.nii.gz';
title_plot = 'Rsoma map Without Partial Fourier';

Vhdr_map = spm_vol(img_path_map);
V_map = spm_read_vols(Vhdr_map);

thr_pve = 0.1:0.1:0.9;
mean_par = [];
Vhdr_pve1 = spm_vol(img_path_pve1);
pve1 = spm_read_vols(Vhdr_pve1);
Vhdr_pve0 = spm_vol(img_path_pve0);
pve0 = spm_read_vols(Vhdr_pve0);
Vhdr_pve2 = spm_vol(img_path_pve2);
pve2 = spm_read_vols(Vhdr_pve2);
for value = thr_pve


    %


    %figure,imagesc(V_pve(:,:,150))

    %make binary mask
    GM=pve1>value;
    WM=pve2>value;
    CSF=pve0>value;

    %figure,imagesc(V_pve(:,:,150))

    %

    V_GM = GM.*V_map;
    V_WM = WM.*V_map;
    V_CSF = CSF.*V_map;
    %figure,imagesc(V_map_masked(:,:,150));

    %calculate mean of parameter value in one region of brain cortex
    %in a second moment, choose a precise region (like the temporal one)
    %
    %if value == 0.1 || value == 0.5 || value == 0.9
        %
        %figure,imagesc(V_WM(1:176,1:256,150));
        %
        %
        %
        %
    %end

    V_GM(V_GM==0)=NaN;
    V_WM(V_WM==0)=NaN;
    V_CSF(V_CSF==0)=NaN;
    Vmean(1) = nanmean(V_GM(:));%Compute the arithmetic mean along the specified axis, ignoring NaNs.
    Vmean(2) = nanmean(V_WM(:));
    Vmean(3) = nanmean(V_CSF(:));
    mean_par(end+1,:) = Vmean;

end

figure, plot(thr_pve, mean_par);
%hold on
%plot(thr_pve, mean_par_wo);
xlabel('Threashold value');
ylabel('Mean value of Rsoma');
title(title_plot);





