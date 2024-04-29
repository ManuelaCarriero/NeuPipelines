%Program to estimate partial volume effect in images with Partial Fourier
%and without Partial Fourier

%Insert the filename if you are in the same folder of the data.
img_path_pve = '/storage/shared/SANDI_240418/NIFTI/t1_mp2rage_sag_p3_iso_fast/Segmented/t1_mp2rage_MRF01UNIMaskedwithINV2_pve_1.nii.gz';
img_path_map = '/storage/shared/SANDI_240418/NIFTI/Rsoma_to_T1.nii.gz';
title_plot = 'Rsoma map Without Partial Fourier';


Vhdr_map = spm_vol(img_path_map);
V_map = spm_read_vols(Vhdr_map);

thr_pve = 0:0.1:1;
mean_par = [];

for value = thr_pve

    Vhdr_pve = spm_vol(img_path_pve);
    V_pve = spm_read_vols(Vhdr_pve);

    %

    size(V_pve)

    %figure,imagesc(V_pve(:,:,150))

    %make binary mask
    V_pve(V_pve<value) = 0;
    V_pve(V_pve>value) = 1;

    %figure,imagesc(V_pve(:,:,150))

    %

    V_map_masked = V_pve.*V_map;

    %figure,imagesc(V_map_masked(:,:,150));

    %calculate mean of parameter value in one region of brain cortex
    %in a second moment, choose a precise region (like the temporal one)

    if value == 0 || value == 0.5 || value == 1

        figure,imagesc(V_map_masked(100:150,120:180,150));




    end

    V_mm_ROI = V_map_masked(100:150,120:180,150);
    Vmean = mean(V_mm_ROI,"all");

    mean_par(end+1) = Vmean;

end

figure, plot(thr_pve, mean_par);
xlabel('Threashold value');
ylabel('Mean value of Rsoma');
title(title_plot);



