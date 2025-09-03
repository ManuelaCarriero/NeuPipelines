%% load data
% load subjs idx
run='run-01';
subjects = importdata(strcat('/media/nas_rete/Vitality/code/subjs_DWI.txt'));
n_subjs=length(subjects);
start_subj=1;
%%
%load pve on subj space
V_pves_dwi={};
V_pves_func={};

for i = start_subj:n_subjs
    
    subj = subjects{i};

    img_path_pve_dwi=strcat('/media/nas_rete/Vitality/maps2SUBJSPACE/dwi/pve_on_b0/',subj,'_',run,'_PVE_1_on_b0.nii.gz');
    img_path_pve_func=strcat('/media/nas_rete/Vitality/maps2SUBJSPACE/func/pve_on_M0/',subj,'_',run,'_PVE_1_on_M0.nii.gz');

    V_vol_pve_dwi = spm_vol(img_path_pve_dwi);
    V_pve_dwi=spm_read_vols(V_vol_pve_dwi);
    V_pves_dwi{end+1}=V_pve_dwi;

    V_vol_pve_func = spm_vol(img_path_pve_func);
    V_pve_func=spm_read_vols(V_vol_pve_func);
    V_pves_func{end+1}=V_pve_func;

end

%% plot pves to check

% V_pves_0_dwi_first=V_pves_0_dwi{1};
% % figure, imagesc(V_pves_0_dwi_first(:,:,45))
% 
% thr=0.1;
% V_pves_0_dwi_first(V_pves_0_dwi_first<thr)=1;
% V_pves_0_dwi_first(V_pves_0_dwi_first<1)=0;
% figure, imagesc(V_pves_0_dwi_first(:,:,45))
% title(num2str(thr))
% 
% V_pves_1_dwi_first=V_pves_1_dwi{1};
% V_pves_1_dwi_first(V_pves_1_dwi_first>0.5)=1;
% V_pves_1_dwi_first(V_pves_1_dwi_first<1)=0;
% % figure, imagesc(V_pves_1_dwi_first(:,:,45))
% 
% V_refined = V_pves_0_dwi_first.*V_pves_1_dwi_first;
% figure, imagesc(V_refined(:,:,45))
%%
%load atlas on subj space
V_atlases_dwi={};
V_atlases_func={};

for i = start_subj:n_subjs
    
    subj = subjects{i};
    
    img_path_atlas_dwi=strcat('/media/nas_rete/Vitality/maps2SUBJSPACE/dwi/atlas_on_b0/AAL3v1_2mm_on_',subj,'_',run,'_acq-dwi_B0_brain_corr.nii.gz');
    img_path_atlas_func=strcat('/media/nas_rete/Vitality/maps2SUBJSPACE/func/atlas_on_M0/AAL3v1_2mm_on_',subj,'_task-bh_',run,'_acq-dexi_M0.nii.gz');

    V_vol_atlas_dwi = spm_vol(img_path_atlas_dwi);
    V_atlas_dwi=spm_read_vols(V_vol_atlas_dwi);
    V_atlases_dwi{end+1}=V_atlas_dwi;

    V_vol_atlas_func = spm_vol(img_path_atlas_func);
    V_atlas_func=spm_read_vols(V_vol_atlas_func);
    V_atlases_func{end+1}=V_atlas_func;

end

%% load parametric maps
% V_CMRO2_maps={};
% V_CBF_maps={};

V_rsoma_maps={};
V_fsoma_maps={};
V_fneurite_maps={};
V_De_maps={};
V_Din_maps={};
V_fextra_maps={};

V_mse_maps={};

for i = start_subj:n_subjs
    
    subj = subjects{i};
    
%     img_path_CMRO2=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/perf/outcome/',subj,'_task-bh_',run,'_acq-dexi_volreg_asl_topup_CMRO2_map.nii.gz');
%     img_path_CBF=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/perf/outcome/',subj,'_task-bh_',run,'_acq-dexi_volreg_asl_topup_CBF_map.nii.gz');

%     img_path_rsoma=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_Rsoma.nii.gz');
%     img_path_fsoma=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_fsoma.nii.gz');
%     img_path_fneurite=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_fneurite.nii.gz');
%     img_path_De=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_De.nii.gz');
%     img_path_Din=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_Din.nii.gz');
%     img_path_fextra=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_fextra.nii.gz');
%     img_path_mse=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_mse.nii.gz');
    img_path_rsoma=strcat('/media/nas_rete/Work_manuela/Vitality_data_SANDInewrelease/',subj,'/SANDI_output/',subj,'_',run,'_SANDI-fit_Rsoma.nii.gz');
    img_path_fsoma=strcat('/media/nas_rete/Work_manuela/Vitality_data_SANDInewrelease/',subj,'/SANDI_output/',subj,'_',run,'_SANDI-fit_fsoma.nii.gz');
    img_path_fneurite=strcat('/media/nas_rete/Work_manuela/Vitality_data_SANDInewrelease/',subj,'/SANDI_output/',subj,'_',run,'_SANDI-fit_fneurite.nii.gz');
    img_path_De=strcat('/media/nas_rete/Work_manuela/Vitality_data_SANDInewrelease/',subj,'/SANDI_output/',subj,'_',run,'_SANDI-fit_De.nii.gz');
    img_path_Din=strcat('/media/nas_rete/Work_manuela/Vitality_data_SANDInewrelease/',subj,'/SANDI_output/',subj,'_',run,'_SANDI-fit_Din.nii.gz');
    img_path_fextra=strcat('/media/nas_rete/Work_manuela/Vitality_data_SANDInewrelease/',subj,'/SANDI_output/',subj,'_',run,'_SANDI-fit_fextra.nii.gz');
    img_path_mse=strcat('/media/nas_rete/Work_manuela/Vitality_data_SANDInewrelease/',subj,'/SANDI_output/',subj,'_',run,'_SANDI-fit_mse.nii.gz');

%     V_vol_CMRO2 = spm_vol(img_path_CMRO2);
%     V_CMRO2=spm_read_vols(V_vol_CMRO2);
%     V_CMRO2_maps{end+1}=V_CMRO2;
% 
%     V_vol_CBF = spm_vol(img_path_CBF);
%     V_CBF=spm_read_vols(V_vol_CBF);
%     V_CBF_maps{end+1}=V_CBF;

    V_vol_rsoma = spm_vol(img_path_rsoma);
    V_rsoma=spm_read_vols(V_vol_rsoma);
    V_rsoma_maps{end+1}=V_rsoma;

    V_vol_fsoma = spm_vol(img_path_fsoma);
    V_fsoma=spm_read_vols(V_vol_fsoma);
    V_fsoma_maps{end+1}=V_fsoma;

    V_vol_fneurite = spm_vol(img_path_fneurite);
    V_fneurite=spm_read_vols(V_vol_fneurite);
    V_fneurite_maps{end+1}=V_fneurite;

    V_vol_De = spm_vol(img_path_De);
    V_De=spm_read_vols(V_vol_De);
    V_De_maps{end+1}=V_De;

    V_vol_Din = spm_vol(img_path_Din);
    V_Din=spm_read_vols(V_vol_Din);
    V_Din_maps{end+1}=V_Din;

    V_vol_fextra = spm_vol(img_path_fextra);
    V_fextra=spm_read_vols(V_vol_fextra);
    V_fextra_maps{end+1}=V_fextra;

    V_vol_mse = spm_vol(img_path_mse);
    V_mse=spm_read_vols(V_vol_mse);
    V_mse_first = V_mse(:,:,:,1);
    V_mse_maps{end+1}=V_mse_first;

end



%% select binary masks thresholds
pve_1_threshold=0.5;
pve_0_threshold=0.4;%threshold beyond which the CSF is set to 0
pve_2_threshold=0.4;%
fsoma_threshold=0.15;

%% number of cells density map (many subjects) 

V_fc_maps={};

for i = 1:n_subjs    
    V_rsoma = V_rsoma_maps{i};
    V_fsoma = V_fsoma_maps{i};
    V_pve = V_pves_dwi{i};
    
    for i = 1:length(V_rsoma_maps(:))
        if V_rsoma(i) < 4
            V_rsoma(i)=0;
        end
    end
    
    V_pve(V_pve>pve_threshold)=1;
    V_pve(V_pve<1)=0;
    V_fsoma_to_mask=V_fsoma;
    V_fsoma_to_mask(V_fsoma_to_mask>fsoma_threshold)=1;
    V_fsoma_to_mask(V_fsoma_to_mask<1)=0;
    V_rsoma_masked = V_rsoma.*V_fsoma_to_mask.*V_pve;
    
    %convert to m^3
    V_rsoma_masked = V_rsoma_masked.*10^-6;
    
    %voxel wise divide fs map over 4/3pir^3
    fc_map = V_fsoma./((4/3)*pi*V_rsoma_masked.^3);
    %figure, imagesc(fc_map(:,:,45));
    for i = 1:length(fc_map(:))
        if fc_map(i)==Inf
            fc_map(i)=NaN;%VALUTA SE METTERE A 0.
    %     elseif isnan(fc_map(i))
    %         fc_map(i)=0;
        end
    end

    V_fc_maps{end+1} = fc_map;

end

% %find a method to remove too much high values at the borders.
% fc_1=V_fc_maps{1}*10^(-9);
% figure, imagesc(rot90(fc_1(:,:,40)))
% title('Numerical soma density')
% clim([0,100000])

%For the remaining unresonable values, you could try MSE threshold.
%% superficial density (many subjects)

V_fsup_maps={};

for i = 1:n_subjs    
    V_rsoma = V_rsoma_maps{i};
    V_fsoma = V_fsoma_maps{i};
    V_pve = V_pves_dwi{i};
    
    for i = 1:length(V_rsoma_maps(:))
        if V_rsoma(i) < 4
            V_rsoma(i)=0;
        end
    end
    
    V_pve(V_pve>pve_threshold)=1;
    V_pve(V_pve<1)=0;
    V_fsoma_to_mask=V_fsoma;
    V_fsoma_to_mask(V_fsoma_to_mask>fsoma_threshold)=1;
    V_fsoma_to_mask(V_fsoma_to_mask<1)=0;
    V_rsoma_masked = V_rsoma.*V_fsoma_to_mask.*V_pve;

    %convert to m^3
    V_rsoma_masked = V_rsoma_masked.*10^-6;    

    %voxel wise divide fs map over 4/3pir^3
    fsup_map = V_fsoma./V_rsoma_masked;
    fsup_map = 3*fsup_map;
    %figure, imagesc(fc_map(:,:,45));
    for i = 1:length(fsup_map(:))
        if fsup_map(i)==Inf
            fsup_map(i)=NaN;%VALUTA SE METTERE A 0.
    %     elseif isnan(fc_map(i))
    %         fc_map(i)=0;
        end
    end

    V_fsup_maps{end+1} = fsup_map;

end

V_fsup_one=V_fsup_maps{1}.*10^(-6);
figure, imagesc(rot90(V_fsup_one(:,:,40)));
title('Superficial Soma Density map')
% % V_fsup_one_array=V_fsup_one(:);
% % figure, hist(V_fsup_one_array);
% % title('Superficial Soma Density Distribution');
% % grid on

%%  Count total number of regions
img_path_atlas='/storage/shared/Atlas/AAL3v1_2mm_resampled.nii.gz';
Vhdr = spm_vol(img_path_atlas);
V_atlas_tot = spm_read_vols(Vhdr);

regions = unique(V_atlas_tot(:));
%background removal
regions(1)=[];
n_regions=numel(regions);

%% in case you don't have the original atlas
n_regions=165;

%% compute medians
%prepare empty matrices with maximum size 
%so to not have problems of different sizze
%then take the common regions (intersection): 
%1. first keep the common regions across subjects;
%2. then keep the common regions across spaces (func and dwi)

labels_func_subjs = zeros(n_subjs,n_regions);
medians_CMRO2_subjs = zeros(n_subjs,n_regions);
medians_CBF_subjs = zeros(n_subjs,n_regions);

labels_dwi_subjs = zeros(n_subjs,n_regions);
medians_rsoma_subjs = zeros(n_subjs,n_regions);
medians_fsoma_subjs = zeros(n_subjs,n_regions);
medians_fc_subjs = zeros(n_subjs,n_regions);
medians_fsup_subjs = zeros(n_subjs,n_regions);

medians_fneurite_subjs = zeros(n_subjs,n_regions);
medians_fextra_subjs = zeros(n_subjs,n_regions);
medians_Din_subjs = zeros(n_subjs,n_regions);
medians_De_subjs = zeros(n_subjs,n_regions);

percentage_removal_CMRO2_subjs = zeros(n_subjs,n_regions);
percentage_removal_CBF_subjs = zeros(n_subjs,n_regions); 
percentage_nans_CMRO2_subjs = zeros(n_subjs,n_regions);
percentage_nans_CBF_subjs = zeros(n_subjs,n_regions); 
percentage_removal_microparameter_subjs = zeros(n_subjs,n_regions); 
percentage_removal_rsoma_subjs = zeros(n_subjs,n_regions);

medians_mse_subjs = zeros(n_subjs,n_regions);
medians_mse_rsoma_subjs = zeros(n_subjs,n_regions);
percentage_nans_rsoma_subjs = zeros(n_subjs,n_regions);%%%%

V_pve_0_func_median_subjs = zeros(n_subjs,n_regions);
V_pve_1_func_median_subjs = zeros(n_subjs,n_regions);
V_pve_2_func_median_subjs = zeros(n_subjs,n_regions);

V_pve_0_dwi_median_subjs = zeros(n_subjs,n_regions);
V_pve_1_dwi_median_subjs = zeros(n_subjs,n_regions);
V_pve_2_dwi_median_subjs = zeros(n_subjs,n_regions);

start_time=tic;
for subj = 1:n_subjs
    tic
    %load atlases
    V_atlas_func = V_atlases_func{subj};
    V_atlas_dwi = V_atlases_dwi{subj};

    %load pve maps
    V_pve_1_func = V_pves_1_func{subj};
    V_pve_1_dwi = V_pves_1_dwi{subj};
    V_pve_0_func = V_pves_0_func{subj};
    V_pve_0_dwi = V_pves_0_dwi{subj};
    V_pve_2_func = V_pves_2_func{subj};
    V_pve_2_dwi = V_pves_2_dwi{subj};

    %load parametric maps
    V_CMRO2 = V_CMRO2_maps{subj};
    V_CBF = V_CBF_maps{subj};

    V_rsoma = V_rsoma_maps{subj};
    V_fsoma = V_fsoma_maps{subj};
    V_fc = V_fc_maps{subj};
    V_fsup = V_fsup_maps{subj};

    V_fneurite = V_fneurite_maps{subj};
    V_Din = V_Din_maps{subj};
    V_De = V_De_maps{subj};
    V_fextra = V_fextra_maps{subj};
    V_MSE = V_mse_maps{subj};

    %where to save medians
    medians_CMRO2_subj = [];
    medians_CBF_subj = [];
    medians_rsoma_subj = [];
    medians_fsoma_subj = [];
    medians_fc_subj = [];
    medians_fsup_subj = [];
    medians_mse_subj = [];%
    medians_mse_rsoma_subj = [];%


    medians_fneurite_subj = [];
    medians_Din_subj = [];
    medians_De_subj = [];
    medians_fextra_subj = [];
    labels_func_subj = [];
    labels_dwi_subj = [];
    percentage_removal_CBF_subj = [];
    percentage_removal_CMRO2_subj = [];
    percentage_nans_CBF_subj = [];
    percentage_nans_CMRO2_subj = [];
    percentage_removal_microparameter_subj = [];
    percentage_removal_rsoma_subj = [];
    percentage_nans_rsoma_subj=[];%%%%

    V_pve_0_func_median_subj = [];
    V_pve_1_func_median_subj = [];
    V_pve_2_func_median_subj = [];

    V_pve_0_dwi_median_subj = [];
    V_pve_1_dwi_median_subj = [];
    V_pve_2_dwi_median_subj = [];

    %%%%FUNC SPACE

    V_pve_0_func_original = V_pve_0_func;
    V_pve_1_func_original = V_pve_1_func;
    V_pve_2_func_original = V_pve_2_func;

    %define binary mask which is different for each subject
    V_pve_1_func(V_pve_1_func>pve_1_threshold)=1;
    V_pve_1_func(V_pve_1_func<1)=0;

    %In case you want to mask by considering not all the WM
    V_pve_2_func(V_pve_2_func<pve_2_threshold)=NaN;    
    V_pve_2_func(V_pve_2_func>=pve_2_threshold)=0;
    V_pve_2_func(isnan(V_pve_2_func))=1;

    % %in case you want all the WM
    % V_pve_2_func(V_pve_2_func==0)=NaN;
    % V_pve_2_func(V_pve_2_func>0)=0;
    % V_pve_2_func(isnan(V_pve_2_func))=1;

    %In case you want to mask by considering not all the CSF
    V_pve_0_func(V_pve_0_func<pve_0_threshold)=NaN;
    V_pve_0_func(V_pve_0_func>=pve_0_threshold)=0;
    V_pve_0_func(isnan(V_pve_0_func))=1;

    % %in case you want all the CSF
    % V_pve_0_func(V_pve_0_func==0)=NaN;
    % V_pve_0_func(V_pve_0_func>0)=0;
    % V_pve_0_func(isnan(V_pve_0_func))=1;

    %loop over regions
    %list of regions (it can vary among subjs and among spaces)
    regions_func = unique(V_atlas_func(:));
    %background removal
    regions_func(1)=[];
    n_regions_func=numel(regions_func);

    V_atlas_mask_func=V_atlas_func; 
    for region = 1:n_regions_func
        %binarize atlas
        for ii = 1:length(V_atlas_func(:))
            if V_atlas_func(ii) == regions_func(region)
                V_atlas_mask_func(ii) = 1;
            else
                V_atlas_mask_func(ii) = 0;
            end
        end

    V_GM = V_pve_1_func.*V_pve_0_func.*V_pve_2_func;
    V_mask_func = V_GM.*V_atlas_mask_func;

    %mask

    V_CMRO2_masked = V_CMRO2.*V_mask_func;
    V_CBF_masked = V_CBF.*V_mask_func;
    
    %remove background
    mask_func_zeros = find(V_mask_func==0);
    V_CMRO2_masked(mask_func_zeros)=[];
    V_CBF_masked(mask_func_zeros)=[];

    %remove zeros inside the regions %CHECK IF CMRO2 AND CBF AND NAN AND
    %NOT ZEROS OR THEY COULD HAVE BOTH TYPES OF UNPHYSICAL VALUES.
    CMRO2_zeros = find(V_CMRO2_masked<=0);    
    %count how many voxels we remove
    n_CMRO2_region_voxels_tot = numel(V_CMRO2_masked);
    n_CMRO2_zeros = numel(CMRO2_zeros);
    percentage_removal_CMRO2 = n_CMRO2_zeros/n_CMRO2_region_voxels_tot;
    percentage_removal_CMRO2_subj(end+1)=percentage_removal_CMRO2;%

%     %remove CMRO2 zeros
%     V_CMRO2_masked(CMRO2_zeros)=[];

    CMRO2_nans=find(isnan(V_CMRO2_masked));
    n_CMRO2_nans=length(CMRO2_nans);
    percentage_nans_CMRO2=n_CMRO2_nans/n_CMRO2_region_voxels_tot;
    percentage_nans_CMRO2_subj(end+1)=percentage_nans_CMRO2;

    %remove zeros inside the regions
    CBF_zeros = find(V_CBF_masked<=0);    
    %count how many voxels we remove
    n_CBF_region_voxels_tot = numel(V_CBF_masked);
    n_CBF_zeros = numel(CBF_zeros);
    percentage_removal_CBF = n_CBF_zeros/n_CBF_region_voxels_tot*100;
    percentage_removal_CBF_subj(end+1)=percentage_removal_CBF;   

%     V_CBF_masked(CBF_zeros)=[];

    CBF_nans=find(isnan(V_CBF_masked));
    n_CBF_nans=length(CBF_nans);
    percentage_nans_CBF=n_CBF_nans/n_CBF_region_voxels_tot;
    percentage_nans_CBF_subj(end+1)=percentage_nans_CBF;

    
    %compute medians
    medians_CMRO2_subj(end+1) = nanmedian(V_CMRO2_masked);
    medians_CBF_subj(end+1) = nanmedian(V_CBF_masked);

    labels_func_subj(end+1) = regions_func(region);

    %calculate PVE regional medians 
    V_pve_0_func_masked = V_pve_0_func_original.*V_mask_func;
    %remove background
    mask_func_zeros = find(V_pve_0_func_masked==0);
    V_pve_0_func_masked(mask_func_zeros)=[];

    V_pve_0_func_median = median(V_pve_0_func_masked(:)); 
    V_pve_0_func_median_subj(end+1) = V_pve_0_func_median;

    V_pve_1_func_masked = V_pve_1_func_original.*V_mask_func;
    %remove background
    mask_func_zeros = find(V_pve_1_func_masked==0);
    V_pve_1_func_masked(mask_func_zeros)=[];

    V_pve_1_func_median = median(V_pve_1_func_masked(:));
    V_pve_1_func_median_subj(end+1) = V_pve_1_func_median;

    V_pve_2_func_masked = V_pve_2_func_original.*V_mask_func;
    %remove background
    mask_func_zeros = find(V_pve_2_func_masked==0);
    V_pve_2_func_masked(mask_func_zeros)=[];

    V_pve_2_func_median = median(V_pve_2_func_masked(:)); 
    V_pve_2_func_median_subj(end+1) = V_pve_2_func_median;

    end
        
    medians_CMRO2_subjs(subj,1:n_regions_func) = medians_CMRO2_subj;
    medians_CBF_subjs(subj,1:n_regions_func) = medians_CBF_subj;
    labels_func_subjs(subj,1:n_regions_func) = labels_func_subj;
    percentage_removal_CMRO2_subjs(subj,1:n_regions_func) = percentage_removal_CMRO2_subj;
    percentage_removal_CBF_subjs(subj,1:n_regions_func) = percentage_removal_CBF_subj;
    percentage_nans_CMRO2_subjs(subj,1:n_regions_func) = percentage_nans_CMRO2_subj;
    percentage_nans_CBF_subjs(subj,1:n_regions_func) = percentage_nans_CBF_subj;

    V_pve_0_func_median_subjs(subj,1:n_regions_func) = V_pve_0_func_median_subj;
    V_pve_1_func_median_subjs(subj,1:n_regions_func) = V_pve_1_func_median_subj;
    V_pve_2_func_median_subjs(subj,1:n_regions_func) = V_pve_2_func_median_subj;

    %%%%DWI space
    regions_dwi = unique(V_atlas_dwi(:));
    %background removal
    regions_dwi(1)=[];
    n_regions_dwi=numel(regions_dwi);

    V_pve_0_dwi_original = V_pve_0_dwi;
    V_pve_1_dwi_original = V_pve_1_dwi;
    V_pve_2_dwi_original = V_pve_2_dwi;

    %define mask
    V_fsoma_to_mask = V_fsoma;
    V_fsoma_to_mask(V_fsoma_to_mask>fsoma_threshold)=1;
    V_fsoma_to_mask(V_fsoma_to_mask<1)=0;

    V_pve_1_dwi(V_pve_1_dwi>pve_1_threshold)=1;
    V_pve_1_dwi(V_pve_1_dwi<1)=0;

    %In case you want to mask by considering not all the WM
    V_pve_2_dwi(V_pve_2_dwi<pve_2_threshold)=NaN;    
    V_pve_2_dwi(V_pve_2_dwi>=pve_2_threshold)=0;
    V_pve_2_dwi(isnan(V_pve_2_dwi))=1;

    % %in case you want all the WM
    % V_pve_2_dwi(V_pve_2_dwi==0)=NaN;
    % V_pve_2_dwi(V_pve_2_dwi>0)=0;
    % V_pve_2_dwi(isnan(V_pve_2_dwi))=1;

    %In case you want to mask by considering not all the CSF
    V_pve_0_dwi(V_pve_0_dwi<pve_0_threshold)=NaN;
    V_pve_0_dwi(V_pve_0_dwi>=pve_0_threshold)=0;
    V_pve_0_dwi(isnan(V_pve_0_dwi))=1;

    % %in case you want all the CSF
    % V_pve_0_dwi(V_pve_0_dwi==0)=NaN;
    % V_pve_0_dwi(V_pve_0_dwi>0)=0;
    % V_pve_0_dwi(isnan(V_pve_0_dwi))=1;


    V_atlas_mask_dwi=V_atlas_dwi;
    %binarize atlas
    for region = 1:n_regions_dwi
    for ii = 1:length(V_atlas_dwi(:))
        if V_atlas_dwi(ii) == region
            V_atlas_mask_dwi(ii) = 1;
        else
            V_atlas_mask_dwi(ii) = 0;
        end
    end



    V_GM = V_pve_1_dwi.*V_pve_2_dwi.*V_pve_0_dwi;
    V_mask_dwi = V_GM.*V_atlas_mask_dwi;
    V_mask_dwi_fs = V_mask_dwi.*V_fsoma_to_mask;

    %mask

    V_rsoma_masked = V_rsoma.*V_mask_dwi_fs;
    V_fsoma_masked = V_fsoma.*V_mask_dwi;
    V_fc_masked = V_fc.*V_mask_dwi;
    V_fsup_masked = V_fsup.*V_mask_dwi;

    V_fneurite_masked = V_fneurite.*V_mask_dwi;  
    V_De_masked = V_De.*V_mask_dwi;
    V_Din_masked = V_Din.*V_mask_dwi;
    V_fextra_masked = V_fextra.*V_mask_dwi;
    
    V_MSE_masked_rsoma = V_MSE.*V_mask_dwi_fs;
    V_MSE_masked = V_MSE.*V_mask_dwi;

    %remove background
    mask_dwi_fs_zeros = find(V_mask_dwi_fs==0);
    mask_dwi_zeros = find(V_mask_dwi==0);

    V_rsoma_masked(mask_dwi_fs_zeros)=[];
    V_fsoma_masked(mask_dwi_zeros)=[];
    V_fc_masked(mask_dwi_zeros)=[];
    V_fsup_masked(mask_dwi_zeros)=[];
    
    V_fneurite_masked(mask_dwi_zeros)=[];
    V_De_masked(mask_dwi_zeros)=[];
    V_Din_masked(mask_dwi_zeros)=[];
    V_fextra_masked(mask_dwi_zeros)=[];
    
    V_MSE_masked_rsoma(mask_dwi_fs_zeros)=[];
    V_MSE_masked(mask_dwi_zeros)=[];

    

    % %remove voxels which have MSE higher than 75th percentile
    idx_high_MSE_micropar=find(V_MSE_masked>prctile(V_MSE_masked,75));
    idx_high_MSE_rsoma=find(V_MSE_masked_rsoma>prctile(V_MSE_masked_rsoma,75));

%FIND INDEX THAT HAVE ERROR HIGHER THAN A CERTAIN DISTRIBUTION
    %PERCENTILE.
    %COUNT THE PERCENTAGE
    percentage_high_MSE_rsoma = numel(idx_high_MSE_rsoma)/numel(V_rsoma_masked);

    percentage_high_MSE_micropar = numel(idx_high_MSE_micropar)/numel(V_fsoma_masked);

    %REMOVE THEM
    V_rsoma_masked(idx_high_MSE_rsoma)=[];
    V_fsoma_masked(idx_high_MSE_micropar)=[];
    V_fneurite_masked(idx_high_MSE_micropar)=[];
    V_De_masked(idx_high_MSE_micropar)=[];
    V_Din_masked(idx_high_MSE_micropar)=[];
    V_fextra_masked(idx_high_MSE_micropar)=[];

    n_rsoma_region_voxels_tot = numel(V_rsoma_masked);%%%%
    rsoma_nans=find(isnan(V_rsoma_masked));%%%%
    n_rsoma_nans=length(rsoma_nans);%%%%
    percentage_nans_rsoma=n_rsoma_nans/n_rsoma_region_voxels_tot;%%%%
    percentage_nans_rsoma_subj(end+1)=percentage_nans_rsoma;%%%%

    %compute medians and save results
    
    medians_rsoma_subj(end+1) = nanmedian(V_rsoma_masked);%
    medians_fsoma_subj(end+1) = nanmedian(V_fsoma_masked);%
    medians_fc_subj(end+1) = nanmedian(V_fc_masked);%
    medians_fsup_subj(end+1) = nanmedian(V_fsup_masked);%

    medians_fneurite_subj(end+1) = nanmedian(V_fneurite_masked);%
    medians_De_subj(end+1) = nanmedian(V_De_masked);%
    medians_Din_subj(end+1) = nanmedian(V_Din_masked);%
    medians_fextra_subj(end+1) = nanmedian(V_fextra_masked);%
    medians_mse_subj(end+1) = nanmedian(V_MSE_masked);%
    medians_mse_rsoma_subj(end+1) = nanmedian(V_MSE_masked_rsoma);%

    labels_dwi_subj(end+1) = regions_dwi(region);%   
    
    percentage_removal_microparameter_subj(end+1) = percentage_high_MSE_micropar;
    percentage_removal_rsoma_subj(end+1) = percentage_high_MSE_rsoma;%%%

    %calculate PVE regional medians %use V_mask_dwi_fs as we are interested
    %in Rsoma
    V_pve_0_dwi_masked = V_pve_0_dwi_original.*V_mask_dwi_fs;
    %remove background
    mask_dwi_zeros = find(V_pve_0_dwi_masked==0);
    V_pve_0_dwi_masked(mask_dwi_zeros)=[];

    V_pve_0_dwi_median = median(V_pve_0_dwi_masked(:)); 
    V_pve_0_dwi_median_subj(end+1) = V_pve_0_dwi_median;

    V_pve_1_dwi_masked = V_pve_1_dwi_original.*V_mask_dwi_fs;
    %remove background
    mask_dwi_zeros = find(V_pve_1_dwi_masked==0);
    V_pve_1_dwi_masked(mask_dwi_zeros)=[];

    V_pve_1_dwi_median = median(V_pve_1_dwi_masked(:));
    V_pve_1_dwi_median_subj(end+1) = V_pve_1_dwi_median;

    V_pve_2_dwi_masked = V_pve_2_dwi_original.*V_mask_dwi_fs;
    %remove background
    mask_dwi_zeros = find(V_pve_2_dwi_masked==0);
    V_pve_2_dwi_masked(mask_dwi_zeros)=[];

    V_pve_2_dwi_median = median(V_pve_2_dwi_masked(:)); 
    V_pve_2_dwi_median_subj(end+1) = V_pve_2_dwi_median;
    end
    
    labels_dwi_subjs(subj,1:n_regions_dwi) = labels_dwi_subj;%
    medians_rsoma_subjs(subj,1:n_regions_dwi) = medians_rsoma_subj;%
    medians_fsoma_subjs(subj,1:n_regions_dwi) = medians_fsoma_subj;%
    medians_fc_subjs(subj,1:n_regions_dwi) = medians_fc_subj;%
    medians_fsup_subjs(subj,1:n_regions_dwi) = medians_fsup_subj;%

    medians_fneurite_subjs(subj,1:n_regions_dwi) = medians_fneurite_subj;%
    medians_Din_subjs(subj,1:n_regions_dwi) = medians_Din_subj;%
    medians_De_subjs(subj,1:n_regions_dwi) = medians_De_subj;%
    medians_fextra_subjs(subj,1:n_regions_dwi) = medians_fextra_subj;%

    medians_mse_rsoma_subjs(subj,1:n_regions_dwi)=medians_mse_rsoma_subj;
    medians_mse_subjs(subj,1:n_regions_dwi)=medians_rsoma_subj;
    percentage_removal_microparameter_subjs(subj,1:n_regions_dwi) = percentage_removal_microparameter_subj;
    percentage_removal_rsoma_subjs(subj,1:n_regions_dwi) = percentage_removal_rsoma_subj;
    percentage_nans_rsoma_subjs(subj,1:n_regions_dwi)=percentage_nans_rsoma_subj;%%%%

    V_pve_0_dwi_median_subjs(subj,1:n_regions_dwi) = V_pve_0_dwi_median_subj;
    V_pve_1_dwi_median_subjs(subj,1:n_regions_dwi) = V_pve_1_dwi_median_subj;
    V_pve_2_dwi_median_subjs(subj,1:n_regions_dwi) = V_pve_2_dwi_median_subj;

    toc
    disp(strcat('Finished subject', num2str(subj),'Starting subject', num2str(subj+1)))

end
timeElapsed = toc(start_time);

%% select variables to examine

energy_parameter = 'CMRO2';
micro_parameter = 'Rsoma';

%% remove uncommon labels between subjects and between spaces
%% select matrices

if strcmp(micro_parameter,'Rsoma')
    medians_dwi_subjs = medians_rsoma_subjs;
elseif strcmp(micro_parameter,'fsoma')
    medians_dwi_subjs = medians_fsoma_subjs;
elseif strcmp(micro_parameter,'fsup')
    medians_dwi_subjs = medians_fsup_subjs;
elseif strcmp(micro_parameter,'fc')
    medians_dwi_subjs = medians_fc_subjs;
end

if strcmp(energy_parameter,'CMRO2')
    medians_func_subjs = medians_CMRO2_subjs;
    percentage_removal_func_subjs = percentage_removal_CMRO2_subjs;
    percentage_nans_func_subjs = percentage_nans_CMRO2_subjs;
elseif strcmp(energy_parameter,'CBF')
    medians_func_subjs = medians_CBF_subjs;
    percentage_removal_func_subjs = percentage_removal_CBF_subjs;
    percentage_nans_func_subjs = percentage_removal_CBF_subjs;
end


%% The following blocks are to select common regions

%% dwi space 

% detect common labels in dwi space across subjs
initial_common = labels_dwi_subjs(1,:);
common = initial_common;
for i = 1:length(labels_dwi_subjs(:,1))
    commonElements = intersect(common, labels_dwi_subjs(i,:));
    common = commonElements;
end

%First, detect indices of common Elements
commonElements(commonElements==0)=[];
lst_idx_matrix_dwi=[];
for row = 1:length(labels_dwi_subjs(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_dwi_subjs(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix_dwi(row,:)=lst_idx_row;
end

%Secondly, select elements corresponding to common indices
%apply this both to the median values and labels
medians_dwi_subjs_final=[];
for row = 1:length(medians_dwi_subjs(:,1))
    medians_dwi_subjs_row=medians_dwi_subjs(row,:);
    lst_idx = lst_idx_matrix_dwi(row,:);
    medians_dwi_subjs_row_final=medians_dwi_subjs_row(lst_idx);
    medians_dwi_subjs_final(row,:)=medians_dwi_subjs_row_final;
end

labels_dwi_subjs_final=[];
for row = 1:length(labels_dwi_subjs(:,1))
    labels_dwi_subjs_row=labels_dwi_subjs(row,:);
    lst_idx = lst_idx_matrix_dwi(row,:);
    labels_dwi_subjs_row_final=labels_dwi_subjs_row(lst_idx);
    labels_dwi_subjs_final(row,:)=labels_dwi_subjs_row_final;
end

V_pve_0_dwi_median_subjs_final=[];
for row = 1:length(V_pve_0_dwi_median_subjs(:,1))
    V_pve_0_dwi_median_subjs_row=V_pve_0_dwi_median_subjs(row,:);
    lst_idx = lst_idx_matrix_dwi(row,:);
    V_pve_0_dwi_median_subjs_row_final=V_pve_0_dwi_median_subjs_row(lst_idx);
    V_pve_0_dwi_median_subjs_final(row,:)=V_pve_0_dwi_median_subjs_row_final;
end

V_pve_1_dwi_median_subjs_final=[];
for row = 1:length(V_pve_1_dwi_median_subjs(:,1))
    V_pve_1_dwi_median_subjs_row=V_pve_1_dwi_median_subjs(row,:);
    lst_idx = lst_idx_matrix_dwi(row,:);
    V_pve_1_dwi_median_subjs_row_final=V_pve_1_dwi_median_subjs_row(lst_idx);
    V_pve_1_dwi_median_subjs_final(row,:)=V_pve_1_dwi_median_subjs_row_final;
end

V_pve_2_dwi_median_subjs_final=[];
for row = 1:length(V_pve_2_dwi_median_subjs(:,1))
    V_pve_2_dwi_median_subjs_row=V_pve_2_dwi_median_subjs(row,:);
    lst_idx = lst_idx_matrix_dwi(row,:);
    V_pve_2_dwi_median_subjs_row_final=V_pve_2_dwi_median_subjs_row(lst_idx);
    V_pve_2_dwi_median_subjs_final(row,:)=V_pve_2_dwi_median_subjs_row_final;
end

if strcmp(micro_parameter,'Rsoma')    
    medians_mse_rsoma_subjs_final=[];
    for row = 1:length(medians_mse_rsoma_subjs(:,1))
        medians_mse_rsoma_subjs_row=medians_mse_rsoma_subjs(row,:);
        lst_idx = lst_idx_matrix_dwi(row,:);
        medians_mse_rsoma_subjs_row_final=medians_mse_rsoma_subjs_row(lst_idx);
        medians_mse_rsoma_subjs_final(row,:)=medians_mse_rsoma_subjs_row_final;
    end
else
    medians_mse_subjs_final=[];
    for row = 1:length(medians_mse_subjs(:,1))
        medians_mse_subjs_row=medians_mse_subjs(row,:);
        lst_idx = lst_idx_matrix_dwi(row,:);
        medians_mse_subjs_row_final=medians_mse_subjs_row(lst_idx);
        medians_mse_subjs_final(row,:)=medians_mse_subjs_row_final;
    end
end
%in case you remove voxels with high MSE values
if strcmp(micro_parameter,'rsoma')    
    percentage_removal_rsoma_subjs_final=[];
    for row = 1:length(percentage_removal_rsoma_subjs(:,1))
        percentage_removal_rsoma_subjs_row=percentage_removal_rsoma_subjs(row,:);
        lst_idx = lst_idx_matrix_dwi(row,:);
        percentage_removal_rsoma_subjs_row_final=percentage_removal_rsoma_subjs_row(lst_idx);
        percentage_removal_rsoma_subjs_final(row,:)=percentage_removal_rsoma_subjs_row_final;
    end
else
    percentage_removal_microparameter_subjs_final=[];
    for row = 1:length(percentage_removal_microparameter_subjs(:,1))
        percentage_removal_microparameter_subjs_row=percentage_removal_microparameter_subjs(row,:);
        lst_idx = lst_idx_matrix_dwi(row,:);
        percentage_removal_microparameter_subjs_row_final=percentage_removal_microparameter_subjs_row(lst_idx);
        percentage_removal_microparameter_subjs_final(row,:)=percentage_removal_microparameter_subjs_row_final;
    end
end

%repeat everything for func space
%% func space
% detect common labels in func space across subjs
initial_common = labels_func_subjs(1,:);
common = initial_common;
for i = 1:length(labels_func_subjs(:,1))
    commonElements = intersect(common, labels_func_subjs(i,:));
    common = commonElements;
end

%First, detect indices of common Elements
commonElements(commonElements==0)=[];
lst_idx_matrix_func=[];
for row = 1:length(labels_func_subjs(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_func_subjs(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix_func(row,:)=lst_idx_row;
end

%Secondly, select elements corresponding to common indices
%apply this both to the median values and labels
medians_func_subjs_final=[];
for row = 1:length(medians_func_subjs(:,1))
    medians_func_subjs_row=medians_func_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    medians_func_subjs_row_final=medians_func_subjs_row(lst_idx);
    medians_func_subjs_final(row,:)=medians_func_subjs_row_final;
end

labels_func_subjs_final=[];
for row = 1:length(labels_func_subjs(:,1))
    labels_func_subjs_row=labels_func_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    labels_func_subjs_row_final=labels_func_subjs_row(lst_idx);
    labels_func_subjs_final(row,:)=labels_func_subjs_row_final;
end

V_pve_0_func_median_subjs_final=[];
for row = 1:length(V_pve_0_func_median_subjs(:,1))
    V_pve_0_func_median_subjs_row=V_pve_0_func_median_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    V_pve_0_func_median_subjs_row_final=V_pve_0_func_median_subjs_row(lst_idx);
    V_pve_0_func_median_subjs_final(row,:)=V_pve_0_func_median_subjs_row_final;
end

V_pve_1_func_median_subjs_final=[];
for row = 1:length(V_pve_1_func_median_subjs(:,1))
    V_pve_1_func_median_subjs_row=V_pve_1_func_median_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    V_pve_1_func_median_subjs_row_final=V_pve_1_func_median_subjs_row(lst_idx);
    V_pve_1_func_median_subjs_final(row,:)=V_pve_1_func_median_subjs_row_final;
end

V_pve_2_func_median_subjs_final=[];
for row = 1:length(V_pve_2_func_median_subjs(:,1))
    V_pve_2_func_median_subjs_row=V_pve_2_func_median_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    V_pve_2_func_median_subjs_row_final=V_pve_2_func_median_subjs_row(lst_idx);
    V_pve_2_func_median_subjs_final(row,:)=V_pve_2_func_median_subjs_row_final;
end

%percentage of zeros in the different regions of different subjects
percentage_removal_func_subjs_final=[];
for row = 1:length(percentage_removal_func_subjs(:,1))
    percentage_removal_func_subjs_row=percentage_removal_func_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    percentage_removal_func_subjs_row_final=percentage_removal_func_subjs_row(lst_idx);
    percentage_removal_func_subjs_final(row,:)=percentage_removal_func_subjs_row_final;
end

%percentage of nans in the different regions of different subjects
percentage_nans_func_subjs_final=[];
for row = 1:length(percentage_nans_func_subjs(:,1))
    percentage_nans_func_subjs_row=percentage_nans_func_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    percentage_nans_func_subjs_row_final=percentage_nans_func_subjs_row(lst_idx);
    percentage_nans_func_subjs_final(row,:)=percentage_nans_func_subjs_row_final;
end

%%
%then for each subj (row of the two matrices),
%find common elements (regions) between a (dwi) and b (func) space,
%find the idx 
%keep only those common both in DWI and func space
%create a new dwi and func matrix which will have same size.

%%
%select common elements between the two spaces

commonElements_lst=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    commonElements = intersect(labels_dwi_subjs_final(row,:),labels_func_subjs_final(row,:));
    if length(commonElements)>1
        commonElements_lst(1,:)=commonElements;
    else
        commonElements_lst(end+1)=commonElements;
    end
end

commonElements_lst=unique(commonElements_lst);
commonElements(commonElements==0)=[];

%% dwi space
%First, detect indices of common Elements

lst_idx_matrix_dwi_spaces=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_dwi_subjs_final(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix_dwi_spaces(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_dwi_subjs_final_spaces=[];
for row = 1:length(medians_dwi_subjs_final(:,1))
    medians_dwi_subjs_row=medians_dwi_subjs_final(row,:);
    lst_idx = lst_idx_matrix_dwi_spaces(row,:);
    medians_dwi_subjs_row_final=medians_dwi_subjs_row(lst_idx);
    medians_dwi_subjs_final_spaces(row,:)=medians_dwi_subjs_row_final;
end

labels_dwi_subjs_final_spaces=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    labels_dwi_subjs_row=labels_dwi_subjs_final(row,:);
    lst_idx = lst_idx_matrix_dwi_spaces(row,:);
    labels_dwi_subjs_row_final=labels_dwi_subjs_row(lst_idx);
    labels_dwi_subjs_final_spaces(row,:)=labels_dwi_subjs_row_final;
end

V_pve_0_dwi_median_subjs_final_spaces=[];
for row = 1:length(V_pve_0_dwi_median_subjs_final(:,1))
    V_pve_0_dwi_median_subjs_row=V_pve_0_dwi_median_subjs_final(row,:);
    lst_idx = lst_idx_matrix_dwi_spaces(row,:);
    V_pve_0_dwi_median_subjs_row_final=V_pve_0_dwi_median_subjs_row(lst_idx);
    V_pve_0_dwi_median_subjs_final_spaces(row,:)=V_pve_0_dwi_median_subjs_row_final;
end

V_pve_1_dwi_median_subjs_final_spaces=[];
for row = 1:length(V_pve_1_dwi_median_subjs_final(:,1))
    V_pve_1_dwi_median_subjs_row=V_pve_1_dwi_median_subjs_final(row,:);
    lst_idx = lst_idx_matrix_dwi_spaces(row,:);
    V_pve_1_dwi_median_subjs_row_final=V_pve_1_dwi_median_subjs_row(lst_idx);
    V_pve_1_dwi_median_subjs_final_spaces(row,:)=V_pve_1_dwi_median_subjs_row_final;
end

V_pve_2_dwi_median_subjs_final_spaces=[];
for row = 1:length(V_pve_2_dwi_median_subjs_final(:,1))
    V_pve_2_dwi_median_subjs_row=V_pve_2_dwi_median_subjs_final(row,:);
    lst_idx = lst_idx_matrix_dwi_spaces(row,:);
    V_pve_2_dwi_median_subjs_row_final=V_pve_2_dwi_median_subjs_row(lst_idx);
    V_pve_2_dwi_median_subjs_final_spaces(row,:)=V_pve_2_dwi_median_subjs_row_final;
end

if strcmp(micro_parameter,'Rsoma')    
    medians_mse_rsoma_subjs_final_spaces=[];
    for row = 1:length(medians_mse_rsoma_subjs_final(:,1))
        medians_mse_rsoma_subjs_row=medians_mse_rsoma_subjs_final(row,:);
        lst_idx = lst_idx_matrix_dwi_spaces(row,:);
        medians_mse_rsoma_subjs_row_final=medians_mse_rsoma_subjs_row(lst_idx);
        medians_mse_rsoma_subjs_final_spaces(row,:)=medians_mse_rsoma_subjs_row_final;
    end
else
    medians_mse_subjs_final_spaces=[];
    for row = 1:length(medians_mse_subjs_final(:,1))
        medians_mse_subjs_row=medians_mse_subjs_final(row,:);
        lst_idx = lst_idx_matrix_dwi_spaces(row,:);
        medians_mse_subjs_row_final=medians_mse_subjs_row(lst_idx);
        medians_mse_subjs_final_spaces(row,:)=medians_mse_subjs_row_final;
    end
end
%In case you remove high MSE voxels
if strcmp(micro_parameter,'rsoma')
    percentage_removal_rsoma_subjs_final_spaces=[];
    for row = 1:length(percentage_removal_rsoma_subjs_final(:,1))
        percentage_removal_rsoma_subjs_row=percentage_removal_rsoma_subjs_final(row,:);
        lst_idx = lst_idx_matrix_dwi_spaces(row,:);
        percentage_removal_rsoma_subjs_row_final=percentage_removal_rsoma_subjs_row(lst_idx);
        percentage_removal_rsoma_subjs_final_spaces(row,:)=percentage_removal_rsoma_subjs_row_final;
    end
else
    percentage_removal_microparameter_subjs_final_spaces=[];
    for row = 1:length(percentage_removal_microparameter_subjs_final(:,1))
        percentage_removal_microparameter_subjs_row=percentage_removal_microparameter_subjs_final(row,:);
        lst_idx = lst_idx_matrix_dwi_spaces(row,:);
        percentage_removal_microparameter_subjs_row_final=percentage_removal_microparameter_subjs_row(lst_idx);
        percentage_removal_microparameter_subjs_final_spaces(row,:)=percentage_removal_microparameter_subjs_row_final;
    end
end
%% func space

%First, detect indices of common Elements

lst_idx_matrix_func_spaces=[];
for row = 1:length(labels_func_subjs_final(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_func_subjs_final(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix_func_spaces(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_func_subjs_final_spaces=[];
for row = 1:length(medians_func_subjs_final(:,1))
    medians_func_subjs_row=medians_func_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    medians_func_subjs_row_final=medians_func_subjs_row(lst_idx);
    medians_func_subjs_final_spaces(row,:)=medians_func_subjs_row_final;
end


labels_func_subjs_final_spaces=[];
for row = 1:length(labels_func_subjs_final(:,1))
    labels_func_subjs_row=labels_func_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    labels_func_subjs_row_final=labels_func_subjs_row(lst_idx);
    labels_func_subjs_final_spaces(row,:)=labels_func_subjs_row_final;
end

V_pve_0_func_median_subjs_final_spaces=[];
for row = 1:length(V_pve_0_func_median_subjs_final(:,1))
    V_pve_0_func_median_subjs_row=V_pve_0_func_median_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    V_pve_0_func_median_subjs_row_final=V_pve_0_func_median_subjs_row(lst_idx);
    V_pve_0_func_median_subjs_final_spaces(row,:)=V_pve_0_func_median_subjs_row_final;
end

V_pve_1_func_median_subjs_final_spaces=[];
for row = 1:length(V_pve_1_func_median_subjs_final(:,1))
    V_pve_1_func_median_subjs_row=V_pve_1_func_median_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    V_pve_1_func_median_subjs_row_final=V_pve_1_func_median_subjs_row(lst_idx);
    V_pve_1_func_median_subjs_final_spaces(row,:)=V_pve_1_func_median_subjs_row_final;
end

V_pve_2_func_median_subjs_final_spaces=[];
for row = 1:length(V_pve_2_func_median_subjs_final(:,1))
    V_pve_2_func_median_subjs_row=V_pve_2_func_median_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    V_pve_2_func_median_subjs_row_final=V_pve_2_func_median_subjs_row(lst_idx);
    V_pve_2_func_median_subjs_final_spaces(row,:)=V_pve_2_func_median_subjs_row_final;
end

%percentage of zeros in the different regions of different subjects
percentage_removal_func_subjs_final_spaces=[];
for row = 1:length(percentage_removal_func_subjs_final(:,1))
    percentage_removal_func_subjs_row=percentage_removal_func_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    percentage_removal_func_subjs_row_final=percentage_removal_func_subjs_row(lst_idx);
    percentage_removal_func_subjs_final_spaces(row,:)=percentage_removal_func_subjs_row_final;
end

%percentage of nans in the different regions of different subjects
percentage_nans_func_subjs_final_spaces=[];
for row = 1:length(percentage_nans_func_subjs_final(:,1))
    percentage_nans_func_subjs_row=percentage_nans_func_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    percentage_nans_func_subjs_row_final=percentage_nans_func_subjs_row(lst_idx);
    percentage_nans_func_subjs_final_spaces(row,:)=percentage_nans_func_subjs_row_final;
end

%% Remove regions based on the number of NaNs

% %Check number of nan values vs subjects for each region
% for region = 1:length(labels_func_subjs_final_spaces(1,:))
%     figure, bar(percentage_nans_func_subjs(:,region))
%     title(strcat('Region',num2str(region)))
%     ylabel(strcat('percentage CMRO_{2} nan values  removed'));
%     xlabel('subjects')
% end

% count how many subjects have number of nans > 50% of total number of
% voxels for each region
percentage=0.5;
n_subjects_tot=[];
for region = 1:length(labels_func_subjs_final_spaces(1,:))
    n_subjects = numel(find(percentage_nans_func_subjs_final_spaces(:,region)>percentage))/n_subjs;
    n_subjects_tot(end+1)=n_subjects;
end

% figure, bar(n_subjects_tot);
% xlabel('Region');
% ylabel('Number of subjects ratio')
% title(strcat('Number of NaNs >',num2str(percentage*100),'%'));
% grid on

nans_threshold=0.5;
regions_idx=find(n_subjects_tot>nans_threshold);

medians_func = medians_func_subjs_final_spaces;
medians_dwi = medians_dwi_subjs_final_spaces;
medians_mse = medians_mse_rsoma_subjs_final_spaces;

medians_pve_0_func = V_pve_0_func_median_subjs_final_spaces;
medians_pve_1_func = V_pve_1_func_median_subjs_final_spaces;
medians_pve_2_func = V_pve_2_func_median_subjs_final_spaces;

medians_pve_0_dwi = V_pve_0_dwi_median_subjs_final_spaces;
medians_pve_1_dwi = V_pve_1_dwi_median_subjs_final_spaces;
medians_pve_2_dwi = V_pve_2_dwi_median_subjs_final_spaces;

medians_func(:,regions_idx)=[];
medians_dwi(:,regions_idx)=[];
medians_mse(:,regions_idx)=[];

medians_pve_0_func(:,regions_idx)=[];
medians_pve_1_func(:,regions_idx)=[];
medians_pve_2_func(:,regions_idx)=[];

medians_pve_0_dwi(:,regions_idx)=[];
medians_pve_1_dwi(:,regions_idx)=[];
medians_pve_2_dwi(:,regions_idx)=[];



labels_final=labels_func_subjs_final_spaces(1,:);
labels_final(regions_idx)=[];
disp('Regions with high number of NaNs successfully removed');

%% Remove regions based on the number of high idx MSE


percentage=0.5;
n_subjects_tot=[];
if strcmp(micro_parameter,'rsoma')
    for region = 1:length(labels_func_subjs_final_spaces(1,:))
        n_subjects = numel(find(percentage_removal_rsoma_subjs(:,region)>percentage))/n_subjs;
        n_subjects_tot(end+1)=n_subjects;
    end
else
    for region = 1:length(labels_func_subjs_final_spaces(1,:))
        n_subjects = numel(find(percentage_removal_microparameter_subjs(:,region)>percentage))/n_subjs;
        n_subjects_tot(end+1)=n_subjects;
    end
end
% figure, bar(n_subjects_tot);
% xlabel('Region');
% ylabel('Number of subjects ratio')
% title(strcat('Number of NaNs >',num2str(percentage*100),'%'));
% grid on

high_mse_threshold=0.5;
regions_idx=find(n_subjects_tot>high_mse_threshold);

medians_func = medians_func_subjs_final_spaces;
medians_dwi = medians_dwi_subjs_final_spaces;
if strcmp(micro_parameter,'Rsoma')
    medians_mse = medians_mse_rsoma_subjs_final_spaces;
else
    medians_mse = medians_mse_subjs_final_spaces;
end
medians_func(:,regions_idx)=[];
medians_dwi(:,regions_idx)=[];
medians_mse(:,regions_idx)=[];
disp('Regions with high number of high MSE values successfully removed');

%% across regions (medians across subjs)

%for median analysis
% medians_func = medians_func_subjs_final_spaces;
% medians_dwi = medians_dwi_subjs_final_spaces;

if strcmp(energy_parameter,'CBF')   
    dependent_parameter='CBF (ml/100g/min)';
elseif strcmp(energy_parameter,'CMRO2')
    dependent_parameter='CMRO_2 (\mu mol/100g/min)';
end

medians_energy_vec = nanmedian(medians_func,1);
SE_energy = nanstd(medians_func,0,1)/sqrt(n_subjs);
medians_micro_parameter_vec = nanmedian(medians_dwi,1);
SE_micro_parameter = nanstd(medians_dwi,0,1)/sqrt(n_subjs);
medians_mse = nanmedian(medians_mse,1);

cv_micro_parameter = SE_micro_parameter./medians_micro_parameter_vec;
medians_pve_0_dwi_vec = nanmedian(medians_pve_0_dwi,1);
SE_pve_0_dwi = nanstd(medians_pve_0_dwi_vec,0,1)/sqrt(n_subjs);
medians_pve_1_dwi_vec = nanmedian(medians_pve_1_dwi,1);
SE_pve_1_dwi = nanstd(medians_pve_1_dwi_vec,0,1)/sqrt(n_subjs);
medians_pve_2_dwi_vec = nanmedian(medians_pve_2_dwi,1);
SE_pve_2_dwi = nanstd(medians_pve_2_dwi_vec,0,1)/sqrt(n_subjs);

medians_pve_0_func_vec = nanmedian(medians_pve_0_func,1);
SE_pve_0_func = nanstd(medians_pve_0_func_vec,0,1)/sqrt(n_subjs);
medians_pve_1_func_vec = nanmedian(medians_pve_1_func,1);
SE_pve_1_func = nanstd(medians_pve_1_func_vec,0,1)/sqrt(n_subjs);
medians_pve_2_func_vec = nanmedian(medians_pve_2_func,1);
SE_pve_2_func = nanstd(medians_pve_2_func_vec,0,1)/sqrt(n_subjs);

%%
%compare std vs mse and cv vs mse

figure, 
bar(SE_micro_parameter);
legend('Standard Error')

figure,
bar(medians_mse);
ylim([0,1.5*10^(-3)])
legend('MSE')

figure,
bar(cv_micro_parameter);
legend('CV')



idx_mse = find(medians_mse>prctile(medians_mse,75));
% isequal(labels_func_subjs_final_spaces,labels_dwi_subjs_final_spaces)
labels_high_mse=labels_final(idx_mse);

%CV
idx_cv = find(cv_micro_parameter>prctile(cv_micro_parameter,75));
labels_high_cv = labels_final(idx_cv);
%size(labels_high_mse)==size(labels_high_cv)

%find labels that have both high cv and high MSE
equal_labels=intersect(labels_high_cv,labels_high_mse);
percentage_equal_labels=numel(equal_labels)/numel(labels_high_mse);% 58%
percentage_equal_labels

%SE
idx_SE = find(SE_micro_parameter>prctile(SE_micro_parameter,75));
labels_high_SE = labels_final(idx_SE);
%size(labels_high_mse)==size(labels_high_SE)

%find labels that have both high SE and MSE
equal_labels=intersect(labels_high_SE,labels_high_mse);
percentage_equal_labels=numel(equal_labels)/numel(labels_high_mse);% 55%
percentage_equal_labels

equal_labels_idx_SE_MSE=[];
for i = 1:numel(equal_labels)
    idx=find(labels_final==equal_labels(i));
    equal_labels_idx_SE_MSE(end+1)=idx;
end

% V_mse_maps_1=V_mse_maps{1};
% median(V_mse_maps_1(:))
% figure, imagesc(V_mse_maps_1(:,:,40))

%% remove regions which have high SE and MSE values
medians_energy_vec(equal_labels_idx_SE_MSE)=[];
medians_micro_parameter_vec(equal_labels_idx_SE_MSE)=[];
SE_energy(equal_labels_idx_SE_MSE)=[];
SE_micro_parameter(equal_labels_idx_SE_MSE)=[];

%% remove regions which have high SE values
medians_energy_vec(idx_SE)=[];
medians_micro_parameter_vec(idx_SE)=[];
SE_energy(idx_SE)=[];
SE_micro_parameter(idx_SE)=[];

%% remove regions which have high CV values
medians_energy_vec(idx_cv)=[];
medians_micro_parameter_vec(idx_cv)=[];
SE_energy(idx_cv)=[];
SE_micro_parameter(idx_cv)=[];

%%

if strcmp(micro_parameter,'Rsoma')
    unit_of_measure='(\mum)';
elseif strcmp(micro_parameter,'fsoma')
    unit_of_measure='';
elseif strcmp(micro_parameter,'fsup')
    unit_of_measure='(m^{-1})';
elseif strcmp(micro_parameter,'fc')
    unit_of_measure='(m^{-3})';
elseif strcmp(micro_parameter,'fneurite')
    unit_of_measure='';
end

%% Run if you to remove outliers identified by eyes
rsoma_limit=11.1;
idx_outlier=find(medians_micro_parameter_vec<rsoma_limit);

medians_energy_vec(idx_outlier)=[];
medians_micro_parameter_vec(idx_outlier)=[];
SE_energy(idx_outlier)=[];
SE_micro_parameter(idx_outlier)=[];

medians_pve_0_dwi_vec(idx_outlier)=[];
medians_pve_1_dwi_vec(idx_outlier)=[];
medians_pve_2_dwi_vec(idx_outlier)=[];

medians_pve_0_func_vec(idx_outlier)=[];
medians_pve_1_func_vec(idx_outlier)=[];
medians_pve_2_func_vec(idx_outlier)=[];

SE_pve_0_func(idx_outlier)=[];
SE_pve_1_func(idx_outlier)=[];
SE_pve_2_func(idx_outlier)=[];

SE_pve_0_dwi(idx_outlier)=[];
SE_pve_1_dwi(idx_outlier)=[];
SE_pve_2_dwi(idx_outlier)=[];
% %find which are low rsoma regions
% disp(strcat('Lower rsoma labels are:',num2str(labels_final(idx_outlier))))
% numel(idx_outlier)
% find(labels_final==46)
%%
[r,p]=corrcoef(medians_energy_vec,medians_micro_parameter_vec,'rows','complete');%0.37%p0.0029
%try other kind of correlation
corr_coef_str=num2str(round(r(2),2));

% remove nan microstructural regions
idx = find(isnan(medians_micro_parameter_vec))
medians_energy_vec(idx)=[];
medians_micro_parameter_vec(idx)=[];
SE_energy(idx)=[];
SE_micro_parameter(idx)=[];

medians_pve_0_dwi_vec(idx)=[];
medians_pve_1_dwi_vec(idx)=[];
medians_pve_2_dwi_vec(idx)=[];

medians_pve_0_func_vec(idx)=[];
medians_pve_1_func_vec(idx)=[];
medians_pve_2_func_vec(idx)=[];

SE_pve_0_func(idx)=[];
SE_pve_1_func(idx)=[];
SE_pve_2_func(idx)=[];

SE_pve_0_dwi(idx)=[];
SE_pve_1_dwi(idx)=[];
SE_pve_2_dwi(idx)=[];

% 
idx = find(isnan(medians_energy_vec))
medians_energy_vec(idx)=[];
medians_micro_parameter_vec(idx)=[];
SE_energy(idx)=[];
SE_micro_parameter(idx)=[];

medians_pve_0_dwi_vec(idx)=[];
medians_pve_1_dwi_vec(idx)=[];
medians_pve_2_dwi_vec(idx)=[];

medians_pve_0_func_vec(idx)=[];
medians_pve_1_func_vec(idx)=[];
medians_pve_2_func_vec(idx)=[];

SE_pve_0_func(idx)=[];
SE_pve_1_func(idx)=[];
SE_pve_2_func(idx)=[];

SE_pve_0_dwi(idx)=[];
SE_pve_1_dwi(idx)=[];
SE_pve_2_dwi(idx)=[];

y=medians_energy_vec;
x=medians_micro_parameter_vec;
% Fit a quadratic equation
pol = polyfit(x, y, 2);
% Evaluate the fitted polynomial
y_fit = polyval(pol, x);
% Calculate the R-squared value
Rsquared = 1 - sum((y - y_fit).^2) / sum((y - nanmean(y)).^2);
Rsquared

[x_sorted,index] = sortrows(x');
y_fit=y_fit';
y_fit_sorted = y_fit(index);

%%
%medians_micro_parameter_vec(equal_labels_idx)=
%%
figure, 
s = errorbar(medians_micro_parameter_vec, medians_energy_vec, SE_energy, SE_energy, SE_micro_parameter, SE_micro_parameter,'o');
% hold on
% plot(x_sorted,y_fit_sorted,'--','LineWidth',3,'Color',"#000000");
% diff=setdiff(medians_micro_parameter_vec_ALLpves,medians_micro_parameter_vec_GMpves);
% idx_diff=find(medians_micro_parameter_vec==diff(1));
% hold on
% h = errorbar(medians_micro_parameter_vec(idx_diff), medians_energy_vec(idx_diff), SE_energy(idx_diff), SE_energy(idx_diff), SE_micro_parameter(idx_diff), SE_micro_parameter(idx_diff),'o');
xlabel(strcat(micro_parameter,unit_of_measure),'FontSize',15,'FontWeight','bold');
ylabel(dependent_parameter,'FontSize',15,'FontWeight','bold');
s.LineWidth = 0.6;
% h.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
% h.MarkerEdgeColor = 'r';
s.MarkerFaceColor = [0 0.5 0.5];
% h.MarkerFaceColor = [0 0.5 0.5];
% if p(2)<0.05 && p(2)>0.01    
%     txt = {strcat('r = ',corr_coef_str,'*')};
% elseif p(2)<0.01 && p(2)>0.001   
%     txt = {strcat('r = ',corr_coef_str,'**')};
% elseif p(2)<0.001
%     txt = {strcat('r = ',corr_coef_str,'***')};
% elseif p(2)>0.05
%         txt = {strcat('r = ',corr_coef_str,'')};
% end
% text(60,12,txt,'FontWeight', 'Bold');
% annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold','FontSize',15);
% if strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'Rsoma')
%     text(13.8,60,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fsoma')
%     text(0.39,60,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'Rsoma')
%     text(13.8,140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fsoma')
%     text(0.39,140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fc')
%     text(3*10^(13),140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fsup')
%     text(8*10^4,140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fneurite')
%     text(0.3,140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fsup')
%     text(8*10^4,60,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fc')
%     text(3*10^13,60,txt, 'FontWeight', 'bold','FontSize',12);
% 
% end
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'box','off')
x0=400;
y0=400;
width=550;
height=450;
set(gcf,'position',[x0,y0,width,height]);
grid on

%% check relationship between parameters and PVE medians

%select parameters
pve_tissue = 'GM';
parameter = 'dwi';

%select data
if strcmp(parameter,'dwi')
    medians_parameter = medians_micro_parameter_vec;
    SE_parameter = SE_micro_parameter;
    if strcmp(pve_tissue,'CSF')
        medians_pve = medians_pve_0_dwi_vec;
        SE_pve = SE_pve_0_dwi;
    elseif strcmp(pve_tissue,'WM')
        medians_pve = medians_pve_2_dwi_vec;
        SE_pve = SE_pve_2_dwi;
    elseif strcmp(pve_tissue,'GM')
        medians_pve = medians_pve_1_dwi_vec;
        SE_pve = SE_pve_1_dwi;
    end
elseif strcmp(parameter,'func')
    medians_parameter = medians_energy_vec;
    SE_parameter = SE_energy;
    if strcmp(pve_tissue,'CSF')
        medians_pve = medians_pve_0_func_vec;
        SE_pve = SE_pve_0_func;
    elseif strcmp(pve_tissue,'WM')
        medians_pve = medians_pve_2_func_vec;
        SE_pve = SE_pve_2_func;
    elseif strcmp(pve_tissue,'GM')
        medians_pve = medians_pve_1_func_vec;
        SE_pve = SE_pve_1_func;
        
    end
end



[r,p]=corrcoef(medians_parameter,medians_pve,'rows','complete');
%try other kind of correlation
corr_coef_str=num2str(round(r(2),2));

%plot (manually select y label)
figure, 
s = errorbar(medians_pve, medians_parameter, SE_parameter, SE_parameter, SE_pve, SE_pve,'o');
% hold on
% plot(x_sorted,y_fit_sorted,'--','LineWidth',3,'Color',"#000000");
xlabel(strcat(pve_tissue,'pve'),'FontSize',15,'FontWeight','bold');
ylabel(strcat(micro_parameter,unit_of_measure),'FontSize',15,'FontWeight','bold');
s.LineWidth = 0.6;
% h.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
% h.MarkerEdgeColor = 'r';
s.MarkerFaceColor = [0 0.5 0.5];
% if p(2)<0.05 && p(2)>0.01    
%     txt = {strcat('r = ',corr_coef_str,'*')};
% elseif p(2)<0.01 && p(2)>0.001   
%     txt = {strcat('r = ',corr_coef_str,'**')};
% elseif p(2)<0.001
%     txt = {strcat('r = ',corr_coef_str,'***')};
% elseif p(2)>0.05
%         txt = {strcat('r = ',corr_coef_str,'')};
% end
% title('Regional medians')
% text(60,12,txt,'FontWeight', 'Bold');
% annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold','FontSize',15);
% if strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'Rsoma')
%     text(13.8,60,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fsoma')
%     text(0.39,60,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'Rsoma')
%     text(13.8,140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fsoma')
%     text(0.39,140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fc')
%     text(3*10^(13),140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fsup')
%     text(8*10^4,140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fneurite')
%     text(0.3,140,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fsup')
%     text(8*10^4,60,txt, 'FontWeight', 'bold','FontSize',12);
% elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fc')
%     text(3*10^13,60,txt, 'FontWeight', 'bold','FontSize',12);
% 
% end
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'box','off')
%where the plot will be displayed
x0=400;
y0=400;
width=550;
height=450;
set(gcf,'position',[x0,y0,width,height]);
grid on
%% across regions (for each subject)

corr_for_each_subj = [];
for i = 1:n_subjs
    [r,p] = corrcoef(medians_dwi(i,:),medians_func(i,:),'rows','complete');
    %'rows','complete', Omit any rows of the input containing NaN values
    corr_for_each_subj(end+1)=r(2);
end

mean_with_subjs_corr=nanmean(corr_for_each_subj);%with fsup we could have NaNs
[h,p,ci,stats]=ttest(atanh(corr_for_each_subj));


mean_corr=round(mean_with_subjs_corr,2);
mean_corr=num2str(mean_corr);

pvalue=num2str(round(p,2));


figure, 
s=histogram(corr_for_each_subj,'FaceAlpha',1,'BinWidth',0.07);
s.FaceColor="b";
xlabel('correlation coefficient, r','FontWeight','bold','FontSize',15);
ylabel('Counts (# subjects)','FontWeight','bold','FontSize',15);
ylim([0,16.5]);
xline(0,'--','LineWidth',3);
if p<0.05 && p>0.01    
    txt = {strcat('\mu_r = ',mean_corr,'*')};
elseif p<0.01 && p>0.001   
    txt = {strcat('\mu_r = ',mean_corr,'**')};
else
    txt = {strcat('\mu_r = ',mean_corr,'***')};
end

title(strcat(energy_parameter, 'vs',micro_parameter));
text(-0.5,7,txt, 'FontWeight', 'bold','FontSize',15);
grid on






%% across subjs (GM median)