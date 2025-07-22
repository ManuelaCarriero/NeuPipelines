%% load data
% load subjs idx
run='run-01';
subjects = importdata(strcat('/media/nas_rete/Vitality/code/subjs_DWI.txt'));
n_subjs=length(subjects);

%load pve on subj space
V_pves_dwi={};
V_pves_func={};

for i = 1:n_subjs
    
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

%load atlas on subj space
V_atlases_dwi={};
V_atlases_func={};

for i = 1:n_subjs
    
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

%load parametric maps
V_CMRO2_maps={};
V_CBF_maps={};

V_rsoma_maps={};
V_fsoma_maps={};
V_fneurite_maps={};
V_De_maps={};
V_Din_maps={};
V_fextra_maps={};

for i = 1:n_subjs
    
    subj = subjects{i};
    
    img_path_CMRO2=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/perf/outcome/',subj,'_task-bh_',run,'_acq-dexi_volreg_asl_topup_CMRO2_map.nii.gz');
    img_path_CBF=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/perf/outcome/',subj,'_task-bh_',run,'_acq-dexi_volreg_asl_topup_CBF_map.nii.gz');

    img_path_rsoma=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_Rsoma.nii.gz');
    img_path_fsoma=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_fsoma.nii.gz');
    img_path_fneurite=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_fneurite.nii.gz');
    img_path_De=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_De.nii.gz');
    img_path_Din=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_Din.nii.gz');
    img_path_fextra=strcat('/media/nas_rete/Vitality/derivatives/',subj,'/dwi/SANDI_Output/',subj,'_',run,'_SANDI-fit_fextra.nii.gz');

    V_vol_CMRO2 = spm_vol(img_path_CMRO2);
    V_CMRO2=spm_read_vols(V_vol_CMRO2);
    V_CMRO2_maps{end+1}=V_CMRO2;

    V_vol_CBF = spm_vol(img_path_CBF);
    V_CBF=spm_read_vols(V_vol_CBF);
    V_CBF_maps{end+1}=V_CBF;

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

end

%%  Count total number of regions
img_path_atlas='/storage/shared/Atlas/AAL3v1_2mm_resampled.nii.gz';
Vhdr = spm_vol(img_path_atlas);
V_atlas_tot = spm_read_vols(Vhdr);

regions = unique(V_atlas_tot(:));
%background removal
regions(1)=[];
n_regions=numel(regions);

%% select binary masks thresholds
pve_threshold=0.5;
fsoma_threshold=0.15;

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
medians_fneurite_subjs = zeros(n_subjs,n_regions);
medians_fextra_subjs = zeros(n_subjs,n_regions);
medians_Din_subjs = zeros(n_subjs,n_regions);
medians_De_subjs = zeros(n_subjs,n_regions);

percentage_removal_CMRO2_subjs = zeros(n_subjs,n_regions);
percentage_removal_CBF_subjs = zeros(n_subjs,n_regions); 

start_time=tic;
for subj = 1:n_subjs
    tic
    %load atlases
    V_atlas_func = V_atlases_func{subj};
    V_atlas_dwi = V_atlases_dwi{subj};

    %load pve maps
    V_pve_func = V_pves_func{subj};
    V_pve_dwi = V_pves_dwi{subj};

    %load parametric maps
    V_CMRO2 = V_CMRO2_maps{subj};
    V_CBF = V_CBF_maps{subj};

    V_rsoma = V_rsoma_maps{subj};
    V_fsoma = V_fsoma_maps{subj};
    V_fneurite = V_fneurite_maps{subj};
    V_Din = V_Din_maps{subj};
    V_De = V_De_maps{subj};
    V_fextra = V_fextra_maps{subj};
    %V_MSE = V_MSE_maps{subj};

    %where to save medians
    medians_CMRO2_subj = [];
    medians_CBF_subj = [];
    medians_rsoma_subj = [];
    medians_fsoma_subj = [];
    medians_fneurite_subj = [];
    medians_Din_subj = [];
    medians_De_subj = [];
    medians_fextra_subj = [];
    labels_func_subj = [];
    labels_dwi_subj = [];
    percentage_removal_CBF_subj = [];
    percentage_removal_CMRO2_subj = [];

    %FUNC SPACE
    %list of regions (it can vary among subjs and among spaces)
    regions_func = unique(V_atlas_func(:));
    %background removal
    regions_func(1)=[];
    n_regions_func=numel(regions_func);
    
    V_atlas_mask_func=V_atlas_func;
    %binarize atlas
    for region = 1:n_regions_func
    for ii = 1:length(V_atlas_func(:))
        if V_atlas_func(ii) == regions_func(region)
            V_atlas_mask_func(ii) = 1;
        else
            V_atlas_mask_func(ii) = 0;
        end
    end
    %define mask
    V_pve_func(V_pve_func>pve_threshold)=1;
    V_pve_func(V_pve_func<1)=0;
    V_mask_func = V_pve_func.*V_atlas_mask_func;

    %mask

    V_CMRO2_masked = V_CMRO2.*V_mask_func;
    V_CBF_masked = V_CBF.*V_mask_func;
    
    %remove background
    mask_func_zeros = find(V_mask_func==0);
    V_CMRO2_masked(mask_func_zeros)=[];
    V_CBF_masked(mask_func_zeros)=[];

    %remove zeros inside the regions %CHECK IF CMRO2 AND CBF AND NAN AND
    %NOT ZEROS OR THEY COULD HAVE BOTH TYPES OF UNPHYSICAL VALUES.
    CMRO2_zeros = find(V_CMRO2_masked<=0);    %CMRO2_nans=find(isnan(V_CMRO2_masked));
    %count how many voxels we remove
    n_CMRO2_region_voxels_tot = numel(V_CMRO2_masked);
    n_CMRO2_zeros = numel(CMRO2_zeros);
    percentage_removal_CMRO2 = n_CMRO2_zeros/n_CMRO2_region_voxels_tot*100;
    percentage_removal_CMRO2_subj(end+1)=percentage_removal_CMRO2;%

    %remove CMRO2 zeros
    V_CMRO2_masked(CMRO2_zeros)=[];

    %remove zeros inside the regions
    CBF_zeros = find(V_CBF_masked<=0);    
    %count how many voxels we remove
    n_CBF_region_voxels_tot = numel(V_CBF_masked);
    n_CBF_zeros = numel(CBF_zeros);
    percentage_removal_CBF = n_CBF_zeros/n_CBF_region_voxels_tot*100;
    percentage_removal_CBF_subj(end+1)=percentage_removal_CBF;   

    V_CBF_masked(CBF_zeros)=[];
    
    %compute medians
    medians_CMRO2_subj(end+1) = median(V_CMRO2_masked);
    medians_CBF_subj(end+1) = median(V_CBF_masked);
    labels_func_subj(end+1) = regions_func(region);

    end
        
    medians_CMRO2_subjs(subj,1:n_regions_func) = medians_CMRO2_subj;
    medians_CBF_subjs(subj,1:n_regions_func) = medians_CBF_subj;
    labels_func_subjs(subj,1:n_regions_func) = labels_func_subj;
    percentage_removal_CMRO2_subjs(subj,1:n_regions_func) = percentage_removal_CMRO2_subj;
    percentage_removal_CBF_subjs(subj,1:n_regions_func) = percentage_removal_CBF_subj;

    %DWI space
    regions_dwi = unique(V_atlas_dwi(:));
    %background removal
    regions_dwi(1)=[];
    n_regions_dwi=numel(regions_dwi);

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

    %define mask
    V_fsoma_to_mask = V_fsoma;
    V_fsoma_to_mask(V_fsoma_to_mask>fsoma_threshold);

    V_pve_dwi(V_pve_dwi>pve_threshold)=1;
    V_pve_dwi(V_pve_dwi<1)=0;
    V_mask_dwi = V_pve_dwi.*V_atlas_mask_dwi;
    V_mask_dwi_fs = V_mask_dwi.*V_fsoma_to_mask;

    %mask

    V_rsoma_masked = V_rsoma.*V_mask_dwi_fs;
    V_fsoma_masked = V_fsoma.*V_mask_dwi;
    V_fneurite_masked = V_fneurite.*V_mask_dwi;  
    V_De_masked = V_De.*V_mask_dwi;
    V_Din_masked = V_Din.*V_mask_dwi;
    V_fextra_masked = V_fextra.*V_mask_dwi;
    
    % V_MSE_masked_rsoma = V_MSE.*V_mask_dwi_fs;
    % V_MSE_masked = V_MSE.*V_mask_dwi;

    %remove background
    mask_dwi_fs_zeros = find(V_mask_dwi_fs==0);
    mask_dwi_zeros = find(V_mask_dwi==0);

    V_rsoma_masked(mask_dwi_fs_zeros)=[];
    V_fsoma_masked(mask_dwi_zeros)=[];
    V_fneurite_masked(mask_dwi_zeros)=[];
    V_De_masked(mask_dwi_zeros)=[];
    V_Din_masked(mask_dwi_zeros)=[];
    V_fextra_masked(mask_dwi_zeros)=[];
    
    % V_MSE_masked_rsoma(mask_dwi_fs_zeros)=[];
    % V_MSE_masked(mask_dwi_zeros)=[];

    % %remove voxels which have MSE higher than 75th percentile
    %idx_high_MSE=find(mean_n_voxels_tot<prctile(V_MSE_masked,40));
    %idx_high_MSE_rsoma=find(mean_n_voxels_tot<prctile(V_MSE_masked_rsoma,40));
    % %FIND INDEX THAT HAVE ERROR HIGHER THAN A CERTAIN DISTRIBUTION
    % %PERCENTILE.
    % %COUNT THE PERCENTAGE
    %percentage_high_MSE_rsoma =
    %numel(idx_high_MSE_rsoma)/numel(V_rsoma_masked);

    %percentage_high_MSE_micropar =
    %numel(idx_high_MSE_micropar)/numel(V_fsoma_masked);

    % %REMOVE THEM
    % V_rsoma_masked(idx_high_MSE_rsoma)=[];
    % V_fsoma_masked(idx_high_MSE_micropar)=[];
    % V_fneurite_masked(idx_high_MSE_micropar)=[];
    % V_De_masked(idx_high_MSE_micropar)=[];
    % V_Din_masked(idx_high_MSE_micropar)=[];
    % V_fextra_masked(idx_high_MSE_micropar)=[];

    %compute medians and save results
    
    medians_rsoma_subj(end+1) = median(V_rsoma_masked);%
    medians_fsoma_subj(end+1) = median(V_fsoma_masked);%
    medians_fneurite_subj(end+1) = median(V_fneurite_masked);%
    medians_De_subj(end+1) = median(V_De_masked);%
    medians_Din_subj(end+1) = median(V_Din_masked);%
    medians_fextra_subj(end+1) = median(V_fextra_masked);%

    labels_dwi_subj(end+1) = regions_dwi(region);%   
    end
    
    labels_dwi_subjs(subj,1:n_regions_dwi) = labels_dwi_subj;%
    medians_rsoma_subjs(subj,1:n_regions_dwi) = medians_rsoma_subj;%
    medians_fsoma_subjs(subj,1:n_regions_dwi) = medians_fsoma_subj;%
    medians_fneurite_subjs(subj,1:n_regions_dwi) = medians_fneurite_subj;%
    medians_Din_subjs(subj,1:n_regions_dwi) = medians_Din_subj;%
    medians_De_subjs(subj,1:n_regions_dwi) = medians_De_subj;%
    medians_fextra_subjs(subj,1:n_regions_dwi) = medians_fextra_subj;%
    toc
    disp(strcat('Finished subject', num2str(subj),'Starting subject', num2str(subj+1)))

end
timeElapsed = toc(start_time);

%% apply thresholds
%a=labels_dwi b=medians_dwi c=labels_func d=medians_func


%%

initial_common = labels_dwi_subjs(1,:);
common = initial_common;
for i = 1:length(labels_dwi_subjs(:,1))
    commonElements = intersect(common, labels_dwi_subjs(i,:));
    common = commonElements;
end

%First, detect indices of common Elements
commonElements(commonElements==0)=[];
lst_idx_matrix=[];
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
    lst_idx_matrix(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_rsoma_subjs_final=[];
for row = 1:length(medians_rsoma_subjs(:,1))
    medians_rsoma_subjs_row=medians_rsoma_subjs(row,:);
    lst_idx = lst_idx_matrix(row,:);
    medians_rsoma_subjs_row_final=medians_rsoma_subjs_row(lst_idx);
    medians_rsoma_subjs_final(row,:)=medians_rsoma_subjs_row_final;
end

labels_dwi_subjs_final=[];
for row = 1:length(labels_dwi_subjs(:,1))
    labels_dwi_subjs_row=labels_dwi_subjs(row,:);
    lst_idx = lst_idx_matrix(row,:);
    labels_dwi_subjs_row_final=labels_dwi_subjs_row(lst_idx);
    labels_dwi_subjs_final(row,:)=labels_dwi_subjs_row_final;
end
%apply everything both to the dwi space and func space, so also c and d

initial_common = labels_func_subjs(1,:);
common = initial_common;
for i = 1:length(labels_func_subjs(:,1))
    commonElements = intersect(common, labels_func_subjs(i,:));
    common = commonElements;
end

%First, detect indices of common Elements
commonElements(commonElements==0)=[];
lst_idx_matrix=[];
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
    lst_idx_matrix(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_CMRO2_subjs_final=[];
for row = 1:length(medians_CMRO2_subjs(:,1))
    medians_CMRO2_subjs_row=medians_CMRO2_subjs(row,:);
    lst_idx = lst_idx_matrix(row,:);
    medians_CMRO2_subjs_row_final=medians_CMRO2_subjs_row(lst_idx);
    medians_CMRO2_subjs_final(row,:)=medians_CMRO2_subjs_row_final;
end

labels_func_subjs_final=[];
for row = 1:length(labels_func_subjs(:,1))
    labels_func_subjs_row=labels_func_subjs(row,:);
    lst_idx = lst_idx_matrix(row,:);
    labels_func_subjs_row_final=labels_func_subjs_row(lst_idx);
    labels_func_subjs_final(row,:)=labels_func_subjs_row_final;
end

%then for each subj (row of the two matrices),
%find common elements (regions) between a (dwi) and b (func) space,
%find the idx 
%keep only those common both in DWI and func space
%create a new dwi and func matrix which will have same size.
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

%First, detect indices of common Elements

lst_idx_matrix=[];
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
    lst_idx_matrix(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_rsoma_subjs_final_spaces=[];
for row = 1:length(medians_rsoma_subjs_final(:,1))
    medians_rsoma_subjs_row=medians_rsoma_subjs_final(row,:);
    lst_idx = lst_idx_matrix(row,:);
    medians_rsoma_subjs_row_final=medians_rsoma_subjs_row(lst_idx);
    medians_rsoma_subjs_final_spaces(row,:)=medians_rsoma_subjs_row_final;
end

labels_dwi_subjs_final_spaces=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    labels_dwi_subjs_row=labels_dwi_subjs_final(row,:);
    lst_idx = lst_idx_matrix(row,:);
    labels_dwi_subjs_row_final=labels_dwi_subjs_row(lst_idx);
    labels_dwi_subjs_final_spaces(row,:)=labels_dwi_subjs_row_final;
end

%do the analogous for c and d

%First, detect indices of common Elements

lst_idx_matrix=[];
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
    lst_idx_matrix(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_CMRO2_subjs_final_spaces=[];
for row = 1:length(medians_CMRO2_subjs_final(:,1))
    medians_CMRO2_subjs_row=medians_CMRO2_subjs_final(row,:);
    lst_idx = lst_idx_matrix(row,:);
    medians_CMRO2_subjs_row_final=medians_CMRO2_subjs_row(lst_idx);
    medians_CMRO2_subjs_final_spaces(row,:)=medians_CMRO2_subjs_row_final;
end

labels_func_subjs_final_spaces=[];
for row = 1:length(labels_func_subjs_final(:,1))
    labels_func_subjs_row=labels_func_subjs_final(row,:);
    lst_idx = lst_idx_matrix(row,:);
    labels_func_subjs_row_final=labels_func_subjs_row(lst_idx);
    labels_func_subjs_final_spaces(row,:)=labels_func_subjs_row_final;
end

%%

medians_CMRO2 = medians_CMRO2_subjs_final_spaces;
medians_rsoma = medians_rsoma_subjs_final_spaces;

medians_CMRO2_vec = nanmedian(medians_CMRO2,1);
medians_rsoma_vec = nanmedian(medians_rsoma,1);

idx_outlier=find(medians_CMRO2_vec>270);

medians_CMRO2_vec(idx_outlier)=[];
medians_rsoma_vec(idx_outlier)=[];

figure, scatter(medians_rsoma_vec,medians_CMRO2_vec);

[r,p]=corrcoef(medians_CMRO2_vec,medians_rsoma_vec,'rows','complete');%0.37%p0.0029
