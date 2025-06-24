%Load GM 
img_path_GM = '/storage/shared/Atlas/atlas_GM_on_MNI152_T1_2mm.nii.gz'; %GM resampled on T1 2mm
%atlas_GM resampled on SANDI and CBF dimensions
Vhdr_GM = spm_vol(img_path_GM);
V_GM_tot = spm_read_vols(Vhdr_GM); 

%Loading atlas
img_path_atlas='/storage/shared/Atlas/AAL3v1_2mm_resampled.nii.gz';
Vhdr = spm_vol(img_path_atlas);
V_atlas_tot = spm_read_vols(Vhdr);

%% Import data 

run='run-01';%CHANGE
%%%%%%%%%%%%%

subjects = importdata(strcat('/media/nas_rete/Vitality/code/subjs_DWI.txt'));

%run2
%subjects([11,27])=[];

n_subjs=length(subjects);

%%%%load data
V_Din_tots={};
V_De_tots={};
V_rsoma_tots={};
V_fsoma_tots={};

V_fextra_tots={};
V_fneurite_tots={};
V_CMRO2_tots={};
V_CBF_tots={};





for i = 1:1:n_subjs%lst



    subj=num2str(subjects{i});
%     img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/250101/CBF0toMNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');
%     img_path_CMRO2 = strcat('/media/nas_rete/Vitality/maps2MNI/250101/CMRO20toMNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CMRO2_map_2MNI2mm.nii.gz');
% 
%     img_path_rsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_Rsoma_2MNI2mm.nii.gz');
%     img_path_fsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_fsoma_2MNI2mm.nii.gz');
%     %img_path_fc = strcat('/media/nas_rete/Vitality/maps2MNI/SANDItoMNI/sub-00',i,'_run-01_SANDI-fit_fc_2MNI2mm.nii.gz');
%     img_path_De = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_De_2MNI2mm.nii.gz');
%     img_path_Din = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_Din_2MNI2mm.nii.gz');
%     img_path_fneurite = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_fneurite_2MNI2mm.nii.gz');
%     img_path_fextra = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_fextra_2MNI2mm.nii.gz');

     img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/250418_corrected_for_RIM/CBF02MNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CBF_map_2MNI.nii.gz');
    

    img_path_CMRO2 = strcat('/media/nas_rete/Vitality/registered/perf/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CMRO2_map.nii.gz'); %_2MNI2mm
    
    img_path_rsoma = strcat('/media/nas_rete/Vitality/registered/SANDI/',subj,'_',run,'_SANDI-fit_Rsoma_2mm.nii.gz');
    img_path_fsoma = strcat('/media/nas_rete/Vitality/registered/SANDI/',subj,'_',run,'_SANDI-fit_fsoma_2mm.nii.gz');


    %img_path_fc = strcat('/media/nas_rete/Vitality/maps2MNI/SANDItoMNI/sub-00',i,'_run-01_SANDI-fit_fc_2MNI2mm.nii.gz');
    img_path_De = strcat('/media/nas_rete/Vitality/registered/SANDI/',subj,'_',run,'_SANDI-fit_De_2mm.nii.gz');
    img_path_Din = strcat('/media/nas_rete/Vitality/registered/SANDI/',subj,'_',run,'_SANDI-fit_Din_2mm.nii.gz');
    img_path_fneurite = strcat('/media/nas_rete/Vitality/registered/SANDI/',subj,'_',run,'_SANDI-fit_fneurite_2mm.nii.gz');
    img_path_fextra = strcat('/media/nas_rete/Vitality/registered/SANDI/',subj,'_',run,'_SANDI-fit_fextra_2mm.nii.gz');

 

    Vhdr = spm_vol(img_path_CBF);
    V_CBF_tot = spm_read_vols(Vhdr);
%     figure, imagesc(rot90(squeeze(V_CBF_tot(30,:,:)))) 
%     imshow(rot90(squeeze(V(30,:,:))),[])
    V_CBF_tots{end+1} = V_CBF_tot;

    Vhdr = spm_vol(img_path_CMRO2);
    V_CMRO2_tot = spm_read_vols(Vhdr);
%     figure, imagesc(rot90(squeeze(V_CBF_tot(30,:,:)))) 
%     imshow(rot90(squeeze(V(30,:,:))),[])
    V_CMRO2_tots{end+1} = V_CMRO2_tot;

    Vhdr = spm_vol(img_path_rsoma);
    V_rsoma_tot = spm_read_vols(Vhdr);
    V_rsoma_tots{end+1} = V_rsoma_tot;

    Vhdr = spm_vol(img_path_fsoma);
    V_fsoma_tot = spm_read_vols(Vhdr);
    V_fsoma_tots{end+1} = V_fsoma_tot;

    Vhdr = spm_vol(img_path_fneurite);
    V_fneurite_tot = spm_read_vols(Vhdr);
    V_fneurite_tots{end+1} = V_fneurite_tot;

    Vhdr = spm_vol(img_path_fextra);
    V_fextra_tot = spm_read_vols(Vhdr);
    V_fextra_tots{end+1} = V_fextra_tot;

    Vhdr = spm_vol(img_path_Din);
    V_Din_tot = spm_read_vols(Vhdr);
    V_Din_tots{end+1} = V_Din_tot;

    Vhdr = spm_vol(img_path_De);
    V_De_tot = spm_read_vols(Vhdr);
    V_De_tots{end+1} = V_De_tot;

end

%% number of cells density map (many subjects) 

V_fc_tots={};

for i = 1:n_subjs    
    V_rsoma = V_rsoma_tots{i};
    V_fsoma = V_fsoma_tots{i};

    
    for i = 1:length(V_rsoma(:))
        if V_rsoma(i) < 5 %check
            V_rsoma(i)=0;
        end
    end
    
    
    %convert to m^3
    V_rsoma = V_rsoma.*10^-6;
    
    %voxel wise divide fs map over 4/3pir^3
    fc_map = V_fsoma./((4/3)*pi*V_rsoma.^3);
    %figure, imagesc(fc_map(:,:,45));
    for i = 1:length(fc_map(:))
        if fc_map(i)==Inf
            fc_map(i)=NaN;%VALUTA SE METTERE A 0.
    %     elseif isnan(fc_map(i))
    %         fc_map(i)=0;
        end
    end

    V_fc_tots{end+1} = fc_map;

end

%% superficial density (many subjects)

V_fsup_tots={};

for i = 1:n_subjs    
    V_rsoma = V_rsoma_tots{i};
    V_fsoma = V_fsoma_tots{i};
    
%     V_fsoma_mask=V_fsoma;
%     V_fsoma_mask(V_fsoma_mask>0.35)=1;
%     V_fsoma_mask(V_fsoma_mask<1)=0;
% 
%     V_rsoma = V_rsoma.*V_fsoma_mask;

    
    for i = 1:length(V_rsoma(:))
        if V_rsoma(i) < 5 %check
            V_rsoma(i)=0;
        end
    end
%     
    
    %convert to m^3
    V_rsoma = V_rsoma.*10^-6;%Should it be 10^-6? Why is it measured as mm?
    
    %voxel wise divide fs map over 4/3pir^3
    fsup_map = V_fsoma./V_rsoma;
    fsup_map = 3*fsup_map;
    %figure, imagesc(fc_map(:,:,45));
    for i = 1:length(fsup_map(:))
        if fsup_map(i)==Inf
            fsup_map(i)=NaN;%VALUTA SE METTERE A 0.
    %     elseif isnan(fc_map(i))
    %         fc_map(i)=0;
        end
    end

    V_fsup_tots{end+1} = fsup_map;

end

% V_fsup_one=V_fsup_tots{1};
% figure, imagesc(rot90(V_fsup_one(:,:,45)));
% title('Superficial Soma Density map')
% V_fsup_one_array=V_fsup_one(:);
% 
% figure, hist(V_fsup_one_array);
% title('Superficial Soma Density Distribution');
% grid on

%% Method IV (remove zeros of binary maps in all maps and CBF<=0, CMRO2<=0 in all maps)
V_GM = V_GM_tot;
V_GM(V_GM>0.5)=1;
V_GM(V_GM<1)=0;


regions = unique(V_atlas_tot(:));
n_regions = numel(regions);

medians_CBF = [];
medians_CMRO2 = [];
medians_fc= [];
medians_fextra= [];
medians_fneurite = [];
medians_fsoma= [];
medians_Din = [];
medians_De =[];
medians_rsoma = [];
medians_fsup = [];


reshaped_medians_CBF=[];
reshaped_medians_CMRO2=[];
reshaped_medians_De = [];
reshaped_medians_rsoma = [];
reshaped_medians_fsoma = [];
reshaped_medians_fneurite = [];
reshaped_medians_fextra = [];
reshaped_medians_Din = [];
reshaped_medians_fc = [];
reshaped_medians_fsup = [];

n_voxels_lst = [];
n_voxels_lst_allsubjs=[];

corr_for_each_subj_rsoma_CBF=[];
pvalue_for_each_subj_rsoma_CBF=[];
corr_for_each_subj_fsoma_CBF=[];
pvalue_for_each_subj_fsoma_CBF=[];
corr_for_each_subj_fsup_CBF=[];
pvalue_for_each_subj_fsup_CBF=[];

corr_for_each_subj_rsoma_CMRO2=[];
pvalue_for_each_subj_rsoma_CMRO2=[];
corr_for_each_subj_fsoma_CMRO2=[];
pvalue_for_each_subj_fsoma_CMRO2=[];
corr_for_each_subj_fsup_CMRO2=[];
pvalue_for_each_subj_fsup_CMRO2=[];

CBF_zeros_brain_subjs=[];
CMRO2_zeros_brain_subjS=[];

tic
for i = 1:1:n_subjs %1:1:length(lst)
    
    V_CBF_tot = V_CBF_tots{i};
    V_CMRO2_tot = V_CMRO2_tots{i};
    V_fc_tot = V_fc_tots{i};
    V_fsoma_tot = V_fsoma_tots{i};
    V_fneurite_tot = V_fneurite_tots{i};
    V_fextra_tot = V_fextra_tots{i};
    V_Din_tot = V_Din_tots{i};
    V_De_tot = V_De_tots{i};
    V_rsoma_tot = V_rsoma_tots{i};
    V_fsup_tot = V_fsup_tots{i};

    medians_CBF=[];
    medians_CMRO2=[];
    medians_fsoma=[];
    medians_fc=[];
    medians_fneurite=[];
    medians_fextra=[];
    medians_Din=[];
    medians_De=[];
    medians_rsoma=[];
    medians_fsup=[];

    n_voxels_lst=[];

    for k = 2:numel(regions)
        tic


        V_atlas = V_atlas_tot;

        for ii = 1:length(V_atlas_tot(:))
            if V_atlas_tot(ii)==regions(k)
                V_atlas(ii)=1;
            else
                V_atlas(ii)=0;
            end
        end



        
        V_CBF = V_CBF_tot;
        V_CMRO2 = V_CMRO2_tot;
       
        V_fc = V_fc_tot;
        V_fextra = V_fextra_tot;
        V_fneurite = V_fneurite_tot;
        V_rsoma = V_rsoma_tot;
        V_fsoma = V_fsoma_tot;
        V_Din = V_Din_tot;
        V_De = V_De_tot;
        V_fsup = V_fsup_tot;
        V_fsoma_to_mask=V_fsoma_tot;

        
        V_fsoma_to_mask(V_fsoma_to_mask>0.15)=1;
        V_fsoma_to_mask(V_fsoma_to_mask<1)=0;


        V_mask = V_GM.*V_fsoma_to_mask.*V_atlas;

        V_CBF_withoutzerosatlas=V_CBF;
        V_CMRO2_withoutzerosatlas=V_CMRO2;

        mask_zeros=find(V_mask==0);
        
        %check
        V_CBF_withoutzerosatlas(mask_zeros)=[];
        V_CMRO2_withoutzerosatlas(mask_zeros)=[];
        %figure,hist(V_CBF_withoutzerosatlas(:));
        %figure, hist(V_CMRO2_withoutzerosatlas(:));

        V_CBF_masked = V_CBF.*V_mask;
        V_CMRO2_masked = V_CMRO2.*V_mask;
        V_fc_masked = V_fc.*V_mask;
        V_fextra_masked = V_fextra.*V_mask;
        V_fneurite_masked = V_fneurite.*V_mask;
        V_rsoma_masked = V_rsoma.*V_mask;
        V_fsoma_masked = V_fsoma.*V_mask;
        V_Din_masked = V_Din.*V_mask;
        V_De_masked = V_De.*V_mask;
        V_fsup_masked = V_fsup.*V_mask;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp(strcat('region', num2str(k)))

        v_CBF_masked=V_CBF_masked;
        v_CMRO2_masked=V_CMRO2_masked;
        v_fc_masked=V_fc_masked;
        v_fsoma_masked=V_fsoma_masked;
        v_fsup_masked=V_fsup_masked;
        v_rsoma_masked=V_rsoma_masked;
        v_fneurite_masked=V_fneurite_masked;
        v_fextra_masked=V_fextra_masked;
        v_Din_masked=V_Din_masked;
        v_De_masked=V_De_masked;

        v_CBF_masked(mask_zeros)=[];
        v_CMRO2_masked(mask_zeros)=[];
        v_fc_masked(mask_zeros)=[];
        v_fsoma_masked(mask_zeros)=[];
        v_rsoma_masked(mask_zeros)=[];
        v_fsup_masked(mask_zeros)=[];
        v_fneurite_masked(mask_zeros)=[];
        v_fextra_masked(mask_zeros)=[];
        v_Din_masked(mask_zeros)=[];
        v_De_masked(mask_zeros)=[];

        indices_CBF = [];
        for ii = 1:numel(v_CBF_masked)
            if v_CBF_masked(ii)<=0 %|| V_energy_masked(ii)>up_thr %50 200
                indices_CBF(end+1)=ii;
            end
        end

        indices_CMRO2 = [];
        for ii = 1:numel(v_CMRO2_masked)
            if v_CMRO2_masked(ii)<=0 %|| V_energy_masked(ii)>up_thr %50 200
                indices_CMRO2(end+1)=ii;
            end
        end

        disp(strcat('region', num2str(k)))

        indices = cat(2, indices_CBF, indices_CMRO2);
        indices_unique = unique(indices);

        v_CBF_masked(indices_unique)=[];
        v_CMRO2_masked(indices_unique)=[];
        v_fc_masked(indices_unique)=[];
        v_fsoma_masked(indices_unique)=[];
        v_rsoma_masked(indices_unique)=[];
        v_fsup_masked(indices_unique)=[];
        v_fneurite_masked(indices_unique)=[];
        v_fextra_masked(indices_unique)=[];
        v_Din_masked(indices_unique)=[];
        v_De_masked(indices_unique)=[];



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        n_voxels = numel(v_CBF_masked);%In this case it is not the same for all parameters
      
        n_voxels_lst(end+1) = n_voxels;

        median_CBF = median(v_CBF_masked,'omitnan');
        median_CMRO2 = median(v_CMRO2_masked,'omitnan');
        median_fc = median(v_fc_masked,'omitnan');
        median_fextra = median(v_fextra_masked,'omitnan');
        median_fneurite = median(v_fneurite_masked,'omitnan');
        median_fsoma = median(v_fsoma_masked,'omitnan');
        median_rsoma = median(v_rsoma_masked,'omitnan');
        median_Din = median(v_Din_masked,'omitnan');
        median_De = median(v_De_masked,'omitnan');
        median_fsup=median(v_fsup_masked,'omitnan');

        medians_CBF(end+1) = median_CBF;
        medians_CMRO2(end+1) = median_CMRO2;
        medians_fc(end+1) = median_fc;
        medians_fextra(end+1) = median_fextra;
        medians_fneurite(end+1) = median_fneurite;
        medians_fsoma(end+1) = median_fsoma;
        medians_Din(end+1) = median_Din;
        medians_De(end+1) = median_De;
        medians_rsoma(end+1) = median_rsoma;
        medians_fsup(end+1) = median_fsup;

        
    toc
    end

    n_voxels_lst_allsubjs(i,:) = n_voxels_lst;

    reshaped_medians_CBF(i,:) = medians_CBF;
    reshaped_medians_CMRO2(i,:) = medians_CMRO2;
    reshaped_medians_fc(i,:) = medians_fc;
    reshaped_medians_fsoma(i,:) = medians_fsoma;
    reshaped_medians_rsoma(i,:) = medians_rsoma;
    reshaped_medians_fextra(i,:) = medians_fextra;
    reshaped_medians_fneurite(i,:) = medians_fneurite;
    reshaped_medians_Din(i,:) = medians_Din;
    reshaped_medians_De(i,:) = medians_De;
    reshaped_medians_fsup(i,:) = medians_fsup;

    [r_rsoma_CBF,p_rsoma_CBF] = corrcoef(medians_CBF, medians_rsoma, 'rows','complete');
    [r_fsoma_CBF,p_fsoma_CBF] = corrcoef(medians_CBF, medians_fsoma, 'rows','complete');
    [r_fsup_CBF,p_fsup_CBF] = corrcoef(medians_CBF, medians_fsup, 'rows','complete');

    [r_rsoma_CMRO2,p_rsoma_CMRO2] = corrcoef(medians_CMRO2, medians_rsoma, 'rows','complete');
    [r_fsoma_CMRO2,p_fsoma_CMRO2] = corrcoef(medians_CMRO2, medians_fsoma, 'rows','complete');
    [r_fsup_CMRO2,p_fsup_CMRO2] = corrcoef(medians_CMRO2, medians_fsup, 'rows','complete');

    corr_for_each_subj_rsoma_CBF(i)=r_rsoma_CBF(2);
    pvalue_for_each_subj_rsoma_CBF(i)=p_rsoma_CBF(2);

    corr_for_each_subj_fsoma_CBF(i)=r_fsoma_CBF(2);
    pvalue_for_each_subj_fsoma_CBF(i)=p_fsoma_CBF(2);

    corr_for_each_subj_fsup_CBF(i)=r_fsup_CBF(2);
    pvalue_for_each_subj_fsup_CBF(i)=p_fsup_CBF(2);

    corr_for_each_subj_rsoma_CMRO2(i)=r_rsoma_CMRO2(2);
    pvalue_for_each_subj_rsoma_CMRO2(i)=p_rsoma_CMRO2(2);

    corr_for_each_subj_fsoma_CMRO2(i)=r_fsoma_CMRO2(2);
    pvalue_for_each_subj_fsoma_CMRO2(i)=p_fsoma_CMRO2(2);

    corr_for_each_subj_fsup_CMRO2(i)=r_fsup_CMRO2(2);
    pvalue_for_each_subj_fsup_CMRO2(i)=p_fsup_CMRO2(2);

    disp(strcat('Finished subject', num2str(i),'Starting subject', num2str(i+1)))

end
timeElapsed = toc;

%% select parameters to investigate
par='Rsoma'; 
%the energy consumption parameter has to be selected for a plotting reason
%(both CBF and CMRO2 are processed anyway)
energy_par='CMRO2';

%% for each subj 


if strcmp(par,'Rsoma') && strcmp(energy_par,'CBF')
    corr_for_each_subj=corr_for_each_subj_rsoma_CBF;
elseif strcmp(par,'fsoma') && strcmp(energy_par,'CBF')
    corr_for_each_subj=corr_for_each_subj_fsoma_CBF;
elseif strcmp(par,'fsup') && strcmp(energy_par,'CBF')
    corr_for_each_subj=corr_for_each_subj_fsup_CBF;
elseif strcmp(par,'Rsoma') && strcmp(energy_par,'CMRO2')
    corr_for_each_subj=corr_for_each_subj_rsoma_CMRO2;
elseif strcmp(par,'fsoma') && strcmp(energy_par,'CMRO2')
    corr_for_each_subj=corr_for_each_subj_fsup_CMRO2;
elseif strcmp(par,'fsup') && strcmp(energy_par,'CMRO2')
    corr_for_each_subj=corr_for_each_subj_fsup_CMRO2;
end

mean_with_subjs_corr=nanmean(corr_for_each_subj);%why do we have NaNs?
[h,p,ci,stats]=ttest(atanh(corr_for_each_subj));


mean_corr=round(mean_with_subjs_corr,2);
mean_corr=num2str(mean_corr);

pvalue=num2str(round(p,2));


figure, 
s=histogram(corr_for_each_subj,'FaceAlpha',1,'BinWidth',0.07);
s.FaceColor="b";
xlabel('correlation coefficient, r','FontWeight','bold','FontSize',15);
ylabel('Counts (# subjects)','FontWeight','bold','FontSize',15);
ylim([0,14]);
xline(0,'--','LineWidth',3);
if p<0.05 && p>0.01    
    txt = {strcat('\mu_r = ',mean_corr,'*')};
elseif p<0.01 && p>0.001   
    txt = {strcat('\mu_r = ',mean_corr,'**')};
else
    txt = {strcat('\mu_r = ',mean_corr,'***')};
end

title(strcat(energy_par, 'vs',par));
text(0.5,7,txt, 'FontWeight', 'bold','FontSize',15);
grid on
%% average across subjects

mean_CMRO2_tot = nanmean(reshaped_medians_CMRO2,1);
mean_CBF_tot = nanmean(reshaped_medians_CBF,1);
mean_fc_tot = nanmean(reshaped_medians_fc,1);
mean_fextra_tot = nanmean(reshaped_medians_fextra,1);
mean_fneurite_tot = nanmean(reshaped_medians_fneurite,1);
mean_fsoma_tot = nanmean(reshaped_medians_fsoma,1);
mean_rsoma_tot = nanmean(reshaped_medians_rsoma,1);
mean_Din_tot = nanmean(reshaped_medians_Din,1);
mean_De_tot = nanmean(reshaped_medians_De,1);
mean_fsup_tot = nanmean(reshaped_medians_fsup,1);

n=numel(reshaped_medians_CMRO2(:,1));%you can choose any parameter

SE_CMRO2_tot = nanstd(reshaped_medians_CMRO2,0,1)/sqrt(n);
SE_CBF_tot = nanstd(reshaped_medians_CBF,0,1)/sqrt(n);
SE_fc_tot = nanstd(reshaped_medians_fc,0,1)/sqrt(n);
SE_fextra_tot = nanstd(reshaped_medians_fextra,0,1)/sqrt(n);
SE_fneurite_tot = nanstd(reshaped_medians_fneurite,0,1)/sqrt(n);
SE_fsoma_tot = nanstd(reshaped_medians_fsoma,0,1)/sqrt(n);
SE_rsoma_tot = nanstd(reshaped_medians_rsoma,0,1)/sqrt(n);
SE_Din_tot = nanstd(reshaped_medians_Din,0,1)/sqrt(n);
SE_De_tot = nanstd(reshaped_medians_De,0,1)/sqrt(n);
SE_fsup_tot = nanstd(reshaped_medians_fsup,0,1)/sqrt(n);

var_rsoma=SE_rsoma_tot./mean_rsoma_tot;
var_fsoma=SE_fsoma_tot./mean_fsoma_tot;
var_fc=SE_fc_tot./mean_fc_tot;
var_fneurite=SE_fneurite_tot./mean_fneurite_tot;
var_fextra=SE_fextra_tot./mean_fextra_tot;
var_Din=SE_Din_tot./mean_Din_tot;
var_De=SE_De_tot./mean_De_tot;
var_fsup=SE_fsup_tot./mean_fsup_tot;

mean_n_voxels_tot = mean(n_voxels_lst_allsubjs,1);
regions=1:numel(mean_CMRO2_tot);
mean_n_voxels=mean_n_voxels_tot;
% 
% figure, bar(regions,mean_n_voxels_tot);
% xlabel('regions','FontSize',15);
% ylabel('mean number of voxels','FontSize',15);


idx_low_n_voxels=[];
for i=1:numel(regions)
    if mean_n_voxels_tot(i)<prctile(mean_n_voxels_tot,40)%mean(mean_n_voxels)/5 %CHECK
        idx_low_n_voxels(end+1)=i;
    end
end


labels_tot=unique(V_atlas_tot(:));
labels=labels_tot;
labels(idx_low_n_voxels)=[];
mean_n_voxels=mean_n_voxels_tot;
mean_n_voxels(idx_low_n_voxels)=[];

mean_CBF_tot(idx_low_n_voxels)=[];
mean_CMRO2_tot(idx_low_n_voxels)=[];
mean_fc_tot(idx_low_n_voxels)=[];
mean_fextra_tot(idx_low_n_voxels)=[];
mean_fneurite_tot(idx_low_n_voxels)=[];
mean_fsoma_tot(idx_low_n_voxels)=[];
mean_rsoma_tot(idx_low_n_voxels)=[];
mean_Din_tot(idx_low_n_voxels)=[];
mean_De_tot(idx_low_n_voxels)=[];
mean_fsup_tot(idx_low_n_voxels)=[];

%this is to plot errors on both x and y values
SE_rsoma_tot(idx_low_n_voxels)=[];
SE_fsoma_tot(idx_low_n_voxels)=[];
SE_fsup_tot(idx_low_n_voxels)=[];
SE_CBF_tot(idx_low_n_voxels)=[];
SE_CMRO2_tot(idx_low_n_voxels)=[];

%this is to calculate coefficient of variation 
var_rsoma(idx_low_n_voxels)=[];
var_fsoma(idx_low_n_voxels)=[];
var_fextra(idx_low_n_voxels)=[];
var_fneurite(idx_low_n_voxels)=[];
var_fc(idx_low_n_voxels)=[];
var_De(idx_low_n_voxels)=[];
var_Din(idx_low_n_voxels)=[];
var_fsup(idx_low_n_voxels)=[];

%idx_high_SE=find(SE_rsoma_tot>prctile(SE_rsoma_tot,95));
%For now we used GM parameters
idx_high_var_rsoma=find(var_rsoma>prctile(var_rsoma,75));
idx_high_var_fsoma=find(var_fsoma>prctile(var_fsoma,75));
idx_high_var_fsup=find(var_fsup>prctile(var_fsup,75));

%% in order to do Energy vs one Microparameter analysis (run if corr analysis)

if strcmp(par,'Rsoma')
    mean_CBF_tot(idx_high_var_rsoma)=[];
    mean_CMRO2_tot(idx_high_var_rsoma)=[];
    mean_rsoma_tot(idx_high_var_rsoma) = [];
    SE_rsoma_tot(idx_high_var_rsoma)=[];
    SE_CBF_tot(idx_high_var_rsoma)=[];
    SE_CMRO2_tot(idx_high_var_rsoma)=[];
elseif strcmp(par,'fsoma')
    mean_CBF_tot(idx_high_var_fsoma)=[];
    mean_CMRO2_tot(idx_high_var_fsoma)=[];
    mean_fsoma_tot(idx_high_var_fsoma)=[];
    SE_fsoma_tot(idx_high_var_fsoma)=[];
    SE_CBF_tot(idx_high_var_fsoma)=[];
    SE_CMRO2_tot(idx_high_var_fsoma)=[];
elseif strcmp(par,'fsup')
    mean_CBF_tot(idx_high_var_fsup)=[];
    mean_CMRO2_tot(idx_high_var_fsup)=[];
    mean_fsup_tot(idx_high_var_fsup)=[];
    SE_fsup_tot(idx_high_var_fsup)=[];
    SE_CBF_tot(idx_high_var_fsup)=[];  
    SE_CMRO2_tot(idx_high_var_fsup)=[]; 
end

%% in order to do Energy vs many Microparameters analysis (run if GLM analysis)
%it must be so in order to run this analysis

idx_tot=unique([idx_high_var_rsoma,idx_high_var_fsoma,idx_high_var_fsup]);

mean_CBF_tot(idx_tot)=[];
mean_CMRO2_tot(idx_tot)=[];
mean_rsoma_tot(idx_tot)=[];
mean_fsoma_tot(idx_tot)=[];
mean_fsup_tot(idx_tot)=[];

%% Matrix for correlation across subjs analysis (region by region)

V_atlas=V_atlas_tot;
labels_tot=unique(V_atlas);
labels_tot(labels_tot==0)=[];

reshaped_medians_CBF_tot=reshaped_medians_CBF;
reshaped_medians_CMRO2_tot=reshaped_medians_CMRO2;
reshaped_medians_fsoma_tot=reshaped_medians_fsoma;
reshaped_medians_fsup_tot=reshaped_medians_fsup;
reshaped_medians_rsoma_tot=reshaped_medians_rsoma;

reshaped_medians_CBF_tot(:,idx_low_n_voxels)=[];
reshaped_medians_CMRO2_tot(:,idx_low_n_voxels)=[];
labels_tot(idx_low_n_voxels)=[];

if strcmp(par,'Rsoma')
    reshaped_medians_rsoma_tot(:,idx_low_n_voxels)=[];
    labels_tot(idx_high_var_rsoma)=[];
    reshaped_medians_CBF_tot(:,idx_high_var_rsoma)=[];
    reshaped_medians_CMRO2_tot(:,idx_high_var_rsoma)=[];
    reshaped_medians_rsoma_tot(:,idx_high_var_rsoma)=[];
elseif strcmp(par,'fsoma')
    reshaped_medians_fsoma_tot(:,idx_low_n_voxels)=[];
    labels_tot(idx_high_var_fsoma)=[];
    reshaped_medians_CBF_tot(:,idx_high_var_fsoma)=[];
    reshaped_medians_CMRO2_tot(:,idx_high_var_fsoma)=[];
    reshaped_medians_fsoma_tot(:,idx_high_var_fsoma)=[];
elseif strcmp(par,'fsup')
    reshaped_medians_fsup_tot(:,idx_low_n_voxels)=[];
    labels_tot(idx_high_var_fsup)=[];
    reshaped_medians_CBF_tot(:,idx_high_var_fsup)=[];
    reshaped_medians_CMRO2_tot(:,idx_high_var_fsup)=[];
    reshaped_medians_fsup_tot(:,idx_high_var_fsup)=[];
end

%% plot energy vs microparameter 

if strcmp(par,'Rsoma')
    mean_SANDI=mean_rsoma_tot;
    SE_SANDI_tot=SE_rsoma_tot;
elseif strcmp(par,'fsoma')
    mean_SANDI=mean_fsoma_tot;
    SE_SANDI_tot=SE_fsoma_tot;
elseif strcmp(par,'fsup')
    mean_SANDI=mean_fsup_tot;
    SE_SANDI_tot=SE_fsup_tot;
end

if strcmp(energy_par,'CBF')
    mean_energy_tot=mean_CBF_tot;
    dependent_parameter='CBF (ml/100g/min)';
    SE_energy_tot=SE_CBF_tot;
elseif strcmp(energy_par,'CMRO2')
    mean_energy_tot=mean_CMRO2_tot;
    dependent_parameter='CMRO_2 (\mu mol/100g/min)';
    SE_energy_tot=SE_CMRO2_tot;
end

% if any(isnan(mean_SANDI))%fsup can have some NaN values
%     idx=find(isnan(mean_SANDI));
%     mean_SANDI(idx)=[];
%     mean_energy_nanparremoved=mean_energy_tot;
%     mean_energy_nanparremoved(idx)=[];
%     SE_energy_tot(idx)=[];
%     SE_SANDI_tot(idx)=[];
% else
%     mean_energy_nanparremoved=mean_energy_tot;
% end
% % mean_fsoma=nanmean(reshaped_medians_fsoma,1);
% % mean_energy=nanmean(reshaped_medians_energy,1);
% % 
% % n
% % SE_energy = nanstd(reshaped_medians_energy,0,1)/sqrt(numel(reshaped_medians_energy(:,1)));
% % SE_SANDI = nanstd(reshaped_medians_fsoma,0,1)/sqrt(numel(reshaped_medians_energy(:,1)));

[r,p] = corrcoef(mean_energy_tot, mean_SANDI, 'rows','complete');

corr_coef = round(r(2),2);
p_value = p(2);
corr_coef_str = num2str(corr_coef);
p_value_str = num2str(p_value);

P = polyfit(mean_SANDI,mean_energy_tot,1);
yfit_energy = P(1)*mean_SANDI+P(2);




figure, 
s = errorbar(mean_SANDI, mean_energy_tot, SE_energy_tot, SE_energy_tot, SE_SANDI_tot, SE_SANDI_tot,'o');
hold on
plot(mean_SANDI,yfit_energy,'--','LineWidth',3,'Color',"#000000");
xlabel(strcat(par),'FontSize',15,'FontWeight','bold');
ylabel(dependent_parameter,'FontSize',15,'FontWeight','bold');
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
%n_subjs_str=num2str(n_subjs);
if p(2)<0.05 && p(2)>0.01    
    txt = {strcat('r = ',corr_coef_str,'*')};
elseif p(2)<0.01 && p(2)>0.001   
    txt = {strcat('r = ',corr_coef_str,'**')};
elseif p(2)<0.001
    txt = {strcat('r = ',corr_coef_str,'***')};
end
% text(60,12,txt,'FontWeight', 'Bold');
%annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold','FontSize',15);
if strcmp(energy_par,'CBF') && strcmp(par,'Rsoma')
    text(13.8,60,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_par,'CBF') && strcmp(par,'fsoma')
    text(0.39,60,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_par,'CMRO2') && strcmp(par,'Rsoma')
    text(13.8,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_par,'CMRO2') && strcmp(par,'fsoma')
    text(0.39,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_par,'CMRO2') && strcmp(par,'fsup')
    text(8*10^4,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_par,'CBF') && strcmp(par,'fsup')
    text(8*10^4,60,txt, 'FontWeight', 'bold','FontSize',12);
end
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'box','off')
x0=400;
y0=400;
width=550;
height=450;
set(gcf,'position',[x0,y0,width,height]);
grid on

%% check correlation distribution 

%correlation between CBF and SANDI in each region

% %this is done for all the "original" regions.
% % you could remove columns corresponding to regions removed in reshaped_medians_CBF
% %and reshaped_medians_SANDI
% for i = 1:(numel(labels_tot))
%     [r,p] = corrcoef(reshaped_medians_energy_tot(:,i), reshaped_medians_SANDI_tot(:,i),'rows','complete');
%     z(i)=r(2);
%     pvalue(i)=p(2);
% end_tot
micro_parameter='Rsoma';
energy='CMRO2';
if strcmp(energy,'CMRO2')
    if strcmp(micro_parameter,'Rsoma')
        reshaped_medians_energy_tot=reshaped_medians_CMRO2_tot;
        reshaped_medians_SANDI_tot=reshaped_medians_rsoma_tot;
    elseif strcmp(micro_parameter,'fsoma')
        reshaped_medians_energy_tot=reshaped_medians_CMRO2_tot;
        reshaped_medians_SANDI_tot=reshaped_medians_fsoma_tot;
    elseif strcmp(micro_parameter,'fsup')
        reshaped_medians_energy_tot=reshaped_medians_CMRO2_tot;
        reshaped_medians_SANDI_tot=reshaped_medians_fsup_tot;
    end
elseif strcmp(energy,'CBF')
    if strcmp(micro_parameter,'Rsoma')
        reshaped_medians_energy_tot=reshaped_medians_CBF_tot;
        reshaped_medians_SANDI_tot=reshaped_medians_rsoma_tot;
    elseif strcmp(micro_parameter,'fsoma')
        reshaped_medians_energy_tot=reshaped_medians_CBF_tot;
        reshaped_medians_SANDI_tot=reshaped_medians_fsoma_tot;
    elseif strcmp(micro_parameter,'fsup')
        reshaped_medians_energy_tot=reshaped_medians_CBF_tot;
        reshaped_medians_SANDI_tot=reshaped_medians_fsup_tot;
    end
end
% 
labels=labels_tot;
z=[];
pvalue=[];
for i = 1:(numel(labels))
    [r,p] = corrcoef(reshaped_medians_energy_tot(:,i), reshaped_medians_SANDI_tot(:,i),'rows','complete');
    z(end+1)=r(2);
    pvalue(end+1)=p(2);
end

% z=[];
% pvalue=[];
% for i = 1:(numel(labels))
%     [r,p] = corrcoef(reshaped_medians_energy(:,i), reshaped_medians_SANDI(:,i),'rows','complete');
%     z(i)=r(2);
%     pvalue(i)=p(2);
% end


% figure, hist(z);
% xlabel('correlation, r');
% ylabel('# regions');

z_accepted=[];
labels_accepted=[];
pvalue_accepted=[];
for i=1:numel(z)
    if pvalue(i)<0.05
        z_accepted(end+1)=z(i);
        labels_accepted(end+1)=labels(i);
        pvalue_accepted(end+1)=pvalue(i);
    end
end

mean_regions_corr=nanmean(z_accepted);%why do we have NaNs?
[h,p,ci,stats]=ttest(atanh(z_accepted));

mean_corr=round(mean_regions_corr,2);
mean_corr=num2str(mean_corr);

figure, hist(z_accepted);
xlabel('correlation, r','FontWeight','bold','FontSize',15);
ylabel('# regions','FontWeight','bold','FontSize',15);
txt=strcat('\mu_r=',mean_corr);
text(0.4,14,txt, 'FontWeight', 'bold','FontSize',15);
grid on 

% z_accepted=z;
% labels_accepted=labels;
if strcmp(energy,'CMRO2')
    if strcmp(micro_parameter,'Rsoma')
        figure, scatter(reshaped_medians_rsoma_tot,reshaped_medians_CMRO2_tot,'filled')
        [r,p]=corrcoef(reshaped_medians_rsoma_tot,reshaped_medians_CMRO2_tot,'rows','complete');
        corr_coef_str=num2str(round(r(2),2));
        txt = {strcat('r = ',corr_coef_str,'***')};
        text(12.5,120,txt,'FontSize',15, 'FontWeight', 'bold','FontSize',12);
        xlabel('Rsoma(\mum)','FontSize',15,'FontWeight','bold');
        ylabel('CMRO_2(\mumol/100g/min)','FontSize',15,'FontWeight','bold');
        grid on

        medians_rsoma_tot_accepted=reshaped_medians_rsoma_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);
        medians_CMRO2_tot_accepted=reshaped_medians_CMRO2_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);

        figure, scatter(medians_rsoma_tot_accepted, medians_CMRO2_tot_accepted,'filled')
        lgd=legend('7','9','30','57','94','95');
        title(lgd,'Labels')
        grid on
        xlabel('Rsoma(\mum)','FontSize',15,'FontWeight','bold');
        ylabel('CMRO_2(\mumol/100g/min)','FontSize',15,'FontWeight','bold');

    elseif strcmp(micro_parameter,'fsoma')
        figure, scatter(reshaped_medians_fsoma_tot,reshaped_medians_CMRO2_tot)

        medians_fsoma_tot_accepted=reshaped_medians_fsoma_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);
        medians_CMRO2_tot_accepted=reshaped_medians_CMRO2_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);

        figure, scatter(medians_fsoma_tot_accepted, medians_CMRO2_tot_accepted)

    elseif strcmp(micro_parameter,'fsup')
        figure, scatter(reshaped_medians_fsup_tot,reshaped_medians_CMRO2_tot)

        medians_fsup_tot_accepted=reshaped_medians_fsup_tot(:,[find(labels==42)]);
        medians_CMRO2_tot_accepted=reshaped_medians_CMRO2_tot(:,[find(labels==42)]);

        figure, scatter(medians_fsup_tot_accepted, medians_CMRO2_tot_accepted)

    end
elseif strcmp(energy,'CBF')
    if strcmp(micro_parameter,'Rsoma')
        figure, scatter(reshaped_medians_rsoma_tot,reshaped_medians_CBF_tot)
        
        medians_rsoma_tot_accepted=reshaped_medians_rsoma_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);
        medians_CBF_tot_accepted=reshaped_medians_CBF_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);

        figure, scatter(medians_rsoma_tot_accepted, medians_CBF_tot_accepted)
    elseif strcmp(micro_parameter,'fsoma')
        figure, scatter(reshaped_medians_fsoma_tot,reshaped_medians_CBF_tot)

        medians_fsoma_tot_accepted=reshaped_medians_fsoma_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);
        medians_CBF_tot_accepted=reshaped_medians_CBF_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);

        figure, scatter(medians_fsoma_tot_accepted, medians_CBF_tot_accepted)
    elseif strcmp(micro_parameter,'fsup')
        figure, scatter(reshaped_medians_fsup_tot,reshaped_medians_CBF_tot)

        medians_fsup_tot_accepted=reshaped_medians_fsup_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);
        medians_CBF_tot_accepted=reshaped_medians_CBF_tot(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);

        figure, scatter(medians_fsup_tot_accepted, medians_CBF_tot_accepted)
    end
end

%%
corr_and_regions = [z_accepted;labels_accepted].';
diff=setxor(labels_accepted,unique(V_atlas_tot));

%% run this section if
%you want to check where pos and neg corr regions are
binary_z=[];
for i=1:length(z_accepted)
    if z_accepted(i)>0
        binary_z(i)=1;
    else
        binary_z(i)=-1;
    end
end
corr_and_regions = [binary_z;labels_accepted].';
%%
%
V_atlas = V_atlas_tot;
for ii = 1:length(V_atlas_tot(:))%lo fa per tutti i valori dell'immagine
    if any(0==V_atlas_tot(ii))
        V_atlas(ii)=0;
    elseif any(diff==V_atlas_tot(ii))
        V_atlas(ii)=0;
    else
        idx=find(corr_and_regions(:,2)==V_atlas_tot(ii));
        corr=corr_and_regions(:,1);
        V_atlas(ii)=corr(idx);
    end
end


%
%plot only regions with significant p-values
fig=figure;
a=28:4:68;%12:4:72;
for i = 1:length(a)
    subplot(4,4,i)
    imagesc(rot90(V_atlas(:,:,a(i))))%[0,8]
    axis equal
    axis off
    caxis([-1,+1])
end
h=axes(fig,'visible','off');
colormap Gray
colorbar(h,'orientation','horizontal','Location','SouthOutside','FontSize',12);
sgtitle('Regional correlation map thresholded for p<0.05')
caxis([-1,+1])

%check for symmetry between positive and negative correlation.


% Save maps
V_corr_map=V_atlas;
if strcmp(energy,'CMRO2')
    if strcmp(micro_parameter,'fsoma')
        V_corr_CMRO2vsfsoma_map=V_corr_map;
        labels_accepted_CMRO2vsfsoma_map=labels_accepted;
        z_accepted_CMRO2vsfsoma_map=z_accepted;
    elseif strcmp(micro_parameter,'Rsoma')
        V_corr_CMRO2vsrsoma_map=V_corr_map;
        labels_accepted_CMRO2vsrsoma_map=labels_accepted;
        z_accepted_CMRO2vsrsoma_map=z_accepted;
    elseif strcmp(micro_parameter,'fsup')
        V_corr_CMRO2vsfsup_map=V_corr_map;
        labels_accepted_CMRO2vsfsup_map=labels_accepted;
        z_accepted_CMRO2vsfsup_map=z_accepted;
    end
elseif strcmp(energy,'CBF')
    if strcmp(micro_parameter,'fsoma')
        V_corr_CBFvsfsoma_map=V_corr_map;
        labels_accepted_CBFvsfsoma_map=labels_accepted;
        z_accepted_CBFvsfsoma_map=z_accepted;
    elseif strcmp(micro_parameter,'Rsoma')
        V_corr_CBFvsrsoma_map=V_corr_map;
        labels_accepted_CBFvsrsoma_map=labels_accepted;
        z_accepted_CBFvsrsoma_map=z_accepted;
    elseif strcmp(micro_parameter,'fsup')
        V_corr_CBFvsfsup_map=V_corr_map;
        labels_accepted_CBFvsfsup_map=labels_accepted;
        z_accepted_CBFvsfsup_map=z_accepted;
    end
end
% V_atlas=V_mean_energy_map;

%% plot spatially mean parametric maps
micro_parameter='Rsoma';

labels=unique(V_atlas_tot);
labels(labels==0)=[];
diff=setxor(labels,unique(V_atlas_tot));
diff=0;

V_atlas = V_atlas_tot;
labels_and_par = [labels';mean_rsoma_tot].';
for ii = 1:length(V_atlas_tot(:))%lo fa per tutti i valori dell'immagine
    if any(0==V_atlas_tot(ii))
        V_atlas(ii)=0;
    elseif any(diff==V_atlas_tot(ii))
        V_atlas(ii)=0;
    else
        idx=find(labels_and_par(:,1)==V_atlas_tot(ii));
        par=labels_and_par(:,2);
        V_atlas(ii)=par(idx);
    end
end

fig=figure;
a=28:4:68;%12:4:72;
for i = 1:length(a)
    subplot(4,4,i)
    imagesc(rot90(V_atlas(:,:,a(i))))%[0,8]
    axis equal
    axis off
    %caxis([-1,+1])
end
h=axes(fig,'visible','off');
%colormap Gray
colorbar(h,'orientation','horizontal','Location','SouthOutside','FontSize',12);
sgtitle('Regional correlation map thresholded for p<0.05')
%caxis([-1,+1])

%V_mean_energy_map=V_atlas;

if strcmp(micro_parameter,'fsoma')
    V_mean_fsoma_map=V_atlas;
else
    V_mean_rsoma_map=V_atlas;
end

%% GM across subjects analysis
medians_CBF = [];
medians_CMRO2 = [];
medians_fsoma = [];
medians_rsoma = [];
medians_fc = [];
medians_fsup = [];

medians_fextra = [];
medians_fneurite = [];
medians_Din = [];
medians_De = [];

V_GM=V_GM_tot;
V_GM(V_GM>0.5)=1;
V_GM(V_GM<1)=0;

for i = 1:1:n_subjs %here we loose information about the real number of subj
    tic
    V_CBF_tot = V_CBF_tots{i};
    V_CMRO2_tot = V_CMRO2_tots{i};
    V_fsoma_tot = V_fsoma_tots{i};
    V_rsoma_tot = V_rsoma_tots{i};
    V_fc_tot = V_fc_tots{i}; 
    V_fsup_tot = V_fsup_tots{i}; 

    V_fextra_tot = V_fextra_tots{i}; 
    V_fneurite_tot = V_fneurite_tots{i}; 
    V_Din_tot = V_Din_tots{i};
    V_De_tot = V_De_tots{i};

    V_fsoma_to_mask=V_fsoma_tot;
    V_fsoma_to_mask(V_fsoma_to_mask>0.15)=1;
    V_fsoma_to_mask(V_fsoma_to_mask<1)=0;

    V_CBF = V_CBF_tot;
    V_CMRO2 = V_CMRO2_tot;
    V_fsoma = V_fsoma_tot;
    V_rsoma = V_rsoma_tot;
    V_fc = V_fc_tot;
    V_fsup = V_fsup_tot;

    V_fextra = V_fextra_tot;
    V_fneurite = V_fneurite_tot;
    V_Din = V_Din_tot;
    V_De = V_De_tot;

    V_CBF_masked = V_CBF.*V_fsoma_to_mask.*V_GM;
    V_CMRO2_masked = V_CMRO2.*V_fsoma_to_mask.*V_GM;
    V_fsoma_masked = V_fsoma.*V_fsoma_to_mask.*V_GM;
    V_rsoma_masked = V_rsoma.*V_fsoma_to_mask.*V_GM;
    V_fc_masked = V_fc.*V_fsoma_to_mask.*V_GM;
    V_fsup_masked = V_fsup.*V_fsoma_to_mask.*V_GM;

    V_fextra_masked = V_fextra.*V_fsoma_to_mask.*V_GM;
    V_fneurite_masked = V_fneurite.*V_fsoma_to_mask.*V_GM;
    V_Din_masked = V_Din.*V_fsoma_to_mask.*V_GM;
    V_De_masked = V_De.*V_fsoma_to_mask.*V_GM;    

    V_CBF_masked_vec = V_CBF_masked(:);
    V_CMRO2_masked_vec = V_CMRO2_masked(:);
    V_fsoma_masked_vec = V_fsoma_masked(:);
    V_rsoma_masked_vec = V_rsoma_masked(:);
    V_fc_masked_vec = V_fc_masked(:);
    V_fsup_masked_vec = V_fsup_masked(:);

    V_fextra_masked_vec = V_fextra_masked(:);
    V_fneurite_masked_vec = V_fneurite_masked(:);
    V_Din_masked_vec = V_Din_masked(:);
    V_De_masked_vec = V_De_masked(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    indices_CBF = [];
    for ii = 1:numel(V_CBF_masked_vec)
        if V_CBF_masked_vec(ii)==0 %|| V_energy_masked(ii)>up_thr
            indices_CBF(end+1)=ii;
        end
    end

    indices_CMRO2 = [];
    for ii = 1:numel(V_CMRO2_masked_vec)
        if V_CMRO2_masked_vec(ii)==0 %|| V_energy_masked(ii)>up_thr
            indices_CMRO2(end+1)=ii;
        end
    end


    disp(strcat('Processing subj',num2str(i)))

    indices = cat(2, indices_CBF, indices_CMRO2);
    indices_unique = unique(indices);


    V_CBF_masked_vec(indices_unique)=[];
    V_CMRO2_masked_vec(indices_unique)=[];
    V_fsoma_masked_vec(indices_unique)=[];
    V_rsoma_masked_vec(indices_unique)=[];
    V_fc_masked_vec(indices_unique)=[];
    V_fsup_masked_vec(indices_unique)=[];

    V_fextra_masked_vec(indices_unique)=[];
    V_fneurite_masked_vec(indices_unique)=[];
    V_Din_masked_vec(indices_unique)=[];
    V_De_masked_vec(indices_unique)=[];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    median_CBF = nanmedian(V_CBF_masked_vec);
    median_CMRO2 = nanmedian(V_CMRO2_masked_vec);
    median_fsoma = nanmedian(V_fsoma_masked_vec);
    median_rsoma = nanmedian(V_rsoma_masked_vec);
    median_fc = nanmedian(V_fc_masked_vec);
    median_fsup = nanmedian(V_fsup_masked_vec);

    median_fextra = nanmedian(V_fextra_masked_vec);
    median_fneurite = nanmedian(V_fneurite_masked_vec);
    median_Din = nanmedian(V_Din_masked_vec);
    median_De = nanmedian(V_De_masked_vec);

    medians_CBF(end+1) = median_CBF;
    medians_CMRO2(end+1) = median_CMRO2;
    medians_fsoma(end+1) = median_fsoma;
    medians_rsoma(end+1) = median_rsoma;
    medians_fc(end+1) = median_fc;
    medians_fsup(end+1) = median_fsup;

    medians_fextra(end+1) = median_fextra;
    medians_fneurite(end+1) = median_fneurite;
    medians_Din(end+1) = median_Din;
    medians_De(end+1) = median_De;

    toc
end

%% Choose the couple of parameters for GM EMI model 

energy_parameter='CMRO2';%CBF
micro_parameter='Rsoma';%other GM parameters: fsoma fc fsup 

%%
if strcmp(energy_parameter,'CMRO2')
    dependent_parameter_label='CMRO_2(\mumol/100g/min)';
    medians_energy=medians_CMRO2;
elseif strcmp(energy_parameter,'CBF')
    dependent_parameter_label='CBF(ml/100g/min)';
    medians_energy=medians_CBF;
end

if strcmp(micro_parameter,'Rsoma')
    micro_parameter_label='Rsoma(\mum)';
    medians_SANDI=medians_rsoma;
elseif strcmp(micro_parameter,'fsoma')
    micro_parameter_label='fsoma';
    medians_SANDI=medians_fsoma;
elseif strcmp(micro_parameter,'fsup')
    micro_parameter_label='fsup(m^{-1})';
    medians_SANDI=medians_fsup;
elseif strcmp(micro_parameter,'fc')
    micro_parameter_label='fc(m^{-3})';
    medians_SANDI=medians_fc;
end

%%
[r,p] = corrcoef(medians_energy, medians_SANDI, 'rows','complete');
corr_coef = round(r(2),2);
p_value = p(2);
corr_coef_str = num2str(corr_coef);
p_value_str = num2str(p_value);




[P,S] = polyfit(medians_SANDI,medians_energy,1);
%yfit = P(1)*medians_SANDI+P(2);

[yfit,delta] = polyconf(P,medians_SANDI,S,'alpha',0.05);

fit_low=cat(1,medians_SANDI,yfit-delta);
fit_low=fit_low';
fit_low_sorted=sortrows(fit_low,1,'descend');

fit_high=cat(1,medians_SANDI,yfit+delta);
fit_high=fit_high';
fit_high_sorted=sortrows(fit_high,1,'descend');



figure, 
s=scatter(medians_SANDI,medians_energy);
hold on
h=plot(fit_low_sorted(:,1),fit_low_sorted(:,2),'r--',fit_low_sorted(:,1),fit_high_sorted(:,2),'r--','LineWidth',1.5);
%s = errorbar(medians_SANDI, medians_energy, SEs_energy, SEs_energy, SEs_SANDI, SEs_SANDI,'o');
%hold on
%plot(medians_SANDI,yfit,'--','LineWidth',3,'Color',"#000000");
xlabel(micro_parameter_label,'FontSize',15,'FontWeight','bold');
ylabel(dependent_parameter_label,'FontSize',15,'FontWeight','bold');
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
txt = {strcat('r = ',corr_coef_str,'')};%,strcat('p-value = ',p_value_str)
if strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'Rsoma')
    text(13.5,60,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fsoma')
    text(0.39,60,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'Rsoma')
    text(13.5,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fsoma')
    text(0.39,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fsup')
    text(8*10^4,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fsup')
    text(8*10^4,60,txt, 'FontWeight', 'bold','FontSize',12);
end
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'box','off')
x0=400;
y0=400;
width=550;
height=450;
set(gcf,'position',[x0,y0,width,height])
grid on



%%
%Test if CBF and CMRO2 distribution
%in the GM is correct (save hist fig) and if
%there are zeros in the brain regions (save into energy_zeros_brain_subjs)

for i = 1:1:n_subjs %1:1:length(lst)
    
    V_CBF_tot = V_CBF_tots{i};
    V_CMRO2_tot = V_CMRO2_tots{i};
    V_fc_tot = V_fc_tots{i};
    V_fsoma_tot = V_fsoma_tots{i};
    V_fneurite_tot = V_fneurite_tots{i};
    V_fextra_tot = V_fextra_tots{i};
    V_Din_tot = V_Din_tots{i};
    V_De_tot = V_De_tots{i};
    V_rsoma_tot = V_rsoma_tots{i};
    V_fsup_tot = V_fsup_tots{i};

    V_CBF_GM=V_CBF.*V_GM;
    V_CMRO2_GM=V_CMRO2.*V_GM;
    V_CBF_GM(GM_zeros)=[];
    V_CMRO2_GM(GM_zeros)=[];
    %fig_CBF=figure('visible','off');
    figure,
    hist(V_CBF_GM);
    saveas(fig_CBF, strcat('CBF_subj',num2str(i),'.png'));
    %fig_CMRO2=figure('visible','off');
    figure,
    hist(V_CMRO2_GM);
    saveas(fig_CMRO2, strcat('CMRO2_subj',num2str(i),'.png'));
    CBF_zeros_brain_subj=numel(find(V_CBF_GM==0));
    CMRO2_zeros_brain_subj=numel(find(V_CMRO2_GM==0));
    CBF_zeros_brain_subjs(i) = CBF_zeros_brain_subj;
    CMRO2_zeros_brain_subj(i) = CMRO2_zeros_brain_subj;
end



