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

%     %old motion correction
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

    %Linear Interpolation
    img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/250418_corrected_for_RIM/CBF02MNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CBF_map_2MNI.nii.gz');
    img_path_CMRO2 = strcat('/media/nas_rete/Vitality/registered/perf/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CMRO2_map.nii.gz'); %_2MNI2mm
    
%     %Nearest Neighbor Interpolation
%     img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/250627_NN_interpolation/CBF22MNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CBF_map_2MNI.nii.gz');
%     img_path_CMRO2 = strcat('/media/nas_rete/Vitality/maps2MNI/250627_NN_interpolation/CMRO22MNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CMRO2_map_2MNI.nii.gz'); %_2MNI2mm

%     %Nearest Neighbor
%     img_path_rsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250627_NN_interpolation/SANDI2MNI/',subj,'_',run,'_SANDI-fit_Rsoma_2MNI.nii.gz');
%     img_path_fsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250627_NN_interpolation/SANDI2MNI/',subj,'_',run,'_SANDI-fit_fsoma_2MNI.nii.gz'); %_2MNI2mm
%     img_path_De = strcat('/media/nas_rete/Vitality/maps2MNI/250627_NN_interpolation/SANDI2MNI/',subj,'_',run,'_SANDI-fit_fextra_2MNI.nii.gz');
%     img_path_Din = strcat('/media/nas_rete/Vitality/maps2MNI/250627_NN_interpolation/SANDI2MNI/',subj,'_',run,'_SANDI-fit_fneurite_2MNI.nii.gz'); %_2MNI2mm
%     img_path_fneurite = strcat('/media/nas_rete/Vitality/maps2MNI/250627_NN_interpolation/SANDI2MNI/',subj,'_',run,'_SANDI-fit_De_2MNI.nii.gz');
%     img_path_fextra = strcat('/media/nas_rete/Vitality/maps2MNI/250627_NN_interpolation/SANDI2MNI/',subj,'_',run,'_SANDI-fit_Din_2MNI.nii.gz'); %_2MNI2mm

    %Linear Interpolation
    img_path_rsoma = strcat('/media/nas_rete/Vitality/registered/SANDI/',subj,'_',run,'_SANDI-fit_Rsoma_2mm.nii.gz');
    img_path_fsoma = strcat('/media/nas_rete/Vitality/registered/SANDI/',subj,'_',run,'_SANDI-fit_fsoma_2mm.nii.gz');
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
        if V_rsoma(i) < 5
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

% first_fc=V_fc_tots{1};
% figure, imagesc(first_fc(:,:,45))
% clim([0,10^(14)]);
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

%
% V_fsup_one=V_fsup_tots{1};
% figure, imagesc(rot90(V_fsup_one(:,:,45)));
% title('Superficial Soma Density map')
% V_fsup_one_array=V_fsup_one(:);
% figure, hist(V_fsup_one_array);
% title('Superficial Soma Density Distribution');
% grid on

%%
%linear interpolation
load('/media/nas_rete/Work_manuela/DWI_En_modeling_results/250625_IVMethod_withCBFandCMRO2separated/load_data.mat');

%% The energy vs microstructural parameters relationship
% in GM is first explored to exclude outliers 

%choose the energy parameter: either CMRO2 or CBF
energy_parameter='CMRO2';

%% GM across subjects analysis
medians_energy=[];

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
    if strcmp(energy_parameter,'CBF')
        V_energy_tot = V_CBF_tots{i};
    elseif strcmp(energy_parameter,'CMRO2')
        V_energy_tot = V_CMRO2_tots{i};
    end

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

    V_energy = V_energy_tot;
    V_fsoma = V_fsoma_tot;
    V_rsoma = V_rsoma_tot;
    V_fc = V_fc_tot;
    V_fsup = V_fsup_tot;

    V_fextra = V_fextra_tot;
    V_fneurite = V_fneurite_tot;
    V_Din = V_Din_tot;
    V_De = V_De_tot;

    V_energy_masked = V_energy.*V_fsoma_to_mask.*V_GM;

    V_fsoma_masked = V_fsoma.*V_fsoma_to_mask.*V_GM;
    V_rsoma_masked = V_rsoma.*V_fsoma_to_mask.*V_GM;
    V_fc_masked = V_fc.*V_fsoma_to_mask.*V_GM;
    V_fsup_masked = V_fsup.*V_fsoma_to_mask.*V_GM;
 
    V_fextra_masked = V_fextra.*V_fsoma_to_mask.*V_GM;
    V_fneurite_masked = V_fneurite.*V_fsoma_to_mask.*V_GM;
    V_Din_masked = V_Din.*V_fsoma_to_mask.*V_GM;
    V_De_masked = V_De.*V_fsoma_to_mask.*V_GM;    

    V_energy_masked_vec = V_energy_masked(:);

    V_fsoma_masked_vec = V_fsoma_masked(:);
    V_rsoma_masked_vec = V_rsoma_masked(:);
    V_fc_masked_vec = V_fc_masked(:);
    V_fsup_masked_vec = V_fsup_masked(:);
   
    V_fextra_masked_vec = V_fextra_masked(:);
    V_fneurite_masked_vec = V_fneurite_masked(:);
    V_Din_masked_vec = V_Din_masked(:);
    V_De_masked_vec = V_De_masked(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    indices_energy = [];
    indices_energy=find(V_energy_masked_vec<=0);

    disp(strcat('Processing subj',num2str(i)))

    V_energy_masked_vec(indices_energy)=[];
    V_fsoma_masked_vec(indices_energy)=[];
    V_rsoma_masked_vec(indices_energy)=[];
    V_fc_masked_vec(indices_energy)=[];
    V_fsup_masked_vec(indices_energy)=[];

    V_fextra_masked_vec(indices_energy)=[];
    V_fneurite_masked_vec(indices_energy)=[];
    V_Din_masked_vec(indices_energy)=[];
    V_De_masked_vec(indices_energy)=[];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    median_energy = nanmedian(V_energy_masked_vec);
    median_fsoma = nanmedian(V_fsoma_masked_vec);
    median_rsoma = nanmedian(V_rsoma_masked_vec);
    median_fc = nanmedian(V_fc_masked_vec);
    median_fsup = nanmedian(V_fsup_masked_vec);

    median_fextra = nanmedian(V_fextra_masked_vec);
    median_fneurite = nanmedian(V_fneurite_masked_vec);
    median_Din = nanmedian(V_Din_masked_vec);
    median_De = nanmedian(V_De_masked_vec);

    medians_energy(end+1) = median_energy;
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

%% 
%if you want to run without for loop
%micro_parameter='Rsoma';%other GM parameters: fsoma fsup 

%% detect outlier
% note: you can run this section after the plot section

outlier_idx=find(medians_energy>120);%200 FOR CMRO2


medians_energy(outlier_idx)=[];

medians_fsoma(outlier_idx)=[];
medians_fsup(outlier_idx)=[];
medians_rsoma(outlier_idx)=[];
medians_fc(outlier_idx)=[];
disp('outlier successfully removed')

%% plot

micro_parameters={'Rsoma' 'fsoma' 'fsup' 'fc'};
for i=1:numel(micro_parameters)
    micro_parameter=micro_parameters{i};


    
    if strcmp(energy_parameter,'CMRO2')
        dependent_parameter_label='CMRO_2(\mumol/100g/min)';
    elseif strcmp(energy_parameter,'CBF')
        dependent_parameter_label='CBF(ml/100g/min)';
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
    
    [r,p] = corrcoef(medians_energy, medians_SANDI, 'rows','complete');
    corr_coef = round(r(2),2);
    p_value = p(2);
    corr_coef_str = num2str(corr_coef);
    p_value_str = num2str(p_value);
%
    [P,S] = polyfit(medians_SANDI,medians_energy,1);
    %yfit = P(1)*medians_SANDI+P(2);
    %P polynomial coefficients of grade n=1;
    %S structure which contains information needed to estimate 
    % the error through polyconf.


    [yfit,delta] = polyconf(P,medians_SANDI,S,'alpha',0.05);

    fit_low=cat(1,medians_SANDI,yfit-delta);
    fit_low=fit_low';
    fit_low_sorted=sortrows(fit_low,1,'descend');

    fit_high=cat(1,medians_SANDI,yfit+delta);
    fit_high=fit_high';
    fit_high_sorted=sortrows(fit_high,1,'descend');
    
    array_low=yfit-delta;
    low=sort(array_low,'descend');
%
    figure,
    s=scatter(medians_SANDI,medians_energy);
%     hold on
%     l=plot(medians_SANDI,low,'r-');
    hold on
    %fit_low_sorted(:,1) is equal to fit_high_sorted(:,1)
    h=plot(fit_low_sorted(:,1),fit_low_sorted(:,2),'r--',fit_high_sorted(:,1),fit_high_sorted(:,2),'r--','LineWidth',1.5);
    xlabel(micro_parameter_label,'FontSize',15,'FontWeight','bold');
    ylabel(dependent_parameter_label,'FontSize',15,'FontWeight','bold');
    s.LineWidth = 0.6;
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = [0 0.5 0.5];

    if p(2)<0.05 && p(2)>0.01
        txt = {strcat('r = ',corr_coef_str,'*')};
    elseif p(2)<0.01 && p(2)>0.001
        txt = {strcat('r = ',corr_coef_str,'**')};
    elseif p(2)<0.001
        txt = {strcat('r = ',corr_coef_str,'***')};
    elseif p(2)>0.05
        txt = {strcat('r = ',corr_coef_str,'')};
    end

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
    elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fc')
        text(3.55*10^(13),60,txt, 'FontWeight', 'bold','FontSize',12);
    elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fc')
        text(3.55*10^(13),120,txt, 'FontWeight', 'bold','FontSize',12);
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

end





%% Correlation across regions
GM_pve_threshold=0.5;
fsoma_threshold=0.15;

%set GM partial volume estimate map
V_GM = V_GM_tot;
V_GM(V_GM>GM_pve_threshold)=1;
V_GM(V_GM<1)=0;

regions = unique(V_atlas_tot(:));
%background removal
regions(1)=[];
n_regions=numel(regions);

%% Choose energy parameter to be considered
energy_parameter='CMRO2';

%%


medians_energy_subj=[];
medians_fc_subj= [];
medians_fextra_subj= [];
medians_fneurite_subj = [];
medians_fsoma_subj= [];
medians_Din_subj = [];
medians_De_subj =[];
medians_rsoma_subj = [];
medians_fsup_subj = [];

medians_energy_subjs=[];
medians_De_subjs = [];
medians_rsoma_subjs = [];
medians_fsoma_subjs = [];
medians_fneurite_subjs = [];
medians_fextra_subjs = [];
medians_Din_subjs = [];
medians_fc_subjs = [];
medians_fsup_subjs = [];

n_voxels_lst = [];
n_voxels_lst_allsubjs=[];

corr_for_each_subj_energy_rsoma=[];
pvalue_for_each_subj_energy_rsoma=[];
corr_for_each_subj_energy_fsoma=[];
pvalue_for_each_subj_energy_fsoma=[];
corr_for_each_subj_energy_fsup=[];
pvalue_for_each_subj_energy_fsup=[];
corr_for_each_subj_energy_fc=[];
pvalue_for_each_subj_energy_fc=[];
corr_for_each_subj_energy_fneurite=[];
pvalue_for_each_subj_energy_fneurite=[];

start_time=tic;
for i = 1:1:n_subjs %1:1:length(lst)

    if strcmp(energy_parameter,'CBF')
        V_energy_tot = V_CBF_tots{i};
    elseif strcmp(energy_parameter,'CMRO2')
        V_energy_tot = V_CMRO2_tots{i};
    end
    V_rsoma_tot = V_rsoma_tots{i};
    V_fsup_tot = V_fsup_tots{i};
    V_fc_tot = V_fc_tots{i};
    V_fsoma_tot = V_fsoma_tots{i};
    V_fneurite_tot = V_fneurite_tots{i};
    V_fextra_tot = V_fextra_tots{i};
    V_Din_tot = V_Din_tots{i};
    V_De_tot = V_De_tots{i};
   
    medians_energy_subj=[];
    medians_rsoma_subj=[];
    medians_fsup_subj=[];
    medians_fsoma_subj=[];
    medians_fc_subj=[];
    medians_fneurite_subj=[];
    medians_fextra_subj=[];
    medians_Din_subj=[];
    medians_De_subj=[];

    n_voxels_lst=[];

    for k = 1:n_regions
        tic

        V_atlas = V_atlas_tot;

        for ii = 1:length(V_atlas_tot(:))
            if V_atlas_tot(ii)==regions(k)
                V_atlas(ii)=1;
            else
                V_atlas(ii)=0;
            end
        end
      
        V_energy = V_energy_tot;
        V_fc = V_fc_tot;
        V_fextra = V_fextra_tot;
        V_fneurite = V_fneurite_tot;
        V_rsoma = V_rsoma_tot;
        V_fsoma = V_fsoma_tot;
        V_Din = V_Din_tot;
        V_De = V_De_tot;
        V_fsup = V_fsup_tot;
        V_fsoma_to_mask=V_fsoma_tot;

        V_fsoma_to_mask(V_fsoma_to_mask>fsoma_threshold)=1;
        V_fsoma_to_mask(V_fsoma_to_mask<1)=0;

        V_mask = V_atlas.*V_GM;%.*V_fsoma_to_mask


        
        %%check
        %V_energy_withoutzerosatlas=V_energy;
        %V_energy_withoutzerosatlas(mask_zeros)=[];
        %figure,hist(V_energy_withoutzerosatlas(:));

        V_energy_masked = V_energy.*V_mask;
        V_fc_masked = V_fc.*V_mask;
        V_fextra_masked = V_fextra.*V_mask;
        V_fneurite_masked = V_fneurite.*V_mask;
        V_rsoma_masked = V_rsoma.*V_mask;
        V_fsoma_masked = V_fsoma.*V_mask;
        V_Din_masked = V_Din.*V_mask;
        V_De_masked = V_De.*V_mask;
        V_fsup_masked = V_fsup.*V_mask;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        disp(strcat('region', num2str(k), 'subject', num2str(i)))
        
        %find zeros of the mask (so outside brain)
        mask_zeros=find(V_mask==0);

        V_energy_masked(mask_zeros)=[];
        V_fc_masked(mask_zeros)=[];
        V_fsoma_masked(mask_zeros)=[];
        V_rsoma_masked(mask_zeros)=[];
        V_fsup_masked(mask_zeros)=[];
        V_fneurite_masked(mask_zeros)=[];
        V_fextra_masked(mask_zeros)=[];
        V_Din_masked(mask_zeros)=[];
        V_De_masked(mask_zeros)=[];
        
%         %check distribution here
%         v_energy_masked_withoutremovinginnerzeros=v_energy_masked;     
%         figure, hist(v_energy_masked_withoutremovinginnerzeros(:));
%       
        %find zeros inside brain
        indices_energy = [];
        indices_energy=find(V_energy_masked<=0);
        
        disp(strcat('region', num2str(k)))

        V_energy_masked(indices_energy)=[];
        %remove the corresponding voxels in microstructural maps
        %otherwise you can have imbalance data problem
        V_fc_masked(indices_energy)=[];
        V_fsoma_masked(indices_energy)=[];
        V_rsoma_masked(indices_energy)=[];
        V_fsup_masked(indices_energy)=[];
        V_fneurite_masked(indices_energy)=[];
        V_fextra_masked(indices_energy)=[];
        V_Din_masked(indices_energy)=[];
        V_De_masked(indices_energy)=[];

%        %check distribution here
%        v_energy_masked_afterremovinginnerzeros=v_energy_masked;     
%        figure, hist(v_energy_masked_afterremovinginnerzeros(:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        n_voxels = numel(V_energy_masked); 
        %n voxels is the same for all paramaters
        %and different for each subj
        n_voxels_lst(end+1) = n_voxels;

        median_energy = median(V_energy_masked,'omitnan');
        median_fc = median(V_fc_masked,'omitnan');
        median_fextra = median(V_fextra_masked,'omitnan');
        median_fneurite = median(V_fneurite_masked,'omitnan');
        median_fsoma = median(V_fsoma_masked,'omitnan');
        median_rsoma = median(V_rsoma_masked,'omitnan');
        median_Din = median(V_Din_masked,'omitnan');
        median_De = median(V_De_masked,'omitnan');
        median_fsup=median(V_fsup_masked,'omitnan');

        medians_energy_subj(end+1) = median_energy;
        medians_rsoma_subj(end+1) = median_rsoma;
        medians_fsoma_subj(end+1) = median_fsoma;
        medians_fsup_subj(end+1) = median_fsup;
        medians_fc_subj(end+1) = median_fc;
        medians_fextra_subj(end+1) = median_fextra;
        medians_fneurite_subj(end+1) = median_fneurite;
        medians_Din_subj(end+1) = median_Din;
        medians_De_subj(end+1) = median_De;
      
    toc
    end

    n_voxels_lst_allsubjs(i,:) = n_voxels_lst;

    medians_energy_subjs(i,:) = medians_energy_subj;
    medians_fc_subjs(i,:) = medians_fc_subj;
    medians_fsoma_subjs(i,:) = medians_fsoma_subj;
    medians_rsoma_subjs(i,:) = medians_rsoma_subj;
    medians_fextra_subjs(i,:) = medians_fextra_subj;
    medians_fneurite_subjs(i,:) = medians_fneurite_subj;
    medians_Din_subjs(i,:) = medians_Din_subj;
    medians_De_subjs(i,:) = medians_De_subj;
    medians_fsup_subjs(i,:) = medians_fsup_subj;

    [r_energy_rsoma,p_energy_rsoma] = corrcoef(medians_energy_subj, medians_rsoma_subj, 'rows','complete');
    [r_energy_fsoma,p_energy_fsoma] = corrcoef(medians_energy_subj, medians_fsoma_subj, 'rows','complete');
    [r_energy_fsup,p_energy_fsup] = corrcoef(medians_energy_subj, medians_fsup_subj, 'rows','complete');
    [r_energy_fc,p_energy_fc] = corrcoef(medians_energy_subj, medians_fc_subj, 'rows','complete');
    [r_energy_fneurite,p_energy_fneurite] = corrcoef(medians_energy_subj, medians_fneurite_subj, 'rows','complete');

    corr_for_each_subj_energy_rsoma(i)=r_energy_rsoma(2);
    pvalue_for_each_subj_energy_rsoma(i)=p_energy_rsoma(2);

    corr_for_each_subj_energy_fsoma(i)=r_energy_fsoma(2);
    pvalue_for_each_subj_energy_fsoma(i)=p_energy_fsoma(2);

    corr_for_each_subj_energy_fsup(i)=r_energy_fsup(2);
    pvalue_for_each_subj_energy_fsup(i)=p_energy_fsup(2);

    corr_for_each_subj_energy_fc(i)=r_energy_fc(2);
    pvalue_for_each_subj_energy_fc(i)=p_energy_fc(2);

    corr_for_each_subj_energy_fneurite(i)=r_energy_fneurite(2);
    pvalue_for_each_subj_energy_fneurite(i)=p_energy_fneurite(2);


    disp(strcat('Finished subject', num2str(i),'Starting subject', num2str(i+1)))

end
timeElapsed = toc(start_time);

disp(strcat('Total computing time =',num2str(round(timeElapsed/60,2)),'min'));

%% select parameters to investigate
micro_parameter='fneurite'; 

%% remove outlier detected in GM corr across subjs section
outlier_idx=26;
corr_for_each_subj_energy_rsoma(outlier_idx)=[];
corr_for_each_subj_energy_fsoma(outlier_idx)=[];
corr_for_each_subj_energy_fsup(outlier_idx)=[];
corr_for_each_subj_energy_fc(outlier_idx)=[];

%% for each subj 

if strcmp(micro_parameter,'Rsoma') 
    corr_for_each_subj=corr_for_each_subj_energy_rsoma;
elseif strcmp(micro_parameter,'fsoma') 
    corr_for_each_subj=corr_for_each_subj_energy_fsoma;
elseif strcmp(micro_parameter,'fsup') 
    corr_for_each_subj=corr_for_each_subj_energy_fsup;
elseif strcmp(micro_parameter,'fc')
    corr_for_each_subj=corr_for_each_subj_energy_fc;
elseif strcmp(micro_parameter,'fneurite')
    corr_for_each_subj=corr_for_each_subj_energy_fneurite;
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

%% remove outlier detected in GM section
outlier_idx=26;
medians_energy_subjs(outlier_idx,:)=[];
medians_rsoma_subjs(outlier_idx,:)=[];
medians_fsoma_subjs(outlier_idx,:)=[];
medians_fsup_subjs(outlier_idx,:)=[];
medians_fc_subjs(outlier_idx,:)=[];

medians_fextra_subjs(outlier_idx,:)=[];
medians_fneurite_subjs(outlier_idx,:)=[];
medians_Din_subjs(outlier_idx,:)=[];
medians_De_subjs(outlier_idx,:)=[];

disp('Outlier successfully removed');

%% average across subjects

mean_energy_tot = nanmean(medians_energy_subjs,1);
mean_fc_tot = nanmean(medians_fc_subjs,1);
mean_fextra_tot = nanmean(medians_fextra_subjs,1);
mean_fneurite_tot = nanmean(medians_fneurite_subjs,1);
mean_fsoma_tot = nanmean(medians_fsoma_subjs,1);
mean_rsoma_tot = nanmean(medians_rsoma_subjs,1);
mean_Din_tot = nanmean(medians_Din_subjs,1);
mean_De_tot = nanmean(medians_De_subjs,1);
mean_fsup_tot = nanmean(medians_fsup_subjs,1);

mean_energy = mean_energy_tot;
mean_fc = mean_fc_tot;
mean_fextra = mean_fextra_tot;
mean_fneurite = mean_fneurite_tot;
mean_fsoma = mean_fsoma_tot;
mean_rsoma = mean_rsoma_tot;
mean_Din = mean_Din_tot;
mean_De = mean_De_tot;
mean_fsup = mean_fsup_tot;

SE_energy = nanstd(medians_energy_subjs,0,1)/sqrt(n_subjs);%0 flag -> normalization by n.
SE_fc = nanstd(medians_fc_subjs,0,1)/sqrt(n_subjs);
SE_fextra = nanstd(medians_fextra_subjs,0,1)/sqrt(n_subjs);
SE_fneurite = nanstd(medians_fneurite_subjs,0,1)/sqrt(n_subjs);
SE_fsoma = nanstd(medians_fsoma_subjs,0,1)/sqrt(n_subjs);
SE_rsoma = nanstd(medians_rsoma_subjs,0,1)/sqrt(n_subjs);
SE_Din = nanstd(medians_Din_subjs,0,1)/sqrt(n_subjs);
SE_De = nanstd(medians_De_subjs,0,1)/sqrt(n_subjs);
SE_fsup = nanstd(medians_fsup_subjs,0,1)/sqrt(n_subjs);

% in order to put a threshold based on coefficient of variation

var_rsoma=SE_rsoma./mean_rsoma;
var_fsoma=SE_fsoma./mean_fsoma;
var_fc=SE_fc./mean_fc;
var_fsup=SE_fsup./mean_fsup;

var_fneurite=SE_fneurite./mean_fneurite;
var_fextra=SE_fextra./mean_fextra;
var_Din=SE_Din./mean_Din;
var_De=SE_De./mean_De;

mean_n_voxels_tot = mean(n_voxels_lst_allsubjs,1);

% 
% figure, bar(regions,mean_n_voxels_tot);
% xlabel('regions','FontSize',15);
% ylabel('mean number of voxels','FontSize',15);

idx_low_n_voxels=[];
for i=1:n_regions
    if mean_n_voxels_tot(i)<prctile(mean_n_voxels_tot,40)%mean(mean_n_voxels)/5 %CHECK
        idx_low_n_voxels(end+1)=i;
    end
end

%%%remove idx with low n voxels
labels_tot=regions;
labels=labels_tot;
labels(idx_low_n_voxels)=[];

mean_n_voxels=mean_n_voxels_tot;
mean_n_voxels(idx_low_n_voxels)=[];

mean_energy(idx_low_n_voxels)=[];
mean_fc(idx_low_n_voxels)=[];
mean_fextra(idx_low_n_voxels)=[];
mean_fneurite(idx_low_n_voxels)=[];
mean_fsoma(idx_low_n_voxels)=[];
mean_rsoma(idx_low_n_voxels)=[];
mean_Din(idx_low_n_voxels)=[];
mean_De(idx_low_n_voxels)=[];
mean_fsup(idx_low_n_voxels)=[];

%this array is needed for the MLR analysis
mean_energy_high_n_voxels=mean_energy;
mean_rsoma_high_n_voxels=mean_rsoma;
mean_fsoma_high_n_voxels=mean_fsoma;
mean_fc_high_n_voxels=mean_fc;
mean_fsup_high_n_voxels=mean_fsup;
mean_fextra_high_n_voxels=mean_fextra;
mean_fneurite_high_n_voxels=mean_fneurite;
mean_Din_high_n_voxels=mean_Din;
mean_De_high_n_voxels=mean_De;

%this is to plot errors on both x and y values
SE_energy(idx_low_n_voxels)=[];
SE_rsoma(idx_low_n_voxels)=[];
SE_fsoma(idx_low_n_voxels)=[];
SE_fc(idx_low_n_voxels)=[];
SE_fsup(idx_low_n_voxels)=[];

SE_fneurite(idx_low_n_voxels)=[];
SE_fextra(idx_low_n_voxels)=[];
SE_Din(idx_low_n_voxels)=[];
SE_De(idx_low_n_voxels)=[];




%this is to calculate coefficient of variation 
% of microstructural parameter 
var_rsoma(idx_low_n_voxels)=[];
var_fsoma(idx_low_n_voxels)=[];
var_fc(idx_low_n_voxels)=[];
var_fsup(idx_low_n_voxels)=[];

var_fneurite(idx_low_n_voxels)=[];
var_fextra(idx_low_n_voxels)=[];
var_Din(idx_low_n_voxels)=[];
var_De(idx_low_n_voxels)=[];

disp('First thresholding step successfully done');

%% Find indices corresponding to high coefficient of variation of microstructural parameters

idx_high_var_rsoma=find(var_rsoma>prctile(var_rsoma,75));
idx_high_var_fsoma=find(var_fsoma>prctile(var_fsoma,75));
idx_high_var_fsup=find(var_fsup>prctile(var_fsup,75));
idx_high_var_fc=find(var_fc>prctile(var_fc,75));
idx_high_var_fneurite=find(var_fneurite>prctile(var_fneurite,75));

%if based on SE instead of var
%idx_high_SE=find(SE_rsoma_tot>prctile(SE_rsoma_tot,95));

%% in order to do Energy vs one Microparameter analysis (run if corr analysis)

if strcmp(micro_parameter,'Rsoma')
    mean_energy(idx_high_var_rsoma)=[];
    SE_energy(idx_high_var_rsoma)=[];
    mean_rsoma(idx_high_var_rsoma) = [];
    SE_rsoma(idx_high_var_rsoma)=[];   
elseif strcmp(micro_parameter,'fsoma')
    mean_energy(idx_high_var_fsoma)=[];
    SE_energy(idx_high_var_fsoma)=[];
    mean_fsoma(idx_high_var_fsoma)=[];
    SE_fsoma(idx_high_var_fsoma)=[];
elseif strcmp(micro_parameter,'fsup')
    mean_energy(idx_high_var_fsup)=[];
    SE_energy(idx_high_var_fsup)=[]; 
    mean_fsup(idx_high_var_fsup)=[];
    SE_fsup(idx_high_var_fsup)=[]; 
elseif strcmp(micro_parameter,'fc')
    mean_energy(idx_high_var_fc)=[];
    SE_energy(idx_high_var_fc)=[];
    mean_fc(idx_high_var_fc)=[];
    SE_fc(idx_high_var_fc)=[];
elseif strcmp(micro_parameter,'fneurite')
    mean_energy(idx_high_var_fneurite)=[];
    SE_energy(idx_high_var_fneurite)=[];
    mean_fneurite(idx_high_var_fneurite)=[];
    SE_fneurite(idx_high_var_fneurite)=[];
end

disp('Second thresholding step successfully done');

%% plot energy vs microparameter 
micro_parameter='Rsoma';
if strcmp(micro_parameter,'Rsoma')
    mean_SANDI=mean_rsoma;
    SE_SANDI=SE_rsoma;
    unit_of_measure='(\mum)';
elseif strcmp(micro_parameter,'fsoma')
    mean_SANDI=mean_fsoma;
    SE_SANDI=SE_fsoma;
    unit_of_measure='';
elseif strcmp(micro_parameter,'fsup')
    mean_SANDI=mean_fsup;
    SE_SANDI=SE_fsup;
    unit_of_measure='(m^{-1})';
elseif strcmp(micro_parameter,'fc')
    mean_SANDI=mean_fc;
    SE_SANDI=SE_fc;
    unit_of_measure='(m^{-3})';
elseif strcmp(micro_parameter,'fneurite')
    mean_SANDI=mean_fneurite;
    SE_SANDI=SE_fneurite;
    unit_of_measure='';
end

if strcmp(energy_parameter,'CBF')   
    dependent_parameter='CBF (ml/100g/min)';
elseif strcmp(energy_parameter,'CMRO2')
    dependent_parameter='CMRO_2 (\mu mol/100g/min)';
end

[r,p] = corrcoef(mean_energy, mean_SANDI, 'rows','complete');

corr_coef = round(r(2),2);
p_value = p(2);
corr_coef_str = num2str(corr_coef);
p_value_str = num2str(p_value);

P = polyfit(mean_SANDI,mean_energy,1);
yfit_energy = P(1)*mean_SANDI+P(2);

figure, 
s = errorbar(mean_SANDI, mean_energy, SE_energy, SE_energy, SE_SANDI, SE_SANDI,'o');
hold on
plot(mean_SANDI,yfit_energy,'--','LineWidth',3,'Color',"#000000");
xlabel(strcat(micro_parameter,unit_of_measure),'FontSize',15,'FontWeight','bold');
ylabel(dependent_parameter,'FontSize',15,'FontWeight','bold');
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
%ylim([40,180]);
%n_subjs_str=num2str(n_subjs);
if p(2)<0.05 && p(2)>0.01    
    txt = {strcat('r = ',corr_coef_str,'*')};
elseif p(2)<0.01 && p(2)>0.001   
    txt = {strcat('r = ',corr_coef_str,'**')};
elseif p(2)<0.001
    txt = {strcat('r = ',corr_coef_str,'***')};
elseif p(2)>0.05
        txt = {strcat('r = ',corr_coef_str,'')};
end
% text(60,12,txt,'FontWeight', 'Bold');
%annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold','FontSize',15);
if strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'Rsoma')
    text(13.8,60,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fsoma')
    text(0.39,60,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'Rsoma')
    text(13.8,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fsoma')
    text(0.39,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fc')
    text(3*10^(13),140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fsup')
    text(8*10^4,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CMRO2') && strcmp(micro_parameter,'fneurite')
    text(0.3,140,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fsup')
    text(8*10^4,60,txt, 'FontWeight', 'bold','FontSize',12);
elseif strcmp(energy_parameter,'CBF') && strcmp(micro_parameter,'fc')
    text(3*10^13,60,txt, 'FontWeight', 'bold','FontSize',12);

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

%% in order to do Energy vs many Microparameters analysis (run if GLM analysis)
%it must be so in order to run this analysis

%unique is needed, otherwise you delete two times the voxel corresponding
%to the same index and in the second time you will delete
%a different one.
idx_tot=unique([idx_high_var_rsoma,idx_high_var_fsoma,idx_high_var_fsup]);

mean_energy_high_n_voxels(idx_tot)=[];
mean_rsoma_high_n_voxels(idx_tot)=[];
mean_fsoma_high_n_voxels(idx_tot)=[];
mean_fsup_high_n_voxels(idx_tot)=[];
mean_fc_high_n_voxels(idx_tot)=[];

%% Matrix for correlation across subjs analysis (region by region)
%for this analysis you need to use the matrix n_subjsxn_regions

V_atlas=V_atlas_tot;
labels_tot=unique(V_atlas);
labels_tot(labels_tot==0)=[];

medians_energy_thresholded_subjs=medians_energy_subjs;
medians_fsoma_thresholded_subjs=medians_fsoma_subjs;
medians_fsup_thresholded_subjs=medians_fsup_subjs;
medians_rsoma_thresholded_subjs=medians_rsoma_subjs;

medians_energy_thresholded_subjs(:,idx_low_n_voxels)=[];
labels_tot(idx_low_n_voxels)=[];

if strcmp(micro_parameter,'Rsoma')

    medians_rsoma_thresholded_subjs(:,idx_low_n_voxels)=[];
    
    %remove indices corresponding to high var
    labels_tot(idx_high_var_rsoma)=[];

    medians_energy_thresholded_subjs(:,idx_high_var_rsoma)=[];
    medians_rsoma_thresholded_subjs(:,idx_high_var_rsoma)=[];
elseif strcmp(micro_parameter,'fsoma')

    medians_fsoma_thresholded_subjs(:,idx_low_n_voxels)=[];

    labels_tot(idx_high_var_fsoma)=[];
    medians_energy_thresholded_subjs(:,idx_high_var_fsoma)=[];
    medians_fsoma_thresholded_subjs(:,idx_high_var_fsoma)=[];

elseif strcmp(micro_parameter,'fsup')

    medians_fsup_thresholded_subjs(:,idx_low_n_voxels)=[];

    labels_tot(idx_high_var_fsup)=[];
    medians_energy_thresholded_subjs(:,idx_high_var_fsup)=[];
    medians_fsup_thresholded_subjs(:,idx_high_var_fsup)=[];

end

%% check correlation distribution 

if strcmp(micro_parameter,'Rsoma')
    medians_SANDI_thresholded_subjs = medians_rsoma_thresholded_subjs;
elseif strcmp(micro_parameter,'fsoma')
    medians_SANDI_thresholded_subjs = medians_fsoma_thresholded_subjs;
elseif strcmp(micro_parameter,'fsup')
    medians_SANDI_thresholded_subjs = medians_fsup_thresholded_subjs;
end



% 
z=[];
pvalue=[];
for i = 1:n_regions
    [r,p] = corrcoef(medians_energy_thresholded_subjs(:,i), medians_SANDI_thresholded_subjs(:,i),'rows','complete');
    z(end+1)=r(2);
    pvalue(end+1)=p(2);
end

% distribution without thresholding
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

%% if you want to plot regional medians for each subject


  
figure, scatter(medians_SANDI_thresholded_subjs,medians_energy_thresholded_subjs,'filled')
[r,p]=corrcoef(medians_SANDI_thresholded_subjs,medians_energy_thresholded_subjs,'rows','complete');
corr_coef_str=num2str(round(r(2),2));
txt = {strcat('r = ',corr_coef_str,'***')};
text(12.5,120,txt,'FontSize',15, 'FontWeight', 'bold','FontSize',12);
xlabel('Rsoma(\mum)','FontSize',15,'FontWeight','bold');
ylabel('CMRO_2(\mumol/100g/min)','FontSize',15,'FontWeight','bold');
grid on

%plot only samples of regions significantly correlated across subjs
medians_rsoma_tot_accepted=medians_SANDI_thresholded_subjs(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);
medians_energy_tot_accepted=medians_energy_thresholded_subjs(:,[find(labels==7),find(labels==9),find(labels==30),find(labels==57),find(labels==94),find(labels==95)]);

figure, scatter(medians_rsoma_tot_accepted, medians_energy_tot_accepted,'filled')
lgd=legend('7','9','30','57','94','95');
title(lgd,'Labels')
grid on
xlabel(micro_parameter,'FontSize',15,'FontWeight','bold');
ylabel(dependent_parameter,'FontSize',15,'FontWeight','bold');

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
        %substitute the corresponding corr coefficient
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
if strcmp(energy_parameter,'CMRO2')
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
elseif strcmp(energy_parameter,'CBF')
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



