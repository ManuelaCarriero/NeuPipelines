%Loading SANDI map to compare with vascular and metabolic map

%SISTEMA TUTTO IN FUNZIONI

%Each code section contains:
%voxel-wise analysis (one subject):
 %scatter plot 
 %corr curve
 %map varying PVE threashold
 %cbf vs Rsoma (masked with fsoma) 
 %correlation between CBF and fsoma
%voxel-wise analysis (all subjects)

%glm (remove unphysical rsoma and cbf values)


%regional correlation cbf vs rsoma (all subjects)

%clustering
%number of cells density maps calculation
%regional analysis with numerical cellular density
%glm regional analysis

%fsoma vs GM





%Load GM 
img_path_GM = '/storage/shared/Atlas/atlas_GM_on_MNI152_T1_2mm.nii.gz'; %GM resampled on T1 2mm
%atlas_GM resampled on SANDI and CBF dimensions
Vhdr = spm_vol(img_path_GM);
V_GM_tot = spm_read_vols(Vhdr); 

%Loading atlas
img_path_atlas='/storage/shared/Atlas/AAL3v1_2mm.nii.gz';
Vhdr = spm_vol(img_path_atlas);
V_atlas_tot = spm_read_vols(Vhdr);




%% voxel-wise analysis (one subject). PVE simulation. 1

%Load SANDI map
SANDI_parameter='Rsoma';
SANDI_unit_of_measure="(\mu m)";
img_path_SANDI=strcat('/media/nas_rete/Vitality/maps2MNI/sub-001_run-01_SANDI-fit_',SANDI_parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path_SANDI);
V_SANDI_tot = spm_read_vols(Vhdr);



%Loading CBF map
img_path_CBF = '/media/nas_rete/Vitality/maps2MNI/CBF0toMNI/sub-001_task-bh_run-01_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz';
Vhdr = spm_vol(img_path_CBF);
V_CBF_tot = spm_read_vols(Vhdr);
%figure, imagesc(V_CBF(:,:,70))
V_CBF_tot=max(V_CBF_tot,0);%replace negative numbers with zeros
%figure, imagesc(V_CBF_tot(:,:,45))
cbf_unit_of_measure = " (ml/100g/min)";



% min_rsoma=10;
% min_cbf=20;

% min_rsoma_str=num2str(min_rsoma);
% min_cbf_str=num2str(min_cbf);

% path=strcat('/media/nas_rete/Work_manuela/DWI_En_modeling/main/min_rsoma',min_rsoma_str,'andCBF',min_cbf_str);
% if exist([path])==0
% mkdir ([path])
% end

threasholds = 0.1:0.1:1;
index = 1:1:numel(threasholds);

plot_subj_corr(threasholds, index,  V_SANDI_tot, V_CBF_tot, V_GM_tot,SANDI_parameter,SANDI_unit_of_measure,cbf_unit_of_measure);



%% 2
%%%%%fneurite vs CBF (voxel wise considering medians over many subjects and
%%%%%mask with GM and fsoma)

n_subjs=12;

V_CBF_tots={};

for i = 1:1:n_subjs%lst
    if i < 10
        i=num2str(i);
        img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBF0toMNI/sub-00',i,'_task-bh_run-01_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');
    else
        i=num2str(i);
        img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBF0toMNI/sub-0',i,'_task-bh_run-01_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');
    end
    Vhdr = spm_vol(img_path_CBF);
    V_CBF_tot = spm_read_vols(Vhdr);
%     figure, imagesc(rot90(squeeze(V_CBF_tot(30,:,:)))) 
%     imshow(rot90(squeeze(V(30,:,:))),[])
    V_CBF_tots{end+1} = V_CBF_tot;
end



V_fsoma_tots={};

V_fneurite_tots={};

V_fextra_tots={};

V_Din_tots={};

V_De_tots={};

%parameter = 'De';



V_rsoma_tots={};

parameters = {'Rsoma', 'fsoma','fneurite','Din','De','fextra'};

for j = 1:length(parameters)
    parameter = parameters{j};
    for i = 1:1:n_subjs%lst
        if i < 10

            i=num2str(i);
            img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/sub-00',i,'_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
        else

            i=num2str(i);
            img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/sub-0',i,'_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
        end
        Vhdr = spm_vol(img_path_SANDI);
        V_SANDI_tot = spm_read_vols(Vhdr);

        if strcmp(parameter,'Rsoma')

            V_rsoma_tots{end+1} = V_SANDI_tot;

        elseif strcmp(parameter,'fsoma')

            V_fsoma_tots{end+1} = V_SANDI_tot;

        elseif strcmp(parameter,'fneurite')

            V_fneurite_tots{end+1} = V_SANDI_tot;
        
        elseif strcmp(parameter,'Din')

            V_Din_tots{end+1} = V_SANDI_tot;

        elseif strcmp(parameter,'De')

            V_De_tots{end+1} = V_SANDI_tot;

        elseif strcmp(parameter,'fextra')

            V_fextra_tots{end+1} = V_SANDI_tot;
        end
    end
end

threasholds = 0.1:0.1:1;
index = 1:1:numel(threasholds);

plot_subjs_corr(threasholds, 0.5, V_GM_tot, V_fsoma_tots, V_rsoma_tots, V_CBF_tots, n_subjs, 'subsample');






%%  GLM voxel-wise 4

n_subjs=12;

V_CBF_tots={};

for i = 1:1:n_subjs%lst
    if i < 10
        i=num2str(i);
        %img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBF0_to_MNI/sub-00',i,'_run-01_bh_CBF0_2MNI.nii.gz');
        img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBF0toMNI/sub-00',i,'_task-bh_run-01_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBFbhtoMNI/sub-00',i,'_task-bh_run-01_dexi_volreg_asl_topup_CBF_map_2MNI.nii.gz');
    else
        i=num2str(i);
        img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBF0toMNI/sub-0',i,'_task-bh_run-01_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');
    end
    Vhdr = spm_vol(img_path_CBF);
    V_CBF_tot = spm_read_vols(Vhdr);
%     figure, imagesc(rot90(squeeze(V_CBF_tot(30,:,:)))) 
%     imshow(rot90(squeeze(V(30,:,:))),[])
    V_CBF_tots{end+1} = V_CBF_tot;
end



V_fsoma_tots={};

V_fneurite_tots={};

V_fextra_tots={};

V_Din_tots={};

V_De_tots={};

%parameter = 'De';



V_rsoma_tots={};

parameters = {'Rsoma', 'fsoma','fneurite','Din','De','fextra'};

for j = 1:length(parameters)
    parameter = parameters{j};
    for i = 1:1:n_subjs%lst
        if i < 10

            i=num2str(i);
            img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/SANDItoMNI/sub-00',i,'_run-01_SANDI-fit_', parameter,'_2MNI2mm.nii.gz');
        else

            i=num2str(i);
            img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/SANDItoMNI/sub-0',i,'_run-01_SANDI-fit_',parameter,'_2MNI2mm.nii.gz');
        end
        Vhdr = spm_vol(img_path_SANDI);
        V_SANDI_tot = spm_read_vols(Vhdr);

        if strcmp(parameter,'Rsoma')

            V_rsoma_tots{end+1} = V_SANDI_tot;

        elseif strcmp(parameter,'fsoma')

            V_fsoma_tots{end+1} = V_SANDI_tot;

        elseif strcmp(parameter,'fneurite')

            V_fneurite_tots{end+1} = V_SANDI_tot;
        
        elseif strcmp(parameter,'Din')

            V_Din_tots{end+1} = V_SANDI_tot;

        elseif strcmp(parameter,'De')

            V_De_tots{end+1} = V_SANDI_tot;

        elseif strcmp(parameter,'fextra')

            V_fextra_tots{end+1} = V_SANDI_tot;
        end
    end
end



[v_fsoma, v_rsoma, v_fneurite, v_fextra, v_De, v_Din]=vectorize(n_subjs, V_GM_tot, V_CBF_tots, V_fsoma_tots, V_fneurite_tots, V_fextra_tots, V_Din_tots, V_De_tots, V_rsoma_tots);

% glm

v_CBF_tr = v_CBF';
v_fsoma_tr = v_fsoma';
v_rsoma_tr = v_rsoma';
v_fneurite_tr = v_fneurite';
v_fextra_tr = v_fextra';
v_Din_tr = v_Din';
v_De_tr = v_De';


X = [v_fsoma_tr v_rsoma_tr v_fneurite_tr v_fextra_tr v_Din_tr v_De_tr]; %aggiungi fc
y = v_CBF_tr;

% X = [v_fsoma v_rsoma v_fneurite v_fextra v_Din v_De]; %aggiungi fc
% y = v_CBF;
%one subject v_fsoma size 1      239794
%multisubjects v_fsoma 232597           1
%figure, hist(y);

mdl = fitglm(X,y,'linear');

%length(v_CBF)

%mean on all the subjects

% figure, scatter(v_fneurite_tr, v_CBF_tr);
% index_permutation=randperm(length(v_fneurite_tr));
% %you must order CBF and SANDI in the same way
% v_SANDI_permutated=v_fneurite_tr(index_permutation);
% v_CBF_permutated = v_CBF_tr(index_permutation);
% figure,
% scatter(v_SANDI_permutated(1:1:round(numel(v_fneurite_tr)*0.01)),v_CBF_permutated(1:1:round(numel(v_CBF_tr)*0.01)))
% [r,p]=corrcoef(v_SANDI_permutated, v_CBF_permutated);



%% clustering voxel-wise k-means 5
%figure, scatter(v_rsoma, v_CBF);
z_rsoma = zscore(v_rsoma);
z_CBF = zscore(v_CBF);

X=[z_rsoma;z_CBF].';


vars_lst=[];
vars_lst_c1=[];

n_clusters=1:10;

for i=n_clusters
    [idx,C]=kmeans(X,i);
    vars_lst_c1(i)=var(C(:,1));
    vars_lst(i)=var(C(:,2));
    
end


%figure, scatter(n_clusters, vars_lst);
yinterp=interp1(n_clusters,vars_lst, n_clusters);
figure, plot(n_clusters, vars_lst,'o',n_clusters,yinterp,'-');
xlabel('# of clusters','FontSize',15);
ylabel('Variance c2','FontSize',15);
yinterp=interp1(n_clusters,vars_lst_c1, n_clusters);
figure, plot(n_clusters, vars_lst_c1,'o',n_clusters,yinterp,'-');
xlabel('# of clusters','FontSize',15);
ylabel('Variance c1','FontSize',15);

%PLOT THAT SHOWS CENTROIDS AND CLUSTERS: DO THEY CORRESPOND TO PARTICULAR
%REGIONS ? How do you recover spatial information once vectorized ?

%T=cluster(X,"maxclust",2);
%TRY
% for q=1:n_clusters
% C(q)=mean(X_ex(idx_ex==q,:));
% end

%% hierarchical clustering

z_rsoma = zscore(v_rsoma);
z_CBF = zscore(v_CBF);

X=[z_rsoma;z_CBF].';

Y=pdist(X);
Z=linkage(Y);
%%
n_clusters=1:10;
var_1=[];
var_2=[];
for j = n_clusters
    idx=cluster(Z,"maxclust",j);
    C_1=[];
    C_2=[];
    for i = 1:1:j
        %make the mean over y data and calculate ycentroid
        X_1=X(:,1);
        C_1(i)=mean(X_1(idx==i));
        %make the mean over x data and calculate xcentroid        
        X_2=X(:,2);
        C_2(i)=mean(X_2(idx==i));
        
    end
    var_1(i)=var(C_1);
    var_2(i)=var(C_2);
end

figure, scatter(n_clusters, var_1);
xlabel('# of clusters','FontWeight','bold');
ylabel('Variance Rsoma','FontWeight','bold');
figure, scatter(n_clusters, var_2);
xlabel('# of clusters','FontWeight','bold');
ylabel('Variance CBF','FontWeight','bold');

%highlight regions with significant pvalue





%% regional correlation (all subjects) 7
%load data
%%%%%%%%%%%%%
run='run-02';%CHANGE
%%%%%%%%%%%%%

subjects = importdata(strcat('/media/nas_rete/Vitality/code/subjs_DWI.txt'));

V_CBF_tots={};
V_CMRO2_tots={};
V_SANDI_tots={};
V_fsoma_tots={};


%run2
subjects([11,27])=[];

n_subjs=length(subjects);

% lst = 2:1:n_subjs;GLOVE
% lst(lst==8) = [];
% lst = 1:1:n_subjs;PRIN
% lst(lst==2) = [];
for i = 1:1:n_subjs %lst

    %if i < 10

        %i=num2str(i);
        subj=num2str(subjects{i});
        img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/250101/CBF0toMNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');        
        img_path_CMRO2 = strcat('/media/nas_rete/Vitality/maps2MNI/250101/CMRO20toMNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CMRO2_map_2MNI2mm.nii.gz');
        img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_Rsoma_2MNI2mm.nii.gz');
        img_path_fsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_fsoma_2MNI2mm.nii.gz');
         
        %img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBF0toMNI/sub-00',i,'_task-bh_',run,'_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CMRO20/sub-00',i,'_task-bh_run-01_dexi_volreg_asl_topup_CMRO2_map_on_MNI2mm.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBF0toMNIGrubb038/sub-00',i,'_run-01_bh_CBF0_2MNI2mm.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/CBF2MNI-resting/pil00',i,'_resting_task_CBF_map_2MNI.nii.gz');
        %img_path_CBF = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/perf/CBF0_2MNI/sub-0',i,'_ep2d_dexi_pc_v1_rs_CBF0_2MNI.nii.gz');
        %img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/SANDItoMNI/sub-00',i,'_',run,'_SANDI-fit_Rsoma_2MNI2mm.nii.gz');
        %img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/sub-00',i,'_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
        %img_path_SANDI = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil00',i,'_resting_SANDI-fit_Rsoma.nii.gz');
        %img_path_SANDI = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/dwi/SANDI_2MNI/sub-0',i,'_Rsoma_SANDI-fit_2MNI.nii.gz');
        %img_path_fsoma = strcat('/media/nas_rete/Vitality/maps2MNI/SANDItoMNI/sub-00',i,'_',run,'_SANDI-fit_fsoma_2MNI2mm.nii.gz');
 
    %else

        %i=num2str(i);
        %img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBF0toMNI/sub-0',i,'_task-bh_',run,'_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CMRO20/sub-0',i,'_task-bh_run-01_dexi_volreg_asl_topup_CMRO2_map_on_MNI2mm.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/CBF0toMNIGrubb038/sub-0',i,'_run-01_bh_CBF0_2MNI2mm.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/CBF2MNI-resting/pil0',i,'_resting_task_CBF_map_2MNI.nii.gz');
        %img_path_CBF = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/perf/CBF0_2MNI/sub-',i,'_ep2d_dexi_pc_v1_rs_CBF0_2MNI.nii.gz');
        
        %img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/SANDItoMNI/sub-0',i,'_',run,'_SANDI-fit_Rsoma_2MNI2mm.nii.gz');
        %img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/sub-0',i,'_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
        %img_path_SANDI = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil0',i,'_resting_SANDI-fit_Rsoma.nii.gz');
        %img_path_SANDI = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/dwi/SANDI_2MNI/sub-',i,'_Rsoma_SANDI-fit_2MNI.nii.gz');
        %img_path_fsoma = strcat('/media/nas_rete/Vitality/maps2MNI/SANDItoMNI/sub-0',i,'_',run,'_SANDI-fit_fsoma_2MNI2mm.nii.gz');


    %end

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

    Vhdr = spm_vol(img_path_SANDI);
    V_SANDI_tot = spm_read_vols(Vhdr);
    V_SANDI_tots{end+1} = V_SANDI_tot;

    Vhdr = spm_vol(img_path_fsoma);
    V_fsoma_tot = spm_read_vols(Vhdr);
    V_fsoma_tots{end+1} = V_fsoma_tot;


end

%% 
%process data
%%%%%%%%%%%%%
v='CBF';%CHANGE



if strcmp(v,'CBF')
    low_thr=5;
    up_thr=100;
else
    low_thr=50;
    up_thr=200;
end
%%%%%%%%%%%%%



if strcmp(v,'CBF')
    V_energy_tots=V_CBF_tots;
else
    V_energy_tots=V_CMRO2_tots;
end

%plot_regional_subjs_corr
regions = unique(V_atlas_tot(:));
n_regions = numel(regions);

reshaped_medians_energy=[];
reshaped_medians_SANDI=[];


medians_energy=[];
medians_SANDI=[];

V_GM = V_GM_tot;
V_GM(V_GM>0.5)=1;
V_GM(V_GM<1)=0;

n_voxels_lst = [];
n_voxels_lst_allsubjs=[];

corr_for_each_subj=[];
pvalue_for_each_subj=[];

corr_for_each_subj_reduced=[];
pvalue_for_each_subj_reduced=[];





% tic
for i = 1:1:n_subjs %here we loose information about the real number of subj
    
    V_energy_tot = V_energy_tots{i};
    V_SANDI_tot = V_SANDI_tots{i};
    V_fsoma_tot = V_fsoma_tots{i};

    V_fsoma_tot(V_fsoma_tot>0.15)=1;
    V_fsoma_tot(V_fsoma_tot<1)=0;

    medians_energy=[];
    medians_SANDI=[];
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

        
        V_energy = V_energy_tot;
        V_SANDI = V_SANDI_tot;

        V_energy_fsoma_masked = V_energy.*V_fsoma_tot;
        V_SANDI_fsoma_masked = V_SANDI.*V_fsoma_tot;

        V_energy_atlas = V_energy_fsoma_masked.*V_atlas;
        V_SANDI_atlas = V_SANDI_fsoma_masked.*V_atlas;




        V_energy_masked = V_energy_atlas.*V_GM;
        V_SANDI_masked = V_SANDI_atlas.*V_GM;


        
        
        
        V_energy_masked = V_energy_masked(:);
        V_SANDI_masked = V_SANDI_masked(:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        indices_energy = [];
        for ii = 1:numel(V_energy_masked)
            if V_energy_masked(ii)<low_thr || V_energy_masked(ii)>up_thr
                %v_CBF_reduced(ii)=[];
                indices_energy(end+1)=ii;
            end
        end

%         indices_rsoma_zeros = [];
%         for j = 1:numel(V_SANDI_masked)
%             if V_SANDI_masked(j)<5
%                 indices_rsoma_zeros(end+1)=j;
%             end
%         end
        disp(strcat('region', num2str(k),'sub-',num2str(i)))

%         indices = cat(2, indices_rsoma_zeros, indices_energy);
%         indices_unique = unique(indices);
%         %fsoma
         indices_unique = unique(indices_energy);

        v_energy_masked=V_energy_masked;
        v_SANDI_masked=V_SANDI_masked;

        v_energy_masked(indices_unique)=[];
        v_SANDI_masked(indices_unique)=[];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        n_voxels = numel(v_energy_masked);

        n_voxels_lst(end+1) = n_voxels;

        median_energy = median(v_energy_masked);
        median_SANDI = median(v_SANDI_masked);

        medians_energy(end+1) = median_energy;
        medians_SANDI(end+1) = median_SANDI;


        

    toc
    end

    n_voxels_lst_allsubjs(i,:) = n_voxels_lst;

    reshaped_medians_energy(i,:) = medians_energy;
    reshaped_medians_SANDI(i,:) = medians_SANDI;

    [r,p] = corrcoef(medians_energy, medians_SANDI, 'rows','complete');

    corr_for_each_subj(i)=r(2);
    pvalue_for_each_subj(i)=p(2);

%%%%
    idx=find(n_voxels_lst<prctile(n_voxels_lst,40));
    medians_energy_reduced=medians_energy;
    medians_SANDI_reduced=medians_SANDI;
    medians_energy_reduced(idx)=[];
    medians_SANDI_reduced(idx)=[];

    [r_reduced,p_reduced] = corrcoef(medians_energy_reduced, medians_SANDI_reduced, 'rows','complete');

    corr_for_each_subj_reduced(i)=r_reduced(2);
    pvalue_for_each_subj_reduced(i)=p_reduced(2);

%     reshaped_medians_CBF_reduced(i,:) = medians_CBF_reduced;
%     reshaped_medians_SANDI_reduced(i,:) = medians_SANDI_reduced;
%%%%

    disp(strcat('Finished subject', num2str(i),'Starting subject', num2str(i+1)))

end

% disp('Total time for running');
% toc


%%
%regional corr for each subj

mean_with_subjs_corr=mean(corr_for_each_subj);
[h,p,ci,stats]=ttest(atanh(corr_for_each_subj));

% corr_PRIN_for_each_subj_12subjs=corr_for_each_subj;
% pvalues_PRIN_for_each_subj_12subjs=pvalue_for_each_subj;

%plot regional correlation for each subject
figure, 
s=scatter(1:n_subjs,corr_for_each_subj,45);
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
yline(0,'--');
mean_corr=round(mean_with_subjs_corr,2);
mean_corr=num2str(mean_corr);
txt = {strcat('r = ',mean_corr,'*')};%,strcat('p-value = ',p_value_str)
%annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold','FontSize',12);
text(22,0.4,txt,'FontWeight', 'Bold','FontSize',12);
ylabel('r','FontSize',15,'FontWeight','bold');
xlabel('Subjects','FontSize',15,'FontWeight','bold');
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');

%%
% mean regional corr over all subjs

%mean over all subjects for each region

reshaped_medians_energy_tot=reshaped_medians_energy;
reshaped_medians_SANDI_tot=reshaped_medians_SANDI;

mean_energy_tot = nanmean(reshaped_medians_energy_tot,1);%median
mean_SANDI_tot = nanmean(reshaped_medians_SANDI_tot,1);

SE_energy_tot = nanstd(reshaped_medians_energy_tot,0,1)/sqrt(numel(reshaped_medians_energy_tot(:,1)));
SE_SANDI_tot = nanstd(reshaped_medians_SANDI_tot,0,1)/sqrt(numel(reshaped_medians_energy_tot(:,1)));

%other possible types of thresholding

%var_CBF=SE_CBF_tot./mean_CBF_tot;
var_SANDI=SE_SANDI_tot./mean_SANDI_tot;

% idx_lower_than_four=[];
% for i = 1:numel(reshaped_medians_SANDI(1,:))
%     if mean(reshaped_medians_SANDI(:,i)) < 4
%         idx_lower_than_four(end+1)=i;
%     end
% end

% nan_cols=[];
%thresholding for number of nans
% for i = 1:numel(reshaped_medians_SANDI(1,:))         
%         nans=sum(isnan(reshaped_medians_SANDI(:,i)));
%         nan_cols(end+1)=nans;
% end

% idx_outlier=[];
% for i = 1:numel(mean_CBF_tot)
%     if mean_CBF_tot(i)>170
%         idx_outlier(end+1)=i;
%     end
% end
% 

%thresholding step
labels_tot=unique(V_atlas_tot(:));
labels_tot(labels_tot==0)=[];
labels=labels_tot';

mean_energy=mean_energy_tot;
mean_SANDI=mean_SANDI_tot;
SE_energy=SE_energy_tot;
SE_SANDI=SE_SANDI_tot;


%remove low number of voxels regions
mean_n_voxels_tot = mean(n_voxels_lst_allsubjs,1);
regions=1:numel(mean_n_voxels_tot);
mean_n_voxels=mean_n_voxels_tot;

%plot n voxels distribution
% figure, bar(regions,mean_n_voxels_tot);
% xlabel('regions','FontSize',15);
% ylabel('mean number of voxels','FontSize',15);

idx_low_n_voxels=[];
for i=1:numel(regions)
    if mean_n_voxels_tot(i)<prctile(mean_n_voxels_tot,40)%600%prctile(mean_n_voxels_tot,40)%mean(mean_n_voxels)/5
        idx_low_n_voxels(end+1)=i;
    end
end

labels(idx_low_n_voxels)=[];
mean_n_voxels(idx_low_n_voxels)=[];
mean_SANDI(idx_low_n_voxels)=[];
mean_energy(idx_low_n_voxels)=[];
SE_energy(idx_low_n_voxels)=[];
SE_SANDI(idx_low_n_voxels)=[];
var_SANDI(idx_low_n_voxels)=[];



reshaped_medians_energy_tot(:,idx_low_n_voxels)=[];
reshaped_medians_SANDI_tot(:,idx_low_n_voxels)=[];

% 
% idx=find(mean_SANDI<12);
% labels(idx)=[];
% mean_n_voxels(idx)=[];
% mean_SANDI(idx)=[];
% mean_CBF(idx)=[];
% SE_SANDI(idx)=[];
% SE_CBF(idx)=[];
% 
% idx_high_SE=find(SE_SANDI>prctile(SE_SANDI,95));
idx_high_SE=find(var_SANDI>prctile(var_SANDI,75));
labels(idx_high_SE)=[];
mean_SANDI(idx_high_SE)=[];
mean_energy(idx_high_SE)=[];
SE_SANDI(idx_high_SE)=[];
SE_energy(idx_high_SE)=[];

reshaped_medians_energy_tot(:,idx_high_SE)=[];
reshaped_medians_SANDI_tot(:,idx_high_SE)=[];

% idx=find(mean_SANDI<4);
% labels(idx)=[];
% mean_n_voxels(idx)=[];
% mean_SANDI(idx)=[];
% mean_energy(idx)=[];
% SE_SANDI(idx)=[];
% SE_energy(idx)=[];

[r,p] = corrcoef(mean_energy, mean_SANDI, 'rows','complete');
corr_coef = round(r(2),2);
p_value = p(2);
corr_coef_str = num2str(corr_coef);
p_value_str = num2str(p_value);

parameter='Rsoma';
unit_of_measure='(\mum)';

if strcmp(v,'CBF')
    dependent_parameter='CBF(ml/100g/min)';
else
    dependent_parameter='CMRO2(\mumol/100g/min)';
end

P = polyfit(mean_SANDI,mean_energy,1);
yfit = P(1)*mean_SANDI+P(2);

figure, 
%s = scatter(mean_energy,mean_SANDI);
s = errorbar(mean_SANDI, mean_energy, SE_energy, SE_energy, SE_SANDI, SE_SANDI,'o');
hold on
plot(mean_SANDI,yfit,'--','LineWidth',3,'Color',"#000000");
xlabel(strcat(parameter, unit_of_measure),'FontSize',15,'FontWeight','bold');
ylabel(dependent_parameter,'FontSize',15,'FontWeight','bold');
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
txt = {strcat('r = ',corr_coef_str,'*')};%,strcat('p-value = ',p_value_str)
if strcmp(v,'CBF')
    text(12.5,60,txt, 'FontWeight', 'bold','FontSize',12);
else
    text(12.5,140,txt, 'FontWeight', 'bold','FontSize',12);
end
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'box','off')
grid on
% path_main = '/media/nas_rete/Work_manuela/DWI_En_modeling/main';
% path_to_image=strcat(path_main,'/regional_corr/','regional_corr_Rsoma_CBF_multisubjects_Vitality_PVE05_ERR_intervalincrease.png');
% saveas(figure(2), path_to_image);

fitglm(mean_SANDI, mean_energy)

%% bootstrapped subjects
n_boots=42;
raws_mean_energy=[];
raws_mean_SANDI=[];

all_ks=[];%check the random numbers are different at each n boot strap.
for n = 1:n_boots
    raws_energy=[];
    raws_SANDI=[];
    all_ks_one_run=[];
    for j=1:n_subjs
        k=randperm(n_subjs);
        all_ks_one_run(end+1)=k(1);
        raw_energy=reshaped_medians_energy(k(1),:);
        raw_SANDI=reshaped_medians_SANDI(k(1),:);

        raws_energy(j,:)=raw_energy;
        raws_SANDI(j,:)=raw_SANDI;
    end
    all_ks(n,:)=all_ks_one_run;
    mean_energy_tot = nanmean(raws_energy,1);%median
    mean_SANDI_tot = nanmean(raws_SANDI,1);
    
    var_SANDI=SE_SANDI_tot./mean_SANDI_tot;
    
    %thresholding step
    labels_tot=unique(V_atlas_tot(:));
    labels_tot(labels_tot==0)=[];
    labels=labels_tot';
    
    mean_energy=mean_energy_tot;
    mean_SANDI=mean_SANDI_tot;
    
    %remove low number of voxels regions
    mean_n_voxels_tot = mean(n_voxels_lst_allsubjs,1);
    regions=1:numel(mean_n_voxels_tot);
    mean_n_voxels=mean_n_voxels_tot;
    
    idx_low_n_voxels=[];
    for i=1:numel(regions)
        if mean_n_voxels_tot(i)<prctile(mean_n_voxels_tot,40)
            idx_low_n_voxels(end+1)=i;
        end
    end
    
    labels(idx_low_n_voxels)=[];
    mean_SANDI(idx_low_n_voxels)=[];
    mean_energy(idx_low_n_voxels)=[];
    var_SANDI(idx_low_n_voxels)=[];
    
    idx_high_SE=find(var_SANDI>prctile(var_SANDI,75));
    labels(idx_high_SE)=[];
    mean_SANDI(idx_high_SE)=[];
    mean_energy(idx_high_SE)=[];

    raws_mean_energy(n,:)=mean_energy;
    raws_mean_SANDI(n,:)=mean_SANDI;

end




%%
labels_run2_CBF=labels;
mean_SANDI_run2_CBF=mean_SANDI;
mean_energy_run2_CBF=mean_energy;
SE_SANDI_run2_CBF=SE_SANDI;
SE_energy_run2_CBF=SE_energy;
reshaped_medians_energy_run2_CBF=reshaped_medians_energy;
reshaped_medians_SANDI_run2_CBF=reshaped_medians_SANDI;

%% clustering regional scatterplot CMRO2 vs Rsoma 8

z_rsoma = zscore(mean_SANDI);
z_energy = zscore(mean_energy);

X=[z_rsoma;z_energy].';
%%
%clusters made by hand

% 

% %5 clusters for negative correlation
% X(1:35,1)=-10+randn(35,1);
% X(1:35,2)=10+randn(35,1);
% X(36:75,1)=10+randn(40,1);
% X(36:75,2)=-10+randn(40,1);
% X(76:115,1)=randn(40,1);
% X(76:115,2)=randn(40,1);
% X(116:135,1)=-20+randn(20,1);
% X(116:135,2)=20+randn(20,1);
% X(136:155,1)=-30+randn(20,1);
% X(136:155,2)=30+randn(20,1);

% %5 clusters for positive correlation
X(1:35,1)=-10+randn(35,1);
X(1:35,2)=-10+randn(35,1);
X(36:75,1)=10+randn(40,1);
X(36:75,2)=10+randn(40,1);
X(76:115,1)=randn(40,1);
X(76:115,2)=randn(40,1);
X(116:135,1)=20+randn(20,1);
X(116:135,2)=20+randn(20,1);
X(136:155,1)=30+randn(20,1);
X(136:155,2)=30+randn(20,1);

figure, scatter(X(:,1),X(:,2))



%only noise
% X=randn(75,2);

%plot points within a circle of radius r
% maxr = 0.001;
% %create data
% N = 20;
% r = rand(1,N) * maxr;
% theta = rand(1,N) * 2*pi;
% [x, y] = pol2cart(theta, r);
% %plot
% scatter(x, y, 'b*');
% hold on
% viscircles([0 0], maxr);
% hold off
% axis equal


%%
%%%%% Hierarchical clustering
Y=pdist(X);
Z=linkage(Y);

n_clusters=1:10;
var_1=[];
var_2=[];
for j = n_clusters
    idx=cluster(Z,"maxclust",j);
    C_1=[];
    C_2=[];
    for i = 1:1:j
        %make the mean over y data and calculate ycentroid
        X_1=X(:,1);
        C_1(i)=mean(X_1(idx==i));%it's like X_1(find(idx==10))
        %make the mean over x data and calculate xcentroid
        X_2=X(:,2);
        C_2(i)=mean(X_2(idx==i));

    end
    var_1(i)=var(C_1);
    var_2(i)=var(C_2);

end
%% k-means clustering (changing X)
slopes_n_boots=[];
slopes=[];

all_ks=[];
n_boots=42;
for n=1:n_boots
    z_rsoma = zscore(raws_mean_SANDI(n,:));
    z_energy = zscore(raws_mean_energy(n,:));

    X=[z_rsoma;z_energy].';
% % considering hand made samples
%     k_one_run=[];
%     new_Xs=[];
%     rng(n)
%     for j = 1:length(X(:,1))
%     
%     k=randperm(length(X(:,1)),1);
%     k_one_run(end+1)=k;
%     
%     new_X=X(k,:);
% 
% 
%     new_Xs(j,:)=new_X;
%     end
%     all_ks(n,:)=k_one_run;
% %     figure, hist(all_ks(n,:));
% %
    n_clusters=1:length(X);

    slopes=[];
    n
    for i=n_clusters
        %i
        rng(42);
        [idx,C]=kmeans(X,i);
        yfit=polyfit(C(:,1),C(:,2),1);
        slope=yfit(1);
        slopes(end+1)=slope;
        %figure, scatter(C(:,1),C(:,2));
        
    end
    slopes_n_boots(n,:)=slopes;
end

%% k-means clustering changing seed

var_1=[];
var_2=[];

n_clusters=2:length(X);
%idx_for_n_clusters=[]; %in case yoy run it once
n_clusterings=42;

slopes_n_kmeans=[];
slopes=[];
for k=1:n_clusterings
    slopes=[];
    disp(strcat('kmeans seed',num2str(k)))
    for i=n_clusters
        disp(strcat('Cluster',num2str(i)))
        rng(k);%42
        [idx,C]=kmeans(X,i);
        diff=[];
        for j=1:i

            diff(idx==j,:)=X(idx==j,:)-repmat(C(j,:),sum(idx==j),1);
        end
        sse(i)=sum(sum(diff.^2));


        %idx_for_n_clusters(:,i)=idx;
        var_1(i)=var(C(:,1));
        var_2(i)=var(C(:,2));
        C_1=C(:,1);
        C_2=C(:,2);

        P=polyfit(C(:,1),C(:,2),1); %if you have one point, linear fit starting from 0.
        slope=P(1);
        slopes(end+1)=slope;

        yfit = P(1)*C(:,1)+P(2);

%         if i==1 && k==1
% 
%             figure;
%             s=plot(X(idx==1,1),X(idx==1,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b");
%             %s.MarkerEdgeColor = 'b';
%             %s.MarkerFaceColor = [0 0.5 0.5];
%             hold on
%             plot(C_1,C_2,'kx',...
%                 'MarkerSize',15,'LineWidth',3)
%             hold on
%             plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
%             xlabel('Rsoma','FontSize',15,'FontWeight','bold');
%             ylabel('CMRO2','FontSize',15,'FontWeight','bold');
%             legend('Cluster 1','Centroids','Location','NW')
%             title 'Cluster Assignments and Centroids'
%             hold off
% 
%         elseif i==2 && k==1
% 
%             figure;
%             plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',10,'MarkerFaceColor',"r")
%             hold on
%             plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',10,'MarkerFaceColor',"b")
%             hold on
%             plot(C_1,C_2,'kx',...
%                 'MarkerSize',15,'LineWidth',3)
%             hold on
%             plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
%             xlabel('Rsoma','FontWeight','bold','FontSize',15);
%             ylabel('CBF','FontWeight','bold','FontSize',15);
%             legend('Cluster 1','Cluster 2','Centroids','Location','NW')
%             title 'Cluster Assignments and Centroids'
%             hold off
% 
% 
%         elseif i==3 && k==1
%             figure;
%             plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
%             hold on
%             plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
%             hold on
%             plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
%             hold on
%             plot(C_1,C_2,'kx',...
%                 'MarkerSize',15,'LineWidth',3)
%             hold on
%             plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
%             legend('Cluster 1','Cluster 2','Cluster 3','Centroids','Location','NW')
%             title 'Cluster Assignments and Centroids'
%             hold off
% 
% 
%         elseif i==4 && k==1
%             figure;
%             plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
%             hold on
%             plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
%             hold on
%             plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
%             hold on
%             plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',5,'MarkerFaceColor',"c")
%             hold on
%             plot(C_1,C_2,'kx',...
%                 'MarkerSize',15,'LineWidth',3)
%             hold on
%             plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
%             legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Location','NW')
%             title 'Cluster Assignments and Centroids'
%             hold off
% 
% 
%         elseif i==5 && k==1
% 
%             figure;
%             plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
%             hold on
%             plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
%             hold on
%             plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
%             hold on
%             plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',5,'MarkerFaceColor',"c")
%             hold on
%             plot(X(idx==5,1),X(idx==5,2),"o",'Color',"m",'MarkerSize',5,'MarkerFaceColor',"m")
%             hold on
%             plot(C_1,C_2,'kx',...
%                 'MarkerSize',15,'LineWidth',3)
%             hold on
%             plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
%             legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Cluster 5','Centroids',...
%                 'Location','NW')
%             title 'Cluster Assignments and Centroids'
%             hold off
% 
% 
%         elseif i==6 && k==1
% 
%             figure;
%             plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
%             hold on
%             plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
%             hold on
%             plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
%             hold on
%             plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',5,'MarkerFaceColor',"c")
%             hold on
%             plot(X(idx==5,1),X(idx==5,2),"o",'Color',"m",'MarkerSize',5,'MarkerFaceColor',"m")
%             hold on
%             plot(X(idx==6,1),X(idx==6,2),"o",'Color',"y",'MarkerSize',5,'MarkerFaceColor',"y")
%             hold on
%             plot(C_1,C_2,'kx',...
%                 'MarkerSize',15,'LineWidth',3)
%             hold on
%             plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
%             legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Centroids',...
%                 'Location','NW')
%             title 'Cluster Assignments and Centroids'
%             hold off
% 
%         elseif i==7 && k==1
% 
%             figure;
%             plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
%             hold on
%             plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
%             hold on
%             plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
%             hold on
%             plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',5,'MarkerFaceColor',"c")
%             hold on
%             plot(X(idx==5,1),X(idx==5,2),"o",'Color',"m",'MarkerSize',5,'MarkerFaceColor',"m")
%             hold on
%             plot(X(idx==6,1),X(idx==6,2),"o",'Color',"y",'MarkerSize',5,'MarkerFaceColor',"y")
%             hold on
%             plot(X(idx==7,1),X(idx==7,2),"o",'Color',"#0072BD",'MarkerSize',5,'MarkerFaceColor',"#0072BD")
%             hold on
%             plot(C_1,C_2,'kx',...
%                 'MarkerSize',15,'LineWidth',3)
%             hold on
%             plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
%             legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Cluster 7','Centroids',...
%                 'Location','NW')
%             title 'Cluster Assignments and Centroids'
%             hold off
%         end

    end

    slopes_n_kmeans(k,:)=slopes;

    %             if i==1
    %
    %                 figure;
    %                 s=plot(X(idx==1,1),X(idx==1,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b");
    %                 %s.MarkerEdgeColor = 'b';
    %                 %s.MarkerFaceColor = [0 0.5 0.5];
    %                 hold on
    %                 plot(C_1,C_2,'kx',...
    %                     'MarkerSize',15,'LineWidth',3)
    %                 hold on
    %                 plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
    %                 xlabel('Rsoma','FontSize',15,'FontWeight','bold');
    %                 ylabel('CMRO2','FontSize',15,'FontWeight','bold');
    %                 legend('Cluster 1','Centroids','Location','NW')
    %                 title 'Cluster Assignments and Centroids'
    %                 hold off
    %
    %             elseif i==2
    %
    %                 figure;
    %                 plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',10,'MarkerFaceColor',"r")
    %                 hold on
    %                 plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',10,'MarkerFaceColor',"b")
    %                 hold on
    %                 plot(C_1,C_2,'kx',...
    %                     'MarkerSize',15,'LineWidth',3)
    %                 hold on
    %                 plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
    %                 xlabel('Rsoma','FontWeight','bold','FontSize',15);
    %                 ylabel('CBF','FontWeight','bold','FontSize',15);
    %                 legend('Cluster 1','Cluster 2','Centroids','Location','NW')
    %                 title 'Cluster Assignments and Centroids'
    %                 hold off
    %
    %
    %             elseif i==3
    %                 figure;
    %                 plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
    %                 hold on
    %                 plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
    %                 hold on
    %                 plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
    %                 hold on
    %                 plot(C_1,C_2,'kx',...
    %                     'MarkerSize',15,'LineWidth',3)
    %                 hold on
    %                 plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
    %                 legend('Cluster 1','Cluster 2','Cluster 3','Centroids','Location','NW')
    %                 title 'Cluster Assignments and Centroids'
    %                 hold off
    %
    %
    %             elseif i==4
    %                 figure;
    %                 plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
    %                 hold on
    %                 plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
    %                 hold on
    %                 plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
    %                 hold on
    %                 plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',5,'MarkerFaceColor',"c")
    %                 hold on
    %                 plot(C_1,C_2,'kx',...
    %                     'MarkerSize',15,'LineWidth',3)
    %                 hold on
    %                 plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
    %                 legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Location','NW')
    %                 title 'Cluster Assignments and Centroids'
    %                 hold off
    %
    %
    %             elseif i==5
    %
    %                 figure;
    %                 plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
    %                 hold on
    %                 plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
    %                 hold on
    %                 plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
    %                 hold on
    %                 plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',5,'MarkerFaceColor',"c")
    %                 hold on
    %                 plot(X(idx==5,1),X(idx==5,2),"o",'Color',"m",'MarkerSize',5,'MarkerFaceColor',"m")
    %                 hold on
    %                 plot(C_1,C_2,'kx',...
    %                     'MarkerSize',15,'LineWidth',3)
    %                 hold on
    %                 plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
    %                 legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Cluster 5','Centroids',...
    %                     'Location','NW')
    %                 title 'Cluster Assignments and Centroids'
    %                 hold off
    %
    %
    %             elseif i==6
    %
    %                 figure;
    %                 plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
    %                 hold on
    %                 plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
    %                 hold on
    %                 plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
    %                 hold on
    %                 plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',5,'MarkerFaceColor',"c")
    %                 hold on
    %                 plot(X(idx==5,1),X(idx==5,2),"o",'Color',"m",'MarkerSize',5,'MarkerFaceColor',"m")
    %                 hold on
    %                 plot(X(idx==6,1),X(idx==6,2),"o",'Color',"y",'MarkerSize',5,'MarkerFaceColor',"y")
    %                 hold on
    %                 plot(C_1,C_2,'kx',...
    %                     'MarkerSize',15,'LineWidth',3)
    %                 hold on
    %                 plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
    %                 legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Centroids',...
    %                     'Location','NW')
    %                 title 'Cluster Assignments and Centroids'
    %                 hold off
    %
    %             elseif i==7
    %
    %                 figure;
    %                 plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
    %                 hold on
    %                 plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
    %                 hold on
    %                 plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
    %                 hold on
    %                 plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',5,'MarkerFaceColor',"c")
    %                 hold on
    %                 plot(X(idx==5,1),X(idx==5,2),"o",'Color',"m",'MarkerSize',5,'MarkerFaceColor',"m")
    %                 hold on
    %                 plot(X(idx==6,1),X(idx==6,2),"o",'Color',"y",'MarkerSize',5,'MarkerFaceColor',"y")
    %                 hold on
    %                 plot(X(idx==7,1),X(idx==7,2),"o",'Color',"#0072BD",'MarkerSize',5,'MarkerFaceColor',"#0072BD")
    %                 hold on
    %                 plot(C_1,C_2,'kx',...
    %                     'MarkerSize',15,'LineWidth',3)
    %                 hold on
    %                 plot(C(:,1),yfit,'--','LineWidth',3,'Color',"#000000");
    %                 legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Cluster 7','Centroids',...
    %                     'Location','NW')
    %                 title 'Cluster Assignments and Centroids'
    %                 hold off
    %             end
    %
    %     elseif i==8
    %         figure;
    %         plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',5,'MarkerFaceColor',"r")
    %         hold on
    %         plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',5,'MarkerFaceColor',"b")
    %         hold on
    %         plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',5,'MarkerFaceColor',"g")
    %         hold on
    %         plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',5,'MarkerFaceColor',"c")
    %         hold on
    %         plot(X(idx==5,1),X(idx==5,2),"o",'Color',"m",'MarkerSize',5,'MarkerFaceColor',"m")
    %         hold on
    %         plot(X(idx==6,1),X(idx==6,2),"o",'Color',"y",'MarkerSize',5,'MarkerFaceColor',"y")
    %         hold on
    %         plot(X(idx==7,1),X(idx==7,2),"o",'Color',"#0072BD",'MarkerSize',5,'MarkerFaceColor',"#0072BD")
    %         hold on
    %         plot(X(idx==8,1),X(idx==8,2),"o",'Color',"#D95319",'MarkerSize',5,'MarkerFaceColor',"#D95319")
    %         hold on
    %         plot(C_1,C_2,'kx',...
    %             'MarkerSize',15,'LineWidth',3)
    %         legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Cluster 7','Cluster 8','Centroids',...
    %             'Location','NW')
    %         title 'Cluster Assignments and Centroids'
    %         hold off
    %
    %     if i==9
    %         figure;
    %         plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',10,'MarkerFaceColor',"r")
    %         hold on
    %         plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',10,'MarkerFaceColor',"b")
    %         hold on
    %         plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',10,'MarkerFaceColor',"g")
    %         hold on
    %         plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',10,'MarkerFaceColor',"c")
    %         hold on
    %         plot(X(idx==5,1),X(idx==5,2),"o",'Color',"m",'MarkerSize',10,'MarkerFaceColor',"m")
    %         hold on
    %         plot(X(idx==6,1),X(idx==6,2),"o",'Color',"y",'MarkerSize',10,'MarkerFaceColor',"y")
    %         hold on
    %         plot(X(idx==7,1),X(idx==7,2),"o",'Color',"#0072BD",'MarkerSize',10,'MarkerFaceColor',"#0072BD")
    %         hold on
    %         plot(X(idx==8,1),X(idx==8,2),"o",'Color',"#D95319",'MarkerSize',10,'MarkerFaceColor',"#D95319")
    %         hold on
    %         plot(X(idx==9,1),X(idx==9,2),"o",'Color',"#EDB120",'MarkerSize',10,'MarkerFaceColor',"#EDB120")
    %         hold on
    %         plot(C_1,C_2,'kx',...
    %             'MarkerSize',15,'LineWidth',3)
    %         xlabel('Rsoma','FontWeight','bold','FontSize',15);
    %         ylabel('CBF','FontWeight','bold','FontSize',15);
    %         legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Cluster 7','Cluster 8','Cluster 9','Centroids',...
    %             'Location','NW','FontSize',12)
    %         title 'Cluster Assignments and Centroids'
    %         hold off
    %
    %     elseif i==10
    %         figure;
    %         plot(X(idx==1,1),X(idx==1,2),"o",'Color',"r",'MarkerSize',10,'MarkerFaceColor',"r")
    %         hold on
    %         plot(X(idx==2,1),X(idx==2,2),"o",'Color',"b",'MarkerSize',10,'MarkerFaceColor',"b")
    %         hold on
    %         plot(X(idx==3,1),X(idx==3,2),"o",'Color',"g",'MarkerSize',10,'MarkerFaceColor',"g")
    %         hold on
    %         plot(X(idx==4,1),X(idx==4,2),"o",'Color',"c",'MarkerSize',10,'MarkerFaceColor',"c")
    %         hold on
    %         plot(X(idx==5,1),X(idx==5,2),"o",'Color',"m",'MarkerSize',10,'MarkerFaceColor',"m")
    %         hold on
    %         plot(X(idx==6,1),X(idx==6,2),"o",'Color',"y",'MarkerSize',10,'MarkerFaceColor',"y")
    %         hold on
    %         plot(X(idx==7,1),X(idx==7,2),"o",'Color',"#0072BD",'MarkerSize',10,'MarkerFaceColor',"#0072BD")
    %         hold on
    %         plot(X(idx==8,1),X(idx==8,2),"o",'Color',"#D95319",'MarkerSize',10,'MarkerFaceColor',"#D95319")
    %         hold on
    %         plot(X(idx==9,1),X(idx==9,2),"o",'Color',"#EDB120",'MarkerSize',10,'MarkerFaceColor',"#EDB120")
    %         hold on
    %         plot(X(idx==10,1),X(idx==10,2),"o",'Color',"#7E2F8E",'MarkerSize',10,'MarkerFaceColor',"#7E2F8E")
    %         hold on
    %         plot(C_1,C_2,'kx',...
    %             'MarkerSize',15,'LineWidth',3)
    %         legend('Cluster 1','Cluster 2','Cluster 3', 'Cluster 4','Cluster 5','Cluster 6','Cluster 7','Cluster 8','Cluster 9','Cluster 10','Centroids',...
    %             'Location','NW','FontSize',12)
    %         title('Cluster Assignments and Centroids','Fontsize',15);
    %         xlabel('Rsoma','FontWeight','bold','FontSize',15);
    %         ylabel('CMRO2','FontWeight','bold','FontSize',15);
    %         hold off
    %    end

end
%%

%using variance of centroids
% var_1_perc=var_1/sum(var_1);
% var_2_perc=var_2/sum(var_2);
% figure, 
% s=scatter(n_clusters, var_1);
% xlabel('# of clusters','FontWeight','bold','FontSize',15);
% ylabel('1- Variance Rsoma','FontWeight','bold','FontSize',15);
% s.MarkerEdgeColor = 'b';
% s.MarkerFaceColor = [0 0.5 0.5];
% figure, 
% s=scatter(n_clusters, var_2);
% xlabel('# of clusters','FontWeight','bold','FontSize',15);
% ylabel('1-Variance CBF','FontWeight','bold','FontSize',15);
% s.MarkerEdgeColor = 'b';
% s.MarkerFaceColor = [0 0.5 0.5];
% 
% var_tot= (var_1+var_2)/(var_1(end)+var_2(end));
% figure, 
% s=scatter(n_clusters, 1-var_tot);
% xlabel('# of clusters','FontWeight','bold','FontSize',15);
% ylabel('1-Variance tot','FontWeight','bold','FontSize',15);
% s.MarkerEdgeColor = 'b';
% s.MarkerFaceColor = [0 0.5 0.5];

%elbow plot
figure, 
p=plot(sse,'-o');
p.MarkerEdgeColor = 'b';
p.MarkerFaceColor = [0 0.5 0.5];
xline(20,'--','LineWidth',3,'Color','r');
xline(9,'--','LineWidth',3,'Color','r');
ax=gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
grid on

mean_slopes=mean(slopes_n_boots,1);%slopes_n_kmeans

figure,
for i = 1:n_boots%n_clusterings
plot(slopes_n_boots(i,:));%slopes_n_kmeans
xlabel('# clusters','FontSize',15,'FontWeight','bold');
ylabel('GLM slope','FontSize',15,'FontWeight','bold');
% xline(5,'--');
hold on
end
hold on
s=plot(mean_slopes,'LineWidth',3,'Color','r');
grid on
legend(s,'mean');
%try to find a relationship with varying kmeans.





%% GLM analysis of clusters
X_clust=[X,idx];
slopes=[];
p_values=[];
corr_coefs=[];

n_clusters=9;
for n=1:n_clusters

    X_n_clust=X_clust;
    X_n_clust=X_n_clust(X_n_clust(:,3)==n,1:3);
%     figure,
%     scatter(X_n_clust(:,1), X_n_clust(:,2))
    
    %fitglm(X_n_clust(:,1), X_n_clust(:,2))
    yfit=polyfit(X_n_clust(:,1),X_n_clust(:,2),1);
    slope=yfit(1);
    slopes(end+1)=slope;

    [r,p]=corrcoef(X_n_clust(:,1), X_n_clust(:,2));
    p_values(end+1)=p(2);
    corr_coefs(end+1)=r(2);
end

%% Antonello's method to check if the linear relationship is driven by two or more clusters
X_clust=[X,idx];

X_n_clust_subsampled=[];

slopes=[];
p_values=[];
corr_coefs=[];


n_clusters=9;
for n=1:n_clusters
    X_n_clust_subsampled=[];
    
    
    for i = 1:n
    X_n_clust=X_clust;
    X_n_clust=X_n_clust(X_n_clust(:,3)==i,1:3);%| X_n_clust(:,3)==9 | X_n_clust(:,3)==5
    X_n_clust_subsampled=cat(1,X_n_clust_subsampled,X_n_clust);
    end
%     figure,
%     scatter(X_n_clust(:,1), X_n_clust(:,2))
    
    %fitglm(X_n_clust(:,1), X_n_clust(:,2))
    %X_n_clust_subsampled=X_n_clust_subsampled{1:n};
    yfit=polyfit(X_n_clust_subsampled(:,1),X_n_clust_subsampled(:,2),1);
    slope=yfit(1);
    slopes(end+1)=slope;

    [r,p]=corrcoef(X_n_clust_subsampled(:,1), X_n_clust_subsampled(:,2));
    p_values(end+1)=p(2);
    corr_coefs(end+1)=r(2);
end

figure, 
p=plot(1:n_clusters,slopes,'-o');
p.MarkerEdgeColor = 'b';
p.MarkerFaceColor = [0 0.5 0.5];
xlabel('# clusters','FontSize',12,'FontWeight','bold');
ylabel('Beta coefficient','FontSize',12,'FontWeight','bold');
grid on

%% plot spatially clusters found



%I method
% Cluster3D=NaN(size(V_atlas_tot));
% for i=1:size(idx_and_regions,1)
% Cluster3D(V_atlas_tot==idx_and_regions(i,2))=idx_and_regions(i,1);
% dimension(i)=sum(V_atlas_tot(:)==idx_and_regions(i,2));
% end

n=9;%choose clusters to plot
idx_tr=idx_for_n_clusters(:,n);
idx_tr=idx_tr';
%idx_and_regions=[idx_tr;labels];
idx_and_regions_tr=[idx_tr;labels]';



regions=unique(V_atlas_tot);
regions_tr=regions';
diff=setdiff(regions_tr,labels);
%II method

V_atlas = V_atlas_tot;
for ii = 1:length(V_atlas_tot(:))
    if any(0==V_atlas_tot(ii))
        V_atlas(ii)=0;
    elseif any(diff==V_atlas_tot(ii))
        V_atlas(ii)=0;
    else
        index=find(idx_and_regions_tr(:,2)==V_atlas_tot(ii));
        idx_col=idx_and_regions_tr(:,1);
        V_atlas(ii)=idx_col(index);%check! Is it selecting the idx of clusters ?
    end
end



fig=figure;
a=12:4:72;%1:4:44;%
for i = 1:length(a)
    subplot(4,4,i)
    imagesc(rot90(V_atlas(:,:,a(i)))),[0,n+1];
    axis equal
    axis off
    caxis([0,n+1])
end
h=axes(fig,'visible','off');
% map = [0.2 0.1 0.5
%     0.1 0.5 0.8
%     0.2 0.7 0.6
%     0.8 0.7 0.3
%     0.9 1 0];
map=[0.5 0.5 0.5
    1 0 0
    0 0 1
    0.1 1 0
     0 1 1
     1 0 1
     1 1 0
     0 0.4470 0.7410
     0.8500 0.3250 0.0980
     0.9290 0.6940 0.1250];
%     0.4940 0.1840 0.5560];
colormap(map)%hsv
%colorbar(h,'orientation','horizontal','Location','SouthOutside','FontSize',12);
%caxis([0,3])


%FAI PER IL CASO VOXEL WISE.

%% check for CBF, SANDI and n voxels distributions 9

regions=1:numel(mean_energy);
figure, bar(regions,std_energy);
xlabel('regions','FontSize',15);
ylabel('CBF standard deviation','FontSize',15);

regions=1:numel(mean_energy);

mean_of_means_CBF = mean(mean_energy);
mean_of_means_SANDI = mean(mean_SANDI);
mean_of_SE_CBF = mean(SE_energy);
mean_of_SE_SANDI=mean(SE_SANDI);

figure, bar(regions,SE_SANDI);
xlabel('regions','FontSize',15);
ylabel('CBF Standard Error','FontSize',15);
yline(mean_of_SE_SANDI,'--','LineWidth',3)
% [up,lo] = envelope(SE_CBF);
% hold on
% plot(regions,up,'linewidth',1.5)
% legend('SE dist','up','lo')
% hold off

regions=1:numel(mean_energy);
figure, bar(regions,mean_SANDI);
xlabel('regions','FontSize',15);
ylabel('means Rsoma','FontSize',15);



%number of voxels
mean_n_voxels = mean(n_voxels_lst_allsubjs,1);
regions=1:numel(mean_energy_tot);
figure, bar(regions,mean_n_voxels);
xlabel('regions','FontSize',15);
ylabel('mean number of voxels','FontSize',15);

[r,p] = corrcoef(mean_n_voxels, mean_SANDI, 'rows','complete');
corr_coef = r(2);
p_value = p(2);

corr_coef_str = num2str(corr_coef);
p_value_str = num2str(p_value);

figure, plot(mean_n_voxels, mean_SANDI,'o');
xlabel('Mean number of voxels','FontSize',15);
ylabel('Mean Rsoma','FontSize',15);
txt = {strcat('r = ',corr_coef_str),strcat('p-value = ',p_value_str)};
% text(60,12,txt,'FontWeight', 'Bold');
annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold');
% P = polyfit(mean_n_voxels,SE_SANDI,1);
% yfit = P(1)*mean_n_voxels+P(2);
% hold on;
% plot(mean_n_voxels,yfit,'b--');

%% linear model parameters

% mdl = fitlm(mean_CBF,mean_SANDI);
% plot(mdl);
% xlabel('CBF (ml/100g/min)');
% ylabel(strcat(parameter, ' (μm)'));

%% check correlation distribution 8

%correlation between CBF and SANDI in each region

%this is done for all the "original" regions.
% you could remove columns corresponding to regions removed in reshaped_medians_CBF
%and reshaped_medians_SANDI
for i = 1:(numel(labels_tot))
    [r,p] = corrcoef(reshaped_medians_energy_tot(:,i), reshaped_medians_SANDI_tot(:,i),'rows','complete');
    z(i)=r(2);
    pvalue(i)=p(2);
end



for i = 1:(numel(labels))
    [r,p] = corrcoef(reshaped_medians_energy(:,i), reshaped_medians_SANDI(:,i),'rows','complete');
    z(i)=r(2);
    pvalue(i)=p(2);
end

figure, hist(z);
xlabel('correlation, r');
ylabel('# regions');

z_accepted=[];
labels_accepted=[];
for i=1:numel(z)
    if pvalue(i)<0.05
        z_accepted(end+1)=z(i);
        labels_accepted(end+1)=labels(i);
    end
end

figure, hist(z_accepted);
xlabel('correlation, r');
ylabel('# regions');
      


%plot correlation values vs regions
x = 1:numel(z_accepted);
figure, 
s = scatter(x,z_accepted);
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
xlabel('regions');
ylabel('r');
yline(0.29,'-','Little if any');
yline(0.49,'-','Low');
yline(0.69,'-','Moderate');
yline(0.89,'-','High');
yline(-0.29,'-','Little if any');
yline(-0.49,'-','Low');
yline(-0.69,'-','Moderate');
yline(-0.89,'-','High');
% path_to_image=strcat(path_main,'/regional_corr/','regions_classification.png');
% saveas(figure, path_to_image);


% 0.90 to 1.00 Very high correlation
% 0.70 to 0.89 High correlation
% 0.50 to 0.69 Moderate correlation
% 0.30 to 0.49 Low correlation
% 0.00 to 0.29 Little if any correlation
very_high_corr=[];
high_corr=[];
moderate_corr=[];
low_corr=[];
very_low_corr=[];
low_neg_corr=[];
moderate_neg_corr=[];
high_neg_corr=[];
very_high_neg_corr=[];



for j=1:numel(z_accepted)
    if z_accepted(j)>=0.9
        very_high_corr(end+1)=j;
    elseif z_accepted(j)>=0.7 && z_accepted(j)<0.9
        high_corr(end+1)=j;
    elseif z_accepted(j)>=0.5 && z_accepted(j)<0.7
        moderate_corr(end+1)=j;
    elseif z_accepted(j)>=0.30 && z_accepted(j)<0.5
        low_corr(end+1)=j;
    elseif z_accepted(j)>=-0.30 && z_accepted(j)<0.30
        very_low_corr(end+1)=j;
    elseif z_accepted(j)>=-0.50 && z_accepted(j)<-0.30
        low_neg_corr(end+1)=j;
    elseif z_accepted(j)>=-0.7 && z_accepted(j)<-0.50
        moderate_neg_corr(end+1)=j;
    elseif z_accepted(j)>=-0.9 && z_accepted(j)<-0.70
        high_neg_corr(end+1)=j;
    else
        very_high_neg_corr(end+1)=j;
    end
end

n_very_high=numel(very_high_corr);
n_high=numel(high_corr);
n_moderate=numel(moderate_corr);
n_low=numel(low_corr);
n_very_low=numel(very_low_corr);
n_very_high_neg=numel(very_high_neg_corr);
n_high_neg=numel(high_neg_corr);
n_moderate_neg=numel(moderate_neg_corr);
n_low_neg=numel(low_neg_corr);
very_high_neg_corr=[];
n = vertcat(n_very_high, n_high, n_moderate, n_low, n_very_low, n_very_high_neg,n_high_neg,n_moderate_neg,n_low_neg);
   

% plot spacially correlation 8
% for i = 1:(numel(unique(V_atlas_tot(:)))-1)
%     [r,p] = corrcoef(reshaped_medians_CBF(:,i), reshaped_medians_SANDI(:,i),'rows','complete');
%     z(i)=r(2);
%     pvalue(i)=p(2);
% end

%remove columns corresponding to regions removed in reshaped_medians_CBF
%and reshaped_medians_SANDI
% for i = 1:(numel(labels))
%     [r,p] = corrcoef(reshaped_medians_CBF(:,i), reshaped_medians_SANDI(:,i),'rows','complete');
%     z(i)=r(2);
%     pvalue(i)=p(2);
% end

% regions_copy=regions;
% regions_copy=regions_copy';
% regions_copy(regions_copy==0)=[];
% corr_and_regions=[z;regions_copy].';
% 
% indices_ordered=1:numel(z);
% diff=setdiff(indices_ordered,regions);


corr_and_regions = [z_accepted;labels_accepted].';
diff=setxor(labels_accepted,unique(V_atlas_tot));
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
colormap pink
colorbar(h,'orientation','horizontal','Location','SouthOutside','FontSize',12);
sgtitle('Regional correlation map thresholded for p<0.05')
caxis([-1,+1])


%% number of cells density map (one subject) 10

%load fs map
%img_path_fsoma = strcat('/media/nas_rete/Vitality/clustering/sub-00',i,'_run-01_SANDI-fit_fsoma_2MNI.nii.gz');
img_path_fsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/sub-001_run-01_SANDI-fit_fsoma_2MNI2mm.nii.gz');
Vhdr = spm_vol(img_path_fsoma);
V_fsoma_tot = spm_read_vols(Vhdr);

%load rsoma map
img_path_rsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/sub-001_run-01_SANDI-fit_Rsoma_2MNI2mm.nii.gz');
%img_path_rsoma = strcat('/media/nas_rete/Vitality/clustering/sub-00',i,'_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
Vhdr = spm_vol(img_path_rsoma);
V_rsoma_tot = spm_read_vols(Vhdr);

V_rsoma = V_rsoma_tot;

for i = 1:length(V_rsoma(:))
    if V_rsoma(i) < 10
        V_rsoma(i)=0;
    end
end


%convert to m^3
V_rsoma = V_rsoma.*10^-3;

%voxel wise divide fs map over 4/3pir^3
fc_map = V_fsoma_tot./((4/3)*pi*V_rsoma.^3);
%figure, imagesc(fc_map(:,:,45));
for i = 1:length(fc_map(:))
    if fc_map(i)==Inf
        fc_map(i)=NaN;
%     elseif isnan(fc_map(i))
%         fc_map(i)=0;
    end
end

V_fsoma_tot(V_fsoma_tot>0.1)=1;
V_fsoma_tot(V_fsoma_tot<1)=0;
fc_masked = fc_map.*(V_fsoma);

figure; 
imagesc(rot90(fc_masked(:,:,45)));
colorbar(h,'orientation','horizontal','Location','SouthOutside','FontSize',9);
figure, imagesc(rot90(V_fsoma_tot(:,:,45)));
figure, imagesc(rot90(V_rsoma_tot(:,:,45)));





%% GLM regional calculation 12

run='run-01';%CHANGE
%%%%%%%%%%%%%

subjects = importdata(strcat('/media/nas_rete/Vitality/code/subjs_DWI.txt'));

%run2
subjects([11,27])=[];

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
    img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/250101/CBF0toMNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');
    img_path_CMRO2 = strcat('/media/nas_rete/Vitality/maps2MNI/250101/CMRO20toMNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CMRO2_map_2MNI2mm.nii.gz');

    img_path_rsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_Rsoma_2MNI2mm.nii.gz');
    img_path_fsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_fsoma_2MNI2mm.nii.gz');
    %img_path_fc = strcat('/media/nas_rete/Vitality/maps2MNI/SANDItoMNI/sub-00',i,'_run-01_SANDI-fit_fc_2MNI2mm.nii.gz');
    img_path_De = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_De_2MNI2mm.nii.gz');
    img_path_Din = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_Din_2MNI2mm.nii.gz');
    img_path_fneurite = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_fneurite_2MNI2mm.nii.gz');
    img_path_fextra = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_fextra_2MNI2mm.nii.gz');


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

%%%%%%%%%%%%%%%%%Calculate fc maps for many subjects
V_fc_tots={};

for i = 1:n_subjs    
    V_rsoma = V_rsoma_tots{i};
    V_fsoma = V_fsoma_tots{i};

    
    for i = 1:length(V_rsoma(:))
        if V_rsoma(i) < 10 %check
            V_rsoma(i)=0;
        end
    end
    
    
    %convert to m^3
    V_rsoma = V_rsoma.*10^-3;
    
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
% V_fc=V_fc_tots{1};
% figure, imagesc(V_fc(:,:,45))
% maxi=max(fc_map(:));
% caxis([0,maxi])
%%
v='CBF';

if strcmp(v,'CBF')
    V_energy_tots=V_CBF_tots;
    low_thr=5;
    up_thr=100;
else
    V_energy_tots=V_CMRO2_tots;
    low_thr=90;
    up_thr=200;
end

V_energy_matrix=[];
V_rsoma_matrix=[];
V_fsoma_matrix=[];
V_fneurite_matrix=[];
V_Din_matrix=[];
V_De_matrix=[];
V_fextra_matrix=[];
V_fc_matrix=[];

V_GM = V_GM_tot;
threashold = 0.5;

regions = unique(V_atlas_tot(:));
n_regions = numel(regions);

% reshaped_medians_CBF=[];
% reshaped_medians_SANDI=[];


medians_energy = [];
medians_fc= [];
medians_fextra= [];
medians_fneurite = [];
medians_fsoma= [];
medians_Din = [];
medians_De =[];
medians_rsoma = [];

V_GM = V_GM_tot;
V_GM(V_GM>0.5)=1;
V_GM(V_GM<1)=0;

reshaped_medians_energy=[];
reshaped_medians_De = [];
reshaped_medians_rsoma = [];
reshaped_medians_fsoma = [];
reshaped_medians_fneurite = [];
reshaped_medians_fextra = [];
reshaped_medians_Din = [];
reshaped_medians_fc = [];

n_voxels_lst = [];
n_voxels_lst_allsubjs=[];

% corr_for_each_subj=[];
% pvalue_for_each_subj=[];

% tic
for i = 1:1:n_subjs %1:1:length(lst)
    
    V_energy_tot = V_energy_tots{i};
    V_fc_tot = V_fc_tots{i};
    V_fsoma_tot = V_fsoma_tots{i};
    V_fneurite_tot = V_fneurite_tots{i};
    V_fextra_tot = V_fextra_tots{i};
    V_Din_tot = V_Din_tots{i};
    V_De_tot = V_De_tots{i};
    V_rsoma_tot = V_rsoma_tots{i};

    medians_energy=[];
    medians_fsoma=[];
    medians_fc=[];
    medians_fneurite=[];
    medians_fextra=[];
    medians_Din=[];
    medians_De=[];
    medians_rsoma=[];

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

        
        V_energy = V_energy_tot;
        V_fc = V_fc_tot;
        V_fextra = V_fextra_tot;
        V_fneurite = V_fneurite_tot;
        V_rsoma = V_rsoma_tot;
        V_fsoma = V_fsoma_tot;
        V_Din = V_Din_tot;
        V_De = V_De_tot;
        V_fsoma_to_mask=V_fsoma_tot;

        
        V_fsoma_to_mask(V_fsoma_to_mask>0.15)=1;
        V_fsoma_to_mask(V_fsoma_to_mask<1)=0;

        V_energy_atlas = V_energy.*V_atlas.*V_fsoma_to_mask;
        V_fc_atlas = V_fc.*V_atlas;%.*V_fsoma_to_mask;
        V_fextra_atlas = V_fextra.*V_atlas.*V_fsoma_to_mask;
        V_fneurite_atlas = V_fneurite.*V_atlas.*V_fsoma_to_mask;
        V_rsoma_atlas = V_rsoma.*V_atlas.*V_fsoma_to_mask;
        V_fsoma_atlas = V_fsoma.*V_atlas;%.*V_fsoma_to_mask;
        V_Din_atlas = V_Din.*V_atlas.*V_fsoma_to_mask;
        V_De_atlas = V_De.*V_atlas.*V_fsoma_to_mask;




        V_energy_masked = V_energy_atlas.*V_GM;
        V_fc_masked = V_fc_atlas.*V_GM;
        V_fextra_masked = V_fextra_atlas.*V_GM;
        V_fneurite_masked = V_fneurite_atlas.*V_GM;
        V_rsoma_masked = V_rsoma_atlas.*V_GM;
        V_fsoma_masked = V_fsoma_atlas.*V_GM;
        V_Din_masked = V_Din_atlas.*V_GM;
        V_De_masked = V_De_atlas.*V_GM;

        V_energy_masked = V_energy_masked(:);
        V_fc_masked = V_fc_masked(:);
        V_fextra_masked = V_fextra_masked(:);
        V_fneurite_masked = V_fneurite_masked(:);
        V_rsoma_masked = V_rsoma_masked(:);
        V_fsoma_masked = V_fsoma_masked(:);
        V_Din_masked = V_Din_masked(:);
        V_De_masked = V_De_masked(:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        indices_energy = [];
        for ii = 1:numel(V_energy_masked)
            if V_energy_masked(ii)<low_thr || V_energy_masked(ii)>up_thr %50 200
                indices_energy(end+1)=ii;
            end
        end

%         indices_rsoma_zeros = [];
%         for j = 1:numel(V_rsoma_masked)
%             if V_rsoma_masked(j)<5
%                 indices_rsoma_zeros(end+1)=j;
%             end
%         end
        disp(strcat('region', num2str(k)))

%         indices = cat(2, indices_rsoma_zeros, indices_CMRO2);
%         indices_unique = unique(indices);
%         %fsoma
       indices_unique = unique(indices_energy);

        v_energy_masked=V_energy_masked;
        v_fc_masked=V_fc_masked;
        v_fsoma_masked=V_fsoma_masked;
        v_fneurite_masked=V_fneurite_masked;
        v_fextra_masked=V_fextra_masked;
        v_rsoma_masked=V_rsoma_masked;
        v_Din_masked=V_Din_masked;
        v_De_masked=V_De_masked;

        v_energy_masked(indices_unique)=[];
        v_fc_masked(indices_unique)=[];
        v_fsoma_masked(indices_unique)=[];
        v_fneurite_masked(indices_unique)=[];
        v_fextra_masked(indices_unique)=[];
        v_Din_masked(indices_unique)=[];
        v_De_masked(indices_unique)=[];
        v_rsoma_masked(indices_unique)=[];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        n_voxels = numel(v_energy_masked);
      
        n_voxels_lst(end+1) = n_voxels;

        median_energy = median(v_energy_masked);
        median_fc = median(v_fc_masked);
        median_fextra = median(v_fextra_masked);
        median_fneurite = median(v_fneurite_masked);
        median_fsoma = median(v_fsoma_masked);
        median_rsoma = median(v_rsoma_masked);
        median_Din = median(v_Din_masked);
        median_De = median(v_De_masked);


        medians_energy(end+1) = median_energy;
        medians_fc(end+1) = median_fc;
        medians_fextra(end+1) = median_fextra;
        medians_fneurite(end+1) = median_fneurite;
        medians_fsoma(end+1) = median_fsoma;
        medians_Din(end+1) = median_Din;
        medians_De(end+1) = median_De;
        medians_rsoma(end+1) = median_rsoma;


        
    toc
    end

    n_voxels_lst_allsubjs(i,:) = n_voxels_lst;

    reshaped_medians_energy(i,:) = medians_energy;
    reshaped_medians_fc(i,:) = medians_fc;
    reshaped_medians_fextra(i,:) = medians_fextra;
    reshaped_medians_fneurite(i,:) = medians_fneurite;
    reshaped_medians_fsoma(i,:) = medians_fsoma;
    reshaped_medians_rsoma(i,:) = medians_rsoma;
    reshaped_medians_Din(i,:) = medians_Din;
    reshaped_medians_De(i,:) = medians_De;
    %[r,p] = corrcoef(medians_CBF, medians_SANDI, 'rows','complete');

%     corr_for_each_subj(i)=r(2);
%     pvalue_for_each_subj(i)=p(2);

    disp(strcat('Finished subject', num2str(i),'Starting subject', num2str(i+1)))

end

mean_energy_tot = nanmean(reshaped_medians_energy,1);%median
mean_fc_tot = nanmean(reshaped_medians_fc,1);
mean_fextra_tot = nanmean(reshaped_medians_fextra,1);
mean_fneurite_tot = nanmean(reshaped_medians_fneurite,1);
mean_fsoma_tot = nanmean(reshaped_medians_fsoma,1);
mean_rsoma_tot = nanmean(reshaped_medians_rsoma,1);
mean_Din_tot = nanmean(reshaped_medians_Din,1);
mean_De_tot = nanmean(reshaped_medians_De,1);

n=numel(reshaped_medians_energy(:,1));

SE_energy_tot = nanstd(reshaped_medians_energy,0,1)/sqrt(n);
SE_fc_tot = nanstd(reshaped_medians_fc,0,1)/sqrt(n);
SE_fextra_tot = nanstd(reshaped_medians_fextra,0,1)/sqrt(n);
SE_fneurite_tot = nanstd(reshaped_medians_fneurite,0,1)/sqrt(n);
SE_fsoma_tot = nanstd(reshaped_medians_fsoma,0,1)/sqrt(n);
SE_rsoma_tot = nanstd(reshaped_medians_rsoma,0,1)/sqrt(n);
SE_Din_tot = nanstd(reshaped_medians_Din,0,1)/sqrt(n);
SE_De_tot = nanstd(reshaped_medians_De,0,1)/sqrt(n);

var=SE_rsoma_tot./mean_rsoma_tot;

mean_n_voxels_tot = mean(n_voxels_lst_allsubjs,1);
regions=1:numel(mean_energy_tot);
mean_n_voxels=mean_n_voxels_tot;
% 
% figure, bar(regions,mean_n_voxels_tot);
% xlabel('regions','FontSize',15);
% ylabel('mean number of voxels','FontSize',15);


idx_low_n_voxels=[];
for i=1:numel(regions)
    if mean_n_voxels_tot(i)<prctile(mean_n_voxels_tot,40)%mean(mean_n_voxels)/5
        idx_low_n_voxels(end+1)=i;
    end
end


labels_tot=unique(V_atlas_tot(:));
labels=labels_tot;
labels(idx_low_n_voxels)=[];
mean_n_voxels=mean_n_voxels_tot;
mean_n_voxels(idx_low_n_voxels)=[];

mean_energy_tot(idx_low_n_voxels)=[];
mean_fc_tot(idx_low_n_voxels)=[];
mean_fextra_tot(idx_low_n_voxels)=[];
mean_fneurite_tot(idx_low_n_voxels)=[];
mean_fsoma_tot(idx_low_n_voxels)=[];
mean_rsoma_tot(idx_low_n_voxels)=[];
mean_Din_tot(idx_low_n_voxels)=[];
mean_De_tot(idx_low_n_voxels)=[];

SE_rsoma_tot(idx_low_n_voxels)=[];
SE_fsoma_tot(idx_low_n_voxels)=[];
SE_energy_tot(idx_low_n_voxels)=[];
%non so se serve fare la trasposta
var(idx_low_n_voxels)=[];
% idx=find(mean_rsoma_tot<12);
% 
% mean_energy_tot(idx)=[];
% mean_fc_tot(idx)=[];
% mean_fextra_tot(idx)=[];
% mean_fneurite_tot(idx)=[];
% mean_fsoma_tot(idx)=[];
% mean_rsoma_tot(idx)=[];
% mean_Din_tot(idx)=[];
% mean_De_tot(idx)=[];
% 
% SE_rsoma_tot(idx)=[];
% SE_fsoma_tot(idx)=[];
% SE_energy_tot(idx)=[];



%idx_high_SE=find(SE_rsoma_tot>prctile(SE_rsoma_tot,95));
idx_high_SE=find(var>prctile(var,75));
mean_energy_tot(idx_high_SE)=[];
mean_rsoma_tot(idx_high_SE) = [];
mean_fc_tot(idx_high_SE)=[];
mean_fextra_tot(idx_high_SE)=[];
mean_fneurite_tot(idx_high_SE)=[];
mean_fsoma_tot(idx_high_SE)=[];
%mean_rsoma_tot(idx_high_SE)=[];
mean_Din_tot(idx_high_SE)=[];
mean_De_tot(idx_high_SE)=[];

SE_rsoma_tot(idx_high_SE)=[];
SE_fsoma_tot(idx_high_SE)=[];
SE_energy_tot(idx_high_SE)=[];

%var(idx_high_SE)=[];

v_energy_tr = mean_energy_tot';
v_fc_tr = mean_fc_tot';
v_rsoma_tr = mean_rsoma_tot';
v_fneurite_tr = mean_fneurite_tot';
v_fextra_tr = mean_fextra_tot';
v_Din_tr = mean_Din_tot';
v_De_tr = mean_De_tot';
v_fsoma_tr = mean_fsoma_tot';

X = [v_rsoma_tr v_fneurite_tr v_fextra_tr v_Din_tr v_De_tr v_fsoma_tr]; %v_fsoma_tr
y = v_energy_tr;

X = [v_rsoma_tr v_fsoma_tr]; %aggiungi fc
y = v_energy_tr;

%one subject v_fsoma size 1      239794
%multisubjects v_fsoma 232597           1
%figure, hist(y);

mdl = fitglm(X,y,'linear');

%% correlation matrix
A=[v_rsoma_tr v_fneurite_tr v_fextra_tr v_Din_tr v_De_tr v_fsoma_tr]; 
corrcoef(A)

%% plot energy vs fsoma

% mean_fsoma=nanmean(reshaped_medians_fsoma,1);
% mean_energy=nanmean(reshaped_medians_energy,1);
% 
% n=
% SE_energy = nanstd(reshaped_medians_energy,0,1)/sqrt(numel(reshaped_medians_energy(:,1)));
% SE_SANDI = nanstd(reshaped_medians_fsoma,0,1)/sqrt(numel(reshaped_medians_energy(:,1)));

[r,p] = corrcoef(mean_energy_tot, mean_fsoma_tot, 'rows','complete');

corr_coef = round(r(2),2);
p_value = p(2);
corr_coef_str = num2str(corr_coef);
p_value_str = num2str(p_value);

P = polyfit(mean_fsoma_tot,mean_energy_tot,1);
yfit_energy = P(1)*mean_fsoma_tot+P(2);

parameter='fsoma';
%unit_of_measure='(%)';
if strcmp(v,'CBF')
    dependent_parameter='CBF (ml/100g/min)';%'CMRO2(\mumol/100g/min)';%'CBF(ml/100g/min)';%'CMRO2(\mumol/100g/min)';
else
    dependent_parameter='CMRO_2 (\mu mol/100g/min)';
end

figure, 
%s = scatter(mean_CBF,mean_SANDI);
s = errorbar(mean_fsoma_tot, mean_energy_tot, SE_energy_tot, SE_energy_tot, SE_fsoma_tot, SE_fsoma_tot,'o');
hold on
plot(mean_fsoma_tot,yfit_energy,'--','LineWidth',3,'Color',"#000000");
xlabel(strcat(parameter),'FontSize',15,'FontWeight','bold');
ylabel(dependent_parameter,'FontSize',15,'FontWeight','bold');
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
n_subjs_str=num2str(n_subjs);
txt = {strcat('r = ',corr_coef_str,'**')};%,strcat('p-value = ',p_value_str)
% text(60,12,txt,'FontWeight', 'Bold');
%annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold','FontSize',15);
text(0.40,67,txt, 'FontWeight', 'bold','FontSize',12);
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'box','off')
grid on













