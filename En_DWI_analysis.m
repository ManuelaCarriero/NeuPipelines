%Loading SANDI map to compare with vascular and metabolic map


%Each code section contains:
%voxel-wise analysis (one subject):
 %scatter plot 
 %corr curve
 %map varying PVE threashold
 %Rsoma (masked with fsoma) vs cbf
 %correlation between fsoma and CBF
%glm (remove unphysical rsoma and cbf values)
%regional correlation (one subject)
%regional correlation (all subjects)
%Check Rsoma distribution in Glove study
%number of cells density maps calculation
%regional analysis with numerical cellular density
%voxel-wise analysis (all subjects)
%fsoma vs GM

parameter='Rsoma';%'Rsoma';'fsoma';'De'; insert parameter name compatible with SANDI files.
unit_of_measure = ' (\mu m)';%'(\mu m^{2}/ms)';%'(\mu m)';(%);
label = strcat(parameter,unit_of_measure);

%Load GM 
img_path = '/storage/shared/Atlas/atlas_GM_on_MNI152_T1_2mm.nii.gz'; %presa da internet (prova a fare un FAST sulla MNI)
Vhdr = spm_vol(img_path);
V_GM_tot = spm_read_vols(Vhdr); 

%Loading hcp atlas
img_path='/storage/shared/Atlas/HCPMMP1_on_MNI152_ICBM2009a_nlin_res_on_MNI_T1_2mm.nii.gz';
Vhdr = spm_vol(img_path);
V_hcp_tot = spm_read_vols(Vhdr);

%% voxel-wise analysis (one subject)

%Load SANDI map

img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_SANDI_tot = spm_read_vols(Vhdr);

%Loading CBF map
img_path = '/media/nas_rete/Vitality/clustering/perf/sub-001_run-01_bh_CBF0_2MNI.nii.gz';
Vhdr = spm_vol(img_path);
V_CBF_tot = spm_read_vols(Vhdr);
%figure, imagesc(V_CBF(:,:,70))
V_CBF_tot=max(V_CBF_tot,0);%replace negative numbers with zeros
%figure, imagesc(V_CBF_tot(:,:,45))
cbf_unit_of_measure = ' (ml/g/min)';

threasholds = 0.1:0.1:1;
index = 1:1:numel(threasholds);

%%

min_rsoma=10;
min_cbf=20;

min_rsoma_str=num2str(min_rsoma);
min_cbf_str=num2str(min_cbf);
path=strcat('/media/nas_rete/Work_manuela/DWI_En_modeling/main/min_rsoma',min_rsoma_str,'andCBF',min_cbf_str);

if exist([path])==0
mkdir ([path])
end

figure(1),
for i = index

    V_GM = V_GM_tot;

    V_GM(V_GM>=threasholds(i))=1;
    V_GM(V_GM<threasholds(i))=0;%you could write V_GM(V_GM<1)=0.
    
    V_SANDI = V_SANDI_tot;
    %Mask SANDI map with grey matter
    V_SANDI_masked=V_SANDI.*V_GM;
    
    slices_SANDI{i} = V_SANDI_masked(:,:,45);

    V_CBF = V_CBF_tot;
    V_CBF_masked = V_CBF.*V_GM;

    %%%%
    v_SANDI=V_SANDI_masked(:);
    v_CBF=V_CBF_masked(:);

    indices_rsoma_zeros=[];
    for j = 1:numel(v_SANDI)
        if v_SANDI(j)<min_rsoma
            indices_rsoma_zeros(end+1)=j;
        end
    end
    
%     V_SANDI_masked(V_SANDI_masked==indices_rsoma_zeros)=0;
%     figure, imagesc(V_SANDI_masked(:,:,45)); in order to
%     show the image without removed points, you should take indices
% in 3d like V_SANDI_masked(1,3,4)

    %prepare vectors to substitute with reduced data points
    v_CBF_reduced = v_CBF;
    v_SANDI_reduced = v_SANDI;
    
    %Remove CBF=0, 95 and 5 prctile values from CBF distribution

    prctile99=prctile(v_CBF,99);
    prctile5=prctile(v_CBF,5);

    indices_cbf = [];
    for ii = 1:length(v_CBF)
        if v_CBF(ii)>prctile99 || v_CBF(ii)<prctile5 || v_CBF(ii)==0 || V_CBF(ii)<20 %|| ii == any(indices_rsoma_zeros)
            %v_CBF_reduced(ii)=[];
            indices_cbf(end+1)=ii;
        end
    end
    %try delete elements in 3d image and plot it again
    
    %Remove corresponding values of CBF to Rsoma=0
    indices_tot=cat(2,indices_rsoma_zeros, indices_cbf);

    indices_tot_unique = unique(indices_tot);

    v_CBF_reduced([indices_tot_unique])=[];
    v_SANDI_reduced([indices_tot_unique])=[];

%     figure, hist(v_CBF_reduced);
%     figure, hist(v_SANDI_reduced); %This has to be used for quality check
%     of SANDI maps!
    
    %%%%
%     figure, hist(v_CBF_reduced);
%     title('CBF distribution')
%     xlabel('CBF (ml/100g/min)');
%     ylabel('counts');
%     figure, hist(v_SANDI_reduced);
%     title('Rsoma distribution')
%     xlabel(label);
%     ylabel('counts');
    r = corrcoef(v_CBF_reduced, v_SANDI_reduced);
    corr_values(i) = r(2);

    

    subplot(4,4,i)
    scatter(v_CBF_reduced,v_SANDI_reduced)
    thr = num2str(threasholds(i));
    title(strcat("thr>=",thr))
    xlabel('CBF (ml/100g/min)');
    %ylabel(strcat(label,' masked with fsoma'));
    ylabel(strcat(label,' Signal'));
end
sgtitle("Rsoma vs CBF increasing GM PVE threashold (thr)");
set(figure(1), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
path_to_image=strcat(path,'/rsoma_vs_cbf_GMmask_min_rsoma',min_rsoma_str,'_cbf',min_cbf_str,'.png');
saveas(figure(1), path_to_image);




%correlation curve
figure(2), 
p = plot(threasholds, corr_values, '-o', 'LineWidth', 1, 'MarkerFaceColor','b', 'MarkerSize',3);
xlabel('PVE threashold');
ylabel('r');
title('corr\_coef Rsoma (masked with GM) vs CBF');
set(figure(2), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
path_to_image=strcat(path,'/corr_coef_rsomavscbf_GMmask_min_rsoma',min_rsoma_str,'_cbf',min_cbf_str,'.png');
saveas(figure(2), path_to_image);

%for the last PVE threashold
figure(3), hist(v_CBF_reduced);
title('CBF distribution')
xlabel('CBF (ml/100g/min)');
ylabel('counts');
path_to_image=strcat(path,'/CBFdistribution_min_rsoma',min_rsoma_str,'_cbf',min_cbf_str,'.png');
saveas(figure(3), path_to_image);

figure(4), hist(v_SANDI_reduced);
title('Rsoma distribution')
xlabel(label);
ylabel('counts');
path_to_image=strcat(path,'/Rsomadistribution_min_rsoma',min_rsoma_str,'_cbf',min_cbf_str,'.png');
saveas(figure(4), path_to_image);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




fig=figure(5);
for i = index
    subplot(4,4,i)
    imagesc(rot90(slices_SANDI{i}))
    thr = num2str(threasholds(i));
    title(strcat("thr>=",thr))
    axis equal
    axis off
end
caxis([0,1000]);
h = axes(fig,'visible','off'); 
sgtitle("Rsoma map masked with GM varying thr")
maxi=max(V_SANDI_masked, [],'all');
dy=maxi/10;
%labels = round(0:dy:maxi);
c = colorbar(h,'Position',[0.93 0.4 0.019 0.4]);  % 'TickLabels',[labels]) attach colorbar to h
% c.Limits=[0 maxi];
set(figure(5), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(5), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/rsoma_maps_threasholdingPVE.png');

%c = colorbar(h,'Location','southoutside','Position',[0.5 0.02 0.022 0.7]);
%c = colorbar(h,'Location','southoutside');
%set(c,'Position',[0.5 0.168 0.9 0.002])


%colormap(c,'jet')


%%
%mask Rsoma map with fsoma
img_path='/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_fsoma_2MNI.nii.gz';
Vhdr = spm_vol(img_path);
V_fsoma = spm_read_vols(Vhdr);

figure(6), 
for i = index
    
    V_GM = V_GM_tot;

    V_GM(V_GM>0)=1;
    V_GM(V_GM<=0)=0;

    
    V_SANDI_masked=V_SANDI.*V_GM;

    V_fsoma_mask = V_fsoma;
    
%     V_fsoma_tot(V_fsoma_tot>threasholds(i))=1;
%     V_fsoma_tot(V_fsoma_tot<threasholds(i))=0;
%     
%     V_fsoma_masked = V_fsoma_tot;

    V_Rsoma_per_fsoma = V_SANDI_masked.*(V_fsoma_mask>threasholds(i));

    slices_SANDI{i} = V_Rsoma_per_fsoma(:,:,45);

    V_CBF = V_CBF_tot;

    V_CBF_masked = V_CBF.*V_GM;
    
    % Plot and calculate correlation
    
    v_CBF=V_CBF_masked(:);
    v_Rsoma_per_fsoma = V_Rsoma_per_fsoma(:);
    
    r = corrcoef(v_CBF, v_Rsoma_per_fsoma); 
    corr_values(i) = r(2);
    
    subplot(4,4,i)
    scatter(v_CBF,v_Rsoma_per_fsoma)
    thr = num2str(threasholds(i));
    title(strcat("thr>=",thr))
    xlabel('CBF (ml/100g/min)');
    %ylabel(strcat(label,' masked with fsoma'));
    ylabel(strcat(label,' Signal'));
end
sgtitle("Rsoma masked with fsoma vs CBF");
set(figure(6), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(6), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/rsoma_vs_cbf_fsomamask.png');



figure(7), 
p = plot(threasholds, corr_values, '-o', 'LineWidth', 1, 'MarkerFaceColor','b', 'MarkerSize',3);
xlabel('PVE threashold');
ylabel('r');
title('corr\_coef Rsoma (masked with fsoma) vs CBF');
set(figure(7), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(7), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/corr_coef_rsomavscbf_fsomamask.png');


   
fig=figure(8);
for i = index
    subplot(4,4,i)
    imagesc(rot90(slices_PVE{i}))
    thr = num2str(threasholds(i));
    title(strcat("thr>=",thr))
    axis equal
    axis off
end
h = axes(fig,'visible','off'); 
sgtitle("Rsoma map masked with fsoma greater or equal than thr")
maxi=max(V_Rsoma_per_fsoma, [],'all');
dy=maxi/10;
labels = round(0:dy:maxi);
c = colorbar(h,'Position',[0.93 0.4 0.019 0.4], 'TickLabels',[0 2 3 5 6 8 9 11 12 14 16]);
set(figure(8), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(8), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/rsoma_maps_threasholdingfsomamap.png');

%%
%correlation between fsoma and CBF
V_GM = V_GM_tot;
V_GM(V_GM>0)=1;
V_GM(V_GM<=0)=0;

V_fsoma_masked = V_fsoma.*V_GM;
V_CBF_masked = V_CBF.*V_GM;

r = corrcoef(V_fsoma_masked(:), V_CBF_masked(:)); 

figure(9), 
scatter(V_CBF_tot(:), V_fsoma_masked(:))
ylabel('fsoma(%)');
%ylabel(strcat(label,' masked with fsoma'));
xlabel('CBF (ml/100g/min)');
title("corr\_coef="+r(2));
saveas(figure(9), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/corr_coef_rsomavscbf.png');

V_SANDI = V_fsoma;

figure(10), 
for i = index
    
    V_GM = V_GM_tot;

    V_GM(V_GM>threasholds(i))=1;
    V_GM(V_GM<threasholds(i))=0;    
    
    %Mask SANDI map with grey matter
    V_SANDI_masked=V_SANDI.*V_GM;

    V_CBF = V_CBF_tot;
    V_CBF_masked = V_CBF.*V_GM;
        
    % Plot and calculate correlation
    
    v_CBF=V_CBF_masked(:);
    v_SANDI=V_SANDI_masked(:);
    %v_Rsoma_per_fsoma = V_Rsoma_per_fsoma(:);
    
    r = corrcoef(v_CBF, v_SANDI); 
    corr_values(i) = r(2);
    
    subplot(4,4,i)
    scatter(v_CBF,v_SANDI)
    thr = num2str(threasholds(i));
    title(strcat("thr>=",thr))
    xlabel('CBF (ml/100g/min)');
    %ylabel(strcat(label,' masked with fsoma'));
    ylabel(strcat('fsoma(%)',' Signal'));
end
sgtitle("fsoma vs CBF increasing GM PVE threashold (thr)");
set(figure(10), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(10), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/fsoma_vs_cbf_GMmask.png');

figure(11), 
p = plot(threasholds, corr_values, '-o', 'LineWidth', 1, 'MarkerFaceColor','b', 'MarkerSize',3);
xlabel('PVE threashold');
ylabel('r');
title('corr\_coef fsoma (masked with GM) vs CBF');
set(figure(11), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(11), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/corr_coef_fsomavscbf_GMmask.png');



%%
%glm (remove unphysical rsoma and cbf values)

parameter='Rsoma';%'Rsoma';'fsoma';'De'; insert parameter name compatible with SANDI files.
% unit_of_measure = '(\mu m)';%'(\mu m^{2}/ms)';%'(\mu m)';(%);
% label = strcat(parameter,unit_of_measure);
img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_rsoma = spm_read_vols(Vhdr);


parameter='fsoma';
img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_fsoma = spm_read_vols(Vhdr);


parameter='De';
img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_De = spm_read_vols(Vhdr);

parameter='Din';
img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_Din = spm_read_vols(Vhdr);

parameter='fneurite';
img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_fneurite = spm_read_vols(Vhdr);

parameter='fextra';
img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_fextra = spm_read_vols(Vhdr);

img_path = '/media/nas_rete/Vitality/clustering/perf/sub-001_run-01_bh_CBF0_2MNI.nii.gz';
Vhdr = spm_vol(img_path);
V_CBF_tot = spm_read_vols(Vhdr);
%figure, imagesc(V_CBF(:,:,45))
V_CBF_tot=max(V_CBF_tot,0);%replace negative numbers with zeros



GM_path='/media/nas_rete/Vitality/clustering/SegSub01/GM2atlas.nii.gz';
Vhdr = spm_vol(GM_path);
V_GM_tot = spm_read_vols(Vhdr);
V_GM = V_GM_tot;

threashold = 0;

V_fsoma = V_fsoma.*(V_GM>threashold);
V_rsoma = V_rsoma.*(V_GM>threashold);
V_CBF_tot = V_CBF_tot.*(V_GM>threashold);

V_fsoma = V_fsoma(:);
V_rsoma = V_rsoma(:);
V_CBF_tot = V_CBF_tot(:);

X = [V_fsoma V_rsoma];
y = V_CBF_tot;

figure, hist(y);

mdl = fitglm(X,y,'linear','Distribution','normal');






%% regional correlation (one subject)
regions = unique(V_hcp_tot(:));
n_regions = numel(regions);

medians_CBF_onesubj=[];
medians_SANDI_onesubj=[];

for k = 2:numel(regions)
    tic
    disp(strcat('region', num2str(k)))
    V_hcp = V_hcp_tot;


    for ii = 1:length(V_hcp_tot(:))
        if V_hcp_tot(ii)==regions(k)
            V_hcp(ii)=1;
        else
            V_hcp(ii)=0;
        end
    end

    %Not modify original data
    V_GM = V_GM_tot;
    V_CBF = V_CBF_tot;
    V_SANDI = V_SANDI_tot;

    V_CBF_atlas = V_CBF.*V_hcp;
    V_SANDI_atlas = V_SANDI.*V_hcp;
    
    V_GM(V_GM > 0)=1;
    V_GM(V_GM==0)=0;

    V_CBF_masked = V_CBF_atlas.*V_GM;
    V_SANDI_masked = V_SANDI_atlas.*V_GM;
    

    V_CBF_masked = V_CBF_masked(:);
    V_SANDI_masked = V_SANDI_masked(:);
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

%     up_prctile=prctile(V_CBF_masked,95);
%     low_prctile=prctile(V_CBF_masked,5);%distribution is strongly peaked
%     around 0, hence calculated percentiles are wrong

    indices_cbf = [];
    for ii = 1:numel(V_CBF_masked)
        if V_CBF_masked(ii)<20 || V_CBF_masked(ii)>90
            %v_CBF_reduced(ii)=[];
            indices_cbf(end+1)=ii;
        end
    end

    indices_rsoma_zeros = [];
    for j = 1:numel(V_SANDI_masked)
        if V_SANDI_masked(j)<10
            indices_rsoma_zeros(end+1)=j;
        end
    end
    


    indices = cat(2, indices_rsoma_zeros, indices_cbf);
    indices_unique = unique(indices);

    v_CBF_masked=V_CBF_masked; 
    v_SANDI_masked=V_SANDI_masked;

    v_CBF_masked([indices_unique])=[];
    v_SANDI_masked([indices_unique])=[];
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    median_CBF = median(v_CBF_masked);
    median_SANDI = median(v_SANDI_masked);

    medians_CBF_onesubj(end+1) = median_CBF;
    medians_SANDI_onesubj(end+1) = median_SANDI;

    toc
    
end

[r,p] = corrcoef(medians_CBF_onesubj, medians_SANDI_onesubj);
corr_coef = r(2);
p_value = p(2);

corr_coef_str = num2str(corr_coef);
p_value_str = num2str(p_value);

figure(10), 
s = scatter(medians_CBF_onesubj,medians_SANDI_onesubj);
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
txt = {"one subject:", strcat('r = ',corr_coef_str),strcat('p-value = ',p_value_str)};
% text(60,12,txt,'FontWeight', 'Bold');
annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold');
ylabel(strcat(parameter, unit_of_measure));
xlabel(strcat("CBF",cbf_unit_of_measure));
% title_str = "one subject";
% title(title_str);
path_main = '/media/nas_rete/Work_manuela/DWI_En_modeling/main';
path_to_image=strcat(path_main,'/regional_corr/','regional_corr_Rsoma_CBF_onesubject.png');
saveas(figure(10), path_to_image);



%Clustering on corr plot ?

%% regional correlation (all subjects)




V_CBF_tots={};
V_SANDI_tots={};

n_subjs=12;

% lst = 2:1:n_subjs;GLOVE
% lst(lst==8) = [];
% lst = 1:1:n_subjs;PRIN
% lst(lst==2) = [];
for i = 1:1:n_subjs%lst

    if i < 10

        i=num2str(i);
        img_path_CBF = strcat('/media/nas_rete/Vitality/clustering/perf/sub-00',i,'_run-01_bh_CBF0_2MNI.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/CBF2MNI-resting/pil00',i,'_resting_task_CBF_map_2MNI.nii.gz');
        %img_path_CBF = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/perf/CBF0_2MNI/sub-0',i,'_ep2d_dexi_pc_v1_rs_CBF0_2MNI.nii.gz');
        img_path_SANDI = strcat('/media/nas_rete/Vitality/clustering/sub-00',i,'_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
        %img_path_SANDI = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil00',i,'_resting_SANDI-fit_Rsoma.nii.gz');
        %img_path_SANDI = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/dwi/SANDI_2MNI/sub-0',i,'_Rsoma_SANDI-fit_2MNI.nii.gz');

    else

        i=num2str(i);
        img_path_CBF = strcat('/media/nas_rete/Vitality/clustering/perf/sub-0',i,'_run-01_bh_CBF0_2MNI.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/CBF2MNI-resting/pil0',i,'_resting_task_CBF_map_2MNI.nii.gz');
        %img_path_CBF = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/perf/CBF0_2MNI/sub-',i,'_ep2d_dexi_pc_v1_rs_CBF0_2MNI.nii.gz');
        img_path_SANDI = strcat('/media/nas_rete/Vitality/clustering/sub-0',i,'_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
        %img_path_SANDI = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil0',i,'_resting_SANDI-fit_Rsoma.nii.gz');
        %img_path_SANDI = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/dwi/SANDI_2MNI/sub-',i,'_Rsoma_SANDI-fit_2MNI.nii.gz');


    end

    Vhdr = spm_vol(img_path_CBF);
    V_CBF_tot = spm_read_vols(Vhdr);
    V_CBF_tots{end+1} = V_CBF_tot;

    Vhdr = spm_vol(img_path_SANDI);
    V_SANDI_tot = spm_read_vols(Vhdr);
    V_SANDI_tots{end+1} = V_SANDI_tot;

end



regions = unique(V_hcp_tot(:));
n_regions = numel(regions);

reshaped_medians_CBF=[];
reshaped_medians_SANDI=[];


medians_CBF=[];
medians_SANDI=[];

V_GM = V_GM_tot;
V_GM(V_GM>0.5)=1;
V_GM(V_GM<1)=0;

for i = 1:1:n_subjs %1:1:length(lst)
    
    V_CBF_tot = V_CBF_tots{i};
    V_SANDI_tot = V_SANDI_tots{i};

    medians_CBF=[];
    medians_SANDI=[];
    for k = 2:numel(regions)
        tic


        V_hcp = V_hcp_tot;

        for ii = 1:length(V_hcp_tot(:))
            if V_hcp_tot(ii)==regions(k)
                V_hcp(ii)=1;
            else
                V_hcp(ii)=0;
            end
        end

        
        V_CBF = V_CBF_tot;
        V_SANDI = V_SANDI_tot;

        V_CBF_atlas = V_CBF.*V_hcp;
        V_SANDI_atlas = V_SANDI.*V_hcp;




        V_CBF_masked = V_CBF_atlas.*V_GM;
        V_SANDI_masked = V_SANDI_atlas.*V_GM;

        V_CBF_masked = V_CBF_masked(:);
        V_SANDI_masked = V_SANDI_masked(:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %     up_prctile=prctile(V_CBF_masked,95);
        %     low_prctile=prctile(V_CBF_masked,5);%La distribuzione ha un superpicco a 0 quindi i percentili calcolati non sono giusti

        indices_cbf = [];
        for ii = 1:numel(V_CBF_masked)
            if V_CBF_masked(ii)<20 || V_CBF_masked(ii)>90
                %v_CBF_reduced(ii)=[];
                indices_cbf(end+1)=ii;
            end
        end

        indices_rsoma_zeros = [];
        for j = 1:numel(V_SANDI_masked)
            if V_SANDI_masked(j)<10
                indices_rsoma_zeros(end+1)=j;
            end
        end
        disp(strcat('region', num2str(k)))

        indices = cat(2, indices_rsoma_zeros, indices_cbf);
        indices_unique = unique(indices);
%         %fsoma
%        indices_unique = unique(indices_cbf);

        v_CBF_masked=V_CBF_masked;
        v_SANDI_masked=V_SANDI_masked;

        v_CBF_masked(indices_unique)=[];
        v_SANDI_masked(indices_unique)=[];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        median_CBF = median(v_CBF_masked);
        median_SANDI = median(v_SANDI_masked);

        medians_CBF(end+1) = median_CBF;
        medians_SANDI(end+1) = median_SANDI;

    toc
    end

    reshaped_medians_CBF(i,:) = medians_CBF;
    reshaped_medians_SANDI(i,:) = medians_SANDI;

    disp(strcat('Finished subject', num2str(i),'Starting subject', num2str(i+1)))

end



%mean over all subjects for each region

mean_CBF = nanmean(reshaped_medians_CBF,1);%median
mean_SANDI = nanmean(reshaped_medians_SANDI,1);

% mean_CBF_PRIN=mean_CBF;
% mean_SANDI_PRIN=mean_SANDI;

[r,p] = corrcoef(mean_CBF, mean_SANDI, 'rows','complete');

corr_coef = r(2);
p_value = p(2);

corr_coef_str = num2str(corr_coef);
p_value_str = num2str(p_value);

figure(11), 
s = scatter(mean_CBF,mean_SANDI);
ylabel(strcat(parameter, unit_of_measure));
%ylabel("fsoma");
xlabel(strcat("CBF", cbf_unit_of_measure));
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
txt = {"12 subjects:", strcat('r = ',corr_coef_str),strcat('p-value = ',p_value_str)};
% text(60,12,txt,'FontWeight', 'Bold');
annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold');
path_main = '/media/nas_rete/Work_manuela/DWI_En_modeling/main';
path_to_image=strcat(path_main,'/regional_corr/','regional_corr_Rsoma_CBF_multisubjects_Vitality_PVE05.png');
saveas(figure(11), path_to_image);


% mdl = fitlm(mean_CBF,mean_SANDI);
% plot(mdl);
% xlabel('CBF (ml/100g/min)');
% ylabel(strcat(parameter, ' (μm)'));


%correlation between CBF and SANDI in each region
for i = 1:180
    r = corrcoef(reshaped_medians_CBF(:,i), reshaped_medians_SANDI(:,i),'rows','complete');
    z(i)=r(2);
end
%ADD P-VALUE and highlight regions with significant pvalue


figure, hist(z);
xlabel('correlation, r');
ylabel('counts');


%plot correlation values vs regions
x = 1:numel(z);
figure(12), 
s = scatter(x,z);
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
path_to_image=strcat(path_main,'/regional_corr/','regions_classification.png');
saveas(figure(12), path_to_image);


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



for j=1:180
    if z(j)>=0.9
        very_high_corr(end+1)=j;
    elseif z(j)>=0.7 && z(j)<0.9
        high_corr(end+1)=j;
    elseif z(j)>=0.5 && z(j)<0.7
        moderate_corr(end+1)=j;
    elseif z(j)>=0.30 && z(j)<0.5
        low_corr(end+1)=j;
    elseif z(j)>=-0.30 && z(j)<0.30
        very_low_corr(end+1)=j;
    elseif z(j)>=-0.50 && z(j)<-0.30
        low_neg_corr(end+1)=j;
    elseif z(j)>=-0.7 && z(j)<-0.50
        moderate_neg_corr(end+1)=j;
    elseif z(j)>=-0.9 && z(j)<-0.70
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
   

  


V_hcp = V_hcp_tot;
for ii = 1:length(V_hcp_tot(:))
    if any(0==V_hcp_tot(ii))
        V_hcp(ii)=0;
    else
        V_hcp(ii)=z(V_hcp_tot(ii));
    end
end


fig=figure;
a=28:4:68;%12:4:72;
for i = 1:length(a)
    subplot(4,4,i)
    imagesc(rot90(V_hcp(:,:,a(i)))),[0,8]
    axis equal
    axis off
end
h=axes(fig,'visible','off');
colormap gray
colorbar(h,'orientation','horizontal','Location','SouthOutside','FontSize',12);
caxis([-1,+1])



%%
% Check Rsoma distribution in Glove study

V_SANDI_maps = {};
n_subjs=12;

lst = 2:1:n_subjs;
lst(lst==8) = [];

for i = lst

    if i < 10

        i=num2str(i);
        img_path=strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil00',i,'_resting_SANDI-fit_Rsoma.nii.gz');
        %img=strcat('/media/nas_rete/GLOVE_STUDY/DDC/derivatives/pil00',i,'/resting/dwi/SANDI_Output_2MNI/pil00',i,'_resting_SANDI-fit_Rsoma.nii.gz');

    else
        i=num2str(i);
        img_path=strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil0',i,'_resting_SANDI-fit_Rsoma.nii.gz');
        %img=strcat('/media/nas_rete/GLOVE_STUDY/DDC/derivatives/pil0',i,'/resting/dwi/SANDI_Output_2MNI/pil0',i,'_resting_SANDI-fit_Rsoma.nii.gz');

    end
    
    Vhdr = spm_vol(img_path);
    V_SANDI = spm_read_vols(Vhdr);
    V_SANDI_maps{end+1} = V_SANDI;
end

for i = 1:1:length(lst)

    V_SANDI = V_SANDI_maps{i};
    
    V_SANDI_vec = V_SANDI(:);
    
    indices_rsoma_zeros = [];
    for j = 1:numel(V_SANDI_vec)
        if V_SANDI_vec(j)<10
            indices_rsoma_zeros(end+1)=j;
        end
    end

    V_SANDI_vec([indices_rsoma_zeros])=[];
    
    path = '/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/rsoma_dist';
    if exist([path])==0
    mkdir ([path])
    end
    
    i = lst(i);
    if i < 10

        i=num2str(i);
    
        path_to_image = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/rsoma_dist/hist_pil00',i,'_Rsoma.jpg');
        f=figure('visible','off');
        hist(V_SANDI_vec);
        saveas(f, path_to_image);

    else
        
        i=num2str(i);

        path_to_image = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/rsoma_dist/hist_pil0',i,'_Rsoma.jpg');
        f=figure('visible','off');
        hist(V_SANDI_vec);
        saveas(f, path_to_image);
    end


end

%% number of cells density map

%load fs map
%img_path_fsoma = strcat('/media/nas_rete/Vitality/clustering/sub-00',i,'_run-01_SANDI-fit_fsoma_2MNI.nii.gz');
img_path_fsoma = strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_fsoma_2MNI.nii.gz');
Vhdr = spm_vol(img_path_fsoma);
V_fsoma_tot = spm_read_vols(Vhdr);

%load rsoma map
img_path_rsoma = strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
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


fc_masked = fc_map.*(V_fsoma_tot>0.1);

figure; 
imagesc(rot90(fc_masked(:,:,45)));
colorbar(h,'orientation','horizontal','Location','SouthOutside','FontSize',9);
figure, imagesc(rot90(V_fsoma_tot(:,:,45)));
figure, imagesc(rot90(V_rsoma_tot(:,:,45)));




img_path_CBF = '/media/nas_rete/Vitality/clustering/perf/sub-001_run-01_bh_CBF0_2MNI.nii.gz';
Vhdr = spm_vol(img_path_CBF);
V_CBF_tot = spm_read_vols(Vhdr);

CBF_masked = V_CBF_tot.*(V_fsoma_tot>0.1);

indices_cbf=[];
for i = 1:length(CBF_masked(:))
    if CBF_masked(i)>90 || CBF_masked(i)<20
        indices_cbf(end+1)=i;
    end
end

indices_fc = [];
for i = 1:length(fc_masked(:))
    if isnan(fc_masked(i))
        indices_fc(end+1)=i;
    end
end

indices = cat(2, indices_cbf, indices_fc);
CBF_masked(indices)=[];
fc_masked(indices)=[];


figure, scatter(fc_masked(:),CBF_masked(:));
xlabel('fc');
ylabel('CBF');
[r,p]=corrcoef(fc_masked(:),CBF_masked(:));

%%
% regional analysis with numerical cellular density
V_CBF_tots={};
V_rsoma_tots={};
V_fsoma_tots={};

n_subjs=12;

% lst = 2:1:n_subjs;
% lst(lst==8) = [];

for i = 1:1:n_subjs%lst

    if i < 10

        i=num2str(i);
        img_path_CBF = strcat('/media/nas_rete/Vitality/clustering/perf/sub-00',i,'_run-01_bh_CBF0_2MNI.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/CBF2MNI-resting/pil00',i,'_resting_task_CBF_map_2MNI.nii.gz');
        img_path_rsoma = strcat('/media/nas_rete/Vitality/clustering/sub-00',i,'_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
        img_path_fsoma = strcat('/media/nas_rete/Vitality/clustering/sub-00',i,'_run-01_SANDI-fit_fsoma_2MNI.nii.gz');
        %img_path_rsoma = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil00',i,'_resting_SANDI-fit_Rsoma.nii.gz')
        %img_path_fsoma = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil00',i,'_resting_SANDI-fit_fsoma.nii.gz');

    else

        i=num2str(i);
        img_path_CBF = strcat('/media/nas_rete/Vitality/clustering/perf/sub-0',i,'_run-01_bh_CBF0_2MNI.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/CBF2MNI-resting/pil0',i,'_resting_task_CBF_map_2MNI.nii.gz');
        img_path_rsoma = strcat('/media/nas_rete/Vitality/clustering/sub-0',i,'_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
        img_path_fsoma = strcat('/media/nas_rete/Vitality/clustering/sub-0',i,'_run-01_SANDI-fit_fsoma_2MNI.nii.gz');
        %img_path_rsoma = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil0',i,'_resting_SANDI-fit_Rsoma.nii.gz')
        %img_path_fsoma = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil0',i,'_resting_SANDI-fit_fsoma.nii.gz');

    end

    Vhdr = spm_vol(img_path_CBF);
    V_CBF_tot = spm_read_vols(Vhdr);
    V_CBF_tots{end+1} = V_CBF_tot;

    Vhdr = spm_vol(img_path_rsoma);
    V_rsoma_tot = spm_read_vols(Vhdr);
    V_rsoma_tots{end+1} = V_rsoma_tot;

    Vhdr = spm_vol(img_path_fsoma);
    V_fsoma_tot = spm_read_vols(Vhdr);
    V_fsoma_tots{end+1} = V_fsoma_tot;
end

regions = unique(V_hcp_tot(:));
n_regions = numel(regions);

reshaped_medians_CBF=[];
reshaped_medians_fc=[];



medians_CBF=[];
medians_fc=[];



for j = 1:1:n_subjs

    V_CBF_tot = V_CBF_tots{j};
    V_rsoma_tot = V_rsoma_tots{j};
    V_fsoma_tot = V_fsoma_tots{j};

    medians_CBF=[];
    medians_fc=[];
    for k = 2:numel(regions)
        %tic


        V_hcp = V_hcp_tot;

        for ii = 1:length(V_hcp_tot(:))
            if V_hcp_tot(ii)==regions(k)
                V_hcp(ii)=1;
            else
                V_hcp(ii)=0;
            end
        end



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

        for i = 1:length(fc_map(:))
            if fc_map(i)==Inf
                fc_map(i)=NaN;
            end
        end
        
        V_fsoma = V_fsoma_tot;
        fc_masked = fc_map.*(V_fsoma_tot>0.1);

        CBF_masked = V_CBF_tot.*(V_fsoma_tot>0.1);

%         figure, hist(fc_masked(:));
%         figure, hist(CBF_masked(:));

        fc_masked = fc_masked.*V_hcp;
        CBF_masked = CBF_masked.*V_hcp;

        indices_cbf=[];
        for i = 1:length(CBF_masked(:))
            if CBF_masked(i)>90 || CBF_masked(i)<20
                indices_cbf(end+1)=i;
            end
        end

        indices_fc = [];
        for i = 1:length(fc_masked(:))
            if isnan(fc_masked(i))
                indices_fc(end+1)=i;
            end
        end

        indices = cat(2, indices_cbf, indices_fc);

        v_CBF_masked=CBF_masked;
        v_fc_masked=fc_masked;

        v_CBF_masked(indices)=[];

        v_fc_masked(indices)=[];
% 
%         figure,hist(v_CBF_masked(:));
%         figure, hist(v_fc_masked(:));



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        regions_processed = [];
        median_CBF = median(v_CBF_masked);
        median_fc = median(v_fc_masked);

        medians_CBF(end+1) = median_CBF;
        medians_fc(end+1) = median_fc;
        regions_processed(k-1)=k;
        disp(regions_processed)
        %toc
    end

    reshaped_medians_CBF(j,:) = medians_CBF;
    reshaped_medians_fc(j,:) = medians_fc;

    disp(strcat('Finished subject', num2str(j),'Starting subject', num2str(j+1)))

end

mean_CBF = nanmean(reshaped_medians_CBF,1);
mean_fc = nanmean(reshaped_medians_fc,1);

[r,p] = corrcoef(mean_CBF, mean_fc, 'rows','complete');

corr_coef = r(2);
p_value = p(2);

corr_coef_str = num2str(corr_coef);
p_value_str = num2str(p_value);

figure(11), 
s = scatter(mean_CBF,mean_fc);
%ylabel(strcat(parameter, unit_of_measure));
%ylabel(parameter);
ylabel("fc");
xlabel(strcat("CBF", cbf_unit_of_measure));
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
subjs=num2str(n_subjs);
txt = {strcat(subjs,"subjects:"), strcat('r = ',corr_coef_str),strcat('p-value = ',p_value_str)};
% text(60,12,txt,'FontWeight', 'Bold');
annotation('textbox',[.6,.2,.30,.13], 'String', txt, 'FontWeight', 'Bold');
% path_main = '/media/nas_rete/Work_manuela/DWI_En_modeling/main';
% path_to_image=strcat(path_main,'/regional_corr/','regional_corr_Rsoma_CBF_multisubjects_GLOVE.png');
% saveas(figure(11), path_to_image);

%% voxel wise analysis (all  subjects)

V_CBF_tots={};
V_rsoma_tots={};
V_fsoma_tots={};

n_subjs=12;

% lst = 2:1:n_subjs;
% lst(lst==8) = [];

for i = 1:1:n_subjs%lst

    if i < 10

        i=num2str(i);
        img_path_CBF = strcat('/media/nas_rete/Vitality/clustering/perf/sub-00',i,'_run-01_bh_CBF0_2MNI.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/CBF2MNI-resting/pil00',i,'_resting_task_CBF_map_2MNI.nii.gz');
        img_path_rsoma = strcat('/media/nas_rete/Vitality/clustering/sub-00',i,'_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
        img_path_fsoma = strcat('/media/nas_rete/Vitality/clustering/sub-00',i,'_run-01_SANDI-fit_fsoma_2MNI.nii.gz');
        %img_path_rsoma = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil00',i,'_resting_SANDI-fit_Rsoma.nii.gz')
        %img_path_fsoma = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil00',i,'_resting_SANDI-fit_fsoma.nii.gz');

    else

        i=num2str(i);
        img_path_CBF = strcat('/media/nas_rete/Vitality/clustering/perf/sub-0',i,'_run-01_bh_CBF0_2MNI.nii.gz');
        %img_path_CBF = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/CBF2MNI-resting/pil0',i,'_resting_task_CBF_map_2MNI.nii.gz');
        img_path_rsoma = strcat('/media/nas_rete/Vitality/clustering/sub-0',i,'_run-01_SANDI-fit_Rsoma_2MNI.nii.gz');
        img_path_fsoma = strcat('/media/nas_rete/Vitality/clustering/sub-0',i,'_run-01_SANDI-fit_fsoma_2MNI.nii.gz');
        %img_path_rsoma = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil0',i,'_resting_SANDI-fit_Rsoma.nii.gz')
        %img_path_fsoma = strcat('/media/nas_rete/GLOVE_STUDY/DDC/registered/SANDI/pil0',i,'_resting_SANDI-fit_fsoma.nii.gz');

    end

    Vhdr = spm_vol(img_path_CBF);
    V_CBF_tot = spm_read_vols(Vhdr);
    V_CBF_tots{end+1} = V_CBF_tot;

    Vhdr = spm_vol(img_path_rsoma);
    V_rsoma_tot = spm_read_vols(Vhdr);
    V_rsoma_tots{end+1} = V_rsoma_tot;

    Vhdr = spm_vol(img_path_fsoma);
    V_fsoma_tot = spm_read_vols(Vhdr);
    V_fsoma_tots{end+1} = V_fsoma_tot;
end

regions = unique(V_hcp_tot(:));
n_regions = numel(regions);

reshaped_medians_CBF=[];
reshaped_medians_fc=[];



medians_CBF=[];
medians_rsoma=[];

%make the mean of all subjects
%mask
%corr
for i = 1:1:n_subjs
    v = V_CBF_tots{i};
    v=v(:);
    medians_CBF(i,:)=v;
end
mean_CBF = nanmean(medians_CBF,1);

for i = 1:1:n_subjs
    medians_rsoma(end+1) = V_rsoma_tots{i};
end



%% fsoma vs GM
