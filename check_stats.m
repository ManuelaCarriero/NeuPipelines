%% Comparison bw UniCh and UniGe 
% edited by Alessandra Caporale, 14/11/2024

run='run-02';%CHANGE
%%%%%%%%%%%%%

subjects = importdata(strcat('/media/nas_rete/Vitality/code/subjs_DWI.txt'));

V_CBF_tots={};
V_CMRO2_tots={};
V_SANDI_tots={};
V_fsoma_tots={};


%run2
subjects([11,13,27])=[];

n_subjs=length(subjects);

for i = 1:1:n_subjs 


    subj=num2str(subjects{i});
    img_path_CBF = strcat('/media/nas_rete/Vitality/maps2MNI/250101/CBF0toMNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CBF_map_2MNI2mm.nii.gz');
    img_path_CMRO2 = strcat('/media/nas_rete/Vitality/maps2MNI/250101/CMRO20toMNI/',subj,'_task-bh_',run,'_dexi_volreg_asl_topup_CMRO2_map_2MNI2mm.nii.gz');
    img_path_SANDI = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_Rsoma_2MNI2mm.nii.gz');
    img_path_fsoma = strcat('/media/nas_rete/Vitality/maps2MNI/250101/SANDItoMNI/',subj,'_',run,'_SANDI-fit_fsoma_2MNI2mm.nii.gz');

    Vhdr = spm_vol(img_path_CBF);
    V_CBF_tot = spm_read_vols(Vhdr);
    V_CBF_tots{end+1} = V_CBF_tot;

    Vhdr = spm_vol(img_path_CMRO2);
    V_CMRO2_tot = spm_read_vols(Vhdr);
    V_CMRO2_tots{end+1} = V_CMRO2_tot;

    Vhdr = spm_vol(img_path_SANDI);
    V_SANDI_tot = spm_read_vols(Vhdr);
    V_SANDI_tots{end+1} = V_SANDI_tot;

    Vhdr = spm_vol(img_path_fsoma);
    V_fsoma_tot = spm_read_vols(Vhdr);
    V_fsoma_tots{end+1} = V_fsoma_tot;


end

if strcmp(run,'run-01')
    V_rsoma_tots_vitality_run1=V_SANDI_tots;
else
    V_rsoma_tots_vitality_run2=V_SANDI_tots;
end

%% PRIN

V_Din_tots_prin={};
V_De_tots_prin={};
V_rsoma_tots_prin={};
V_fsoma_tots_prin={};
V_fextra_tots_prin={};
V_fneurite_tots_prin={};


n_subjs=12;

for i = 1:1:n_subjs

    if i < 10

        
        i=num2str(i);
        img_path_rsoma_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/dwi/SANDI_2MNI/sub-0',i,'_Rsoma_SANDI-fit_2MNI.nii.gz');
        img_path_fsoma_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/dwi/SANDI_2MNI/sub-0',i,'_fsoma_SANDI-fit_2MNI.nii.gz');
        img_path_De_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/dwi/SANDI_2MNI/sub-0',i,'_De_SANDI-fit_2MNI.nii.gz');
        img_path_Din_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/dwi/SANDI_2MNI/sub-0',i,'_Din_SANDI-fit_2MNI.nii.gz');
        img_path_fneurite_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/dwi/SANDI_2MNI/sub-0',i,'_fneurite_SANDI-fit_2MNI.nii.gz');
        img_path_fextra_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-0',i,'/dwi/SANDI_2MNI/sub-0',i,'_fextra_SANDI-fit_2MNI.nii.gz');

    else

        i=num2str(i);
        img_path_rsoma_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/dwi/SANDI_2MNI/sub-',i,'_Rsoma_SANDI-fit_2MNI.nii.gz');
        img_path_fsoma_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/dwi/SANDI_2MNI/sub-',i,'_fsoma_SANDI-fit_2MNI.nii.gz');
        img_path_De_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/dwi/SANDI_2MNI/sub-',i,'_De_SANDI-fit_2MNI.nii.gz');
        img_path_Din_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/dwi/SANDI_2MNI/sub-',i,'_Din_SANDI-fit_2MNI.nii.gz');
        img_path_fneurite_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/dwi/SANDI_2MNI/sub-',i,'_fneurite_SANDI-fit_2MNI.nii.gz');
        img_path_fextra_prin = strcat('/storage/shared/PRINAntonello2022/BIDS/derivatives/sub-',i,'/dwi/SANDI_2MNI/sub-',i,'_fextra_SANDI-fit_2MNI.nii.gz');

    end

    Vhdr = spm_vol(img_path_rsoma_prin);
    V_rsoma_tot = spm_read_vols(Vhdr);
    V_rsoma_tots_prin{end+1} = V_rsoma_tot;

    Vhdr = spm_vol(img_path_fsoma_prin);
    V_fsoma_tot = spm_read_vols(Vhdr);
    V_fsoma_tots_prin{end+1} = V_fsoma_tot;

    Vhdr = spm_vol(img_path_fneurite_prin);
    V_fneurite_tot = spm_read_vols(Vhdr);
    V_fneurite_tots_prin{end+1} = V_fneurite_tot;

    Vhdr = spm_vol(img_path_fextra_prin);
    V_fextra_tot = spm_read_vols(Vhdr);
    V_fextra_tots_prin{end+1} = V_fextra_tot;

    Vhdr = spm_vol(img_path_Din_prin);
    V_Din_tot = spm_read_vols(Vhdr);
    V_Din_tots_prin{end+1} = V_Din_tot;

    Vhdr = spm_vol(img_path_De_prin);
    V_De_tot = spm_read_vols(Vhdr);
    V_De_tots_prin{end+1} = V_De_tot;

end

%% load GM
img_path_GM='/storage/shared/Atlas/atlas_GM_on_MNI152_T1_2mm.nii.gz';
Vhdr = spm_vol(img_path_GM);
V_GM = spm_read_vols(Vhdr);
%%
%FOR EACH SUBJECT:
%MASK WITH GM
%MEAN OVER ALL GM PIXELS
%CHECK IF VITALITY MEAN DISTRIBUTION IS EQUAL TO PRIN MEAN DISTRIBUTION
means_vitality_run1=[];

for i = 1:n_subjs
    V_rsoma=V_rsoma_tots_vitality_run1{i};
    V_rsoma_masked=V_rsoma.*(V_GM>0.5);
    V_rsoma_masked(V_rsoma_masked==0)=NaN;
    means_vitality_run1(end+1)=nanmean(V_rsoma_masked(:));
end

means_prin_tot=[];
n_subjs=12;
for i = 1:n_subjs
    V_rsoma=V_rsoma_tots_prin{i};
    V_rsoma_masked=V_rsoma.*(V_GM>0.5);
    V_rsoma_masked(V_rsoma_masked==0)=NaN;
    means_prin_tot(end+1)=nanmean(V_rsoma_masked(:));
end


means_vitality_run2=[];
for i = 1:n_subjs
    V_rsoma=V_rsoma_tots_vitality_run2{i};
    V_rsoma_masked=V_rsoma.*(V_GM>0.5);
        V_rsoma_masked(V_rsoma_masked==0)=NaN;
    means_vitality_run2(end+1)=nanmean(V_rsoma_masked(:));
end

%%

[h,p]=ttest(means_prin_tot, means_vitality_run1_tot);
[r_corr,p_corr]=corrcoef(means_prin_tot, means_vitality_run1_tot,'rows','complete');

h_str=num2str(h);
p_str=num2str(round(p,2));

r_corr_str=num2str(round(r_corr(2),2));
p_corr_str=num2str(round(p_corr(2),2));

figure, 
h1=histogram(means_vitality_run1,BinWidth=0.3);
hold on
h2=histogram(means_prin,BinWidth=0.3);
txt = {strcat('t = ',h_str,', p-value:',p_str)};
text(11.4,5.5,txt, 'FontWeight', 'bold','FontSize',15);
txt = {strcat('r = ',r_corr_str,', p-value:',p_corr_str)};
text(11.4,6.5,txt, 'FontWeight', 'bold','FontSize',15);
ylim([0,8])
xlabel('Mean soma size in Grey Matter','FontWeight','bold','FontSize',15);
ylabel('Counts','FontWeight','bold','FontSize',15);
legend();
%%
means_prin=means_prin_tot;
means_prin(11)=[];
[h,p]=ttest(means_prin, means_vitality_run2);
[r_corr,p_corr]=corrcoef(means_prin, means_vitality_run2,'rows','complete');

h_str=num2str(h);
p_str=num2str(round(p,2));

r_corr_str=num2str(round(r_corr(2),2));
p_corr_str=num2str(round(p_corr(2),2));

figure, 
h1=histogram(means_vitality_run2,BinWidth=0.3);
hold on
h2=histogram(means_prin,BinWidth=0.3);
txt = {strcat('t = ',h_str,', p-value:',p_str)};
text(11,3,txt, 'FontWeight', 'bold','FontSize',15);
txt = {strcat('r = ',r_corr_str,', p-value:',p_corr_str)};
text(11,7,txt, 'FontWeight', 'bold','FontSize',15);
ylim([0,8])
xlabel('Mean soma size in Grey Matter','FontWeight','bold','FontSize',15);
ylabel('Counts','FontWeight','bold','FontSize',15);
legend();

%%

[h,p]=ttest(means_vitality_run1, means_vitality_run2);
[r_corr,p_corr]=corrcoef(means_vitality_run1, means_vitality_run2,'rows','complete');

h_str=num2str(h);
p_str=num2str(round(p,2));

r_corr_str=num2str(round(r_corr(2),2));
p_corr_str=num2str(round(p_corr(2),2));

figure, 
h1=histogram(means_vitality_run2,BinWidth=0.1,EdgeAlpha=0.7);
hold on
h2=histogram(means_vitality_run1,BinWidth=0.1,EdgeAlpha=0.7);
txt = {strcat('t = ',h_str,', p-value:',p_str)};
text(11.9,9.5,txt, 'FontWeight', 'bold','FontSize',12);
txt = {strcat('r = ',r_corr_str,', p-value:',p_corr_str)};
text(11.9,10,txt, 'FontWeight', 'bold','FontSize',12);
ylim([0,12])
xlabel('Mean soma size in Grey Matter','FontWeight','bold','FontSize',15);
ylabel('Counts','FontWeight','bold','FontSize',15);
legend();

bias12=nanmean(means_vitality_run1,2)-nanmean(means_vitality_run2,2);

%% check if regions with small rsoma and big rsoma correspond between acquisitions

%% check b0 intensities in raw Vitality image

% img_path=strcat('/media/nas_rete/Vitality/sub-001/dwi/sub-001_run-01_dwi.nii.gz');
% Vhdr = spm_vol(img_path);
% V_run1 = spm_read_vols(Vhdr);
% size(V_run1)
% 
% img_path=strcat('/media/nas_rete/Vitality/maps2MNI/SegSub01_t1_mp2rage/sub-001_run-01_PVE_0_on_b0.nii.gz');
% Vhdr = spm_vol(img_path);
% V_PVE_run1 = spm_read_vols(Vhdr);
% size(V_PVE_run1)
% 
% img_path=strcat('/media/nas_rete/Vitality/sub-001/dwi/sub-001_run-02_dwi.nii.gz');
% Vhdr = spm_vol(img_path);
% V_run2 = spm_read_vols(Vhdr);
% size(V_run2)
% 
% img_path=strcat('/media/nas_rete/Vitality/maps2MNI/SegSub01_t1_mp2rage/sub-001_run-02_PVE_0_on_b0.nii.gz');
% Vhdr = spm_vol(img_path);
% V_PVE_run2 = spm_read_vols(Vhdr);
% size(V_PVE_run2)

run='run-01';%CHANGE
%%%%%%%%%%%%%

subjects = importdata(strcat('/media/nas_rete/Vitality/code/subjs_DWI.txt'));

V_pve0_tots={};
V_dwi_tots={};



%run2
subjects([11,13,27])=[];

n_subjs=length(subjects);

for i = 1:1:n_subjs
    i

    subj=num2str(subjects{i});
    img_path_pve0 = strcat('/media/nas_rete/Vitality/maps2MNI/pve_on_b0/',subj,'_',run,'_PVE_0_on_b0.nii.gz');

    Vhdr = spm_vol(img_path_pve0);
    V_pve0_tot = spm_read_vols(Vhdr);
    V_pve0_tots{end+1} = V_pve0_tot;
end

%It's really slow (maybe it depends on images since the raw images are heavier)
for i = 1:1:n_subjs
    i
    img_path_dwi = strcat('/media/nas_rete/Vitality/',subj,'/dwi/',subj,'_',run,'_dwi.nii.gz');


    Vhdr = spm_vol(img_path_dwi);
    V_dwi_tot = spm_read_vols(Vhdr);
    V_dwi_tots{end+1} = V_dwi_tot;

end

if strcmp(run,'run-01')
    means_pve0_run1_subjs=[];
    for j=1:n_subjs
        j
        V_dwi_tots_run1=V_dwi_tots{j};
        V_pve0_tots_run1=V_pve0_tots{j};

        %create dataframe with volumes in one column and bvalues in an other column
        %tell to reorder rows of dataframe according to ordered values of second
        %column
        a=size(V_dwi_tots_run1);
        means_pve0_run1=[];
        for i = 1:a(4)
            V_pve0_tots_run1(V_pve0_tots_run1>0)=1;
            V_pve0_tots_run1(V_pve0_tots_run1<1)=0;
            V_masked_with_PVE=V_dwi_tots_run1(:,:,:,i).*V_pve0_tots_run1;
            mean_pve0=mean(V_masked_with_PVE,"all");
            means_pve0_run1(end+1)=mean_pve0;
        end
        means_pve0_run1_subjs(j,:)=means_pve0_run1;
    end

else
    means_pve0_run2_subjs=[];
    for j=1:n_subjs
        j
        V_dwi_tots_run2=V_dwi_tots{j};
        V_pve0_tots_run2=V_pve0_tots{j};

        b=size(V_dwi_tots_run2);
        means_pve0_run2=[];
        for i = 1:b(4)
            V_pve0_tots_run2(V_pve0_tots_run2>0)=1;
            V_pve0_tots_run2(V_pve0_tots_run2<1)=0;
            V_masked_with_PVE=V_dwi_tots_run2(:,:,:,i).*V_pve0_tots_run2;
            mean_pve0=mean(V_masked_with_PVE,"all");
            means_pve0_run2(end+1)=mean_pve0;
        end
        means_pve0_run2_subjs(j,:)=means_pve0_run2;
    end
end
%%


mean_pve0_all_subjs_run1=mean(means_pve0_run1_subjs,1);
mean_pve0_all_subjs_run2=mean(means_pve0_run2_subjs,1);


% figure, 
% scatter(1:length(means_pve0_run1), means_pve0_run1);
% hold on
% scatter(1:length(means_pve0_run2), means_pve0_run2);
% legend('run1','run2');

%assuming b-values are the same for all runs and subjs
sub001run01dwi=load('/media/nas_rete/Vitality/sub-001/dwi/sub-001_run-01_dwi.bval');
sub001run02dwi=load('/media/nas_rete/Vitality/sub-001/dwi/sub-001_run-02_dwi.bval');

bvalues=unique(sub001run01dwi);

vec_run1=cat(1,mean_pve0_all_subjs_run1, sub001run01dwi);
vec_run1=vec_run1';
mean_b0_run1=mean(vec_run1(vec_run1(:,2)==0,1));

vec_run2 = cat(1, mean_pve0_all_subjs_run2, sub001run02dwi);
vec_run2=vec_run2';
mean_b0_run2=mean(vec_run2(vec_run2(:,2)==0,1));



figure, 
s = scatter(vec_run1(:,2),log10(vec_run1(:,1)/mean_b0_run1),'o');
hold on
h = scatter(vec_run2(:,2), log10(vec_run2(:,1)/mean_b0_run2),'o');
legend("run1", "run2", "FontSize",15, 'FontWeight','bold');
xlabel("b-values", "FontSize",15, 'FontWeight','bold');
ylabel("S/S_0", "FontSize",15, 'FontWeight','bold');
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = 'b';
h.MarkerEdgeColor = 'r';
h.MarkerFaceColor = 'r';
%ylim([-1.2,-0.2])



figure,
for i = 1:length(bvalues)

    bval=bvalues(i);
    %select only b=0 values
    b0_run2=vec_run2(vec_run2(:,2)==bval);
    b0_run1=vec_run1(vec_run1(:,2)==bval);

    j=num2str(bval);

    subplot(3,3,i);
    s = scatter(1:length(b0_run1),b0_run1,'o');%/mean_b0_run1
    hold on
    h = scatter(1:length(b0_run2),b0_run2,'o');%/mean_b0_run2
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = 'b';
    h.MarkerEdgeColor = 'r';
    h.MarkerFaceColor = 'r';    
    xlabel("time point", "FontSize",15, 'FontWeight','bold');
    ylabel("CSF Signal", "FontSize",15, 'FontWeight','bold');
    title(strcat("b-val=",j));
end
lgd=legend("run1", "run2", "FontSize",15, 'FontWeight','bold');
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 15;



b0_run2=vec_run2(vec_run2(:,2)==0);
b0_run1=vec_run1(vec_run1(:,2)==0);
figure,
s = scatter(1:length(b0_run1),b0_run1/mean_b0_run1,'o');
hold on
h = scatter(1:length(b0_run2),b0_run2/mean_b0_run2,'o');
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = 'b';
h.MarkerEdgeColor = 'r';
h.MarkerFaceColor = 'r';    
xlabel("time point", "FontSize",15, 'FontWeight','bold');
ylabel("CSF Signal", "FontSize",15, 'FontWeight','bold');
title(strcat("b-val=0"));
legend("run1","run2");

%%
b0_run2=vec_run2(vec_run2(:,2)==0);
b0_run1=vec_run1(vec_run1(:,2)==0);
%create dataframe and color them.
b0_all_runs=cat(1, b0_run1/mean_b0_run1, b0_run2/mean_b0_run2);
figure,
s = scatter(1:15,b0_all_runs(1:15,1),'o');
hold on
h = scatter(16:30,b0_all_runs(16:30,1),'o');
% hold on
% yline(b0_all_runs(1),'--');
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = 'b';
h.MarkerEdgeColor = 'r';
h.MarkerFaceColor = 'r';    
xlabel("time point", "FontSize",15, 'FontWeight','bold');
ylabel("CSF Signal", "FontSize",15, 'FontWeight','bold');
title(strcat("b-val=0"));
legend("run1","run2");

b0_run2=vec_run2(vec_run2(:,2)==0);
b0_run1=vec_run1(vec_run1(:,2)==0);
%create dataframe and color them.
b0_all_runs=cat(1, b0_run1, b0_run2);
figure,
s = scatter(1:15,b0_all_runs(1:15,1),'o');
hold on
h = scatter(16:30,b0_all_runs(16:30,1),'o');
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = 'b';
h.MarkerEdgeColor = 'r';
h.MarkerFaceColor = 'r';    
xlabel("time point", "FontSize",15, 'FontWeight','bold');
ylabel("CSF Signal", "FontSize",15, 'FontWeight','bold');
title(strcat("b-val=0"));
legend("run1","run2");

%% Check if distributions are significantly different (Chieti vs Genova)

subj='sub-02';
outputpath=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Harmonization/SANDI/Comparison_Chieti_Genova/',subj);

GM_thr=0.5;
fsoma_thr=0.15;

if exist(outputpath)==0
mkdir (outputpath)
end



%import GM pve (80 slices)
img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/BIDS/derivatives/',subj,'/anat/',subj,'_t1_mp2rage_sag_p3_iso_fast_UNI_brain_pve_1_on_SANDI_Genova.nii.gz');%subj,'_t1_mp2rage_sag_p3_iso_fast_UNI_brain_seg_pve_1_res_on_SANDI.nii.gz'subj,'_t1_mp2rage_sag_p3_iso_fast_UNI_brain_pve_1_on_SANDI_resampled.nii.gz
Vhdr = spm_vol(img_path);
V_GM_tot = spm_read_vols(Vhdr);
size(V_GM_tot)
%import atlas (80 slices)
img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/BIDS/derivatives/',subj,'/dwiGenova/',subj,'_AAL3v1_2mm_on_Rsoma_SANDI-fit_Genova.nii.gz');
Vhdr = spm_vol(img_path);
V_atlas_tot = spm_read_vols(Vhdr);
size(V_atlas_tot)
%import rsoma map Chieti (80 slices)
img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/BIDS/derivatives/',subj,'/dwiGenova/SANDI_output/',subj,'_Rsoma_SANDI-fit.nii.gz');
Vhdr = spm_vol(img_path);
V_rsoma_tot = spm_read_vols(Vhdr);

img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/BIDS/derivatives/',subj,'/dwiGenova/SANDI_output/',subj,'_fsoma_SANDI-fit.nii.gz');
Vhdr = spm_vol(img_path);
V_fsoma_tot = spm_read_vols(Vhdr);

V_fsoma = V_fsoma_tot;
V_fsoma(V_fsoma>fsoma_thr)=1;
V_fsoma(V_fsoma<1)=0;

V_GM = V_GM_tot;
V_GM(V_GM>fsoma_thr)=1;
V_GM(V_GM<1)=0;

V_rsoma_GM = V_rsoma_tot.*V_GM.*V_fsoma; 

% figure, imagesc(V_rsoma_GM(:,:,45))
% figure, hist(V_rsoma_GM(:))
% figure, imagesc(V_atlas_tot(:,:,40))



regions=unique(V_atlas_tot(:));
means_Genova=[];
for k = 2:numel(regions)

    V_atlas = V_atlas_tot;

    for ii = 1:length(V_atlas_tot(:))
        if V_atlas_tot(ii)==regions(k)
            V_atlas(ii)=1;
        else
            V_atlas(ii)=0;
        end
    end

    V_rsoma_masked=V_rsoma_GM.*V_atlas;
    V_rsoma_masked(V_rsoma_masked==0)=NaN;

    mean_rsoma=nanmean(V_rsoma_masked(:));

    
    means_Genova(end+1)=mean_rsoma;
end

figure, 
h=histogram(means_chieti_AFNI,50);
h.FaceColor="#4DBEEE";
hold on
k=histogram(means_Genova,50);
k.FaceColor="#FF0000";
title('sub-02','FontWeight','bold','FontSize',15)
xlabel('Mean soma size in Grey Matter','FontWeight','bold','FontSize',15);
ylabel('Counts','FontWeight','bold','FontSize',15);
legend("Chieti2x2x2","Genova2x2x2");
%ylim([0,25]);



%import GM pve Genova (72 slices)
img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/BIDS/derivatives/',subj,'/anat/','T1_MP2RAGE_brain_pve_1_on_sub-Genova_De_SANDI-fit_Genova.nii.gz');
Vhdr = spm_vol(img_path);
V_GM_tot = spm_read_vols(Vhdr);

%import atlas Genova (72 slices)
img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/BIDS/derivatives/',subj,'/dwiGenova/',subj,'_AAL3v1_2mm_on_Rsoma_SANDI-fit_Genova.nii.gz');
Vhdr = spm_vol(img_path);
V_atlas_tot = spm_read_vols(Vhdr);

%import rsoma map Genova (72 slices)
img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/BIDS/derivatives/',subj,'/dwiGenova/SANDI_output/',subj,'_Rsoma_SANDI-fit.nii.gz');
Vhdr = spm_vol(img_path);
V_rsoma_tot = spm_read_vols(Vhdr);

img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/BIDS/derivatives/',subj,'/dwiGenova/SANDI_output/',subj,'_fsoma_SANDI-fit.nii.gz');
Vhdr = spm_vol(img_path);
V_fsoma_tot = spm_read_vols(Vhdr);

V_fsoma = V_fsoma_tot;
V_fsoma(V_fsoma>fsoma_thr)=1;
V_fsoma(V_fsoma<1)=0;

V_GM = V_GM_tot;
V_GM(V_GM>fsoma_thr)=1;
V_GM(V_GM<1)=0;

% DWI Genova ha 72 slices; DWI Chieti ha 80 slices
V_rsoma_GM = V_rsoma_tot.*V_GM;%.*V_fsoma;  % edited by AC, 13-11-24

means_Genova=[];
for k = 2:numel(regions)

    V_atlas = V_atlas_tot;

    for ii = 1:length(V_atlas_tot(:))
        if V_atlas_tot(ii)==regions(k)
            V_atlas(ii)=1;
        else
            V_atlas(ii)=0;
        end
    end

    V_rsoma_masked=V_rsoma_GM.*V_atlas;
    V_rsoma_masked(V_rsoma_masked==0)=NaN;

    mean_rsoma=nanmean(V_rsoma_masked(:));

    
    means_Genova(end+1)=mean_rsoma;
end

% figure
% hist(means_Genova,50)
% title('SANDI protocol Genova''FontWeight','bold','FontSize',15)
% xlabel('Mean soma size in Grey Matter','FontWeight','bold','FontSize',15);
% ylabel('Counts','FontWeight','bold','FontSize',15);

%isequal(means_chieti,means_Genova)

[h1,p1] = ttest2(means_chieti,means_Genova);

% Histogram with Rsoma mean in UniCh and UniGe
% h = figure;
figname = 'Rsoma_distribution';
figure(1),
histogram(means_chieti,50)
hold on
histogram(means_Genova,50)
%title('SANDI protocol Chieti ','FontWeight','bold','FontSize',15)
xlabel('Mean soma size in Grey Matter','FontWeight','bold','FontSize',15);
ylabel('Counts','FontWeight','bold','FontSize',15);
legend()
path_to_image=strcat(outputpath,subj,'_',figname,'_GM_',num2str(GM_thr),'_fsoma_',num2str(fsoma_thr),'.png');
saveas(figure(1), path_to_image);
%saveas(h,strcat(outputpath,figname,'_GM_',num2str(GM_thr),'_fsoma_',num2str(fsoma_thr)),'fig');

% Correlation plot
[r,p]=corrcoef(means_chieti,means_Genova,'rows','complete');
corr=num2str(round(r(2),2));
pvalue=num2str(round(p(2),2));
figname='Corr_plot';
figure(2), 
scatter(means_chieti,means_Genova)
hold on
plot(0:16,0:16)
txt = {strcat('r = ',corr,', p-value:',pvalue)};
text(8,4,txt, 'FontWeight', 'bold','FontSize',15);
xlabel('Means Chieti','FontWeight','bold');
ylabel('Means Genova','FontWeight','bold');
path_to_image=strcat(outputpath,subj,'_',figname,'_GM_',num2str(GM_thr),'_fsoma_',num2str(fsoma_thr),'.png');
saveas(figure(2), path_to_image);

%% Reproducibility test (SANDI maps)

subj='sub-03';

parameters={'Rsoma','fneurite','fsoma','fextra','De','Din'};
outputpath=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/',subj);

for par=1:6
    parameter=parameters{par};
    GM_thr=0.5;
    fsoma_thr=0.15;
    
    if exist(outputpath)==0
    mkdir (outputpath)
    end
    
    runs=1:3;
    
    
    % I PART: load data
    
    V_GM_tots={};
    V_atlas_tots={};
    V_rsoma_tots={};
    V_fsoma_tots={};
    
    for run = runs

        if strcmp(subj,'sub-03') && run==3
            run=2;
        end


        i=num2str(run);    
        img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/derivatives/',subj,'/dwi/run-0',i,'/',subj,'_run-0',i,'_desc-UNI_MP2RAGE_brain_pve_1_on_SANDI.nii.gz');
        Vhdr = spm_vol(img_path);
        V_GM_tot = spm_read_vols(Vhdr);
        V_GM_tots{end+1} = V_GM_tot;
   
    
    
    
    
        i=num2str(run);    
        img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/derivatives/',subj,'/dwi/run-0',i,'/',subj,'_run-0',i,'_AAL3v1_2mm_on_res_De_SANDI-fit.nii.gz');
        Vhdr = spm_vol(img_path);
        V_atlas_tot = spm_read_vols(Vhdr);
        V_atlas_tots{end+1} = V_atlas_tot;

    
    
    
    
        i=num2str(run);    
        img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/derivatives/',subj,'/dwi/run-0',i,'/SANDI_output/',subj,'_run-0',i,'_',parameter,'_SANDI-fit.nii.gz');
        Vhdr = spm_vol(img_path);
        V_rsoma_tot = spm_read_vols(Vhdr);
        V_rsoma_tots{end+1} = V_rsoma_tot;
   
    
    
    
   
        i=num2str(run);    
        img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/derivatives/',subj,'/dwi/run-0',i,'/SANDI_output/',subj,'_run-0',i,'_fsoma_SANDI-fit.nii.gz');
        Vhdr = spm_vol(img_path);
        V_fsoma_tot = spm_read_vols(Vhdr);
        V_fsoma_tots{end+1} = V_fsoma_tot;
    end
    
    % II PART: check difference in labels between atlases in subjects spaces 
    run1=V_atlas_tots{1};
    run2=V_atlas_tots{2};
    run3=V_atlas_tots{3};
    
    n_labels_run1=unique(run1);
    n_labels_run2=unique(run2);
    n_labels_run3=unique(run3);
    
    diff_run23=setxor(n_labels_run2,n_labels_run3);
    diff_run12=setxor(n_labels_run1,n_labels_run2);
    diff_run13=setxor(n_labels_run1,n_labels_run3);
    diff_runs=unique([diff_run23,diff_run12,diff_run13]);
    
    run1_diff=run1;
    for i = 1:length(run1(:))
        if run1(i)==diff_runs
            run1_diff(i)=0;
        end
    end
    
    run2_diff=run2;
    for i = 1:length(run2(:))
        if run2(i)==diff_runs
            run2_diff(i)=0;
        end
    end
    
    run3_diff=run3;
    for i = 1:length(run3(:))
        if run3(i)==diff_runs
            run3_diff(i)=0;
        end
    end
    
    V_atlas_diff_tots={run1_diff, run2_diff, run3_diff}; %new list of atlases
    %test
    % size(unique(run3))
    % III PART: COMPUTE MEANS FOR EACH RUN
    
    
    
    
    means_runs=[];
    
    for i = runs
        V_fsoma = V_fsoma_tots{i};
        V_fsoma(V_fsoma>fsoma_thr)=1;
        V_fsoma(V_fsoma<1)=0;
    
        V_GM = V_GM_tots{i};
        V_GM(V_GM>fsoma_thr)=1;
        V_GM(V_GM<1)=0;
    
        V_rsoma_tot = V_rsoma_tots{i};
        V_atlas_tot = V_atlas_diff_tots{i};
    
        V_rsoma_GM = V_rsoma_tot.*V_GM.*V_fsoma;
    
    
    
        regions=unique(V_atlas_tot(:));
    
        means_run=[];
    
        for k = 2:numel(regions)
    
            V_atlas = V_atlas_tot;
    
            for ii = 1:length(V_atlas_tot(:))
                if V_atlas_tot(ii)==regions(k)
                    V_atlas(ii)=1;
                else
                    V_atlas(ii)=0;
                end
            end
    
            V_rsoma_masked=V_rsoma_GM.*V_atlas;
            V_rsoma_masked(V_rsoma_masked==0)=NaN;
    
            mean_rsoma=nanmean(V_rsoma_masked(:));
    
    
            means_run(end+1)=mean_rsoma;
        end
    
        means_runs(i,:)=means_run;
    
    end
    
    %IV PART: HISTOGRAM
    
    %A two tailed hypothesis test
    [h12,p12,ci,stats]=ttest(means_runs(1,:),means_runs(2,:));
    [h23,p23,ci,stats]=ttest(means_runs(2,:),means_runs(3,:));
    [h13,p13,ci,stats]=ttest(means_runs(1,:),means_runs(3,:));
    
    h=[h12, h23, h13];
    p=[p12, p23, p13];
    
    figure, 
    s=histogram(means_runs(1,:),50);
    s.FaceColor="g";
    hold on
    k=histogram(means_runs(2,:),50);
    k.FaceColor="b";
    hold on
    k=histogram(means_runs(3,:),50);
    k.FaceColor="r";
    title(subj,'FontWeight','bold','FontSize',15)
    xlabel('Mean soma size in Grey Matter','FontWeight','bold','FontSize',15);
    ylabel('Counts','FontWeight','bold','FontSize',15);
    legend("run-01","run-02","run-03");
    h12=num2str(h(1));
    h23=num2str(h(2));
    h13=num2str(h(3));
    p12=num2str(round(p(1),5));
    p23=num2str(round(p(2),5));
    p13=num2str(round(p(3),5));
    txt12 = {strcat('t_{12} = ',h12,', p-value_{12}:',p12)};
    txt23 = {strcat('t_{23} = ',h23,', p-value_{23}:',p23)};
    txt13 = {strcat('t_{13} = ',h13,', p-value_{13}:',p13)};
    text(2.5,16,txt12, 'FontWeight', 'bold','FontSize',15);
    text(2.5,14,txt23, 'FontWeight', 'bold','FontSize',15);
    text(2.5,12,txt13, 'FontWeight', 'bold','FontSize',15);
    %ylim([0,25]);
    
    
    
    
    
    % V PART: matrix for ICC analysis
    
    % k = strfind(subj,'-');
    % subj(k)='';
    % 
    % struct.description = strcat('Matrix for ICC analysis',' ',subj);
    % struct = setfield(struct,[strcat('means_runs_',subj)], means_runs);
    folder = strcat(strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/',subj));
    
    if exist([folder])==0
        mkdir(folder)
    end
    
    cd(folder)
    
    save(strcat(subj,'_means_runs_',parameter,'.mat'),'means_runs');
end

bias12=nanmean(means_runs(1,:),2)-nanmean(means_runs(2,:),2);
bias23=nanmean(means_runs(2,:),2)-nanmean(means_runs(3,:),2);
bias13=nanmean(means_runs(1,:),2)-nanmean(means_runs(3,:),2);

biases=[bias12, bias23, bias13];

f = fopen(strcat('biases_',subj,'.txt'), 'w');
fprintf(f, 'bias12 bias23 bias13\n');
fprintf(f, '\n');
writematrix(biases, strcat('biases_',subj,'_',parameter,'.txt'), 'WriteMode', 'append');
fclose(f);

%% for many subjs

% % I PART: load data
% subjects = importdata(strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/code/subjs_dwi.txt'));
% subjects(3)=[];%remove sub-03 as it has no run3
% n_subjs=length(subjects);
% runs=1:3;
% 
% V_GM_tots_run1={};
% V_GM_tots_run2={};
% V_GM_tots_run3={};
% 
% V_atlas_tots_run1={};
% V_atlas_tots_run2={};
% V_atlas_tots_run3={};
% 
% V_rsoma_tots_run1={};
% V_rsoma_tots_run2={};
% V_rsoma_tots_run3={};
% 
% V_fsoma_tots_run1={};
% V_fsoma_tots_run2={};
% V_fsoma_tots_run3={};
% 
% for run = runs
%     i=num2str(run);
%     for j = 1:n_subjs
%     
%     subj=subjects{j};
% 
%     img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/derivatives/',subj,'/dwi/run-0',i,'/',subj,'_run-0',i,'_desc-UNI_MP2RAGE_brain_pve_1_on_SANDI.nii.gz');
%     Vhdr = spm_vol(img_path);
%     V_GM_tot = spm_read_vols(Vhdr);
% 
%     img_path_atlas=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/derivatives/',subj,'/dwi/run-0',i,'/',subj,'_run-0',i,'_AAL3v1_2mm_on_res_De_SANDI-fit.nii.gz');
%     Vhdr_atlas = spm_vol(img_path_atlas);
%     V_atlas_tot = spm_read_vols(Vhdr_atlas);
% 
%     img_path_rsoma=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/derivatives/',subj,'/dwi/run-0',i,'/SANDI_output/',subj,'_run-0',i,'_Rsoma_SANDI-fit.nii.gz');
%     Vhdr_rsoma = spm_vol(img_path_rsoma);
%     V_rsoma_tot = spm_read_vols(Vhdr_rsoma);
% 
%     img_path_fsoma=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/derivatives/',subj,'/dwi/run-0',i,'/SANDI_output/',subj,'_run-0',i,'_fsoma_SANDI-fit.nii.gz');
%     Vhdr_fsoma = spm_vol(img_path_fsoma);
%     V_fsoma_tot = spm_read_vols(Vhdr_fsoma);
% 
%     if run==1
%         V_GM_tots_run1{end+1} = V_GM_tot;
%         V_atlas_tots_run1{end+1} = V_atlas_tot;
%         V_rsoma_tots_run1{end+1} = V_rsoma_tot;
%         V_fsoma_tots_run1{end+1} = V_fsoma_tot;
%     elseif run==2
%         V_GM_tots_run2{end+1} = V_GM_tot;
%         V_atlas_tots_run2{end+1} = V_atlas_tot;
%         V_rsoma_tots_run2{end+1} = V_rsoma_tot;
%         V_fsoma_tots_run2{end+1} = V_fsoma_tot;
%     elseif run==3
%         V_GM_tots_run3{end+1} = V_GM_tot;
%         V_atlas_tots_run3{end+1} = V_atlas_tot;
%         V_rsoma_tots_run3{end+1} = V_rsoma_tot;
%         V_fsoma_tots_run3{end+1} = V_fsoma_tot;
%     end
% 
%     end
%     
% end

% II PART: check difference in labels between atlases in subjects spaces 
%note: for En_DWI_Analysis I didn't have problems for labels lost from one
%subj to an other one.

% diff between runs and between subjs (like sub-01 run-01 and sub-02
% run-01...)
% III PART: COMPUTE MEANS FOR EACH RUN
%IV PART: HISTOGRAM


%% raw data

subj='sub-01';
runs=1:3;
%load data 
PVE_thr=0.5;
pve_n=0;

if pve_n==0
    tissue='CSF';
elseif pve_n==1
    tissue='GM';
elseif pve_n==2
    tissue='WM';
else
    disp('Wrong number of pve. Stop the execution of the code :) !');
    pause(inf);
end

str_pve_n=num2str(pve_n);


V_pve_tots={};

for run = runs
    i=num2str(run);    
    img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/derivatives/',subj,'/dwi/run-0',i,'/',subj,'_run-0',i,'_desc-UNI_MP2RAGE_brain_pve_',str_pve_n,'_on_SANDI.nii.gz');
    Vhdr = spm_vol(img_path);
    V_pve_tot = spm_read_vols(Vhdr);
    V_pve_tots{end+1} = V_pve_tot;
end

V_raw_tots={};

for run = runs
    i=num2str(run);    
    img_path=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/',subj,'/dwi/',subj,'_run-0',i,'_dwi.nii.gz');
    Vhdr = spm_vol(img_path);
    V_raw_tot = spm_read_vols(Vhdr);
    V_raw_tots{end+1} = V_raw_tot;
end

noise_runs=[];
means_pve_runs=[];
bkgs_pve_runs=[];

%I can normalize data. What about standard errors and background noise.

for i = runs

    V_PVE_run=V_pve_tots{i};
    V_raw=V_raw_tots{i};
    disp(strcat('processing run-0',num2str(i)));
    a=size(V_raw);

    means_pve_run=[];
    SEs_pve_run=[];
    bkgs_pve_run=[];


    for j = 1:a(4)

        disp(strcat('processing run-0',num2str(i),' volume', num2str(j)));

        V_PVE_run(V_PVE_run>PVE_thr)=1;
        V_PVE_run(V_PVE_run<1)=0;
        V_masked_with_PVE=V_raw(:,:,:,j).*V_PVE_run;
        mean_pve0=mean(V_masked_with_PVE,"all");
        SE_pve0=std(V_masked_with_PVE(:))/sqrt(numel(V_masked_with_PVE(:)));

        %select 4 ROIS for noise
        r1=V_raw(1:10,1:10,:,j);
        r2=V_raw(1:10,100:110,:,j);
        r3=V_raw(100:110,100:110,:,j);
        r4=V_raw(100:110,1:10,:,j);
        std_r1=std(r1(:));
        std_r2=std(r2(:));
        std_r3=std(r3(:));
        std_r4=std(r4(:));
        stds=[std_r1,std_r2,std_r3,std_r4];
        noise=mean(stds);
 
        bkgs_pve_run(end+1)=noise;
        means_pve_run(end+1)=mean_pve0;
        SEs_pve_run(end+1)=SE_pve0;

    end
    
    means_pve_runs(i,:)=means_pve_run;
    noise_runs(i,:)=SEs_pve_run;
    bkgs_pve_runs(i,:)=bkgs_pve_run;
end



% plots

PVE_thr_str=num2str(PVE_thr);

for run = runs
    i=num2str(run);
    load(strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/BIDS/',subj,'/dwi/',subj,'_run-0',i,'_dwi.bval'));
end

bvals=sub_01_run_01_dwi;%b-vals shuffled always in the same way



vec_run1=cat(1,means_pve_runs(1,:), bvals, noise_runs(1,:),bkgs_pve_runs(1,:));
vec_run1=vec_run1';
mean_b0_run1=mean(vec_run1(vec_run1(:,2)==0,1));

vec_run2=cat(1,means_pve_runs(2,:), bvals, noise_runs(2,:),bkgs_pve_runs(2,:));
vec_run2=vec_run2';
mean_b0_run2=mean(vec_run2(vec_run2(:,2)==0,1));

vec_run3=cat(1,means_pve_runs(3,:), bvals, noise_runs(3,:),bkgs_pve_runs(3,:));
vec_run3=vec_run3';
mean_b0_run3=mean(vec_run3(vec_run3(:,2)==0,1));

if pve_n==0
    vec_run1_CSF=vec_run1;
    vec_run2_CSF=vec_run2;
    vec_run3_CSF=vec_run3;

    b0_run1_CSF=mean_b0_run1;
    b0_run2_CSF=mean_b0_run2;
    b0_run3_CSF=mean_b0_run3;

elseif pve_n==1

    vec_run1_GM=vec_run1;
    vec_run2_GM=vec_run2;
    vec_run3_GM=vec_run3;

    b0_run1_GM=mean_b0_run1;
    b0_run2_GM=mean_b0_run2;
    b0_run3_GM=mean_b0_run3;

elseif pve_n==2

    vec_run1_WM=vec_run1;
    vec_run2_WM=vec_run2;
    vec_run3_WM=vec_run3;

    b0_run1_WM=mean_b0_run1;
    b0_run2_WM=mean_b0_run2;
    b0_run3_WM=mean_b0_run3;
end

%%
n_tissues=3; 
%if you set n_tissue=3, 
% then you need to run 
% the first part for all tissues 
% (in order to have
%the matrices for all tissues)

% plots (calculate rate of decay?)
outputpath=strcat('/media/nas_rete/PRIN2022PNRR_Tomass/Repeatability/',subj,'/Normalized/',tissue,'/');
if exist([outputpath])==0
mkdir ([outputpath])
end

if n_tissues==3
 
    %non linear fit using least square method
%     abc=2000;
%     fun = @(abc,t) exp(-t/abc);
%     abc0=0.1;
%     abcLSQ = lsqcurvefit( fun, abc0, vec_run1_CSF(:,2), vec_run1_CSF(:,1)/max(vec_run1_CSF(:,1)));

    %Linear fit
%     P_run1_CSF = polyfit(vec_run1_CSF(:,2),log10(vec_run1_CSF(:,1)/max(vec_run1_CSF(:,1))),1);
%     yfit_run1_CSF = P_run1_CSF(1)*vec_run1_CSF(:,2)+P_run1_CSF(2);
% 
%     P_run2_CSF = polyfit(vec_run2_CSF(:,2),log10(vec_run2_CSF(:,1)/max(vec_run2_CSF(:,1))),1);
%     yfit_run2_CSF = P_run2_CSF(1)*vec_run2_CSF(:,2)+P_run2_CSF(2);
% 
%     P_run3_CSF = polyfit(vec_run3_CSF(:,2),log10(vec_run3_CSF(:,1)/max(vec_run3_CSF(:,1))),1);
%     yfit_run3_CSF = P_run3_CSF(1)*vec_run3_CSF(:,2)+P_run3_CSF(2);
% 
%     P_run1_GM = polyfit(vec_run1_GM(:,2),log10(vec_run1_GM(:,1)/max(vec_run1_GM(:,1))),1);
%     yfit_run1_GM = P_run1_GM(1)*vec_run1_GM(:,2)+P_run1_GM(2);
% 
%     P_run2_GM = polyfit(vec_run2_GM(:,2),log10(vec_run2_GM(:,1)/max(vec_run2_GM(:,1))),1);
%     yfit_run2_GM = P_run2_GM(1)*vec_run2_GM(:,2)+P_run2_GM(2);
% 
%     P_run3_GM = polyfit(vec_run3_GM(:,2),log10(vec_run3_GM(:,1)/max(vec_run3_GM(:,1))),1);
%     yfit_run3_GM = P_run3_GM(1)*vec_run3_GM(:,2)+P_run3_GM(2);
% 
%     P_run1_WM = polyfit(vec_run1_WM(:,2),log10(vec_run1_WM(:,1)/max(vec_run1_WM(:,1))),1);
%     yfit_run1_WM = P_run1_WM(1)*vec_run1_WM(:,2)+P_run1_WM(2);
% 
%     P_run2_WM = polyfit(vec_run2_WM(:,2),log10(vec_run2_WM(:,1)/max(vec_run2_WM(:,1))),1);
%     yfit_run2_WM = P_run2_WM(1)*vec_run2_WM(:,2)+P_run2_WM(2);
% 
%     P_run3_WM = polyfit(vec_run3_WM(:,2),log10(vec_run3_WM(:,1)/max(vec_run3_WM(:,1))),1);
%     yfit_run3_WM = P_run3_WM(1)*vec_run3_WM(:,2)+P_run3_WM(2);



    figure(1),
    s = scatter(vec_run1_CSF(:,2), log10(vec_run1_CSF(:,1)/b0_run1_CSF),'o');%/max(vec_run1_CSF(:,1))
    hold on
%    plot(vec_run1_CSF(:,2), log10(vec_run1_CSF(:,1)/max(vec_run1_CSF(:,1))),'-','LineWidth',3,'Color',"b");
%     hold on
    h = scatter(vec_run2_CSF(:,2),log10(vec_run2_CSF(:,1)/b0_run2_CSF),'o');
    hold on
    %plot(vec_run2_CSF(:,2), yfit_run2_CSF,'--','LineWidth',3,'Color',"r");
%     hold on
    z = scatter(vec_run3_CSF(:,2),log10(vec_run3_CSF(:,1)/b0_run3_CSF),'o');
    hold on
    %plot(vec_run3_CSF(:,2), yfit_run3_CSF,'--','LineWidth',3,'Color',"#EDB120");
%     hold on
    l = scatter(vec_run1_GM(:,2), log10(vec_run1_GM(:,1)/b0_run1_GM),'o');
    hold on
    %plot(vec_run1_GM(:,2), yfit_run1_GM,'--','LineWidth',3,'Color',"b");
%     hold on
    k = scatter(vec_run2_GM(:,2),log10(vec_run2_GM(:,1)/b0_run2_GM),'o');
    hold on
    %plot(vec_run2_GM(:,2), yfit_run2_GM,'--','LineWidth',3,'Color',"r");
%     hold on
    j = scatter(vec_run3_GM(:,2),log10(vec_run3_GM(:,1)/b0_run3_GM),'o');
    hold on
    %plot(vec_run3_GM(:,2), yfit_run3_GM,'--','LineWidth',3,'Color',"#EDB120");
%     hold on
    q = scatter(vec_run1_WM(:,2), log10(vec_run1_WM(:,1)/b0_run1_WM),'o');
    hold on
    %plot(vec_run1_WM(:,2), yfit_run1_WM,'--','LineWidth',3,'Color',"b");
%     hold on
    f = scatter(vec_run2_WM(:,2),log10(vec_run2_WM(:,1)/b0_run2_WM),'o');
    hold on
    %plot(vec_run2_WM(:,2), yfit_run2_WM,'--','LineWidth',3,'Color',"r");
%     hold on
    g = scatter(vec_run3_WM(:,2),log10(vec_run3_WM(:,1)/b0_run3_WM),'o');
%     hold on
    %plot(vec_run3_WM(:,2), yfit_run3_WM,'--','LineWidth',3,'Color',"#EDB120");
    legend("run1", "run2", "run3", "FontSize",15, 'FontWeight','bold');
    xlabel("b-values", "FontSize",15, 'FontWeight','bold');
    ylabel("S/S_0", "FontSize",15, 'FontWeight','bold');
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = 'b';
    h.MarkerEdgeColor = 'r';
    h.MarkerFaceColor = 'r';
    z.MarkerEdgeColor = "#EDB120";
    z.MarkerFaceColor = "#EDB120";
    %
    l.MarkerEdgeColor = 'b';
    l.MarkerFaceColor = 'b';
    k.MarkerEdgeColor = 'r';
    k.MarkerFaceColor = 'r';
    j.MarkerEdgeColor = "#EDB120";
    j.MarkerFaceColor = "#EDB120";
    %
    q.MarkerEdgeColor = 'b';
    q.MarkerFaceColor = 'b';
    f.MarkerEdgeColor = 'r';
    f.MarkerFaceColor = 'r';
    g.MarkerEdgeColor = "#EDB120";
    g.MarkerFaceColor = "#EDB120";
    
    str='WM';
    dim=[0.5000 0.5100 .3 .3];
    annotation('textbox',dim,'String',str,'FontWeight','bold',"FontSize",15);
    str='GM';
    dim=[0.5000 0.5100 .3 .3];
    annotation('textbox',dim,'String',str,'FontWeight','bold',"FontSize",15);
    str='CSF';
    dim=[0.5000 0.5100 .3 .3];
    annotation('textbox',dim,'String',str,'FontWeight','bold',"FontSize",15);
    %text(5000,-0.6,'WM','FontWeight','bold');
    %text(5000,-0.8,'GM','FontWeight','bold');
    %text(5000,-1.1,'CSF','FontWeight','bold');
    sgtitle(strcat('sub',subj));

else

    figure(1),
    s = scatter(vec_run1(:,2), log10(vec_run1(:,1)/max(vec_run1(:,1))),'o');
    hold on
    h = scatter(vec_run2(:,2),log10(vec_run2(:,1)/max(vec_run2(:,1))),'o');
    hold on
    z = scatter(vec_run3(:,2),log10(vec_run3(:,1)/max(vec_run3(:,1))),'o');
    legend("run1", "run2", "run3", "FontSize",15, 'FontWeight','bold');
    xlabel("b-values", "FontSize",15, 'FontWeight','bold');
    ylabel("S/S_0", "FontSize",15, 'FontWeight','bold');
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = 'b';
    h.MarkerEdgeColor = 'r';
    h.MarkerFaceColor = 'r';
    z.MarkerEdgeColor = "#EDB120";
    z.MarkerFaceColor = "#EDB120";
    sgtitle(strcat(tissue,' signal (', subj,')'));
    path_to_image=strcat(outputpath,strcat(subj,'_',tissue,'_runs_fordifferentbvalues_normalized.jpg'));
    saveas(figure(1), path_to_image);
end

%%
bvalues=unique(bvals);

figure(2),
for i = 1:length(bvalues)

    bval=bvalues(i);
    %select means and noise values for each b-val
%     mean_run3=vec_run3(vec_run3(:,2)==bval);
%     mean_run2=vec_run2(vec_run2(:,2)==bval);
%     mean_run1=vec_run1(vec_run1(:,2)==bval);

    idx_run1=find(vec_run1(:,2)==bval);
    idx_run2=find(vec_run2(:,2)==bval);
    idx_run3=find(vec_run3(:,2)==bval);

    raws_run1=vec_run1(idx_run1,:);
    raws_run2=vec_run2(idx_run2,:);
    raws_run3=vec_run3(idx_run3,:);

    bkg_run1=mean(vec_run1(:,4));
    bkg_run2=mean(vec_run2(:,4));
    bkg_run3=mean(vec_run3(:,4));

    j=num2str(bval);

    subplot(3,3,i);
    %s = scatter(1:length(b0_run1),b0_run1,'o');
    s = errorbar(1:length(raws_run1(:,1)), raws_run1(:,1), raws_run1(:,3), raws_run1(:,3), 'o');
    hold on
    yline(bkg_run1,'--','LineWidth',1,'Color','b');
    hold on
    %h = scatter(1:length(b_run2),b_run2,'o');/max(vec_run1(:,1))
    h = errorbar(1:length(raws_run2(:,1)), raws_run2(:,1), raws_run2(:,3), raws_run2(:,3), 'o');
    hold on
    yline(bkg_run2,'--','LineWidth',1,'Color','r');
    hold on
    %z = scatter(1:length(b_run3),b_run3,'o');
    z = errorbar(1:length(raws_run3(:,1)), raws_run3(:,1), raws_run3(:,3), raws_run3(:,3), 'o'); 
    hold on
    yline(bkg_run3,'--','LineWidth',1,'Color','g');
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = 'b';
    h.MarkerEdgeColor = 'r';
    h.MarkerFaceColor = 'r'; 
    z.MarkerEdgeColor = "#EDB120";
    z.MarkerFaceColor = "#EDB120";  
    xlabel("time point", "FontSize",15, 'FontWeight','bold');
    ylabel(strcat(tissue," Signal"), "FontSize",15, 'FontWeight','bold');
    title(strcat("b-val=",j));
   
end
lgd=legend("S\_run1", "N\_run1", "S\_run2", "N\_run2", "S\_run3", "N\_run3", "FontSize",15, 'FontWeight','bold','Location','southeast');
set(lgd,'position',[.60 .1 .1 .1])
set(figure(2), 'Units', 'Normalized', 'Outerposition', [0 0 0.7 0.9]);
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 15;
sgtitle(strcat(subj,' PVE thr=',PVE_thr_str),'FontSize',14,'FontWeight','bold');
path_to_image=strcat(outputpath,strcat(subj,'_',tissue,'signal_subplot_PVE',PVE_thr_str,'_errors.jpg'));
saveas(figure(2), path_to_image);




figure(3),
for i = 1:length(bvalues)

    bval=bvalues(i);
    %select means and noise values for each b-val
%     mean_run3=vec_run3(vec_run3(:,2)==bval);
%     mean_run2=vec_run2(vec_run2(:,2)==bval);
%     mean_run1=vec_run1(vec_run1(:,2)==bval);

    idx_run1=find(vec_run1(:,2)==bval);%_WM
    idx_run2=find(vec_run2(:,2)==bval);
    idx_run3=find(vec_run3(:,2)==bval);

    %not normalized
%     raws_run1=vec_run1(idx_run1,:);
%     raws_run2=vec_run2(idx_run2,:);
%     raws_run3=vec_run3(idx_run3,:);

    %normalized
    raws_run1_norm=vec_run1(idx_run1,:);%/mean_b0_run1;%/max(vec_run1(idx_run1,:));
    raws_run2_norm=vec_run2(idx_run2,:);%/mean_b0_run2;%/max(vec_run2(idx_run2,:));
    raws_run3_norm=vec_run3(idx_run3,:);%/mean_b0_run3;%/max(vec_run3(idx_run3,:));

    j=num2str(bval);

    subplot(3,3,i);
    s = scatter(1:length(raws_run1_norm(:,1)),raws_run1_norm(:,1),'o');%add log10
    hold on
    h = scatter(1:length(raws_run2_norm(:,1)),raws_run2_norm(:,1),'o');
    hold on
    z = scatter(1:length(raws_run3_norm(:,1)),raws_run3_norm(:,1),'o');    
    s.MarkerEdgeColor = 'b';
    s.MarkerFaceColor = 'b';
    h.MarkerEdgeColor = 'r';
    h.MarkerFaceColor = 'r'; 
    z.MarkerEdgeColor = "#EDB120";
    z.MarkerFaceColor = "#EDB120";  
    xlabel("time point", "FontSize",15, 'FontWeight','bold');
    ylabel(strcat(tissue," Signal"), "FontSize",15, 'FontWeight','bold');
    title(strcat("b-val=",j));
   
end
lgd=legend("run1", "run2", "run3", 'FontWeight','bold','FontSize',15);
set(lgd,'position',[.60 .1 .1 .1])
set(figure(3), 'Units', 'Normalized', 'Outerposition', [0 0 0.7 0.9]);
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 15;
sgtitle(strcat(subj,' PVE thr=',PVE_thr_str),'FontSize',15,'FontWeight','bold');
path_to_image=strcat(outputpath,strcat(subj,'_',tissue,'signal_subplot_PVE',PVE_thr_str,'.jpg'));
saveas(figure(3), path_to_image);




b0_run3=vec_run3(vec_run3(:,2)==0);
b0_run2=vec_run2(vec_run2(:,2)==0);
b0_run1=vec_run1(vec_run1(:,2)==0);
figure(4),
s = scatter(1:length(b0_run1),b0_run1,'o');
hold on
h = scatter(1:length(b0_run2),b0_run2,'o');
hold on
z = scatter(1:length(b0_run3),b0_run3,'o');
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = 'b';
h.MarkerEdgeColor = 'r';
h.MarkerFaceColor = 'r';    
z.MarkerEdgeColor = "#EDB120";
z.MarkerFaceColor = "#EDB120";  
xlabel("time point", "FontSize",15, 'FontWeight','bold');
ylabel(strcat(tissue," Signal"), "FontSize",15, 'FontWeight','bold');
title(strcat(subj,' b-val=0',' PVE thr=',PVE_thr_str));
legend("run1","run2", "run3",'Location','northwest');
path_to_image=strcat(outputpath,strcat(subj,'_',tissue,'_bval0_PVE',PVE_thr_str,'.jpg'));
saveas(figure(4), path_to_image);

% 
% b0_run3=vec_run3(vec_run3(:,2)==0);
% b0_run2=vec_run2(vec_run2(:,2)==0);
% b0_run1=vec_run1(vec_run1(:,2)==0);

bval=0;
idx_run1=find(vec_run1(:,2)==bval);
idx_run2=find(vec_run2(:,2)==bval);
idx_run3=find(vec_run3(:,2)==bval);
raws_run1=vec_run1(idx_run1,:);
raws_run2=vec_run2(idx_run2,:);
raws_run3=vec_run3(idx_run3,:);

bkg_run1=mean(vec_run1(:,4));
bkg_run2=mean(vec_run2(:,4));
bkg_run3=mean(vec_run3(:,4));

figure(5),
s = scatter(1:length(raws_run1(:,1)),raws_run1(:,1),'o');
%s = errorbar(1:length(raws_run1(:,1)), raws_run1(:,1), raws_run1(:,3), raws_run1(:,3), 'o');
hold on
%yline(bkg_run1,'--','LineWidth',1,'Color','b');
h = scatter(16:30,raws_run2(:,1),'o');
%hold on
%h = errorbar(1:length(raws_run2(:,1)), raws_run2(:,1), raws_run2(:,3), raws_run2(:,3), 'o');
% hold on
% yline(bkg_run2,'--','LineWidth',1,'Color','r');
z = scatter(31:45,raws_run3(:,1),'o');
hold on
%z = errorbar(1:length(raws_run3(:,1)), raws_run3(:,1), raws_run3(:,3), raws_run3(:,3), 'o');
% hold on
% yline(bkg_run2,'--','LineWidth',1,'Color','y');
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = 'b';
h.MarkerEdgeColor = 'r';
h.MarkerFaceColor = 'r';    
z.MarkerEdgeColor = "#EDB120";
z.MarkerFaceColor = "#EDB120";  
xlabel("time point", "FontSize",15, 'FontWeight','bold');
ylabel(strcat(tissue, " Signal"), "FontSize",15, 'FontWeight','bold');
title(strcat(subj,' b-val=0',' PVE thr=',PVE_thr_str));
%legend("S\_run1", "N\_run1", "S\_run2", "N\_run2", "S\_run3", "N\_run3");
legend("run1", "run2", "run3",'Location','northwest');
path_to_image=strcat(outputpath,strcat(subj,'_',tissue,'_bval0_PVE',PVE_thr_str,'_asfunctionoftime.jpg'));
saveas(figure(5), path_to_image);

figure(6),
%s = scatter(1:length(raws_run1(:,1)),raws_run1(:,1),'o');
s = errorbar(1:length(raws_run1(:,1)), raws_run1(:,1), raws_run1(:,3), raws_run1(:,3), 'o');
hold on
%yline(bkg_run1,'--','LineWidth',1,'Color','b');
%h = scatter(16:30,raws_run2(:,1),'o');
%hold on
h = errorbar(1:length(raws_run2(:,1)), raws_run2(:,1), raws_run2(:,3), raws_run2(:,3), 'o');
% hold on
% yline(bkg_run2,'--','LineWidth',1,'Color','r');
%z = scatter(31:45,raws_run3(:,1),'o');
hold on
z = errorbar(1:length(raws_run3(:,1)), raws_run3(:,1), raws_run3(:,3), raws_run3(:,3), 'o');
% hold on
% yline(bkg_run2,'--','LineWidth',1,'Color','y');
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = 'b';
h.MarkerEdgeColor = 'r';
h.MarkerFaceColor = 'r';    
z.MarkerEdgeColor = "#EDB120";
z.MarkerFaceColor = "#EDB120";  
xlabel("time point", "FontSize",15, 'FontWeight','bold');
ylabel(strcat(tissue, " Signal"), "FontSize",15, 'FontWeight','bold');
title(strcat(subj,' b-val=0',' PVE thr=',PVE_thr_str));
%legend("S\_run1", "N\_run1", "S\_run2", "N\_run2", "S\_run3", "N\_run3");
legend("run1", "run2", "run3",'Location','northwest');
path_to_image=strcat(outputpath,strcat(subj,'_',tissue,'_bval0_PVE',PVE_thr_str,'_asfunctionoftime_errors.jpg'));
saveas(figure(6), path_to_image);



