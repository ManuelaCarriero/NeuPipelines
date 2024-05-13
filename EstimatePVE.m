%Program to estimate partial volume effect in images with Partial Fourier
%and without Partial Fourier

%Insert the filename if you are in the same folder of the data.

%Rsoma map to T1
% img_path_pve1 = '/storage/shared/SANDI_240418/NIFTI/t1_mp2rage_sag_p3_iso_fast/Segmented/t1_mp2rage_UNIMaskedwithINV2_pve_1.nii.gz';
% img_path_pve0 = '/storage/shared/SANDI_240418/NIFTI/t1_mp2rage_sag_p3_iso_fast/Segmented/t1_mp2rage_UNIMaskedwithINV2_pve_0.nii.gz';
% img_path_pve2 = '/storage/shared/SANDI_240418/NIFTI/t1_mp2rage_sag_p3_iso_fast/Segmented/t1_mp2rage_UNIMaskedwithINV2_pve_2.nii.gz';
% img_path_map = '/storage/shared/SANDI_240418/NIFTI/RsomawPF_to_T1.nii.gz';
% title_plot = 'Rsoma map With Partial Fourier';

%T1 to diffusion weighted image
img_path_pve1 = '/storage/shared/SANDI_240418/NIFTI/pve1_to_DWIwPF.nii.gz';
img_path_pve0 = '/storage/shared/SANDI_240418/NIFTI/pve0_to_DWIwPF.nii.gz';
img_path_pve2 = '/storage/shared/SANDI_240418/NIFTI/pve2_to_DWIwPF.nii.gz';
img_path_map = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240418tr3000woPF/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_Rsoma.nii.gz';
%img_path_map = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240418tr3000wPF/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_fsoma.nii.gz';
title_plot = 'Rsoma map Without Partial Fourier';
y_label = 'Mean value of Rsoma';

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

Quantiles_list = [];

CI_list=[];
index = 1:1:6;


for value = 1:1:numel(thr_pve)


    %


    %figure,imagesc(V_pve(:,:,150))

    %make binary mask
    GM=pve1>thr_pve(value);
    WM=pve2>thr_pve(value);
    CSF=pve0>thr_pve(value);

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

    V_GM = rmmissing(V_GM(:));
    V_WM = rmmissing(V_WM(:));
    V_CSF = rmmissing(V_CSF(:));
% mean(V_GM,'omitnan');
    Vmean(1) = mean(V_GM);%Compute the arithmetic mean along the specified axis, ignoring NaNs.
    Vmean(2) = mean(V_WM);
    Vmean(3) = mean(V_CSF);
    mean_par(end+1,:) = Vmean;

    %Quantiles
%     Quantiles = zeros(1,9);
%     Q_GM = quantile(V_GM,[0.25, 0.5, 0.75]);
%     Quantiles(1:1,1:3) = Q_GM;
%     Q_WM = quantile(V_WM,[0.25, 0.5, 0.75]);
%     Quantiles(1:1,4:6) = Q_WM;
%     Q_CSF = quantile(V_CSF,[0.25, 0.5, 0.75]);
%     Quantiles(1:1,7:9) = Q_CSF;
% 
%     Quantiles_list(end+1,:) = Quantiles;

% 
% %     % Confidence Intervals using bootstrapping
    V_GM = bootstrp(100,@mean,V_GM);
    V_WM = bootstrp(100,@mean,V_WM);
    V_CSF = bootstrp(100,@mean,V_CSF);

    CI = zeros(1,6);
    
    SEM_GM = std(V_GM)/sqrt(length(V_GM));% Standard Error
    ts = tinv([0.025  0.975],length(V_GM)-1);% T-Score
           
    CI(1:1,1:2) = mean(V_GM) + ts*SEM_GM;
    
    
    
        
    SEM_WM = std(V_WM)/sqrt(length(V_WM));% Standard Error
    ts = tinv([0.025  0.975],length(V_WM)-1);% T-Score
    
    CI(1:1,3:4) = mean(V_WM) + ts*SEM_WM;

    SEM_CSF = std(V_CSF)/sqrt(length(V_CSF));% Standard Error
    ts = tinv([0.025  0.975],length(V_CSF)-1);% T-Score
    
    CI(1:1,5:6) = mean(V_CSF) + ts*SEM_CSF;
    
    CI_list(end+1,:) = CI;
% err_GM(value,1) = SEM_GM;
% err_WM(value,1) = SEM_WM;
% err_CSF(value,1) = SEM_CSF;
    
     %Check for distribution
%      figure,
%      subplot(1,3,1)
%      hist(V_GM)
%      subplot(1,3,2)
%      hist(V_WM)
%      subplot(1,3,3)
%      hist(V_CSF)

%     figure,
%     subplot(1,3,1)
%     qqplot(V_GM)
%     subplot(1,3,2)
%     qqplot(V_WM)
%     subplot(1,3,3)
%     qqplot(V_CSF)
end

figure, plot(thr_pve, mean_par);
%hold on
%plot(thr_pve, mean_par_wo);
xlabel('Threashold value');
%ylabel(y_label);
ylabel(y_label);
title(title_plot);
ylim([10 13.5]);










figure, 
p1 = errorbar(thr_pve, mean_par(:,1), mean_par(:,1)-CI_list(:,1), CI_list(:,2)-mean_par(:,1));
hold on
p2 = errorbar(thr_pve, mean_par(:,2), mean_par(:,2)-CI_list(:,3), CI_list(:,4)-mean_par(:,2));
hold on
p3 = errorbar(thr_pve, mean_par(:,3), mean_par(:,3)-CI_list(:,5), CI_list(:,6)-mean_par(:,3));
p1.LineStyle = '-';
p1.Marker = '.';
p2.LineStyle = '-';
p2.Marker = '.';
p3.LineStyle = '-';
p3.Marker = '.';
ylim([10 13.5]);
legend([p1 p2 p3],{'Grey Matter', 'White Matter', 'CSF'});
xlabel('Threashold value');
ylabel(y_label);
title(title_plot);

%% Figure with SEM as errorbars
% figure, 
% p1 = errorbar(thr_pve, mean_par(:,1), err_GM);
% hold on
% p2 = errorbar(thr_pve, mean_par(:,2), err_WM);
% hold on
% p3 = errorbar(thr_pve, mean_par(:,3), err_CSF);
% p1.LineStyle = '-';
% p1.Marker = '.';
% p2.LineStyle = '-';
% p2.Marker = '.';
% p3.LineStyle = '-';
% p3.Marker = '.';
% %ylim([10 13.5]);
% %ylim([10 13.5]);
% legend([p1 p2 p3],{'Grey Matter', 'White Matter', 'CSF'});
% xlabel('Threashold value');
% %ylabel(ylabel);
% ylabel(ylabel);
% title(title_plot);

%% Figure with quantiles as errorbars
%{
figure, 
p1 = errorbar(thr_pve, Quantiles_list(:,2), Quantiles_list(:,2) - Quantiles_list(:,1), Quantiles_list(:,2) - Quantiles_list(:,3));
% hold on
% p2 = errorbar(thr_pve, Quantiles_list(:,5), Quantiles_list(:,5) - Quantiles_list(:,4), Quantiles_list(:,6) - Quantiles_list(:,5));
% hold on
% p3 = errorbar(thr_pve, Quantiles_list(:,8), Quantiles_list(:,8) - Quantiles_list(:,7), Quantiles_list(:,9) - Quantiles_list(:,8));
p1.LineStyle = '-';
p1.Marker = '.';
% p2.LineStyle = '-';
% p2.Marker = '.';
% p3.LineStyle = '-';
% p3.Marker = '.';
%ylim([0 15]);
%ylim([10 13.5]);
%legend([p1 p2 p3],{'Grey Matter', 'White Matter', 'CSF'});
xlabel('Threashold value');
ylabel(y_label);
title(title_plot);
%}


%cambia tipo di interpolazione come suggerito da Davide.

%% Slope
%deltayslope_woPF = mean_par(7,1) - mean_par(1,1); 0.9128
deltayslope_wPF = mean_par(7,1) - mean_par(1,1); %0.8684
