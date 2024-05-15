%Program to estimate partial volume effect in images with Partial Fourier
%and without Partial Fourier


%T1 to diffusion weighted image
img_path_pve1_wPF = '/storage/shared/SANDI_240418/NIFTI/pve1_to_DWIwPF.nii.gz';
img_path_pve0_wPF = '/storage/shared/SANDI_240418/NIFTI/pve0_to_DWIwPF.nii.gz';
img_path_pve2_wPF = '/storage/shared/SANDI_240418/NIFTI/pve2_to_DWIwPF.nii.gz';
img_path_map_wPF = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240418tr3000wPF/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_Rsoma.nii.gz';

img_path_pve1_woPF = '/storage/shared/SANDI_240418/NIFTI/pve1_to_DWIwoPF.nii.gz';
img_path_pve0_woPF = '/storage/shared/SANDI_240418/NIFTI/pve0_to_DWIwoPF.nii.gz';
img_path_pve2_woPF = '/storage/shared/SANDI_240418/NIFTI/pve2_to_DWIwoPF.nii.gz';
img_path_map_woPF = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240418tr3000woPF/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_Rsoma.nii.gz';

%img_path_map = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240418tr3000wPF/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_fsoma.nii.gz';
title_plot = 'Rsoma map with and without Partial Fourier';
y_label = 'Mean value of Rsoma';

Vhdr_map_woPF = spm_vol(img_path_map_woPF);
V_map_woPF = spm_read_vols(Vhdr_map_woPF);

Vhdr_map_wPF = spm_vol(img_path_map_wPF);
V_map_wPF = spm_read_vols(Vhdr_map_wPF);

V_maps = {V_map_woPF V_map_wPF};

thr_pve = 0.1:0.1:0.9;

Vhdr_pve1_wPF = spm_vol(img_path_pve1_wPF);
pve1_wPF = spm_read_vols(Vhdr_pve1_wPF);
Vhdr_pve0_wPF = spm_vol(img_path_pve0_wPF);
pve0_wPF = spm_read_vols(Vhdr_pve0_wPF);
Vhdr_pve2_wPF = spm_vol(img_path_pve2_wPF);
pve2_wPF = spm_read_vols(Vhdr_pve2_wPF);

Vhdr_pve1_woPF = spm_vol(img_path_pve1_woPF);
pve1_woPF = spm_read_vols(Vhdr_pve1_woPF);
Vhdr_pve0_woPF = spm_vol(img_path_pve0_woPF);
pve0_woPF = spm_read_vols(Vhdr_pve0_woPF);
Vhdr_pve2_woPF = spm_vol(img_path_pve2_woPF);
pve2_woPF = spm_read_vols(Vhdr_pve2_woPF);

Quantiles_list_wPF = [];
Quantiles_list_woPF = [];

CI_list_woPF=[];
mean_par_woPF = [];
CI_list_wPF=[];
mean_par_wPF = [];
V_map = V_map_woPF;
%metti in una funzione che data la Vmap restituisce 
%mean_par ed i corrispettivi CI_list
for element = 1:numel(V_maps)

    for value = 1:1:numel(thr_pve)
    
    
        %
        if element == 1
            pve0 = pve0_woPF;
            pve1 = pve1_woPF;
            pve2 = pve2_woPF;
        else
            pve0 = pve0_wPF;
            pve1 = pve1_wPF;
            pve2 = pve2_wPF;
        end

    
        %figure,imagesc(V_pve(:,:,150))
    
        %make binary mask
        GM=pve1>thr_pve(value);
        WM=pve2>thr_pve(value);
        CSF=pve0>thr_pve(value);
    
        %figure,imagesc(V_pve(:,:,150))
    
        %
    
        V_GM = GM.*V_maps{element};
        V_WM = WM.*V_maps{element};
        V_CSF = CSF.*V_maps{element};
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
        if element == 1
            mean_par_woPF(end+1,:) = Vmean;
        else
            mean_par_wPF(end+1,:) = Vmean;
        end
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

%         Quantiles = zeros(1,9);
%         Q_GM = quantile(V_GM,[0.25, 0.5, 0.75]);
%         Quantiles(1:1,1:3) = Q_GM;
%         Q_WM = quantile(V_WM,[0.25, 0.5, 0.75]);
%         Quantiles(1:1,4:6) = Q_WM;
%         Q_CSF = quantile(V_CSF,[0.25, 0.5, 0.75]);
%         Quantiles(1:1,7:9) = Q_CSF;
%     
%         if element == 1
%             Quantiles_list_woPF(end+1,:) = Quantiles;
%         else
%             Quantiles_list_wPF(end+1,:) = Quantiles;
%         end
    
        CI = zeros(1,6);
        
        SEM_GM = std(V_GM);% Standard Error /sqrt(length(V_GM))
        ts = tinv([0.025  0.975],length(V_GM)-1);% T-Score
               
        CI(1:1,1:2) = mean(V_GM) + ts*SEM_GM;
        
        
        
            
        SEM_WM = std(V_WM);% Standard Error /sqrt(length(V_WM))
        ts = tinv([0.025  0.975],length(V_WM)-1);% T-Score
        
        CI(1:1,3:4) = mean(V_WM) + ts*SEM_WM;
    
        SEM_CSF = std(V_CSF);% Standard Error /sqrt(length(V_CSF))
        ts = tinv([0.025  0.975],length(V_CSF)-1);% T-Score
        
        CI(1:1,5:6) = mean(V_CSF) + ts*SEM_CSF;
        
        if element == 1
            CI_list_woPF(end+1,:) = CI;
        else 
            CI_list_wPF(end+1,:) = CI;
        end 
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
end






%togli la divisione per il numero di campioni ! No però a quel punto la
%media diventa il nuovo campione quindi occorre comunque dividere 
%per il numero di campioni.

mean_par = {mean_par_wPF mean_par_woPF};
CI_list = {CI_list_wPF CI_list_woPF};
title_plots = {'Rsoma map with Partial Fourier' 'Rsoma map without Partial Fourier'};

for index = 1:1:numel(mean_par)

    figure,
    p1 = errorbar(thr_pve, mean_par{index}(:,1), mean_par{index}(:,1)-CI_list{index}(:,1), CI_list{index}(:,2)-mean_par{index}(:,1));
    hold on
    p2 = errorbar(thr_pve, mean_par{index}(:,2), mean_par{index}(:,2)-CI_list{index}(:,3), CI_list{index}(:,4)-mean_par{index}(:,2));
    hold on
    p3 = errorbar(thr_pve, mean_par{index}(:,3), mean_par{index}(:,3)-CI_list{index}(:,5), CI_list{index}(:,6)-mean_par{index}(:,3));
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
    title(title_plots{index});

end



%GM with and without Partial Fourier
figure, 
p1 = errorbar(thr_pve, mean_par_wPF(:,1), mean_par_wPF(:,1)-CI_list_wPF(:,1), CI_list_wPF(:,2)-mean_par_wPF(:,1));
hold on
p2 = errorbar(thr_pve, mean_par_woPF(:,1), mean_par_woPF(:,1)-CI_list_woPF(:,1), CI_list_woPF(:,2)-mean_par_woPF(:,1));
p1.LineStyle = '-';
p1.Marker = '.';
p2.LineStyle = '-';
p2.Marker = '.';
legend([p1 p2],{'GM with PF', 'GM without PF'});
xlabel('Threashold value');
ylabel(y_label);
ylim([10 13.5])
title(title_plot);



%WM with and without Partial Fourier
figure, 
p1 = errorbar(thr_pve, mean_par_wPF(:,2), mean_par_wPF(:,2)-CI_list_wPF(:,3), CI_list_wPF(:,4)-mean_par_wPF(:,2));
hold on
p2 = errorbar(thr_pve, mean_par_woPF(:,2), mean_par_woPF(:,2)-CI_list_woPF(:,3), CI_list_woPF(:,4)-mean_par_woPF(:,2));
p1.LineStyle = '-';
p1.Marker = '.';
p2.LineStyle = '-';
p2.Marker = '.';
legend([p1 p2],{'WM with PF', 'WM without PF'});
xlabel('Threashold value');
ylabel(y_label);
ylim([10 13.5])
title(title_plot);

%CSF with and without Partial Fourier
figure, 
p1 = errorbar(thr_pve, mean_par_wPF(:,3), mean_par_wPF(:,3)-CI_list_wPF(:,5), CI_list_wPF(:,6)-mean_par_wPF(:,3));
hold on
p2 = errorbar(thr_pve, mean_par_woPF(:,3), mean_par_woPF(:,3)-CI_list_woPF(:,5), CI_list_woPF(:,6)-mean_par_woPF(:,3));
p1.LineStyle = '-';
p1.Marker = '.';
p2.LineStyle = '-';
p2.Marker = '.';
legend([p1 p2],{'CSF with PF', 'CSF without PF'});
xlabel('Threashold value');
ylabel(y_label);
ylim([10 13.5])
title(title_plot);

%distribution error with and without Partial Fourier
figure, 
hist(CI_list_wPF(:,2)-CI_list_wPF(:,1));
title('Confidence intervals GM distribution with PF');

figure,
hist(CI_list_woPF(:,2)-CI_list_woPF(:,1));
title('Confidence intervals GM distribution without PF');

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

% figure, 
% p1 = errorbar(thr_pve, Quantiles_list_wPF(:,2), Quantiles_list_wPF(:,2) - Quantiles_list_wPF(:,1), Quantiles_list_wPF(:,2) - Quantiles_list_wPF(:,3));
% hold on
% p2 = errorbar(thr_pve, Quantiles_list_wPF(:,5), Quantiles_list_wPF(:,5) - Quantiles_list_wPF(:,4), Quantiles_list_wPF(:,6) - Quantiles_list_wPF(:,5));
% hold on
% p3 = errorbar(thr_pve, Quantiles_list_wPF(:,8), Quantiles_list_wPF(:,8) - Quantiles_list_wPF(:,7), Quantiles_list_wPF(:,9) - Quantiles_list_wPF(:,8));
% p1.LineStyle = '-';
% p1.Marker = '.';
% p2.LineStyle = '-';
% p2.Marker = '.';
% p3.LineStyle = '-';
% p3.Marker = '.';
% ylim([0 15]);
% ylim([10 13.5]);
% legend([p1 p2 p3],{'Grey Matter', 'White Matter', 'CSF'});
% xlabel('Threashold value');
% ylabel(y_label);
% title(title_plot);



%cambia tipo di interpolazione come suggerito da Davide.

%% Slope
deltayslope_woPF = mean_par_woPF(7,1) - mean_par_woPF(1,1); 0.9128
deltayslope_wPF = mean_par_wPF(7,1) - mean_par_wPF(1,1); %0.8684
