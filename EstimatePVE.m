%Program to estimate partial volume effect in images with Partial Fourier
%and without Partial Fourier

% img_path_pve1_woPF = '/storage/shared/SANDI_240418/NIFTI/pve1_to_DWIwoPF.nii.gz';
% img_path_pve0_woPF = '/storage/shared/SANDI_240418/NIFTI/pve0_to_DWIwoPF.nii.gz';
% img_path_pve2_woPF = '/storage/shared/SANDI_240418/NIFTI/pve2_to_DWIwoPF.nii.gz';
% img_path_map_woPF = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240418tr3000woPF/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_Rsoma.nii.gz';
% title_plot = 'Rsoma map without Partial Fourier Acq240418';

% img_path_pve1_woPF = '/media/nas_rete/Work_manuela/DICOM_PRIN2022Antonello/NIFTI/SANDI_analysis/InterpNN_pve1_to_b0.nii.gz';
% img_path_pve0_woPF = '/media/nas_rete/Work_manuela/DICOM_PRIN2022Antonello/NIFTI/SANDI_analysis/InterpNN_pve0_to_b0.nii.gz';
% img_path_pve2_woPF = '/media/nas_rete/Work_manuela/DICOM_PRIN2022Antonello/NIFTI/SANDI_analysis/InterpNN_pve2_to_b0.nii.gz';
% img_path_map_woPF = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240516/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_Rsoma.nii.gz';
% title_plot = 'Rsoma map without Partial Fourier Acq240516';


img_path_pve1_woPF = '/media/nas_rete/Work_manuela/PVE240523/pve1_to_b0.nii.gz';
img_path_pve0_woPF = '/media/nas_rete/Work_manuela/PVE240523/pve0_to_b0.nii.gz';
img_path_pve2_woPF = '/media/nas_rete/Work_manuela/PVE240523/pve2_to_b0.nii.gz';

img_path_map_woPF = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240523/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_Rsoma.nii.gz';
%img_path_map_woPF = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240523/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_fsoma.nii.gz';

title_plot = 'Rsoma map without Partial Fourier Acq240523';
y_label = 'Mean value of Rsoma';

Vhdr_map_woPF = spm_vol(img_path_map_woPF);
V_map_woPF = spm_read_vols(Vhdr_map_woPF);

%Vhdr_map_wPF = spm_vol(img_path_map_wPF);
%V_map_wPF = spm_read_vols(Vhdr_map_wPF);

%V_maps = {V_map_woPF V_map_wPF};

thr_pve = 0.1:0.1:0.9;

Vhdr_pve1_woPF = spm_vol(img_path_pve1_woPF);
pve1_woPF = spm_read_vols(Vhdr_pve1_woPF);
Vhdr_pve0_woPF = spm_vol(img_path_pve0_woPF);
pve0_woPF = spm_read_vols(Vhdr_pve0_woPF);
Vhdr_pve2_woPF = spm_vol(img_path_pve2_woPF);
pve2_woPF = spm_read_vols(Vhdr_pve2_woPF);

% err_list_woPF = [];
% err_list_wPF = [];

mean_par_woPF = [];

% mean_par_wPF = [];


%mask with fsoma
img_path_map = '/storage/shared/SANDI_240229_MatlabToolbox/SANDI-Matlab-Toolbox-Latest-Release-main/dataset240523/SANDI_MainFolder/derivatives/SANDI_analysis/sub-01/ses-01/SANDI_Output/SANDI-fit_fsoma.nii.gz';
Vhdr_map = spm_vol(img_path_map);
V_mapfsoma = spm_read_vols(Vhdr_map);

%figure, imagesc(V_map(:,:,40));
% figure, histogram(V_map);

V_mapfsoma=V_mapfsoma>0.4;
%deep grey matter putamen region
% figure, imagesc(V_mapfsoma(:,:,40));

V_mapRsoma = spm_read_vols(Vhdr_map_woPF);
Vmasked = V_mapRsoma.*V_mapfsoma;
V_maps = {Vmasked};

for element = 1:numel(V_maps)

    for value = 1:1:numel(thr_pve)

        if element == 1
            pve0 = pve0_woPF;
            pve1 = pve1_woPF;
            pve2 = pve2_woPF;
        else
            pve0 = pve0_wPF;
            pve1 = pve1_wPF;
            pve2 = pve2_wPF;
        end

        %make binary mask
        GM=pve1>thr_pve(value);
        WM=pve2>thr_pve(value);
        CSF=pve0>thr_pve(value);
    
        V_GM = GM.*V_maps{element};
        V_WM = WM.*V_maps{element};
        V_CSF = CSF.*V_maps{element};

        %figure,imagesc(V_map_masked(:,:,150));
    
%         calculate mean of parameter value in one region of brain cortex
%         in a second moment, choose a precise region (like the temporal one)
        
%         if value == 0.1 || value == 0.5 || value == 0.9
%         
%             figure,
%             imagesc(V_WM(:,:,70))
%         
%         
%         
%         
%         end
    
        V_GM(V_GM==0)=NaN;
        V_WM(V_WM==0)=NaN;
        V_CSF(V_CSF==0)=NaN;
    
        V_GM = rmmissing(V_GM(:));
        V_WM = rmmissing(V_WM(:));
        V_CSF = rmmissing(V_CSF(:));

        Vmean(1) = mean(V_GM);
        Vmean(2) = mean(V_WM);
        Vmean(3) = mean(V_CSF);

%         if element == 1
%             mean_par_woPF(end+1,:) = Vmean;
%         else
%             mean_par_wPF(end+1,:) = Vmean;
%         end

        % Confidence Intervals using bootstrapping
        V_GM = bootstrp(100,@mean,V_GM);
        V_WM = bootstrp(100,@mean,V_WM);
        V_CSF = bootstrp(100,@mean,V_CSF);

    

        
        SEM_GM = std(V_GM);
                  
        SEM_WM = std(V_WM);
    
        SEM_CSF = std(V_CSF);
        

    err_GM(value,1) = SEM_GM;
    err_WM(value,1) = SEM_WM;
    err_CSF(value,1) = SEM_CSF;

        if element == 1
            err_list_woPF(end+1,:) = CI;
        else 
            err_list_wPF(end+1,:) = CI;
        end 
        
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

% Figure with SEM as errorbars
figure, 
p1 = errorbar(thr_pve, mean_par_woPF(:,1), err_GM);
hold on
p2 = errorbar(thr_pve, mean_par_woPF(:,2), err_WM);
hold on
p3 = errorbar(thr_pve, mean_par_woPF(:,3), err_CSF);
p1.LineStyle = '-';
p1.Marker = '.';
p2.LineStyle = '-';
p2.Marker = '.';
p3.LineStyle = '-';
p3.Marker = '.';
%ylim([10 13.5]);
%ylim([10 13.5]);
legend([p1 p2 p3],{'Grey Matter', 'White Matter', 'CSF'});
xlabel('Threashold value');
%ylabel(ylabel);
ylabel(y_label);
title(title_plot);








%% Slope
%deltayslope_woPF = mean_par_woPF(7,1) - mean_par_woPF(1,1); 0.9128
%deltayslope_wPF = mean_par_wPF(7,1) - mean_par_wPF(1,1); %0.8684



%%
% With Partial Fourier and Without Partial Fourier 

%{
mean_par = {mean_par_wPF mean_par_woPF};
err_list = {err_list_wPF err_list_woPF};
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


%GM, CSF and WM without Partial Fourier
figure,
p1 = errorbar(thr_pve, mean_par_woPF(:,1), mean_par_woPF(:,1)-err_list_woPF(:,1), err_list_woPF(:,2)-mean_par_woPF(:,1));
hold on
p2 = errorbar(thr_pve, mean_par_woPF(:,2), mean_par_woPF(:,2)-err_list_woPF(:,3), err_list_woPF(:,4)-mean_par_woPF(:,2));
hold on
p3 = errorbar(thr_pve, mean_par_woPF(:,3), mean_par_woPF(:,3)-err_list_woPF(:,5), err_list_woPF(:,6)-mean_par_woPF(:,3));
p1.LineStyle = '-';
p1.Marker = '.';
p2.LineStyle = '-';
p2.Marker = '.';
p3.LineStyle = '-';
p3.Marker = '.';
%ylim([10 13.5]);
legend([p1 p2 p3],{'Grey Matter', 'White Matter', 'CSF'});
xlabel('Threashold value');
ylabel(y_label);


%GM with and without Partial Fourier
figure, 
p1 = errorbar(thr_pve, mean_par_wPF(:,1), mean_par_wPF(:,1)-err_list_wPF(:,1), err_list_wPF(:,2)-mean_par_wPF(:,1));
hold on
p2 = errorbar(thr_pve, mean_par_woPF(:,1), mean_par_woPF(:,1)-err_list_woPF(:,1), err_list_woPF(:,2)-mean_par_woPF(:,1));
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
p1 = errorbar(thr_pve, mean_par_wPF(:,2), mean_par_wPF(:,2)-err_list_wPF(:,3), err_list_wPF(:,4)-mean_par_wPF(:,2));
hold on
p2 = errorbar(thr_pve, mean_par_woPF(:,2), mean_par_woPF(:,2)-err_list_woPF(:,3), err_list_woPF(:,4)-mean_par_woPF(:,2));
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
p1 = errorbar(thr_pve, mean_par_wPF(:,3), mean_par_wPF(:,3)-err_list_wPF(:,5), err_list_wPF(:,6)-mean_par_wPF(:,3));
hold on
p2 = errorbar(thr_pve, mean_par_woPF(:,3), mean_par_woPF(:,3)-err_list_woPF(:,5), err_list_woPF(:,6)-mean_par_woPF(:,3));
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
hist(err_list_wPF(:,2)-err_list_wPF(:,1));
title('Confidence intervals GM distribution with PF');

figure,
hist(err_list_woPF(:,2)-err_list_woPF(:,1));
title('Confidence intervals GM distribution without PF');
%}


