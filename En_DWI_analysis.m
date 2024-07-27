%Loading SANDI map to compare with vascular and metabolic map
%Choose a map

parameter='Rsoma';%'Rsoma';'fsoma';'De'; insert parameter name compatible with SANDI files.
unit_of_measure = '(\mu m)';%'(\mu m^{2}/ms)';%'(\mu m)';(%);
label = strcat(parameter,unit_of_measure);
img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_SANDI = spm_read_vols(Vhdr);



%Loading GM pve map 

GM_path='/media/nas_rete/Vitality/clustering/SegSub01/GM2atlas.nii.gz';
Vhdr = spm_vol(GM_path);
V_GM_tot = spm_read_vols(Vhdr);
% figure, imagesc(V_GM(:,:,45))
% title('GM pve')

%Loading CBF map
img_path = '/media/nas_rete/Vitality/clustering/perf/sub-001_run-01_bh_CBF0_2MNI.nii.gz';
Vhdr = spm_vol(img_path);
V_CBF_tot = spm_read_vols(Vhdr);
%figure, imagesc(V_CBF(:,:,70))
V_CBF_tot=max(V_CBF_tot,0);%replace negative numbers with zeros


threasholds = 0.1:0.1:1;
index = 1:1:numel(threasholds);

figure(1), 
for i = index
    
    V_GM = V_GM_tot;

    V_GM(V_GM>=threasholds(i))=1;
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
    ylabel(strcat(label,' Signal'));
end
sgtitle("Rsoma vs CBF increasing GM PVE threashold (thr)");
set(figure(1), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(1), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/rsoma_vs_cbf_GMmask.png');




figure(2), 
p = plot(threasholds, corr_values, '-o', 'LineWidth', 1, 'MarkerFaceColor','b', 'MarkerSize',3);
xlabel('PVE threashold');
ylabel('r');
title('corr\_coef Rsoma (masked with GM) vs CBF');
set(figure(2), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(2), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/corr_coef_rsomavscbf_GMmask.png');


%Plot PVE maps threasholded
for i = index    
    V_GM = V_GM_tot;

    V_GM(V_GM>=threasholds(i))=1;
    V_GM(V_GM<threasholds(i))=0;
    % figure, imagesc(V_GM(:,:,45))
    % title('GM mask')
    
    PVE = V_GM.*V_GM_tot;
    
    slices_PVE{i} = PVE(:,:,45);
    
    %Mask SANDI map with grey matter
    
    V_SANDI_masked=V_SANDI.*V_GM;
    
    slices_SANDI{i} = V_SANDI_masked(:,:,45);
    

end



   
fig=figure(3);
for i = index
    subplot(4,4,i)
    imagesc(rot90(slices_SANDI{i}))
    thr = num2str(threasholds(i));
    title(strcat("thr>=",thr))
    axis equal
    axis off
end
h = axes(fig,'visible','off'); 
sgtitle("Rsoma map masked with GM varying thr")
maxi=max(V_SANDI_masked, [],'all');
dy=maxi/10;
labels = round(0:dy:maxi);
c = colorbar(h,'Position',[0.93 0.4 0.019 0.4], 'TickLabels',[0 2 3 5 6 8 9 11 12 14 16]);  % attach colorbar to h
% c.Limits=[0 maxi];
set(figure(3), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(3), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/rsoma_maps_threasholdingPVE.png');

%c = colorbar(h,'Location','southoutside','Position',[0.5 0.02 0.022 0.7]);
%c = colorbar(h,'Location','southoutside');
%set(c,'Position',[0.5 0.168 0.9 0.002])


%colormap(c,'jet')



%mask Rsoma map with fsoma
img_path='/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_fsoma_2MNI.nii.gz';
Vhdr = spm_vol(img_path);
V_fsoma = spm_read_vols(Vhdr);

figure(4), 
for i = index
    
    V_GM = V_GM_tot;

    V_GM(V_GM>0)=1;
    V_GM(V_GM<=0)=0;

    
    V_SANDI_masked=V_SANDI.*V_GM;

    V_fsoma_tot = V_fsoma;
    

    


    V_fsoma_tot(V_fsoma_tot>threasholds(i))=1;
    V_fsoma_tot(V_fsoma_tot<threasholds(i))=0;
    
    V_fsoma_masked = V_fsoma_tot;

    V_Rsoma_per_fsoma = V_SANDI_masked.*V_fsoma_masked;
    
    
    
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
set(figure(4), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(4), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/rsoma_vs_cbf_fsomamask.png');



figure(5), 
p = plot(threasholds, corr_values, '-o', 'LineWidth', 1, 'MarkerFaceColor','b', 'MarkerSize',3);
xlabel('PVE threashold');
ylabel('r');
title('corr\_coef Rsoma (masked with fsoma) vs CBF');
set(figure(5), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(5), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/corr_coef_rsomavscbf_fsomamask.png');

%Plot PVE maps threasholded
for i = index
    V_GM = V_GM_tot;

    V_GM(V_GM>0)=1;
    V_GM(V_GM<=0)=0;

    V_SANDI_masked=V_SANDI.*V_GM;

    V_fsoma_tot = V_fsoma;

    V_fsoma_tot(V_fsoma_tot>threasholds(i))=1;
    V_fsoma_tot(V_fsoma_tot<threasholds(i))=0;

    V_fsoma_masked = V_fsoma_tot;

    V_Rsoma_per_fsoma = V_SANDI_masked.*V_fsoma_masked;

    slices_SANDI{i} = V_Rsoma_per_fsoma(:,:,45);


end



   
fig=figure(6);
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
set(figure(6), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(6), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/rsoma_maps_threasholdingfsomamap.png');

%%%%
V_GM = V_GM_tot;
V_GM(V_GM>0)=1;
V_GM(V_GM<=0)=0;

V_fsoma_masked = V_fsoma.*V_GM;
V_CBF_masked = V_CBF.*V_GM;

r = corrcoef(V_fsoma_masked(:), V_CBF_masked(:)); 

figure(7), 
scatter(V_CBF_tot(:), V_fsoma_masked(:))
ylabel('fsoma(%)');
%ylabel(strcat(label,' masked with fsoma'));
xlabel('CBF (ml/100g/min)');
title("corr\_coef="+r(2));
saveas(figure(7), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/corr_coef_rsomavscbf.png');

V_SANDI = V_fsoma;

figure(8), 
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
set(figure(8), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(8), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/fsoma_vs_cbf_GMmask.png');

figure(9), 
p = plot(threasholds, corr_values, '-o', 'LineWidth', 1, 'MarkerFaceColor','b', 'MarkerSize',3);
xlabel('PVE threashold');
ylabel('r');
title('corr\_coef fsoma (masked with GM) vs CBF');
set(figure(9), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
saveas(figure(9), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/corr_coef_fsomavscbf_GMmask.png');


%%

%studying the relationship voxel by voxel: in voxel 1, cbf has value tot, fsoma value tot, De value tot.

%3d plot: CBF, fsoma, De

%3d plot: CBF, fsoma, Rsoma
parameter='De';%'Rsoma';%'fsoma'; insert parameter name compatible with SANDI files.
unit_of_measure = '(\mu m^{2}/ms)';%'(\mu m)';(%);
label_De = strcat(parameter,unit_of_measure);
img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_De = spm_read_vols(Vhdr);

V_De_masked=V_De.*V_GM;

parameter='fsoma';
unit_of_measure = '%';
label_fsoma = strcat(parameter,unit_of_measure);
img_path=strcat('/media/nas_rete/Vitality/clustering/sub-001_run-01_SANDI-fit_',parameter,'_2MNI.nii.gz');
Vhdr = spm_vol(img_path);
V_fsoma = spm_read_vols(Vhdr);

V_fsoma_masked=V_fsoma.*V_GM;

%V_CBF_masked

label_CBF = 'CBF (ml/100g/min)';

figure, plot3(V_De_masked(:),V_CBF_masked(:),V_fsoma_masked(:), '*');
xlabel(label_De);
ylabel(label_CBF);
zlabel(label_fsoma);

figure, hist(V_De(:))

figure, hist(V_fsoma_masked(:),100)

De_value = 0.5; %or any other value is allowed since it then collapses all data
V_De_masked = De_value + zeros(numel(V_De_masked),1);
figure, plot3(V_De_masked(:),V_CBF_masked(:),V_fsoma_masked(:), '*');
xlabel(label_De);
ylabel(label_CBF);
zlabel(label_fsoma);


 





v_De = V_De_masked(:);
v_CBF = V_CBF_masked(:);
v_fsoma = V_fsoma_masked(:);


step = 1000;
start = 1:step:numel(v_De);
last_step = rem(numel(v_De),step)-1;

for i = start
    if i < numel(v_De) - rem(numel(v_De),step)

        figure, plot3(v_De(i:i+step),v_CBF(i:i+step),v_fsoma(i:i+step),'*');
        xlabel(label_De);
        ylabel(label_CBF);
        zlabel(label_fsoma);

    else
        figure, plot3(v_De(i:i+last_step),v_CBF(i:i+last_step),v_fsoma(i:i+last_step),'*');
        xlabel(label_De);
        ylabel(label_CBF);
        zlabel(label_fsoma);
    end

end

%fsoma vs De

r = corrcoef(v_fsoma, v_De); 

figure, scatter(v_fsoma,v_De)
xlabel('fsoma(%)');
%ylabel(strcat(label,' masked with fsoma'));
ylabel('De(\mu m^{2}/ms)');
title("corr\_coef="+r(2));


%parameters correlated by the model


%try cmro2, CBF and SANDI parameter which come from different models

%%
%glm


%%
%Visualize map in sagital view
% img_path = '/media/nas_rete/Work_manuela/SANDI-fit_Rsoma2x2x2.nii';
% Vhdr = spm_vol(img_path);
% V = spm_read_vols(Vhdr);
% figure, imagesc(V(30,:,:)) 
% title('Rsoma map 2x2x2')
% imshow(rot90(squeeze(V(30,:,:))),[])



