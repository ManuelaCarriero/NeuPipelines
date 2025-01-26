%This function makes voxel wise analysis between energy metrics and SANDI
%parameter

%voxel wise analysis
%one subj
%threshold with GM

%OUTPUT:
%Scatterplots of energy metrics vs SANDI parameter for different values of PVE;
%correlation value as function of PVE value;
%SANDI map for different values of PVE.



function []=plot_subj_corr(threasholds, index,  V_SANDI_tot, V_CBF_tot, V_GM_tot,SANDI_parameter,SANDI_unit_of_measure,cbf_unit_of_measure) 



figure,
for i = index

    V_GM = V_GM_tot;

    V_GM(V_GM>=threasholds(i))=1;
    V_GM(V_GM<threasholds(i))=0;%you could write V_GM(V_GM<1)=0.

    V_SANDI = V_SANDI_tot;
    %Mask SANDI map with grey matter
    V_SANDI_masked=V_SANDI.*V_GM;



    V_CBF = V_CBF_tot;
    V_CBF_masked = V_CBF.*V_GM;

    slices_SANDI{i} = V_SANDI_masked(:,:,45);
    slices_CBF{i} = V_CBF_masked(:,:,45);

    %%%%
    v_SANDI=V_SANDI_masked(:);
    v_CBF=V_CBF_masked(:);

    %     indices_rsoma_zeros=[];
    %     for j = 1:numel(v_SANDI)
    %         if v_SANDI(j)<min_rsoma
    %             indices_rsoma_zeros(end+1)=j;
    %         end
    %     end


    %prepare vectors to substitute with reduced data points
    v_CBF_reduced = v_CBF;
    v_SANDI_reduced = v_SANDI;

    %Remove CBF=0, 95 and 5 prctile values from CBF distribution

    %     prctile99=prctile(v_CBF,99);
    %     prctile5=prctile(v_CBF,5);

    indices_cbf = [];
    for ii = 1:length(v_CBF)
        if v_CBF(ii)>100 || v_CBF(ii)<5 %|| ii == any(indices_rsoma_zeros)
            %v_CBF_reduced(ii)=[];
            indices_cbf(end+1)=ii;
        end
    end
    %try delete elements in 3d image and plot it again

    %Remove corresponding values of CBF to Rsoma=0
    %indices_tot=cat(2,indices_rsoma_zeros, indices_cbf);
    indices_tot=indices_cbf;

    indices_tot_unique = unique(indices_tot);

    v_CBF_reduced(indices_tot_unique)=[];
    v_SANDI_reduced(indices_tot_unique)=[];

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
    r = corrcoef(v_SANDI_reduced, v_CBF_reduced);
    corr_values(i) = r(2);

    P = polyfit(v_SANDI_reduced,v_CBF_reduced,1);
    yfit = P(1)*v_SANDI_reduced+P(2);

    %v_CBF_reduced(1:numel(v_CBF_reduced)*0.1)
    %v_CBF_reduced(numel(v_CBF_reduced)*0.1:numel(v_CBF_reduced)*0.1*2)
    subplot(4,4,i)
    scatter(v_SANDI_reduced,v_CBF_reduced)
    %numel(v_CBF_reduced)
    hold on;
    plot(v_SANDI_reduced,yfit,'r--');
    thr = num2str(threasholds(i));
    title(strcat("thr>=",thr))
    ylabel(strcat('CBF', cbf_unit_of_measure));
    %ylabel(strcat(label,' masked with fsoma'));
    label=strcat(SANDI_parameter,SANDI_unit_of_measure);
    xlabel(strcat(label,' Signal'));
end
sgtitle("CBF vs fneurite increasing GM PVE threashold (thr)");
%set(figure, 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);

% path_to_image=strcat(path,'/rsoma_vs_cbf_GMmask_min_rsoma',min_rsoma_str,'_cbf',min_cbf_str,'.png');
% path_to_image=strcat(path,'/fsoma_vs_cbf_GMmask_min_rsoma.png');
% saveas(figure(1), path_to_image);



%correlation curve
figure, 
p = plot(threasholds, corr_values, '-o', 'LineWidth', 1, 'MarkerFaceColor','b', 'MarkerSize',3);
xlabel('PVE threashold', 'FontSize',15);
ylabel('r','FontSize',15);
title('corr\_coef fneurite (masked with GM) vs CBF');
%set(figure, 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
% path_to_image=strcat(path,'/corr_coef_fneuritevscbf_GMmask_min_rsoma',min_rsoma_str,'_cbf',min_cbf_str,'.png');
% saveas(figure(2), path_to_image);

%for the last PVE threashold
% figure(3), hist(v_CBF_reduced);
% title('CBF distribution')
% xlabel('CBF (ml/100g/min)');
% ylabel('counts');
% path_to_image=strcat(path,'/CBFdistribution_min_rsoma',min_rsoma_str,'_cbf',min_cbf_str,'.png');
% saveas(figure(3), path_to_image);

% figure(4), hist(v_SANDI_reduced);
% title('Rsoma distribution')
% xlabel(label);
% ylabel('counts');
% path_to_image=strcat(path,'/Rsomadistribution_min_rsoma',min_rsoma_str,'_cbf',min_cbf_str,'.png');
% saveas(figure(4), path_to_image);



fig=figure;
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
sgtitle(strcat(SANDI_parameter,' map masked with GM varying thr'));
% maxi=max(V_CBF_masked, [],'all');%V_CBF_masked is of the last interaction so it is the maximum value.
% dy=maxi/10;
%labels = round(0:dy:maxi);

c = colorbar(h,'Position',[0.93 0.4 0.019 0.4]);  % 'TickLabels',[labels]) attach colorbar to h
%c = colorbar(h,'Location','southoutside');

% c.Limits=[0 maxi];
% set(figure(5), 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
% saveas(figure(5), '/media/nas_rete/Work_manuela/DWI_En_modeling/main/rsoma_maps_threasholdingPVE.png');
end