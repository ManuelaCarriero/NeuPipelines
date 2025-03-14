%This function makes voxel wise analysis between energy metrics and SANDI
%parameter

%makes the median of many subjects
%considers GM threshold and fsoma threshold > 0.4

%This function allows to choose if calculate with different values of GM
%PVE or just one value
%index argument contains either the value of the single threshold
%or the indices values corresponding to each threshold value



%OUTPUT:
%Scatterplots of energy metrics vs SANDI parameter for different values of PVE;
%correlation value as function of PVE value;
%SANDI map for different values of PVE.

%if one PVE value is chosen:
%scatterplot of energy metric vs SANDI parameter and corr value.

function []=plot_subjs_corr(threasholds, index, V_GM_tot, V_fsoma_tots, V_SANDI_tots, V_energy_tots, n_subjs, low_thr,  subsample, V_fsoma_thr, micro_parameter,micro_parameter_unit_of_measure, energy, energy_unit_of_measure,subsample_perc, fig_scatter,fig_corr, fig_slices)
thr=size(index);

if thr(2)>1
    if strcmp(fig_scatter,'on')
    figure,
    for i = index
    
        V_energy_vecs=[];
        %V_fsoma_vecs=[];
        V_SANDI_vecs=[];
        
        V_energy_matrices = {};
        V_SANDI_matrices = {};
    
        V_GM = V_GM_tot;
    
        %threashold_fsoma = 0.5;
    
        for j=1:n_subjs
            V_GM = V_GM_tot;
            V_GM(V_GM>threasholds(i))=1;
            V_GM(V_GM<threasholds(i))=0;
    
            V_fsoma_tot=V_fsoma_tots{j};
    
            V_energy_tot=V_energy_tots{j};
            V_energy_tot=V_energy_tot.*V_GM.*(V_fsoma_tot>V_fsoma_thr);
            V_energy_vecs(j,:)=V_energy_tot(:)';
        
            V_energy_matrices{j} = V_energy_tot;


    
    
            V_SANDI_tot=V_SANDI_tots{j};%CHANGE 
    
    
            V_SANDI_tot=V_SANDI_tot.*V_GM.*(V_fsoma_tot>V_fsoma_thr);
            V_SANDI_vecs(j,:)=V_SANDI_tot(:)';
            V_SANDI_matrices{j} = V_SANDI_tot;

        end


    
        fourdmatrix_CBF=cat(4,V_energy_matrices{:});
        fourdmatrix_SANDI=cat(4,V_SANDI_matrices{:});
    
        median_matrix_energy=median(fourdmatrix_CBF,4);
        median_matrix_SANDI=median(fourdmatrix_SANDI,4);

        if strcmp(fig_slices,'on')
            slice_GM=V_GM.*(V_GM>threasholds(i));
            slices_GM{i}=slice_GM(:,:,45);
            slices_SANDI{i}=median_matrix_SANDI(:,:,45);
            slices_energy{i}=median_matrix_energy(:,:,45);
        end

        median_energy = median(V_energy_vecs,1);
       
        %median_fsoma = median(V_fsoma_matrix,1);
        median_SANDI = median(V_SANDI_vecs,1);
        
    
    
        %remove zero values of CBF and higher values
    
        indices_energy = [];
        for ii = 1:length(median_energy)
            if  median_energy(ii)<low_thr %median_energy(ii)>up_thr ||
                
                indices_energy(end+1)=ii;
            end
        end
    
    
        %Remove corresponding values of CBF to Rsoma=0
        %indices_tot=cat(2,indices_rsoma_zeros, indices_cbf);
        indices_tot=indices_energy;
    
        indices_tot_unique = unique(indices_tot);
    
    
        v_SANDI = median_SANDI;
        v_energy = median_energy;
    
    
        v_SANDI(indices_tot_unique)=[];
        v_energy(indices_tot_unique)=[];
    
    
        r = corrcoef(v_energy, v_SANDI);
        corr_values(i) = r(2);
%
%         P = polyfit(v_SANDI,v_energy,1);
%         yfit = P(1)*v_SANDI+P(2);


        
        %
        index_permutation=randperm(length(v_SANDI));
        %you must order CBF and SANDI in the same way
        v_SANDI_permutated=v_SANDI(index_permutation);
        v_energy_permutated = v_energy(index_permutation);

        %fit reduced
        P = polyfit(v_SANDI_permutated(1:1:round(numel(v_SANDI)*subsample_perc)),v_energy_permutated(1:1:round(numel(v_SANDI)*subsample_perc)),1);
        yfit = P(1)*v_SANDI_permutated(1:1:round(numel(v_SANDI)*subsample_perc))+P(2);

        %

        subplot(3,3,i)
        scatter(v_SANDI_permutated(1:1:round(numel(v_SANDI)*subsample_perc)),v_energy_permutated(1:1:round(numel(v_SANDI)*subsample_perc)))%v_SANDI,v_CBF
        %v_SANDI_permutated(1:1:round(numel(v_SANDI)*0.001)),v_energy_permutated(1:1:round(numel(v_SANDI)*0.001))
        %numel(v_energy_reduced)
        hold on;
        %p=plot(v_SANDI,yfit,'--','LineWidth',2.5,'Color',"#000000");%'r--',
        p=plot(v_SANDI_permutated(1:1:round(numel(v_SANDI)*subsample_perc)),yfit,'--','LineWidth',2.5,'Color',"#000000");%'r--',
        p.LineStyle='--';
        thr = num2str(threasholds(i));
        title(strcat("thr>=",thr))
        ylabel(strcat(energy,energy_unit_of_measure),'FontWeight','bold');
        %ylabel(strcat(label,' masked with fsoma'));
        xlabel(strcat(micro_parameter,micro_parameter_unit_of_measure),'FontWeight','bold');
    
    end
    %n_subjs=num2str(n_subjs);%if you want to instert this information in the title
    sgtitle(strcat(energy,' vs ',micro_parameter,' increasing GM PVE threashold (thr)'),'fontweight','bold');
    %set(figure, 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
    x0=400;
    y0=400;
    width=750;
    height=651;
    set(gcf,'position',[x0,y0,width,height])
    end
    
    if strcmp(fig_corr,'on')
    %Plot correlation curve
    figure,
    p = plot(threasholds, corr_values, '-o', 'LineWidth', 1, 'MarkerFaceColor','b', 'MarkerSize',7);
    xlabel('PVE threshold', 'FontSize',15,'FontWeight','bold');
    ylabel('r','FontSize',15,'FontWeight','bold');
    set(get(gca, 'XAxis'), 'FontWeight', 'bold');
    set(get(gca, 'YAxis'), 'FontWeight', 'bold');
    set(gca,'box','off')
    grid on
    %title(strcat('corr\_coef', energy, 'vs', micro_parameter));
    %set(figure, 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
    end

    %Plot slices
    if strcmp(fig_slices,'on')

        fig=figure;
        for i = index
            subplot(3,3,i)
            imagesc(rot90(slices_GM{i}))
            thr = num2str(threasholds(i));
            title(strcat("thr>=",thr))
            axis equal
            axis off
        end
        caxis([0,1000]);
        h = axes(fig,'visible','off');
        %     sgtitle(strcat(SANDI_parameter,' map masked with GM varying thr'));
        colorbar(h,'Position',[0.93 0.4 0.019 0.4]);

        x0=400;
        y0=400;
        width=750;
        height=651;
        set(gcf,'position',[x0,y0,width,height])

        fig=figure;
        for i = index
            subplot(3,3,i)
            imagesc(rot90(slices_SANDI{i}))
            thr = num2str(threasholds(i));
            title(strcat("thr>=",thr))
            caxis([5,16]);
            axis equal
            axis off
        end
        caxis([5,16]);
        h = axes(fig,'visible','off');
        %     sgtitle(strcat(SANDI_parameter,' map masked with GM varying thr'));
        colormap('jet')
        colorbar(h,'Position',[0.93 0.4 0.019 0.4]);
        set(h, 'CLim',[5 15])

        x0=400;
        y0=400;
        width=750;
        height=651;
        set(gcf,'position',[x0,y0,width,height])

        fig=figure;
        for i = index
            subplot(3,3,i)
            imagesc(rot90(slices_energy{i}))
            thr = num2str(threasholds(i));
            title(strcat("thr>=",thr))
            caxis([5,130]);
            axis equal
            axis off
        end
        %caxis([0,1000]);
        h = axes(fig,'visible','off');
        %     sgtitle(strcat(SANDI_parameter,' map masked with GM varying thr'));
        colorbar(h,'Position',[0.93 0.4 0.019 0.4]);
        set(h, 'CLim',[5 130])
        
        x0=400;
        y0=400;
        width=750;
        height=651;
        set(gcf,'position',[x0,y0,width,height])


    end

    

else
%in this case: if subsample is "all", it shows only one plot for the GM threashold chosen 
%else it will plot different plots corresponding to different subsampling
%of the total dataset (in order to see better the structure of data).

        V_CBF_vecs=[];
        %V_fsoma_vecs=[];
        V_SANDI_vecs=[];
        
        V_CBF_matrices = {};
        V_SANDI_matrices = {};
    
        V_GM = V_GM_tot;
    
        %threashold_fsoma = 0.5;
    
        for j=1:n_subjs
            V_GM = V_GM_tot;
            V_GM(V_GM>index)=1;
            V_GM(V_GM<index)=0;
    
            V_fsoma_tot=V_fsoma_tots{j};
    
            V_CBF_tot=V_CBF_tots{j};
            V_CBF_tot=V_CBF_tot.*V_GM.*(V_fsoma_tot>V_fsoma_thr);
            V_CBF_vecs(j,:)=V_CBF_tot(:)';
        
            V_CBF_matrices{j} = V_CBF_tot;
    
    
            V_SANDI_tot=V_SANDI_tots{j};%CHANGE 
    
    
            V_SANDI_tot=V_SANDI_tot.*V_GM.*(V_fsoma_tot>V_fsoma_thr);
            V_SANDI_vecs(j,:)=V_SANDI_tot(:)';
            V_SANDI_matrices{j} = V_SANDI_tot;
        end
  
        fourdmatrix_CBF=cat(4,V_CBF_matrices{:});
        fourdmatrix_SANDI=cat(4,V_SANDI_matrices{:});
    
        median_matrix_CBF=median(fourdmatrix_CBF,4);
        median_matrix_SANDI=median(fourdmatrix_SANDI,4);
       
        median_CBF = median(V_CBF_vecs,1);
       
        %median_fsoma = median(V_fsoma_matrix,1);
        median_SANDI = median(V_SANDI_vecs,1);
        
    
    
        %remove zero values of CBF and higher values
    
        indices_cbf = [];
        for ii = 1:length(median_CBF)
            if median_CBF(ii)<low_thr  % median_CBF(ii)>up_thr || 
                
                indices_cbf(end+1)=ii;
            end
        end
    
    
        %Remove corresponding values of CBF to Rsoma=0
        %indices_tot=cat(2,indices_rsoma_zeros, indices_cbf);
        indices_tot=indices_cbf;
    
        indices_tot_unique = unique(indices_tot);
    
    
        v_SANDI = median_SANDI;
        v_CBF = median_CBF;
    
    
        v_SANDI(indices_tot_unique)=[];
        v_CBF(indices_tot_unique)=[];
    
    
        [r,p] = corrcoef(v_CBF, v_SANDI);
        corr_value = r(2);
        p_value = p(2);
        corr_coef_str=num2str(corr_value);
        p_value_str=num2str(p_value);
    
        P = polyfit(v_SANDI,v_CBF,1);
        yfit = P(1)*v_SANDI+P(2);
        
        if strcmp(subsample,'all')

        %v_CBF_reduced(1:numel(v_CBF_reduced)*0.1)
        %v_CBF_reduced(numel(v_CBF_reduced)*0.1:numel(v_CBF_reduced)*0.1*2)
            figure,
            scatter(v_SANDI,v_CBF)
            %numel(v_CBF_reduced)
            hold on;
            plot(v_SANDI,yfit,'r--');
            index=num2str(index);
            title(strcat("thr>=",index))
            ylabel('CBF (ml/100g/min)');
            %ylabel(strcat(label,' masked with fsoma'));
            xlabel('fsoma');
            txt = {strcat('r = ',corr_coef_str,'**')};%,strcat('p-value = ',p_value_str)
            text(max(v_SANDI)/2,max(v_CBF)/2, txt, 'FontWeight', 'bold','FontSize',12)

        else
            
            figure,
            scatter(v_SANDI,v_CBF)
            %numel(v_CBF_reduced)
            hold on;
            plot(v_SANDI,yfit,'r--');
            index=num2str(index);
            title(strcat("thr>=",index))
            ylabel('CBF (ml/100g/min)','FontWeight','bold');
            %ylabel(strcat(label,' masked with fsoma'));
            xlabel(micro_parameter,'FontWeight','bold');
            txt = {strcat('r = ',corr_coef_str,'**')};%,strcat('p-value = ',p_value_str)
            text(max(v_SANDI)/2,max(v_CBF)/2, txt, 'FontWeight', 'bold','FontSize',12)

            for i=1:10
                index_permutation=randperm(length(v_SANDI));
                %you must order CBF and SANDI in the same way
                v_SANDI_permutated=v_SANDI(index_permutation);
                v_CBF_permutated = v_CBF(index_permutation);
    
                figure,
                scatter(v_SANDI_permutated(1:1:round(numel(v_SANDI)*subsample_perc)),v_CBF_permutated(1:1:round(numel(v_SANDI)*subsample_perc)))
                %numel(v_CBF_reduced)
                hold on;
                plot(v_SANDI,yfit,'r--');
                %thr = num2str(threshold);
                %title(strcat("thr>=",thr));
                % thr = num2str(corr_value);
                % text(strcat("r=",corr_value));
                ylabel('CBF (ml/100g/min)','fontweight','bold');
                %ylabel(strcat(label,' masked with fsoma'));
                xlabel(micro_parameter,'fontweight','bold');
            end

        end




end 