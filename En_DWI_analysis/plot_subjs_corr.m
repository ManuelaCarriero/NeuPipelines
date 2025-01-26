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

function []=plot_subjs_corr(threasholds, index, V_GM_tot, V_fsoma_tots, V_rsoma_tots, V_CBF_tots, n_subjs, subsample)
thr=size(index);

if thr(2)>1
    figure,
    for i = index
    
        V_CBF_vecs=[];
        %V_fsoma_vecs=[];
        V_SANDI_vecs=[];
        
        V_CBF_matrices = {};
        V_SANDI_matrices = {};
    
        V_GM = V_GM_tot;
    
        %threashold_fsoma = 0.5;
    
        for j=1:n_subjs
            V_GM = V_GM_tot;
            V_GM(V_GM>threasholds(i))=1;
            V_GM(V_GM<threasholds(i))=0;
    
            V_fsoma_tot=V_fsoma_tots{j};
    
            V_CBF_tot=V_CBF_tots{j};
            V_CBF_tot=V_CBF_tot.*V_GM.*(V_fsoma_tot>0.4);
            V_CBF_vecs(j,:)=V_CBF_tot(:)';
        
            V_CBF_matrices{j} = V_CBF_tot;
    
    
            V_SANDI_tot=V_rsoma_tots{j};%CHANGE 
    
    
            V_SANDI_tot=V_SANDI_tot.*V_GM.*(V_fsoma_tot>0.4);
            V_SANDI_vecs(j,:)=V_SANDI_tot(:)';
            V_SANDI_matrices{j} = V_SANDI_tot;
        end
        slice_GM=V_GM.*(V_GM>threasholds(i));
        slices_GM{i}=slice_GM(:,:,45);
    
        fourdmatrix_CBF=cat(4,V_CBF_matrices{:});
        fourdmatrix_SANDI=cat(4,V_SANDI_matrices{:});
    
        median_matrix_CBF=median(fourdmatrix_CBF,4);
        median_matrix_SANDI=median(fourdmatrix_SANDI,4);
    
        slices_SANDI{i}=median_matrix_SANDI(:,:,45);
        slices_CBF{i}=median_matrix_CBF(:,:,45);
    
        median_CBF = median(V_CBF_vecs,1);
       
        %median_fsoma = median(V_fsoma_matrix,1);
        median_SANDI = median(V_SANDI_vecs,1);
        
    
    
        %remove zero values of CBF and higher values
    
        indices_cbf = [];
        for ii = 1:length(median_CBF)
            if median_CBF(ii)>100 || median_CBF(ii)<5 
                
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
    
    
        r = corrcoef(v_CBF, v_SANDI);
        corr_values(i) = r(2);
    
        P = polyfit(v_SANDI,v_CBF,1);
        yfit = P(1)*v_SANDI+P(2);
    
        %v_CBF_reduced(1:numel(v_CBF_reduced)*0.1)
        %v_CBF_reduced(numel(v_CBF_reduced)*0.1:numel(v_CBF_reduced)*0.1*2)
        subplot(4,4,i)
        scatter(v_SANDI,v_CBF)
        %numel(v_CBF_reduced)
        hold on;
        plot(v_SANDI,yfit,'r--');
        thr = num2str(threasholds(i));
        title(strcat("thr>=",thr))
        ylabel('CBF (ml/100g/min)');
        %ylabel(strcat(label,' masked with fsoma'));
        xlabel('fsoma');
    end
    sgtitle("CBF vs fsoma increasing GM PVE threashold (thr) (12 subjects)");
    %set(figure, 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);
    
    figure, 
    p = plot(threasholds, corr_values, '-o', 'LineWidth', 1, 'MarkerFaceColor','b', 'MarkerSize',3);
    xlabel('PVE threashold', 'FontSize',15);
    ylabel('r','FontSize',15);
    title('corr\_coef CBF (masked with GM) vs fsoma');
    %set(figure, 'Units', 'Normalized', 'Outerposition', [0 0 1 1]);

else
    
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
            V_CBF_tot=V_CBF_tot.*V_GM.*(V_fsoma_tot>0.4);
            V_CBF_vecs(j,:)=V_CBF_tot(:)';
        
            V_CBF_matrices{j} = V_CBF_tot;
    
    
            V_SANDI_tot=V_rsoma_tots{j};%CHANGE 
    
    
            V_SANDI_tot=V_SANDI_tot.*V_GM.*(V_fsoma_tot>0.4);
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
            if median_CBF(ii)>100 || median_CBF(ii)<5 
                
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
            ylabel('CBF (ml/100g/min)');
            %ylabel(strcat(label,' masked with fsoma'));
            xlabel('fsoma');
            txt = {strcat('r = ',corr_coef_str,'**')};%,strcat('p-value = ',p_value_str)
            text(max(v_SANDI)/2,max(v_CBF)/2, txt, 'FontWeight', 'bold','FontSize',12)

            for i=1:10
                index_permutation=randperm(length(v_SANDI));
                %you must order CBF and SANDI in the same way
                v_SANDI_permutated=v_SANDI(index_permutation);
                v_CBF_permutated = v_CBF(index_permutation);
    
                figure,
                scatter(v_SANDI_permutated(1:1:round(numel(v_SANDI)*0.01)),v_CBF_permutated(1:1:round(numel(v_SANDI)*0.01)))
                %numel(v_CBF_reduced)
                hold on;
                plot(v_SANDI,yfit,'r--');
                %thr = num2str(threshold);
                %title(strcat("thr>=",thr));
                % thr = num2str(corr_value);
                % text(strcat("r=",corr_value));
                ylabel('CBF (ml/100g/min)','fontweight','bold');
                %ylabel(strcat(label,' masked with fsoma'));
                xlabel('Rsoma (\mum)','fontweight','bold');
            end

        end




end 