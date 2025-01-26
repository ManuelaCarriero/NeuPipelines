%This function selects all GM areas multiplying by GM PVE = 0.5
% vectorizes all 3D images matrices 
% makes the median through all subjects
% remove unphysical CBF and rsoma values from SANDI parameters
% (you should modify the threshold over rsoma)


function [v_fsoma, v_rsoma, v_fneurite, v_fextra, v_De, v_Din]=vectorize(n_subjs, V_GM_tot, V_CBF_tots, V_fsoma_tots, V_fneurite_tots, V_fextra_tots, V_Din_tots, V_De_tots, V_rsoma_tots)

V_CBF_matrix=[];
V_rsoma_matrix=[];
V_fsoma_matrix=[];
V_fneurite_matrix=[];
V_Din_matrix=[];
V_De_matrix=[];
V_fextra_matrix=[];

V_GM = V_GM_tot;
threashold = 0.5;

for i=1:n_subjs
    V_fsoma_tot=V_fsoma_tots{i};
    V_fsoma_tot=V_fsoma_tot.*(V_GM>threashold);
    V_fsoma_matrix(i,:)=V_fsoma_tot(:)';
    V_CBF_tot=V_CBF_tots{i};
    V_CBF_tot=V_CBF_tot.*(V_GM>threashold);%.*(V_fsoma_tot>0.5);
    V_CBF_matrix(i,:)=V_CBF_tot(:)';
    V_rsoma_tot=V_rsoma_tots{i};
    V_rsoma_tot=V_rsoma_tot.*(V_GM>threashold);%.*(V_fsoma_tot>0.5);
    V_rsoma_matrix(i,:)=V_rsoma_tot(:)';
    V_fneurite_tot=V_fneurite_tots{i};
    V_fneurite_tot=V_fneurite_tot.*(V_GM>threashold);
    V_fneurite_matrix(i,:)=V_fneurite_tot(:)';
    V_Din_tot=V_Din_tots{i};
    V_Din_tot=V_Din_tot.*(V_GM>threashold);
    V_Din_matrix(i,:)=V_Din_tot(:)';
    V_De_tot=V_De_tots{i};
    V_De_tot=V_De_tot.*(V_GM>threashold);
    V_De_matrix(i,:)=V_De_tot(:)';
    V_fextra_tot=V_fextra_tots{i};
    V_fextra_tot=V_fextra_tot.*(V_GM>threashold);
    V_fextra_matrix(i,:)=V_fextra_tot(:)';
end

median_CBF = median(V_CBF_matrix,1);
median_rsoma = median(V_rsoma_matrix,1);
median_fsoma = median(V_fsoma_matrix,1);
median_fneurite = median(V_fneurite_matrix,1);
median_Din = median(V_Din_matrix,1);
median_De = median(V_De_matrix,1);
median_fextra = median(V_fextra_matrix,1);%median


%remove zero values of median CBF and higher values

indices_cbf = [];
for ii = 1:length(median_CBF)%V_CBF_tot
    if median_CBF(ii)>100 || median_CBF(ii)<5 %|| ii == any(indices_rsoma_zeros)
        %v_CBF_reduced(ii)=[];
        indices_cbf(end+1)=ii;
    end
end
% 
indices_rsoma = [];
for ii = 1:length(median_rsoma)%V_CBF_tot
    if median_rsoma(ii)<5 
        %v_CBF_reduced(ii)=[];
        indices_rsoma(end+1)=ii;
    end
end
%try delete elements in 3d image and plot it again

%Remove corresponding values of CBF to Rsoma=0
%indices_tot=indices_cbf;
indices_tot=cat(2,indices_rsoma,indices_cbf);

indices_tot_unique = unique(indices_tot);

v_fsoma = median_fsoma;
v_rsoma = median_rsoma;
v_fneurite = median_fneurite;
v_fextra = median_fextra;
v_Din = median_Din;
v_De = median_De;
v_CBF = median_CBF;


v_fsoma(indices_tot_unique)=[];
v_rsoma(indices_tot_unique)=[];
v_fneurite(indices_tot_unique)=[];
v_fextra(indices_tot_unique)=[];
v_Din(indices_tot_unique)=[];
v_De(indices_tot_unique)=[];
v_CBF(indices_tot_unique)=[];

end