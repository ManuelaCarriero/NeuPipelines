mean_energy=mean_energy_high_n_voxels;
mean_fsoma=mean_fsoma_high_n_voxels;
mean_rsoma=mean_rsoma_high_n_voxels;
mean_fsup=mean_fsup_high_n_voxels;

% mean_De=mean_De_high_n_voxels;
% mean_fneurite=mean_fneurite_high_n_voxels;
% mean_fextra=mean_fextra_high_n_voxels;
% mean_fc=mean_fc_high_n_voxels;

energy_scored=zscore(mean_energy);
fsoma_scored=zscore(mean_fsoma);
rsoma_scored=zscore(mean_rsoma);
fsup_scored=zscore(mean_fsup);

% De_scored=zscore(mean_De);
% fneurite_scored=zscore(mean_fneurite);
% fextra_scored=zscore(mean_fextra);
% fc_scored=zscore(mean_fc);

energy_tr=energy_scored';
fsoma_tr=fsoma_scored';
rsoma_tr=rsoma_scored';
fsup_tr=fsup_scored';

% De_tr=De_scored';
% fneurite_tr=fneurite_scored';
% fextra_tr=fextra_scored';
% fc_tr=fc_scored';
%% GLM and PLS models (single run)
%GLM

%one with everything
% y=energy_tr;
% X=[fsoma_tr rsoma_tr fsup_tr De_tr fneurite_tr fextra_tr fc_tr];

%one with only soma characteristics
y=energy_tr;
X=[fsoma_tr rsoma_tr fsup_tr ];

fitglm(X,y,'linear')

%PLS with cross validation
s_x=size(X);
ncomp=s_x(2);

[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X,y,ncomp); %with n_comp=1 fsoma is positive as we expect from corr.

%nested cross validation
%k-fold cross validation to select the right number of components
[XL,yl,XS,YS,beta,PCTVAR,PLSmsep] = plsregress(X,y,ncomp,'CV',10);

min_err=min(PLSmsep(2,:));%min MSE for the predictor variable
idx=find(PLSmsep(2,:)==min_err);
new_ncomp=idx-1;

[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X,y,new_ncomp);


figure, plot(1:ncomp+1,PLSmsep(2,:))
%% Bootstrapping to obtain a distribution of beta values and the error
%This is to create the random data stream for reproducibility
%s=RandStream('mlfg6331_64','Seed',i);
for i=1:1000
    disp(i)    
    %try to sample from 1 to the length and then use the index 
    % for each parameter

    idx_samples=datasample(1:length(energy_tr),length(energy_tr),'Replace',true);

    energy_sampled = energy_tr(idx_samples);
    rsoma_sampled = rsoma_tr(idx_samples);
    fsoma_sampled = fsoma_tr(idx_samples);
    fsup_sampled = fsup_tr(idx_samples);
% 
%     De_sampled = De_tr(idx_samples);
%     fneurite_sampled = fneurite_tr(idx_samples);
%     fextra_sampled = fextra_tr(idx_samples);
%     fc_sampled = fc_tr(idx_samples);
    
    
    y = energy_sampled;
    X = [rsoma_sampled fsoma_sampled fsup_sampled];%De_tr fneurite_tr fextra_tr fc_tr

    
    size_x=size(X);
    ncomp=size_x(2);%ncomp<=nvariables

    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X,y,ncomp); %with n_comp=1 fsoma is positive as we expect from corr.

    %nested cross validation
    %k-fold cross validation to select the right number of components
    [XL,yl,XS,YS,beta,PCTVAR,PLSmsep] = plsregress(X,y,ncomp,'CV',10);

    min_err=min(PLSmsep(2,:));%min MSE for the predictor variable
    idx=find(PLSmsep(2,:)==min_err);
    new_ncomp=idx-1;

    if new_ncomp<=0
        beta_rsoma(i)=NaN;
        beta_fsoma(i)=NaN;
        beta_fsup(i)=NaN;
    else
    %new_ncomp=2;
    [XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X,y,new_ncomp);

    beta_rsoma(i)=beta(2);
    beta_fsoma(i)=beta(3);
    beta_fsup(i)=beta(4);
    end

end

beta_fsoma(isnan(beta_fsoma))=[];
beta_rsoma(isnan(beta_rsoma))=[];
beta_fsup(isnan(beta_fsup))=[];

figure, hist(beta_fsoma);
xlabel('fsoma beta coefficients','FontWeight','bold','FontSize',12);
ylabel('Counts','FontWeight','bold','FontSize',12);
title(strcat(energy,'vs Soma microparameters'));
grid on
figure, hist(beta_rsoma);
xlabel('Rsoma beta coefficients','FontWeight','bold','FontSize',12);
ylabel('Counts','FontWeight','bold','FontSize',12);
title(strcat(energy,'vs Soma microparameters'));
grid on
figure, hist(beta_fsup);
xlabel('fsup beta coefficients','FontWeight','bold','FontSize',12);
ylabel('Counts','FontWeight','bold','FontSize',12);
title(strcat(energy,'vs Soma microparameters'));
grid on




%For the central limit theorem, the mean has to converge to zero
mean_beta_fsoma=mean(beta_fsoma);
mean_beta_rsoma=mean(beta_rsoma);
mean_beta_fsup=mean(beta_fsup);

%Calculate the expectation value
micro_parameter='Rsoma';
if strcmp(micro_parameter,'fsoma')
    beta=beta_fsoma;
elseif strcmp(micro_parameter,'Rsoma')
    beta=beta_rsoma;
elseif strcmp(micro_parameter,'fsup')
    beta=beta_fsup;
end

[N,edge]=histcounts(beta);
length(edge)
means=[];
for i=1:length(edge)
    if i < length(edge)
     mean=(edge(i)+edge(i+1))/2;
     means(i)=mean;
    else
        break
    end
end

if strcmp(micro_parameter,'fsoma')   
    weighted_fsoma=means.*N;
    E_fsoma=sum(weighted_fsoma)./sum(N);%-0.0054
elseif strcmp(micro_parameter,'Rsoma')
    weighted_rsoma=means.*N;
    E_Rsoma=sum(weighted_rsoma)./sum(N);%-0.0207
elseif strcmp(micro_parameter,'fsup')
    weighted_fsup=means.*N;
    E_fsup=sum(weighted_fsup)./sum(N);%0.0038
end

std_beta_fsoma=std(beta_fsoma);
std_beta_rsoma=std(beta_rsoma);
std_beta_fsup=std(beta_fsup);

%%
figure,
plot(0:ncomp,PLSmsep(2,:),'b-o');
xlabel('Number of components','FontWeight','bold');
ylabel('Estimated Mean Squared Prediction Error','FontWeight','bold');
legend({'PLSR'},'location','NE');
grid on
%

figure,
plot(1:ncomp,cumsum(100*PCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in y');

figure, bar(1:ncomp,PCTVAR(2,:))
xlabel('Components');
ylabel('Variance Explained (%)');
grid on

%% Corr coef matrix

A=[fsoma_scored' rsoma_scored' fsup_scored'];
corr_matr=corrcoef(A);

% Plot the data using imagesc() for later comparison
clrLim = [-1,1];

N = 255;
cmap = uint8([N*ones(1,N), N:-1:0; 0:N-1, N:-1:0; 0:N, N*ones(1,N)].');

figure()
imagesc(corr_matr)
xticklabels({'','fsoma','','Rsoma','','fsup'})
xtickangle(0);
yticklabels({'','fsoma','','Rsoma','','fsup'})
ax=gca;
colormap(gca,cmap);
colorbar();
ax.FontWeight='bold';
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
caxis(clrLim);
axis equal
axis tight




%% Other kind of analysis 
% (if you want to study the property of latent space 
% or check the goodness of fit) 

yfit = [ones(size(X,1),1) X]*beta;
figure,
plot(y,yfit,'o')

TSS = sum((y-mean(y)).^2);
RSS = sum((y-yfit).^2);
Rsquared = 1 - RSS/TSS;

figure,
plot(1:6,stats.W,'o-')
legend({'c1','c2','c3','c4','c5','c6'},'Location','best')
xlabel('Predictor')
ylabel('Weight')


figure,
scatter(XS(:,1),YS(:,1)) 

figure, 
%s = scatter(mean_energy,mean_SANDI);
s = plot(XS(:,1), YS(:,1),'o');
xlabel('XS PLS1','FontSize',15,'FontWeight','bold');
ylabel('YS PLS1','FontSize',15,'FontWeight','bold');
s.LineWidth = 0.6;
s.MarkerEdgeColor = 'b';
s.MarkerFaceColor = [0 0.5 0.5];
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(gca,'box','off')
grid on

figure,
scatter(XS(:,1),XS(:,4))

figure,
bar(1:numel(X(:,1)),XS(:,6))

%% Print correlation coefficients

disp(beta)