mean_energy=mean_energy_high_n_voxels;
mean_fsoma=mean_fsoma_high_n_voxels;
mean_rsoma=mean_rsoma_high_n_voxels;
mean_fsup=mean_fsup_high_n_voxels;

energy_scored=zscore(mean_energy);
fsoma_scored=zscore(mean_fsoma);
rsoma_scored=zscore(mean_rsoma);
fsup_scored=zscore(mean_fsup);

energy_tr=energy_scored';
fsoma_tr=fsoma_scored';
rsoma_tr=rsoma_scored';
fsup_tr=fsup_scored';



y = energy_tr;
X = [rsoma_tr fsoma_tr fsup_tr]; %
%X=[x_rsoma  x_fsoma x_Din x_fextra x_fneurite x_De]; 

%nested cross validation 
size_x=size(X);
ncomp=size_x(2);%ncomp<=nvariables

[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X,y,ncomp); %with n_comp=1 fsoma is positive as we expect from corr.

%k-fold cross validation to select the right number of components
[XL,yl,XS,YS,beta,PCTVAR,PLSmsep] = plsregress(X,y,ncomp,'CV',10);





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




%%
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

disp(beta);