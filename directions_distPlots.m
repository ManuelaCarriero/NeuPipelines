%This script is for the analysis 
% of diffusion gradient directions

% we consider four files of diffusion gradient directions:
% three files used in Genova research centre and one used in Chieti
% we analyze both input files given to the Scanner and output file
% the results show that input and output bvals files are equal
% input and output bvecs are different

%Dataset

%Genova (GE)

%a=bvec file after scanner session
%b=the first file GE
%c=the second file GE
%supfile=the data file in the supplementary material
%d=output of the second session

%b = importdata(path_to_file/1_SANDI_scheme_b6000minnesota.dvs);
%c = importdata(path_to_file/2_sandi2024.dvs);
%supfile = importdata(path_to_file/supfile.txt);
%a = importdata(path_to_file/cmrr_mbep2d_diff_7shell_20240126153140_4.bvec);
%d = importdata(path_to_file/cmrr_mbep2d_diff_7shell_20240229163117_2.bvec);

%input_bvecs = importdata(path_to_file/cmrrdiff_206dir_Chieti.dvs);
%output_bvecs = importdata(path_to_file/sub-01_cmrr_mbep2d_diff_7shell_tr3000.bvec);
%output_bval = importdata(path_to_file/sub-01_cmrr_mbep2d_diff_7shell_tr3000.bval);


%sphere
[xx, yy, zz] = sphere(20); %n=20 is the number of faces
figure, hold on
mesh(xx, yy, zz) %creates a mesh plot, which is a three-dimensional surface 
% that has solid edge colors and no face colors. 
% The function plots the values in matrix Z as heights above 
% a grid in the x-y plane defined by X and Y. 
% The edge colors vary according to the heights specified by Z.
plot3(a(:,1), a(:,2), a(:,3), '-*')

figure, plot3(a(:,1), a(:,2), a(:,3), 'o')
title('bvec file after first scanner session')




[xx, yy, zz] = sphere(20); %n=20 is the number of faces
figure, hold on
mesh(xx, yy, zz) 
plot3(c(:,1), c(:,2), c(:,3))
title('bvec input file second session')


%normalize
GE_rescaled_inputvecs=[];
%input_bvecs = input_bvecs';
for i = 1:1:size(c,1)
    j = c(i,:)./norm(c(i,:));
    GE_rescaled_inputvecs(end+1,:) = j;
end

GE_rescaled_inputvecs(isnan(GE_rescaled_inputvecs))=0;

[xx, yy, zz] = sphere(20); %n=20 is the number of faces
figure, hold on
mesh(xx, yy, zz) 
plot3(GE_rescaled_inputvecs(1,:), GE_rescaled_inputvecs(2,:), GE_rescaled_inputvecs(3,:), '*')








%Try
% [xx, yy, zz] = sphere(20);
% figure, hold on
% mesh(xx, yy, zz) 
% figure, plot3(a(:,1), a(:,2), a(:,3), '*')
% [U,V,W] = surfnorm(a(:,1), a(:,2), a(:,3)); %surfnorm(X,Y,Z) creates a three-dimensional 
% % surface plot and displays its surface normals. A surface normal is the imaginary 
% % line perpendicular to a flat surface, or perpendicular to the tangent plane at 
% % a point on a non-flat surface.surfnorm(X,Y,Z) creates a three-dimensional surface plot 
% % and displays its surface normals. A surface normal is 
% % the imaginary line perpendicular to a flat surface, 
% % or perpendicular to the tangent plane at a point on a non-flat surface.
% figure, quiver3(a(:,1), a(:,2), a(:,3), U, V, W,0)
% %surfnorm(X,Y,Z) creates a three-dimensional 
% % surface plot and displays its surface normals. 
% % A surface normal is the imaginary line perpendicular
% % to a flat surface, or perpendicular to the tangent 
% % plane at a point on a non-flat surface.



[xx, yy, zz] = sphere(20);
figure, hold on
mesh(xx, yy, zz) 
plot3(b(:,4), b(:,5), b(:,6), 'ro')
figure, plot3(b(:,4), b(:,5), b(:,6), 'ro')
title('first input file')

figure, plot3(c(:,4),c(:,5),c(:,6),'ko')
title('second input bvecs')

%d = d';
figure, plot3(d(:,1), d(:,2), d(:,3), 'co')
title('second output bvecs')



%normalized input bvecs
c1=c(:,4:6);
c1=c1./repmat(sqrt(sum(c1.^2'))',1,3);%repmat matrix with one raw and one column
% % it contains values of modules
% c2= (sum(c(:,4:6).^2'))'*6000;%sqrt(sum(c(:,4:6).^2'))'*6000

[xx, yy, zz] = sphere(20); 
figure, hold on
mesh(xx, yy, zz)
plot3(c1(:,1), c1(:,2), c1(:,3), 'bo')
title('Normalized input bvecs')

figure, plot3(c1(:,1), c1(:,2), c1(:,3), 'bo')
hold on, plot3(d(:,1), d(:,2), d(:,3), 'co')
title('input (blue) vs output bvecs data (cyan)')

[xx, yy, zz] = sphere(20); 
figure, hold on
mesh(xx, yy, zz)
plot3(supfile(:,2), supfile(:,3), supfile(:,4), 'mo')

figure, plot3(supfile(:,2), supfile(:,3), supfile(:,4), 'mo')
title('Supplementary material')



%figure, scatter(c2,bval(2:end))


%Check for phi and theta distributions

% Transform Cartesian coordinates to spherical
[phi_b, theta_b] = cart2sph(b(:,4), b(:,5), b(:,6));
[phi_a, theta_a] = cart2sph(a(:,1), a(:,2), a(:,3));
[phi_c, theta_c] = cart2sph(c(:,4), c(:,5), c(:,6));
[phi_supfile, theta_supfile] = cart2sph(supfile(:,2), supfile(:,3), supfile(:,4));

figure, hist(phi_a, 25)
title('bvec file after scanner session (phi)')

figure, hist(phi_b, 25)
title('first file (phi)')

figure, hist(theta_a, 25)
title('bvec file after scanner session (theta)')

figure, hist(theta_b, 25)
title('first file (theta)')

figure, hist(phi_c, 25)
title('second file (phi)')

figure, hist(theta_c, 25)
title('second file (theta)')

figure, hist(phi_supfile, 25)
title('supplementary file (phi)')

figure, hist(theta_supfile, 25)
title('supplementary file (theta)')



%Chieti

%Compare input vs output (official acquisition)

%bvals

%Check b_vals are equal for input and output
%you obtain input bvals making the quadratic sum times the maximum bval
input_bvecs = input_bvecs(:,4:6);

input_bvecs_tr = input_bvecs';
input_bvecs_square = input_bvecs_tr.^2;
input_bvals_calculated = sum(input_bvecs_square)*max(output_bvals);

figure, 

scatter(output_bvals(2:207), input_bvals_calculated);




%bvecs 3dplot

[xx, yy, zz] = sphere(20);
figure, hold on
mesh(xx, yy, zz)
plot3(output_bvecs(:,1), output_bvecs(:,2), output_bvecs(:,3), '*');%unit norm
title('Output bvecs');



%normalize inputbvecs

rescaled_inputvecs=[];
%input_bvecs = input_bvecs';
for i = 1:1:size(input_bvecs,1)
    j = input_bvecs(i,:)./norm(input_bvecs(i,:));
    rescaled_inputvecs(end+1,:) = j;
end

rescaled_inputvecs(isnan(rescaled_inputvecs))=0;



[xx, yy, zz] = sphere(20);
figure, hold on
mesh(xx, yy, zz)
plot3(rescaled_inputvecs(:,1), rescaled_inputvecs(:,2), rescaled_inputvecs(:,3), '*');
title('Input bvecs');






%figure, histogram(diff);%Further step: chack at which b-values are
%associated small and big difference using angles
%title('Difference distribution');
%xlabel('Diff');
%ylabel('Counts');
% Bval in uscita vs sum(bvalingresso)*6000 (I bval non cambiano, cambiano solo I bvec).
% Plot su sfera unitaria uscita norma 1; ingress single righe dividi per somma quadratica delle direzioni (x,y,z).

edges = unique(output_bvals);
counts = histc(output_bvals, edges);%11	10	20	30	40	40	56





%b-values = 0/500/1000/2000/3000/4000/6000 s/mm2 with 15/6/32/40/40/40/40 measurements per shell





%Testing all the output bvecs directions are equal

%1.we plot different output bvecs on the same sphere
%each one with a different color
%the last plot should cotain only dots of the same color.

for i = 1:1:7
    my_field = strcat('sub',num2str(i));
    s = "/storage/shared/PRINAntonello2022/BIDS/sub-0" + i + "/dwi/sub-0" + i + "_cmrr_mbep2d_diff_7shell_tr3000.bvec";
    data.(my_field) = struct([]);
    data.(my_field) = importdata(s);
 
end 

fields_cell = fieldnames(data);

color = {'ro' 'go' 'bo' 'co' 'mo' 'yo' 'ko'};

for i = 1:1:numel(fields_cell)
    output_bvecs = data.(fields_cell{i});
    c = color{i};
    plot3(output_bvecs(1,:), output_bvecs(2,:), output_bvecs(3,:),c);
    hold on
end

%Create cell of matrices (equivalent to a list of matrices in Python)

% mat_cell = {};
% for i = 1:1:numel(fields_cell)
%     output_bvecs = data.(fields_cell{i});
%     mat_cell{i} = output_bvecs;
% end

%2. Check if all matrices of the cell of matrices
%are equal to one matrix of reference
%This test is more accurate

for i = 1:1:numel(fields_cell)
    output_bvecs2 = data.(fields_cell{i});
    tf = isequal(data.(fields_cell{1}),output_bvecs2);
    disp(tf);
end


