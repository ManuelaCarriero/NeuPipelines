img_path='/media/nas_rete/GLOVE_STUDY/DDC/derivatives/pil002/resting/dwi/pil002_resting_dwi.nii.gz';
Vhdr=spm_vol(img_path);
V=spm_read_vols(Vhdr);

disp('tot volumes');

tot_fourth_length=length(V(1,1,1,:));
%unique(bval)
%% load matrix

load('/media/nas_rete/GLOVE_STUDY/DDC/derivatives/pil002/resting/dwi/bval');
load('/media/nas_rete/GLOVE_STUDY/DDC/derivatives/pil002/resting/dwi/bvec');
bval_reduced=bval;
bval_tot=bval;

bvec_reduced=bvec;
bvec_tot=bvec;

disp('tot matrix')
%% load indices
idx_6000=find(bval_tot>5990);
idx_200=find(bval_tot==200);

disp('laod indices')

%%
bval_reduced(idx_6000)=[];
disp('bval removed 6000')

%%
bval_reduced(idx_200)=[];
disp('bval removed 6000')


%informal testing
%bval_tot(idx_200)
%unique(bval_reduced)
%%
bvec_reduced(:,idx_200)=[];
disp('bvec removed 200')

%%
bvec_reduced(:,idx_6000)=[];
disp('bvec removed 6000')

%%
V(:,:,:,idx_200)=[];
tot_fourth_length_reduced200=length(V(1,1,1,:));
disp('Vols removed 200')

%%
V(:,:,:,idx_6000)=[];
tot_fourth_length_reduced6000=length(V(1,1,1,:));
disp('Vols removed 6000')

%%

hdr = niftiinfo(img_path);
hdr.Datatype = 'double';

subj='pil002';
ses='resting';


hdr.ImageSize = size(V);
niftiwrite(V,strcat('/media/nas_rete/GLOVE_STUDY/DDC/derivatives/pil002/resting/dwi/',subj,'_',ses,'_dwi_reduced200.nii'),hdr,"Compressed",true);


