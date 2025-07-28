%test on toy dataset, which is easly interpretable

a = zeros(3,8);
b = rand(3,8);
e = rand(3,8);
c = zeros(3,8);
d = rand(3,8);
f = rand(3,8);

lst = [1,2,3,4,5];
a(1,1:length(lst)) = lst;
a(:,4)=4;
a(2,4)=0;
a(2,5)=4;
%a(:,7)=5;%This case is not present in our data (a row with a label appearing more than once)

lst = [10,5,6,7,8];
c(1,1:length(lst)) = lst;
c(:,4)=4;
c(2,4)=0;
c(2,5)=4;
%c(:,7)=5;

%% select matrices

labels_dwi_subjs = a;
medians_dwi_subjs = b;
SE_dwi_subjs = e;
labels_func_subjs = c;
medians_func_subjs = d;
SE_func_subjs = f;

%%

% detect common labels in dwi space across subjs
initial_common = labels_dwi_subjs(1,:);
common = initial_common;
for i = 1:length(labels_dwi_subjs(:,1))
    commonElements = intersect(common, labels_dwi_subjs(i,:));
    common = commonElements;
end

%First, detect indices of common Elements
commonElements(commonElements==0)=[];
lst_idx_matrix_dwi=[];
for row = 1:length(labels_dwi_subjs(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_dwi_subjs(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix_dwi(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median and SE values and labels
medians_dwi_subjs_final=[];
for row = 1:length(medians_dwi_subjs(:,1))
    medians_dwi_subjs_row=medians_dwi_subjs(row,:);
    lst_idx = lst_idx_matrix_dwi(row,:);
    medians_dwi_subjs_row_final=medians_dwi_subjs_row(lst_idx);
    medians_dwi_subjs_final(row,:)=medians_dwi_subjs_row_final;
end

SE_dwi_subjs_final=[];
for row = 1:length(SE_dwi_subjs(:,1))
    SE_dwi_subjs_row=SE_dwi_subjs(row,:);
    lst_idx = lst_idx_matrix_dwi(row,:);
    SE_dwi_subjs_row_final=SE_dwi_subjs_row(lst_idx);
    SE_dwi_subjs_final(row,:)=SE_dwi_subjs_row_final;
end

labels_dwi_subjs_final=[];
for row = 1:length(labels_dwi_subjs(:,1))
    labels_dwi_subjs_row=labels_dwi_subjs(row,:);
    lst_idx = lst_idx_matrix_dwi(row,:);
    labels_dwi_subjs_row_final=labels_dwi_subjs_row(lst_idx);
    labels_dwi_subjs_final(row,:)=labels_dwi_subjs_row_final;
end

%repeat everything for func space

initial_common = labels_func_subjs(1,:);
common = initial_common;
for i = 1:length(labels_func_subjs(:,1))
    commonElements = intersect(common, labels_func_subjs(i,:));
    common = commonElements;
end

%First, detect indices of common Elements
commonElements(commonElements==0)=[];
lst_idx_matrix_func=[];
for row = 1:length(labels_func_subjs(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_func_subjs(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix_func(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median and SE values and labels
medians_func_subjs_final=[];
for row = 1:length(medians_func_subjs(:,1))
    medians_func_subjs_row=medians_func_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    medians_func_subjs_row_final=medians_func_subjs_row(lst_idx);
    medians_func_subjs_final(row,:)=medians_func_subjs_row_final;
end

SE_func_subjs_final=[];
for row = 1:length(SE_func_subjs(:,1))
    SE_func_subjs_row=SE_func_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    SE_func_subjs_row_final=SE_func_subjs_row(lst_idx);
    SE_func_subjs_final(row,:)=SE_func_subjs_row_final;
end

labels_func_subjs_final=[];
for row = 1:length(labels_func_subjs(:,1))
    labels_func_subjs_row=labels_func_subjs(row,:);
    lst_idx = lst_idx_matrix_func(row,:);
    labels_func_subjs_row_final=labels_func_subjs_row(lst_idx);
    labels_func_subjs_final(row,:)=labels_func_subjs_row_final;
end

%%
%then for each subj (row of the two matrices),
%find common elements (regions) between a (dwi) and b (func) space,
%find the idx 
%keep only those common both in DWI and func space
%create a new dwi and func matrix which will have same size.

%%
%select common elements between the two spaces

commonElements_lst=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    commonElements = intersect(labels_dwi_subjs_final(row,:),labels_func_subjs_final(row,:));
    if length(commonElements)>1
        commonElements_lst(1,:)=commonElements;
    else
        commonElements_lst(end+1)=commonElements;
    end
end

commonElements_lst=unique(commonElements_lst);
commonElements(commonElements==0)=[];
%%
%First, detect indices of common Elements

lst_idx_matrix_dwi_spaces=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_dwi_subjs_final(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix_dwi_spaces(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_dwi_subjs_final_spaces=[];
for row = 1:length(medians_dwi_subjs_final(:,1))
    medians_dwi_subjs_row=medians_dwi_subjs_final(row,:);
    lst_idx = lst_idx_matrix_dwi_spaces(row,:);
    medians_dwi_subjs_row_final=medians_dwi_subjs_row(lst_idx);
    medians_dwi_subjs_final_spaces(row,:)=medians_dwi_subjs_row_final;
end

SE_dwi_subjs_final_spaces=[];
for row = 1:length(SE_dwi_subjs_final(:,1))
    SE_dwi_subjs_row=SE_dwi_subjs_final(row,:);
    lst_idx = lst_idx_matrix_dwi_spaces(row,:);
    SE_dwi_subjs_row_final=SE_dwi_subjs_row(lst_idx);
    SE_dwi_subjs_final_spaces(row,:)=SE_dwi_subjs_row_final;
end

labels_dwi_subjs_final_spaces=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    labels_dwi_subjs_row=labels_dwi_subjs_final(row,:);
    lst_idx = lst_idx_matrix_dwi_spaces(row,:);
    labels_dwi_subjs_row_final=labels_dwi_subjs_row(lst_idx);
    labels_dwi_subjs_final_spaces(row,:)=labels_dwi_subjs_row_final;
end

%%
%do the analogous for c and d

%First, detect indices of common Elements

lst_idx_matrix_func_spaces=[];
for row = 1:length(labels_func_subjs_final(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_func_subjs_final(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix_func_spaces(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_func_subjs_final_spaces=[];
for row = 1:length(medians_func_subjs_final(:,1))
    medians_func_subjs_row=medians_func_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    medians_func_subjs_row_final=medians_func_subjs_row(lst_idx);
    medians_func_subjs_final_spaces(row,:)=medians_func_subjs_row_final;
end

SE_func_subjs_final_spaces=[];
for row = 1:length(SE_func_subjs_final(:,1))
    SE_func_subjs_row=SE_func_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    SE_func_subjs_row_final=SE_func_subjs_row(lst_idx);
    SE_func_subjs_final_spaces(row,:)=SE_func_subjs_row_final;
end

labels_func_subjs_final_spaces=[];
for row = 1:length(labels_func_subjs_final(:,1))
    labels_func_subjs_row=labels_func_subjs_final(row,:);
    lst_idx = lst_idx_matrix_func_spaces(row,:);
    labels_func_subjs_row_final=labels_func_subjs_row(lst_idx);
    labels_func_subjs_final_spaces(row,:)=labels_func_subjs_row_final;
end

%% Testing
% % 1. check if the dimensions (number of rows) of the final matrices are
% correct
% % 2. check if the median and SE values in the final matrices correspond
% to the correct labels.
% % 3. check that the labels are unique

% Test 1
assert(length(medians_func_subjs_final_spaces(:,1)) == length(a(:,1)))

% Test 2
idx=find(a==4);
assert(b(idx(1))==medians_dwi_subjs_final_spaces(1))

% Test 3
assert(length(labels_dwi_subjs_final_spaces(1,:))==length(unique(labels_dwi_subjs_final_spaces)))

%% version code in case one label appears more than once in a row 
% (you need to create a cell instead of a matrix)
%otherwise code gives you error as a matrix with rows with different
%lengths cannot be created.

initial_common = labels_dwi_subjs(1,:);
common = initial_common;
for i = 1:length(labels_dwi_subjs(:,1))
    commonElements = intersect(common, labels_dwi_subjs(i,:));
    common = commonElements;
end

%First, detect indices of common Elements
commonElements(commonElements==0)=[];
lst_idx_matrix = cell(1, length(labels_dwi_subjs(:,1)));
for row = 1:length(labels_dwi_subjs(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_dwi_subjs(row,:) == commonElement);
        for j = 1:length(idx)
            lst_idx_row(end+1)=idx(j);
        end
    end
    lst_idx_matrix{1,row}=lst_idx_row;
end

%%
%Secondly, select common indices
%apply this both to the median values and labels
medians_rsoma_subjs_final=cell(1, length(labels_dwi_subjs(:,1)));
for row = 1:length(medians_rsoma_subjs(:,1))
    medians_rsoma_subjs_row=medians_rsoma_subjs(row,:);
    lst_idx = lst_idx_matrix{1,row};
    medians_rsoma_subjs_row_final=medians_rsoma_subjs_row(lst_idx);
    medians_rsoma_subjs_final{1,row}=medians_rsoma_subjs_row_final;
end
%%
labels_dwi_subjs_final=cell(1, length(labels_dwi_subjs(:,1)));
for row = 1:length(labels_dwi_subjs(:,1))
    labels_dwi_subjs_row=labels_dwi_subjs(row,:);
    lst_idx = lst_idx_matrix{1,row};
    labels_dwi_subjs_row_final=labels_dwi_subjs_row(lst_idx);
    labels_dwi_subjs_final{1,row}=labels_dwi_subjs_row_final;
end
%%
%apply everything both to the dwi space and func space, so also c and d

initial_common = labels_func_subjs(1,:);
common = initial_common;
for i = 1:length(labels_func_subjs(:,1))
    commonElements = intersect(common, labels_func_subjs(i,:));
    common = commonElements;
end

%First, detect indices of common Elements
commonElements(commonElements==0)=[];
lst_idx_matrix=[];
for row = 1:length(labels_func_subjs(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_func_subjs(row,:) == commonElement);
        for j = 1:length(idx)
            lst_idx_row(end+1)=idx(j);
        end
    end
    lst_idx_matrix(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_CMRO2_subjs_final=[];
for row = 1:length(medians_CMRO2_subjs(:,1))
    medians_CMRO2_subjs_row=medians_CMRO2_subjs(row,:);
    lst_idx = lst_idx_matrix(row,:);
    medians_CMRO2_subjs_row_final=medians_CMRO2_subjs_row(lst_idx);
    medians_CMRO2_subjs_final(row,:)=medians_CMRO2_subjs_row_final;
end

labels_func_subjs_final=[];
for row = 1:length(labels_func_subjs(:,1))
    labels_func_subjs_row=labels_func_subjs(row,:);
    lst_idx = lst_idx_matrix(row,:);
    labels_func_subjs_row_final=labels_func_subjs_row(lst_idx);
    labels_func_subjs_final(row,:)=labels_func_subjs_row_final;
end

%%
%then for each subj (row of the two matrices),
%find common elements (regions) between a (dwi) and b (func) space,
%find the idx 
%keep only those common both in DWI and func space
%create a new dwi and func matrix which will have same size.
commonElements_lst=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    commonElements = intersect(labels_dwi_subjs_final(row,:),labels_func_subjs_final(row,:));
    if length(commonElements)>1
        commonElements_lst(1,:)=commonElements;
    else
        commonElements_lst(end+1)=commonElements;
    end

end

commonElements_lst=unique(commonElements_lst);
commonElements(commonElements==0)=[];

%First, detect indices of common Elements

lst_idx_matrix=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_dwi_subjs_final(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_rsoma_subjs_final_spaces=[];
for row = 1:length(medians_rsoma_subjs_final(:,1))
    medians_rsoma_subjs_row=medians_rsoma_subjs_final(row,:);
    lst_idx = lst_idx_matrix(row,:);
    medians_rsoma_subjs_row_final=medians_rsoma_subjs_row(lst_idx);
    medians_rsoma_subjs_final_spaces(row,:)=medians_rsoma_subjs_row_final;
end

labels_dwi_subjs_final_spaces=[];
for row = 1:length(labels_dwi_subjs_final(:,1))
    labels_dwi_subjs_row=labels_dwi_subjs_final(row,:);
    lst_idx = lst_idx_matrix(row,:);
    labels_dwi_subjs_row_final=labels_dwi_subjs_row(lst_idx);
    labels_dwi_subjs_final_spaces(row,:)=labels_dwi_subjs_row_final;
end

%do the analogous for c and d

%First, detect indices of common Elements

lst_idx_matrix=[];
for row = 1:length(labels_func_subjs_final(:,1))
    lst_idx_row=[];
    for i = 1:length(commonElements)
        commonElement = commonElements(i);            
        idx = find(labels_func_subjs_final(row,:) == commonElement);
        if length(idx)>1
            lst_idx_row(1,:)=idx;
        else
            lst_idx_row(end+1)=idx;
        end
    end
    lst_idx_matrix(row,:)=lst_idx_row;
end

%Secondly, select common indices
%apply this both to the median values and labels
medians_CMRO2_subjs_final_spaces=[];
for row = 1:length(medians_CMRO2_subjs_final(:,1))
    medians_CMRO2_subjs_row=medians_CMRO2_subjs_final(row,:);
    lst_idx = lst_idx_matrix(row,:);
    medians_CMRO2_subjs_row_final=medians_CMRO2_subjs_row(lst_idx);
    medians_CMRO2_subjs_final_spaces(row,:)=medians_CMRO2_subjs_row_final;
end

labels_func_subjs_final_spaces=[];
for row = 1:length(labels_func_subjs_final(:,1))
    labels_func_subjs_row=labels_func_subjs_final(row,:);
    lst_idx = lst_idx_matrix(row,:);
    labels_func_subjs_row_final=labels_func_subjs_row(lst_idx);
    labels_func_subjs_final_spaces(row,:)=labels_func_subjs_row_final;
end
