% Analysis_Of_Cancer_DiffGenes.mat
% Tomer Zohar, 2018

%% Structure of 'part_1_analysis_v1.1.mat' Data File
clear all 
close all 
clc

var_names = {'S', 'nGenes', 'geneIDs'};
filePath = 'data/matlab_io/part_1_analysis_v1.2.mat';
load(filePath, var_names{:});
can_type = fieldnames(S);

for i = 1:length(geneIDs)
    v_num = find(geneIDs{i} == '.');
    geneIDs{i}(v_num:end) = [];
end

val_ref = struct2cell(S); 
isSignif_mat = val_ref{1}.isSignif; 
for i = 2:length(can_type)
    isSignif_mat = horzcat(isSignif_mat,val_ref{i}.isSignif);
end

isSignif_Uni = zeros(length(isSignif_mat),1);
for i = 1:length(isSignif_mat)
    if nnz(isSignif_mat(i,:)) > 0
        isSignif_Uni(i) = 1;
    else
        isSignif_Uni(i) = 0;
    end
end

%log2fc matrix
isSignif_Um = repmat(isSignif_Uni,1,length(can_type));
logfc_mat = val_ref{1}.logfc; 
for i = 2:length(can_type)
    logfc_mat = horzcat(logfc_mat,val_ref{i}.logfc);
end
Sign_logfc = (logfc_mat).*(isSignif_Um);

% Delete all rows with zeros
Sign_logfc = Sign_logfc(any(Sign_logfc,2),:);

clustergram(Sign_logfc,'Cluster','column','ColumnLabels',can_type','Colormap','parula','Linkage','complete','Dendrogram',6)

[coeff,score,latent,tsquared,explained] = pca(Sign_logfc);
figure(1)
bar(explained)

[coeff,score,latent,tsquared,explained] = pca(Sign_logfc,'NumComponents',2);
figure(2) 
biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'BLCA','BRCA','HNSC','KICH','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','STAD','THCA'});

figure(3)
rng default % for reproducibility
mapped = tsne(Sign_logfc);
gscatter(mapped(:,1),mapped(:,2));

figure(4)
[cidx,ctrs] = kmeans(Sign_logfc,10,'dist','corr','rep',5,'disp','final');
hold on
for k = 1:2
    scatter(1:length(Sign_logfc),Sign_logfc(:,k),'MarkerFaceColor',[rand(1,3)])
end
hold off
