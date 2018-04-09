% Analysis_Of_Cancer_DiffGenes.m
% Tomer Zohar, 2018

%% Structure of 'part_1_analysis_v1.1.mat' Data File
clear all 
close all 
clc

% var_names = {'S', 'nGenes', 'geneIDs'};
filePath = 'data/matlab_io/part_1_analysis_v2.0f.mat';
load(filePath);
can_type = fieldnames(S);

for i = 1:length(geneIDs)
    v_num = find(geneIDs{i} == '.');
    geneIDs{i}(v_num:end) = [];
end

figure(1)
rng('default') 
mapped = tsne(rot90(combData),'Perplexity',10,'Algorithm','exact','Distance','euclidean');
gscatter(mapped(:,1),mapped(:,2),combLabels',hsv(12),'.',3);
disp('done')
