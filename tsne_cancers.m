% tsne_cancers.m
% Tomer Zohar, 2018

%% TSNE Clustering Across Cancer Types
clear all 
close all 
clc

% var_names = {'S', 'nGenes', 'geneIDs'};
filePath = 'data/matlab_io/combined_analysis_v3.0.mat';
load(filePath);

fh = figure(1);
rng('default') 
mapped = tsne(combData','Perplexity',10,'Algorithm','exact','Distance','euclidean');
gscatter(mapped(:,1),mapped(:,2),combLabels',hsv(12),'.',3);
disp('done')
% savefig(fh, 'figures/tnse_cluster_pooled.fig')
