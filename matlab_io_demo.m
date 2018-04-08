% matlab_io_demo.mat
% Corban Swain, 2018

%% Structure of 'part_1_analysis_v1.2.mat' Data File

% there are three variables in the data file
var_names = {'S', 'nGenes', 'geneIDs'};
load('data/matlab_io/part_1_analysis_v1.2.mat', var_names{:});

% 'S': a struct of structs with a field for each cancer type, named 
% according to its four letter reference code
fprintf('Field names contained in S (the cancer types):\n')
disp(fieldnames(S));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each field (cancer type) in 'S' contains another struct with th following
% fields:

% 	'tStat': a column vector of with shape [nGenes, 1] for the t statistic
% 	comparing normal and tumor samples

% 	'pValue': a column vector with shape [nGenes, 1] for the p value for
% 	each t statistic

%	'fc': a column vector with shape [nGenes, 1] for the fold change values
%	between the normal and tumor conditions. A value greater than 1
%	represents an increase in the gene for cancer samples.

%	'logfc': a column vector with shape [nGenes, 1] for the log2(fold
%	change) values between the normal and tumor conditions. A value greater
%	than 0 represents an increase in the gene for cancer samples.

%	'isValid': a column vector of logicals with shape [nGenes, 1] of
%	logicals for selecting the genes with analyzable expression levels

%	'isSignif': a column vector of logicals with shape [nGenes, 1] for
%	selecting the genes with significant differential expression between
%	normal and tumor samples

%	'isNormal': a row vector of logicals with shape [1, nSamples] for
%	selecting samples from normal (non-tumor) sources and/or patients

%	'isTumor': a row vector of logicals with shape [1, nSamples] for
%	selecting samples from tumor sources

%	'expression': a matrix with shape [nGenes, nSamples] with the original
%	RPKM expression data

% 	'nSamples': a scalar with the number of samples

%	'sampleNames': a cell array with shape [1, nSamples] containing the
%	names of each sample

fprintf('Field names contained in each S.(''cancer_type''):\n')
disp(fieldnames(S.BLCA));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 'nGenes': a scalar indicating the number of genes, this is the same for 
% all datasets
fprintf('nGenes: ')
disp(nGenes);

% 'geneIDs': a cell array of with shape [nGenes, 1] containing the ID codes
%  for each lncRNA gene
fprintf('One gene ID: ')
disp(geneIDs{1});
