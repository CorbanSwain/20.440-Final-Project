% multi_analysis_demo.m
% Corban Swain, 2018

% get the list of files
multiAnalysisDir = fullfile('data', 'matlab_io', 'multi_analysis');
fileListPath = fullfile(multiAnalysisDir, 'file_list.mat');
varNames1 = {'fileNames'};
load(fileListPath, varNames1{:});

% go through each file
nFiles = length(fileNames);
varNames2 = {'values', 'geneNames', 'sampleNames', 'analysisMetadata'};
for i = 1:nFiles 
   dataFilePath = fullfile(multiAnalysisDir ,fileNames{i});
   load(dataFilePath, varNames2{:})
   
   % INFO Printing
   fprintf(['\n%d - ANALYSIS INFO:\n\tFilter by: %s\n\tMetric: %s\n\t', ...
      'Listed Samples: %s\n\n'], ...
      i, ...
      analysisMetadata.filter_method, ...
      analysisMetadata.metric, ...
      analysisMetadata.samples)
   disp('Value Matrix Dimensions:')
   disp(size(values))
   disp('Gene Name List Dimensions:')
   disp(size(geneNames))
   disp('Sample Name List Dimensions:')
   disp(size(sampleNames))
 
end
