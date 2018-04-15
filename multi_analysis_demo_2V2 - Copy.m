% multi_analysis_demo.m
% Corban Swain, 2018
%%
close all

% get the list of files
multiAnalysisDir = fullfile('data', 'matlab_io', 'multi_analysis');
fileListPath = fullfile(multiAnalysisDir, 'file_list.mat');
varNames1 = {'fileNames'};
load(fileListPath, varNames1{:});

% go through each file
nFiles = length(fileNames);
for i = 1:nFiles
    dataFilePath = fullfile(multiAnalysisDir, fileNames{i});
    load(dataFilePath)
    
    % INFO Printing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    % ANALYSIS HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sampleGroupNumbers = sampleGroupNumbers + 1;
    values = zscore(values');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Name',[fileNames{i} 'Explained Varience'],'NumberTitle','off','Color','w');
    [coeff,score,latent,tsquared,explained] = pca(values);
    bar(explained)
    axis([0 10 0 100])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Name',[fileNames{i} 'Bi-Plot'],'NumberTitle','off','Color','w');
    hbi = biplot(coeff(:,1:2),'scores',score(:,1:2),...
        'ObsLabels',num2str((1:size(score,1))'));
    
    cval = hex2rgb(['#e6194b';'#3cb44b';'#ffe119';'#f58231';...
        '#911eb4';'#46f0f0';'#fabebe';'#008080';...
        '#aa6e28';'#aaffc3';'#000080';'#000000']);
    
    ref = get(hbi(:),'Tag');
    ind = find(contains(ref,'obsmarker'));
    Norm_ind = find(contains(sampleNames,'Normal'));
    Can_ind = find(~contains(sampleNames,'Normal'));
    
    if ~isempty(Norm_ind)
    set(hbi(Norm_ind),'Marker', 'p', 'MarkerSize', 12);
    cval = repelem(cval,2,1);
    end 
    
    set(hbi(Can_ind),'Marker','.','MarkerSize',12);
    
    for ii = 1:length(ind)
        set(hbi(ind(ii)),'Color',cval(sampleGroupNumbers(ii),:));
    end
    
        values = zscore(values');
    [coeff,score,latent] = pca(values,'NumComponents',2,'Algorithm','als');
    
    opts = statset('Display','off');
    clusterNum = 12;
    [idscore,C] = kmeans(score,clusterNum,'Distance','correlation',...
        'Replicates',50,'Options',opts);
    
    figure('Name',[fileNames{i} ' K-means plot'],'NumberTitle','off');
    x1 = min(score(:,1)):0.01:max(score(:,1));
    x2 = min(score(:,2)):0.01:max(score(:,2));
    [x1G,x2G] = meshgrid(x1,x2);
    XGrid = [x1G(:),x2G(:)];
    idx2Region = kmeans(XGrid,clusterNum,'MaxIter',1,'Start',C);
    gscatter(XGrid(:,1),XGrid(:,2),idx2Region,hsv(clusterNum),'..');
    cval = hsv(clusterNum);
    for v = 1:clusterNum
        scatter(score(idscore==v,1),score(idscore==v,2),12,cval(v,:),'filled')
        hold on
    end
    plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3)
    hold off
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %    figValues = figValues + 1;
    %
    %    figure(figValues)
    %    cval = hsv(length(sampleList));
    %    for i = 1:length(sampleList)
    %        class = find(sampleGroupNumbers == i);
    %        scatter3(score(class,1),score(class,2),score(class,3),10,'MarkerFaceColor',cval(i,:))
    %        hold on
    %    end
    %    legend(sampleList);
    %    hold off
    %    axis equal
    %    view(45,45)
    %    figValues = figValues + 1;
    
   figure('Name',[fileNames{i} 'Gscatter Plot'],'NumberTitle','off','Color','w');
   gscatter(score(:,1),score(:,2),sampleNames',cval,'.',20)
%     for i = 1:length(sampleNames)
%         %        class = find(sampleGroupNumbers == i);
%         %        scatter(score(class,1),score(class,2),10,'MarkerFaceColor',cval(i,:))
%         %        hold on
%     end
    %    legend(sampleNames);
    %    axis equal
    %    figValues = figValues + 1;
    
    %
    %    figure(4)
    %    clustergram(data,'Cluster','column','Colormap','jet','Linkage','complete','Dendrogram',6)
    % Pairwise t-test liked no filter
    % make gscatter plots with different symbol healthy cells
    % k-means
    % choose a centriod pt, dot product all samples with that piont then choose
    % the samples with the max inner dot product, choose the pca that those
    % samples are loaded on most strongly, then choose the genes that
    % contribute most to those components
    % MA plot to PCA space
    % gene space clustering with labeled significance and gene expression
    % clustering
    
    
    
    
end
