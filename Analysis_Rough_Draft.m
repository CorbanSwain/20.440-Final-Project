% Analysis_Rough_Draft.m
% Tomer Zohar
clc
% get the list of files
fileListPath = fullfile('file_list.mat');
load(fileListPath);

% go through each file
nFiles = length(fileNames);
for i = 1
    dataFilePath = fullfile(fileNames{i});
    load(dataFilePath)
    
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
    

    % ANALYSIS HERE
    sampleGroupNumbers = sampleGroupNumbers + 1;
    valuesper = zscore(values,0,2)';
    
    % PCA PLOT
    figure('Name',[fileNames{i} 'gscatter'] ,'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
    [coeff,score1] = pca(valuesper,'Algorithm','eig');
    cval = hex2rgb(['#e6194b';'#3cb44b';'#ffe119';'#f58231';...
        '#911eb4';'#46f0f0';'#fabebe';'#008080';...
        '#aa6e28';'#aaffc3';'#000080';'#000000']);
    zscor1 = zscore(score1,0,2);
    gscatter(zscor1(:,1),zscor1(:,2),sampleNames',cval,'.',25);
    title('Score Plot');xlabel('PCA-1');ylabel('PCA-2');
    axis([-60 60 -40 40])

    % HISTOGRAM COEFF PLOT
    figure('Name',[fileNames{i} 'histogram'],'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600])
    for k = 1:12
        subplot(4,3,k)
        histogram(coeff(:,k))
        title(num2str(k));xlabel('Coefficient Value');ylabel('Occurance');
    end

    % PCA OF JUST SINGIF GENES
    figure('Name',[fileNames{i} 'PCA OF JUST SIGNIF GENES'],'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600])
    valuesSignif = valuesper(:,geneNumSignif > 0);
    [coeff,score2] = pca(valuesSignif,'Algorithm','eig');
    zscor2 = zscore(score2,0,2);
    
    hold on
    for k = 1:12
        row = find(sampleGroupNumbers == k);
        scatter(zscor2(row,1),zscor2(row,2),25,cval(k,:),'filled',...
            'MarkerFaceAlpha',0.9);
    end
    hold off
    legend(sampleList);title('Score Plot');xlabel('PCA-1');ylabel('PCA-2');
    
    %%%%%
    figure('Name',[fileNames{i} 'Cancer PCA 1-3'],'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600])
    hold on
    for k = 1:12
        row = find(sampleGroupNumbers == k);
        scatter(zscor2(row,1),zscor2(row,3),25,cval(k,:),'filled',...
            'MarkerFaceAlpha',0.9);
    end
    hold off
    legend(sampleList);title('Score Plot');xlabel('PCA-1');ylabel('PCA-3');
    
    %%%%%
    figure('Name',[fileNames{i} 'Cancer PCA 1-2-3'],'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600])
    hold on
    for k = 1:12
        row = find(sampleGroupNumbers == k);
        scatter3(zscor2(row,1),zscor2(row,2),zscor2(row,3),25,cval(k,:),'filled',...
            'MarkerFaceAlpha',0.7);
    end
    hold off
    legend(sampleList);title('Score Plot');xlabel('PCA-1');ylabel('PCA-2');zlabel('PCA-3');
    view([45 45])

    % PCA Removals GLOBAL
    figure('Name',[fileNames{i} 'PCA Removals GLOBAL'] ,'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
    sigcol = find(geneNumSignif > 0);
    valsigr = valuesper;
    valsigr(:,sigcol) = [];
    [coeff,score3] = pca(valsigr,'Algorithm','eig');
    zscor3 = zscore(score3,0,2);
    gscatter(zscor3(:,1),zscor3(:,2),sampleNames',cval,'.',25);
    title('Score Plot');xlabel('PCA-1');ylabel('PCA-2');
    axis([-60 60 -40 40])
    
    % PCA Removals GLOBAL RAND
    figure('Name',[fileNames{i} 'PCA Removals GLOBAL RAND'] ,'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
    samp = datasample(1:length(geneNumSignif),length(find(geneNumSignif > 0)));
    valsigr = valuesper;
    valsigr(:,samp) = [];
    [coeff,score4] = pca(valsigr,'Algorithm','eig');
    zscor4 = zscore(score4,0,2);
    gscatter(zscor4(:,1),zscor4(:,2),sampleNames',cval,'.',25);
    title('Score Plot');xlabel('PCA-1');ylabel('PCA-2');
    axis([-60 60 -40 40])

    % PCA Global Groups
    figure('Name',[fileNames{i} 'PCA Locals'] ,'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
    subGroup = {[5],[1 3 7 8 9],[3 8 11],[2 10],[7 10 11],[3 8 9],[3 7 11],[2 7 10 12],[4 6 12]};
    for k = 1:9
        subtemp = subGroup{k};
        subplot(3,3,k)
        hold on
        for ij = 1:length(subtemp)
            tru = ismember(sampleGroupNumbers,subtemp(ij));
            scatter(zscor1(tru,1),zscor1(tru,2),25,cval(subtemp(ij),:),'filled')
        end
        hold off
        title(['Sub Group ' num2str(k)]);
        legend(sampleList(subGroup{k}))
        axis([-60 60 -40 40])
        if k > 6
            xlabel('PCA-1');ylabel('PCA-2')
        end
    end
    
    figure('Name',[fileNames{i} 'PCA Locals 2'] ,'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
    for k = 1:9
        subtemp = subGroup{k};
        subplot(3,3,k)
        hold on
        for ij = 1:length(subtemp)
            tru = ismember(sampleGroupNumbers,subtemp(ij));
            scatter(zscor1(tru,1),zscor1(tru,3),25,cval(subtemp(ij),:),'filled')
        end
        hold off
        title(['Sub Group ' num2str(k)]);
        legend(sampleList(subGroup{k}))
        axis([-60 60 -40 40])
        if k > 6
            xlabel('PCA-1');ylabel('PCA-3')
        end
    end
    
% RAND EVALUATE REMOVALS
    % PCA Global Groups
    figure('Name',[fileNames{i} 'PCA Locals Removal'] ,'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
    subGroup = {[5],[1 3 7 8 9],[3 8 11],[2 10],[7 10 11],[3 8 9],[3 7 11],[2 7 10 12],[4 6 12]};
    for k = 1:9
        subtemp = subGroup{k};
        subplot(3,3,k)
        hold on
        for ij = 1:length(subtemp)
            tru = ismember(sampleGroupNumbers,subtemp(ij));
            scatter(zscor3(tru,1),zscor3(tru,2),25,cval(subtemp(ij),:),'filled')
        end
        hold off
        title(['Sub Group ' num2str(k)]);
        legend(sampleList(subGroup{k}))
        axis([-60 60 -40 40])
        if k > 6
            xlabel('PCA-1');ylabel('PCA-2')
        end
    end
    
    figure('Name',[fileNames{i} 'PCA Locals 2 removal'] ,'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
    for k = 1:9
        subtemp = subGroup{k};
        subplot(3,3,k)
        hold on
        for ij = 1:length(subtemp)
            tru = ismember(sampleGroupNumbers,subtemp(ij));
            scatter(zscor3(tru,1),zscor3(tru,3),25,cval(subtemp(ij),:),'filled')
        end
        hold off
        title(['Sub Group ' num2str(k)]);
        legend(sampleList(subGroup{k}))
        axis([-60 60 -40 40])
        if k > 6
            xlabel('PCA-1');ylabel('PCA-3')
        end
    end
% RAND EVALUATE REMOVAL RANDOM 
    % PCA Global Groups
    figure('Name',[fileNames{i} 'PCA Locals Removal RAND'] ,'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
    subGroup = {[5],[1 3 7 8 9],[3 8 11],[2 10],[7 10 11],[3 8 9],[3 7 11],[2 7 10 12],[4 6 12]};
    for k = 1:9
        subtemp = subGroup{k};
        subplot(3,3,k)
        hold on
        for ij = 1:length(subtemp)
            tru = ismember(sampleGroupNumbers,subtemp(ij));
            scatter(zscor4(tru,1),zscor4(tru,2),25,cval(subtemp(ij),:),'filled')
        end
        hold off
        title(['Sub Group ' num2str(k)]);
        legend(sampleList(subGroup{k}))
        axis([-60 60 -40 40])
        if k > 6
            xlabel('PCA-1');ylabel('PCA-2')
        end
    end
    
    figure('Name',[fileNames{i} 'PCA Locals 2 removal RAND'] ,'NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
    for k = 1:9
        subtemp = subGroup{k};
        subplot(3,3,k)
        hold on
        for ij = 1:length(subtemp)
            tru = ismember(sampleGroupNumbers,subtemp(ij));
            scatter(zscor4(tru,1),zscor4(tru,3),25,cval(subtemp(ij),:),'filled')
        end
        hold off
        title(['Sub Group ' num2str(k)]);
        legend(sampleList(subGroup{k}))
        axis([-60 60 -40 40])
        if k > 6
            xlabel('PCA-1');ylabel('PCA-3')
        end
    end
    
    s1 = silhouette(zscor1,sampleGroupNumbers);
    s2 = silhouette(zscor2,sampleGroupNumbers);
    s3 = silhouette(zscor3,sampleGroupNumbers);
    s4 = silhouette(zscor4,sampleGroupNumbers);
    
for i =1:12
    v1(i) = mean(s1(sampleGroupNumbers == i))
end 

for i =1:12
    v2(i) = mean(s2(sampleGroupNumbers == i))
end 

for i =1:12
    v3(i) = mean(s3(sampleGroupNumbers == i))
end 

for i =1:12
    v4(i) = mean(s4(sampleGroupNumbers == i))
end

figure('Name','Mean of the Silhouette Scores for Every Cancer Type','NumberTitle','off','Color','w',...
        'rend','painters','pos',[10 10 700 600]);
bg1 = vertcat(v1,v2,v3)';
bar(bg1)
title('Mean of the Silhouette Scores for Every Cancer Type')
legend('No Removal','Removal of Significant Genes','Removal of Random Genes')
set(gca,'XTickLabel',sampleList');
xtickangle(45)

% for jk = 1:12
%     subplot(3,4,jk)
%     bar(s1(sampleGroupNumbers == jk))
%     title(num2str(jk))
% end 


end
