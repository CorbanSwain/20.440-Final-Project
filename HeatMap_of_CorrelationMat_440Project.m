% Analysis_Rough_Draft.m
% Tomer Zohar

clc
close all
% get the list of files
fileListPath = fullfile('file_list.mat');
load(fileListPath);

% go through each file
nFiles = length(fileNames);
for k = 1
    dataFilePath = fullfile(fileNames{k});
    load('no_filter-fold_change_pairwise-tumor_v7.1.mat')
    
        'Listed Samples: %s\n\n'], ...
        k, ...
        analysisMetadata.filter_method, ...
        analysisMetadata.metric, ...
        analysisMetadata.samples)
    disp('Value Matrix Dimensions:')
    disp(size(values))
    disp('Gene Name List Dimensions:')
    disp(size(geneNames))
    disp('Sample Name List Dimensions:')
    disp(size(sampleNames))
    
    cval = hex2rgb(['#e6194b';'#3cb44b';'#ffe119';'#f58231';...
        '#911eb4';'#46f0f0';'#fabebe';'#008080';...
        '#aa6e28';'#aaffc3']);
    cmap = ones(100,3);
    cmap2 = ones(100,3);
    cmap2(51:end,:) = repelem(hex2rgb(['#AAABB8';'#2E9CCA';'#29648A';'#464866';...
        '#25274D']),10,1);
    cmap(51:end,:) = flip(repelem(hex2rgb(['#78281F';'#943126';'#B03A2E';'#CB4335';'#E74C3C';...
        '#EC7063';'#F1948A';'#F5B7B1';'#FADBD8';'#FDEDEC']),5,1));
    cmap(1:50,:) = flip(repelem(hex2rgb(['#EBEDEF';'#D6DBDF';'#AEB6BF';'#85929E';'#5D6D7E';...
        '#34495E';'#2E4053';'#283747';'#212F3C';'#1B2631']),5,1));
    
    % ANALYSIS HERE
    sampleGroupNumbers = sampleGroupNumbers + 1;
    valuesper = zscore(values,0,2);
    SMx = reshape(randperm(numel(values)), size(values,1), size(values,2));
    shufVal = zscore(values(SMx),0,2);
    
    for i = 1:length(sampleList)
        gidx(i) = max(find(sampleGroupNumbers == i));
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%% RAW COR MAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     rm = struct('GroupNumber',gidx,'Annotation',sampleList,...
%      'Color',cval);
%     cm = struct('GroupNumber',{4,5},'Annotation',{'Time1','Time2'},...
%      'Color',{[1 1 0],[0.6 0.6 1]});
 
    corMat = corrcoef(values);
    h0 = HeatMap(corMat,'Colormap','jet','Title','All Cancers lncRNA Fold Change Correlation Coefficient');
    set(h0,'Colormap',cmap)
    
    %%%%%%%%%%%%%%%%%%%%%%%% SIGNIF ONLY COR MAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     sigcol = find(geneNumSignif ~= 1);
%     valsigr = values;
%     valsigr(sigcol,:) = [];
%     corMatrem = corrcoef(valsigr);
%     h1 = HeatMap(corMatrem,'Colormap','jet','Title','Overlaping Differentially Expressed Genes Removed');
%     set(h1,'Colormap',cmap)
%     
    
    %%%%%%%%%%%%%%%%%%%%%%%% MEAN COR MAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v = zeros(length(values),length(sampleList));
    for j =1:length(sampleList)
        v(:,j) = mean(values(:,sampleGroupNumbers == j),2);
    end
    corMean= corrcoef(v);
    h2 = clustergram(corMean,'Cluster','all','Colormap','hsv',...
        'RowLabels',sampleList,...
        'ColumnLabels',sampleList);
    addTitle(h2,'Correlation Coefficient of Mean lncRNA Fold Change in Cancers');
    set(h2,'Colormap',cmap2);
    
    %%%%%%%%%%%%%%%%%%%%%%%% MEAN COR MAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigcol = find(geneNumSignif == 0);
    valsigr = values;
    valsigr(sigcol,:) = [];
    corMatrem = corrcoef(valsigr);
    h3 = HeatMap(corMatrem,'Colormap','jet','Title','Correlation Coefficient Heat Map of Differentially Expressed Genes Only');
    set(h3,'Colormap',cmap);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%% MEAN SIGNIF ONLY COR MAT %%%%%%%%%%%%%%%%%%%%%%%%%
    v = zeros(size(valsigr,1),length(sampleList));
    for j =1:length(sampleList)
        v(:,j) = mean(valsigr(:,sampleGroupNumbers == j),2);
    end
    corMean= corrcoef(v);
    h4 = clustergram(corMean,'Cluster','all','Colormap','hsv',...
        'RowLabels',sampleList,...
        'ColumnLabels',sampleList);
    addTitle(h4,'Correlation Coefficient Cluster of The Mean of Differentially Expressed lncRNA Expression in Cancers');
    set(h4,'Colormap',cmap2);
    
    
    Y = pdist(corMean);
    Z = linkage(Y);
    svals = zeros(length(sampleList),length(sampleList));
    cidx = zeros(length(sampleList),length(sampleList));
    for k = 1:length(sampleList)
        cidx(k,:) = cluster(Z,'Maxclust',k);
        subplot(4,3,k)
        [s,h] = silhouette(corMean,cidx(k,:));
        for jk = 1:k
            svals(k,jk) = mean(s(cidx(k,:) == jk));
        end
        title(num2str(k))
    end
    
    
    svals(isnan(svals)) = 0;
    figure('Name','silhoette means for clusters')
    hold on
    for k = 2:length(sampleList)
        selec = 1:k;
        y = svals(k,selec);
        x = ones(k, 1) * k;
        scatter(x,y,50,cval(k,:),'filled')
    end
    hold off
    
    figure('Name','silhoutte means')
    for k = 1:length(sampleList)
        selec = 1:k;
        y = svals(k,selec);
        ymean(k) = mean(y);
    end
    bar(1:10,ymean)

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXPORT GENE NAMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    isGroup = zeros(6,502);
    ggSignif = zeros(6,length(geneNames)); 
    for i = 1:6
    GroupCanNames{i} = sampleList(find(cidx(6,:) == i));
    isGroup(i,:) = ismember(sampleGroupNumbers,find(cidx(6,:) == i));
    ggSignif(i,:) = all(signifMatrix(:,find(cidx(6,:) == i)) > 0,2);
    end
    ggSignif = logical(ggSignif);
    isGroup = logical(isGroup);

    
%     
%     fid = fopen(strcat(GroupCanNames{4},'.txt'),'wt');
% fprintf(fid,'%s\n', geneNames(ggSignif()));
% fclose(fid)
%     
    
    %     c = cluster(,'Cutoff',3);
    
    %
    %     s = [1    1    1    1    1    2    2    2    2    3    3    3    5    8   ];
    %     t = [2    3    8    10   11   3    8    10   11   8    9    11   6    9   ];
    %     w = [0.62 0.64 0.59 0.62 0.66 0.60 0.59 0.59 0.57 0.57 0.56 0.64 0.79 0.74];
    %
    %     for k = 1:length(s)
    %         cvalE(k,:) = mean(vertcat(cval(s(k),:),cval(t(k),:)));
    %     end
    %
    %
    %     f = figure('Name','','NumberTitle','off','Color','w',...
    %         'rend','painters','pos',[10 10 700 600]);
    %     G = graph(s,t,w,12);
    %     G.Nodes.Name = sampleList';
    %     h = plot(G,'linewidth',((w.^5)*30),'NodeColor',cval,...
    %         'MarkerSize',10,'EdgeColor',cvalE);
    
    % Set the 'visible' property 'off'
    %     ax = gca;
    %     ax.Visible = 'off';
    
    %     subGroup = {[5 6],[8 9],[2 8],[1 2 3 8 11],[1 3 11],[1 2 3]};
    %     for ij = 1:6
    %         subtemp = subGroup{ij};
    %         SignGeneGroup  = v(:,subtemp);
    %         [y,i] = sort(std(SignGeneGroup,0,2));
    %         GroupSignNames{ij} = geneNames(i(1:10));
    %     end
    
    
    % The axes is removed from PDF file, too.
    %     print(gcf,'-dpdf','-r300');
    
    
end