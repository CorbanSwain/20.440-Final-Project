
clear all
clc

% get the list of files
% load('t_test-fold_change_pairwise-tumor_v5.0.mat');
% 
% 
% % INFO Printing
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf(['\n%d - ANALYSIS INFO:\n\tFilter by: %s\n\tMetric: %s\n\t', ...
%     'Listed Samples: %s\n\n'], ...
%     analysisMetadata.filter_method, ...
%     analysisMetadata.metric, ...
%     analysisMetadata.samples)
% disp('Value Matrix Dimensions:')
% disp(size(values))
% disp('Gene Name List Dimensions:')
% disp(size(geneNames))
% disp('Sample Name List Dimensions:')
% disp(size(sampleNames))
% 
% % ANALYSIS HERE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampleGroupNumbers = sampleGroupNumbers + 1;
% valuesper = zscore(values,0,2);
% FOR PCAs 1 and 2
%figure('Name','PCAs 1-2','NumberTitle','off','Color','w');
% % clf;
% % low = min(geneNumSignif);
% % high = max(geneNumSignif);
% % for j = low:high
% %     values = valuesper(geneNumSignif == j ,:)';
% %     [coeff,score,latent,tsquared,explained] = pca(values);
% %      score = zscore(score,0,2);
% %     cval = hex2rgb(['#e6194b';'#3cb44b';'#ffe119';'#f58231';...
% %         '#911eb4';'#46f0f0';'#fabebe';'#008080';...
% %         '#aa6e28';'#aaffc3';'#000080';'#000000']);
% %     try
% %         subplot(3,4,double(j)+double(1))
% %         gscatter(score(:,1),score(:,2),sampleNames',cval,'.',10);
% %         title(num2str(j))
% %         legend off
% %     catch
% %         print('error')
% %     end
% %     
% % end
% % 
% FOR PCAs 3 and 4
% % figure('Name','PCAs 3-4','NumberTitle','off','Color','w');
% % clf;
% % low = min(geneNumSignif);
% % high = max(geneNumSignif);
% % for j = low:high
% %     values = valuesper(geneNumSignif == j ,:)';
% %     [coeff,score,latent,tsquared,explained] = pca(values);
% %     score = zscore(score,0,2);
% %     cval = hex2rgb(['#e6194b';'#3cb44b';'#ffe119';'#f58231';...
% %         '#911eb4';'#46f0f0';'#fabebe';'#008080';...
% %         '#aa6e28';'#aaffc3';'#000080';'#000000']);
% %     try
% %         subplot(3,4,double(j)+double(1))
% %         gscatter(score(:,3),score(:,4),sampleNames',cval,'.',10);
% %         title(num2str(j))
% %         legend off
% %     catch
% %         print('error')
% %     end
% %     
% % end
% % 
% 
% % FOR PCAs 1 and 3
% figure('Name','PCAs 1-3','NumberTitle','off','Color','w');
% clf;
% low = min(geneNumSignif);
% high = max(geneNumSignif);
% for j = low:high
%     values = valuesper(geneNumSignif == j ,:)';
%     [coeff,score,latent,tsquared,explained] = pca(values);
%     score = zscore(score,0,2);
%     cval = hex2rgb(['#e6194b';'#3cb44b';'#ffe119';'#f58231';...
%         '#911eb4';'#46f0f0';'#fabebe';'#008080';...
%         '#aa6e28';'#aaffc3';'#000080';'#000000']);
%     try
%         subplot(3,4,double(j)+double(1))
%         gscatter(score(:,1),score(:,2),sampleNames',cval,'.',10);
%         title(num2str(j))
%         axis([-5 5 -5 5])
%         legend off
%     catch
%         print('error')
%     end
%     
% end
% 
% %SHUFFLING THE VALUES AND COMPARING
% load('t_test-fold_change_pairwise-tumor_v5.0.mat');
% SMx = reshape(randperm(622*537),622,537);
% values = values(SMx); 
% valuesper = zscore(values,0,2);
% figure('Name','PCAs','NumberTitle','off','Color','w');
% clf;
% low = min(geneNumSignif);
% high = max(geneNumSignif);
% for j = low:high
%     values = valuesper(geneNumSignif == j ,:)';
%     [coeff,score,latent,tsquared,explained] = pca(values);
%     score = zscore(score,0,2);
%     cval = hex2rgb(['#e6194b';'#3cb44b';'#ffe119';'#f58231';...
%         '#911eb4';'#46f0f0';'#fabebe';'#008080';...
%         '#aa6e28';'#aaffc3';'#000080';'#000000']);
%     try
%         subplot(3,4,double(j)+double(1))
%         gscatter(score(:,1),score(:,2),sampleNames',cval,'.',10);
%         title(num2str(j))
%         legend off
%         axis([-5 5 -5 5])
%     catch
%         print('error')
%     end
%     
% end


% FOR REDUCING THE SIGNIF GENES I DONT KNOW IF THIS IS QUITE RIGHT
load('data/matlab_io/multi_analysis/no_filter-fold_change_mean-tumor_v5.0.mat');
valuesper = zscore(values,0,2);
%figure('Name','w/o','NumberTitle','off','Color','w');
figure(1)
clf;

axs = gobjects(10, 1);
num = nnz(geneNumSignif);
rng = linspace(.1,1,10);
rng = round(num.*rng);
for j = 1:10
    val = valuesper';
    elim = datasample(find(geneNumSignif > 0), rng(j), 'Replace', false);
    val(:,elim) = [];
    [coeff,score,latent,tsquared,explained] = pca(val);
%     score = zscore(score,0,2);
    cval = hex2rgb(['#e6194b';'#3cb44b';'#ffe119';'#f58231';...
        '#911eb4';'#46f0f0';'#fabebe';'#008080';...
        '#aa6e28';'#aaffc3';'#000080';'#000000']);
    disp(num2str(length(score(:,1))))
    try
        axs(j) = subplot(3,4,double(j));
        gscatter(score(:,3),score(:,4),sampleNames',cval,'.',10);
        title(['# Signifgenes Reduced ' num2str(rng(j))])
        legend off
        axis([-100 100 -60 60])
    catch
        print('error')
    end
    
end
linkaxes(axs)
% FOR REDUCING THE SIGNIF GENES I DONT KNOW IF THIS IS QUITE RIGHT 
% RANDOM TESTING
load('data/matlab_io/multi_analysis/threshold-fold_change_pairwise-tumor_v5.0.mat');
valuesper = zscore(values,0,2);
%figure('Name','w/o shuffled','NumberTitle','off','Color','w');
figure(2);
clf;

num = nnz(geneNumSignif);
rng = linspace(.1,1,10);
rng = round(num.*rng);
for j = 1:10
    val = valuesper';
    elim = datasample(1:length(val), rng(j), 'replace', false);
    val(:,elim) = [];
    [coeff,score,latent,tsquared,explained] = pca(val);
%     score = zscore(score,0,2);
    cval = hex2rgb(['#e6194b';'#3cb44b';'#ffe119';'#f58231';...
        '#911eb4';'#46f0f0';'#fabebe';'#008080';...
        '#aa6e28';'#aaffc3';'#000080';'#000000']);
    disp(num2str(length(score(:,1))))
    try
        subplot(3,4,double(j))
        gscatter(score(:,1),score(:,2),sampleNames',cval,'.',10);
        title(['# randGene reduced ' num2str(rng(j))])
        legend off
        axis([-100 100 -60 60])
    catch
        print('error')
    end
    
end
