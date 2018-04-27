load('no_filter-fold_change_mean-tumor_v6.2.mat')
load('signif_matrix.mat')
%  
% 
% for i = 1:12727
% if nnz(isSignif(i,:)) < 7
% isSignif(i,:) = 0;
% end
% end
% 
% Sign7 = geneNames(any(isSignif,2));


subGroup = {[3 8 11],[2 10],[7 10 11],[3 8 9],[3 7 11],[4 6 12]};
sampleGroupNumbers = sampleGroupNumbers + 1;

ov3 = isSignif(:,subGroup{1});
ov4 = isSignif(:,subGroup{2});
ov5 = isSignif(:,subGroup{3});
ov6 = isSignif(:,subGroup{4});
ov7 = isSignif(:,subGroup{5});
ov9 = isSignif(:,subGroup{6});

[i3,j] = find(ov3 == 0);
id3 = ~ismember(1:12727,i3);
N3 = geneNames(id3);

[i4,j] = find(ov4 == 0);
id4 = ~ismember(1:12727,i4);
N4 = geneNames(id4);


[i5,j] = find(ov5 == 0);
id5 = ~ismember(1:12727,i5);
N5 = geneNames(id5);

[i6,j] = find(ov6 == 0);
id6 = ~ismember(1:12727,i6);
N6 = geneNames(id6);

[i7,j] = find(ov7 == 0);
id7 = ~ismember(1:12727,i7);
N7 = geneNames(id7);

[i9,j] = find(ov9 == 0);
id9 = ~ismember(1:12727,i9);
N9 = geneNames(id9);


fid = fopen('N3.txt','wt');
fprintf(fid,'%s\n', N3{:});
fclose(fid)
fid = fopen('N4.txt','wt');
fprintf(fid,'%s\n', N4{:});
fclose(fid)
fid = fopen('N5.txt','wt');
fprintf(fid,'%s\n', N5{:});
fclose(fid)
fid = fopen('N6.txt','wt');
fprintf(fid,'%s\n', N6{:});
fclose(fid)
fid = fopen('N7.txt','wt');
fprintf(fid,'%s\n', N7{:});
fclose(fid)
fid = fopen('N9.txt','wt');
fprintf(fid,'%s\n', N9{:});
fclose(fid)
% fid = fopen('Signif7.txt','wt');
% fprintf(fid,'%s\n', Sign7{:});
% fclose(fid)


