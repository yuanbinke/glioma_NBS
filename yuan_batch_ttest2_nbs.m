clc
clear
PorQ=0.01;
CPThrd=0.05;
M=1000;
covdir='I:\ZhangNan_Pre';
cd(covdir)
load HGG_4_26_2.txt
load LGG77_cov_2.txt
data1dir='I:\ZhangNan_Pre\FC_LGG_HGG\FC_LGG\subatlas_LP_BNSL080\GSRFB';
data2dir='I:\ZhangNan_Pre\FC_LGG_HGG\FC_HGG4\subatlas_LP_BNSL080\GSRFB';
OutputDir='I:\ZhangNan_Pre\FC_LGG_HGG\FC_HGG4\t2_BNSL080_001_GSR_LGG_HGG4_code_nocov';mkdir(OutputDir)
cd(data1dir)
sublist=dir([data1dir filesep '*.txt']);
for s=1:length(sublist)
    Group1Cells(s)={[data1dir filesep sublist(s).name]};
end
[MatrixGroup1, AliasList1]=yuan_GetGroupData(Group1Cells');
% group2
cd(data2dir)
sublist=dir([data2dir filesep '*.txt']);
for s=1:length(sublist)
    Group2Cells(s)={[data2dir filesep sublist(s).name]};
end
[MatrixGroup2, AliasList2]=yuan_GetGroupData(Group2Cells');

fprintf('Group2:\n');
for i=1:numel(AliasList2)
    fprintf('\t%s\n', AliasList2{i});
end
AllMatrix1=zeros([size(MatrixGroup1{1}),numel(MatrixGroup1)]);
AllMatrix2=zeros([size(MatrixGroup2{1}),numel(MatrixGroup2)]);
%Group1
for i=1:numel(MatrixGroup1)
    Matrix=MatrixGroup1{i};
    Matrix=Matrix-diag(diag(Matrix));
    AllMatrix1(:,:,i)=Matrix;
end
%Group2
for i=1:numel(MatrixGroup2)
    Matrix=MatrixGroup2{i};
    Matrix=Matrix-diag(diag(Matrix));
    AllMatrix2(:,:,i)=Matrix;
end
[n11, n12, n13]=size(AllMatrix1);
[n21, n22, n23]=size(AllMatrix2);
%Group1
AllMatrix1=reshape(AllMatrix1, [], n13);
MIndex=triu(true(n11, n12), 1);
AllMatrix1=AllMatrix1(MIndex(:), :);
%Group2
AllMatrix2=reshape(AllMatrix2, [], n23);
MIndex=triu(true(n21, n22), 1);
AllMatrix2=AllMatrix2(MIndex(:), :);

GroupMatrix=cell(2, 1);
GroupMatrix{1, 1}=AllMatrix1';
GroupMatrix{2, 1}=AllMatrix2';

% CovCells{1}=LGG77_cov_2;
% CovCells{2}=HGG_4_26_2;
CovCells=[];
[T, P]=gretna_TTest2(GroupMatrix, CovCells);
PThrd=PorQ;
Index=find(P<PThrd);
if isempty(Index)
    msgbox('No Edge Left');
    fprintf('\n\tApplying Edge P: Done.\n');
    return
end
TThrd=min(abs(T(Index)));
% NBS
NMsk=true(n11, n12);
TMap=zeros(n11*n12, 1);
TMap(MIndex(:))=T;
TMap=reshape(TMap, [n11, n12]);
TMap=TMap+TMap';
PMap=zeros(n11*n12, 1);
PMap(MIndex(:))=P;
PMap=reshape(PMap, [n11, n12]);
PMap=PMap+PMap';
[Comnet, Comnet_P]=gretna_TTest2_NBS(GroupMatrix,...
    NMsk, CovCells, PThrd, CPThrd, TMap, PMap, M);
if isempty(Comnet)
    msgbox(sprintf('No Comnet Left under Edge P < %g, Comnet P < %g', PThrd, CPThrd));
else
    save(fullfile(OutputDir, [Prefix, '_ComnetMat.mat']), 'Comnet', 'Comnet_P');
    save(fullfile(OutputDir, [Prefix, '_ComnetP.txt']),  'Comnet_P', '-ASCII',  '-DOUBLE', '-TABS');
    for c=1:size(Comnet, 1)
        comnet=Comnet{c, 2};
        save(fullfile(OutputDir, sprintf('%s_%s_Comnet_%s.txt', Prefix, Comnet{c, 1}, Comnet{c, 3})),  'comnet', '-ASCII',  '-DOUBLE', '-TABS');
    end
end
fprintf('\n\tNBS: Done.\n');

TMap=zeros(n11*n12, 1);
PMap=zeros(n11*n12, 1);
TMap(MIndex(:))=T;
PMap(MIndex(:))=P;
TMap=reshape(TMap, [n11, n12]);
PMap=reshape(PMap, [n11, n12]);
TMap=TMap+TMap';
PMap=PMap+PMap';
Prefix='Edge';
save(fullfile(OutputDir, [Prefix, '_TNet.txt']),  'TMap', '-ASCII',  '-DOUBLE', '-TABS');
save(fullfile(OutputDir, [Prefix, '_PNet.txt']),  'PMap', '-ASCII',  '-DOUBLE', '-TABS');
save(fullfile(OutputDir, [Prefix, '_TThrd.txt']), 'TThrd', '-ASCII', '-DOUBLE', '-TABS');
save(fullfile(OutputDir, [Prefix, '_PThrd.txt']), 'PThrd', '-ASCII', '-DOUBLE', '-TABS');
fprintf('\nFinished.\n');