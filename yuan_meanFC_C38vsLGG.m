clc
clear
Cfcdir='I:\ZhangNan_Pre\FC_LGG_HGG\FC_C38\BNSL080_FunImgARWSDCF_GSRB';
%Group1
cd(Cfcdir)
sublist=dir([Cfcdir filesep '*.txt']);
for s=1:length(sublist)
    Group1Cells(s)={[Cfcdir filesep sublist(s).name]};
end
[MatrixGroup1, AliasList1]=yuan_GetGroupData(Group1Cells');
fprintf('Group1:\n');
for i=1:numel(AliasList1)
    fprintf('\t%s\n', AliasList1{i});
end

AllMatrix1=zeros(numel(MatrixGroup1),1);
% ClockC=zeros(numel(MatrixGroup1),1);
for i=1:numel(MatrixGroup1)
    Matrix=MatrixGroup1{i};
    Matrix=Matrix-diag(diag(Matrix));
    Matrix(Matrix<0.2)=0;
    AllMatrix1(i)=sum(sum(Matrix));
end

Pfcdir='I:\ZhangNan_Pre\FC_LGG_HGG\FC_LGG\subatlas_LP_BNSL080\GSRFB';
%Group2
cd(Pfcdir)
sublist=dir([Pfcdir filesep '*.txt']);
for s=1:length(sublist)
    Group2Cells(s)={[Pfcdir filesep sublist(s).name]};
end
[MatrixGroup2, AliasList2]=yuan_GetGroupData(Group2Cells');
fprintf('Group2:\n');
for i=1:numel(AliasList2)
    fprintf('\t%s\n', AliasList2{i});
end

AllMatrix2=zeros(numel(MatrixGroup2),1);
for i=1:numel(MatrixGroup2)
    Matrix=MatrixGroup2{i};
    Matrix=Matrix-diag(diag(Matrix));
    Matrix(Matrix<0.2)=0;
    AllMatrix2(i)=sum(sum(Matrix));
end

Pfcdir='I:\ZhangNan_Pre\FC_LGG_HGG\FC_HGG\subatlas_LP_BNSL080\GSRFB';
%Group3
cd(Pfcdir)
sublist=dir([Pfcdir filesep '*.txt']);
for s=1:length(sublist)
    Group3Cells(s)={[Pfcdir filesep sublist(s).name]};
end
[MatrixGroup3, AliasList3]=yuan_GetGroupData(Group3Cells');
fprintf('Group3:\n');
for i=1:numel(AliasList3)
    fprintf('\t%s\n', AliasList3{i});
end

AllMatrix3=zeros(numel(MatrixGroup3),1);
for i=1:numel(MatrixGroup3)
    Matrix=MatrixGroup3{i};
    Matrix=Matrix-diag(diag(Matrix));
    Matrix(Matrix<0.2)=0;
    AllMatrix3(i)=sum(sum(Matrix));
end
% [h,p,ci,stats] = ttest2(AllMatrix1,AllMatrix2);
% [h1,p1,ci1,stats1] = ttest2(AllMatrix1,AllMatrix3);
% [h2,p2,ci2,stats2] = ttest2(AllMatrix2,AllMatrix3);
AllMatrix1=AllMatrix1/2;
AllMatrix2=AllMatrix2/2;
AllMatrix3=AllMatrix3/2;

behadir='I:\ZhangNan_Pre';
cd(behadir)
load Beha_HGG_53.mat
load Beha_LGG_77.mat
load C38_cov_FD_2.txt
%one-way ANCOVA
DependentVariable=[AllMatrix1;AllMatrix2;AllMatrix3];
GroupLabel=[ones(length(AllMatrix1),1);ones(length(AllMatrix2),1)*2;ones(length(AllMatrix3),1)*3];
Covariates=[C38_cov_FD_2;Beha_LGG_77(:,[1 2 3 6 8]);Beha_HGG_53(:,[1 2 3 6 8])];
[GCF, GCP]=y_ancova1(DependentVariable,GroupLabel,Covariates);
[GF, GP]=y_ancova1(DependentVariable,GroupLabel);
%% ttest2 with cov
group_A=AllMatrix1;
group_B=AllMatrix2;
Cov_group_A=C38_cov_FD_2;
Cov_group_B=Beha_LGG_77(:,[1 2 3 6 8]);
DependentVariable=[group_A;group_B];
%定义协变量
grouplable=[ones(size(group_A))*1;ones(size(group_B))*(-1)];
Regressors = [grouplable,ones((size(group_A,1)+size(group_B,1)),1),];
CovVariable=[Cov_group_A;Cov_group_B];% Cov_group_A和Cov_group_A需要事先load,每一组的协变量可以是多列
Regressors=[Regressors,CovVariable];
Contrast = zeros(1,size(Regressors,2));
Contrast(1) = 1;
TF_Flag='T';
%ss
[b,r,SSE,SSR, T12, TF_ForContrast]=y_regress_ss(DependentVariable,Regressors,Contrast,TF_Flag);
%如果是双尾检验
p12_cov = 2 * tcdf(-abs(T12(1)),(size(group_A,1)+size(group_B,1))-2);
%
group_A=AllMatrix1;
group_B=AllMatrix3;
Cov_group_A=C38_cov_FD_2;
Cov_group_B=Beha_HGG_53(:,[1 2 3 6 8]);
DependentVariable=[group_A;group_B];
%定义协变量
grouplable=[ones(size(group_A))*1;ones(size(group_B))*(-1)];
Regressors = [grouplable,ones((size(group_A,1)+size(group_B,1)),1),];
CovVariable=[Cov_group_A;Cov_group_B];% Cov_group_A和Cov_group_A需要事先load,每一组的协变量可以是多列
Regressors=[Regressors,CovVariable];
Contrast = zeros(1,size(Regressors,2));
Contrast(1) = 1;
TF_Flag='T';
%ss
[b,r,SSE,SSR, T13, TF_ForContrast]=y_regress_ss(DependentVariable,Regressors,Contrast,TF_Flag);
%如果是双尾检验
p13_cov = 2 * tcdf(-abs(T13(1)),(size(group_A,1)+size(group_B,1))-2);
%
group_A=AllMatrix2;
group_B=AllMatrix3;
Cov_group_A=Beha_LGG_77(:,[1 2 3 6 8]);
Cov_group_B=Beha_HGG_53(:,[1 2 3 6 8]);
DependentVariable=[group_A;group_B];
%定义协变量
grouplable=[ones(size(group_A))*1;ones(size(group_B))*(-1)];
Regressors = [grouplable,ones((size(group_A,1)+size(group_B,1)),1),];
CovVariable=[Cov_group_A;Cov_group_B];% Cov_group_A和Cov_group_A需要事先load,每一组的协变量可以是多列
Regressors=[Regressors,CovVariable];
Contrast = zeros(1,size(Regressors,2));
Contrast(1) = 1;
TF_Flag='T';
%ss
[b,r,SSE,SSR, T23, TF_ForContrast]=y_regress_ss(DependentVariable,Regressors,Contrast,TF_Flag);
%如果是双尾检验
p23_cov = 2 * tcdf(-abs(T23(1)),(size(group_A,1)+size(group_B,1))-2);

Data{1,1}=AllMatrix1;
Data{1,2}=AllMatrix2;
Data{1,3}=AllMatrix3;


Gname={'HC','LGG','HGG'};
Lname={''};
% yuan_gretna_plot_bar2(Data, Gname, Lname);
yuan_gretna_plot_dot_3groups(Data, Gname, Lname)
hold on
p=roundn(p12_cov,-3);
p1=roundn(p13_cov,-3);
p2=roundn(p23_cov,-3);

set(gca, 'FontName','Arial','FontSize',36,'LineWidth', 3.5);
set(gcf, 'WindowStyle','normal');
set(gca,'tickdir','in');
set(gca,'Box','on');
% xlabel('Naming (res)','FontName','Arial','FontSize',24);
% legend(['P (HC vs HGG) = ' num2str(p1)],'Location','NorthEast');
legend(['P (HC vs LGG) = ' num2str(p) '*'], ['P (HC vs HGG) < 0.001 *'],['P (LGG vs HGG) = ' num2str(p2)],'Location','NorthEast');
% legend(['P (HC vs HGG) = ' num2str(p1)],['P (LGG vs HGG) = ' num2str(p2)],'Location','NorthEast');
% legend(['P (HC vs HGG) = ' num2str(p1)],'Location','NorthEast');
ylabel('Global strength','FontName','Arial','FontSize',40);
% xlim([min(res(:,1)), max(res(:,1))]);
% ylim([200,750]);
ylim([80,580]);
% ylim([500,4500]);
% ylim([300,2300]);
% ylim([min([mean(Data{1,1}),mean(Data{1,2}),mean(Data{1,3})])*0.8, max([mean(Data{1,1})+std(Data{1,1}),mean(Data{1,2})+std(Data{1,1}),mean(Data{1,3})+std(Data{1,3})])*1.1]);
cd('I:\ZhangNan_Pre\FC_LGG_HGG')
filename=['MeanFC_BNSL080_C38_LGG77_HGG53_115_cov'];
print(1,'-dtiff',filename);