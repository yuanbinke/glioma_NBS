clc
clear
behadir='L:\ZN_YJ_CP\ZN_YJ_Pre';
cd(behadir)
load YJ_beha.mat
load ZN_beha.mat

ExcSubs_YJ=[6,8,20,31,43];
ExcSubs_ZN=[16,50];
YJ_beha(ExcSubs_YJ,:)=[];
ZN_beha(ExcSubs_ZN,:)=[];

LGGIndex_ZN=find(ZN_beha(:,4)<3);
LGGIndex_YJ=find(YJ_beha(:,4)<3);
HGGIndex_ZN=find(ZN_beha(:,4)>2);
HGGIndex_YJ=find(YJ_beha(:,4)>2);

Beha_LGG_77=[ZN_beha(LGGIndex_ZN,1:15);YJ_beha(LGGIndex_YJ,1:15)];
Beha_HGG_53=[ZN_beha(HGGIndex_ZN,1:15);YJ_beha(HGGIndex_YJ,1:15)];

NN=80;
Pfcdir='I:\ZhangNan_Pre\FC_LGG_HGG\FC_HGG\subatlas_LP_BNSL080\GSRFB';
cd(Pfcdir)
sublist1=dir('pa*.txt');
for s=1:length(sublist1)
    Group1Cells1(s)={[Pfcdir filesep sublist1(s).name]};
end
[MatrixGroup1, AliasList1]=yuan_GetGroupData(Group1Cells1');
fprintf('Group1:\n');
for i=1:numel(AliasList1)
    fprintf('\t%s\n', AliasList1{i});
end
AllMatrix1=zeros([size(MatrixGroup1{1}),numel(MatrixGroup1)]);
for i=1:numel(MatrixGroup1)
    Matrix=MatrixGroup1{i};
    Matrix=Matrix-diag(diag(Matrix));
    AllMatrix1(:,:,i)=Matrix;
end

cd(Pfcdir)
sublist2=dir('tem*.txt');
for s=1:length(sublist2)
    Group1Cells2(s)={[Pfcdir filesep sublist2(s).name]};
end
[MatrixGroup2, AliasList2]=yuan_GetGroupData(Group1Cells2');
fprintf('Group2:\n');
for i=1:numel(AliasList2)
    fprintf('\t%s\n', AliasList2{i});
end
AllMatrix2=zeros([size(MatrixGroup2{1}),numel(MatrixGroup2)]);
for i=1:numel(MatrixGroup2)
    Matrix=MatrixGroup2{i};
    Matrix=Matrix-diag(diag(Matrix));
    AllMatrix2(:,:,i)=Matrix;
end

cd(Pfcdir)
sublist3=dir('YJ_*.txt');
for s=1:length(sublist3)
    Group1Cells3(s)={[Pfcdir filesep sublist3(s).name]};
end
[MatrixGroup3, AliasList3]=yuan_GetGroupData(Group1Cells3');
fprintf('Group3:\n');
for i=1:numel(AliasList3)
    fprintf('\t%s\n', AliasList3{i});
end
AllMatrix3=zeros([size(MatrixGroup3{1}),numel(MatrixGroup3)]);
for i=1:numel(MatrixGroup3)
    Matrix=MatrixGroup3{i};
    Matrix=Matrix-diag(diag(Matrix));
    AllMatrix3(:,:,i)=Matrix;
end
AllMatrix=cat(3,AllMatrix1,AllMatrix2,AllMatrix3);

PartialR=zeros(NN,NN);
PartialP=zeros(NN,NN);
for i=1:NN
    for j=1:NN
        tmp=AllMatrix(i,j,:);
        tmpp=reshape(tmp,[size(tmp,1)*size(tmp,2)*size(tmp,3),1]);
        [PartialR(i,j),PartialP(i,j)]=partialcorr(tmpp,Beha_HGG_53(:,15),Beha_HGG_53(:,[1:6,8]));
    end
end
PartialR(isnan(PartialR))=0;
PartialP(isnan(PartialP))=0;
FCtdir='I:\ZhangNan_Pre\FC_LGG_HGG\FC_HGG\t2_BNSL080_001_GSR';
cd(FCtdir)
load Edge_Pos_Comnet_P_0_000999001.mat
Edge_Pos_Comnet=Edge_Pos_Comnet_P_0_000999001;
Edge_Pos_Comnet(Edge_Pos_Comnet>0)=1;
PartialR_Pos=PartialR;
PartialR_Pos_AQ=PartialR_Pos.*Edge_Pos_Comnet;

Index_R_Pos_001=zeros(2,4);
clock=0;
for i=1:NN-1
    for j=(i+1):NN
        if PartialR_Pos_AQ(i,j)~=0 && PartialP(i,j)<0.01;
        clock=clock+1;
        Index_R_Pos_001(clock,1)=i;
        Index_R_Pos_001(clock,2)=j;
        Index_R_Pos_001(clock,3)=PartialR_Pos_AQ(i,j);
        Index_R_Pos_001(clock,4)=PartialP(i,j);
        end
    end
end

PartialR_Pos_AQ_FDR=PartialR_Pos_AQ;
Pvec1=zeros(10,1);
clock=0;
for i=1:NN
    for j=(i+1):NN
        if PartialR_Pos_AQ(i,j)~=0;
            clock=clock+1;
            Pvec1(clock)=PartialP(i,j);
        end
    end
end
[pID1,pN1] = gretna_FDR(Pvec1,0.05);
if ~isnan(pID1)
    PartialR_Pos_AQ_FDR(PartialP>pID1)=0;
end

FCtdir='I:\ZhangNan_Pre\FC_LGG_HGG\FC_HGG\t2_BNSL080_001_GSR';
cd(FCtdir)
load Edge_Neg_Comnet_P_0_014985.mat
Edge_Neg_Comnet=Edge_Neg_Comnet_P_0_014985;
Edge_Neg_Comnet(Edge_Neg_Comnet<0)=1;
PartialR_Neg=PartialR;
PartialR_Neg_AQ=PartialR_Neg.*Edge_Neg_Comnet;
PartialR_Neg_AQ_FDR=PartialR_Neg_AQ;

Index_R_Neg_001=zeros(2,4);
clock=0;
for i=1:NN-1
    for j=(i+1):NN
        if PartialR_Neg_AQ(i,j)~=0 && PartialP(i,j)<0.01;
        clock=clock+1;
        Index_R_Neg_001(clock,1)=i;
        Index_R_Neg_001(clock,2)=j;
        Index_R_Neg_001(clock,3)=PartialR_Neg_AQ(i,j);
        Index_R_Neg_001(clock,4)=PartialP(i,j);
        end
    end
end

Pvec=zeros(10,1);
clock=0;
for i=1:NN
    for j=(i+1):NN
        if PartialR_Neg_AQ(i,j)~=0;
            clock=clock+1;
            Pvec(clock)=PartialP(i,j);
        end
    end
end
[pID,pN] = gretna_FDR(Pvec,0.05);
if ~isnan(pID)
    PartialR_Neg_AQ_FDR(PartialP>pID)=0;
end

    
