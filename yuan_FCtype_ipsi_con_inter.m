clc
clear
labeldir='G:\Parcellations_atlas\BNSL';
cd(labeldir)
Atmp=fopen('Example_BNSL_80_3mm_name.node');
AAtmp=textscan(Atmp,'%f%f%f%f%f%s');
BN_cor_MNI=[AAtmp{1,1},AAtmp{1,2},AAtmp{1,3}];

behadir='I:\ZhangNan_Pre\FC_LGG_HGG_PCA\T2_C_LGG_BNSL_PCA';
cd(behadir)
load Edge_Pos_Comnet_P_0_000999001.edge
FC_mean=Edge_Pos_Comnet_P_0_000999001;
clockn =0;
NN=80;
PosHemlabel1=zeros(1,4);
for i=1:NN-1
    for j=i+1:NN
        clockn=clockn+1;
        if BN_cor_MNI(i,1)<0 && BN_cor_MNI(j,1)<0 % left
            if FC_mean(i,j)~=0
                PosHemlabel1(1)=PosHemlabel1(1)+1;
            end
        end
        if BN_cor_MNI(i,1)>0 && BN_cor_MNI(j,1)<0 % inter
            if FC_mean(i,j)~=0
                PosHemlabel1(2)=PosHemlabel1(2)+1;
            end
        end
        if BN_cor_MNI(i,1)<0 && BN_cor_MNI(j,1)>0 % inter
            if FC_mean(i,j)~=0
                PosHemlabel1(2)=PosHemlabel1(2)+1;
            end
        end
        if BN_cor_MNI(i,1)>0 && BN_cor_MNI(j,1)>0 % right
            if FC_mean(i,j)~=0
                PosHemlabel1(3)=PosHemlabel1(3)+1;
            end
        end
            if FC_mean(i,j)~=0
                PosHemlabel1(4)=PosHemlabel1(4)+1;
            end
        
    end   
end


cd(behadir)
load Edge_Neg_Comnet_P_0_000999001.edge
FC_mean=Edge_Neg_Comnet_P_0_000999001;
clockn =0;
NN=80;
PosHemlabel2=zeros(1,4);
for i=1:NN-1
    for j=i+1:NN
        clockn=clockn+1;
        if BN_cor_MNI(i,1)<0 && BN_cor_MNI(j,1)<0 % left
            if FC_mean(i,j)~=0
                PosHemlabel2(1)=PosHemlabel2(1)+1;
            end
        end
        if BN_cor_MNI(i,1)>0 && BN_cor_MNI(j,1)<0 % inter
            if FC_mean(i,j)~=0
                PosHemlabel2(2)=PosHemlabel2(2)+1;
            end
        end
        if BN_cor_MNI(i,1)<0 && BN_cor_MNI(j,1)>0 % inter
            if FC_mean(i,j)~=0
                PosHemlabel2(2)=PosHemlabel2(2)+1;
            end
        end
        if BN_cor_MNI(i,1)>0 && BN_cor_MNI(j,1)>0 % right
            if FC_mean(i,j)~=0
                PosHemlabel2(3)=PosHemlabel2(3)+1;
            end
        end
            if FC_mean(i,j)~=0
                PosHemlabel2(4)=PosHemlabel2(4)+1;
            end
        
    end   
end


LGGedge=PosHemlabel1;
HGGedge=PosHemlabel2;
H=bar([LGGedge;HGGedge]');figure(gcf);
set(H(1), 'FaceColor', [1 0 0]); 
% set(H(2), 'FaceColor', [0.7 0 0]); 
set(H(2), 'FaceColor', [0 1 1]); 
set(gcf, 'WindowStyle','normal');
set(gca,'tickdir','in');
set(gca,'Box','on');
ylabel('FC number','FontSize',30);
maxtmp=max(max(PosHemlabel1),max(PosHemlabel2));
axis([0.5 4.5 0 maxtmp*1.05]);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Ipsi', 'Inter','Contra','Total'})
hold on 
% ll=legend('LGG > WHO III','LGG > WHO IV','Location','Best');
ll=legend('HC > WHO IV','HC < WHO IV','Location','NorthWest');
set(ll,'FontName','Arial','FontSize',24)
set(gca, 'FontName','Arial','FontSize',24,'LineWidth', 3);