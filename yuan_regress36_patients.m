
clc
clear
FunWdir='I:\ZhangNan_Pre\FunImgARWD';
cd(FunWdir)
Funlist=dir('*p*');
wtumordir='I:\ZhangNan_Pre\wtumormask_3mm';
cd(wtumordir)
wlist=dir('w*.nii');
rpdir='I:\ZhangNan_Pre\RealignParameter';
cd(rpdir)
rplist=dir('*p*');
Nvol=230;
FunCdir='I:\ZhangNan_Pre\FunImgARWD_NoGSR';mkdir(FunCdir)
for i=1:length(Funlist)
    fprintf('RegreeOut %s\n',  Funlist(i).name);
    cd([FunWdir filesep Funlist(i).name])
    Inputtmp=dir('D*.nii');
    Inputfile={[FunWdir filesep Funlist(i).name filesep Inputtmp(1).name]}
    TMMsk=[wtumordir filesep wlist(i).name]
    ttfile=which('dpabi');
    [pathstr, ~] = fileparts(ttfile);
    GSMsk=[];
    WMMsk=fullfile(pathstr,'Templates','WhiteMask_09_61x73x61.hdr');
    CSFMsk=fullfile(pathstr,'Templates','CsfMask_07_61x73x61.hdr');
    HMInd=4;
    cd([rpdir filesep rplist(i).name])
    rptmp=dir('rp*.txt');
    HMFile={[rpdir filesep rplist(i).name filesep rptmp(1).name]};
    Run_RegressOut(Inputfile, TMMsk, GSMsk, WMMsk, CSFMsk, HMInd, HMFile)
    FunCsubdir=[FunCdir filesep Funlist(i).name];mkdir(FunCsubdir)
    cd([FunWdir filesep Funlist(i).name])  
    movefile('cNGS*.nii',FunCsubdir)
%     Run_RegressOut({fullfile('FunImgARW', D{i}, 'wrarest.nii')},...
%         fullfile('T1Img', D{i}, ['wReslice_Tumor_Mask_', D{i}, '.nii']),...
%         [],...
%         'WhiteMask_09_61x73x61.hdr',...
%         'CsfMask_07_61x73x61.hdr',...
%         4, {fullfile('RealignParameter', D{i}, 'rp_arest.txt')});
end


FunCdir='I:\ZhangNan_Pre\FunImgARWD_GSR';mkdir(FunCdir)
for i=84:length(Funlist)
    fprintf('RegreeOut %s\n',  Funlist(i).name);
    cd([FunWdir filesep Funlist(i).name])
    Inputtmp=dir('D*.nii');
    Inputfile={[FunWdir filesep Funlist(i).name filesep Inputtmp(1).name]}
    TMMsk=[wtumordir filesep wlist(i).name]
    ttfile=which('dpabi');
    [pathstr, ~] = fileparts(ttfile);
    GSMsk=fullfile(pathstr,'Templates','BrainMask_05_61x73x61.hdr');
    WMMsk=fullfile(pathstr,'Templates','WhiteMask_09_61x73x61.hdr');
    CSFMsk=fullfile(pathstr,'Templates','CsfMask_07_61x73x61.hdr');
    HMInd=4;
    cd([rpdir filesep rplist(i).name])
    rptmp=dir('rp*.txt');
    HMFile={[rpdir filesep rplist(i).name filesep rptmp(1).name]};
    Run_RegressOut(Inputfile, TMMsk, GSMsk, WMMsk, CSFMsk, HMInd, HMFile)
    FunCsubdir=[FunCdir filesep Funlist(i).name];mkdir(FunCsubdir)
    cd([FunWdir filesep Funlist(i).name])  
    movefile('cWGS*.nii',FunCsubdir)
%     Run_RegressOut({fullfile('FunImgARW', D{i}, 'wrarest.nii')},...
%         fullfile('T1Img', D{i}, ['wReslice_Tumor_Mask_', D{i}, '.nii']),...
%         [],...
%         'WhiteMask_09_61x73x61.hdr',...
%         'CsfMask_07_61x73x61.hdr',...
%         4, {fullfile('RealignParameter', D{i}, 'rp_arest.txt')});
end

