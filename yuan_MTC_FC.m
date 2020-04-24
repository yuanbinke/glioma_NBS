function yuan_MTC_FC(subname,subfile,tmpfile,RoiIndex,MTCdir,FCdir)
% read the template file
Vtem = spm_vol(tmpfile);
[Ytem, ~] = spm_read_vols(Vtem);
Ytem(isnan(Ytem)) = 0;
Ytem=round(Ytem);
MNI_coord = cell(length(RoiIndex),1);
for j = 1:length(RoiIndex)
    Region = RoiIndex(j);
    ind = find(Region == Ytem(:));
    
%     if ~isempty(ind)
        [I,J,K] = ind2sub(size(Ytem),ind);
        XYZ = [I J K]';
        XYZ(4,:) = 1;
        MNI_coord{j,1} = XYZ;
%     else
%         error (['There are no voxels in ROI' blanks(1) num2str(RoiIndex(j)) ', please specify ROIs again']);
%     end
end

fprintf('Extracting time series for %s\n', subname);

File_filter='';
    
    cd (subfile)
    File_name = spm_select('List',pwd, ['^' File_filter '.*\.img$']);
    if isempty(File_name)
        File_name = spm_select('List',pwd, ['^' File_filter '.*\.nii$']);
    end
    
    Vin = spm_vol(File_name);
    MTC = zeros(size(Vin,1),length(RoiIndex));
    
    for j = 1:length(RoiIndex)
            VY = spm_get_data(Vin,MNI_coord{j,1});
            MTC(:,j) = mean(VY,2);
    end
MTC(isnan(MTC))=0;    
% % desMTCname=[MTCdir filesep subname '.txt'];
% save(desMTCname, 'MTC','-ASCII');
% fprintf('Extracting time series for %s ...... is done\n', subname);

% generate the FC matrix
    [FC_r, FC_p] = corrcoef(MTC);
    
    FC_r = FC_r - diag(diag(FC_r));
%     FC_p = FC_p - diag(diag(FC_p));
% Fisher z
   FC_Fisher=0.5*log((1+FC_r)./(1-FC_r));
   FC_Fisher(isnan(FC_Fisher))=0;
   FC_Fisher(isinf(FC_Fisher))=0;
%    FC_Fisher(abs(FC_Fisher)>2)=0;
%    FC_Fisher2=reshape(FC_Fisher, [length(RoiIndex)*length(RoiIndex),1]);
%    MIndex=triu(true(length(RoiIndex), length(RoiIndex)), 1);
%    FC_Fisher2=FC_Fisher2(MIndex(:), :);
   cd(FCdir)
%    mkdir(subname)
%    cd([FCdir filesep subname])
%    save FC_Fisher2.mat FC_Fisher2
%    save([subname '_R.txt'], 'FC_r','-ASCII');
   save([subname '_zR.txt'], 'FC_Fisher','-ASCII');
%    save([subname '_p.txt'], 'FC_p','-ASCII');
fprintf('Calculateing Pearson correlation for %s ...... is done\n', subname);  
end