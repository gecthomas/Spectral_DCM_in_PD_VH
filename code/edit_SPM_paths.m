%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script created by George Thomas 05/24
%
% For SPM practical session May 2024
%
% Run this script once after downloading and unpacking the data. It will
% edit the metadata (namely paths to files) within subjects' SPM.mat files
% so that they match up with the local directory.
%
% Subsequent scripts will not work correctly unless this is done.
%
% The input should be the full path to the local DCM batch directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edit_SPM_paths(batchDir)

subjects = importdata('VH_non_VH_list.txt');

for ii = 1:length(subjects)
   
    thisSub = subjects{ii};
    disp(thisSub);
    subGlmDir = fullfile(batchDir, '..', 'GLM', thisSub);
    subFuncDir = fullfile(batchDir, '..', 'func', thisSub);
    
    load(fullfile(subGlmDir, 'SPM.mat'),'SPM');
    
    SPM.swd = subGlmDir;
    
    nscan = SPM.nscan;
    
    SPM.xY = rmfield(SPM.xY,'P');
    
    for jj = 1:nscan
    
        SPM.xY.P(jj,:) = [fullfile(subFuncDir, sprintf('rsfmri_ica-aroma%04d.nii',jj-1)) ',1'];
        SPM.xY.VY(jj).fname = fullfile(subFuncDir, sprintf('rsfmri_ica-aroma%04d.nii',jj-1));
    
    end
    
    save(fullfile(subGlmDir, 'SPM.mat'),'SPM');
    
end

end