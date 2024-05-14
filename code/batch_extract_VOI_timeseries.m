%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script created by George Thomas 10/21
%
% This scipt will extract timeseries from the following regions in each
% hemisphere:
%
% LGN 
% Hippocampus
% PFC (lateral) 
%
% And the following midline regions:
%
% Medial thalamus
% V1
%
%
% Assumes the following directory structure:
% DCM/
%   batch/
%       <this script>
%       <subject list> .txt file
%   GLM/
%       subject_XXX/
%           SPM.mat
%   func/
%       subject_XXX/
%           rsfmri*1.img
%               .
%               .
%               .
%           rsfmri*N.img
%
% N.B. run from within directory DCM/batch/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear matlabbatch
tic
% get current directory (/path/to/DCM/batch/)
d0 = pwd;

% load subject list
subjects = importdata('VH_non_VH_list.txt');

% specify global parameters
% name of main subdirectories relative to d0 
glmdir = [d0 '/../GLM/'];

% define VOIs (assuming sphere with movement of centre fixed)
% names of VOIs
all_VOInames = {'LGN_L', 'LGN_R','MTh_mid','V1_mid',...
    'Hippo_L','Hippo_R','PFCl_L','PFCl_R'};

% VOI centres [x,y,z]mm
all_VOIcentres = [-22 -29 -4;... LGN_L c.o.g seg
    21 -27 -4;... % LGN_R 
    0 -15 6;... % MTh_mid c.o.g seg
    0 -83 2;... % V1_mid from Razi (2017)
    -25 -15 -20;... % Hippo_L c.o.g. seg
    25 -13 -21;... % Hippo_R 
    -38 33 16;... % PFCl_L from Yeo (2011)
    38 33 16]; % PFCl_R 
% radius in mm for each VOI (LGN / MTh smaller)
all_VOIradii = [4 4 8 10 8 8 8 8];

% initialise SPM
spm('defaults', 'FMRI');

%%

% loop through subjects
for n = 1:length(subjects)
    
    disp(['>>>>>>>>>> ', subjects{n} ' <<<<<<<<<<'])
    
    % specify local (per subject) parameters
    
    % path to the subject's GLM directory
    sub_glmdir = [glmdir, subjects{n}];
    % get the SPM.mat file from GLM estimation
    SPM_file = {strcat(sub_glmdir,'/','SPM.mat')};
    % get the brain mask image to use with ROI 
    bmask_file = {strcat(sub_glmdir,'/','mask.nii',',1')};
    
    % loop through VOIs
    % approx. 2x faster to make a batch with all jobs and run all together, 
    % rather than running lots of 1 batch jobs separately
    % preallocation of matlabbatch has minimal effect on speed
    for v = 1:length(all_VOInames)
        
        % define local (per VOI) parameters
        VOIname = all_VOInames{v};
        VOIcentre = all_VOIcentres(v,:);
        VOIradius = all_VOIradii(v);
        
        % make the SPM job
        % SPM file
        matlabbatch{v}.spm.util.voi.spmmat = SPM_file; 
        % no adjustment
        matlabbatch{v}.spm.util.voi.adjust = NaN;
        matlabbatch{v}.spm.util.voi.session = 1; 
        % name of VOI
        matlabbatch{v}.spm.util.voi.name = VOIname; 
        % centre, radius, movement (i1)
        matlabbatch{v}.spm.util.voi.roi{1}.sphere.centre = VOIcentre; 
        matlabbatch{v}.spm.util.voi.roi{1}.sphere.radius = VOIradius; 
        matlabbatch{v}.spm.util.voi.roi{1}.sphere.move.fixed = 1; 
        % brain mask image (i3)
        matlabbatch{v}.spm.util.voi.roi{2}.mask.image = bmask_file; 
        matlabbatch{v}.spm.util.voi.roi{2}.mask.threshold = 0.25; 
        % use voxels in left + right that are within the mask
        matlabbatch{v}.spm.util.voi.expression = 'i1&i2'; 
        
    end
    
    % run the jobs for this subject
    spm_jobman('run', matlabbatch);
    
end

cd(d0)

toc