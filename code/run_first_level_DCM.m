%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script created by George Thomas 10/21
%
% This script will specify and estimate first-level dynamic causal models
% (DCMs) for all subjects.
%
% Assumes the following directory structure:
% DCM/
%   batch/
%       <this script>
%       <subject list> .txt file
%   GLM/
%       subject_XXX/
%           SPM.mat
%
% N.B. run from within directory DCM/batch/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Global settings

clear
close all

% get current directory (/path/to/DCM/batch/)
d0 = pwd;

% network name
name = 'VN_LR';

% load subject list
all_subjects = importdata('VH_non_VH_list.txt');


% specify global parameters
% name of main subdirectories relative to d0
glmdir = [d0 '/../GLM/'];
outdir = [d0 '/../analyses/'];

% MRI scanner settings
TR = 3.36;   % Repetition time (secs)
TE = 0.03;  % Echo time (secs)

% Experiment settings
nconditions = 1;

%% Specify DCMs (one per subject)

% regions 
regions = {'LGN_L','LGN_R','MTh_mid','V1_mid',...
    'Hippo_L','Hippo_R','PFCl_L','PFCl_R'};
nregions = length(regions); 

% A-matrix (on / off) 
a = ones(nregions,nregions);
apfx = 'full';

% B-matrix (won't modify anything as resting-state)
b = zeros(nregions); % Rest

% C-matrix (again no modification due to task)
c = zeros(nregions,nconditions);

% D-matrix (disabled)
d = zeros(nregions,nregions,0);

% loop through subjects
for n = 1:length(all_subjects)
    
    sub_name = all_subjects{n};
    sub_glmdir = [glmdir, sub_name];
    % Load SPM
    SPM     = load(fullfile(sub_glmdir,'SPM.mat'));
    SPM     = SPM.SPM;
    
    % Load ROIs
    for r = 1:nregions
        f = fullfile(sub_glmdir, ['VOI_' regions{r} '_1.mat']);
        XY = load(f);
        xY(r) = XY.xY;
    end
    
    % Move to output directory
    cd(sub_glmdir)
    
    % Specify a fully connected DCM
    s = struct();
    s.name       = [apfx '_' name];
    s.u          = 1;                 % Include rest condition
    s.delays     = repmat(TR/2,1,nregions); % Slice timing for each region
    s.TE         = TE;      % echo time
    s.nonlinear  = false;   % bilinear or nonlinear 
    s.two_state  = false;   % one-state or two state
    s.stochastic = false;   % spectral DCM renders model deterministic
    s.centre     = true;    % mean centre input
    s.induced    = 1;       % specifies CSDs rather than timeseries
    s.a          = a;       % conn matrix A (baseline EC)
    s.b          = b;       % conn matrix B (modulation due to input)
    s.c          = c;       % conn matrix C (sensitivity to driving input)
    s.d          = d;
    
    % specify DCM using (SPM.mat,timeseries,forward model)
    DCM = spm_dcm_specify(SPM,xY,s);
    
    % Return to script directory
    cd(d0);
    
end

% Collate into a GCM file and estimate

% Find all fully connected DCM files
dcms = spm_select('FPListRec',glmdir,['DCM_' apfx '_' name '.mat']);
% sort in the same way as the list
[~,b] = ismember(all_subjects, sort(all_subjects));
dcms = dcms(b,:);

% Check if a GCM file already exists in outdir
if exist(fullfile(outdir,[name '_GCM_' apfx '_all.mat']),'file')
    opts.Default = 'No';
    opts.Interpreter = 'none';
    % overwrite it if so
    f = questdlg('Overwrite existing GCM?','Overwrite?','Yes','No',opts);
    tf = strcmp(f,'Yes');
else
    tf = true;
end

% Collate individual DCMs & estimate
if tf
    % Character array -> cell array
    GCM_spec = cellstr(dcms);
    
    % Load all DCM structures
    GCM_spec = spm_dcm_load(GCM_spec);
    
    % give the DCMs names according to subjects (makes it easier to check
    % the group DCMs are correct afterwards)
    mismatch = 0;
    for s = 1:length(GCM_spec)
        if ~isempty(strfind(dcms(s,:),all_subjects{s}))
            GCM_spec{s,1}.name = sprintf(['%s_' apfx], all_subjects{s});
        else 
            mismatch = mismatch + 1;
        end
    end
    % if there are mismatches then the GCM may be badly ordered 
    disp([num2str(mismatch), ' GCM-DCM name mismatches'])
    
    % Estimate DCMs (this won't effect original DCM files)
    use_parfor = true; % switch for parallel processing

    % ESTIMATE DCMS ----------------
    
    % estimate GCM with default priors
    GCM_est_all = spm_dcm_fit(GCM_spec,use_parfor);
    
    % save (has to be called 'GCM' to work with the GUI as well)
    GCM = GCM_est_all;
    save(['../analyses/' name '_GCM_' apfx '_VH.mat'],'GCM');
end

%% delete temp DCM files to keep directory clean
delete( ['sub-*_' apfx '.mat'] )