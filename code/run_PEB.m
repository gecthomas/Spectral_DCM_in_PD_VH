clear

%% Load PEB prerequisites

tic;
 
% network name
net = 'VN_LR';

% Load design matrix 

% VH vs non-VH
load('../analyses/dmx_VH.mat');

% Load GCM
GCM_VH = load(['../analyses/' net '_GCM_full_VH.mat']);
GCM_VH = GCM_VH.(char(fieldnames(GCM_VH)));

% PEB global settings
M = struct();
% M.Q is the choice of precision components to use. By setting this to 
% ’all’, the between-subject variability for each DCM connection will be 
% individually estimated.
M.Q      = 'all';
M.maxit  = 128;

%% Build PEB
% Build PEB using parameters from matrix A, the baseline connectivty mtx
% Always the case for rs-fMRI (can ommit 'A' from name for brevity)

% For hallucinators vs non-hallucinators
M.X      = X_VH;
M.Xnames = labels_VH;
[PEB_VH,RCM_VH] = spm_dcm_peb(GCM_VH,M,{'A'});
save(['../analyses/' net '_PEB_VH.mat'],'PEB_VH','RCM_VH');

%% Automatic search over parameters

% Visual hallucinators / non-hallucinators
BMA_auto_VH = spm_dcm_peb_bmc(PEB_VH);
save(['../analyses/' net '_BMA_auto_VH.mat'],'BMA_auto_VH');

toc