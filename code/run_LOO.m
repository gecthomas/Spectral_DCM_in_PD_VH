%% leave-one-out cross validation 

tic 
clear
% define regions
regions = {'LGN_L','LGN_R','MTh_mid','V1_mid','Hippo_L','Hippo_R',...
    'PFCl_L','PFCl_R'};
%% load DM + GCM

% load GCM
GCM_VH = load('../analyses/VN_LR_no_regr_GCM_full_VH.mat');
GCM_VH = GCM_VH.(char(fieldnames(GCM_VH)));

% and design matrix
load('../analyses/dmx_VH.mat');

%% PEB global settings
M = struct();
% M.Q is the choice of precision components to use. By setting this to 
% ’all’, the between-subject variability for each DCM connection will be 
% individually estimated.
M.Q      = 'all';
M.maxit  = 512;

%% run LOO for VH vs non-VH using largest PEB parameter

% find largest PEB parameter
load('../analyses/VN_LR_no_regr_BMA_A_179models.mat');
[VH_maxVal, VH_maxIdx] = max(abs(BMA.Ep(...
    (length(regions)^2)+1:(length(regions)^2)*2)));

% apply analysis specific PEB settings
M.X      = X_VH;
M.Xnames = labels_VH;
% define max EP connection to try and predict
[maxTo,maxFrom] = ind2sub([length(regions) length(regions)],VH_maxIdx);
% run LOO
disp(['run LOO with parameter: ' regions{maxFrom} ' to ' regions{maxTo}])
[qE_MP,qC_MP,Q_MP] = spm_dcm_loo(GCM_VH, M,...
    {['A(' num2str(maxTo) ',' num2str(maxFrom) ')']});
save('../analyses/VN_LR_no_regr_LOO_MP_VH.mat','qE_MP','qC_MP','Q_MP');

%% run LOO with some combinations of parameters

% LGN L   = 1
% LGN R   = 2
% MTh mid = 3
% V1 mid  = 4
% Hippo L = 5
% Hippo R = 6
% PFC L   = 7
% PFC R   = 8

% for the top 5 rather than the top 1
[~,sortEp] = sort(abs(BMA.Ep(...
    (length(regions)^2)+1:(length(regions)^2)*2)),'descend');
top5 = sortEp(1:5);
[top5to,top5from] = ind2sub([length(regions) length(regions)],top5);
% and for the top 10 parameters
top10 = sortEp(1:10);
[top10to,top10from] = ind2sub([length(regions) length(regions)],top10);

% M settings same as before

%% run LOO for top 5 params
disp(['run LOO with parameters: '...
    regions{top5from(1)} ' to ' regions{top5to(1)} ', '...
    regions{top5from(2)} ' to ' regions{top5to(2)} ', '...
    regions{top5from(3)} ' to ' regions{top5to(3)} ', '...
    regions{top5from(4)} ' to ' regions{top5to(4)} ', '...
    regions{top5from(5)} ' to ' regions{top5to(5)}])
[qE_top5P,qC_top5P,Q_top5P] = spm_dcm_loo(GCM_VH, M,...
    {['A(' num2str(top5to(1)) ',' num2str(top5from(1)) ')'],...
     ['A(' num2str(top5to(2)) ',' num2str(top5from(2)) ')'],...
     ['A(' num2str(top5to(3)) ',' num2str(top5from(3)) ')'],...
     ['A(' num2str(top5to(4)) ',' num2str(top5from(4)) ')'],...
     ['A(' num2str(top5to(5)) ',' num2str(top5from(5)) ')']});
save('../analyses/VN_LR_no_regr_LOO_top5P_VH.mat',...
    'qE_top5P','qC_top5P','Q_top5P');

%% run LOO for top 10 params
disp(['run LOO with parameters: '...
    regions{top10from(1)} ' to ' regions{top10to(1)} ', '...
    regions{top10from(2)} ' to ' regions{top10to(2)} ', '...
    regions{top10from(3)} ' to ' regions{top10to(3)} ', '...
    regions{top10from(4)} ' to ' regions{top10to(4)} ', '...
    regions{top10from(5)} ' to ' regions{top10to(5)} ', '...
    regions{top10from(6)} ' to ' regions{top10to(6)} ', '...
    regions{top10from(7)} ' to ' regions{top10to(7)} ', '...
    regions{top10from(8)} ' to ' regions{top10to(8)} ', '...
    regions{top10from(9)} ' to ' regions{top10to(9)} ', '...
    regions{top10from(10)} ' to ' regions{top10to(10)}])
[qE_top10P,qC_top10P,Q_top10P] = spm_dcm_loo(GCM_VH, M,...
    {['A(' num2str(top10to(1)) ',' num2str(top10from(1)) ')'],...
     ['A(' num2str(top10to(2)) ',' num2str(top10from(2)) ')'],...
     ['A(' num2str(top10to(3)) ',' num2str(top10from(3)) ')'],...
     ['A(' num2str(top10to(4)) ',' num2str(top10from(4)) ')'],...
     ['A(' num2str(top10to(5)) ',' num2str(top10from(5)) ')'],...
     ['A(' num2str(top10to(6)) ',' num2str(top10from(6)) ')'],...
     ['A(' num2str(top10to(7)) ',' num2str(top10from(7)) ')'],...
     ['A(' num2str(top10to(8)) ',' num2str(top10from(8)) ')'],...
     ['A(' num2str(top10to(9)) ',' num2str(top10from(9)) ')'],...
     ['A(' num2str(top10to(10)) ',' num2str(top10from(10)) ')']});
save('../analyses/VN_LR_no_regr_LOO_top10P_VH.mat',...
    'qE_top10P','qC_top10P','Q_top10P');

%% get some LOO stats

% load LOOs
mydir = [pwd '\..\analyses\'];
files = fullfile(mydir,'VN_LR_no_regr_LOO*.mat');
theFiles = dir(files);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(mydir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  load(fullFileName);
end

% and dmx
load('../analyses/dmx_VH.mat');

% get subject order
hasVH = X_VH(:,2) > 0;
%k = 1:size(X_VH,1);
%k = [k(~hasVH), k(hasVH)];

% get 90% ci
ci = 0.90;
ci = 1 - (1-ci)/2; 
ci = spm_invNcdf(ci);

%% for MP
disp(' --- Using_MP --- ')

% subjects within 90% ci
c_MP  = ci*sqrt(qC_MP);
lower_MP  = qE_MP - c_MP;
higher_MP = qE_MP + c_MP;
ci90_MP = X_VH(:,2)' >= lower_MP & X_VH(:,2)' <= higher_MP;
disp('>>> within 90% CI:')
disp(['PD-no-VH ' num2str(sum(ci90_MP & ~hasVH')) '/' num2str(sum(~hasVH)) ])
disp(['PD-VH ' num2str(sum(ci90_MP & hasVH')) '/' num2str(sum(hasVH)) ])

% correctly predicted
disp('>>> correctly predicted subjects')
disp('PD-no-VH:')
disp(['correct ' num2str(sum(Q_MP(1,~hasVH) > 0.95)) '/' num2str(sum(~hasVH))])
disp(['incorrect ' num2str(sum(Q_MP(2,~hasVH) > 0.95)) '/' num2str(sum(~hasVH))])
disp('PD-VH:')
disp(['correct ' num2str(sum(Q_MP(2,hasVH) > 0.95)) '/' num2str(sum(hasVH))])
disp(['incorrect ' num2str(sum(Q_MP(1,hasVH) > 0.95)) '/' num2str(sum(hasVH))])

% get correlations between predicted and actual group membership
% N.B. better to use point biserial correlation
[~,df_MP] = spm_ancova(X_VH(:,1:2),[],qE_MP(:),[0;1]);
[r_MP,~,p_MP,~] = pointbiserial(double(hasVH),full(qE_MP),0.05,'both');

disp('>>> correlation between predicted and actual group effect')
fprintf('corr(df:%-2.0f) = %-0.2f: p = %-0.5f',df_MP(2),r_MP,p_MP)
disp(' ')

%% for top5P
disp(' --- Using_top5P --- ')

% subjects within 90% ci
c_top5P  = ci*sqrt(qC_top5P);
lower_top5P  = qE_top5P - c_top5P;
higher_top5P = qE_top5P + c_top5P;
ci90_top5P = X_VH(:,2)' >= lower_top5P & X_VH(:,2)' <= higher_top5P;
disp('>>> within 90% CI:')
disp(['PD-no-VH ' num2str(sum(ci90_top5P & ~hasVH')) '/' num2str(sum(~hasVH)) ])
disp(['PD-VH ' num2str(sum(ci90_top5P & hasVH')) '/' num2str(sum(hasVH)) ])

disp('>>> correctly predicted')
disp('PD-no-VH:')
disp(['correct ' num2str(sum(Q_top5P(1,~hasVH) > 0.95)) '/' num2str(sum(~hasVH))])
disp(['incorrect ' num2str(sum(Q_top5P(2,~hasVH) > 0.95)) '/' num2str(sum(~hasVH))])
disp('PD-VH:')
disp(['correct ' num2str(sum(Q_top5P(2,hasVH) > 0.95)) '/' num2str(sum(hasVH))])
disp(['incorrect ' num2str(sum(Q_top5P(1,hasVH) > 0.95)) '/' num2str(sum(hasVH))])

% get correlations between predicted and actual group membership
[~,df_top5P] = spm_ancova(X_VH(:,1:2),[],qE_top5P(:),[0;1]);
[r_top5P,~,p_top5P,~] = pointbiserial(double(hasVH),full(qE_top5P),0.05,'both');

disp('>>> correlation between predicted and actual group effect')
fprintf('corr(df:%-2.0f) = %-0.2f: p = %-0.5f',df_top5P(2),r_top5P,p_top5P)
disp(' ')

%% for top10P
disp(' --- Using_top10P --- ')

% subjects within 90% ci
c_top10P  = ci*sqrt(qC_top10P);
lower_top10P  = qE_top10P - c_top10P;
higher_top10P = qE_top10P + c_top10P;
ci90_top10P = X_VH(:,2)' >= lower_top10P & X_VH(:,2)' <= higher_top10P;
disp('>>> within 90% CI:')
disp(['PD-no-VH ' num2str(sum(ci90_top10P & ~hasVH')) '/' num2str(sum(~hasVH)) ])
disp(['PD-VH ' num2str(sum(ci90_top10P & hasVH')) '/' num2str(sum(hasVH)) ])

disp('>>> correctly predicted')
disp('PD-no-VH:')
disp(['correct ' num2str(sum(Q_top10P(1,~hasVH) > 0.95)) '/' num2str(sum(~hasVH))])
disp(['incorrect ' num2str(sum(Q_top10P(2,~hasVH) > 0.95)) '/' num2str(sum(~hasVH))])
disp('PD-VH:')
disp(['correct ' num2str(sum(Q_top10P(2,hasVH) > 0.95)) '/' num2str(sum(hasVH))])
disp(['incorrect ' num2str(sum(Q_top10P(1,hasVH) > 0.95)) '/' num2str(sum(hasVH))])

% correlation between predicted and actual group membership
[~,df_top10P] = spm_ancova(X_VH(:,1:2),[],qE_top10P(:),[0;1]);
[r_top10P,~,p_top10P,~] = pointbiserial(double(hasVH),full(qE_top10P),0.05,'both');

disp('>>> correlation between predicted and actual group effect')
fprintf('corr(df:%-2.0f) = %-0.2f: p = %-0.5f',df_top10P(2),r_top10P,p_top10P)
disp(' ')
%%

toc