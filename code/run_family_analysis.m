% First and second level family analysis for VN_LR_no_regr VH vs non-VH
clear

%% Define hypotheses at first level

% region_names = {'LGN_L','LGN_R','MTh_mid','V1_mid','Hippo_L','Hippo_R',...
%    'PFCl_L','PFCl_R'};

% Define A matrix for top-down vs bottom-up
%--------------------------------------------------------------------------

% Top-down
a_tdbu_fam{1} = [1 1 1 1 1 1 1 1;
                 1 1 1 1 1 1 1 1;
                 0 0 1 0 0 0 1 1;
                 0 0 1 1 1 1 1 1;
                 0 0 1 0 1 1 1 1;
                 0 0 1 0 1 1 1 1;
                 0 0 0 0 0 0 1 1;
                 0 0 0 0 0 0 1 1];
                 
 
% Bottom-up
a_tdbu_fam{2} = a_tdbu_fam{1}';

% Both
a_tdbu_fam{3} = ones(8);

% Names 
a_tdbu_fam_names = {'Top-down','Bottom-up','Both'};

% Define A matrix for inter vs intra-hemispheric
% -------------------------------------------------------------------------

% inter-hemi
a_hemi_fam{1} = [1 1 0 0 0 1 0 1;
                 1 1 0 0 1 0 1 0;
                 0 0 1 0 0 0 0 0;
                 0 0 0 1 0 0 0 0;
                 0 1 0 0 1 1 0 1;
                 1 0 0 0 1 1 1 0;
                 0 1 0 0 0 1 1 1;
                 1 0 0 0 1 0 1 1];
% intra-hemi
a_hemi_fam{2} = [1 0 1 1 1 0 1 0;
                 0 1 1 1 0 1 0 1;
                 1 1 1 1 1 1 1 1;
                 1 1 1 1 1 1 1 1;
                 1 0 1 1 1 0 1 0;
                 0 1 1 1 0 1 0 1;
                 1 0 1 1 1 0 1 0;
                 0 1 1 1 0 1 0 1];

% both
a_hemi_fam{3} = ones(8);

% names
a_hemi_fam_names = {'Interhemi','Intrahemi','Both'};

% Define A matrix for regions
% -------------------------------------------------------------------------

% get full subsets of 5 regional combinations
S5 = dec2bin(0:2^5-1);
S5num = zeros(size(S5));
for ii = 1:numel(S5)
    S5num(ii) = str2num(S5(ii));
end
    
% region numbers
nums = 1:5;
% switch all combinations of regions + add to base matrix
for ii = 1:size(S5num,1)-1
    temp = zeros(5);
    % switch on regions
    temp(:,[logical(S5num(ii+1,:))]) = 1;
    temp([logical(S5num(ii+1,:))],:) = 1;
    % duplicate rows + columns (for L + R)
    temp = repmat(temp, 1, 2);
    temp = reshape(temp', [], 10)';
    temp = kron(temp,ones(1,2));
    % delete duplciate midline regions (MTh and V1)
    temp([4,6],:) = [];
    temp(:,[4,6]) = [];
    % switch on intrinsic connections
    temp(logical(eye(8))) = 1;
    % put into matrix
    a_reg_fam{ii} = temp;
    a_reg_fam_names{ii} = ['regions: '...
        num2str(nums(logical(S5num(ii+1,:))))];
end

% take only unique conn matrices
% (some combos of 4 regions are equivalent to fully connected model)
a_reg_fam_mat = cat(3,a_reg_fam{:});
x = reshape(a_reg_fam_mat,8,[],1);
y = reshape(x(:),8*8,[])';
[~,uni] = unique(y,'rows','stable');
a_reg_fam = a_reg_fam(uni);
a_reg_fam_names = a_reg_fam_names(uni);

for ii = 1:length(a_reg_fam)
    % give the one with all connections on a better name
   if isequal(a_reg_fam{ii}, ones(8))
       a_reg_fam_names{ii} = 'regions: all';
   end
end

%% make a connectivity matrix for each of these combos
% -------------------------------------------------------------------------
m = 1;
for aa = 1:length(a_tdbu_fam)
    for bb = 1:length(a_hemi_fam)
        for cc = 1:length(a_reg_fam)
            all_A{m} = double(a_tdbu_fam{aa} &...
                              a_hemi_fam{bb} &...
                              a_reg_fam{cc});
            % also save name and family assignment
            all_A_names{m} = [a_tdbu_fam_names{aa} ', '...
                              a_hemi_fam_names{bb} ', '...
                              a_reg_fam_names{cc}];
            tdbu_family(m) = aa;
            hemi_family(m) = bb;
            reg_family(m) = cc;
            m = m+1;
        end
    end
end

% take only unique conn matrices 
all_Amat = cat(3,all_A{:});
x = reshape(all_Amat,8,[],1);
y = reshape(x(:),8*8,[])';
[~,uni2] = unique(y,'rows','stable');
all_A = all_A(uni2);
all_A_names = all_A_names(uni2);
tdbu_family = tdbu_family(uni2); 
hemi_family = hemi_family(uni2);
reg_family = reg_family(uni2);

%% define a DCM for each of these models
% -------------------------------------------------------------------------

% load and unpack example DCM
GCM_VH = load('../analyses/VN_LR_no_regr_GCM_full_VH.mat');
GCM_VH = spm_dcm_load(GCM_VH.(char(fieldnames(GCM_VH))));
DCM_template = GCM_VH{1,1};
b = DCM_template.b;
c = DCM_template.c;
d = DCM_template.d;
options = DCM_template.options;

% Output cell array for new models
GCM_templates = {};

for ii = 1:length(all_A)
    
    % fill this empty cell with templates
    DCM = struct();
    DCM.a = all_A{ii};
    DCM.b = b;
    DCM.c = c;
    DCM.d = d;
    DCM.options = options;
    DCM.name = all_A_names{ii};
    GCM_templates{1,ii} = DCM;
    
end

% make a null model with no connections switched on
% -------------------------------------------------------------------------

DCM.a = zeros(8);
DCM.b = b;
DCM.c = c;
DCM.d = d;
DCM.options = options;
DCM.name = 'No connections';


% add in to template space 
GCM_templates{1,length(all_A)+1} = DCM;
tdbu_family(length(all_A)+1) = length(a_tdbu_fam)+1;
hemi_family(length(all_A)+1) = length(a_hemi_fam)+1;
reg_family(length(all_A)+1) = length(a_reg_fam)+1;

% save name
all_A_names{end+1} = 'Null';
a_tdbu_fam_names{end+1} = 'Null';
a_hemi_fam_names{end+1} = 'Null';
a_reg_fam_names{end+1} = 'Null';

% run another check in case any of the templates are the same as the null
% model
remove = [];
for ii = 1:length(GCM_templates) - 1
    if isequal(GCM_templates{ii}.a, zeros(8))
        disp(['template ' num2str(ii) ' is the same as null'])
        remove = [remove, ii];
    elseif isequal(GCM_templates{ii}.a, ones(8))
        disp(['template ' num2str(ii) ' is the full model'])
        all_A_names{ii} = 'Fully connected';
    end
end
all_A_names(remove) = [];
tdbu_family(remove) = [];
hemi_family(remove) = [];
reg_family(remove) = [];
GCM_templates(remove) = [];

save(['../analyses/VN_LR_no_regr_GCM_' num2str(length(GCM_templates)) 'templates.mat'],...
    'GCM_templates','tdbu_family','hemi_family','reg_family');

%% run second level

% load required PEB and design matrix
load('../analyses/VN_LR_no_regr_PEB_VH.mat')
load('../analyses/dmx_VH.mat');

% delete some variables to clear RAM
clear GCM_VH RCM_VH DCM_template all_Amat x y a_reg_fam_mat DCM

%% run hypothesis based analysis 
% -------------------------------------------------------------------------
[BMA,BMR] = spm_dcm_peb_bmc(PEB_VH, GCM_templates);
save(['../analyses/VN_LR_no_regr_BMA_A_' num2str(length(GCM_templates)) 'models.mat'],...
    'BMA','-v7.3');
save(['../analyses/VN_LR_no_regr_BMR_A_' num2str(length(GCM_templates)) 'models.mat'],...
    'BMR','-v7.3');
%% run family based analysis
% -------------------------------------------------------------------------

% Compare families

%% top-down, bottom-up or both
[BMA_fam_tdbu,fam_tdbu] = spm_dcm_peb_bmc_fam(BMA, BMR, tdbu_family, 'ALL');

%% interhemi, intrahemi or both
[BMA_fam_hemi,fam_hemi] = spm_dcm_peb_bmc_fam(BMA, BMR, hemi_family, 'ALL');

%% regional involvement
[BMA_fam_reg,fam_reg] = spm_dcm_peb_bmc_fam(BMA, BMR, reg_family, 'ALL');

%%
save(['../analyses/VN_LR_no_regr_BMA_' num2str(length(GCM_templates)) 'fam_tdbu.mat'],...
    'BMA_fam_tdbu','fam_tdbu');
save(['../analyses/VN_LR_no_regr_BMA_' num2str(length(GCM_templates)) 'fam_hemi.mat'],...
    'BMA_fam_hemi','fam_hemi');
save(['../analyses/VN_LR_no_regr_BMA_' num2str(length(GCM_templates)) 'fam_reg.mat'],...
    'BMA_fam_reg','fam_reg');

%% report results
% -------------------------------------------------------------------------

% report highest probability models
[max1, idx1] = max(BMA.P(:));
[comms1,diffs1] = ind2sub(size(BMA.P),idx1);
disp(['Highest joint probability across template models: ' num2str(max1)])
disp(['Commonalities: model ' num2str(comms1) ', ' all_A_names{comms1}])
disp(['Differences: model ' num2str(diffs1) ', ' all_A_names{diffs1}])

% report families with best explanatory power

[max2, idx2] = max(fam_tdbu.family.post(:));
[comms2,diffs2] = ind2sub(size(fam_tdbu.family.post),idx2);
disp(['Highest joint probability across tdbu families: ' num2str(max2)])
disp(['Commonalities: family ' num2str(comms2) ', ' a_tdbu_fam_names{comms2}])
disp(['Differences: family ' num2str(diffs2) ', ' a_tdbu_fam_names{diffs2}])

[max3, idx3] = max(fam_hemi.family.post(:));
[comms3,diffs3] = ind2sub(size(fam_hemi.family.post),idx3);
disp(['Highest joint probability across hemi families: ' num2str(max3)])
disp(['Commonalities: family ' num2str(comms3) ', ' a_hemi_fam_names{comms3}])
disp(['Differences: family ' num2str(diffs3) ', ' a_hemi_fam_names{diffs3}])

[max4, idx4] = max(fam_reg.family.post(:));
[comms4,diffs4] = ind2sub(size(fam_reg.family.post),idx4);
disp(['Highest joint probability across reg families: ' num2str(max4)])
disp(['Commonalities: family ' num2str(comms4) ', ' a_reg_fam_names{comms4}])
disp(['Differences: family ' num2str(diffs4) ', ' a_reg_fam_names{diffs4}])

%% thresholded parameter averaging

% Get comms and diffs from template BMA models
Wcomms = full(BMA.Ep(1:(8^2)));
Wdiffs = full(BMA.Ep(((8^2)+1):2*(8^2)));
% Thresholding
Wcomms(BMA.Pw < 0.95) = 0; 
Wdiffs(BMA.Px < 0.95) = 0;
Wcomms = reshape(Wcomms,8,8);
Wdiffs = reshape(Wdiffs,8,8);