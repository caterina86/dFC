

%% FC comparing controls vs patients

baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
DirData = fullfile(projectDir, 'derivatives/HMM');

mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR
addpath(genpath(fullfile(projectDir, 'code')))


K = (2:15); % no. states to be included in the analysis
k = 6; % IMPORTANT: choose the specific model (number of states)

N = length(subjects); % no. subjects
reduced = 1; % 1 -> reduced number of volumes; 0 -> upsampling of dataset 2017
pca_i = 10;
zeromean = 0;

% number of volumes per subject
if reduced == 1
    T = 464; 
else
    T = 480;
end

if reduced == 1
    if pca_i == 0
        resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], ['zeromean' num2str(zeromean)]);
    else
        resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca_i), ['zeromean' num2str(zeromean)]);
    end        
else
    if pca_i == 0
        resultsDir = fullfile(DirData, ['PCA_UP_' num2str(max(K))], ['zeromean' num2str(zeromean)]);
    else
        resultsDir = fullfile(DirData, ['PCA_UP_' num2str(max(K))], num2str(pca_i),['zeromean' num2str(zeromean)]);
    end        
end

load(fullfile(resultsDir, 'HMMs_GROUP.mat'));

% Find the minimum value and its indices, based on the Free Energy
[minValue_fe, minIndex_fe] = min(cell2mat(fe(:,:)));
[rowIndex_fe, colIndex_fe] = ind2sub(size(fe), minIndex_fe);
rep_i = rowIndex_fe(1,k); % repetition chosen based on the lowest free energy

disp(['repetition chosen ' num2str(rep_i)])
disp(['number states chosen ' num2str(K(k))])

load(fullfile(resultsDir, ['HMMs_r' num2str(rep_i) '_SUBJECTS.mat']));
mkdir(fullfile(resultsDir, 'figures/'));
    
state_i = K(k); 
outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)]);



%% NBS parameters

% options -> t-test; permutations = 10000; significance = 0.05; method = NBS; 
% Components size = Extent

% Create the contrast file
contrast = [-1 1 0 0 0];
contrastFile1 = fullfile(outputDir, 'NBS/contrast_EBvsSC.txt');
writematrix(contrast, contrastFile1, 'Delimiter','tab');  % tab-delimited

contrast = [1 -1 0 0 0];
contrastFile2 = fullfile(outputDir, 'NBS/contrast_SCvsEB.txt');
writematrix(contrast, contrastFile2, 'Delimiter','tab');  % tab-delimited

thresholds = 3.5;
sizes = {'Extent'}; %{'Intensity'}; %{'Extent'};

% Number of subjects
nSubjects = 83;
% Create S struct
S = struct();
S.Nsub = nSubjects;
S.OUT.waitbar = [];   % optional

num_perm = '10000';

outputStatDir = fullfile(outputDir, 'NBS', sizes{1}, 'staticFC');
mkdir(outputStatDir);



%% Assign the mean FC of controls to patients

shuffle_state = 3;

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFC_state' num2str(shuffle_state) '_Mean']);
mkdir(outputShuffle);

% =====================================================
% PRECOMPUTE CONTENT THAT DOES NOT DEPEND ON PERMUTATION
% =====================================================
load(fullfile(outputStatDir, 'matrices_staticFC.mat'))

input = static_FCz; % your real static FC (Fisher-z)
[m, n, Sub] = size(input);
triangle = triu(true(m,n),1);
E = sum(triangle(:));

patients = 1:33;
controls = 34:83;
nP = numel(patients);
nC = numel(controls);

% ---------------------------------------------------------
% Compute FC for all states
% ---------------------------------------------------------
FC_allStates = cell(Sub, state_i);
for st = 1:state_i
    for sub_i = 1:Sub
        [~, FC_allStates{sub_i, st}] = getFuncConn(HMMs_dualregr{k, sub_i}, st);
    end
end

% ---------------------------------------------------------
% 1) ORIGINAL FO COLUMN (unchanged)
% ---------------------------------------------------------
FO_original = FOdual{rep_i,k}(:, shuffle_state);

% ---------------------------------------------------------
% 2) GROUP-SWAP OF FC ONLY (take mean FC of the opposite group)
% ---------------------------------------------------------
FC_perm = FC_allStates;  % start from original FC

if shuffle_state == 1 || shuffle_state == 6
    % Patients receive mean FC of controls
    mean_FC_controls(:, :) = mean(cat(3, FC_allStates{controls, shuffle_state}), 3);
    for sub_i = patients
        FC_perm{sub_i, shuffle_state} = mean_FC_controls(:, :);
    end
    % Controls remain unchanged
else
    mean_FC_patients(:, :) = mean(cat(3, FC_allStates{patients, shuffle_state}), 3);
    for sub_i = controls
        FC_perm{sub_i, shuffle_state} = mean_FC_patients(:, :);
    end
    % Patients remain unchanged
end

% ---------------------------------------------------------
% 3) Reconstruct static FC with original FO
% ---------------------------------------------------------
full_FC_perm = zeros(m, m, Sub);


for sub_i = 1:Sub
    FC_static_sub = zeros(m, m);
    FO_use = zeros(1,state_i);
    
    for st = 1:state_i
        FO_use(st) = FO_original(sub_i);  % keep original FO
        FC_static_sub = FC_static_sub + FO_use(st) * FC_perm{sub_i, st};
    end
    
    % normalize / bound / zero diagonal
    FC_static_sub = min(max(FC_static_sub, -0.999999), 0.999999);
    FC_static_sub(eye(m) == 1) = 0;
    
    % Fisher Z
    full_FC_perm(:,:,sub_i) = atanh(FC_static_sub);
end

% ---------------------------------------------------------
% 4) SAVE MATRIX FOR NBS
% ---------------------------------------------------------
tmpfile = fullfile(outputShuffle, 'tmp_FCperm_mean.mat');
full_FO_perm(isnan(full_FC_perm)) = 0;
save(tmpfile, 'full_FC_perm');


% ========= RUN NBS ==========
if shuffle_state == 1 || shuffle_state == 6
    contrastFile = contrastFile1;
    label = 'EBvsSC';
elseif shuffle_state == 3 || shuffle_state == 7
    contrastFile = contrastFile2;
    label = 'SCvsEB';
end

clear global nbs
UI = [];
UI.test.ui = 't-test';
UI.size.ui = sizes{1};
UI.matrices.ui = tmpfile;
UI.design.ui = fullfile(outputDir, 'NBS/designMatrix.mat');
UI.contrast.ui = contrastFile;
UI.exchange.ui = '';
UI.node_label.ui = fullfile(outputDir, 'NBS/nodeLabels.txt');
UI.node_coor.ui = fullfile(outputDir, 'NBS/COG.txt');
UI.perms.ui = num_perm;
UI.method.ui = 'Run NBS';
UI.thresh.ui = num2str(thresholds(1));
UI.alpha.ui = '0.05';

NBSrun(UI, S);
close all;

global nbs

if ~isempty(nbs.NBS.con_mat)
    tmat = nbs.NBS.test_stat;

    merged_mask = false(size(tmat));
    for ci = 1:numel(nbs.NBS.con_mat)
        merged_mask = merged_mask | full(nbs.NBS.con_mat{ci});
    end

    t_masked = tmat .* merged_mask;
    t_masked(1:size(t_masked,1)+1:end) = 0;
    t_masked(isnan(t_masked)) = 0;

    tvals_perm(:, 1) = t_masked(triangle);
else
    tvals_perm(:, 1) = 0;
end

save(fullfile(outputShuffle, 'tvals.mat'), 'tvals_perm')


% load the output from NBS
if shuffle_state == 1 | shuffle_state == 6
    load(fullfile(outputStatDir, 'Results_staticFC_EBvsSC.mat'))
else
    load(fullfile(outputStatDir, 'Results_staticFC_SCvsEB.mat'))
end

tvalues = nbs_outputs.NBS.test_stat;    % t-stat matrix        
tvalues(tvalues < thresholds) = 0; % remove all values of t < thr
numComponents = numel(nbs_outputs.NBS.con_mat);
merged_mask = false(size(tvalues));      % initialize empty logical mask
    
for i = 1:numComponents
    comp = full(nbs_outputs.NBS.con_mat{i});    % convert to full
    merged_mask = merged_mask | (comp == 1);    % logical OR to merge
end
t_masked_static = tvalues .* merged_mask;           % masked matrix
t_masked_static(1:size(t_masked_static,1)+1:end) = 0;      % zero diagonal
t_masked_static(isnan(t_masked_static)) = 0;               % replace NaNs
tvals_masked = t_masked_static(triangle);

% Choose the tvals within the connections that are significant in the empirical static FC
tvals_perm_nan = tvals_perm;      % copy original
sig_rows = tvals_masked ~= 0;     % rows to keep
n_nonzero_staticFC = nnz(tvals_masked);

% Set non-significant rows to NaN
tvals_perm_nan(~sig_rows, :) = NaN;
tvals_perm_nan(tvals_perm_nan == 0) = NaN;
tvals_perm_nan(tvals_perm_nan < thresholds) = NaN; % remove all values of t < thr

mean(tvals_perm_nan, 'omitnan');

n_nonNaN_shuffling = sum(~isnan(tvals_perm_nan), 1);
perc_edges = (max(n_nonNaN_shuffling)/n_nonzero_staticFC)*100;
disp(['Percentage significant connections = ' num2str(perc_edges)])



