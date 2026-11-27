
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



%% Shuffle FO values of controls to patients, and estimate the model based FC

% =====================================================
% FO SHUFFLING – WITHOUT RESIDUALS
% =====================================================

shuffle_state = 6;

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFO_state' num2str(shuffle_state)]);
mkdir(outputShuffle);

perm_N = 100;

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


% =====================================================
% PRECOMPUTE STATE-LEVEL FC FOR ALL STATES
% =====================================================

FC_allStates = cell(Sub, state_i);
for st = 1:state_i
    for sub_i = 1:Sub
        [~, FC_allStates{sub_i, st}] = getFuncConn( ...
            HMMs_dualregr{k, sub_i}, st);
    end
end

% Pattern (any number after tmp_FOperm_)
pattern = fullfile(outputShuffle, 'tmp_FOperm_*.mat');
files = dir(pattern);

% Extract numbers from filenames
nums = zeros(numel(files),1);
for i = 1:numel(files)
    % Extract any number in the filename (before extension)
    token = regexp(files(i).name, '\d+', 'match');
    if ~isempty(token)
        nums(i) = str2double(token{end});   % use last number if multiple
    else
        nums(i) = NaN;                      % no number found
    end
end

% Sort by the extracted numbers
[~, idx] = sort(nums);
files = files(idx);

if ~isempty(files)
    % Extract numbers from filenames
    nums = zeros(length(files),1);
    for ii = 1:length(files)
        fname = files(ii).name;
        % Remove prefix and suffix
        numStr = erase(fname, {'tmp_FOperm_', '.mat'});
        nums(k) = str2double(numStr);
    end
    
    % Find the highest completed permutation
    lastPerm = max(nums);
    startPerm = lastPerm + 1;
    
    % Load existing t-values and mean values
    tvalsFile = fullfile(outputShuffle, 'tvals.mat');
    meanvalsFile = fullfile(outputShuffle, 'mean_vals.mat');

    if exist(tvalsFile, 'file')
        load(tvalsFile, 'tvals_perm');     % assumes variable is tvalues
    end
    
    if exist(meanvalsFile, 'file')
        load(meanvalsFile, 'value_mean');   % assumes variable is meanFC
    end
    
else
    startPerm = 1; % no previous files
end

% Determine how many previous permutations exist
if exist('value_mean','var') && ~isempty(value_mean)
else
    value_mean = [];   % initialize empty cell array
end

if exist('tvals_perm','var') && ~isempty(tvals_perm)
else
    tvals_perm = [];
end

fprintf('Starting permutation loop from %d\n', startPerm);


% =====================================================
% PERMUTATION LOOP
% =====================================================

for perm_i = startPerm:perm_N

    fprintf('Permutation %d / %d\n', perm_i, perm_N);

    % ---------------------------------------------------------
    % 1) ORIGINAL FO COLUMN
    % ---------------------------------------------------------
    FO_original = FOdual{rep_i,k}(:, shuffle_state);

    % ---------------------------------------------------------
    % 2) GROUP-SWAP OF FO ONLY FOR shuffle_state
    % ---------------------------------------------------------
    FO_perm = FO_original;  % start from original
    
    if shuffle_state == 1 | shuffle_state == 6
        % Patients receive FO from random controls
        rand_controls = controls(randi(nC, nP, 1));
        FO_perm(patients) = FO_original(rand_controls);
        % Controls remain unchanged
        FO_perm(controls) = FO_original(controls);  
    else
        % Controls receive FO from random patients
        rand_patients = patients(randi(nP, nC, 1));
        FO_perm(controls) = FO_original(rand_patients);        
        % Controls remain unchanged
        FO_perm(patients) = FO_original(patients);
    end

    value_mean(perm_i, 1) = mean(FO_perm(patients));
    value_mean(perm_i, 2) = mean(FO_perm(controls));    

    full_FO_perm = zeros(m, m, Sub);

    for sub_i = 1:Sub

        % start with zero
        FC_static_sub = zeros(m, m);
        FO_use = zeros(1,7);

        for st = 1:state_i

            if st == shuffle_state
                FO_use(1,st) = FO_perm(sub_i); % SHUFFLED FO
            else
                FO_use(1,st) = FOdual{rep_i,k}(sub_i, st); % original FO
            end

            % accumulate weighted FC
            FC_static_sub = FC_static_sub + ...
                FO_use(1,st) * FC_allStates{sub_i, st};

        end

        % Normalize by the total proportion of time spent at each state
        %total_proportion_time = sum(FO_use(1,:));
        %FC_static_sub = FC_static_sub / total_proportion_time;
        FC_static_sub = min(max(FC_static_sub, -0.999999), 0.999999);
        FC_static_sub(logical(eye(100))) = 0;

        % Fisher Z
        full_FO_perm(:,:,sub_i) = atanh(FC_static_sub);
    end

    % ---------------------------------------------------------
    % 4) SAVE MATRIX FOR NBS
    % ---------------------------------------------------------
    tmpfile = fullfile(outputShuffle, sprintf('tmp_FOperm_%d.mat', perm_i));
    full_FO_perm(isnan(full_FO_perm)) = 0;
    save(tmpfile, 'full_FO_perm');

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

    % ---------------------------------------------------------
    % 6) EXTRACT SIGNIFICANCE AND T-VALUES
    % ---------------------------------------------------------
    if ~isempty(nbs.NBS.con_mat)
        tmat = nbs.NBS.test_stat;

        merged_mask = false(size(tmat));
        for ci = 1:numel(nbs.NBS.con_mat)
            merged_mask = merged_mask | full(nbs.NBS.con_mat{ci});
        end

        t_masked = tmat .* merged_mask;
        t_masked(1:size(t_masked,1)+1:end) = 0;
        t_masked(isnan(t_masked)) = 0;

        tvals_perm(:, perm_i) = t_masked(triangle);
    else
        tvals_perm(:, perm_i) = 0;
    end

    save(fullfile(outputShuffle, 'tvals.mat'), 'tvals_perm')
    save(fullfile(outputShuffle, 'mean_vals.mat'), 'value_mean')

end




%% Shuffle the Functional Connectome of controls to patients, and estimate the model based FC

% =====================================================
% FUNCTIONAL CONNECTOME SHUFFLING
% =====================================================

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFC_state' num2str(shuffle_state)]);
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


% =====================================================
% PRECOMPUTE STATE-LEVEL FC FOR ALL STATES
% =====================================================

FC_allStates = cell(Sub, state_i);
for st = 1:state_i
    for sub_i = 1:Sub
        [~, FC_allStates{sub_i, st}] = getFuncConn( ...
            HMMs_dualregr{k, sub_i}, st);
    end
end

% Pattern (any number after tmp_FOperm_)
pattern = fullfile(outputShuffle, 'tmp_FCperm_*.mat');
files = dir(pattern);

% Extract numbers from filenames
nums = zeros(numel(files),1);
for i = 1:numel(files)
    % Extract any number in the filename (before extension)
    token = regexp(files(i).name, '\d+', 'match');
    if ~isempty(token)
        nums(i) = str2double(token{end});   % use last number if multiple
    else
        nums(i) = NaN;                      % no number found
    end
end

% Sort by the extracted numbers
[~, idx] = sort(nums);
files = files(idx);


if ~isempty(files)
    % Extract numbers from filenames
    nums = zeros(length(files),1);
    for ii = 1:length(files)
        fname = files(ii).name;
        % Remove prefix and suffix
        numStr = erase(fname, {'tmp_FCperm_', '.mat'});
        nums(k) = str2double(numStr);
    end
    
    % Find the highest completed permutation
    lastPerm = max(nums);
    startPerm = lastPerm + 1;
    
    % Load existing t-values and mean values
    tvalsFile = fullfile(outputShuffle, 'tvals.mat');
    meanvalsFile = fullfile(outputShuffle, 'mean_vals.mat');

    if exist(tvalsFile, 'file')
        load(tvalsFile, 'tvals_perm');     % assumes variable is tvalues
    end
    
    if exist(meanvalsFile, 'file')
        load(meanvalsFile, 'meanFC');   % assumes variable is meanFC
    end
    
else
    startPerm = 1; % no previous files
end

% Determine how many previous permutations exist
if exist('meanFC','var') && ~isempty(meanFC)
    prevPerms = size(meanFC,1);
else
    prevPerms = 0;
    meanFC = {};   % initialize empty cell array
end

if exist('tvals_perm','var') && ~isempty(tvals_perm)
    prevPerms_tvals = size(tvals_perm,2);
else
    prevPerms_tvals = 0;
    tvals_perm = [];
end

fprintf('Starting permutation loop from %d\n', startPerm);


% =====================================================
% PERMUTATION LOOP
% =====================================================
for perm_i = startPerm:perm_N

    fprintf('Permutation %d / %d\n', perm_i, perm_N);

    % ---------------------------------------------------------
    % 1) GROUP-SWAP OF FC ONLY FOR shuffle_state
    % ---------------------------------------------------------
    
    FC_perm = FC_allStates;
    
    % Patients receive FC from random controls
    rand_controls = controls(randi(nC, nP, 1));
    rand_patients = patients(randi(nP, nC, 1));

    if shuffle_state == 1 | shuffle_state == 6
        for ii = 1:nP
            FC_perm{patients(ii), shuffle_state} = ...
                FC_allStates{rand_controls(ii), shuffle_state};
        end 
    else
        for ii = 1:nC
            FC_perm{controls(ii), shuffle_state} = ...
                FC_allStates{rand_patients(ii), shuffle_state};
        end 
    end

    M = size(FC_perm{1,1},1); 

    % Initialize new entry in meanFC
    meanFC{perm_i,1} = zeros(M);   % patients
    meanFC{perm_i,2} = zeros(M);   % controls

    % patients
    for s = patients
        meanFC{perm_i,1} = meanFC{perm_i,1} + FC_perm{s,shuffle_state};
    end
    meanFC{perm_i,1} = meanFC{perm_i,1} / nP;
    
    % controls
    for s = controls
        meanFC{perm_i,2} = meanFC{perm_i,2} + FC_perm{s,shuffle_state};
    end
    meanFC{perm_i,2} = meanFC{perm_i,2} / nC;    

    m = 100;

    full_FC_perm = zeros(m, m, Sub);
    
    for sub_i = 1:Sub
    
        FC_static_sub = zeros(m, m);
        FOw = zeros(1,state_i);
    
        for st = 1:state_i
            % Use ORIGINAL FOdual
            FOw(1,st) = FOdual{rep_i,k}(sub_i, st);
    
            % Use SHUFFLED FC_perm
            FC_use = FC_perm{sub_i, st};
    
            % accumulate weighted FC
            FC_static_sub = FC_static_sub + ...
                FOw(1,st) * FC_use;

        end
    
        % Normalize by the total proportion of time spent at each state
        %total_proportion_time = sum(FOw(1,:));
        %FC_static_sub = FC_static_sub / total_proportion_time;

        FC_static_sub(logical(eye(100))) = 0;
        FC_static_sub = min(max(FC_static_sub, -0.999999), 0.999999);

        % Fisher Z
        full_FC_perm(:,:,sub_i) = atanh(FC_static_sub);

    end


    % ---------------------------------------------------------
    % 4) SAVE MATRIX FOR NBS
    % ---------------------------------------------------------
    tmpfile = fullfile(outputShuffle, sprintf('tmp_FCperm_%d.mat', perm_i));
    full_FC_perm(isnan(full_FC_perm)) = 0;
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

    % ---------------------------------------------------------
    % 6) EXTRACT SIGNIFICANCE AND T-VALUES
    % ---------------------------------------------------------

    tvals_perm(:, perm_i) = zeros(4950,1);

    if ~isempty(nbs.NBS.con_mat)
        tmat = nbs.NBS.test_stat;
        merged_mask = false(size(tmat));
        for ci = 1:numel(nbs.NBS.con_mat)
            merged_mask = merged_mask | full(nbs.NBS.con_mat{ci});
        end

        t_masked = tmat .* merged_mask;
        t_masked(1:size(t_masked,1)+1:end) = 0;
        t_masked(isnan(t_masked)) = 0;

        tvals_perm(:, perm_i) = t_masked(triangle);
    else
        tvals_perm(:, perm_i) = 0;
    end

    save(fullfile(outputShuffle, 'tvals.mat'), 'tvals_perm')
    save(fullfile(outputShuffle, 'mean_vals.mat'), 'meanFC')

end



