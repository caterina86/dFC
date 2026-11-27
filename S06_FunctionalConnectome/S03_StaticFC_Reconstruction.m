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

m = 100;
Sub = 83;

mkdir(fullfile(outputStatDir, 'Reconstruction_noFO'));

FC_allStates = cell(Sub, state_i);
for st = 1:state_i
    for sub_i = 1:Sub
        [~, FC_allStates{sub_i, st}] = getFuncConn( ...
            HMMs_dualregr{k, sub_i}, st);
    end
end




%% Exclude one state at a time
for excl_st = 1:state_i % loop over each state to exclude
    
    for sub_i = 1:Sub
        
        FC_static_sub = zeros(m, m);

        for st = 1:state_i

            if st == excl_st
                % Skip this state in this permutation
                continue;
            end

            FOw = FOdual{rep_i,k}(sub_i, st);
            FC_static_sub = FC_static_sub + FOw * FC_allStates{sub_i, st};
        end
        
        % --- post-processing ---
        FC_static_sub(logical(eye(m))) = 0;
        FC_static_sub = min(max(FC_static_sub, -0.999999), 0.999999);

        % Fisher Z
        full_FC_perm_excl(:,:,sub_i, excl_st) = atanh(FC_static_sub);

    end
end


for excl_st = 1:state_i
    
    % Create NBS-format cell array
    matrices = cell(Sub,1);
    for s = 1:Sub
        matrices{s} = full_FC_perm_excl(:,:,s,excl_st);
    end

    % Save Temporary file
    tmpfile = fullfile(outputStatDir, 'Reconstruction', sprintf('FC_excluding_state_%02d.mat', excl_st));
    save(tmpfile, 'matrices')

    for cont_i = 1:2
    
        if cont_i == 1
            contrastFile = contrastFile1;
            label = 'EBvsSC';
        else
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

        save(fullfile(outputStatDir, 'Reconstruction', ['Results_NOState' num2str(excl_st) '_' label '.mat']), 'nbs');      
        clear nbs;  

    end
end




%% ===== SANITY CHECK: per-state group contribution =====

sanity = struct();
sanity.delta_mean = zeros(state_i,1);
sanity.delta_pos_frac = zeros(state_i,1);
sanity.delta_neg_frac = zeros(state_i,1);

patients_idx = 1:33;
controls_idx = 34:83;

for st = 1:state_i
    
    % accumulate subject-wise weighted FC
    FCw_pat = [];
    FCw_con = [];
    
    for sub_i = 1:Sub
        
        FCw = FOdual{rep_i,k}(sub_i,st) * FC_allStates{sub_i,st};
        
        if ismember(sub_i, patients_idx)
            FCw_pat(:,:,end+1) = FCw;
        elseif ismember(sub_i, controls_idx)
            FCw_con(:,:,end+1) = FCw;
        end
    end
    
    % group means
    FCw_pat_mean = median(FCw_pat, 3, 'omitnan');
    FCw_con_mean = median(FCw_con, 3, 'omitnan');
    
    % difference
    Delta = FCw_pat_mean - FCw_con_mean;
    
    % ignore diagonal
    Delta(logical(eye(m))) = NaN;
    
    % summary stats
    sanity.delta_mean(st)     = mean(Delta(:), 'omitnan');
    sanity.delta_pos_frac(st) = mean(Delta(:) > 0, 'omitnan');
    sanity.delta_neg_frac(st) = mean(Delta(:) < 0, 'omitnan');
    
end

save(fullfile(outputStatDir, 'Reconstruction', 'SanityCheck_StateContributions.mat'), 'sanity');




%% Exclude a pair of states

% Define states to exclude
states_to_exclude = [1, 6];

% Preallocate
full_FC_perm_excl_pair = zeros(m, m, Sub);

for sub_i = 1:Sub
    FC_static_sub = zeros(m, m);

    for st = 1:state_i
        if ismember(st, states_to_exclude)
            continue;  % skip excluded states
        end

        FOw = FOdual{rep_i,k}(sub_i, st);
        FC_static_sub = FC_static_sub + FOw * FC_allStates{sub_i, st};
    end

    % --- post-processing ---
    FC_static_sub(logical(eye(m))) = 0;
    FC_static_sub = min(max(FC_static_sub, -0.999999), 0.999999);

    % Fisher Z
    full_FC_perm_excl_pair(:,:,sub_i) = atanh(FC_static_sub);
end

    
% Create NBS-format cell array
matrices = cell(Sub,1);
for s = 1:Sub
    matrices{s} = full_FC_perm_excl_pair(:,:,s);
end

% Save Temporary file
tmpfile = fullfile(outputStatDir, 'Reconstruction', ['FC_excluding_state' num2str(states_to_exclude(1,1)) '_state' num2str(states_to_exclude(1,2)) '.mat']);
save(tmpfile, 'matrices')

for cont_i = 1:2

    if cont_i == 1
        contrastFile = contrastFile1;
        label = 'EBvsSC';
    else
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

    save(fullfile(outputStatDir, 'Reconstruction', ['Results_NOState' num2str(states_to_exclude(1,1)) '_State' num2str(states_to_exclude(1,2)) '_' label '.mat']), 'nbs');      
    clear nbs;  

end




%% Create figure for the reconstruction (after create the .edge file)

label = 'EBvsSC'; %'SCvsEB'; %'EBvsSC';
label_matrix = 'TimeAveragedFC'; %'TimeAveragedFC'; % 'staticFC'

% load the model-based static FC
staticFC = load(fullfile(outputStatDir, ['Results_' label_matrix '_' label '.edge']));
staticFC (staticFC  < 4) = 0;
staticFC = double(staticFC ~= 0);

for k_i = 7
    % load the specific recontruction
    ReconFC = load(fullfile(outputStatDir, 'Reconstruction', ['Results_NOState' num2str(k_i) '_' label '.edge']));
    % Zero out values < 4
    ReconFC(ReconFC < 4) = 0;
    ReconFC = double(ReconFC ~= 0);
    
    % ------------------------------------------------
    % Encode edge presence
    % 1 = static only
    % 2 = both
    % 3 = recon only
    % ------------------------------------------------
    c = zeros(size(staticFC));

    c(staticFC == 1 & ReconFC == 0) = 1;
    c(staticFC == 1 & ReconFC == 1) = 2;
    c(staticFC == 0 & ReconFC == 1) = 3;

    % Flatten and remove zeros (no connection)
    c_vec = c(:);
    c_vec(c_vec == 0) = [];
    
    N = numel(c_vec);  % total number of connections
    
    perc_static_only = 100 * sum(c_vec == 1) / N;
    perc_both  = 100 * sum(c_vec == 2) / N;
    perc_recon_only = 100 * sum(c_vec == 3) / N;

    disp(['State ' num2str(k_i) ' Perc Static Only ' num2str(perc_static_only) ' Perc Recon only ' num2str(perc_recon_only) ' Both ' num2str(perc_both)])

    dlmwrite(fullfile(outputStatDir, 'Reconstruction', [label_matrix '_NOState' num2str(k_i) '_' label '.edge']), c, 'delimiter', ' ');

end

