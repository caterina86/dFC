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
            FC_static_sub = FC_static_sub + FC_allStates{sub_i, st};
        end
        
        FC_static_sub = FC_static_sub / (state_i - 1);

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
    tmpfile = fullfile(outputStatDir, 'Reconstruction_noFO', sprintf('FC_excluding_state_%02d.mat', excl_st));
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

        save(fullfile(outputStatDir, 'Reconstruction_noFO', ['Results_NOState' num2str(excl_st) '_' label '.mat']), 'nbs');      
        clear nbs;  

    end
end




