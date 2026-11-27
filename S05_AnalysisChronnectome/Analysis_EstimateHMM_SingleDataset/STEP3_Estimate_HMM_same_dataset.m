% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Script to perform temporal filtering, extraction of parcels and
% timecourses, preparation of the input and estimation of the HMM

% Run this script only if you need to prepare the data and estimate the HMM;
% In our case, we have already the matrices of HMM, so we don't need to run it again.

% Dependencies:
% HMM-MAR


%% Define paths and parameters

clearvars;

% Select the directory of the dFC folder
baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
cd(fullfile(projectDir, 'code/'));

dataset = 'dataset2023';
parcels = 1:100; % number of parcels of the atlas used
subjects = importdata(fullfile(projectDir, 'code', ['list_subjects_' dataset '.txt']));

N = length(subjects); % no. subjects
f = cell(N,1);
T = cell(N,1);

% dependency
mydir = '/Volumes/data_cate/GitHub';
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR

derivativesDir = fullfile(projectDir, 'derivatives');
DirData = fullfile(derivativesDir, 'HMM');



%% Estimate the HMM

for j=1:N
    my_string = ['sub-' subjects{j}];
    f{j} = fullfile(DirData, ['sub-' subjects{j}], ['sub-' subjects{j} '_data_Atlas_reduced.txt']);
    load(f{j}); T{j} = 464;
end

K = (2:15); % no. states to be included in the analysis
covtype = 'full'; % type of covariance matrix
zeromean = 0; % if 1 -> the mean of the timeseries will not be used to drive the states
repetitions = 20; % to run it multiple times (keeping all the results)
pca = 10;

% Define the options
options = struct();
options.order = 0; % no autoregressive components
options.covtype = covtype;
options.zeromean = zeromean;
options.verbose = 1;
% standardinzation of the concatenated data (all subjects) for each component
options.standardise = 1; % in order to avoid the segmentation to be dominated by differences in terms of signal scaling or shifting (mean =0, std = 1 for all subjects and components)
options.inittype = 'HMM-MAR';
options.cyc = 500; % maximum number of variational inference cycles
options.initcyc = 10;
options.initrep = 3;

% Options of the single-subject estimation
options_singlesubj = struct();
options_singlesubj.order = 0; % no autoregressive components
options_singlesubj.zeromean = zeromean; % don't model the mean
options_singlesubj.covtype = covtype;
options_singlesubj.standardise = 1;
options_singlesubj.cyc = 100;
options_singlesubj.verbose = 0;


%% Train the model

for pca_i = 1:length(pca)

    DirOut = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K)) '_train_test_' dataset], num2str(pca(pca_i)), ['zeromean' num2str(zeromean)]);
    mkdir(DirOut)

    % Reduce Dimensionality
    % apply the PCA on the data, and then apply the HMM on the clean data
    if numel(pca) > 0
        options.pca = pca(pca_i);
        options_singlesubj.pca = pca(pca_i);
    elseif numel(lowrank) > 0
        % probabilistic PCA observational model (reduce dimensionality and encodes
        % a pattern of functional connectivity)
        options.lowrank= lowrank(pca_i); % can't be run with zeromean = 0
        options_singlesubj.lowrank = lowrank(pca_i);
    end

    for r = 1:repetitions

        fileName = [DirOut '/HMMs_r' num2str(r) '_GROUP.mat'];
        fileInfo = dir(fileName);

        if ~exist(fileName, 'file') || isempty(fileInfo)

            for state_i = 1:length(K)

                k = K(state_i);
                disp(['estimation for state ' num2str(k)])

                options.K = k; % number of states

                % Run the HMM at the group level and get some statistics
                [hmm{r,state_i},Gamma{r,state_i}, Xi{r,state_i}, vpath{r,state_i}] = hmmmar(f,T,options);

                FOgroup{r,state_i}= zeros(N,k); % Fractional occupancy
                meanActivations{r,state_i}= zeros(length(parcels),k); % maps

                for j=1:N
                    ind = (1:T{1}) + T{1}*(j-1);
                    FOgroup{r,state_i}(j,:) = mean(Gamma{r,state_i}(ind,:));
                    X = zscore(dlmread(f{j}));
                    meanActivations{r,state_i}(:,k) = meanActivations{r,state_i}(:,k) + ...
                        sum(X .* repmat(Gamma{r,state_i}(ind,k),1,length(parcels)))';
                end
                % mean activation for all the original components
                meanActivations{r,state_i}= meanActivations{r,state_i} ./ repmat(sum(Gamma{r,state_i}),length(parcels),1);
                % measure of stability per subject (state switching rate for each subject); vector with a value of stability per
                % each subject (rows = subjects; 1 column)
                switchingRate{r,state_i} = getSwitchingRate(Gamma{r,state_i},T,options);
                % temporal stability of the states (vector of life time per state-> number of time points per visit)
                % this matrix has as many columns as the states and the sum of the elements in all columns has to be (subjects * time)
                lifeTime{r,state_i} = getStateLifeTimes(vpath{r,state_i},T{j},options);
                % value for each subject; it is used to assess whether the states are assigned to describe entire subjects (poor job)
                maxFO{r,state_i} = getMaxFractionalOccupancy(Gamma{r,state_i},T,options);
                FO{r,state_i} = getFractionalOccupancy(Gamma{r,state_i},T,options, 2);

                % mean of all states
                if zeromean == 0
                    figure;
                    hmm{r,state_i}.state(k).W.Mu_W
                    hold on;
                    plot(hmm{r,state_i}.state(state_i).W.Mu_W)
                    xlabel('Components', 'FontSize', 18)
                    title('Mean', 'FontSize', 18)
                end

                % get the transition probability matrix
                % TP{r,state_i} = getTransProbs(hmm{r,state_i});
                % visualize the transition probability matrix without the
                % diagonal values
                test =hmm{r,state_i}.P;
                for i = 1:size(hmm{r,state_i}.P,1)
                    test(i,i) = NaN;
                end
                figure; imagesc(test); colorbar

                % compute the free energy and log likelihood (log likelihood of the observed data, per trial -> avLL)
                [fe{r,state_i}, ll] = hmmfe(f, T, hmm{r,state_i});
                ll_average{r,state_i} = mean(ll,1);

                % Subject-level estimation
                options_singlesubj.hmm = hmm;
                FOdual{r,state_i} = zeros(N,k);

                for j = 1:N
                    X = dlmread(f{j});
                    % dual-estimation
                    [HMMs_dual_test{state_i,j},Gamma_test{state_i,j}, vpath_test{state_i,j}] = hmmdual(X,T{j},hmm{r,state_i});
                    for k_i = 1:k
                        HMMs_dual_test{state_i,j}.state(k_i).prior = [];
                    end
                    FO_test{r,state_i}(j,:) = mean(Gamma_test{state_i,j});
                end

                save([DirOut '/HMMs_r' num2str(r) '_SUBJECTS.mat'],...
                    'HMMs_dual_test','Gamma_test','vpath_test', 'FO_test', '-v7.3')
            end

        else
            disp(['The estimation ' fileName ' exists and is not empty'])
        end

        close all;

    end

    save([DirOut '/HMMs_GROUP.mat'],...
        'hmm','FOgroup','meanActivations','switchingRate','maxFO', 'FO', 'fe', 'Gamma', 'vpath', 'lifeTime', 'll_average', '-v7.3')
    save([DirOut '/HMMs_options.mat'], 'options')

end



%% Dual Estimation

subjects_test = importdata(fullfile(projectDir, 'code/Review_October2026/DualEstimation', ['list_subjects_' dataset '.txt']));

N_test = length(subjects_test); % no. subjects
f_test = cell(N_test,1); 
T_test = cell(N_test,1);

for j=1:N_test
    my_string = ['sub-' subjects_test{j}];
    f_test{j} = fullfile(DirData, ['sub-' subjects_test{j}], ['sub-' subjects_test{j} '_data_Atlas_reduced.txt']);           
    load(f_test{j}); T_test{j} = 464;
end 

% Load the selected HMM trained on Dataset 1
load(fullfile(DirOut, '/HMMs_GROUP.mat'));

[~, bestRep] = min(cell2mat(fe));
hmm_train = hmm{bestRep};

% Load Dataset 2 participants and their time series
% Define f_test and T_test for these participants

HMMs_dual_test = cell(N_test, 1);
Gamma_test = cell(N_test, 1);
vpath_test = cell(N_test, 1);
FO_test = zeros(N_test, K);

for j = 1:N_test
    X = dlmread(f_test{j});

    [HMMs_dual_test{j}, Gamma_test{j}, vpath_test{j}] = ...
        hmmdual(X, T_test{j}, hmm_train);

    FO_test(j,:) = mean(Gamma_test{j}, 1);
end

save(fullfile(DirOut, 'HMMs_SUBJECTS'), ...
    'HMMs_dual_test', 'Gamma_test', 'vpath_test', 'FO_test', '-v7.3');