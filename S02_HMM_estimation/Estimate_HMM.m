% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Script to perform temporal filtering, extraction of parcels and
% timecourses, preparation of the input and estimation of the HMM

% Run this script only if you need to prepare the data and estimate the HMM;
% In our case, we have already the matrices of HMM, so we don't need to run it again.

% Dependencies:
% HMM-MAR


%% Define paths and parameters

% Select the directory of the dFC folder
baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
cd(fullfile(projectDir, 'code/'));

parcels = 1:100; % number of parcels of the atlas used
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
filtering = 1; % 1 -> bandpass filtering; 2 -> highpass filtering (f=1/2000HZ); 3 -> highpass filtering (f=0.01HZ)

N = length(subjects); % no. subjects
f = cell(N,1); 
T = cell(N,1);

% dependency
mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR



%% Preprocessing steps after fmriprep and ICA denoising

% To run this code, you need to have the:
% 1. ica denoised data
% 2. the atlas in the derivatives folder

% extract parcels from the atlas
system(['bash S01_extract_parcels.sh ' baseDir]);

% Temporal filtering -> bandpass filtering 0.01-0.1HZ
filter='0.01_0.1';
system(['bash S02_filtering_bandpass.sh ' baseDir]);
derivativesDir = fullfile(projectDir, 'derivatives');
atlasDir = (fullfile(derivativesDir, 'Atlas'));
mkdir(atlasDir)
% extract timecourses from the filtered and clean functional data
system(['bash S03_extract_timecourses.sh ' baseDir ' ' atlasDir ' ' projectDir ' ' filter]);

DirData = fullfile(derivativesDir, 'HMM');



%% Prepare the data for the HMM

% extract the timecourses based on the atlas parcellation

disp('Extract Input based on Atlas')
for sub_i = 1:length(subjects) % for all subjects
    
    for ii = 1:length(parcels) % for all parcels           
        atlasDir = fullfile(derivativesDir, 'Atlas/', ['sub-' subjects{sub_i}], 'timecourses/');
        textFilename = sprintf(fullfile(atlasDir, 'timecourse_parcel%d.txt'), ii);
        data(:,ii) = load(textFilename, '-ascii');
    end
    
    HMMDir = fullfile(derivativesDir, 'HMM/', ['sub-' subjects{sub_i}]);
    mkdir(HMMDir)
    
    writematrix(data, fullfile(HMMDir, ['sub-' subjects{sub_i} '_data_Atlas.txt']))   
    clear data;
    
end



%% Estimate the HMM

% remove the first volumes
for j=1:N
    my_string = ['sub-' subjects{j}];

    if contains(my_string, 'blnd') || contains(my_string, 'ctrl')
       values = load(fullfile(DirData, ['sub-' subjects{j}], ['sub-' subjects{j} '_data_Atlas.txt']));
       values(1:8,:) = [];

       % Save upsampled_values_spline to a text file
       save(fullfile(fullfile(DirData, ['sub-' subjects{j}], ['sub-' subjects{j} '_data_Atlas_reduced.txt'])), 'values', '-ascii');

    else
       values = load(fullfile(DirData, ['sub-' subjects{j}], ['sub-' subjects{j} '_data_Atlas.txt']));
       values(1:16,:) = [];

       % Save upsampled_values_spline to a text file
       save(fullfile(fullfile(DirData, ['sub-' subjects{j}], ['sub-' subjects{j} '_data_Atlas_reduced.txt'])), 'values', '-ascii');

    end
end

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

for pca_i = 1:length(pca)
        
    DirOut = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca(pca_i)), ['zeromean' num2str(zeromean)]);
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
                
                % cross-validated log-likelihood
                %[mcv{r,state_i}, cv{r,state_i}] = cvhmmmar(f,T,options);

                % Subject-level estimation
                options_singlesubj.hmm = hmm;
                FOdual{r,state_i} = zeros(N,k);

                for j = 1:N
                    X = dlmread(f{j});
                    % dual-estimation
                    [HMMs_dualregr{state_i,j},Gammaj{state_i,j}, vpathj{state_i,j}] = hmmdual(X,T{j},hmm{r,state_i}); 
                    for k_i = 1:k
                        HMMs_dualregr{state_i,j}.state(k_i).prior = [];
                    end
                    FOdual{r,state_i}(j,:) = mean(Gammaj{state_i,j});
                end
  
            end
            
                save([DirOut '/HMMs_r' num2str(r) '_SUBJECTS.mat'],...
            'HMMs_dualregr','Gammaj','FOdual', '-v7.3')

               
        else
            disp(['The estimation ' fileName ' exists and is not empty'])
        end
            
        close all;
  
    end
    
    save([DirOut '/HMMs_GROUP.mat'],...
        'hmm','FOgroup','meanActivations','switchingRate','maxFO', 'FO', 'fe', 'Gamma', 'vpath', 'lifeTime', 'll_average', '-v7.3')
    save([DirOut '/HMMs_options.mat'], 'options')
                
end
    
