% Caterina Pedersini, October 2026, cpedersini86@gmail.com

% Plot single-subjects trend for each model order

% Dependencies:
% HMM-MAR

clearvars; close all;

% load paths
baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
DirData = fullfile(projectDir, 'derivatives/HMM');

mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR
addpath(genpath(fullfile(projectDir, 'code')))

zeromean = 0;
K = 2:15; % no. states to be included in the analysis
pca = 10;
repetitions = 20;
T = 464; % number of timepoints

%resultsDir = fullfile(DirData, ['REDUCED_NOpca_' num2str(max(K))], '0', ['zeromean' num2str(zeromean)]);
resultsDir = fullfile(DirData, 'PCA_REDUCED_15', num2str(pca), ['zeromean' num2str(zeromean)]);
load(fullfile(resultsDir, 'HMMs_GROUP.mat'));
load(fullfile(resultsDir, 'HMMs_options.mat'));
N = length(subjects);
  



%% Choose the repetition with the lowest free energy, for each model order

for num_states_i = 1:length(K)
    
    state_i = K(num_states_i);
    k = state_i-1; 

    % Find the minimum value and its indices, based on the Free Energy
    [minValue_fe, minIndex_fe] = min(cell2mat(fe(:,:)));
    [rowIndex_fe, colIndex_fe] = ind2sub(size(fe), minIndex_fe);

    % Find the minimum value and its indices, based on the Free Energy
    [minValue_ll, minIndex_ll] = min(cell2mat(ll_average(:,:)));
    [rowIndex_ll, colIndex_ll] = ind2sub(size(ll_average), minIndex_ll);

    rep_i = rowIndex_fe(1,k); % repetition chosen based on the lowest free energy
    disp(['repetition chosen ' num2str(rep_i)])
    disp(['number states chosen ' num2str(K(k))])

    % Plot single-subjects dynamics
    mkdir(fullfile(resultsDir, 'figures/'));
    outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)]);
    mkdir(outputDir);

    % Viterbi Path
    vector = vpath{rep_i,k};
    numSegments = numel(vector) / T;
    
    for i = 1:numSegments
    
        startIndex = (i-1)*T + 1;
        endIndex = i*T;
    
        v = vector(startIndex:endIndex);
    
        figure('Color','w','Position',[100 100 900 300]);
    
        stairs(1:T, v, 'LineWidth', 1.5);
    
        ylim([0 hmm{rep_i,k}.K+1]);
        xlim([1 T]);
    
        yticks(1:hmm{rep_i,k}.K);
        xlabel('Time points');
        ylabel('Hidden state');
        title(['Viterbi path – subject ' subjects{i}]);
    
        set(gca,'FontSize',13,'LineWidth',1.2, ...
            'TickDir','out','Box','off');
    
        exportgraphics(gcf, ...
            fullfile(outputDir, ...
            ['sub-' subjects{i} '_ViterbiPath.png']), ...
            'Resolution',300);
    
        close(gcf);
    end
    
    % Gamma
    for sub_i = 1:N
    
        % Extract subject-specific Gamma
        startIndex = (sub_i-1)*T + 1;
        endIndex = sub_i*T;
    
        G = Gamma{rep_i,k}(startIndex:endIndex,:);
        
        % Plot Gamma
        figure('Color','w');
        plot_Gamma(G,T,0);
        sgtitle(['Subject ' subjects{sub_i} ' – Group-level Gamma']);
    
        saveas(gcf,fullfile(outputDir, ...
            ['sub-' subjects{sub_i} '_GroupGamma.png']));
    
        close(gcf);
    
   
    end


end

