% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Compare models with different number of states.
% For each model order, we assessed the distance between FC matrices when 
% each state is active and calculate the mean distance between FC of
% different states, within the same model order. We do the same for the
% mean activation maps, to test which model order is composed of states
% that are more dissimilar among them.

% To run this script you need to run script STEP1_SingleSubejcts_Path.m for
% each state between the range specified as K (in this case, 2:15)

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
K = (2:15); % no. states to be included in the analysis
pca = 10;
repetitions = 20;
T = 464; % number of timepoints
num_networks = 7;

resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca), ['zeromean' num2str(zeromean)]);
mkdir(fullfile(resultsDir, 'figures/Distances'));



%% Geodesic distance and pearson dissimilarity between FC of models with different number of states

for ii = 1:length(K)
    
    states = K(ii);
    disp(states)
    
    % Create the path to the NumStates_XX directory
    numStatesDir = fullfile(resultsDir, ['figures/NumStates_' num2str(states)]);

    % List all directories in the NumStates_XX directory
    repetitionsDirs = dir(numStatesDir);
    repetitionsDirs = repetitionsDirs([repetitionsDirs.isdir]); % Filter out non-directory entries
    repetitionsDirs = repetitionsDirs(~ismember({repetitionsDirs.name}, {'.', '..'})); % Exclude '.' and '..'

    % Get the first directory found (assuming there's only one Repetitions_XX directory)
    if ~isempty(repetitionsDirs)
        firstRepetitionDir = repetitionsDirs(1).name;
        matrixDir = fullfile(numStatesDir, firstRepetitionDir);
    else
        error('No Repetitions_XX directory found.');
    end

    % Initialize arrays to store distances and pairs
    num_matrices = states;
    d_geodesic = zeros(num_matrices*(num_matrices-1)/2, 1);
    d_pearson = zeros(num_matrices*(num_matrices-1)/2, 1);
    comparisons = cell(num_matrices*(num_matrices-1)/2, 2); % Cell array to store pairs

    index = 1;
    for i = 1:num_matrices
        for j = (i+1):num_matrices
            % Load matrices
            matrix1 = importdata(fullfile(matrixDir, ['FC_state' num2str(i) '.mat']));
            matrix2 = importdata(fullfile(matrixDir, ['FC_state' num2str(j) '.mat']));

            % Save matrices to .mat files
            save(fullfile(projectDir, 'code/matrix1.mat'), 'matrix1');
            save(fullfile(projectDir, 'code/matrix2.mat'), 'matrix2');

            cd(fullfile(projectDir, 'code/'));
            % Run Python code and capture output
            [status, output] = system(['python3 ' fullfile(projectDir, 'code/distance_metrics.py ') ' matrix1.mat matrix2.mat']);

            % Parse output to extract distances
            distances = str2double(strsplit(output));

            % Store distances
            d_geodesic(index) = distances(3);
            d_pearson(index) = distances(6);

            % Store pair of matrices
            comparisons{index, 1} = ['FC_state' num2str(i) '.mat'];
            comparisons{index, 2} = ['FC_state' num2str(j) '.mat'];

            index = index + 1;
        end
    end

    if length(d_geodesic) > 1 && length(d_pearson) > 1 
        d_geodesic_all(ii,1) = mean(d_geodesic,1);
        d_pearson_all(ii,1) = mean(d_pearson,1);
    else
        d_geodesic_all(ii,1) = d_geodesic;
        d_pearson_all(ii,1) = d_pearson;
        
    end
        
    % Save distances and comparisons to a MAT file
    %save(fullfile(matrixDir, 'distances.mat'), 'd_geodesic', 'd_pearson', 'comparisons');

end

system('rm matrix1.mat matrix2.mat')

mkdir(fullfile(resultsDir, 'figures/Distances'));
save(fullfile(resultsDir, 'figures/Distances/distances_all.mat'), 'd_geodesic_all', 'd_pearson_all');




% Plot the results for pearson dissimilarity and geodesic distance
x_values = K; % Define x-axis values from 2 to 15
y_values = d_pearson_all(1:(num_networks*2)); 

figure; 
plot(x_values, y_values, 'LineWidth', 8)
xlim([min(K) max(K)]); 
xticks(min(K):1:max(K))
yticks(round(min(y_values),2):0.03:round(max(y_values),2))
title('Pearson Dissimilarity', 'FontSize', 20)
ylabel('Mean Pearson Dissimilarity', 'FontSize', 20)
xlabel('Model States', 'FontSize', 20)
ax = gca;
ax.FontSize = 16; % Set tick label font size

saveas(gca, fullfile(resultsDir, 'figures/Distances/Pearson_Dissimilarity.bmp'));

x_values = K; % Define x-axis values from 2 to 15
y_values = d_geodesic_all(1:(num_networks*2)); 

figure; 
plot(x_values, y_values, 'LineWidth', 8)
xlim([min(K) max(K)]); 
xticks(min(K):1:max(K))
yticks(round(min(y_values),2):0.5:round(max(y_values),2))
title('Geodesic Distance', 'FontSize', 20)
ylabel('Mean Geodesic Distance', 'FontSize', 20); 
xlabel('Model States', 'FontSize', 20)
ax = gca;
ax.FontSize = 16; % Set tick label font size
saveas(gca, fullfile(resultsDir, 'figures/Distances/Geodesic_Distance.bmp'));


%% Correlation between vectors of mean activation

for ii = 1:length(K)
    
    states = K(ii);
    disp(states)
    
    % Create the path to the NumStates_XX directory
    numStatesDir = fullfile(resultsDir, ['figures/NumStates_' num2str(states)]);

    % List all directories in the NumStates_XX directory
    repetitionsDirs = dir(numStatesDir);
    repetitionsDirs = repetitionsDirs([repetitionsDirs.isdir]); % Filter out non-directory entries
    repetitionsDirs = repetitionsDirs(~ismember({repetitionsDirs.name}, {'.', '..'})); % Exclude '.' and '..'

    % Get the first directory found (assuming there's only one Repetitions_XX directory)
    if ~isempty(repetitionsDirs)
        firstRepetitionDir = repetitionsDirs(1).name;
        projectDir = fullfile(numStatesDir, firstRepetitionDir);
    else
        error('No Repetitions_XX directory found.');
    end
    
    load(fullfile(projectDir, 'MeanActivation.mat'))
    % Initialize arrays to store distances and pairs
    num_matrices = states;
    corr_mean = zeros(num_matrices*(num_matrices-1)/2, 1);
    p_mean = zeros(num_matrices*(num_matrices-1)/2, 1);
    comparisons = cell(num_matrices*(num_matrices-1)/2, 2); % Cell array to store pairs

    index = 1;
    for i = 1:num_matrices
        for j = (i+1):num_matrices
            
            [corr_mean(index), p_mean(index)] = corr(net_mean(:,i), net_mean(:,j));
            index = index + 1;
           
        end
    end
    
    figure; histogram(corr_mean,20); title(['Histogram correlations state ' num2str(states)])
    xline(mean(corr_mean,1), 'r--', 'LineWidth', 2);
    
    if length(corr_mean) > 1
        corr_mean_all(ii,1) = mean(abs(corr_mean),1);
        p_mean_all(ii,1) = mean(p_mean,1);
    else
        corr_mean_all(ii,1) = abs(corr_mean);
        p_mean_all(ii,1) = p_mean;
        
    end

    
    % Save comparisons to a MAT file
    %save(fullfile(projectDir, 'corr_mean.mat'), 'corr_mean', 'p_mean');

end
save(fullfile(resultsDir, 'figures/Distances/corr_all.mat'), 'corr_mean_all');

% Plot the results
x_values = K; % Define x-axis values from 2 to 12
y_values = corr_mean_all(1:(max(K)-1)); % Extract corresponding y-values from your data

figure; 
plot(x_values, y_values, 'LineWidth', 8)
xlim([min(K) max(K)]);
xticks(min(K):1:max(K))
yticks(round(min(y_values),2):0.2:round(max(y_values),2))
title('Correlation between mean activations', 'FontSize', 20)
ylabel('Absolute Mean Correlation', 'FontSize', 20); 
xlabel('Model States', 'FontSize', 20)
ax = gca;
ax.FontSize = 16; % Set tick label font size

saveas(gca, fullfile(resultsDir, 'figures/Distances/Correlations.bmp'));

close all;


