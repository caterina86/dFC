% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Script to extract the highest correlation between states belonging to
% different model orders.
% The calculation is based on the correlation between maps of mean
% activation and the pairs of states showing the highest correlations are
% indicated in comparison_max.mat
 
clearvars; close all;

baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
DirData = fullfile(projectDir, 'derivatives/HMM');
resultsDir = fullfile(DirData, 'PCA_REDUCED_15/10/zeromean0/figures');

mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR
addpath(genpath(fullfile(projectDir, 'code')))

% Define states to be compared
states_1 = 8; 
states_2 = 7;

% Get the first repetition directory for both states
[firstRepetitionDir1, net_mean1] = getFirstRepetition(resultsDir, states_1);
[firstRepetitionDir2, net_mean2] = getFirstRepetition(resultsDir, states_2);

% Load the mean activation data for both states
if exist(net_mean1, 'file')
    data1 = load(net_mean1);
else
    error(['MeanActivation.mat not found for state ' num2str(states_1)]);
end

if exist(net_mean2, 'file')
    data2 = load(net_mean2);
else
    error(['MeanActivation.mat not found for state ' num2str(states_2)]);
end


%% Distance between maps of mean activation

index = 1;

for ii = 1:size(data2.net_mean,2)
    
    for jj = 1:size(data1.net_mean,2)
        
        [r(index),p(index)] = corr(data2.net_mean(:,ii), data1.net_mean(:,jj));
        comparisons{index,1} = [num2str(ii)];
        comparisons{index,2} = [num2str(jj)];
        
        index = index + 1;
    end
    
    [r_max(ii), ind_max(ii)] = max(r);
    % extract the indeces of states more similar in terms of mean activation (higher correlation)
    comparison_max(ii,1) = str2num(comparisons{ind_max(ii),1});
    comparison_max(ii,2) = str2num(comparisons{ind_max(ii),2});

    index = 1;

end

% Plot the scatterplots of the mean activations of 2 states
ii = 3; 
jj = 1;

figure; 
scatter(data1.net_mean(:,ii), data2.net_mean(:,jj), 'filled')
lsline
xlabel(['Mean of state ' num2str(ii) ' Model ' num2str(states_1) ' states'], 'FontSize', 18)
ylabel(['Mean of state ' num2str(jj) ' Model ' num2str(states_2) ' states'], 'FontSize', 18)
text(-0.4, 0.4, ['r = ' num2str(r_max(1,ii))])



%% Distance between FC matrices -> NOT used in the paper

% index = 1;
% 
% for ii = 1:states_2
% 
%     for jj = 1:states_1
%     
%         matrix1 = importdata(fullfile(resultsDir, ['NumStates_' num2str(states_2)], firstRepetitionDir2, ['FC_state' num2str(ii) '.mat']));
%         matrix2 = importdata(fullfile(resultsDir, ['NumStates_' num2str(states_1)], firstRepetitionDir1, ['FC_state' num2str(jj) '.mat']));
%           
%         % Save matrices to .mat files
%         save(fullfile(projectDir, 'code/matrix1.mat'), 'matrix1');
%         save(fullfile(projectDir, 'code/matrix2.mat'), 'matrix2');
% 
%         cd(fullfile(projectDir, 'code/'));
%         % Run Python code and capture output
%         [status, output] = system(['python3 ' fullfile(projectDir, 'code/distance_metrics.py ') ' matrix1.mat matrix2.mat']);
% 
%         % Parse output to extract distances
%         distances = str2double(strsplit(output));
% 
%         % Store distances
%         d_geodesic(index) = distances(3);
%         d_pearson(index) = distances(6);
% 
%         % Store pair of matrices
%         comparisons{index,1} = [num2str(ii)];
%         comparisons{index,2} = [num2str(jj)];
% 
%         index = index + 1;
%     end
%     
%     [geodesic_min(ii), ind_geod_min(ii)] = min(d_geodesic);
%     [pearson_min(ii), ind_pear_min(ii)] = min(d_pearson);
%     
%     % extract the indeces of states more similar in terms of mean activation (higher correlation)
%     comparison_geodesic_min(ii,1) = str2double(comparisons{ind_geod_min(ii),1});
%     comparison_geodesic_min(ii,2) = str2double(comparisons{ind_geod_min(ii),2});
%     
%     comparison_pearson_min(ii,1) = str2double(comparisons{ind_pear_min(ii),1});
%     comparison_pearson_min(ii,2) = str2double(comparisons{ind_pear_min(ii),2});
%     
%     clear d_geodesic d_pearson
%     index = 1;
%     
% end
