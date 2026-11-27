% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Script to extract the ratio Coefficient of Variation, for all models
% tested (from the repetition with the lowest FE).


clearvars; close all;

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
num_states = 7; % choose the specific model (number of states)
k = num_states-1; 

resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca), ['zeromean' num2str(zeromean)]);
load(fullfile(resultsDir, 'HMMs_GROUP.mat'));



%% Extract the coefficient of uncertainty -> lowest values indicate less uncertainty
% in the estimation of the mean activation maps

% Find the minimum value and its indices, based on the Free Energy
for k_i = 1:max(K-1)
    [minValue_fe, minIndex_fe] = min(cell2mat(fe(:,k_i)));
    [rowIndex_fe(k_i,1), colIndex_fe] = ind2sub(size(fe(:,k_i)), minIndex_fe);
end

state_i = 1;

% Extract the CV coefficient -> lower ratio indicate lower uncertainty
for repet_i = 1:size(rowIndex_fe,1)
    
        rep_i = rowIndex_fe(repet_i);
            
        K_ii(repet_i,1) = hmm{rep_i,state_i}.K;

        for k_i = 1:K_ii(repet_i,1)

            omega = getFuncConn(hmm{rep_i,state_i},k_i); % variance/covariance matrix
            J = 1./omega; 

            % Update the formula
            CV_values(:, k_i) = 100 ./ sqrt(diag(J));
            mean_CV_by_state_all(state_i,k_i) = mean(CV_values(:,k_i));
            std_CV_by_state_all(state_i,k_i) = std(CV_values(:,k_i));

        end
    
    for ii = 1:size(mean_CV_by_state_all,1)
        ind = mean_CV_by_state_all(ii,:) ~= 0;
        ratio(ii,1) = mean(std_CV_by_state_all(ii,ind),2)/mean(mean_CV_by_state_all(ii,ind),2);
    end
    
    state_i = state_i + 1;
end

% Plot the results
x_values = K; 
y_values = ratio(1:(max(K)-1)); % Extract corresponding y-values from your data

figure; 
plot(x_values, y_values, 'LineWidth', 8)
xlim([min(K) max(K)]) ; 
xticks(min(K):1:max(K))
yticks(round(min(y_values),2):0.005:round(max(y_values),2))
title('Correlation of Variation', 'FontSize', 20)
ylabel('Ratio CV', 'FontSize', 20); 
xlabel('Model States', 'FontSize', 20)
ax = gca;
ax.FontSize = 16; % Set tick label font size

saveas(gca, fullfile(resultsDir, 'figures/Distances/RatioCV.bmp'));


% calculate AIC and BIC for HMM (https://rdrr.io/cran/HMMpa/man/AIC_HMM.html)
% m = k; % number of states
% z = 2; % number of parameters -> for the normal distribution (mean and standard deviation))    
% p = m^2 + z*m - 1; % penalty
% AIC = -2 * cell2mat(ll_average(repet_i,m)) + 2*p;
% BIC = -2 * cell2mat(ll_average(repet_i,m)) + p * log(464);
 