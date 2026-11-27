% Caterina Pedersini, October 2024, cpedersini86@gmail.com

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
K = (2:15); % no. states to be included in the analysis
pca = 10;
repetitions = 20;
T = 464; % number of timepoints

resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca), ['zeromean' num2str(zeromean)]);
load(fullfile(resultsDir, 'HMMs_GROUP.mat'));
load(fullfile(resultsDir, 'HMMs_options.mat'));
N = length(subjects);
  



%% Choose the repetition with the lowest free energy, for each model order

num_states_i = 6;

state_i = K(num_states_i);
k = state_i-1;


figure('Color','w','Position',[100 100 750 550]);
hold on;

% Extract Free Energy values
FE = cell2mat(fe);

% Mean and SD across repetitions
meanFE = mean(FE, 1, 'omitnan');
sdFE   = std(FE, 0, 1, 'omitnan');
lineColor = 0.5 * [0.2 0.45 0.7] + (1-0.5) * [1 1 1];

for repet_i = 1:repetitions
    plot(K, FE(repet_i,:), 'Color', lineColor, 'LineWidth', 0.9);
end

% Mean Free Energy
plot(K, meanFE, '-o', ...
    'Color', [0.15 0.15 0.15], ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', [0.2 0.45 0.7], ...
    'MarkerSize', 7);

% Formatting
ax = gca;
ax.FontSize = 24;
ax.FontName = 'Arial';
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.Box = 'off';
xticks(K);
xlabel('Number of states (K)', 'FontSize', 24);
ylabel('Free Energy', 'FontSize', 24);
xlim([min(K)-0.3 max(K)+0.3]);
grid off;

exportgraphics(gcf, ...
    fullfile('/Volumes/data_cate/Projects/NYUAD/dFC/temporal_filtering/code/Review_October2026/NumberOfStates/PredictivePerformance/PlotFE.jpg'), ...
    'Resolution',600);


% Plot Average Log-Likelihood across repetitions
figure('Color','w','Position',[100 100 750 550]);
hold on;

LL = cell2mat(ll_average);
% Use absolute values only if you want to reproduce the original figure
LL = abs(LL);

meanLL = mean(LL,1,'omitnan');
sdLL = std(LL,0,1,'omitnan');
lineColor = 0.5 * [0.2 0.45 0.7] + (1-0.5) * [1 1 1];

for repet_i = 1:repetitions
    plot(K,LL(repet_i,:), 'Color', lineColor, 'LineWidth', 0.9);
end

% Mean Free Energy
plot(K, meanLL, '-o', ...
    'Color', [0.15 0.15 0.15], ...
    'LineWidth', 2.5, ...
    'MarkerFaceColor', [0.2 0.45 0.7], ...
    'MarkerSize', 7);

% Axes formatting
ax = gca;
ax.FontSize = 24;
ax.FontName = 'Arial';
ax.LineWidth = 1.2;
ax.TickDir = 'out';
ax.Box = 'off';
xticks(K);
xlim([min(K)-0.3 max(K)+0.3]);
xlabel('Number of states (K)','FontSize',24);
ylabel('Average Log-Likelihood','FontSize',24);

grid off;

% Export
exportgraphics(gcf, ...
    fullfile('/Volumes/data_cate/Projects/NYUAD/dFC/temporal_filtering/code/Review_October2026/NumberOfStates/PredictivePerformance/PlotAvLL.jpg'), ...
    'Resolution',600);

% Mean avLL across model orders for each repetition
mean_ll = mean(LL,2,'omitnan');

for repet_i = 1:repetitions
    disp(['Mean log-likelihood repetition ' num2str(repet_i) ...
        ' = ' num2str(mean_ll(repet_i))]);
end

close all;



