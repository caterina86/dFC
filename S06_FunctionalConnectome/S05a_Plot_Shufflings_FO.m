

% Code to plot the results based on the shuffling of values of FO or FC
% across participants (patients vs controls)

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

thresholds = 4;
sizes = {'Extent'}; %{'Intensity'}; %{'Extent'};

outputStatDir = fullfile(outputDir, 'NBS', sizes{1}, 'staticFC');
mkdir(outputStatDir);



%% PARAMETERS
perm_N = 500;
shuffle_state = 1;

load(fullfile(outputStatDir, 'matrices_TimeAveragedFC.mat'))
input = matrices_static; % your real static FC (Fisher-z)
[m, n, Sub] = size(input);
triangle = triu(true(m,n),1);

% load the output from NBS
if shuffle_state == 1 | shuffle_state == 6
    load(fullfile(outputStatDir, 'Results_TimeAveragedFC_EBvsSC.mat'))
else
    load(fullfile(outputStatDir, 'Results_TimeAveragedFC_SCvsEB.mat'))
end
tvalues = nbs_outputs.NBS.test_stat;    % t-stat matrix        
tvalues(tvalues < thresholds) = 0; % remove all values of t < thr
numComponents = numel(nbs_outputs.NBS.con_mat);
merged_mask = false(size(tvalues));      % initialize empty logical mask
    
for i = 1:numComponents
    comp = full(nbs_outputs.NBS.con_mat{i});    % convert to full
    merged_mask = merged_mask | (comp == 1);    % logical OR to merge
end
t_masked_static = tvalues .* merged_mask; % masked matrix
t_masked_static(1:size(t_masked_static,1)+1:end) = 0; % zero diagonal
t_masked_static(isnan(t_masked_static)) = 0; % replace NaNs
tvals_masked = t_masked_static(triangle);





%% --------------------------------------------------------------
% % FRACTIONAL OCCUPANCY
% % --------------------------------------------------------------

FO_original = FOdual{rep_i,k}(:, shuffle_state);
patients_FO = mean(FO_original(1:33,1));

if shuffle_state == 1 | shuffle_state == 6
    others_FO = mean(FO_original(1:33,1));
else
    others_FO = mean(FO_original(34:end,1));
end

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFO_state' num2str(shuffle_state)]);

% load the data
load(fullfile(outputShuffle, 'tvals.mat'));
load(fullfile(outputShuffle, 'mean_vals.mat'));

% Choose the tvals within the connections that are significant in the empirical static FC
tvals_perm_nan = tvals_perm;      % copy original
sig_rows = tvals_masked ~= 0;     % rows to keep
n_nonzero_staticFC = nnz(tvals_masked);

% Set non-significant rows to NaN
tvals_perm_nan(~sig_rows, :) = NaN;
tvals_perm_nan(tvals_perm_nan == 0) = NaN;
tvals_perm_nan(tvals_perm_nan < thresholds) = NaN; % remove all values of t < thr

for ii = 1:perm_N
    summary_tvals(1,ii) = mean(tvals_perm_nan(:, ii), 'omitnan');
end

n_nonNaN_shuffling = sum(~isnan(tvals_perm_nan), 1);
perc_edges = (max(n_nonNaN_shuffling)/n_nonzero_staticFC)*100;
disp(['Percentage significant connections = ' num2str(perc_edges)])

% Create a binary mask for the NaNs
nan_mask = isnan(summary_tvals);   % 1 where NaN, 0 where valid
perc_nan = (sum(nan_mask)/perm_N)*100;
disp(['Percentage of NaN = ' num2str(perc_nan)])

if shuffle_state == 1 | shuffle_state == 6
    indx1 = 1;
    indx2 = 2;
else
    indx1 = 2;
    indx2 = 1;
end
valid_vals = value_mean(nan_mask == 0, indx1);
nan_vals   = value_mean(nan_mask == 1, indx1);
col2_val = value_mean(1, indx2);
col3_val = others_FO;


if ~isempty(nan_vals)

    figure('Units','normalized','Position',[0.1 0.1 0.55 0.55], 'Color','w');
    hold on;    
    % --- Histograms ---
    h1 = histogram(valid_vals, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha',0.6, 'EdgeColor','none'); % EB Sign
    h2 = histogram(nan_vals, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha',0.6, 'EdgeColor','none'); % EB NoSign
    % Controls and Patients at x = 3: we'll represent them as vertical lines
    % Convert their mean values to y-axis limits
    ymax = max([h1.Values, h2.Values]) * 1.1;
    
    % --- Lines for Controls and Patients ---
    plot([col2_val col2_val], [0 ymax], 'k--', 'LineWidth', 3); % Controls mean
    plot([col3_val col3_val], [0 ymax], 'k--', 'LineWidth', 3); % Patients mean
        
    % --- X-axis labels ---
    xticks(0:0.02:0.2);
    xlabel('FO Value', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Counts', 'FontSize', 20, 'FontWeight', 'bold'); xlim([0 0.2])
    
    % --- Make axes thicker ---
    set(gca, 'LineWidth', 2, 'FontSize', 18, 'Box','off', 'TickDir','out');
    % Optional: add legend
    legend({'EB Sign','EB NoSign'}, 'Location','northwest');
    saveas(gca, fullfile(outputShuffle, 'Histogram_Shuffle.jpg'))

else

    figure('Units','normalized','Position',[0.1 0.1 0.55 0.55], 'Color','w');
    hold on;    
    % --- Histograms ---
    h1 = histogram(nan_vals, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha',0.6, 'EdgeColor','none'); % EB Sign
    % Controls and Patients at x = 3: we'll represent them as vertical lines
    % Convert their mean values to y-axis limits
    ymax = max([h1.Values]) * 1.1;
    
    % --- Lines for Controls and Patients ---
    plot([col2_val col2_val], [0 ymax], 'k--', 'LineWidth', 3); % Controls mean
    plot([col3_val col3_val], [0 ymax], 'k--', 'LineWidth', 3); % Patients mean
        
    % --- X-axis labels ---
    xticks(0:0.02:0.2);
    xlabel('FO Value', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Counts', 'FontSize', 20, 'FontWeight', 'bold'); xlim([0 0.2])
    
    % --- Make axes thicker ---
    set(gca, 'LineWidth', 2, 'FontSize', 18, 'Box','off', 'TickDir','out');
    % Optional: add legend
    legend({'EB No Sign'}, 'Location','northwest');
    saveas(gca, fullfile(outputShuffle, 'Histogram_Shuffle.jpg'))

end

% Plot the difference of mean values
if shuffle_state == 1 | shuffle_state == 6
    mean_X = value_mean(:,indx1);   % patients for each permutation
    mean_Y = mean(FOdual{rep_i,k}(34:end, shuffle_state));   % controls (constant)
    mean_diff = mean_X - mean_Y;
else
    mean_X = value_mean(:,indx1);   
    mean_Y = mean(FOdual{rep_i,k}(1:33, shuffle_state));   % controls (constant)
    mean_diff = mean_X - mean_Y;
end

% Scatterplot showing the mean difference as a function of valid or nan values
figure('Units','normalized','Position',[0.1 0.1 0.6 0.6], 'Color','w');
scatter(1:length(mean_diff), mean_diff, 60, ~nan_mask, 'filled');
colormap([0 0 1; 1 0 0]); % red=valid, blue=NaN
xlabel('Permutation', 'FontSize', 20);
ylabel('Mean difference FO', 'FontSize', 20, 'FontWeight', 'bold');
colorbar('Ticks',[0 1],'TickLabels',{'EB NoSign', 'EB Sign'}, 'FontSize', 14);
yline(0,'--k');
% --- Make axes thicker ---
set(gca, 'LineWidth', 2, 'FontSize', 24, 'Box','off', 'TickDir','out');
saveas(gca, fullfile(outputShuffle, 'Scatteplot_Shuffle.jpg'))


x = value_mean(:, indx1);      % mean FO or FC
y = summary_tvals';            % t-values
valid_idx = ~isnan(x) & ~isnan(y);
r = corr(x(valid_idx), y(valid_idx));

figure('Units','normalized','Position',[0.1 0.1 0.6 0.6], 'Color','w');
scatter(value_mean(:,indx1), summary_tvals, 80, 'b', 'filled');
% --- Add correlation value inside plot ---
text(0.05, 0.95, sprintf('r = %.3f', r), ...
    'Units', 'normalized', ...   % position in axes units (0–1)
    'FontSize', 30, ...
    'FontWeight', 'bold', ...
    'Color', 'k');
xlabel('Mean EB FO (permutation)');
ylabel('Mean NBS t-value', 'FontSize', 20, 'FontWeight', 'bold'); 
title('Contribution of FO to network differences');
% --- Make axes thicker ---
set(gca, 'LineWidth', 2, 'FontSize', 24, 'Box','off', 'TickDir','out');
% --- Make boxplot lines thicker ---
saveas(gca, fullfile(outputShuffle, 'Scatteplot_tValues_meanFO.jpg'))


% % Statistics
% group1 = mean_diff(~nan_mask);    % Valid group
% group2 = mean_diff(nan_mask);   % NaN group
% group3 = FOdual{rep_i,k}(34:end, shuffle_state); % controls
% 
% % Significant difference between the mean values of valid and Nan
% p = ranksum(group2, group1);
% 
% % Effect size (Cohen’s d): measure of the effect size
% pooled_std = std([group1; group2]);  % SD across both groups together
% d = (mean(group1) - mean(group2)) / pooled_std;


close all;


[maxVal, idx] = max(n_nonNaN_shuffling);

a = tvals_perm_nan(:,idx);

A = nan(m, n);
A(triangle) = a;
A(1:m+1:end) = 0;   % set diagonal to 0
A(isnan(A)) = 0;
nConn = nnz(A);

B = A + A';

dlmwrite(fullfile(outputShuffle, 'Results_meanTValues.edge'), B, 'delimiter', ' ');

