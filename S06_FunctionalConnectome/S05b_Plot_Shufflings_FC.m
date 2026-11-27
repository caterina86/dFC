

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
shuffle_state = 6;

load(fullfile(outputStatDir, 'matrices_TimeAveragedFC.mat'))
input = matrices_static; % your real static FC (Fisher-z)
[m, n, Sub] = size(input);
triangle = triu(true(m,n),1);

% load the output from NBS
if shuffle_state == 1 || shuffle_state == 6
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
% % SHUFFLING OF THE FUNCTIONAL CONNECTOME
% % --------------------------------------------------------------

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFC_state' num2str(shuffle_state)]);

% load the data
load(fullfile(outputShuffle, 'tvals.mat'));
load(fullfile(outputShuffle, 'mean_vals.mat'));

tvals_perm_nan = tvals_perm;      % copy original
sig_rows = tvals_masked ~= 0;     % rows to keep
n_nonzero_staticFC = nnz(tvals_masked);

% Set non-significant rows to NaN
%tvals_perm_nan(~sig_rows, :) = NaN;
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

valid_vals = meanFC(nan_mask == 0, 1);
nan_vals = meanFC(nan_mask == 1, 1);
col2_val   = meanFC(:, 2);   % constant value
ref_mat = col2_val{1,1};       % your reference (matrix or vector)

% If it's a matrix, extract triangle
if ismatrix(ref_mat)
    N = size(ref_mat,1);
    triangle = triu(true(N), 1);
    ref_vec = ref_mat(triangle);
else
    ref_vec = ref_mat;    % already a vector
end

% extract the mean of patients
FC_allStates = cell(Sub, state_i);
for st = 1:state_i
    for sub_i = 1:Sub
        [~, FC_allStates{sub_i, st}] = getFuncConn( ...
            HMMs_dualregr{k, sub_i}, st);
    end
end

meanFC_test{1,1} = zeros(100);   % patients

if shuffle_state == 1 || shuffle_state == 6
    for s = 1:33
        meanFC_test{1,1} = meanFC_test{1,1} + FC_allStates{s,shuffle_state};
    end
    meanFC_test{1,1} = meanFC_test{1,1} / 33;
else
    for s = 34:83
        meanFC_test{1,1} = meanFC_test{1,1} + FC_allStates{s,shuffle_state};
    end
    meanFC_test{1,1} = meanFC_test{1,1} / 50;
end

% If it's a matrix, extract triangle
if ismatrix(meanFC_test)
    N = size(meanFC_test{1,1},1);
    triangle = triu(true(N), 1);
    ref2_vec = meanFC_test{1,1}(triangle);
else
    ref2_vec = meanFC_test;    % already a vector
end

% Plot the mean values
mean_valid = zeros(numel(valid_vals),1);
for i = 1:numel(valid_vals)
    M = valid_vals{i};
    mean_valid(i) = mean(M(triangle));   % mean of upper triangle
end
mean_other = mean(ref_vec);
mean_other2 = mean(ref2_vec);
group1 = mean_valid;          

% plot the histogram
if ~isempty(nan_vals)

    mean_nan_valid = zeros(numel(nan_vals),1);
    for i = 1:numel(mean_nan_valid)
        M = nan_vals{i};
        mean_nan_valid(i) = mean(M(triangle));   % mean of upper triangle
    end
    group2 = mean_nan_valid;          % vector of mean FCs (EB NoSign)
    
    color_sign   = [0.20 0.60 0.55];  % Teal: EB Sign
    color_nosign = [0.85 0.65 0.25];  % Gold: EB NoSign
    color_sign_dark   = color_sign * 0.65;
    color_nosign_dark = color_nosign * 0.65;

    figure('Units','normalized','Position',[0.1 0.1 0.55 0.55], 'Color','w');
    hold on;    
    % --- Histograms ---
    h1 = histogram(group1, ...
        'FaceColor', color_sign, ...
        'FaceAlpha', 0.6, ...
        'EdgeColor', color_sign_dark, ...
        'LineWidth', 1.0);  % EB Sign
    
    h2 = histogram(group2, ...
        'FaceColor', color_nosign, ...
        'FaceAlpha', 0.6, ...
        'EdgeColor', color_nosign_dark, ...
        'LineWidth', 1.0);  % EB NoSign
    
    ymax = max([h1.Values, h2.Values]) * 1.1;
        
    % --- Lines for Controls and Patients ---
    plot([mean_other mean_other], [0 ymax], 'k--', 'LineWidth', 3); % Controls mean
    plot([mean_other2 mean_other2], [0 ymax], 'k--', 'LineWidth', 3); % Patients mean
        
    % --- X-axis labels ---
    xticks(0.2:0.02:0.4);
    xlabel('Mean FC', 'FontSize', 20, 'FontWeight', 'bold');
    
    % --- Make axes thicker ---
    set(gca, 'LineWidth', 2, 'FontSize', 22, 'Box','off', 'TickDir','out');
    %legend({'EB Sign','EB NoSign'}, 'Location','northwest');
    saveas(gca, fullfile(outputShuffle, 'Histogram_Shuffle.jpg'))

else

    figure('Units','normalized','Position',[0.1 0.1 0.55 0.55], 'Color','w');
    hold on;    
    % --- Histograms ---
    h1 = histogram(group1, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha',0.6, 'EdgeColor','none'); % EB Sign
    % Controls and Patients at x = 3: we'll represent them as vertical lines
    % Convert their mean values to y-axis limits
    ymax = max(h1.Values) * 1.1;
    
    % --- Lines for Controls and Patients ---
    plot([mean_other mean_other], [0 ymax], 'k--', 'LineWidth', 3); % Controls mean
    plot([mean_other2 mean_other2], [0 ymax], 'k--', 'LineWidth', 3); % Patients mean
        
    % --- X-axis labels ---
    xticks(0.2:0.02:0.4);
    xlabel('FC Value', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Counts', 'FontSize', 20, 'FontWeight', 'bold');
    
    % --- Make axes thicker ---
    set(gca, 'LineWidth', 2, 'FontSize', 18, 'Box','off', 'TickDir','out');
    % Optional: add legend
    legend({'EB Sign'}, 'Location','northwest');
    saveas(gca, fullfile(outputShuffle, 'Histogram_Shuffle.jpg'))

end

%close all;


[maxVal, idx] = max(n_nonNaN_shuffling);

a = tvals_perm_nan(:,idx);

A = nan(m, n);
A(triangle) = a;
A(1:m+1:end) = 0;   % set diagonal to 0
A(isnan(A)) = 0;
nConn = nnz(A);

B = A + A';

dlmwrite(fullfile(outputStatDir, ['Results_meanTValues_State' num2str(shuffle_state) '_FC.edge']), B, 'delimiter', ' ');

