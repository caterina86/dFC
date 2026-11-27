
% Caterina Pedersini,  cpedersini86@gmail.com

% Test to assess 
% 1. mean activation map,
% 2. correlation between state Gamma and head motion

clearvars; close all;

baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
DirData = fullfile(projectDir, 'derivatives/HMM');

mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR
addpath(genpath(fullfile(projectDir, 'code')))

% define parameters
zeromean = 0;
K = (2:15); % range of the number of states tested
pca = 10;
repetitions = 20;
T = 464; % number of timepoints

resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca), ['zeromean' num2str(zeromean)]);
load(fullfile(resultsDir, 'HMMs_GROUP.mat'));




%% Find the minimum value and its indices, based on the Free Energy
nStates = 7; % IMPORTANT: choose the specific model (number of states)
k = nStates - 1; 

[minValue_fe, minIndex_fe] = min(cell2mat(fe(:,:)));
[rowIndex_fe, colIndex_fe] = ind2sub(size(fe), minIndex_fe);
rep_i = rowIndex_fe(1,k); % repetition chosen based on the lowest free energy

disp(['repetition chosen ' num2str(rep_i)])
disp(['number states chosen ' num2str(K(k))])


for num_states_i = nStates
    
    
    resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca), ['zeromean' num2str(zeromean)]);
    load(fullfile(resultsDir, 'HMMs_GROUP.mat'));
    N = length(subjects); % no. subjects
  

    %% Find the minimum value and its indices, based on the Free Energy

    [minValue_fe, minIndex_fe] = min(cell2mat(fe(:,:)));
    [rowIndex_fe, colIndex_fe] = ind2sub(size(fe), minIndex_fe);

    rep_i = rowIndex_fe(1,k); % repetition chosen based on the lowest free energy
    disp(['repetition chosen ' num2str(rep_i)])
    disp(['number states chosen ' num2str(K(k))])

    load(fullfile(resultsDir, ['HMMs_r' num2str(rep_i) '_SUBJECTS.mat']));
    mkdir(fullfile(resultsDir, 'figures/'));


    %% Plot the mean activations for each state and model

    % group results
    if zeromean == 0
        net_mean = getMean(hmm{rep_i,k}); % get the mean activity of each state -> 100 components x 12 states
        figure;
        for k_i = 1:K(k)
            hold on;
            plot(net_mean(:,k_i))
        end
        hold off; % Release the hold on the plot
        % Create a legend based on the DisplayNames set in the loop
        legend('Location', 'best'); % You can adjust the location as needed
        %saveas(gca, fullfile(outputDir, 'MeanActivation.bmp'))
    end

    %save(fullfile(outputDir, 'MeanActivation.mat'), 'net_mean')


    mean_abs = mean(abs(net_mean), 1);
    mean_signed = mean(net_mean, 1);
    mean_norm = vecnorm(net_mean, 2, 1);
    
    results = table((1:nStates)', mean_signed', mean_abs', mean_norm', ...
        'VariableNames', {'State', 'SignedMean', ...
        'MeanAbsoluteValue', 'MeanVectorNorm'});
    
    disp(results)


    K = nStates;
    cov_state = cell(K,1);
    corr_state = cell(K,1);
    all_values = [];
    
    for s = 1:K
        [cov_state{s}, corr_state{s}] = getFuncConn(hmm{rep_i,k}, s);
    
        % Exclude diagonal for calculating common color limits
        tmp = cov_state{s};
        tmp(1:size(tmp,1)+1:end) = NaN;
        all_values = [all_values; tmp(:)];
    end
    
    % Common color scale across all states
    all_values = all_values(~isnan(all_values));
    clim_common = max(abs(all_values));
   

    meanAbsCov  = zeros(K,1);
    meanAbsCorr = zeros(K,1);
    normCov     = zeros(K,1);
    normCorr    = zeros(K,1);
    
    for s = 1:K
        C = cov_state{s};
        R = corr_state{s};
    
        % Upper triangle: unique off-diagonal connections
        idx = triu(true(size(C)), 1);
    
        meanAbsCov(s)  = mean(abs(C(idx)), 'omitnan');
        meanAbsCorr(s) = mean(abs(R(idx)), 'omitnan');
    
        % Overall magnitude of off-diagonal connections
        normCov(s)  = norm(C(idx), 2);
        normCorr(s) = norm(R(idx), 2);
    end
    
    results = table((1:K)', meanAbsCov, meanAbsCorr, normCov, normCorr, ...
        'VariableNames', {'State','MeanAbsCovariance','MeanAbsCorrelation', ...
                          'CovarianceNorm','CorrelationNorm'});
    
    disp(results);

    meanVariance = zeros(K,1);
    
    for s = 1:K
        meanVariance(s) = mean(diag(cov_state{s}), 'omitnan');
    end
    
    table((1:K)', meanVariance, ...
        'VariableNames', {'State','MeanParcelVariance'})


    
end



%% Correlation between FD and Gamma

projectDir = '/Volumes/data_cate/Projects/NYUAD/dFC/temporal_filtering';

subjects = importdata(fullfile(projectDir, ...
    'code/S07_AnalysisChronnectome/list_subjects_dataset2023_2017.txt'));

nSub = 83;
nT = 464;



%% 1. Load FD time series once

FD = cell(nSub,1);
motion = nan(nSub,1);

for sub_i = 1:nSub

    file1 = fullfile(projectDir, 'derivatives/A1/fmriprep', ...
        ['sub-' subjects{sub_i} ...
        '_task-rest_run-1_desc-confounds_timeseries.tsv']);

    file2 = fullfile(projectDir, 'derivatives/A1/fmriprep', ...
        ['sub-' subjects{sub_i} ...
        '_task-resting_desc-confounds_timeseries.tsv']);

    if exist(file1, 'file')
        tmpregressor = readtable(file1, 'FileType', 'delimitedtext');
    elseif exist(file2, 'file')
        tmpregressor = readtable(file2, 'FileType', 'delimitedtext');
    else
        error('No confounds file found for participant %s', subjects{sub_i});
    end

    fd = tmpregressor.framewise_displacement;
    nVol = numel(fd);

    if nVol == 472
        fd = fd(9:end);       % Remove first 8 volumes
    elseif nVol == 480
        fd = fd(17:end);      % Remove first 16 volumes
    else
        error('Unexpected number of volumes (%d) for participant %s', ...
            nVol, subjects{sub_i});
    end

    if numel(fd) ~= nT
        error(['FD length mismatch for participant %s: ' ...
            '%d volumes after trimming, expected %d'], ...
            subjects{sub_i}, numel(fd), nT);
    end

    FD{sub_i} = fd(:);
    motion(sub_i) = mean(fd, 'omitnan');

end



%% 2. Correlate FD with posterior probability for each state

rho_all = nan(nSub,nStates);
pval_all = nan(nSub,nStates);

pGroup = nan(nStates,1);
medianRho = nan(nStates,1);
nValid = zeros(nStates,1);

for ii = 1:nStates

    Gamma_selected = Gamma{rep_i, k};

    % Expected dimensions: (464 * 83) x 7
    assert(size(Gamma_selected,1) == nT*nSub, ...
        'Unexpected number of rows in Gamma for state model %d', ii);

    assert(size(Gamma_selected,2) == nStates, ...
        'Unexpected number of states in Gamma for model %d', ii);

    % Each column is a participant; each row is a time point
    Gamma_state = reshape(Gamma_selected(:,ii), nT, nSub);

    for s = 1:nSub

        gamma_s = Gamma_state(:,s);
        fd_s = FD{s};

        valid = isfinite(gamma_s) & isfinite(fd_s);

        if sum(valid) >= 3 && ...
                std(gamma_s(valid)) > 0 && std(fd_s(valid)) > 0

            [rho_all(s,ii), pval_all(s,ii)] = corr( ...
                gamma_s(valid), fd_s(valid), 'Type', 'Spearman');

        end
    end

    validSub = isfinite(rho_all(:,ii));
    nValid(ii) = sum(validSub);

    if nValid(ii) > 0
        medianRho(ii) = median(rho_all(validSub,ii));
        pGroup(ii) = signrank(rho_all(validSub,ii), 0);
    end

    fprintf('State %d: N = %d, median rho = %.3f, p = %.4f\n', ...
        ii, nValid(ii), medianRho(ii), pGroup(ii));

end


figure;
tiledlayout(1, nStates, 'TileSpacing', 'compact');
% Common symmetric y-axis limits across all states
maxRho = max(abs(rho_all(:)), [], 'omitnan');
yl = 1.1 * maxRho;
if isempty(yl) || ~isfinite(yl) || yl == 0
    yl = 1;
end

for ii = 1:nStates
    nexttile;

    rho = rho_all(:,ii);
    rho = rho(isfinite(rho));

    scatter(ones(size(rho)), rho, 22, 'filled', ...
        'MarkerFaceAlpha', 0.45);
    hold on;

    yline(0, '--k');
    yline(median(rho, 'omitnan'), '-r', 'LineWidth', 2.5);

    xlim([0.7 1.3]);
    ylim([-yl yl]);  % Identical scale and zero position
    xticks([]);
    title(sprintf('State %d', ii), 'FontSize', 16);
    grid on;

    if ii == 1
        ylabel('Spearman \rho with FD');
    end
end






%% 3. Correct for multiple comparisons across seven states

pHolm = nan(size(pGroup));
validP = isfinite(pGroup);

if any(validP)
    p = pGroup(validP);
    [pSorted, order] = sort(p);

    m = numel(pSorted);
    adjSorted = min(1, cummax((m - (1:m)' + 1) .* pSorted));

    adj = nan(size(p));
    adj(order) = adjSorted;
    pHolm(validP) = adj;
end

results = table((1:nStates)', nValid, medianRho, pGroup, pHolm, ...
    'VariableNames', {'State','N','MedianRho','RawP','HolmAdjustedP'});
disp(results);


% Median absolute correlation for each state
medianAbsRho = median(abs(rho_all), 1, 'omitnan');
for ii = 1:nStates
    fprintf('State %d: median |rho| = %.4f\n', ...
        ii, medianAbsRho(ii));
end



figure('Color','w');
% Absolute participant-level correlations
data = abs(rho_all);
boxplot(data, 'Labels', compose('%d', 1:nStates), ...
    'Symbol', '');
hold on;

% Overlay individual participants
rng(1);
for ii = 1:nStates
    vals = data(:,ii);
    vals = vals(isfinite(vals));
    x = ii + (rand(size(vals)) - 0.5) * 0.25;
    scatter(x, vals, 40, 'filled', 'MarkerFaceAlpha', 0.4);
end

% Highlight the median lines in red
h = findobj(gca, 'Tag', 'Median');
set(h, 'Color', 'r', 'LineWidth', 1.5);
xlabel('HMM state', 'FontSize', 18);
ylabel('Absolute Spearman correlation', 'FontSize', 18);

% Increase tick-number size and weight
set(gca, 'FontSize', 18);

ylim([0, max(data(:), [], 'omitnan') * 1.1]);
box on;

exportgraphics(gcf, ...
    fullfile('/Volumes/data_cate/Projects/NYUAD/dFC/temporal_filtering/code/Review_October2026/Correlation_Gamma_FD/Correlation_Mean_FD.jpg'), ...
    'Resolution',600);

