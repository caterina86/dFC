%% Cross-validated HMM model selection


% Select the directory of the dFC folder
baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
cd(fullfile(projectDir, 'code/'));

parcels = 1:100; % number of parcels of the atlas used
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));

N = length(subjects); % no. subjects
f = cell(N,1); 
T = cell(N,1);

% dependency
mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR

derivativesDir = fullfile(projectDir, 'derivatives');
DirData = fullfile(derivativesDir, 'HMM');


for j=1:N
    my_string = ['sub-' subjects{j}];
    f{j} = fullfile(DirData, ['sub-' subjects{j}], ['sub-' subjects{j} '_data_Atlas_reduced.txt']);           
    load(f{j}); T{j} = 464;
end 

nSubjects = numel(T);
nFolds = 5;
repetitions = 20;

pca = 10; 
covtype = 'full';
zeromean = 0;

% HMM-MAR options
options = struct();
options.order = 0;
options.covtype = covtype;
options.zeromean = zeromean;
options.verbose = 1;
options.standardise = 1;
options.inittype = 'HMM-MAR';
options.cyc = 500;
options.initcyc = 10;
options.initrep = 3;
options.pca = pca(1);



%% Cross-validation for HMM model order

K_values = 2:15;

% Subject-level cross-validation folds
rng(1234);
cv = cvpartition(nSubjects, 'KFold', nFolds);

% Store results
cvLL = nan(numel(K_values), nFolds);
meanCVLL = nan(numel(K_values), 1);
stdCVLL = nan(numel(K_values), 1);


% Cross-validate each model order
for k_idx = 1:numel(K_values)

    K = K_values(k_idx);

    fprintf('\nCross-validating K = %d\n', K);
    options_cv = options;
    options_cv.K = K;
    options_cv.covtype = 'full';
    options_cv.cvfolds = cv;

    [meanCVLL(k_idx), cvLL(k_idx,:)] = cvhmmmar(f, T, options_cv);

    stdCVLL(k_idx) = std(cvLL(k_idx,:), 0, 2);

end

% Results table
resultsCV = table(K_values(:), meanCVLL, stdCVLL, 'VariableNames', {'K', 'MeanTestLogLikelihood', 'SDAcrossFolds'});
disp(resultsCV);
% Save results
save('HMM_crossvalidation_K2_15.mat', ...
     'resultsCV', 'cvLL', 'K_values', 'cv');


% Plot
figure('Color', 'w', 'Position', [100 100 750 500]);
errorbar(K_values, meanCVLL, stdCVLL, ...
    '-o', 'LineWidth', 1.5, 'MarkerSize', 7);
xlabel('Number of hidden states (K)', 'FontSize', 22);
ylabel('Mean test log-likelihood per time point', 'FontSize', 22);
xticks(K_values);
xlim([1.5 15.5]);
box off;
grid on;





%% Multiple repetitions for each fold


nCVrepeats = 50;

% Results: K x fold x repetition
cvLL_all = nan(numel(K_values), nFolds, nCVrepeats);
meanCVLL_all = nan(numel(K_values), nCVrepeats);

for r = 1:nCVrepeats

    fprintf('\nCV repetition %d/%d\n', r, nCVrepeats);
    % New subject-level partition for each repetition
    rng(1234 + r);
    cv = cvpartition(nSubjects, 'KFold', nFolds);

    for k_idx = 1:numel(K_values)
        K = K_values(k_idx);
        fprintf('  K = %d\n', K);
        options_cv = options;
        options_cv.K = K;
        options_cv.cvfolds = cv;
        [meanCVLL, foldLL] = cvhmmmar(f, T, options_cv);
        cvLL_all(k_idx,:,r) = foldLL(:)';
        meanCVLL_all(k_idx,r) = meanCVLL;
    end
end

% Summarize performance across repetitions
meanAcrossRepeats = mean(meanCVLL_all, 2, 'omitnan');
sdAcrossRepeats = std(meanCVLL_all, 0, 2, 'omitnan');

% How often each K achieved the highest mean test log-likelihood
bestK = nan(nCVrepeats,1);

for r = 1:nCVrepeats
    scores = meanCVLL_all(:,r);
    valid = isfinite(scores);
    if any(valid)
        validIdx = find(valid);
        [~, idx] = max(scores(valid));
        bestK(r) = K_values(validIdx(idx));
    end
end

resultsCV = table(K_values(:), meanAcrossRepeats, ...
    sdAcrossRepeats, 'VariableNames', ...
    {'K','MeanTestLogLikelihood','SDAcrossRepeats'});

disp(resultsCV);
disp('Frequency with which each K achieved the highest score:');
tabulate(bestK);

save('HMM_repeatedCV_K2_15.mat', ...
    'resultsCV','cvLL_all','meanCVLL_all', ...
    'K_values','bestK');


% Plot
figure('Color', 'w', 'Position', [100 100 1000 650]);
hold on;

errorbar(K_values, meanAcrossRepeats, sdAcrossRepeats, ...
    '-o', ...
    'LineWidth', 2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.2 0.45 0.7], ...
    'MarkerEdgeColor', [0.2 0.45 0.7], ...
    'Color', [0.2 0.45 0.7], ...
    'CapSize', 10);

xlabel('Number of hidden states (K)', ...
    'FontSize', 24, 'FontWeight', 'bold');
ylabel('Mean test log-likelihood', ...
    'FontSize', 24, 'FontWeight', 'bold');

xticks(K_values);
xlim([1.5 15.5]);
set(gca, 'FontSize', 24, 'LineWidth', 1.2, ...
    'TickDir', 'out', 'TickLength', [0.008 0.008], ...
    'FontName', 'Arial');
box off;
grid off;

ax = gca;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';
ax.Layer = 'top';

hold off;

exportgraphics(gcf, ...
    fullfile('/Volumes/data_cate/Projects/NYUAD/dFC/temporal_filtering/code/Review_October2026/NumberOfStates/PredictivePerformance/Crossvalidation_2_15_50permutations.jpg'), ...
    'Resolution',600);

%%

for k_idx = 1:numel(K_values)

    K = K_values(k_idx);
    scores = reshape(cvLL_all(k_idx,:,:), [], 1);

    fprintf(['K=%2d | finite folds: %d/%d | ' ...
             'mean fold LL: %.5f | SD: %.5f\n'], ...
        K, sum(isfinite(scores)), numel(scores), ...
        mean(scores, 'omitnan'), std(scores, 'omitnan'));
end


% Mean and SD across CV repetitions
meanAcrossRepeats = mean(meanCVLL_all, 2, 'omitnan');
sdAcrossRepeats   = std(meanCVLL_all, 0, 2, 'omitnan');

% Standard error (descriptive)
nValid = sum(isfinite(meanCVLL_all), 2);
seAcrossRepeats = sdAcrossRepeats ./ sqrt(nValid);

% Summary table
summaryTable = table(K_values(:), meanAcrossRepeats, ...
    sdAcrossRepeats, seAcrossRepeats, ...
    'VariableNames', {'K','MeanCVLL','SD','SE'});

disp(summaryTable);


nK = numel(K_values);

deltaMean = nan(nK-1,1);
deltaSD   = nan(nK-1,1);
pValues   = nan(nK-1,1);
nPairs    = zeros(nK-1,1);

for k = 1:nK-1

    % Paired values from the same CV repetitions
    a = meanCVLL_all(k,:)';
    b = meanCVLL_all(k+1,:)';

    valid = isfinite(a) & isfinite(b);
    delta = b(valid) - a(valid);

    nPairs(k) = numel(delta);
    deltaMean(k) = mean(delta);
    deltaSD(k) = std(delta);

    % Paired nonparametric test
    pValues(k) = signrank(b(valid), a(valid));
end

% Holm correction for the 13 adjacent comparisons
pHolm = nan(size(pValues));
validP = isfinite(pValues);
p = pValues(validP);
m = numel(p);

[pSorted, order] = sort(p);
adjustedSorted = min(1, cummax((m-(1:m)'+1) .* pSorted));
adjusted = nan(size(p));
adjusted(order) = adjustedSorted;
pHolm(validP) = adjusted;

% Results table
comparisonTable = table( ...
    K_values(1:end-1)', K_values(2:end)', ...
    deltaMean, deltaSD, nPairs, pValues, pHolm, ...
    'VariableNames', {'K1','K2','MeanDifference', ...
    'SDDifference','NRepetitions','pRaw','pHolm'});

disp(comparisonTable);