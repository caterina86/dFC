

% Statistical analysis to assess differences between patients and controls
% at the level of the chronnectome

% load paths
baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
DirData = fullfile(projectDir, 'derivatives/HMM');

mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR
addpath(genpath(fullfile(projectDir, 'code')))

K = (2:15); % no. states to be included in the analysis
num_states = 2:15; % IMPORTANT: choose the specific model (number of states)
k = num_states - 1; 
N = length(subjects); % no. subjects
zeromean = 0;


%% Load the data and Extract the switching rate and maxFO

DirPCA10 =  fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], '10', ['zeromean' num2str(zeromean)]);

% Load group-level results
PCA10 = load(fullfile(DirPCA10, 'HMMs_GROUP.mat'), ...
           'switchingRate', 'maxFO', 'lifeTime');
nRep = 20;
Nsub = 83;
K_values = 2:15;
nK = numel(K_values);

bestRep = [15, 8, 6, 19, 15, 2, 13, 14, 17, 13, 5, 20, 16, 20];

% Dimensions: subjects x model sizes
SR_PCA10 = nan(Nsub, nK);
MFO_PCA10 = nan(Nsub, nK);

for k_i = 1:nK

    % Select the repetition with the lowest Free Energy for this K
    r = bestRep(k_i);

    SR_PCA10(:, k_i) = ...
        PCA10.switchingRate{r, k_i}(:);

    MFO_PCA10(:, k_i) = ...
        PCA10.maxFO{r, k_i}(:);

end


%% Plot

figure('Color','w','Position',[100 100 1000 420]);

% Switching rate across model orders K = 2:15
subplot(1,2,1);
hold on;
for s = 1:Nsub
    plot(K_values, SR_PCA10(s,:), ...
        '-', 'Color', [0.8 0.8 0.8]);
end
for k_i = 1:numel(K_values)
    boxchart(K_values(k_i)*ones(Nsub,1), ...
        SR_PCA10(:,k_i), ...
        'BoxFaceColor', [0.20 0.55 0.85], ...
        'MarkerStyle', 'none');
end
xticks(K_values);
xlim([min(K_values)-0.5, max(K_values)+0.5]);
xlabel('Number of hidden states (K)');
ylabel('Switching rate');
set(gca, 'FontSize', 20);
box off;

subplot(1,2,2);
hold on;
for s = 1:Nsub
    plot(K_values, MFO_PCA10(s,:), ...
        '-', 'Color', [0.8 0.8 0.8]);
end
for k_i = 1:numel(K_values)
    boxchart(K_values(k_i)*ones(Nsub,1), ...
        MFO_PCA10(:,k_i), ...
        'BoxFaceColor', [0.20 0.55 0.85], ...
        'MarkerStyle', 'none');
end
xticks(K_values);
xlim([min(K_values)-0.5, max(K_values)+0.5]);
xlabel('Number of hidden states (K)');
ylabel('Max FO');
set(gca, 'FontSize', 20);
box off;





%% Significant differences between states

K_values = 2:15;
nK = numel(K_values);

% Store results for all pairwise comparisons
nComparisons = nK*(nK-1)/2;
K1 = zeros(nComparisons,1);
K2 = zeros(nComparisons,1);
MeanDiff = zeros(nComparisons,1);
tStat = zeros(nComparisons,1);
pRaw = zeros(nComparisons,1);

idx = 0;

for i = 1:nK-1
    for j = i+1:nK

        % Paired subject-level observations
        x =  MFO_PCA10(:,i);
        y =  MFO_PCA10(:,j);

        % Remove subjects with missing values in either model
        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);

        % Paired t-test
        [~, p, ~, stats] = ttest(x,y);

        idx = idx + 1;

        K1(idx) = K_values(i);
        K2(idx) = K_values(j);
        MeanDiff(idx) = mean(y-x);
        tStat(idx) = stats.tstat;
        pRaw(idx) = p;
    end
end

% Holm correction across all pairwise comparisons
pHolm = nan(size(pRaw));
[sortedP, order] = sort(pRaw);
m = numel(pRaw);

adjusted = (m - (1:m)' + 1).*sortedP;
adjusted = cummax(adjusted);
adjusted = min(adjusted,1);

pHolm(order) = adjusted;

% Results table
Results = table(K1,K2,MeanDiff,tStat,pRaw,pHolm, ...
    'VariableNames', ...
    {'K1','K2','MeanDifference','tStatistic','pRaw','pHolm'});

Results.Significant = Results.pHolm < 0.01;

disp(Results);


% %% Calculate variability across subjects using CV
% 
% CV_SR  = nan(1, nK);
% CV_MFO = nan(1, nK);
% 
% for k_i = 1:nK
% 
%     % Switching rate
%     x = SR_PCA10(:,k_i);
%     x = x(isfinite(x));
% 
%     CV_SR(k_i) = std(x, 0) / abs(mean(x));
% 
%     % Maximum fractional occupancy
%     x = MFO_PCA10(:,k_i);
%     x = x(isfinite(x));
% 
%     CV_MFO(k_i) = std(x, 0) / abs(mean(x));
% 
% end
% 
% %% Plot CV across model orders
% 
% figure('Color','w','Position',[100 100 900 400]);
% 
% subplot(1,2,1);
% plot(K_values, CV_SR, '-o', ...
%     'LineWidth', 2, 'MarkerSize', 6);
% hold on;
% xline(7, '--r', 'K = 7');
% xlabel('Number of hidden states (K)');
% ylabel('Coefficient of variation');
% title('Switching rate');
% xticks(K_values);
% box off;
% set(gca,'FontSize',12);
% 
% subplot(1,2,2);
% plot(K_values, CV_MFO, '-o', ...
%     'LineWidth', 2, 'MarkerSize', 6);
% hold on;
% xline(7, '--r', 'K = 7');
% xlabel('Number of hidden states (K)');
% ylabel('Coefficient of variation');
% title('Maximum fractional occupancy');
% xticks(K_values);
% box off;
% set(gca,'FontSize',12);



%% Switching-rate plateau analysis

K = K_values(:)';
SR = MFO_PCA10; %SR_PCA10; %MFO_PCA10;

% Group-level summary
meanSR = mean(SR, 1, 'omitnan');
nValid = sum(isfinite(SR), 1);
sdSR = std(SR, 0, 1, 'omitnan');
seSR = sdSR ./ sqrt(nValid);


%% Fit continuous piecewise linear model
% Require at least 3 model orders on either side
candidates = K(4:end-3);
SSE = nan(size(candidates));
betaAll = cell(size(candidates));

for i = 1:numel(candidates)

    bp = candidates(i);

    X = [ones(numel(K),1)', ...
         K(:)', ...
         max(0, K(:)-bp)'];

    % Correct dimensions
    X = [ones(numel(K),1), K(:), ...
         max(0, K(:)-bp)];

    y = meanSR(:);
    valid = isfinite(y);

    beta = X(valid,:) \ y(valid);
    residuals = y(valid) - X(valid,:)*beta;

    SSE(i) = sum(residuals.^2);
    betaAll{i} = beta;
end

[~, bestIdx] = min(SSE);
bestBP = candidates(bestIdx);
beta = betaAll{bestIdx};

slopeBefore = beta(2);
slopeAfter  = beta(2) + beta(3);

fprintf('Estimated breakpoint: K = %d\n', bestBP);
fprintf('Slope before breakpoint: %.5f per state\n', ...
    slopeBefore);
fprintf('Slope after breakpoint: %.5f per state\n', ...
    slopeAfter);



%% Bootstrap participants to estimate uncertainty
rng(1234);

nBoot = 2000;
nSub = size(SR,1);

bootBP = nan(nBoot,1);
bootSlopeBefore = nan(nBoot,1);
bootSlopeAfter = nan(nBoot,1);

for b = 1:nBoot

    % Resample participants, preserving each participant's
    % measurements across all model orders
    idx = randi(nSub, nSub, 1);

    bootMean = mean(SR(idx,:), 1, 'omitnan')';

    bootSSE = nan(size(candidates));
    bootBeta = cell(size(candidates));

    for i = 1:numel(candidates)

        bp = candidates(i);

        X = [ones(numel(K),1), K(:), ...
             max(0, K(:)-bp)];

        valid = isfinite(bootMean);

        if sum(valid) < 4
            continue
        end

        bb = X(valid,:) \ bootMean(valid);
        res = bootMean(valid) - X(valid,:)*bb;

        bootSSE(i) = sum(res.^2);
        bootBeta{i} = bb;
    end

    [minSSE, j] = min(bootSSE);

    if isfinite(minSSE)
        bootBP(b) = candidates(j);

        bb = bootBeta{j};
        bootSlopeBefore(b) = bb(2);
        bootSlopeAfter(b) = bb(2) + bb(3);
    end
end

% 95% percentile intervals
validBP = bootBP(isfinite(bootBP));
validAfter = bootSlopeAfter(isfinite(bootSlopeAfter));

fprintf('95%% bootstrap CI for breakpoint: K = %.1f to %.1f\n', ...
    prctile(validBP, 2.5), prctile(validBP, 97.5));

fprintf('95%% bootstrap CI for post-break slope: %.5f to %.5f\n', ...
    prctile(validAfter, 2.5), prctile(validAfter, 97.5));



%% Plot group mean and fitted segmented trend
Kfine = linspace(min(K), max(K), 200)';

Xfine = [ones(numel(Kfine),1), Kfine, ...
         max(0, Kfine-bestBP)];

fitFine = Xfine * beta;

figure('Color','w','Position',[100 100 1000 420]);
hold on;

% Mean switching rate and standard error
errorbar(K, meanSR, seSR, 'o', ...
    'Color', [0.20 0.55 0.85], ...
    'MarkerFaceColor', [0.20 0.55 0.85], ...
    'MarkerEdgeColor', [0.20 0.55 0.85], ...
    'MarkerSize', 8, ...
    'LineWidth', 1.5, ...
    'CapSize', 8);

% Piecewise linear fit
plot(Kfine, fitFine, '-', ...
    'Color', [0.20 0.55 0.85], ...
    'LineWidth', 2.5);

% Breakpoint
xline(bestBP, '--', ...
    'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.5, ...
    'FontSize', 18, ...
    'LabelVerticalAlignment', 'middle', ...
    'LabelHorizontalAlignment', 'left');

% Axes
xticks(K_values);
xlim([min(K_values)-0.5, max(K_values)+0.5]);

xlabel('Number of hidden states (K)', 'FontSize', 22);
ylabel('Max FO', 'FontSize', 22);
set(gca, 'FontSize', 22, ...
    'LineWidth', 1.2, ...
    'TickDir', 'out');

box off;
legend('Mean \pm SE', 'Piecewise linear fit', ...
    'Location', 'best', 'FontSize', 20);