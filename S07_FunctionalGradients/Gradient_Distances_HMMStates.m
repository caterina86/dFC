
%% Gradient analysis

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
zeromean = 0;
pca = 10; % 10 -> pca reduction; 0 -> no pca reduction
repetitions = 20;
reduced = 1; % 1 -> reduced number of volumes; 0 -> upsampling of dataset 2017

% number of volumes per subject
if reduced == 1
    T = 464; 
else
    T = 480;
end

if reduced == 1
    if pca == 0
        resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], ['zeromean' num2str(zeromean)]);
    else
        resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca), ['zeromean' num2str(zeromean)]);
    end        
else
    if pca == 0
        resultsDir = fullfile(DirData, ['PCA_UP_' num2str(max(K))], ['zeromean' num2str(zeromean)]);
    else
        resultsDir = fullfile(DirData, ['PCA_UP_' num2str(max(K))], num2str(pca),['zeromean' num2str(zeromean)]);
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
load(fullfile(resultsDir, ['HMMs_r' num2str(rep_i) '_SUBJECTS.mat']));
mkdir(fullfile(resultsDir, 'figures/'));
    
state_i = K(k); 
outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)]);
mkdir(outputDir);

model_i = 1; % 1 -> model-based FC; 2 -> empirical stationary FC

if model_i == 1
    label = 'Gradients_Model-based';
else
    label = 'Gradients_Empirical';
end

mkdir(fullfile(outputDir, label));


%%
cmap=zeros(100,3);
cmap(1:9,:) = repmat([0 0 1], 9, 1);  % Blue, visual 
cmap(51:58,:) = repmat([0 0 1], 8, 1);  % Blue, visual 
cmap(10:15,:) = repmat([0 1 0], 6, 1);  % Green, SomMo
cmap(59:66,:) = repmat([0 1 0], 8, 1);  % Green, SomMo
cmap(16:23,:) = repmat([1 0 0], 8, 1);  % Red, DorsAttn
cmap(67:73,:) = repmat([1 0 0], 7, 1);  % Red, DorsAttn
cmap(24:30,:) = repmat([1 1 0], 7, 1);  % Yellow, SalVentAttn
cmap(74:78,:) = repmat([1 1 0], 5, 1);  % Yellow, SalVentAttn
cmap(31:33,:) = repmat([1 0 1], 3, 1);  % Magenta, Limbic
cmap(79:80,:) = repmat([1 0 1], 2, 1);  % Magenta, Limbic
cmap(34:37,:) = repmat([0 1 1], 4, 1);  % Cyan, Cont
cmap(81:89,:) = repmat([0 1 1], 9, 1);  % Cyan, Cont
cmap(38:50,:) = repmat([0.75 0.25 0.75], 13, 1);  % Gray, DMN
cmap(90:100,:) = repmat([0.75 0.25 0.75], 11, 1);  % Gray, DMN

visual_idx = [1:9 51:58];
som_idx = [10:15 59:66];
dors_idx = [16:23 67:73];
sal_idx = [24:30 74:78];
limbic_idx = [31:33 79:80];
cont_idx = [34:37 81:89];
dmn_idx = [38:50 90:100];

q = 0.025;   % conservative threshold




%% Load the reference gradients based on Margulies (2016)

nROIs = 100;

atlasFile = fullfile(projectDir, 'Atlas/MNI/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz');
V = double(niftiread(atlasFile));  % parcel labels 1..100

gradFiles = { ...
    fullfile(outputDir, 'Gradients', 'volume.all.0.nii.gz'), ... % Gradient 1
    fullfile(outputDir, 'Gradients', 'volume.all.1.nii.gz')};    % Gradient 2

Gref_parcel = zeros(nROIs, 2);

for g = 1:2
    Gvol = double(niftiread(gradFiles{g}));
    
    for p = 1:nROIs
        mask = (V == p) & ~isnan(Gvol) & Gvol ~= 0;
        Gref_parcel(p,g) = mean(Gvol(mask));
    end
end



%% -------------------------------
% HMM-based FC + Gradient Analysis
% -------------------------------

% Parameters
patients_idx = 1:33;
controls_idx = 34:83;
groups = {patients_idx, controls_idx};
Nsub = numel(patients_idx) + numel(controls_idx);
nROIs = 100;         % number of ROIs
density_threshold = 0.10;  % top 10%
nGrad = 2;           % number of gradients to analyze

% Storage
Galign = cell(1,2);          % group gradients
FCz_all = zeros(nROIs, nROIs, Nsub);  % Fisher-Z FC for each subject



%% State-dependent

% -------------------------------
% State-dependent gradient analysis
% -------------------------------

% Network indices
visual_idx = [1:9 51:58];
som_idx    = [10:15 59:66];
dors_idx   = [16:23 67:73];
sal_idx    = [24:30 74:78];
limbic_idx = [31:33 79:80];
cont_idx   = [34:37 81:89];
dmn_idx    = [38:50 90:100];

networks   = {som_idx, dors_idx, sal_idx, limbic_idx, cont_idx, dmn_idx};
net_names  = {'SM','Dorsal','Salience','Limbic','Control','DMN'};

nStates = 7;   % number of HMM states
nGrad   = 2;   % number of gradients
Nsub    = 83;  % total subjects
groups   = {1:33, 34:83}; % patients / controls

% Pre-allocate
Galign_state = cell(2,nStates);


%% -------------------------------
% 1. Compute group gradients for each state
%% -------------------------------
for ii = 1:2   % group: patients / controls
    range = groups{ii};
    
    for k_i = 1:nStates
        nSub = numel(range);
        FCz_all = zeros(nROIs, nROIs, nSub); % Fisher-Z
        
        for s_idx = 1:nSub
            sub_i = range(s_idx);
            
            % Extract FC for this state
            [~, FC] = getFuncConn(HMMs_dualregr{k,sub_i}, k_i);
            FC(1:nROIs+1:end) = 0;       % zero diagonal
            FCz_all(:,:,s_idx) = atanh(FC);
        end
        
        % Group average FC
        FC_group = mean(FCz_all, 3);
        
        % Row-wise thresholding (top 10%)
        for r = 1:nROIs
            row = FC_group(r,:);
            thresh = prctile(row, 100*(1-density_threshold));
            row(row < thresh) = 0;
            FC_group(r,:) = row;
        end
        
        % Cosine similarity
        norms = sqrt(sum(FC_group.^2,2));
        Aff = (FC_group * FC_group') ./ (norms*norms' + eps);
        
        % Avoid disconnected nodes
        zero_rows = all(Aff==0,2);
        if any(zero_rows)
            Aff(zero_rows,:) = eps;
            Aff(:,zero_rows) = eps;
        end
        
        % Gradient decomposition
        G = GradientMaps('kernel','cs','approach','dm','alignment','pa');
        G = G.fit(Aff, 'reference', Gref_parcel);
        Galign_state{ii,k_i} = G;
    end
end


%% -------------------------------
% 2. Project subject FC onto state gradients & compute visual–network distances
%% -------------------------------
pvals_state = cell(nGrad,1); % store p-values per gradient

for k_i = 1:nStates
    Grad_subj = zeros(nROIs,nGrad,Nsub);
    
    % Compute FC per subject
    FCz_all = zeros(nROIs,nROIs,Nsub);
    for s = 1:Nsub
        [~, FC] = getFuncConn(HMMs_dualregr{k,s},k_i);
        FC(1:nROIs+1:end) = 0;
        FCz_all(:,:,s) = atanh(FC);
    end
    
    % Project onto the state specific group gradient.
    % The fact of using the state specific group gradient instead of the
    % stationary grouo gradient, allows to respond to this question: “What is the intrinsic functional geometry of this state?”
    for s = 1:Nsub
        FC = FCz_all(:,:,s);
        FCn = FC ./ (sqrt(sum(FC.^2,2)) + eps);
        
        if s <= numel(groups{1})
            Ggrp = Galign_state{1,k_i}.aligned{1};
        else
            Ggrp = Galign_state{2,k_i}.aligned{1};
        end

        for g = 1:nGrad
            Grad_subj(:,g,s) = FCn * Ggrp(:,g);
        end
    end
    
    % Compute visual–network distance
    for g_i = 1:nGrad
        grad_norm = zscore(squeeze(Grad_subj(:,g_i,:)),0,1); % normalize across ROIs per subject
        dist_all = zeros(Nsub, numel(networks));
        
        for net_j = 1:numel(networks)
            for s = 1:Nsub
                dist_all(s,net_j) = abs(mean(grad_norm(visual_idx,s)) - mean(grad_norm(networks{net_j},s)));
            end
            
            % Compare patients vs controls
            dist_pat = dist_all(1:numel(groups{1}),net_j);
            dist_ctrl = dist_all(numel(groups{1})+1:end,net_j);
            [pvals_state{g_i}(k_i,net_j),~,stats_states{g_i}(k_i,net_j)] = ranksum(dist_pat, dist_ctrl);
        end
        
        % plot
        figure('Color','w'); hold on;

        cmap = [0 0 0.8; 0.8 0 0];
        width = 0.6;
        nPairs = numel(networks);

        all_data = [];
        all_group = [];

        for j = 1:nPairs

            data_pat  = dist_all(1:numel(groups{1}), j);
            data_ctrl = dist_all(numel(groups{1})+1:end, j);

            pos_pat  = 2*(j-1) + 1;
            pos_ctrl = 2*(j-1) + 2;

            all_data  = [all_data; data_pat; data_ctrl];
            all_group = [all_group; ...
                         ones(numel(data_pat),1)*pos_pat; ...
                         ones(numel(data_ctrl),1)*pos_ctrl];
        end

        % --- SINGLE boxplot ---
        h = boxplot(all_data, all_group, ...
            'Positions', 1:(2*nPairs), ...
            'Colors', repmat(cmap, nPairs, 1), ...
            'Widths', width, ...
            'Symbol','');

        % Make boxes thicker
        boxes = findobj(h,'Tag','Box');
        set(boxes,'LineWidth',2);

        hold on;

        % --- Scatter points (separate loop, overlay only) ---
        for j = 1:nPairs

            data_pat  = dist_all(1:numel(groups{1}), j);
            data_ctrl = dist_all(numel(groups{1})+1:end, j);

            pos_pat  = 2*(j-1) + 1;
            pos_ctrl = 2*(j-1) + 2;

            scatter(pos_pat + randn(size(data_pat))*0.07, data_pat, 30, cmap(1,:), ...
                'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor','none');

            scatter(pos_ctrl + randn(size(data_ctrl))*0.07, data_ctrl, 30, cmap(2,:), ...
                'filled', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor','none');
        end

        % --- Axis formatting ---
        xticks(1.5:2:2*nPairs-0.5);
        xticklabels(net_names);

        ylabel(['Visual–Network distance (Gradient ' num2str(g_i) ')']);

        set(gca, 'LineWidth', 2, 'FontSize', 18, 'Box', 'off', 'TickDir', 'out');
        saveas(gca,fullfile(outputDir, label, ['Boxplot_Gradient' num2str(g_i) '_State' num2str(k_i) '.jpg']))

    end
end


close all;




%% FDR on all pvalues - Distances between visual network and other networks

all_p = [pvals_state{1}(:); pvals_state{2}(:)];

p_fdr_all = mafdr(all_p(:), 'BHFDR', true);

pvals_fdr = cell(size(pvals_state));

idx = 1;

% --- pvals{1} ---
n = numel(pvals_state{1}); % 42
pvals_fdr{1} = reshape(p_fdr_all(idx:idx+n-1), size(pvals_state{1}));
idx = idx + n;

% --- pvals{2} ---
n = numel(pvals_state{2}); % 42
pvals_fdr{2} = reshape(p_fdr_all(idx:idx+n-1), size(pvals_state{2}));
idx = idx + n;


for state_i = 1:7

    for ii = 1:6
    
        r(state_i,ii) = stats_states{1,2}(state_i,ii).zval / sqrt(N);
    
    end
end

