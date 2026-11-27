

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



%% -------------------------------
% 1. Compute static FC for each subject
%% -------------------------------
for sub_i = 1:Nsub
    % Determine subject group
    if sub_i <= numel(patients_idx)
        grp = 1; 
    else
        grp = 2; 
    end

    % Extract number of states
    num_states = size(HMMs_dualregr{grp,sub_i},1);  % or get from HMM structure
    FC_state = zeros(nROIs,nROIs,num_states);
    
    if model_i == 1
        % --- Extract model-based stationary FC ---
        for st = 1:num_states
            [~, FC_state(:,:,st)] = getFuncConn(HMMs_dualregr{grp,sub_i}, st); 
        end
        
        % FO weighting
        FO = FOdual{rep_i,k}(sub_i, 1:num_states);   % ensure correct indexing
        FC_weighted = zeros(nROIs);
        for st = 1:num_states
            FC_weighted = FC_weighted + FO(st) * FC_state(:,:,st);
        end
        
        % Zero diagonal and Fisher Z
        FC_weighted(1:nROIs+1:end) = 0;
        FCz_all(:,:,sub_i) = atanh(FC_weighted);

    else
        % --- Extract empirical stationary FC ---
        f = fullfile(DirData, ['sub-' subjects{sub_i}], ['sub-' subjects{sub_i} '_data_Atlas_reduced.txt']);
        timecourses = load(f);          % nTime x nROIs 
        % Stationary FC
        FC = corr(timecourses);         % ROI x ROI
        FC(1:nROIs+1:end) = 0;          % zero diagonal
        
        % Fisher Z
        FCz_all(:,:,sub_i) = atanh(FC);
    end    

end



%% -------------------------------
% 2. Compute group-average FC and gradients
%% -------------------------------
for ii = 1:2
    idx = groups{ii};
    FC_group = mean(FCz_all(:,:,idx),3);  % group mean
    
    % Row-wise thresholding (top 10%)
    for r = 1:nROIs
        row = FC_group(r,:);
        thresh = prctile(row, 100*(1-density_threshold));
        row(row < thresh) = 0;
        FC_group(r,:) = row;
    end
    
    % Cosine affinity
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
    G = G.fit(Aff, 'reference', Gref_parcel);  % reference alignment
    Galign{ii} = G;
end

Galign_patients = Galign{1};
Galign_controls = Galign{2};



%% -------------------------------
% Visualize the gradients on the surface
%% -------------------------------

labeling = load_parcellation('schaefer', nROIs);
labeling = labeling.schaefer_100;
[surf_lh, surf_rh] = load_conte69();

h = plot_hemispheres(labeling, {surf_lh,surf_rh});

gradient_in_euclidean(Gref_parcel(:,1:2),{surf_lh,surf_rh},labeling,cmap,0);
gradient_in_euclidean(Galign_patients.aligned{1}(:,1:2),{surf_lh,surf_rh},labeling,cmap,0);
gradient_in_euclidean(Galign_controls.aligned{1}(:,1:2),{surf_lh,surf_rh},labeling,cmap,0);
close all;

% Visualize the gradient on the surface
grad_parcel = Galign_patients.aligned{1}(:,2); %Gref_parcel(:,2);  % 100 × 1
parc = load_parcellation('schaefer', 100);
% Convert parcel values → vertex values
grad_vertex = parcel2full(grad_parcel, parc.schaefer_100);

% --- Color limits ---
cmin = -30; %min(grad_parcel);
cmax = 30; %max(grad_parcel);
nColors = 100;
half = nColors/2;
% Negative side: dark blue at 0 -> light azzurro at min
neg = [linspace(0,0.23,half)', linspace(0,0.51,half)', linspace(0,0.92,half)'];
neg = flipud(neg);  % dark at 0, light at negative extreme
% Positive side: dark red at 0 -> light yellow at max
pos = [linspace(0.5,1,half)', linspace(0,1,half)', linspace(0,0,half)']; % dark red -> yellow
% Combine
cmap = [neg; pos];

% --- Plot ---
figure('Color','w','Position',[100 100 900 350])
% LEFT hemisphere
subplot(1,2,1)
trisurf(surf_lh.tri, surf_lh.coord(1,:), surf_lh.coord(2,:), surf_lh.coord(3,:), ...
        grad_vertex(1:size(surf_lh.coord,2)), 'EdgeColor','none');
view([-90 0]); axis off; axis equal
lighting gouraud; camlight
title('Left hemisphere')

% RIGHT hemisphere
subplot(1,2,2)
trisurf(surf_rh.tri, surf_rh.coord(1,:), surf_rh.coord(2,:), surf_rh.coord(3,:), ...
        grad_vertex(size(surf_lh.coord,2)+1:end), 'EdgeColor','none');
view([90 0]); axis off; axis equal
lighting gouraud; camlight
title('Right hemisphere')

% --- Apply colormap and color limits ---
colormap(cmap)
colorbar
clim([cmin cmax])



%% -------------------------------
% 3. Project subject FC onto group gradients
%% -------------------------------
Grad_subj = zeros(nROIs, nGrad, Nsub);

for sub_i = 1:Nsub
    FC = FCz_all(:,:,sub_i);  % subject FC
    FCn = FC ./ (sqrt(sum(FC.^2,2)) + eps);  % row-normalized
    
    % Select group gradient
    if sub_i <= numel(patients_idx)
        Ggrp = Galign_patients.aligned{1};
    else
        Ggrp = Galign_controls.aligned{1};
    end
    
    % Projection onto first nGrad gradients
    for g = 1:nGrad
        Grad_subj(:,g,sub_i) = FCn * Ggrp(:,g);  
    end
end


%% -------------------------------
% 4. Gradient range (segregation) and stats
%% -------------------------------
range_subj = zeros(Nsub, nGrad);

for sub_i = 1:Nsub
    for g = 1:nGrad
        Gz = zscore(Grad_subj(:,g,sub_i));  % normalize across ROIs
        range_subj(sub_i,g) = max(Gz) - min(Gz);
    end
end

% Compare groups
for g = 1:nGrad
    patients_range = range_subj(patients_idx,g);
    controls_range = range_subj(controls_idx,g);
    
    [pval_range(g,1),h,stats] = ranksum(patients_range, controls_range);
    fprintf('Gradient %d: p = %.4f, z = %.2f\n', g, pval_range(g,1), stats.zval);
    
    % Effect size
    Ntot = numel(patients_range) + numel(controls_range);
    r = stats.zval / sqrt(Ntot);
    fprintf('Effect size r = %.3f\n', r);
    
    % Boxplot
    figure
    boxplot([patients_range; controls_range], [ones(numel(patients_range),1); 2*ones(numel(controls_range),1)], ...
            'Labels', {'Patients','Controls'})
    ylabel(sprintf('Gradient %d range', g))
end



%% Distance between visual and all the other netowrks, for each gradient

for g_i = 1:2

    Gradient = g_i;

    % --- Network indices ---
    networks = {som_idx, dors_idx, sal_idx, limbic_idx, cont_idx, dmn_idx};
    net_names = {'SM','Dorsal','Salience','Limbic','Control','DMN'};
    visual = visual_idx;
    
    % --- Gradient of interest ---
    grad1_subj = squeeze(Grad_subj(:,Gradient,:));  % 100 × nSub
    nSub = size(grad1_subj,2);
    
    % --- Normalize gradients per subject ---
    grad_norm = zeros(size(grad1_subj));
    for s = 1:nSub
        grad_norm(:,s) = zscore(grad1_subj(:,s));
    end

    % --- Compare visual to each other network ---
    for j = 1:numel(networks)
        dist_subj = zeros(nSub,1);
        for s = 1:nSub
            mean_vis = mean(grad_norm(visual, s));
            mean_net = mean(grad_norm(networks{j}, s));
            dist_subj(s) = abs(mean_vis - mean_net);
        end
        dist_all(:,j) = dist_subj;
        % Compare patients vs controls
        [pvals_empirical{g_i}(j),h,stats_empirical{g_i}(j)] = ranksum(dist_subj(1:33), dist_subj(34:83));
    end
    
    figure('Units','normalized','Position',[0.1 0.1 0.55 0.55], 'Color','w');
    hold on;

    nPairs = numel(networks);
    width = 0.6;
    cmap = [0 0 0.8; 0.8 0 0];

    all_data = [];
    all_group = [];
    all_pos = [];

    for j = 1:nPairs

        % Data
        data_pat  = dist_all(1:33,j);
        data_ctrl = dist_all(34:83,j);

        % Box positions (shifted globally)
        pos_pat  = 2*(j-1) + 1;
        pos_ctrl = 2*(j-1) + 2;

        % Append to global arrays
        all_data  = [all_data; data_pat; data_ctrl];
        all_group = [all_group; ...
                     ones(numel(data_pat),1)*pos_pat; ...
                     ones(numel(data_ctrl),1)*pos_ctrl];
    end

    % --- Boxplots ---
    boxplot(all_data, all_group, ...
        'Positions', 1:(2*nPairs), ...
        'Widths', width, ...
        'Colors', repmat(cmap, nPairs, 1), ...
        'Symbol','');
    
    % Make boxplots more visible
    h = findobj(gca, 'Tag', 'Box');
    for k = 1:numel(h)
        patch(get(h(k),'XData'), get(h(k),'YData'), ...
            get(h(k),'Color'), ...
            'FaceAlpha', 0.20, ...
            'EdgeColor', get(h(k),'Color'), ...
            'LineWidth', 2);
    end
    
    % Thicker median lines
    set(findobj(gca,'Tag','Median'), 'LineWidth', 2.5);
    % Thicker whiskers and caps
    set(findobj(gca,'Tag','Whisker'), 'LineWidth', 1.8);
    set(findobj(gca,'Tag','Cap'), 'LineWidth', 1.8);
    
    % --- Scatter points ---
    for j = 1:nPairs
    
        data_pat  = dist_all(1:33,j);
        data_ctrl = dist_all(34:83,j);
    
        pos_pat  = 2*(j-1) + 1;
        pos_ctrl = 2*(j-1) + 2;
    
        scatter(pos_pat + randn(size(data_pat))*0.07, data_pat, ...
            24, cmap(1,:), 'filled', ...
            'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor','none');
    
        scatter(pos_ctrl + randn(size(data_ctrl))*0.07, data_ctrl, ...
            24, cmap(2,:), 'filled', ...
            'MarkerFaceAlpha', 0.25, 'MarkerEdgeColor','none');
    end

    % --- Legend with coloured squares ---
    h_pat = patch(nan, nan, cmap(1,:), ...
        'FaceAlpha', 0.20, ...
        'EdgeColor', cmap(1,:), ...
        'LineWidth', 2);
    
    h_ctrl = patch(nan, nan, cmap(2,:), ...
        'FaceAlpha', 0.20, ...
        'EdgeColor', cmap(2,:), ...
        'LineWidth', 2);
    
    lgd = legend([h_pat h_ctrl], ...
        {'Early blind', 'Sighted controls'}, ...
        'Location', 'best', ...
        'FontSize', 20, ...
        'Box', 'off');

    % --- X-axis formatting ---
    xticks(1.5:2:2*nPairs-0.5);
    xticklabels(net_names);
    
    ylabel(['Visual–Network distance (Gradient ' num2str(g_i) ')'], ...
        'FontSize', 26, 'FontWeight', 'bold'); 
    title('Visual cortex segregation vs other networks', ...
        'FontSize', 27, 'FontWeight', 'bold');  
    set(gca, ...
        'LineWidth', 2, ...
        'FontSize', 24, ...          % Tick numbers and network labels
        'FontWeight', 'bold', ...    % Make tick labels more visible
        'Box', 'off', ...
        'TickDir', 'out');
    
    % Improve tick appearance
    ax = gca;
    ax.TickLength = [0.015 0.015];
    ax.XAxis.FontSize = 23;
    ax.YAxis.FontSize = 23;

    saveas(gca,fullfile(outputDir, label, ['Boxplot_Gradient' num2str(g_i) '.jpg']))

end




%% Plot the variance explained
lambda_pat = Galign_patients.lambda{1};
lambda_ctrl = Galign_controls.lambda{1};

% Variance explained (%)
varExp_pat = lambda_pat ./ sum(lambda_pat) * 100;
varExp_ctrl = lambda_ctrl ./ sum(lambda_ctrl) * 100;
nGradToPlot = min(9, numel(lambda_pat));  % safety

cmap = [0 0 0.8; 0.8 0 0];

figure('Units','normalized','Position',[0.1 0.1 0.55 0.55], 'Color','w'); 
hold on;

plot(1:nGradToPlot, varExp_pat(1:nGradToPlot), '-^', ...
    'LineWidth',6,'MarkerSize',10, 'Color',cmap(1,:));

plot(1:nGradToPlot, varExp_ctrl(1:nGradToPlot), '-o', ...
    'LineWidth',6,'MarkerSize',10, 'Color',cmap(2,:));

xlabel('Gradient index','FontSize',28);
ylabel('Embedding Contribution','FontSize',28);
legend({'EB','SC'},'Location','northeast');
set(gca, 'LineWidth', 2, 'FontSize', 28, 'Box','off', 'TickDir','out');
saveas(gca,fullfile(outputDir, label, 'VarianceExplained.jpg'))


% Correlation between reference and my data
Gref1 = Gref_parcel(:,1);
Gref2 = Gref_parcel(:,2);
Gmy1  = Galign_patients.aligned{1}(:,1);
Gmy2  = Galign_patients.aligned{1}(:,2);

% Combine into matrices (ROIs x gradients)
GrefMat = [Gref1, Gref2];   % 100 x 2
GmyMat  = [Gmy1,  Gmy2];    % 100 x 2
% Compute 2x2 correlation matrix across ROIs
corrMat = corr(GmyMat, GrefMat);
disp('Correlation between my gradients and reference gradients:');
disp(corrMat);



%% FDR on all pvalues - Distances between visual network and other networks

all_p = [pval_range; pvals_empirical{1}(:); pvals_empirical{2}(:)];

p_fdr_all = mafdr(all_p(:), 'BHFDR', true);


for ii = 1:6

    r(ii,1) = abs(stats_empirical{1,1}(ii).zval / sqrt(N));

end

close all;