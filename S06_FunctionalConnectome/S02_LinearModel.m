
% Do a linear model of the weighted FC of all states, using the empirical
% stationary FC as dependent variable, and extract the variance explained

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

sizes = {'Extent'}; %{'Intensity'}; %{'Extent'};

outputStatDir = fullfile(outputDir, 'NBS', sizes{1}, 'staticFC');
Nsub = 83;

% load the design matrix to add covariates of no interest
load(fullfile(outputDir, 'NBS/designMatrix.mat'))
X_cov = [Y(:,3), Y(:,4), Y(:,5)];




%% Linear model

for ii = 1

    if ii == 1
        % load the model based stationary FC
        load(fullfile(outputStatDir, 'matrices_TimeAveragedFC.mat'));
        input = matrices_static;
    else
        % load the empirical stationary FC
        load(fullfile(outputStatDir, 'matrices_staticFC.mat'));
        input = static_FCz;
    end

    [m, n, Sub] = size(input);  
    triangle = triu(true(m, n), 1);
    E = sum(triangle(:)); % number of edges   
    % Preallocate
    stat_FC_triangle = zeros(Sub, E);   
    % Extract upper triangle for each subject
    for s = 1:Sub
        tmp = input(:,:,s);    
        stat_FC_triangle(s,:) = tmp(triangle);
    end

    % -----------------------------------------------------------
    % --- NEW: BUILD dyn_FC_triangle FROM FO-WEIGHTED FC --------
    % -----------------------------------------------------------
    
    K = k+1;  % number of HMM states
    dyn_FC_triangle = zeros(Sub, E, K);
    dyn_FCNoW_triangle = zeros(Sub, E, K);
    
    % Compute FO-weighted FC for every subject × state
    FCweighted = zeros(m, n, K, Sub);
    FCweighted_z = zeros(m, n, K, Sub);
    FC_NoWeighted = zeros(m, n, K, Sub);
    FC_NoWeighted_z = zeros(m, n, K, Sub);
    
    for sub_i = 1:Sub
        for s_i = 1:K
            [~, FC_state] = getFuncConn(HMMs_dualregr{k, sub_i}, s_i);
            % FO-weighted connectivity weighted by FO
            FCweighted(:,:,s_i,sub_i) = FC_state * FOdual{rep_i,k}(sub_i,s_i);
            FC_NoWeighted(:,:,s_i,sub_i) = FC_state;
            
            % Fisher z-transform
            FCweighted_z(:,:,s_i,sub_i) = atanh(FCweighted(:,:,s_i,sub_i)); % size: 100 x 100 x 7 x 83
            FC_NoWeighted_z(:,:,s_i,sub_i) = atanh(FC_NoWeighted(:,:,s_i,sub_i)); % size: 100 x 100 x 7 x 83
        end
    end
    
    % Convert dynamic FC into edge format
    for sub_i = 1:Sub
        for s_i = 1:K
            tmp = FCweighted_z(:,:,s_i,sub_i);
            dyn_FC_triangle(sub_i,:,s_i) = tmp(triangle); % 83 x 4950 x 7
            tmp2 = FC_NoWeighted_z(:,:,s_i,sub_i);
            dyn_FCNoW_triangle(sub_i,:,s_i) = tmp2(triangle); % 83 x 4950 x 7
        end
    end
    
    E = size(stat_FC_triangle,2);
    STATES = size(dyn_FC_triangle,3);
    n_cov = size(X_cov,2);         
    Sub = size(stat_FC_triangle,1);
            
    % Preallocate
    y_hat_edge_W = zeros(Sub, E);   % predicted values for full model
    y_hat_edge_NW = zeros(Sub, E);    % predicted values for covariates-only
    R2_adj_Weighted = zeros(E,1);
    R2_adj_NoWeighted = zeros(E,1);
    Delta_R2 = zeros(E,1);
    
    for e = 1:E
        y = stat_FC_triangle(:, e);                  % Sub x 1
        x_dyn = squeeze(dyn_FC_triangle(:, e, :));   % Sub x K
        x_dyn2 = squeeze(dyn_FCNoW_triangle(:, e, :));   % Sub x K

        % ---- Covariates-only model ----
        X_cov_edge = [ones(Sub,1), X_cov, x_dyn];          % Sub x (1 + n_cov)
        b_cov = X_cov_edge \ y;
        y_hat_cov = X_cov_edge * b_cov;
        y_hat_edge_W(:, e) = y_hat_cov;
    
        SS_tot = sum((y - mean(y)).^2);
        SS_res_cov = sum((y - y_hat_cov).^2);
        R2_cov = 1 - SS_res_cov / SS_tot;
        k_cov = size(X_cov_edge,2) - 1;
        R2_adj_Weighted(e) = 1 - (1 - R2_cov) * (Sub - 1) / (Sub - k_cov - 1);
    
        % ---- Covariates + dynamic FC model ----
        X_full_edge = [ones(Sub,1), X_cov, x_dyn2];  % Sub x (1 + n_cov + K)
        b_full = X_full_edge \ y;
        y_hat_full = X_full_edge * b_full;
        y_hat_edge_NW(:, e) = y_hat_full;
    
        SS_res_full = sum((y - y_hat_full).^2);
        R2_full = 1 - SS_res_full / SS_tot;
        k_full = size(X_full_edge,2) - 1;
        R2_adj_NoWeighted(e) = 1 - (1 - R2_full) * (Sub - 1) / (Sub - k_full - 1);
    
        % ---- Delta R2 ----
        Delta_R2(e) = R2_adj_Weighted(e) - R2_adj_NoWeighted(e);
    end
    
    % ---- Optional summaries ----
    median_R2_Weighted = median(R2_adj_Weighted);
    median_R2_NoWeighted = median(R2_adj_NoWeighted);
    median_DeltaR2 = median(Delta_R2);

    m = 100; % number of parcels (same as in original FC)
    triangle = triu(true(m,m),1);  % logical mask for upper triangle
    tmp_weighted = zeros(m,m);
    tmp_noweighted = zeros(m,m);
    tmp_delta = zeros(m,m);
    
    % Map edge-wise R2 values to full matrix
    tmp_weighted(triangle) = R2_adj_Weighted;        % 1 x E vector
    tmp_noweighted(triangle) = R2_adj_NoWeighted;
    tmp_delta(triangle) = Delta_R2;
    
    % Mirror upper -> lower triangle
    tmp_weighted = tmp_weighted + tmp_weighted';
    tmp_noweighted = tmp_noweighted + tmp_noweighted';
    tmp_delta = tmp_delta + tmp_delta';
    
    % Store in output (overwrite or create new variable)
    R2_W = tmp_weighted;
    R2_NW = tmp_noweighted;
    Delta_R2_full = tmp_delta;
    
    % Choose label for plotting
    if ii == 1
        x_label = 'Model based FC';
    else
        x_label = 'Empirical FC';
    end
    
    % ---- Histograms ----
    figure('Color', 'w');  % white background
    h = histogram(R2_adj_Weighted, ...
        'FaceColor', [0.3 0.3 0.9], ...     % nicer blue
        'EdgeColor', 'none', ...
        'BinWidth', 0.02);                  % adjust as needed 
    hold on;  
    % Median line
    med_val = median(R2_adj_Weighted);
    xline(med_val, 'r--', 'LineWidth', 2);
    % Labels and fonts
    xlabel('Adjusted R²', 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Count', 'FontSize', 20, 'FontWeight', 'bold');
    title('Distribution of Adjusted R² Values', 'FontSize', 22, 'FontWeight', 'bold');    
    % Ticks and axes styling
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');
    if ii == 1
        saveas(gca, fullfile(outputDir, 'NBS/Extent/Figures/AdjR2_Weigheted_ModelBased.jpg'));
    else
        saveas(gca, fullfile(outputDir, 'NBS/Extent/Figures/AdjR2_Weigheted_Empirical.jpg'));
    end
    
    % figure; 
    % histogram(R2_adj_NoWeighted);
    % xline(median(R2_adj_NoWeighted), 'r--', 'LineWidth', 1.5)
    % xlabel('Adjusted R2 (full model)', 'FontSize', 18)
    % if ii == 1
    %     saveas(gca, fullfile(outputDir, 'NBS/Extent/Figures/AdjR2_NoWeigheted_ModelBased.jpg'));
    % else
    %     saveas(gca, fullfile(outputDir, 'NBS/Extent/Figures/AdjR2_NoWeigheted_Empirical.jpg'));
    % end
    % 
    % figure; 
    % histogram(Delta_R2);
    % xline(median(Delta_R2), 'r--', 'LineWidth', 1.5)
    % xlabel('Delta Adjusted R2 (full - covariates)', 'FontSize', 18)
    % if ii == 1
    %     saveas(gca, fullfile(outputDir, 'NBS/Extent/Figures/DeltaR2_ModelBased.jpg'));
    % else
    %     saveas(gca, fullfile(outputDir, 'NBS/Extent/Figures/DeltaR2_Empirical.jpg'));
    % end
    
    % ---- Scatterplot (optional: full model predictions) ----
    figure('Color', 'w'); 
    x = stat_FC_triangle(:);
    y = y_hat_edge_W(:);
    
    % Scatter with slight transparency (better for dense point clouds)
    scatter(x, y, 10, 'MarkerFaceColor', [0.3 0.3 0.9], ...
                       'MarkerEdgeColor', 'none', ...
                       'MarkerFaceAlpha', 0.4);
    hold on;
    % Identity line
    lims = [min([x; y])  max([x; y])];
    plot(lims, lims, 'r--', 'LineWidth', 2);   
    % Labels & style
    xlabel(x_label, 'FontSize', 20, 'FontWeight', 'bold');
    ylabel('Predicted FC', 'FontSize', 20, 'FontWeight', 'bold');   
    axis equal; xlim(lims); ylim(lims);   
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out');    
    title('Real vs Predicted Static FC', 'FontSize', 22, 'FontWeight', 'bold');   
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out'); 
    hold off;
   
    if ii == 1
        saveas(gca, fullfile(outputDir, 'NBS/Extent/Figures/Scatterplot_WeightedFC_Model.jpg'))
    else
        saveas(gca, fullfile(outputDir, 'NBS/Extent/Figures/Scatterplot_WeightedFC_Empirical.jpg'))
    end
    
    
    % ---- Median edge-wise correlation ----
    r_edges = zeros(E,1);
    for e = 1:E
        r_edges(e) = corr(stat_FC_triangle(:,e), y_hat_edge_W(:,e));
    end
    median_r = median(r_edges);
    
    disp(['Overall median correlation (full model): r = ' num2str(median_r)]);
    disp(['Median Adjusted R2 (weighted FC): ' num2str(median(R2_adj_Weighted ))]);
    disp(['Median Adjusted R2 (No-weighted FC): ' num2str(median(R2_adj_NoWeighted ))]);
    disp(['Median Delta R2: ' num2str(median(Delta_R2))]);


    

    %% Permutations to test whether the variance explained is different than chance
    nPerm = 1000;
    R2_perm = zeros(E, nPerm);
    
    disp('Permutations')

    for perm = 1:nPerm
        perm_idx = randperm(Sub);
        for e = 1:E
            y = stat_FC_triangle(:, e);
            x_dyn= squeeze(dyn_FC_triangle(:, e, :));
            x_dyn_perm = (x_dyn(perm_idx,:));  
            X_edge_perm = [ones(Sub,1), X_cov, x_dyn_perm];
            b_perm = X_edge_perm \ y;
            y_hat_perm = X_edge_perm * b_perm;
    
            SS_res_perm = sum((y - y_hat_perm).^2);
            SS_tot = sum((y - mean(y)).^2);
            R2 = 1 - SS_res_perm / SS_tot;
            k = size(X_edge_perm,2) - 1;
            R2_perm(e,perm) = 1 - (1 - R2) * (Sub - 1) / (Sub - k - 1);
    
        end
    end
    
    % --- COMPUTE P-VALUES ---
    p_vals = zeros(E,1);
    for e = 1:E
        p_vals(e) = mean(R2_perm(e,:) >=  R2_adj_Weighted(e));
    end
    
    % Optional: FDR correction
    [p_fdr] = mafdr(p_vals,'BHFDR',true);
    
    % --- SIGNIFICANT EDGES ---
    sig_edges = p_fdr < 0.05;
    perc_sign_edges = (sum(sig_edges)/length(p_fdr))*100;

    disp(['Number of edges significant at FDR<0.05: ' num2str(sum(sig_edges)) ' for ' x_label]);
    disp(['Percentage of edges significant at FDR<0.05: ' num2str(sum(perc_sign_edges)) ' for ' x_label]);
    
    % Map significant edges to full matrix
    sig_matrix = zeros(m,m);
    sig_matrix(triangle) = sig_edges;
    sig_matrix = sig_matrix + sig_matrix';


end


% % Same type of analysis, but using the function lm in matlab
% for e = 1:E
%     % ----------------------------------
%     % Prepare data for this edge
%     % ----------------------------------
%     y = stat_FC_triangle(:,e);               % static FC (Sub × 1)
%     x = squeeze(dyn_FC_triangle(:,e,:));     % dynamic FC (Sub × K)
% 
%     % Convert to table for fitlm
%     T = array2table([y x X_cov], ...
%         'VariableNames', ['y', strcat('x', string(1:size(x,2))), 'x8', 'x9', 'x10']);
% 
%     % ----------------------------------
%     % Fit linear model
%     % ----------------------------------
%     formula = "y ~ 1 + " + strjoin("x" + string(1:K+3), " + ");
%     lm = fitlm(T, formula);
% 
%     % ----------------------------------
%     % Extract predictions
%     % ----------------------------------
%     y_hat = predict(lm, T);
% 
%     % store predicted values like before
%     y_hat_all((1:Sub) + (e-1)*Sub) = y_hat;
% 
%     % ----------------------------------
%     % Extract R2, adjusted R2, correlation
%     % ----------------------------------
%     R2(e)      = lm.Rsquared.Ordinary;
%     R2_adj(e)  = lm.Rsquared.Adjusted;
%     r_pred_real(e) = corr(y, y_hat);
% end
