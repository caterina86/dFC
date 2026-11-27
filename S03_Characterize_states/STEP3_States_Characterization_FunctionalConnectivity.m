% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Extract functional connectome from each state and test for significant
% differences across states

% Dependencies:
% HMM-MAR
% osl needs to be downloaded in the github folder -> https://ohba-analysis.github.io/osl-docs/

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


%% Estimate the FC for each state

for num_states_i = 1:length(K)
    
    state_i = K(num_states_i);
    k = state_i-1; 
    
    resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca), ['zeromean' num2str(zeromean)]);
    load(fullfile(resultsDir, 'HMMs_GROUP.mat'));
    N = length(subjects); % no. subjects

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
    mkdir(outputDir);


    %% Define the nodes of each resting state network

    visual_network = [(1:9)'; (51:58)'];
    somot_network = [(10:15)'; (59:66)'];
    dorsatt_network = [(16:23)'; (67:73)'];
    salatt_network = [(24:30)'; (74:78)'];
    limbic_network = [(31:33)'; (79:80)'];
    control_network = [(34:37)'; (81:89)'];
    default_network = [(38:50)'; (90:100)'];



    %% Extract the Functional connectivity matrices when a state is active

    % Differences in FC between states
    for k_i = 1:(k+1)

        [~, correlation_matrix] = getFuncConn(hmm{rep_i,k},k_i);
        save(fullfile(outputDir, ['FC_state' num2str(k_i) '.mat']), 'correlation_matrix')

        correlation_matrix(logical(eye(size(correlation_matrix)))) = 0;
        omega_z = atanh(correlation_matrix); % from r to z-fisher

        thr_correlations(:,:,k_i) = omega_z;

        figure; imagesc(omega_z); colorbar; caxis([-1 1])
        title(['All Subjects State ' num2str(k_i)], 'FontSize', 18)
        saveas(gca, fullfile(outputDir, ['FC_state' num2str(k_i) '.bmp'])) 

    end

    close all;

    for i = 1:(k+1)
        for j = i+1:(k+1)       
            % Get the size of the matrix
            [m, n] = size(thr_correlations(:,:,i));
            upper_triangle = triu(true(m, n), 1);

            x = thr_correlations(:,:,i);
            y = thr_correlations(:,:,j);
            x = x(upper_triangle);
            y = y(upper_triangle);

            % Perform Wilcoxon ranksum test
            [pvalue(i,j), h] = signrank(x, y);

        end
    end

    % Flatten p-values matrix
    p_flattened = pvalue(triu(true(size(pvalue)), 1));

    % Apply FDR correction only to non-zero p-values
    adjustedPValues = mafdr(p_flattened, 'BHFDR', true);

    % Reconstruct the matrix structure for the adjustedPValues
    adjustedPMatrix = zeros(size(pvalue));
    adjustedPMatrix(triu(true(size(pvalue)), 1)) = adjustedPValues;

    % Define network labels
    network_labels = cell(100, 1); % Initialize cell array for labels
    network_labels(visual_network) = repmat({'visual'}, length(visual_network), 1); % Set labels for visual network nodes
    network_labels(somot_network) = repmat({'somot'}, length(somot_network), 1); % Set labels for somot network nodes
    network_labels(dorsatt_network) = repmat({'dorsatt'}, length(dorsatt_network), 1); % Set labels for dorsatt network nodes
    network_labels(salatt_network) = repmat({'salatt'}, length(salatt_network), 1); % Set labels for salatt network nodes
    network_labels(limbic_network) = repmat({'limbic'}, length(limbic_network), 1); % Set labels for limbic network nodes
    network_labels(control_network) = repmat({'control'}, length(control_network), 1); % Set labels for control network nodes
    network_labels(default_network) = repmat({'def'}, length(default_network), 1); % Set labels for default network nodes

    close all;
    
end



