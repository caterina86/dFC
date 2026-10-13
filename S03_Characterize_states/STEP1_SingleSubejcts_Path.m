% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Plot single-subjects trend for each model order

% Dependencies:
% HMM-MAR

clearvars; close all;

% load paths
baseDir = uigetdir(pwd, 'Select the base directory'); % select the directory where of the project
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
  
% check the parameters of the HMM
if options.pca ~= 10
    disp(['PCA is different from ' num2str(pca)]);
end
if options.zeromean ~= 0
    disp(['zeromean is different from ' num2str(zeromean)]);
end
if options.K ~= 15
    disp(['max num states is different from ' num2str(max(K))]);
end



%% Choose the repetition with the lowest free energy, for each model order

for num_states_i = 1:length(K)
    
    state_i = K(num_states_i);
    k = state_i-1; 
    
    % Plot FE and Average log-likelihood
    figure;
    for repet_i = 1:repetitions
        hold on;
        plot(K,cell2mat(fe(repet_i,:)))
        xticks(K)
        xlabel('states', 'FontSize', 18); ylabel('free energy', 'FontSize', 18)
        mean_fe(repet_i,1) = mean(cell2mat(fe(repet_i)),1);
        disp(['mean free energy repetition ' num2str(repet_i) ' = ' num2str(mean_fe(repet_i))])
    end
    saveas(gca, fullfile(resultsDir, 'PlotFE.bmp'))

    % Find the minimum value and its indices, based on the Free Energy
    [minValue_fe, minIndex_fe] = min(cell2mat(fe(:,:)));
    [rowIndex_fe, colIndex_fe] = ind2sub(size(fe), minIndex_fe);

    figure;
    for repet_i = 1:repetitions
        hold on;
        plot(K,abs(cell2mat(ll_average(repet_i,:))))
        xticks(K)
        xlabel('states', 'FontSize', 18); ylabel('absolute avLL', 'FontSize', 18)
        mean_ll(repet_i,1) = mean(cell2mat(ll_average(repet_i)),1);
        disp(['mean log-likelihood repetition ' num2str(repet_i) ' = ' num2str(mean_ll(repet_i))])
    end
    saveas(gca, fullfile(resultsDir, 'PlotavLL.bmp'))

    % Find the minimum value and its indices, based on the Free Energy
    [minValue_ll, minIndex_ll] = min(cell2mat(ll_average(:,:)));
    [rowIndex_ll, colIndex_ll] = ind2sub(size(ll_average), minIndex_ll);

    rep_i = rowIndex_fe(1,k); % repetition chosen based on the lowest free energy
    disp(['repetition chosen ' num2str(rep_i)])
    disp(['number states chosen ' num2str(K(k))])


    % Plot single-subjects dynamics
    load(fullfile(resultsDir, ['HMMs_r' num2str(rep_i) '_SUBJECTS.mat']));
    mkdir(fullfile(resultsDir, 'figures/'));
    outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)]);
    mkdir(outputDir);

    for sub_i = 1:N

        % Time courses of the states probabilities (concatenated for all subjects)
        figure; 
        plot_Gamma(Gammaj{k,sub_i}(:,1:k),T,0);
        set(gca,'ylim',[-0.2 1.2]); ylabel('state #');
        set(gca,'xlim',[0 T]); ylabel('state #');
        saveas(gca, fullfile(outputDir, ['sub-' subjects{sub_i} '_TimeCourseStates.bmp']))

        % transition probability matrix of each subject, without the diagonal values  
        figure; 
        tpm_subjects{sub_i} = getTransProbs(HMMs_dualregr{k,sub_i}); % transition probability without the diagonal
        imagesc(tpm_subjects{sub_i}); colormap(jet);
        caxis([min(min(tpm_subjects{sub_i})) max(max(tpm_subjects{sub_i}))]);
        colorbar;
        title(['Transition Probability ' subjects{sub_i}]);
        saveas(gca, fullfile(outputDir, ['sub-' subjects{sub_i} '_TP.bmp']))

        close all;

    end

    % Viterbi Path for each subject
    vector = vpath{rep_i,k};
    % Calculate the number of segments
    numSegments = numel(vector) / T;
    sub_i=1;

    % Plot the viterbi path for each subject
    for i = 1:numSegments
        figure;
        startIndex = (i - 1) * T + 1;
        endIndex = i * T;
        % Adjust x-axis values
        xValues = startIndex:endIndex;
        xValues = xValues - min(xValues) + 1; % Set the minimum x-value to 1

        plot(xValues, vector(startIndex:endIndex), 'LineWidth', 2), set(gca,'Title',text('String','True state path'))
        set(gca,'ylim',[0 hmm{rep_i,k}.K+1]); ylabel('state #')
        title(['Viterbi Path subject ' subjects{sub_i}]);   
        xlim([0 T]);
        saveas(gca, fullfile(outputDir, ['sub-' subjects{sub_i} '_ViterbiPath.bmp']))
        sub_i=sub_i + 1;

    end

    % Plot Fractional occupancy of the group
    correlation = corr(FOgroup{rep_i,k});
    % Perform hierarchical clustering using the 'average' linkage method
    pairwise_distances = pdist(correlation);
    tree = linkage(pairwise_distances, 'ward');

    % Visualize the dendrogram
    figure; 
    subplot(1, 2, 1); % 1 row, 2 columns, first subplot
    dendrogram(tree);
    title('Dendrogram', 'FontSize', 18);
    xlabel('State', 'FontSize', 18);
    ylabel('Distance', 'FontSize', 18);
    % Plot the crrelation of FO
    subplot(1, 2, 2); % 1 row, 2 columns, first subplot
    correlation(eye(size(correlation))==1) = nan;
    imagesc(correlation); colorbar; colormap(jet); caxis([-1 1])
    title('Fractional Occupancy', 'FontSize', 18);    
    pbaspect([size(correlation, 2) size(correlation, 1) 1]);
    saveas(gca, fullfile(outputDir, 'Group_FO.bmp'))

    % Plot Transition Probability of the group
    out = getTransProbs(hmm{rep_i, k});
    pairwise_distances = pdist(out);
    tree = linkage(pairwise_distances, 'ward');

    figure;
    subplot(1, 2, 1); % 1 row, 2 columns, first subplot
    dendrogram(tree);
    title('Dendrogram', 'FontSize', 18);
    xlabel('State', 'FontSize', 18);
    ylabel('Distance', 'FontSize', 18);
    subplot(1, 2, 2); % 1 row, 2 columns, first subplot
    imagesc(out); colormap(jet);
    caxis([min(min(out)) max(max(out))]);
    colorbar;
    title('Transition Probability', 'FontSize', 18);
    pbaspect([size(out, 2) size(out, 1) 1]);
    saveas(gca, fullfile(outputDir, 'Group_TP.bmp'))

    G = digraph(out);
    % Remove edges with NaN weights
    nanEdges = isnan(G.Edges.Weight);
    G = rmedge(G, find(nanEdges));
    % Plot the directed graph
    base_line_width = 1;
    suggested_scaling_factor = 150;

    % Plot the directed graph with adjusted arrow thickness and larger arrows
    figure;
    h = plot(G, 'LineWidth', base_line_width + suggested_scaling_factor * (G.Edges.Weight - min(G.Edges.Weight)), 'MarkerSize', 8, 'ArrowSize', 15);
    labelNode = arrayfun(@(x) {num2str(x)}, 1:numnodes(G))';
    h.NodeLabel = labelNode;
    h.NodeFontSize = 18;  

    close all;

end

