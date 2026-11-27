% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Extract mean activation and plot it as a function of parcel (multiple
% visualization of the results)

% Dependencies:
% HMM-MAR
% osl needs to be downloaded in the github folder -> https://ohba-analysis.github.io/osl-docs/

clearvars; close all;

% Get the FSL directory path from the environment variable
fslPath = getenv('FSLDIR');
% Check if FSLDIR is set
if isempty(fslPath)
    error('FSLDIR environment variable is not set. Please set it before running this script.');
end


% load paths
baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
DirData = fullfile(projectDir, 'derivatives/HMM');

mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR
oslDir = fullfile(mydir, 'osl/osl-core');
addpath(genpath(fullfile(projectDir, 'code')))

% define parameters
zeromean = 0;
K = (2:15); % range of the number of states tested
pca = 10;
repetitions = 20;
T = 464; % number of timepoints
num_networks = 7;



%% Extract the maps of mean activation from each model order

for num_states_i = 1:length(K)
    
    state_i = K(num_states_i);
    k = state_i-1; 
    
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

    outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)]);
    mkdir(outputDir);


    %% Plot the mean activations for each state and model

    cd(oslDir);
    osl_startup;

    visual_network = [(1:9)'; (51:58)'];
    somot_network = [(10:15)'; (59:66)'];
    dorsatt_network = [(16:23)'; (67:73)'];
    salatt_network = [(24:30)'; (74:78)'];
    limbic_network = [(31:33)'; (79:80)'];
    control_network = [(34:37)'; (81:89)'];
    default_network = [(38:50)'; (90:100)'];

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
        saveas(gca, fullfile(outputDir, 'MeanActivation.bmp'))
    end

    save(fullfile(outputDir, 'MeanActivation.mat'), 'net_mean')


    % extract the mean activation of different resting state networks
    mean_networks(1,:) = mean((net_mean(visual_network,:)),1); % visual
    mean_networks(2,:) = mean((net_mean(somot_network,:)),1); % somot
    mean_networks(3,:) = mean((net_mean(dorsatt_network,:)),1); % dorsatt
    mean_networks(4,:) = mean((net_mean(salatt_network,:)),1); % salatt
    mean_networks(5,:) = mean((net_mean(limbic_network,:)),1); % limbic
    mean_networks(6,:) = mean((net_mean(control_network,:)),1); % control
    mean_networks(7,:) = mean((net_mean(default_network,:)),1); % default

    figure;
    for ii = 1:num_networks
        hold on;
        scatter(1:(k+1), (mean_networks(ii,:)), 'filled')
    end
    xlim([0 (k+2)]); xlabel('states', 'FontSize', 20); xticks(1:1:(k+1))
    ylabel('Mean BOLD', 'FontSize', 20)
    yline(0, 'r')
    legend('visual', 'somot', 'dorsatt', 'salatt', 'limbic', 'control', ...
        'defaul', 'Location','Best')
    saveas(gca, fullfile(outputDir, 'Plot_MeanActivation_networks.bmp'))


    % Calculate separatelly the mean for left and right hemisphere
    visual_network_lh = (1:9)';
    visual_network_rh = (51:58)';
    somot_network_lh = (10:15)';
    somot_network_rh = (59:66)';
    dorsatt_network_lh = (16:23)';
    dorsatt_network_rh = (67:73)';
    salatt_network_lh = (24:30)';
    salatt_network_rh = (74:78)';
    limbic_network_lh = (31:33)';
    limbic_network_rh = (79:80)';
    control_network_lh = (34:37)';
    control_network_rh = (81:89)';
    default_network_lh = (38:50)';
    default_network_rh = (90:100)';

    mean_networks(1,:) = mean((net_mean(visual_network_lh,:)),1); % visual
    mean_networks(2,:) = mean((net_mean(somot_network_lh,:)),1); % somot
    mean_networks(3,:) = mean((net_mean(dorsatt_network_lh,:)),1); % dorsatt
    mean_networks(4,:) = mean((net_mean(salatt_network_lh,:)),1); % salatt
    mean_networks(5,:) = mean((net_mean(limbic_network_lh,:)),1); % limbic
    mean_networks(6,:) = mean((net_mean(control_network_lh,:)),1); % control
    mean_networks(7,:) = mean((net_mean(default_network_lh,:)),1); % default
    mean_networks(8,:) = mean((net_mean(visual_network_rh,:)),1); % visual
    mean_networks(9,:) = mean((net_mean(somot_network_rh,:)),1); % somot
    mean_networks(10,:) = mean((net_mean(dorsatt_network_rh,:)),1); % dorsatt
    mean_networks(11,:) = mean((net_mean(salatt_network_rh,:)),1); % salatt
    mean_networks(12,:) = mean((net_mean(limbic_network_rh,:)),1); % limbic
    mean_networks(13,:) = mean((net_mean(control_network_rh,:)),1); % control
    mean_networks(14,:) = mean((net_mean(default_network_rh,:)),1); % default

    % Define colors for each row
    colors = lines(14); % or any other colormap you prefer

    figure;
    hold on;
    for ii = 1:num_networks
        scatter(1:(k+1), mean_networks(ii,:), 'filled', 'MarkerFaceColor', colors(ii,:), 'Marker', 'o')
    end
    for ii = num_networks+1:(num_networks*2)
        scatter(1:(k+1), mean_networks(ii,:), 'filled', 'MarkerFaceColor', colors(ii,:), 'Marker', 'd')
    end
    hold off;

    xlim([0 (k+2)]); 
    xlabel('states', 'FontSize', 20); 
    xticks(1:1:(k+1))
    ylabel('Mean BOLD', 'FontSize', 20)
    yline(0, 'r')
    legend('visual lh', 'somot lh', 'dorsatt lh', 'salatt lh', 'limbic lh', 'control lh', ...
        'default lh', 'visual rh', 'somot rh', 'dorsatt rh', 'salatt rh', 'limbic rh', 'control rh', ...
        'default rh', 'Location','Best')
    saveas(gca, fullfile(outputDir, 'Plot_MeanActivation_networks_Lateralization.bmp'))


    % Calculate the mean level of the data
    mean_level = mean(mean(mean_networks));
    % Compute the absolute deviation
    deviation = mean_networks - mean_level;

    % Define custom colormap
    n = 64; % Number of colors
    colors_hot = hot(n); % Get hot colormap
    colors_cool = jet(n); % Get cool colormap
    zero_index = round(n / 2); % Index for zero value
    colors_custom = [(colors_cool(1:zero_index-1,:)); [1 1 1]; [1 1 1]; flipud(colors_hot(10:(zero_index-1)+9,:))]; % Custom colormap

    % Plot the data
    figure;
    
    for ii = 1:size(deviation,2)
        
        subplot(size(deviation,2),1,ii)
        imagesc((deviation(:,ii))');
        colormap(colors_custom);

        % Define labels for y-axis
        xlabels = {'visual lh', 'somot lh', 'dorsatt lh', 'salatt lh', 'limbic lh', 'control lh', ...
            'default lh', 'visual rh', 'somot rh', 'dorsatt rh', 'salatt rh', 'limbic rh', 'control rh', ...
            'default rh'};

        xticks(1:numel(xlabels));
        xticklabels(xlabels);
        xtickangle(45); % Rotate x-axis labels to 45 degrees for oblique orientation
        set(gca, 'XTickLabel', xlabels, 'FontSize', 14); % Set font size
        yticks(1:(state_i));
        colorbar; caxis([-0.6 0.6]);
    
    end
    
    saveas(gca, fullfile(outputDir, 'Matrix_MeanActivation.bmp'))



    %% Extract the activation from each parcel and plot it

    for k_i = 1:(k+1)

        % Plot
        figure;
        h = scatter(1:length(net_mean), net_mean(:,k_i), 'filled', 'o');  % Plot dots
        hold on;

        % Draw lines connecting dots to x-axis (xtick)
        for i = 1:length(net_mean)
            line([i, i], [0, net_mean(i,k_i)], 'Color', 'k', 'LineStyle', '-');
        end

        % Plot mean line
        yline(0, 'r');  % Red line for mean
        hold off;

        % Add colored background rectangles for each network
        labels = {'visual network', 'somot network', 'dorsatt network', 'salatt network', 'limbic network', 'control network', 'default network', ...
            'visual network', 'somot network', 'dorsatt network', 'salatt network', 'limbic network', 'control network', 'default network'};
        networks = {visual_network, somot_network, dorsatt_network, salatt_network, limbic_network, control_network, default_network};
        colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];  % Specify colors

        for netIdx = 1:length(networks)
            network_indices = networks{netIdx};
            color = colors(netIdx);

            for i = 1:length(network_indices)
                x = [network_indices(i)-0.5, network_indices(i)+0.5, network_indices(i)+0.5, network_indices(i)-0.5];
                y = [min(h.YData)-0.01, min(h.YData)-0.01, max(h.YData)+0.01, max(h.YData)+0.01];
                patch(x, y, color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

            end
        end

        xticks([5:7:(length(net_mean))]);
        xticklabels(labels);
        xtickangle(45);  % Rotate x-axis labels for better visibility

        % Add legend for colors
        set(h, 'HandleVisibility', 'off');
        xlim([0.5, length(net_mean) + 0.5]);  % Set xlim to align with xticks
        ylim([min(h.YData)-0.01, max(h.YData)+0.01]);  % Set ylim
        xlabel('Networks', 'FontSize', 20, 'FontName', 'helvetica');
        ylabel('Mean BOLD', 'FontSize', 20, 'FontName', 'helvetica');
        title(['Mean Values state ' num2str(k_i)], 'FontSize', 20, 'FontName', 'helvetica');

        saveas(gca, fullfile(outputDir, ['Plot_MeanActivation_parcels_state' num2str(k_i) '.bmp']))

    end

    % Create the map of mean activation for each state, as a nifti file
    cd(resultsDir)
    parcellationfile = fullfile(projectDir, 'derivatives/Atlas/MNI/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz');
    maskfile = '/Users/cp3488/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz';
    outputfile = fullfile(outputDir, 'map_mean_networks');
    onconnectivity = 0; % mean activity
    map_mean = makeMap(net_mean,1:k+1,parcellationfile,maskfile, onconnectivity, [], [], outputfile); % size -> voxels * states

    % split the mean map
    system(['fslsplit ' [outputfile '.nii.gz'] ' ' fullfile(outputDir, 'mean') ' -t']);

    % Visualize the results with OSL:
    % repeat values to end up with values for each parcel, and visualize it
    net_mean(visual_network,:) = repmat(mean_networks(1,:),length(visual_network),1);
    net_mean(somot_network,:) = repmat(mean_networks(2,:),length(somot_network),1);
    net_mean(dorsatt_network,:) = repmat(mean_networks(3,:),length(dorsatt_network),1);
    net_mean(salatt_network,:) = repmat(mean_networks(4,:),length(salatt_network),1);
    net_mean(limbic_network,:) = repmat(mean_networks(5,:),length(limbic_network),1);
    net_mean(control_network,:) = repmat(mean_networks(6,:),length(control_network),1);
    net_mean(default_network,:) = repmat(mean_networks(7,:),length(default_network),1);

    atlasfile = 'Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz';
    p = parcellation(atlasfile); % load the parcellation
    p.osleyes(net_mean); % call osleyes



    %% Calculate the uncertainty of the mean estimates: coefficient of variations

    for k_i = 1:(k+1)

        omega = getFuncConn(hmm{rep_i,k},k_i); % covariance matrix

        % Ensure positive diagonal elements
        if any(diag(omega) <= 0)
            error('Negative values found in the diagonal of omega_k_sample.');
        end

        % Assuming J is your matrix
        J = 1./omega; 

        % Update the formula
        CV_values(:, k_i) = 100 ./ sqrt(diag(J));

        mean_CV_by_state(1,k_i) = mean(CV_values(:,k_i));
        std_CV_by_state(1,k_i) = std(CV_values(:,k_i));

        % Plot
        figure;
        scatter(1:length(CV_values), CV_values(:,k_i), 'filled', 'o');  % Plot dots
        hold on;

        % Draw lines connecting dots to x-axis (xtick)
        for i = 1:length(CV_values)
            line([i, i], [0, CV_values(i,k_i)], 'Color', 'k', 'LineStyle', '-');
        end

        % Plot mean line
        yline(mean_CV_by_state(1,k_i), 'r');  % Red line for mean
        hold off;

        % Add colored background rectangles for each network
        labels = {'visual network', 'somot network', 'dorsatt network', 'salatt network', 'limbic network', 'control network', 'default network', ...
            'visual network', 'somot network', 'dorsatt network', 'salatt network', 'limbic network', 'control network', 'default network'};
        networks = {visual_network, somot_network, dorsatt_network, salatt_network, limbic_network, control_network, default_network};
        colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k'];  % Specify colors

        for netIdx = 1:length(networks)
        network_indices = networks{netIdx};
        color = colors(netIdx);

            for i = 1:length(network_indices)
                x = [network_indices(i)-0.5, network_indices(i)+0.5, network_indices(i)+0.5, network_indices(i)-0.5];
                y = [0, 0, 120, 120];
                patch(x, y, color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

            end

            % Add legend entry with color
            %legend_entries{netIdx} = ['Network ' labels{netIdx}];
        end

        % Customize x-axis ticks and labels -> check it!
        xticks([5:7:(length(net_mean))]);
        xticklabels(labels);
        xtickangle(45);  % Rotate x-axis labels for better visibility
        xlim([0.5, length(CV_values) + 0.5]);  % Set xlim to align with xticks
        ylim([0, 120]);  % Set ylim
        xlabel('Data Points', 'FontSize', 20, 'FontName', 'helvetica');
        ylabel('CV Values', 'FontSize', 20, 'FontName', 'helvetica');
        title(['CV for state ' num2str(k_i)], 'FontSize', 20, 'FontName', 'helvetica');
        saveas(gca, fullfile(outputDir, ['Plot_CV_parcels_state' num2str(k_i) '.bmp']))

    end

    for k_i = 1:(k+1)

        delta_CV(:, k_i) = CV_values(:, k_i) - mean_CV_by_state(1,k_i);

    end

    close all;

    % Save the map of coefficient of variations as a nifti image
    cd(resultsDir)
    maskfile = fullfile(fslPath, 'data', 'standard', 'MNI152_T1_2mm_brain.nii.gz');
    outputfile = fullfile(outputDir, 'map_CV');
    onconnectivity = 0; % mean activity
    maps_CV = makeMap(delta_CV,1:k+1,parcellationfile,maskfile, onconnectivity, [], [], outputfile); % size -> voxels * states

    close all;
    clear mean_networks
    
end

