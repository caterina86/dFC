% Compare significant states

clearvars;

baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
DirData = fullfile(projectDir, 'derivatives/HMM');

mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR
addpath(genpath(fullfile(projectDir, 'code')))


% Define networks
visual_network   = [(1:9)'; (51:58)'];
somot_network    = [(10:15)'; (59:66)'];
dorsatt_network  = [(16:23)'; (67:73)'];
salatt_network   = [(24:30)'; (74:78)'];
limbic_network   = [(31:33)'; (79:80)'];
control_network  = [(34:37)'; (81:89)'];
default_network  = [(38:50)'; (90:100)'];

network_indices = {
    visual_network
    somot_network
    dorsatt_network
    salatt_network
    limbic_network
    control_network
    default_network
};

network_names = {
    'Visual'
    'Somatomotor'
    'Dorsal attention'
    'Salience/ventral attention'
    'Limbic'
    'Control'
    'Default mode'
};

% Assign each parcel to a network
network_id = zeros(100,1);
for n = 1:numel(network_indices)
    network_id(network_indices{n}) = n;
end
% Colors: one per network
network_colors = [
    0.00 0.45 0.70;  % Visual - blue
    0.90 0.40 0.10;  % Somatomotor - orange
    0.00 0.60 0.50;  % Dorsal attention - green
    0.80 0.60 0.00;  % Salience/ventral attention - yellow
    0.80 0.40 0.60;  % Limbic - pink
    0.40 0.40 0.40;  % Control - grey
    0.35 0.35 0.80   % Default mode - purple
];


% Check number of parcels assigned to each network
for n = 1:7
    fprintf('%s: %d parcels\n', ...
        network_names{n}, sum(network_id == n));
end

% Check for duplicate assignments
all_indices = vertcat(network_indices{:});
fprintf('Total assigned entries: %d\n', numel(all_indices));
fprintf('Unique parcels: %d\n', numel(unique(all_indices)));
fprintf('Unassigned parcels: %d\n', sum(network_id == 0));



% Load the data
map_2017_all = load(fullfile(projectDir, 'derivatives/HMM/PCA_REDUCED_7_train_test_dataset2017/10/zeromean0/figures/NumStates_7/Repetition_5/MeanActivation.mat'));
map_pooled = load(fullfile(projectDir, 'derivatives/HMM/PCA_REDUCED_15/10/zeromean0/figures/NumStates_7/Repetition_2_original/MeanActivation.mat'));
map_2023_all = load(fullfile(projectDir, 'derivatives/HMM/PCA_REDUCED_7_train_test_dataset2023/10/zeromean0/figures/NumStates_7/Repetition_15/MeanActivation.mat'));




%% Montreal - State 7

for ii = 1:7

    % Ensure column vectors
    map_pooled_1 = map_pooled.net_mean(:,ii);
    map_2023 = map_2023_all.net_mean(:,ii);
    map_2017   = map_2017_all.net_mean(:,7);
    
    % Signed spatial correlation across the 100 parcels
    [r_pooled(ii,1), p_pooled(ii,1)] = corr(map_pooled_1, map_2017, 'Type', 'Pearson');
    fprintf('Signed correlation: r = %.3f, p = %.4g\n', r_pooled(ii,1), p_pooled(ii,1));
    [r_2023(ii,1), p_2023(ii,1)] = corr(map_2023, map_2017, 'Type', 'Pearson');
    fprintf('Signed correlation: r = %.3f, p = %.4g\n', r_2023(ii,1), p_2023(ii,1));

end

% Replace these with your actual map vectors
map2017 = map_2017_all.net_mean(:,7);
map3    = map_pooled.net_mean(:,3);

% Plot comparisons

maps_to_compare = {map3};
state_names = {'Pooled Dataset - State 3'};

figure;

for s = 1

    %subplot(1,2,s);
    hold on;

    x = maps_to_compare{s}(:);
    y = map2017(:);

    % Store scatter handles for the legend
    h_scatter = gobjects(7,1);

    % Plot parcels by resting-state network
    for n = 1:7

        idx = network_id == n;

        h_scatter(n) = scatter(x(idx), y(idx), 90, ...
            network_colors(n,:), 'filled', ...
            'DisplayName', network_names{n});

    end

    % Linear regression line
    valid = isfinite(x) & isfinite(y);

    b = polyfit(x(valid), y(valid), 1);

    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(b, x_fit);

    % Regression line (visible, but excluded from legend)
    plot(x_fit, y_fit, 'k-', 'LineWidth', 2.5);

    % Pearson correlation
    r = corr(x(valid), y(valid), 'Type', 'Pearson');

    text(0.05, 0.95, sprintf('r = %.2f', r), ...
        'Units', 'normalized', ...
        'FontSize', 22, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w');

    % Labels
    xlabel(char(string(state_names{s})), 'FontSize', 18);
    ylabel('Montreal Dataset - State 7', 'FontSize', 18);

    % Axes
    axis square;
    grid on;

    set(gca, ...
        'FontSize', 14, ...
        'LineWidth', 1);

    % Legend: only the seven networks
    legend(h_scatter, network_names, ...
        'Location', 'best', ...
        'FontSize', 14);

end



map2017 = map_2017_all.net_mean(:,7);
map3    = map_pooled.net_mean(:,3);
map7    = map_pooled.net_mean(:,7);

maps_to_compare = [map2017, map3, map7];

nNetworks = numel(network_indices);
mean_activation = zeros(nNetworks, size(maps_to_compare,2));

for s = 1:size(maps_to_compare,2)
    for n = 1:nNetworks
        mean_activation(n,s) = ...
            mean(maps_to_compare(network_indices{n},s), 'omitnan');
    end
end

% correlation considering the mean value within a network
r_network1 = corr(mean_activation(:,1), mean_activation(:,2), ...
                  'Type', 'Pearson');
r_network2 = corr(mean_activation(:,1), mean_activation(:,3), ...
                  'Type', 'Pearson');

fprintf('2017 vs pooled state 3: r = %.3f\n', r_network1);
fprintf('2017 vs pooled state 7: r = %.3f\n', r_network2);






%% Montreal - State 2
for ii = 1:7

    % Ensure column vectors
    map_pooled_1 = map_pooled.net_mean(:,ii);
    map_2023 = map_2023_all.net_mean(:,ii);
    map_2017   = map_2017_all.net_mean(:,2);
    
    % Signed spatial correlation across the 100 parcels
    [r_pooled(ii,1), p_pooled(ii,1)] = corr(map_pooled_1, map_2017, 'Type', 'Pearson');
    fprintf('Signed correlation: r = %.3f, p = %.4g\n',r_pooled(ii,1), p_pooled(ii,1));
    [r_2023(ii,1), p_2023(ii,1)] = corr(map_2023, map_2017, 'Type', 'Pearson');
    fprintf('Signed correlation: r = %.3f, p = %.4g\n', r_2023(ii,1), p_2023(ii,1));

end

% Replace these with your actual map vectors
map2017 = map_2017_all.net_mean(:,2);
map1    = map_pooled.net_mean(:,1);
map7    = map_pooled.net_mean(:,7);
map2    = map_2023_all.net_mean(:,2);

% Assign each parcel to a network
network_id = zeros(100,1);
for n = 1:numel(network_indices)
    network_id(network_indices{n}) = n;
end
% Colors: one per network
network_colors = lines(7);

% Plot comparisons
figure;
maps_to_compare = {map1, map7};
state_names = {'Pooled Dataset - State 1'};

for s = 1

    %subplot(1,2,s);
    hold on;

    x = maps_to_compare{s}(:);
    y = map2017(:);

    % Store scatter handles for the legend
    h_scatter = gobjects(7,1);

    % Plot parcels by resting-state network
    for n = 1:7

        idx = network_id == n;

        h_scatter(n) = scatter(x(idx), y(idx), 90, ...
            network_colors(n,:), 'filled', ...
            'DisplayName', network_names{n});

    end

    % Linear regression line
    valid = isfinite(x) & isfinite(y);

    b = polyfit(x(valid), y(valid), 1);

    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(b, x_fit);

    % Regression line (visible, but excluded from legend)
    plot(x_fit, y_fit, 'k-', 'LineWidth', 2.5);

    % Pearson correlation
    r = corr(x(valid), y(valid), 'Type', 'Pearson');

    text(0.05, 0.95, sprintf('r = %.2f', r), ...
        'Units', 'normalized', ...
        'FontSize', 22, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w');

    % Labels
    xlabel(char(string(state_names{s})), 'FontSize', 18);
    ylabel('Montreal Dataset - State 2', 'FontSize', 18);

    % Axes
    axis square;
    grid on;

    set(gca, ...
        'FontSize', 14, ...
        'LineWidth', 1);

    % Legend: only the seven networks
    legend(h_scatter, network_names, ...
        'Location', 'best', ...
        'FontSize', 14);

end



% Choose between state 1 and 7
map2017 = map_2017_all.net_mean(:,2);
map1    = map_pooled.net_mean(:,1);
map7    = map_pooled.net_mean(:,7);

maps_to_compare = [map2017, map1, map7];

nNetworks = numel(network_indices);
mean_activation = zeros(nNetworks, size(maps_to_compare,2));

for s = 1:size(maps_to_compare,2)
    for n = 1:nNetworks
        mean_activation(n,s) = ...
            mean(maps_to_compare(network_indices{n},s), 'omitnan');
    end
end

% correlation considering the mean value within a network
r_network1 = corr(mean_activation(:,1), mean_activation(:,2), ...
                  'Type', 'Pearson');

r_network7 = corr(mean_activation(:,1), mean_activation(:,3), ...
                  'Type', 'Pearson');

fprintf('2017 vs pooled state 1: r = %.3f\n', r_network1);
fprintf('2017 vs pooled state 7: r = %.3f\n', r_network7);






%% Trento - State 2

for ii = 1:7

    % Ensure column vectors
    map_pooled_1 = map_pooled.net_mean(:,ii);
    map_2017   = map_2017_all.net_mean(:,ii);
    map_2023 = map_2023_all.net_mean(:,2);

    % Signed spatial correlation across the 100 parcels
    [r_pooled(ii,1), p_pooled(ii,1)] = corr(map_pooled_1, map_2023, 'Type', 'Pearson');
    fprintf('Signed correlation: r = %.3f, p = %.4g\n', r_pooled(ii,1), p_pooled(ii,1));
    [r_2017(ii,1), p_2017(ii,1)] = corr(map_2017, map_2023, 'Type', 'Pearson');
    fprintf('Signed correlation: r = %.3f, p = %.4g\n', r_2017(ii,1), p_2017(ii,1));

end

% Replace these with your actual map vectors
map2023 = map_2023_all.net_mean(:,2);
map1    = map_pooled.net_mean(:,1);
map2    = map_2017_all.net_mean(:,2);

% Plot comparisons
figure;
maps_to_compare = {map1};
state_names = {'Pooled Dataset - State 1'};


for s = 1

    %subplot(1,2,s);
    hold on;

    x = maps_to_compare{s}(:);
    y = map2023(:);

    % Store scatter handles for the legend
    h_scatter = gobjects(7,1);

    % Plot parcels by resting-state network
    for n = 1:7

        idx = network_id == n;

        h_scatter(n) = scatter(x(idx), y(idx), 90, ...
            network_colors(n,:), 'filled', ...
            'DisplayName', network_names{n});

    end

    % Linear regression line
    valid = isfinite(x) & isfinite(y);

    b = polyfit(x(valid), y(valid), 1);

    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(b, x_fit);

    % Regression line (visible, but excluded from legend)
    plot(x_fit, y_fit, 'k-', 'LineWidth', 2.5);

    % Pearson correlation
    r = corr(x(valid), y(valid), 'Type', 'Pearson');

    text(0.05, 0.95, sprintf('r = %.2f', r), ...
        'Units', 'normalized', ...
        'FontSize', 22, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w');

    % Labels
    xlabel(char(string(state_names{s})), 'FontSize', 18);
    ylabel('Trento Dataset - State 2', 'FontSize', 18);

    % Axes
    axis square;
    grid on;

    set(gca, ...
        'FontSize', 14, ...
        'LineWidth', 1);

    % Legend: only the seven networks
    legend(h_scatter, network_names, ...
        'Location', 'best', ...
        'FontSize', 14);

end



%% Trento - State 6
for ii = 1:7

    % Ensure column vectors
    map_pooled_1 = map_pooled.net_mean(:,ii);
    map_2017   = map_2017_all.net_mean(:,ii);
    map_2023 = map_2023_all.net_mean(:,6);

    % Signed spatial correlation across the 100 parcels
    [r_pooled(ii,1), p_pooled(ii,1)] = corr(map_pooled_1, map_2023, 'Type', 'Pearson');
    fprintf('Signed correlation: r = %.3f, p = %.4g\n', r_pooled(ii,1), p_pooled(ii,1));
    [r_2017(ii,1), p_2017(ii,1)] = corr(map_2017, map_2023, 'Type', 'Pearson');
    fprintf('Signed correlation: r = %.3f, p = %.4g\n', r_2017(ii,1), p_2017(ii,1));
end

% Replace these with your actual map vectors
map2023 = map_2023_all.net_mean(:,6);
map6    = map_pooled.net_mean(:,6);
map7    = map_2017_all.net_mean(:,7);

% Plot comparisons
figure;
maps_to_compare = {map6};
state_names = {'Pooled Dataset - State 6'};

for s = 1

    %subplot(1,2,s);
    hold on;

    x = maps_to_compare{s}(:);
    y = map2023(:);

    % Store scatter handles for the legend
    h_scatter = gobjects(7,1);

    % Plot parcels by resting-state network
    for n = 1:7

        idx = network_id == n;

        h_scatter(n) = scatter(x(idx), y(idx), 90, ...
            network_colors(n,:), 'filled', ...
            'DisplayName', network_names{n});

    end

    % Linear regression line
    valid = isfinite(x) & isfinite(y);

    b = polyfit(x(valid), y(valid), 1);

    x_fit = linspace(min(x(valid)), max(x(valid)), 100);
    y_fit = polyval(b, x_fit);

    % Regression line (visible, but excluded from legend)
    plot(x_fit, y_fit, 'k-', 'LineWidth', 2.5);

    % Pearson correlation
    r = corr(x(valid), y(valid), 'Type', 'Pearson');

    text(0.05, 0.95, sprintf('r = %.2f', r), ...
        'Units', 'normalized', ...
        'FontSize', 22, ...
        'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w');

    % Labels
    xlabel(char(string(state_names{s})), 'FontSize', 18);
    ylabel('Trento Dataset - State 6', 'FontSize', 18);

    % Axes
    axis square;
    grid on;

    set(gca, ...
        'FontSize', 14, ...
        'LineWidth', 1);

    % Legend: only the seven networks
    legend(h_scatter, network_names, ...
        'Location', 'best', ...
        'FontSize', 14);

end


