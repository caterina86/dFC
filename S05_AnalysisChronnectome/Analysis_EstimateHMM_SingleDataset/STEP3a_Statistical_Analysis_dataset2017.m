% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Statistical analysis to assess differences between patients and controls
% at the level of the chronnectome

clearvars; close all;

baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
DirData = fullfile(projectDir, 'derivatives/HMM');

mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR
addpath(genpath(fullfile(projectDir, 'code')))

K = 7;
num_states = 7; % IMPORTANT: choose the specific model (number of states)
k = num_states - 1; 

N = length(subjects); % no. subjects
zeromean = 0;
pca = 10; % 10 -> pca reduction; 0 -> no pca reduction
repetitions = 20;
% number of volumes per subject
T = 464; 

resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K)) '_train_test_dataset2017'], num2str(pca), ['zeromean' num2str(zeromean)]);
load(fullfile(resultsDir, 'HMMs_GROUP.mat'));



%% Find the minimum value and its indices, based on the Free Energy

[minValue_fe, minIndex_fe] = min(cell2mat(fe(:,:)));
[rowIndex_fe, colIndex_fe] = ind2sub(size(fe), minIndex_fe);
rep_i = rowIndex_fe(1,1); % repetition chosen based on the lowest free energy

disp(['repetition chosen ' num2str(rep_i)])
disp(['number states chosen ' num2str(K)])

state_i = K; % number of states

patients_ind = 17; 
controls_ind = 18;
subjects_list=1:35;
outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)], ['Number_Subjects_' num2str(length(subjects))]);
mkdir(outputDir);  



%% Statistical analysis: define parameters and confounds

Tsubject = repmat(T,1,length(subjects))';
controls = ones(controls_ind,1);
patients = zeros(patients_ind,1);
Y = [ones(length(subjects),1) [patients; controls]];
thr = 0.025; % two-tailed test
options.Nperm = 1000;
options.testP = 0; % to include the transition probability testing
options.grouplevel = 1; 
options.subjectlevel = 0; 

% define confounds for each subject
for sub_i = 1:length(subjects)
    
    if ~exist(fullfile(projectDir, 'derivatives/fmriprep', ['sub-' subjects{sub_i} '_task-rest_run-1_desc-confounds_timeseries.tsv']),'file')
        tmpregressor = readtable(fullfile(projectDir, 'derivatives/fmriprep', ['sub-' subjects{sub_i} '_task-resting_desc-confounds_timeseries.tsv']), 'FileType', 'delimitedtext');
    else
        tmpregressor = readtable(fullfile(projectDir, 'derivatives/fmriprep',['sub-' subjects{sub_i} '_task-rest_run-1_desc-confounds_timeseries.tsv']), 'FileType', 'delimitedtext');
    end        
    motion(sub_i,1) = nanmean(tmpregressor.framewise_displacement);
    clear tmpregressor    
end

regressor_age = readtable(fullfile(projectDir, 'code/subjects_age_2017.xlsx'));
regressor_age = regressor_age.age;

options.confounds = [motion regressor_age];


%% Run the statistical analysis

% test on the group estimates
tests = hmmtest_modified_dataset(Gamma{rep_i, 1},repmat(T,length(subjects),1),Tsubject,Y,options,hmm{rep_i, 1}, [], vpath{rep_i,1}, 1);

if options.testP > 0
    % correct for multiple comparisons
    p_values = [tests.grouplevel.p_fractional_occupancy; tests.grouplevel.p_aggr_fractional_occupancy; tests.grouplevel.p_switching_rate; tests.grouplevel.lifetime'; reshape(tests.grouplevel.p_transProbMat, [],1)];
else
    % correct for multiple comparisons
    p_values = [tests.grouplevel.p_fractional_occupancy; tests.grouplevel.p_aggr_fractional_occupancy; tests.grouplevel.p_switching_rate; tests.grouplevel.lifetime'];
    observed_values_fo = [tests.grouplevel.stat_fractional_occupancy];
    observed_values_lifetime = [tests.grouplevel.stat_lifetime];
end

[fdr_q] = mafdr(p_values, 'BHFDR', true);
fractional_occupancy = fdr_q(1:state_i);
fractional_occupancy_stat = observed_values_fo(2, 1:state_i);
lifetime_stat = observed_values_lifetime(2, 1:state_i);

aggr_fractional_occupancy = fdr_q(state_i+1);
switching_rate = fdr_q(state_i+2);
disp(['p-value switching-rate : ' num2str(switching_rate)])
lifetime = fdr_q(state_i+3:state_i+(k+3));
if options.testP > 0
    transProbMat = reshape(fdr_q(state_i+k+4:end),state_i,state_i);
end





%% Plot the results
switching = switchingRate{rep_i,1};
switching = switching / 0.785;

% 2. Switching rate
figure; 
plot(switching, 'LineWidth', 6);
xlabel('subjects', 'FontSize', 18); ylabel('Switching Rate (HZ)', 'FontSize', 18)
ylim([0 0.3]); xlim([0 83])
ax = gca;
ax.FontSize = 16; % Set tick label font size
saveas(gca, fullfile(outputDir, 'SwitchingRate_Allsubejcts.bmp'))

figure;
switching_group = [zeros(length(patients),1); ones(length(controls),1)];
boxchart(switchingRate{rep_i,1}, 'GroupByColor', switching_group);
ylim([0 0.3]);
ylabel('Switching Rate (Hz)', 'FontSize', 18, 'FontName', 'Helvetica');
ax = gca;
ax.FontSize = 16; % Set tick label font size
saveas(gca, fullfile(outputDir, 'SwitchingRate_boxplot.bmp'))

patients_mean = (mean(switchingRate{rep_i,1}(1:patients_ind),1));
controls_mean = (mean(switchingRate{rep_i,1}(patients_ind+1:end),1));

figure; 
histogram(switchingRate{rep_i,1}(1:patients_ind), 20);
hold on;
histogram(switchingRate{rep_i,1}(patients_ind+1:end), 20);
xline(patients_mean, 'b'); hold on;
xline(controls_mean, 'r')
legend('patients', 'controls')


% 3. Lifetime - vpath estimates on the group
vector = vpath{rep_i,1};
% Calculate the number of segments
numSegments = numel(vector) / T;

for i = 1:numSegments
    startIndex = (i - 1) * T + 1;
    endIndex = i * T;
    lifeTime_mean = getStateLifeTimes(vector(startIndex:endIndex),T, hmm{rep_i, 1}.train);
    for ii = 1:k+1
        current_mean = mean(lifeTime_mean{ii}, 1);
        if isempty(current_mean)
            test(i,ii) = NaN;  % or any other appropriate value
        else
            test(i,ii) = mean(current_mean,2);
        end
    end    
end

% Initialize an empty vector to store all values
allValues_patients = [];
allValues_controls = [];
allRowNumbers_patients = [];
allRowNumbers_controls = [];

% Loop through each cell in the cell array
for columnNumber = 1:(k+1)    
    % Extract values from the current cell
    values_patients = test(1:patients_ind, columnNumber);
    values_controls = test(patients_ind+1:end, columnNumber);    
    % Append values to the overall vector
    allValues_patients = [allValues_patients; values_patients];
    allValues_controls = [allValues_controls; values_controls];    
    % Create a vector of row numbers corresponding to the current values
    rowNumbers_patients = repmat(columnNumber, size(values_patients));
    rowNumbers_controls = repmat(columnNumber, size(values_controls));    
    % Append row numbers to the overall vector
    allRowNumbers_patients = [allRowNumbers_patients; rowNumbers_patients];
    allRowNumbers_controls = [allRowNumbers_controls; rowNumbers_controls];  
end

table.subjects = [allValues_patients; allValues_controls];
table.group = [ones(length(allValues_patients),1); repmat(2,length(allValues_controls),1)];
table.state = [allRowNumbers_patients; allRowNumbers_controls];

figure; 
subplot(1,2,1)
scatter(1:numel(lifetime), lifetime, 'b', 'filled');  % 'b' for blue dots
hold on;
yline(thr, 'r', 'Threshold');
xlim([0 13]); ylim([0 1]);
xticks(1:state_i);
xlabel('states #')
ylabel('FDR p-values', 'FontSize', 18); 

subplot(1,2,2);
boxchart(table.state, table.subjects,'GroupByColor',table.group)
xticks(1:1:(k+1)); xlabel('states', 'FontSize', 18); xlim([0 k+2])
ylabel('Lifetime (number of timepoints)', 'FontSize', 20);%ylim([0 50])
lgd = legend('blind', 'sighted', 'Location','northoutside','Orientation','horizontal');
lgd.FontSize = 16;
ax = gca;
ax.FontSize = 16; % Set tick label font size
saveas(gca, fullfile(outputDir, 'Lifetime_GroupAnalysis.bmp'))


% 4. Fractional Occupancy
figure; 
scatter(1:numel(fractional_occupancy), fractional_occupancy, 'b', 'filled');  % 'b' for blue dots
hold on;
yline(thr, 'r', 'Threshold');
xlabel('states #')
ylabel('FDR p-values', 'FontSize', 18)
xlim([0 k+1]); ylim([0 1])
xticks(1:state_i);

table.subjects = reshape(FO{rep_i,1}(:,:),[],1);
table.group = repmat([ones(patients_ind,1); 2*ones(controls_ind,1)],(k+1),1);
table.state = [];
for state = 1:(state_i)
    table.state = [table.state; repmat(state, length(subjects), 1)];
end

% BoxPlot
table.subjects = reshape(FO{rep_i,1}(:,:),[],1);
table.group = repmat([ones(patients_ind,1); ...
                      2*ones(controls_ind,1)],(k+1),1);
table.state = [];
for state = 1:(state_i)
    table.state = [table.state; repmat(state,length(subjects),1)];
end

% Colours
col_EB = [0.2 0.45 0.75];
col_SC = [0.85 0.3 0.25];

% Figure
figure('Color','w','Position',[100 100 1200 550]);
hold on;

for s = 1:(k+1)

    % Extract data for each group and state
    data_EB = double(table.subjects( ...
        table.state == s & table.group == 1));

    data_SC = double(table.subjects( ...
        table.state == s & table.group == 2));

    % Boxplot positions
    pos_EB = s - 0.25;
    pos_SC = s + 0.25;

    % Blind group
    boxchart(repmat(pos_EB,size(data_EB)), data_EB, 'BoxWidth',0.3, ...
        'BoxFaceColor',col_EB, 'BoxFaceAlpha',0.4, ...
        'WhiskerLineColor',col_EB, 'MarkerColor',col_EB, ...
        'LineWidth',1.2);
    % Individual observations - Blind
    scatter(pos_EB + 0.08*(rand(size(data_EB))-0.5), ...
        data_EB, 25, col_EB, 'filled', ...
        'MarkerFaceAlpha',0.6, ...
        'MarkerEdgeColor','none');

    % Sighted group
    boxchart(repmat(pos_SC,size(data_SC)), data_SC, ...
        'BoxWidth',0.3, 'BoxFaceColor',col_SC, ...
        'BoxFaceAlpha',0.4, 'WhiskerLineColor',col_SC, ...
        'MarkerColor',col_SC, 'LineWidth',1.2);
    % Individual observations - Sighted
    scatter(pos_SC + 0.08*(rand(size(data_SC))-0.5), ...
        data_SC, 25, col_SC, 'filled', ...
        'MarkerFaceAlpha',0.6, ...
        'MarkerEdgeColor','none');
end

% Axes
xlabel('States','FontSize',26);
ylabel('Fractional occupancy (%)','FontSize',26);
xticks(1:(k+1));
xlim([0.5 (k+1)+0.5]);
ylim([0 1]);
ax = gca;
ax.FontSize = 26;
ax.FontName = 'Helvetica';
ax.LineWidth = 1.2;
box off;

% Legend
h1 = patch(NaN,NaN,col_EB, ...
    'FaceAlpha',0.4,'EdgeColor',col_EB);

h2 = patch(NaN,NaN,col_SC, ...
    'FaceAlpha',0.4,'EdgeColor',col_SC);

lgd = legend([h1 h2],{'Blind','Sighted'}, ...
    'Location','northoutside', ...
    'Orientation','horizontal');

lgd.FontSize = 26;
saveas(gca, fullfile(outputDir, 'FO_GroupAnalysis.bmp'))

close all;