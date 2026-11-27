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

K = (2:15); % no. states to be included in the analysis
num_states = 7; %7; % IMPORTANT: choose the specific model (number of states)
k = num_states - 1; 

N = length(subjects); % no. subjects
zeromean = 0;
pca = 10; % 10 -> pca reduction; 0 -> no pca reduction
repetitions = 20;
% number of volumes per subject
T = 464; 

resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca), ['zeromean' num2str(zeromean)]);
load(fullfile(resultsDir, 'HMMs_GROUP.mat'));



%% Find the minimum value and its indices, based on the Free Energy

[minValue_fe, minIndex_fe] = min(cell2mat(fe(:,:)));
[rowIndex_fe, colIndex_fe] = ind2sub(size(fe), minIndex_fe);
rep_i = rowIndex_fe(1,k); % repetition chosen based on the lowest free energy

disp(['repetition chosen ' num2str(rep_i)])
disp(['number states chosen ' num2str(K(k))])

load(fullfile(resultsDir, ['HMMs_r' num2str(rep_i) '_SUBJECTS.mat']));
state_i = K(k); % number of states



%% Select subjects

patients_ind = 33; 
controls_ind = 50;
scanner = [ones(16,1); zeros(17,1); ones(32,1); zeros(18,1)];
subjects_list=1:83;
outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)], ['Number_Subjects_' num2str(length(subjects))]);
mkdir(outputDir);  


%% Statistical analysis: define parameters and confounds

Tsubject = repmat(T,1,length(subjects))';
controls = ones(controls_ind,1);
patients = zeros(patients_ind,1);
Y = [ones(length(subjects),1) [patients; controls]];
thr = 0.05; % two-tailed test
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

regressor_age = readtable(fullfile(projectDir, 'code/subjects_age.xlsx'));
regressor_age = regressor_age.age;

options.confounds = [motion scanner regressor_age];




%% Run the statistical analysis

% test on the group estimates
tests = hmmtest_modified(Gamma{rep_i, k},repmat(T,length(subjects),1),Tsubject,Y,options,hmm{rep_i, k}, [], vpath{rep_i,k});

if options.testP > 0
    % correct for multiple comparisons
    p_values = [tests.grouplevel.p_fractional_occupancy; tests.grouplevel.p_aggr_fractional_occupancy; tests.grouplevel.p_switching_rate; tests.grouplevel.lifetime'; reshape(tests.grouplevel.p_transProbMat, [],1)];
else
    % correct for multiple comparisons
    p_values = [tests.grouplevel.p_fractional_occupancy; tests.grouplevel.p_aggr_fractional_occupancy; tests.grouplevel.p_switching_rate; tests.grouplevel.lifetime'];
    observed_values = [tests.grouplevel.stat_fractional_occupancy];
end

[fdr_q] = mafdr(p_values, 'BHFDR', true);
fractional_occupancy = fdr_q(1:state_i);
fractional_occupancy_stat = observed_values(2,1:state_i);

aggr_fractional_occupancy = fdr_q(state_i+1);
switching_rate = fdr_q(state_i+2);
disp(['p-value switching-rate : ' num2str(switching_rate)])
lifetime = fdr_q(state_i+3:state_i+(k+3));
if options.testP > 0
    transProbMat = reshape(fdr_q(state_i+k+4:end),state_i,state_i);
end



%% Plot the results

% 1. Gamma
figure;
imagesc(Gamma{rep_i, k});
xlabel('states', 'FontSize', 18);
ylabel('trials','FontSize', 18);
saveas(gca, fullfile(outputDir, 'Timecourses_States.bmp'))




%% Switching rate
switching = switchingRate{rep_i,k};
% Correct for differences in TR
switching(17:33) = switching(17:33) / 0.785;
switching(66:83) = switching(66:83) / 0.785;

% Colours
col_EB = [0.2 0.45 0.75];
col_SC = [0.85 0.3 0.25];

% Plot individual switching rates
figure('Color','w','Position',[100 100 1100 500]);
hold on;

% Dataset 1: circles
scatter(1:16, switching(1:16), 150, 'o', 'MarkerFaceColor',col_EB, 'MarkerEdgeColor',col_EB);
scatter(34:65, switching(34:65), 150, 'o', 'MarkerFaceColor',col_SC, 'MarkerEdgeColor',col_SC);
% Dataset 2: stars
scatter(17:33, switching(17:33), 150, '*', 'MarkerEdgeColor',col_EB, 'LineWidth',1.5);
scatter(66:83, switching(66:83), 150, '*', 'MarkerEdgeColor',col_SC, 'LineWidth',1.5);

xlabel('Subjects','FontSize',26);
ylabel('Switching rate (s^{-1})','FontSize',26);
xlim([0 84]);
ylim([0 0.3]);
ax = gca;
ax.FontSize = 26;
ax.LineWidth = 1.2;
box off;

% Legend
h1 = plot(NaN,NaN,'o','MarkerSize',12, 'MarkerFaceColor',col_EB,'MarkerEdgeColor',col_EB, 'LineStyle','none');
h2 = plot(NaN,NaN,'o','MarkerSize',12, 'MarkerFaceColor',col_SC,'MarkerEdgeColor',col_SC, 'LineStyle','none');
h3 = plot(NaN,NaN,'o','MarkerSize',12, 'MarkerFaceColor','k','MarkerEdgeColor','k', 'LineStyle','none');
h4 = plot(NaN,NaN,'*','MarkerSize',16, 'Color','k','LineWidth',2, 'LineStyle','none');
lgd = legend([h1 h2 h3 h4], {'Early Blind','Sighted','Trento','Montreal'}, 'Location','best','FontSize',26);

lgd.ItemTokenSize = [50 30];
saveas(gca, fullfile(outputDir, 'Individual_SR.bmp'));



% ViolinPlot of the switching rate
patients_mean_sr = (mean(switching(1:patients_ind),1));
controls_mean_sr = (mean(switching(patients_ind+1:end),1));

% Data
data_EB = double(switching(1:patients_ind));
data_SC = double(switching(patients_ind+1:end));

% Colours
% col_EB = [0.2 0.45 0.75];
% col_SC = [0.85 0.3 0.25];
% 
% figure('Color','w','Position',[100 100 750 550]);
% hold on;
% % Density estimation
% [y_EB, x_EB] = ksdensity(data_EB);
% [y_SC, x_SC] = ksdensity(data_SC);
% % Normalize violin widths
% violin_width = 0.35;
% y_EB = y_EB / max(y_EB) * violin_width;
% y_SC = y_SC / max(y_SC) * violin_width;
% % Draw violins
% patch([1-y_EB, fliplr(1+y_EB)], [x_EB, fliplr(x_EB)], col_EB, ...
%       'FaceAlpha',0.3,'EdgeColor',col_EB,'LineWidth',1.5);
% 
% patch([2-y_SC, fliplr(2+y_SC)], [x_SC, fliplr(x_SC)], col_SC, ...
%       'FaceAlpha',0.3,'EdgeColor',col_SC,'LineWidth',1.5);
% 
% % Individual observations
% scatter(1 + 0.12*(rand(size(data_EB))-0.5), data_EB, ...
%     60, col_EB,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
% scatter(2 + 0.12*(rand(size(data_SC))-0.5), data_SC, ...
%     60, col_SC,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
% 
% % Median lines
% line([1-0.12 1+0.12], [median(data_EB) median(data_EB)], ...
%     'Color','w','LineWidth',2.5);
% line([2-0.12 2+0.12], [median(data_SC) median(data_SC)], ...
%     'Color','w','LineWidth',2.5);
% 
% % Axes
% xlim([0.4 2.6]);
% ylim([0 0.3]);
% xticks([1 2]);
% xticklabels({'Early blind','Sighted controls'});
% ylabel('Switching rate (s^{-1})','FontSize',26);
% ax = gca;
% ax.FontSize = 26;
% ax.FontName = 'Helvetica';
% ax.LineWidth = 1.2;
% box off;

% Colours
col_EB = [0.2 0.45 0.75];
col_SC = [0.85 0.3 0.25];

figure('Color','w','Position',[100 100 750 550]);
hold on;

% Early blind group
boxchart(ones(size(data_EB)), data_EB, ...
    'BoxWidth',0.45, ...
    'BoxFaceColor',col_EB, ...
    'BoxFaceAlpha',0.4, ...
    'WhiskerLineColor',col_EB, ...
    'MarkerStyle','none', ...
    'LineWidth',1.5);

% Individual observations - Blind
scatter(1 + 0.12*(rand(size(data_EB))-0.5), ...
    data_EB, 35, col_EB, 'filled', ...
    'MarkerFaceAlpha',0.6, ...
    'MarkerEdgeColor','none');

% Sighted controls
boxchart(2*ones(size(data_SC)), data_SC, ...
    'BoxWidth',0.45, ...
    'BoxFaceColor',col_SC, ...
    'BoxFaceAlpha',0.4, ...
    'WhiskerLineColor',col_SC, ...
    'MarkerStyle','none', ...
    'LineWidth',1.5);

% Individual observations - Sighted
scatter(2 + 0.12*(rand(size(data_SC))-0.5), ...
    data_SC, 35, col_SC, 'filled', ...
    'MarkerFaceAlpha',0.6, ...
    'MarkerEdgeColor','none');
% Axes
xlim([0.4 2.6]);
ylim([0 0.3]);
xticks([1 2]);
xticklabels({'Early blind','Sighted controls'});

ylabel('Switching rate (s^{-1})','FontSize',26);

ax = gca;
ax.FontSize = 26;
ax.FontName = 'Helvetica';
ax.LineWidth = 1.2;
box off;

saveas(gca, fullfile(outputDir, 'Group_SR.bmp'))



%% Lifetime - Viterbi estimates on the group
vector = vpath{rep_i,k};
% Number of segments (subjects)
numSegments = numel(vector) / T;
% Initialize
test = NaN(numSegments, k+1);

for i = 1:numSegments
    startIndex = (i - 1) * T + 1;
    endIndex = i * T;
    lifeTime_mean = getStateLifeTimes( ...
        vector(startIndex:endIndex), T, hmm{rep_i,k}.train);

    for ii = 1:k+1
        current_mean = mean(lifeTime_mean{ii}, 1);
        if isempty(current_mean)
            test(i,ii) = NaN;
        else
            test(i,ii) = mean(current_mean,2);
        end
    end
end

% Convert lifetime from samples to seconds
TR = ones(numSegments,1);
TR([17:33 66:end]) = 0.785;
test = test .* TR;

% Separate groups
allValues_patients = [];
allValues_controls = [];
allRowNumbers_patients = [];
allRowNumbers_controls = [];

for columnNumber = 1:(k+1)
    values_patients = test(1:patients_ind, columnNumber);
    values_controls = test(patients_ind+1:end, columnNumber);
    allValues_patients = [allValues_patients; values_patients];
    allValues_controls = [allValues_controls; values_controls];
    rowNumbers_patients = repmat(columnNumber, size(values_patients));
    rowNumbers_controls = repmat(columnNumber, size(values_controls));
    allRowNumbers_patients = [allRowNumbers_patients; rowNumbers_patients];
    allRowNumbers_controls = [allRowNumbers_controls; rowNumbers_controls];
end

% Construct table
table.subjects = [allValues_patients; allValues_controls];
table.group = [ones(length(allValues_patients),1); ...
               repmat(2,length(allValues_controls),1)];
table.state = [allRowNumbers_patients; allRowNumbers_controls];


figure; 
scatter(1:numel(lifetime), lifetime, 'b', 'filled');  % 'b' for blue dots
hold on;
yline(thr, 'r', 'Threshold');
xlim([0 13]); ylim([0 1]);
xticks(1:state_i);
xlabel('states #')
ylabel('FDR p-values', 'FontSize', 18); 

% % Colours
% col_EB = [0.2 0.45 0.75];
% col_SC = [0.85 0.3 0.25];
% 
% figure('Color','w','Position',[100 100 1200 550]);
% hold on;
% 
% for s = 1:(k+1)
%     % Extract data for each group and state
%     data_EB = double(table.subjects(table.state == s & table.group == 1));
%     data_SC = double(table.subjects(table.state == s & table.group == 2));
%     % Violin positions
%     pos_EB = s - 0.18;
%     pos_SC = s + 0.18;
%     % Kernel density estimation
%     [d_EB,y_EB] = ksdensity(data_EB);
%     [d_SC,y_SC] = ksdensity(data_SC);
%     % Normalize violin widths
%     width = 0.15;
%     d_EB = d_EB/max(d_EB)*width;
%     d_SC = d_SC/max(d_SC)*width;
% 
%     % Draw violins
%     patch([pos_EB-d_EB fliplr(pos_EB+d_EB)], ...
%           [y_EB fliplr(y_EB)], col_EB, ...
%           'FaceAlpha',0.3,'EdgeColor',col_EB,'LineWidth',1);
% 
%     patch([pos_SC-d_SC fliplr(pos_SC+d_SC)], ...
%           [y_SC fliplr(y_SC)], col_SC, ...
%           'FaceAlpha',0.3,'EdgeColor',col_SC,'LineWidth',1);
% 
%     % Individual observations
%     scatter(pos_EB+0.08*(rand(size(data_EB))-0.5), ...
%         data_EB,25,col_EB,'filled', ...
%         'MarkerFaceAlpha',0.6,'MarkerEdgeColor','none');
% 
%     scatter(pos_SC+0.08*(rand(size(data_SC))-0.5), ...
%         data_SC,25,col_SC,'filled', ...
%         'MarkerFaceAlpha',0.6,'MarkerEdgeColor','none');
% 
%     % Median lines
%     line([pos_EB-0.08 pos_EB+0.08], ...
%         [median(data_EB,'omitnan') median(data_EB,'omitnan')], ...
%         'Color','w','LineWidth',2);
% 
%     line([pos_SC-0.08 pos_SC+0.08], ...
%         [median(data_SC,'omitnan') median(data_SC,'omitnan')], ...
%         'Color','w','LineWidth',2);
% 
% end
% 
% % Axes
% xticks(1:k+1);
% xlabel('States','FontSize',26);
% xlim([0.5 k+1.5]);
% ylabel('Mean lifetime (s)','FontSize',26);
% ax = gca;
% ax.FontSize = 26;
% ax.FontName = 'Helvetica';
% ax.LineWidth = 1.2;
% box off;
% 
% % Legend
% h1 = patch(NaN,NaN,col_EB,'FaceAlpha',0.3,'EdgeColor',col_EB);
% h2 = patch(NaN,NaN,col_SC,'FaceAlpha',0.3,'EdgeColor',col_SC);
% lgd = legend([h1 h2],{'Blind','Sighted'}, ...
%     'Location','northoutside','Orientation','horizontal');
% lgd.FontSize = 26;


% Colours
col_EB = [0.2 0.45 0.75];
col_SC = [0.85 0.3 0.25];

figure('Color','w','Position',[100 100 1200 550]);
hold on;

for s = 1:(k+1)

    % Extract data for each group and state
    data_EB = double(table.subjects( ...
        table.state == s & table.group == 1));

    data_SC = double(table.subjects( ...
        table.state == s & table.group == 2));

    % Boxchart positions
    pos_EB = s - 0.25;
    pos_SC = s + 0.25;

    % Blind group
    boxchart(repmat(pos_EB,size(data_EB)), data_EB, ...
        'BoxWidth',0.3, ...
        'BoxFaceColor',col_EB, ...
        'BoxFaceAlpha',0.4, ...
        'WhiskerLineColor',col_EB, ...
        'MarkerColor',col_EB, ...
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
xticks(1:k+1);
xlabel('States','FontSize',26);
xlim([0.5 k+1.5]);
ylabel('Mean lifetime (s)','FontSize',26);
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



saveas(gca, fullfile(outputDir, 'Lifetime_GroupAnalysis.bmp'))




%% Fractional Occupancy
figure; 
scatter(1:numel(fractional_occupancy), fractional_occupancy, 'b', 'filled');  % 'b' for blue dots
hold on;
yline(thr, 'r', 'Threshold');
xlabel('states #')
ylabel('FDR p-values', 'FontSize', 18)
xlim([0 k+1]); ylim([0 1])
xticks(1:state_i);

table.subjects = reshape(FO{rep_i,k}(:,:),[],1);
table.group = repmat([ones(patients_ind,1); 2*ones(controls_ind,1)],(k+1),1);
table.state = [];
for state = 1:(state_i)
    table.state = [table.state; repmat(state, length(subjects), 1)];
end


% ViolinPlot
% col_EB = [0.2 0.45 0.75];
% col_SC = [0.85 0.3 0.25];
% 
% figure('Color','w','Position',[100 100 1200 550]);
% hold on;
% 
% for s = 1:(k+1)
% 
%     % Extract data for each group and state
%     data_EB = double(table.subjects(table.state == s & table.group == 1));
%     data_SC = double(table.subjects(table.state == s & table.group == 2));
%     % Violin positions
%     pos_EB = s - 0.18;
%     pos_SC = s + 0.18;
%     % Kernel density estimation
%     [d_EB,y_EB] = ksdensity(data_EB, 'Support', [0 1]);
%     [d_SC,y_SC] = ksdensity(data_SC, 'Support', [0 1]);
% 
%     % Normalize widths
%     width = 0.15;
%     d_EB = d_EB/max(d_EB)*width;
%     d_SC = d_SC/max(d_SC)*width;
% 
%     % Draw violins
%     patch([pos_EB-d_EB fliplr(pos_EB+d_EB)], [y_EB fliplr(y_EB)], col_EB, ...
%           'FaceAlpha',0.3,'EdgeColor',col_EB,'LineWidth',1);
% 
%     patch([pos_SC-d_SC fliplr(pos_SC+d_SC)], [y_SC fliplr(y_SC)], col_SC, ...
%           'FaceAlpha',0.3,'EdgeColor',col_SC,'LineWidth',1);
% 
%     % Individual observations
%     scatter(pos_EB+0.08*(rand(size(data_EB))-0.5), ...
%         data_EB,25,col_EB,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeColor','none');
% 
%     scatter(pos_SC+0.08*(rand(size(data_SC))-0.5), ...
%         data_SC,25,col_SC,'filled','MarkerFaceAlpha',0.6,'MarkerEdgeColor','none');
% 
%     % Median lines
%     line([pos_EB-0.08 pos_EB+0.08], ...
%         [median(data_EB) median(data_EB)], ...
%         'Color','w','LineWidth',2);
%     line([pos_SC-0.08 pos_SC+0.08], ...
%         [median(data_SC) median(data_SC)], ...
%         'Color','w','LineWidth',2);
% 
% end
% 
% % Axes
% xlabel('States','FontSize',26);
% ylabel('Fractional occupancy (%)','FontSize',26);
% xticks(1:(k+1));
% xlim([0.5 (k+1)+0.5]);
% ylim([0 1]);
% ax = gca;
% ax.FontSize = 26;
% ax.FontName = 'Helvetica';
% ax.LineWidth = 1.2;
% box off;
% 
% % Legend
% h1 = patch(NaN,NaN,col_EB,'FaceAlpha',0.3,'EdgeColor',col_EB);
% h2 = patch(NaN,NaN,col_SC,'FaceAlpha',0.3,'EdgeColor',col_SC);
% lgd = legend([h1 h2],{'Blind','Sighted'}, ...
%     'Location','northoutside','Orientation','horizontal');
% lgd.FontSize = 26;
% saveas(gca, fullfile(outputDir, 'FO_GroupAnalysis.bmp'))



% BoxPlot
table.subjects = reshape(FO{rep_i,k}(:,:),[],1);
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

% Save figure
saveas(gca, fullfile(outputDir,'FO_GroupAnalysis.bmp'));




close all;




%% Extract the Transition Probabilites
% diagonal values = nan to extract the topological plot without the diagonal values

if options.testP > 0

    % Plot the significant results of the transition probability
    data = transProbMat; % FDR-corrected p-values
    data(data > thr) = NaN;
    [nr,nc] = size(data);
    figure;
    pcolor([data nan(nr,1); nan(1,nc+1)]);
    colormap(flip(jet))
    % Calculate the middle positions for x-axis and y-axis ticks
    x_middle_positions = 1.5:1:(state_i+0.5);
    y_middle_positions = 1.5:1:(state_i+0.5);
    % Set x-axis and y-axis ticks in the middle of each cell
    xticks(x_middle_positions);
    yticks(y_middle_positions);
    xticklabels(1:state_i);  
    yticklabels(1:state_i);  
    xlabel('To state', 'FontSize', 16); ylabel('From state', 'FontSize', 16)
    set(gca, 'ydir', 'reverse');
    colorbar;
    title('FDR corrected TPM', 'FontSize', 16);
    saveas(gca, fullfile(outputDir, 'TPM_GroupAnalysis.bmp'))

end

for sub_i = 1:length(subjects)
    % transition probability matrix of the group, without the diagonal values  
    figure; 
    tpm_subjects(:,:,sub_i) = getTransProbs(HMMs_dualregr{k,subjects_list(sub_i)});
    close all;

end


%% Extract the patient TPM
data = tpm_subjects(:,:,1:patients_ind);
mean_matrix_patients = median(data,3);

figure; 
subplot(1,2,1)
pairwise_distances = pdist(mean_matrix_patients);
tree = linkage(pairwise_distances, 'average');
% Plot dendrogram and capture its handle
[h, ~,~] = dendrogram(tree);
set(h, 'LineWidth', 2); % Adjust line thickness
set(gca, 'FontSize', 14); % Adjust the size of the numbers on both axes

subplot(1,2,2)
mean_matrix_patients(1:length(mean_matrix_patients)+1:end) = NaN;
data = mean_matrix_patients;
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
colormap(jet)
% Calculate the middle positions for x-axis and y-axis ticks
x_middle_positions = 1.5:1:(state_i+0.5);
y_middle_positions = 1.5:1:(state_i+0.5);
% Set x-axis and y-axis ticks in the middle of each cell
xticks(x_middle_positions);
yticks(y_middle_positions);
xticklabels(1:state_i);  
yticklabels(1:state_i);  
set(gca, 'ydir', 'reverse');
xlabel('To state', 'FontSize', 20); ylabel('From state', 'FontSize', 20); 
colorbar; caxis([0 0.35]);
title('TPM Patients', ['Number subjects ' num2str(length(subjects))], 'FontSize', 20);
pbaspect([size(data, 2) size(data, 1) 1]);
saveas(gca, fullfile(outputDir, 'MeanTPM_patients.bmp'))


% Set the percentile threshold (e.g., 20%)
percentile_threshold = 80;
% Calculate the threshold value
threshold_value = prctile(mean_matrix_patients(:), percentile_threshold);
% Threshold the matrix to keep values above the threshold
thresholded_matrix = mean_matrix_patients;
thresholded_matrix(thresholded_matrix <= threshold_value) = NaN;

figure;
data = thresholded_matrix;
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
colormap(jet)
% Calculate the middle positions for x-axis and y-axis ticks
x_middle_positions = 1.5:1:(state_i+0.5);
y_middle_positions = 1.5:1:(state_i+0.5);
% Set x-axis and y-axis ticks in the middle of each cell
xticks(x_middle_positions);
yticks(y_middle_positions);
xticklabels(1:state_i);  
yticklabels(1:state_i);  
set(gca, 'ydir', 'reverse');
caxis([0 0.07]);
colorbar; xlabel('To state', 'FontSize', 20); ylabel('From state', 'FontSize', 20); 
title('Top 20% probabilities Patients', ['Number subjects ' num2str(length(subjects))], 'FontSize', 20);
saveas(gca, fullfile(outputDir, '20%Prob_Matrix_patients.bmp'))


% Create a directed graph object from the matrix
G = digraph(data);
% Remove edges with NaN weights
nanEdges = isnan(G.Edges.Weight);
G = rmedge(G, find(nanEdges));
% Plot the directed graph
base_line_width = 1;
suggested_scaling_factor = 150;
% Plot the directed graph with adjusted arrow thickness and larger arrows
figure;
h = plot(G, 'LineWidth', base_line_width + suggested_scaling_factor * (G.Edges.Weight - min(G.Edges.Weight)), 'MarkerSize', 8, 'ArrowSize', 15);
% Customize node labels directly in the plot function
labelNode = arrayfun(@(x) {num2str(x)}, 1:numnodes(G))';
h.NodeLabel = labelNode;
h.NodeFontSize = 18;  % Adjust the font size as needed
saveas(gca, fullfile(outputDir, '20%Prob_Matrix_patients_graph.bmp'))



%% Extract the sighted TPM

data = tpm_subjects(:,:,patients_ind+1:end);
mean_matrix_controls = median(data,3);

figure; 
subplot(1,2,1)
pairwise_distances = pdist(mean_matrix_controls);
tree = linkage(pairwise_distances, 'average');
% Plot dendrogram and capture its handle
[h,~,~] = dendrogram(tree);
set(h, 'LineWidth', 2); % Adjust line thickness
set(gca, 'FontSize', 14); % Adjust the size of the numbers on both axes

subplot(1,2,2)
mean_matrix_controls(1:length(mean_matrix_controls)+1:end) = NaN; % delete the diagonal
data = mean_matrix_controls;
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
colormap(jet)
% Calculate the middle positions for x-axis and y-axis ticks
x_middle_positions = 1.5:1:(state_i+0.5);
y_middle_positions = 1.5:1:(state_i+0.5);
% Set x-axis and y-axis ticks in the middle of each cell
xticks(x_middle_positions);
yticks(y_middle_positions);
xticklabels(1:state_i);  
yticklabels(1:state_i);  
set(gca, 'ydir', 'reverse');
colorbar; xlabel('To state', 'FontSize', 20); ylabel('From state', 'FontSize', 20); 
caxis([0 0.35]);
title('TPM Controls', ['Number subjects ' num2str(length(subjects))], 'FontSize', 20);
pbaspect([size(data, 2) size(data, 1) 1]);
saveas(gca, fullfile(outputDir, 'MeanTPM_controls.bmp'))

% Calculate the threshold value
threshold_value = prctile(mean_matrix_controls(:), percentile_threshold);
% Threshold the matrix to keep values above the threshold
thresholded_matrix = mean_matrix_controls;
thresholded_matrix(thresholded_matrix <= threshold_value) = NaN;

figure;
data = thresholded_matrix;
[nr,nc] = size(data);
pcolor([data nan(nr,1); nan(1,nc+1)]);
colormap(jet)
% Calculate the middle positions for x-axis and y-axis ticks
x_middle_positions = 1.5:1:(state_i+0.5);
y_middle_positions = 1.5:1:(state_i+0.5);
% Set x-axis and y-axis ticks in the middle of each cell
xticks(x_middle_positions);
yticks(y_middle_positions);
xticklabels(1:state_i);  
yticklabels(1:state_i);  
set(gca, 'ydir', 'reverse');
colorbar; caxis([0 0.07]);
xlabel('To state', 'FontSize', 20); ylabel('From state', 'FontSize', 20); 
title('Top 20% probabilities Controls', ['Number subjects ' num2str(length(subjects))], 'FontSize', 20);
saveas(gca, fullfile(outputDir, '20%Prob_Matrix_controls.bmp'))

% Create a directed graph object from the matrix
G = digraph(data);
% Remove edges with NaN weights
nanEdges = isnan(G.Edges.Weight);
G = rmedge(G, find(nanEdges));
% Plot the directed graph
base_line_width = 1;
suggested_scaling_factor = 150;
% Plot the directed graph with adjusted arrow thickness and larger arrows
figure;
h = plot(G, 'LineWidth', base_line_width + suggested_scaling_factor * (G.Edges.Weight - min(G.Edges.Weight)), 'MarkerSize', 8, 'ArrowSize', 15);
% Customize node labels directly in the plot function
labelNode = arrayfun(@(x) {num2str(x)}, 1:numnodes(G))';
h.NodeLabel = labelNode;
h.NodeFontSize = 18;  % Adjust the font size as needed
saveas(gca, fullfile(outputDir, '20%Prob_Matrix_controls_graph.bmp'))


%% Difference between groups in the TPM

if options.testP > 0

    for sub_i = 1:length(subjects)
        % transition probability matrix of the group, without the diagonal values  
        figure; 
        tpm_subjects{sub_i} = HMMs_dualregr{k,subjects_list(sub_i)}.P; % transition probability without the diagonal
        close all;

    end
    % 
    % Extract the patient TPM
    data = tpm_subjects(1:patients_ind);
    % Ensure all matrices have the same size
    matrix_sizes = cellfun(@size, data, 'UniformOutput', false);
    assert(all(cellfun(@(s) isequal(s, matrix_sizes{1}), matrix_sizes)), 'All matrices must have the same size.');
    % Initialize the mean matrix
    mean_matrix = zeros(matrix_sizes{1});
    % Sum up all matrices
    for i = 1:numel(data)
        mean_matrix = mean_matrix + data{i};
    end
    % Calculate the mean by dividing by the number of matrices
    mean_matrix_patients = mean_matrix/ numel(data);

    % Extract the controls TPM
    data = tpm_subjects(patients_ind+1:end);
    % Ensure all matrices have the same size
    matrix_sizes = cellfun(@size, data, 'UniformOutput', false);
    assert(all(cellfun(@(s) isequal(s, matrix_sizes{1}), matrix_sizes)), 'All matrices must have the same size.');
    % Initialize the mean matrix
    mean_matrix = zeros(matrix_sizes{1});
    % Sum up all matrices
    for i = 1:numel(data)
        mean_matrix = mean_matrix + data{i};
    end
    % Calculate the mean by dividing by the number of matrices
    mean_matrix_controls = mean_matrix / numel(data);

    % plot the results
    data = transProbMat;
    data(data > thr) = NaN;
    diff_TPM = (mean_matrix_controls - mean_matrix_patients);
    diff_TPM(isnan(data))=nan;
    diff_TPM(diff_TPM==0)=nan;
    % Create a mask for NaN values
    nan_mask = isnan(diff_TPM);

    % Plot the matrix using imagesc with NaN values set to white
    figure;
    imagesc(diff_TPM, 'AlphaData', ~nan_mask); caxis([-0.05 0.05])
    % Set NaN values to transparent
    colormap jet;  % Set the colormap to 'jet'
    colorbar;  % Display colorbar
    set(gca, 'XTick', 1:size(diff_TPM, 2) + 0.5, 'XTickLabel', 1:size(diff_TPM, 2), ...
             'YTick', 1:size(diff_TPM, 1) + 0.5, 'YTickLabel', 1:size(diff_TPM, 1));

    % Display numerical values in the middle of each row and column
    for i = 1:size(diff_TPM, 1)
        for j = 1:size(diff_TPM, 2)
            if ~isnan(diff_TPM(i, j))
                text(j, i, num2str(round(diff_TPM(i, j),3)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
        end
    end

    % Add labels and title if necessary
    xlabel(' to states', 'FontSize', 18);
    ylabel('from states', 'FontSize', 18);
    title('TPM of difference controls - patients', ['Number subjects ' num2str(length(subjects))], 'FontSize', 20);
    saveas(gca, fullfile(outputDir, 'Difference_TPM.bmp'))

end

close all;
