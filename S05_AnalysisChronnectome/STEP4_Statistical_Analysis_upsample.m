% Statistical analysis to assess differences between patients and controls
% Caterina Pedersini, March 2024, caterina.pedersini@nyu.edu

clearvars; close all;

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

if pca == 0
    resultsDir = fullfile(DirData, 'PCA_UPSAMPLE', ['zeromean' num2str(zeromean)]);
else
    resultsDir = fullfile(DirData, 'PCA_UPSAMPLE', num2str(pca),['zeromean' num2str(zeromean)]);
end        

load(fullfile(resultsDir, 'HMMs_GROUP.mat'));

% Find the minimum value and its indices, based on the Free Energy
[minValue_fe, minIndex_fe] = min(cell2mat(fe(:,:)));
[rowIndex_fe, colIndex_fe] = ind2sub(size(fe), minIndex_fe);
rep_i = rowIndex_fe(1,k); % repetition chosen based on the lowest free energy

disp(['repetition chosen ' num2str(rep_i)])
disp(['number states chosen ' num2str(K(k))])

load(fullfile(resultsDir, ['HMMs_r' num2str(rep_i) '_SUBJECTS.mat']));
state_i = K(k); % number of states


%% Select subjects

if k == 4 || k == 5 || k == 6 || k == 7 || k == 8 % don't remove anyone
    patients_ind = 33; 
    controls_ind = 50;
    scanner = [ones(16,1); zeros(17,1); ones(32,1); zeros(18,1)];
    T = [repmat(480,16,1); repmat(370,17,1); repmat(480,32,1); repmat(370,18,1)];
    subjects_list=1:83;
    outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)], ['Number_Subjects_' num2str(length(subjects))]);
    mkdir(outputDir);     
elseif k == 10 % remove blnd05 and eb23
    Gamma{rep_i, k}((T*20)+1:(T*21),:)=[]; 
    Gamma{rep_i, k}((T*15)+1:(T*16),:)=[]; 
    FO{rep_i,k}(21,:) = [];    
    FO{rep_i,k}(16,:) = [];    
    vpath{rep_i,k}((T*20)+1:(T*21),:) = [];
    vpath{rep_i,k}((T*15)+1:(T*16),:) = [];
    switchingRate{rep_i,k}(21,:) = [];
    switchingRate{rep_i,k}(16,:) = [];
    patients_ind = 31; 
    controls_ind = 50;
    subjects(21) = [];
    subjects(16) = [];
    scanner = [ones(15,1); zeros(16,1); ones(32,1); zeros(18,1)];
    subjects_list=[1:15 17:20 22:83];
    outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)], ['Number_Subjects_' num2str(length(subjects))]);
    mkdir(outputDir);     
elseif k == 9 || k == 11 || k == 12 || k == 13 || k == 14 % remove blnd05, ctrl08 and ctrl10
    Gamma{rep_i, k}((T*74)+1:(T*75),:)=[]; 
    Gamma{rep_i, k}((T*72)+1:(T*73),:)=[]; 
    Gamma{rep_i, k}((T*20)+1:(T*21),:)=[]; 
    FO{rep_i,k}(75,:) = [];
    FO{rep_i,k}(73,:) = [];    
    FO{rep_i,k}(21,:) = [];    
    vpath{rep_i,k}((T*74)+1:(T*75),:) = [];
    vpath{rep_i,k}((T*72)+1:(T*73),:) = [];
    vpath{rep_i,k}((T*20)+1:(T*21),:) = [];
    switchingRate{rep_i,k}(75,:)= [];
    switchingRate{rep_i,k}(73,:) = [];
    switchingRate{rep_i,k}(21,:) = [];
    patients_ind = 32; 
    controls_ind = 48;
    subjects(75) = [];
    subjects(73) = [];
    subjects(21) = [];
    scanner = [ones(16,1); zeros(16,1); ones(32,1); zeros(16,1)];
    subjects_list=[1:20 22:72 74 76:83];
    outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)], ['Number_Subjects_' num2str(length(subjects))]);
    mkdir(outputDir);  
end
   

%% Statistical analysis

Tsubject = T;
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
    
    if ~exist(fullfile(baseDir, 'dataset2017/derivatives/fmriprep', ['sub-' subjects{sub_i}], 'func', ['sub-' subjects{sub_i} '_task-rest_run-1_desc-confounds_timeseries.tsv']),'file')
        tmpregressor = readtable(fullfile(baseDir, 'dataset2023/derivatives/fmriprep', ['sub-' subjects{sub_i}], 'func', ['sub-' subjects{sub_i} '_task-resting_desc-confounds_timeseries.tsv']), 'FileType', 'delimitedtext');
    else
        tmpregressor = readtable(fullfile(baseDir, 'dataset2017/derivatives/fmriprep', ['sub-' subjects{sub_i}], 'func', ['sub-' subjects{sub_i} '_task-rest_run-1_desc-confounds_timeseries.tsv']), 'FileType', 'delimitedtext');
    end        
    motion(sub_i,1) = nanmean(tmpregressor.framewise_displacement);
    clear tmpregressor    
end

regressor_age = readtable(fullfile(baseDir, 'dataset2023_2017/subjects_age.xlsx'));
regressor_age = regressor_age.age;

options.confounds = [motion scanner regressor_age];


%% test on the group estimates
tests = hmmtest_stats_upsample(Gamma{rep_i, k},T,Tsubject,Y,options,hmm{rep_i, k}, []);

if options.testP > 0
    % correct for multiple comparisons
    p_values = [tests.grouplevel.p_fractional_occupancy; tests.grouplevel.p_aggr_fractional_occupancy; tests.grouplevel.p_switching_rate; reshape(tests.grouplevel.p_transProbMat, [],1)];
else
    % correct for multiple comparisons
    p_values = [tests.grouplevel.p_fractional_occupancy; tests.grouplevel.p_aggr_fractional_occupancy; tests.grouplevel.p_switching_rate];
    observed_values_fo = [tests.grouplevel.stat_fractional_occupancy];
end

[fdr_q] = mafdr(p_values, 'BHFDR', true);
fractional_occupancy = fdr_q(1:state_i);
fractional_occupancy_stat = observed_values_fo(2, 1:state_i);

aggr_fractional_occupancy = fdr_q(state_i+1);
switching_rate = fdr_q(state_i+2);
disp(['p-value switching-rate : ' num2str(switching_rate)])


%% Plot the results

% Switching rate
switching = switchingRate{rep_i,k};

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
col_EB = [0.2 0.45 0.75];
col_SC = [0.85 0.3 0.25];

figure('Color','w','Position',[100 100 750 550]);
hold on;
% Density estimation
[y_EB, x_EB] = ksdensity(data_EB);
[y_SC, x_SC] = ksdensity(data_SC);
% Normalize violin widths
violin_width = 0.35;
y_EB = y_EB / max(y_EB) * violin_width;
y_SC = y_SC / max(y_SC) * violin_width;
% Draw violins
patch([1-y_EB, fliplr(1+y_EB)], [x_EB, fliplr(x_EB)], col_EB, ...
      'FaceAlpha',0.3,'EdgeColor',col_EB,'LineWidth',1.5);

patch([2-y_SC, fliplr(2+y_SC)], [x_SC, fliplr(x_SC)], col_SC, ...
      'FaceAlpha',0.3,'EdgeColor',col_SC,'LineWidth',1.5);

% Individual observations
scatter(1 + 0.12*(rand(size(data_EB))-0.5), data_EB, ...
    60, col_EB,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');
scatter(2 + 0.12*(rand(size(data_SC))-0.5), data_SC, ...
    60, col_SC,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none');

% Median lines
line([1-0.12 1+0.12], [median(data_EB) median(data_EB)], ...
    'Color','w','LineWidth',2.5);
line([2-0.12 2+0.12], [median(data_SC) median(data_SC)], ...
    'Color','w','LineWidth',2.5);

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

