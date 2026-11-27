% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Statistical analysis to assess differences between patients and controls
% separatelly for the 2 datasets

clearvars; close all;

baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');
subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
DirData = fullfile(projectDir, 'derivatives/HMM');

mydir = uigetdir(pwd, 'Select the GitHub directory');
addpath(genpath(fullfile(mydir, 'HMM-MAR'))) % https://github.com/OHBA-analysis/HMM-MAR
addpath(genpath(fullfile(projectDir, 'code')))

K = (2:15); % no. states to be included in the analysis
num_states = 7; % IMPORTANT: choose the specific model (number of states)
k = num_states - 1; 
group = 2; % 1: Montreal

N = length(subjects); % no. subjects
zeromean = 0;
pca = 10; % 10 -> pca reduction; 0 -> no pca reduction
repetitions = 20;
reduced = 1; % 1 -> reduced number of volumes; 0 -> upsampling of dataset 2017

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



%% Select subgroup of subjects

if group == 2 % 2023
    Gamma{rep_i, k}((T*65)+1:(T*83),:)=[]; 
    Gamma{rep_i, k}((T*16)+1:(T*33),:)=[]; 
    FO{rep_i,k}(66:83,:) = [];    
    FO{rep_i,k}(17:33,:) = [];    
    vpath{rep_i,k}((T*65)+1:(T*83),:) = [];
    vpath{rep_i,k}((T*16)+1:(T*33),:) = [];
    switchingRate{rep_i,k}(66:83,:) = [];    
    switchingRate{rep_i,k}(17:33,:) = [];
    patients_ind = 16; 
    controls_ind = 32;
    subjects(66:83) = [];
    subjects(17:33) = [];
    subjects_list=[1:16 34:65];
    outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)], 'Dataset2023');
    mkdir(outputDir);  
    
elseif group == 1 % dataset2017
    Gamma{rep_i, k}((T*33)+1:(T*65),:)=[]; 
    Gamma{rep_i, k}(1:(T*16),:)=[]; 
    FO{rep_i,k}(34:65,:) = [];    
    FO{rep_i,k}(1:16,:) = [];    
    vpath{rep_i,k}((T*33)+1:(T*65),:) = [];
    vpath{rep_i,k}(1:(T*16),:) = [];
    switchingRate{rep_i,k}(34:65,:) = [];    
    switchingRate{rep_i,k}(1:16,:) = [];
    patients_ind = 17; 
    controls_ind = 18;
    subjects(34:65) = [];
    subjects(1:16) = [];
    subjects_list=[17:33 66:83];
    outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)], 'Dataset2017');
    mkdir(outputDir);    
  
elseif group == 3
    Gamma{rep_i, k}((T*65)+1:(T*83),:)=[]; 
    Gamma{rep_i, k}((T*55)+1:(T*56),:)=[]; 
    Gamma{rep_i, k}((T*47)+1:(T*52),:)=[]; 
    Gamma{rep_i, k}((T*44)+1:(T*45),:)=[]; 
    Gamma{rep_i, k}((T*16)+1:(T*42),:)=[]; 
    FO{rep_i,k}(66:83,:) = [];    
    FO{rep_i,k}(56,:) = [];    
    FO{rep_i,k}(48:52,:) = [];    
    FO{rep_i,k}(45,:) = [];    
    FO{rep_i,k}(17:42,:) = [];    
    vpath{rep_i,k}((T*65)+1:(T*83),:) = [];
    vpath{rep_i, k}((T*55)+1:(T*56),:)=[]; 
    vpath{rep_i, k}((T*47)+1:(T*52),:)=[]; 
    vpath{rep_i, k}((T*44)+1:(T*45),:)=[]; 
    vpath{rep_i,k}((T*16)+1:(T*42),:) = [];
    switchingRate{rep_i,k}(66:83,:) = [];    
    switchingRate{rep_i,k}(56,:) = [];    
    switchingRate{rep_i,k}(48:52,:) = [];    
    switchingRate{rep_i,k}(45,:) = [];
    switchingRate{rep_i,k}(17:42,:) = [];
    patients_ind = 16; 
    controls_ind = 16;
    subjects(66:83) = [];
    subjects(56,:) = [];    
    subjects(48:52,:) = [];    
    subjects(45,:) = []; 
    subjects(17:42) = [];
    subjects_list=[1:16 43:44 46:47 53:55 57:65];
    outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)], 'Dataset2023_AgeMatched');
    mkdir(outputDir);    
end  


    
%% Statistical analysis: define parameters and confounds

Tsubject = repmat(T,1,length(subjects))';
controls = ones(controls_ind,1);
patients = zeros(patients_ind,1);
Y = [ones(length(subjects),1) [patients; controls]];
thr = 0.025;
options.Nperm = 1000;
options.testP = 0; % to include the transition probability testing
options.grouplevel = 1; 
options.subjectlevel = 0; 

for sub_i = 1:length(subjects)
    
    if ~exist(fullfile(projectDir, 'derivatives/fmriprep', ['sub-' subjects{sub_i} '_task-rest_run-1_desc-confounds_timeseries.tsv']),'file')
        tmpregressor = readtable(fullfile(projectDir, 'derivatives/fmriprep', ['sub-' subjects{sub_i} '_task-resting_desc-confounds_timeseries.tsv']), 'FileType', 'delimitedtext');
    else
        tmpregressor = readtable(fullfile(projectDir, 'derivatives/fmriprep',['sub-' subjects{sub_i} '_task-rest_run-1_desc-confounds_timeseries.tsv']), 'FileType', 'delimitedtext');
    end        
    motion(sub_i,1) = nanmean(tmpregressor.framewise_displacement);
    clear tmpregressor    
end

if group == 2 % dataset 2023
    regressor_age = readtable(fullfile(projectDir, 'code/subjects_age_2023.xlsx'));
elseif group == 1
    regressor_age = readtable(fullfile(projectDir, 'code/subjects_age_2017.xlsx'));
end
regressor_age = regressor_age.age;
  
options.confounds = [motion regressor_age];



%% Run the statistical analysis

tests = hmmtest_modified_dataset(Gamma{rep_i, k},repmat(T,length(subjects),1),Tsubject,Y,options,hmm{rep_i, k}, [], vpath{rep_i,k}, group);

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



%% Plot the results

% 1. Gamma
figure;
imagesc(Gamma{rep_i, k});
xlabel('states', 'FontSize', 18);
ylabel('trials','FontSize', 18);
saveas(gca, fullfile(outputDir, 'Timecourses_States.bmp'))




%% Switching rate

if group == 2
    switching = switchingRate{rep_i,k};
elseif group == 1
        switching = switchingRate{rep_i,k};
        switching = switching / 0.785;
end

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

saveas(gca, fullfile(outputDir, 'Switching_Rate.bmp'))






%% Lifetime - vpath estimates on the group
lifetime = fdr_q(state_i+3:state_i+(k+3));

vector = vpath{rep_i,k};
% Calculate the number of segments
numSegments = numel(vector) / T;

for i = 1:numSegments
    startIndex = (i - 1) * T + 1;
    endIndex = i * T;
    lifeTime_mean = getStateLifeTimes(vector(startIndex:endIndex),T, hmm{rep_i, k}.train);
    for ii = 1:k+1
        current_mean = mean(lifeTime_mean{ii}, 1);
        if isempty(current_mean)
            test(i,ii) = NaN;  % or any other appropriate value
        else
            test(i,ii) = mean(current_mean,2);
        end
    end    
end

% Convert lifetime from samples to seconds
TR = ones(numSegments,1);
if group == 1
    TR(:,1) = 0.785;
end
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
xlim([0.5 k+1.5]); ylim([0 50]);
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
saveas(gca, fullfile(outputDir, 'FO_GroupAnalysis.bmp'))



close all;

