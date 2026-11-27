

% Plot FO and FC in the same figure

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

thresholds = 4; %3.5;
sizes = {'Extent'}; %{'Intensity'}; %{'Extent'};

outputStatDir = fullfile(outputDir, 'NBS', sizes{1}, 'staticFC');
mkdir(outputStatDir);
FigureDir = fullfile(outputDir, 'NBS', sizes{1}, 'Figures');


%% Define the shuffle states and load the data

shuffle_state = 1;
outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFC_state' num2str(shuffle_state)]);

% load the data
load(fullfile(outputShuffle, 'tvals.mat'));
load(fullfile(outputShuffle, 'mean_vals.mat'));
tvals_perm_FC = tvals_perm;
perm_N = size(tvals_perm_FC,2);
tvals_perm_FC(tvals_perm_FC < thresholds) = 0; % remove all values of t < thr
perm_survives_FC = any(tvals_perm_FC ~= 0, 1);   % [1 x perm_N]
frac_survive_FC  = sum(perm_survives_FC)/perm_N;    

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFO_state' num2str(shuffle_state)]);

% load the data
load(fullfile(outputShuffle, 'tvals.mat'));
tvals_perm_FO = tvals_perm;
perm_N = size(tvals_perm_FO,2);
load(fullfile(outputShuffle, 'mean_vals.mat'));

tvals_perm_FO(tvals_perm_FO < thresholds) = 0; % remove all values of t < thr
perm_survives_FO = any(tvals_perm_FO ~= 0, 1);   % should be all false
frac_survive_FO  = (sum(perm_survives_FO)/perm_N)*100;    


%%
shuffle_state = 6;

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFC_state' num2str(shuffle_state)]);

% load the data
load(fullfile(outputShuffle, 'tvals.mat'));
load(fullfile(outputShuffle, 'mean_vals.mat'));
tvals_perm_FC_6 = tvals_perm;
perm_N = size(tvals_perm_FC_6,2);
tvals_perm_FC_6(tvals_perm_FC_6 < thresholds) = 0; % remove all values of t < thr
perm_survives_FC_6 = any(tvals_perm_FC_6 ~= 0, 1);   % [1 x perm_N]
frac_survive_FC_6  = sum(perm_survives_FC_6)/perm_N;    

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFO_state' num2str(shuffle_state)]);

% load the data
load(fullfile(outputShuffle, 'tvals.mat'));
tvals_perm_FO_6 = tvals_perm;
perm_N = size(tvals_perm_FO_6,2);
load(fullfile(outputShuffle, 'mean_vals.mat'));

tvals_perm_FO_6(tvals_perm_FO_6 < thresholds) = 0; % remove all values of t < thr
perm_survives_FO_6 = any(tvals_perm_FO_6 ~= 0, 1);   % should be all false
frac_survive_FO_6  = (sum(perm_survives_FO_6)/perm_N)*100;    


%%
shuffle_state = 7;

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFC_state' num2str(shuffle_state)]);

% load the data
load(fullfile(outputShuffle, 'tvals.mat'));
load(fullfile(outputShuffle, 'mean_vals.mat'));
tvals_perm_FC_7 = tvals_perm;
perm_N = size(tvals_perm_FC_7,2);
tvals_perm_FC_7(tvals_perm_FC_7 < thresholds) = 0; % remove all values of t < thr
perm_survives_FC_7 = any(tvals_perm_FC_7 ~= 0, 1);   % [1 x perm_N]
frac_survive_FC_7  = sum(perm_survives_FC_7)/perm_N;    

outputShuffle = fullfile(outputStatDir, ...
    ['Shuffle_GroupSwapFO_state' num2str(shuffle_state)]);

% load the data
load(fullfile(outputShuffle, 'tvals.mat'));
tvals_perm_FO_7 = tvals_perm;
perm_N = size(tvals_perm_FO_7,2);
load(fullfile(outputShuffle, 'mean_vals.mat'));

tvals_perm_FO_7(tvals_perm_FO_7 < thresholds) = 0; % remove all values of t < thr
perm_survives_FO_7 = any(tvals_perm_FO_7 ~= 0, 1);   % should be all false
frac_survive_FO_7  = (sum(perm_survives_FO_7)/perm_N)*100;    




%% Plot

figure('Color','w','Position',[300 300 400 550]);
h = bar([frac_survive_FO, frac_survive_FC, frac_survive_FO_6, frac_survive_FC_6, ...
    frac_survive_FO_7, frac_survive_FC_7]);

% Set colors: State 1 = red, State 6 = blue
h.FaceColor = 'flat';
h.CData(1,:) = [0.8500 0.3250 0.0980];  % FO State 1 -> red
h.CData(2,:) = [0.8500 0.3250 0.0980];  % FC State 1 -> red
h.CData(3,:) = [0.8500 0.3250 0.0980];  % FO State 6 -> red
h.CData(4,:) = [0.8500 0.3250 0.0980];  % FC State 6 -> red
h.CData(5,:) = [0.8500 0.3250 0.0980];  % FO State 6 -> red
h.CData(6,:) = [0.8500 0.3250 0.0980];  % FC State 6 -> red

set(gca,'XTickLabel',{'FO','FC'}, ...
        'FontSize',16,'LineWidth',2);
ylabel('% of permutations with Sigifcant results','FontSize',16);
ylim([0 1]);
box off;

saveas(gca, fullfile(FigureDir, 'Results_FO_FC_Shuffling_AllStates.jpg'))


% HeatMap

% Arrange data: rows = states, columns = FO and FC
data = [
    frac_survive_FO,   frac_survive_FC;
    frac_survive_FO_6, frac_survive_FC_6;
    frac_survive_FO_7, frac_survive_FC_7
];

figure('Color','w','Position',[300 300 450 320]);

h = heatmap({'FO','FC'}, ...
            {'State 1','State 6','State 7'}, data);

% Labels and appearance
h.Title = 'Proportion of permutations with significant results';
h.XLabel = '';
h.YLabel = '';
h.FontSize = 36;
h.CellLabelFormat = '%.2f';

% Common scale
h.ColorLimits = [0 1];

% Purple sequential colormap
n = 256;
cmap = interp1([0 0.5 1], ...
    [0.98 0.97 0.99; ...
     0.72 0.60 0.82; ...
     0.35 0.18 0.50], ...
    linspace(0,1,n));

h.Colormap = cmap;
h.ColorbarVisible = 'on';

saveas(gca, fullfile(FigureDir, 'Results_FO_FC_Shuffling_AllStates_Heatmap.jpg'))
