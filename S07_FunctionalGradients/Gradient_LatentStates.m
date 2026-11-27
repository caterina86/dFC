%% Gradient analysis

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
reduced = 1; % 1 -> reduced number of volumes; 0 -> upsampling of dataset 2017

% number of volumes per subject
if reduced == 1
    T = 464; 
else
    T = 480;
end

if reduced == 1
    if pca == 0
        resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], ['zeromean' num2str(zeromean)]);
    else
        resultsDir = fullfile(DirData, ['PCA_REDUCED_' num2str(max(K))], num2str(pca), ['zeromean' num2str(zeromean)]);
    end        
else
    if pca == 0
        resultsDir = fullfile(DirData, ['PCA_UP_' num2str(max(K))], ['zeromean' num2str(zeromean)]);
    else
        resultsDir = fullfile(DirData, ['PCA_UP_' num2str(max(K))], num2str(pca),['zeromean' num2str(zeromean)]);
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
load(fullfile(resultsDir, ['HMMs_r' num2str(rep_i) '_SUBJECTS.mat']));
mkdir(fullfile(resultsDir, 'figures/'));
    
state_i = K(k); 
outputDir = fullfile(resultsDir, 'figures/', ['NumStates_' num2str(state_i)], ['Repetition_' num2str(rep_i)]);
mkdir(outputDir);



%% Load the reference gradients based on Margulies (2016)

nROIs = 100;

atlasFile = '/Users/caterina/Documents/Projects/dFC/temporal_filtering/Atlas/MNI/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii.gz';
V = double(niftiread(atlasFile));  % parcel labels 1..100

gradFiles = { ...
    fullfile(outputDir, 'Gradients', 'volume.all.0.nii.gz'), ... % Gradient 1
    fullfile(outputDir, 'Gradients', 'volume.all.1.nii.gz')};    % Gradient 2

Gref_parcel = zeros(nROIs, 2);

for g = 1:2
    Gvol = double(niftiread(gradFiles{g}));
    
    for p = 1:nROIs
        mask = (V == p) & ~isnan(Gvol) & Gvol ~= 0;
        Gref_parcel(p,g) = mean(Gvol(mask));
    end
end




%% -------------------------------
% Project the activation mean maps on the low-dimensional space of the gradients
% -------------------------------
nROIs = 100;         % number of ROIs
kStates = k+1;       % number of HMM states
nGrad = 2;           % gradients to use (1 and 2)
patients_idx = 34:83;  % adjust according to your subject indices
controls_idx = 1:33;

for subj = 1:N
    % Get mean activation for state s, subject subj
    meanAct = getMean(HMMs_dualregr{k, subj});  % nROIs x 1
    StateCoord_all(:, :, subj) = (Gref_parcel' * meanAct)/100;
end

StateCoord_controls = mean(StateCoord_all(:,:,controls_idx),3); % nGrad x K
StateCoord_patients = mean(StateCoord_all(:,:,patients_idx),3);
StateCoord_all = mean(StateCoord_all(:,:,:),3);

% Flatten and correlation
vec_controls = reshape(StateCoord_controls, [], 1);
vec_patients = reshape(StateCoord_patients, [], 1);
r_all = corr(vec_controls, vec_patients);
% Euclidean distance per state
dist_per_state = sqrt(sum((StateCoord_controls - StateCoord_patients).^2,1));

disp(['Pearson r across all states: ', num2str(r_all)]);
disp('Euclidean distance per state:');
disp(dist_per_state);




figure('Color', 'w');  hold on;% white background
% Define markers for groups
marker_controls = 'o';   % circle
marker_patients = '^';   % triangle

% Define colors for states
state_colors = [...
    1 0 0;      % state 1 - red
    1 1 0;      % state 2 - yellow
    0.5 0 0.5;  % state 3 - purple
    0 0 1;      % state 4 - blue
    0 0 1;      % state 5 - blue
    0 0 1;      % state 6 - blue
    0.5 0 0.5]; % state 7 - purple


for s = 1:kStates
    % Plot line connecting states
    plot([StateCoord_controls(1,s), StateCoord_patients(1,s)], ...
         [StateCoord_controls(2,s), StateCoord_patients(2,s)], 'k--');
    
    % Plot controls
    scatter(StateCoord_controls(1,s), StateCoord_controls(2,s), 200, ...
        'Marker', marker_controls, 'MarkerFaceColor', state_colors(s,:), 'MarkerEdgeColor','k');
    
    % Plot patients
    scatter(StateCoord_patients(1,s), StateCoord_patients(2,s), 200, ...
        'Marker', marker_patients, 'MarkerFaceColor', state_colors(s,:), 'MarkerEdgeColor','k');
    
    % Add text labels
    % text(StateCoord_controls(1,s)+0.05, StateCoord_controls(2,s)+0.05, num2str(s), 'Color','k');
    % text(StateCoord_patients(1,s)+0.05, StateCoord_patients(2,s)+0.05, num2str(s), 'Color','k');
end

xlabel('Gradient 1', 'FontSize', 20);
ylabel('Gradient 2', 'FontSize', 20);
%title('HMM states in gradient space');
set(gca, 'LineWidth', 2, 'FontSize', 22, 'Box','off', 'TickDir','out');
axis equal;

saveas(gca,fullfile(outputDir, 'Gradients/LowDimensionalManifold.jpg'))



%% Surrogate data
nPerm = 1000;

derivativesDir = fullfile(projectDir, 'derivatives/A1');
DirOut = fullfile(derivativesDir, 'HMM/PCA_REDUCED_15/10/Surrogates');

load(fullfile(DirOut, 'surrogate_hmm.mat'))

% Mean activation from the surrogate data
net_mean_surr = cell(nPerm,1);
for perm_i = 1:nPerm
    net_mean_surr{perm_i} = StateMeans_surr{perm_i}';   % convert to nROIs × nStates
end
save(fullfile(DirOut,'surrogate_net_mean.mat'),'net_mean_surr','-v7.3');


% 1. Latent states distance (1 measure)
surr_meanDist = zeros(nPerm,1);
for perm_i = 1:nPerm
    StateCoord_surr = (Gref_parcel' * net_mean_surr{perm_i}') ./ size(Gref_parcel,1);
    coords2D_surr = StateCoord_surr(1:2,:);  
    D = squareform(pdist(coords2D_surr'));      % 7×7
    surr_meanDist(perm_i) = mean(D(triu(true(size(D)),1)));
end

% --- EMPIRICAL MEAN DISTANCE ---
D_emp = squareform(pdist(StateCoord_controls(1:2,:)'));
emp_meanDist = mean(D_emp(triu(true(size(D_emp)),1)));

% --- TWO-TAILED PERMUTATION TEST ---
p_meanDist = mean(surr_meanDist >= emp_meanDist | surr_meanDist <= emp_meanDist);
fprintf('Empirical mean distance: %.3f, p-value: %.4f\n', emp_meanDist, p_meanDist);


% 2. Gradient extremity analysis
surr_gradpos = zeros(nPerm,2,7);
for perm_i = 1:nPerm
    StateCoord_surr = (Gref_parcel' * net_mean_surr{perm_i}') ./ size(Gref_parcel,1);
    surr_gradpos(perm_i,:,:) = StateCoord_surr(1:2,:); % (perm × grad × state)
end

% State-wise gradient significance (two axes × nStates)
p_grad = zeros(2,7);
for g = 1:2
    for s = 1:7
        emp_val = StateCoord_controls(g,s);
        null_vals = squeeze(surr_gradpos(:,g,s));
        
        p_grad(g,s) = mean(abs(null_vals - mean(null_vals)) >= ...
                          abs(emp_val - mean(null_vals)));
    end
end

p_FDR = mafdr(p_grad(:), 'BHFDR', true);
disp('Gradient extremity analysis');
p_grad_FDR = reshape(p_FDR , size(p_grad));

% significant values = position of the latent state along the gradient
% significantly different than the position it would assume by chance. 






