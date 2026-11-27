%% FC comparing controls vs patients

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

sizes = {'Extent'}; %{'Intensity'}; %{'Extent'};


%% Run NBS for the dynamic FC

% options -> t-test; permutations = 10000; significance = 0.05; method = NBS; 
% Components size = Extent

% Create the contrast file
contrast = [-1 1 0 0 0];
contrastFile1 = fullfile(outputDir, 'NBS/contrast_EBvsSC.txt');
writematrix(contrast, contrastFile1, 'Delimiter','tab');  % tab-delimited

contrast = [1 -1 0 0 0];
contrastFile2 = fullfile(outputDir, 'NBS/contrast_SCvsEB.txt');
writematrix(contrast, contrastFile2, 'Delimiter','tab');  % tab-delimited

thresholds = 3.5;

% Number of subjects
nSubjects = 83;
% Create S struct
S = struct();
S.Nsub = nSubjects;
S.OUT.waitbar = [];   % optional

num_perm = '10000';




%% Whole brain FC for each state - Dynamic FC

outputDynDir = fullfile(outputDir, 'NBS', sizes{1}, 'dynamicFC');
mkdir(outputDynDir);

for k_i = 1:(k+1)
    
    for sub_i = 1:N
        [~,omega(:,:,sub_i)] = getFuncConn(HMMs_dualregr{k,sub_i},k_i);
        subject = omega(:,:,sub_i);
        subject_z = atanh(subject); % from r to z-fisher       
        matrices(:,:,sub_i) = subject_z;     
    end
    
    save(fullfile(outputDynDir, ['matrices_state' num2str(k_i) '.mat']), 'matrices');

end

for k_i = 1:7
    
    clear global nbs
    
    for cont_i = 1:2
    
        if cont_i == 1
            contrastFile = contrastFile1;
            label = 'EBvsSC';
        else
            contrastFile = contrastFile2;
            label = 'SCvsEB';
        end
    
        % settings nbs 
        UI = [];
        UI.test.ui = 't-test';
        UI.size.ui=sizes{1};
        UI.design.ui = fullfile(outputDir, 'NBS/designMatrix.mat');
        UI.matrices.ui = fullfile(outputDynDir, ['matrices_state' num2str(k_i) '.mat']);
        UI.contrast.ui = contrastFile;
        UI.exchange.ui = '';
        UI.node_label.ui = fullfile(outputDir, 'NBS/nodeLabels.txt');
        UI.node_coor.ui = fullfile(outputDir, 'NBS/COG.txt');
        UI.perms.ui = num_perm;
        UI.method.ui = 'Run NBS';
        UI.thresh.ui = num2str(thresholds(1));
        UI.alpha.ui = '0.05';
        
        % --- Run NBS ---
        NBSrun(UI, S);
        global nbs
        nbs_outputs = nbs;

        % Pause briefly to ensure figure is created
        pause(1);  % adjust if needed

        % Find all figure windows
        figHandles = findall(0,'Type','figure');

        for f = 1:length(figHandles)
            fig = figHandles(f);
            % Save only real NBS result figures
            saveas(fig, fullfile(outputDynDir, ...
                ['Results_State' num2str(k_i) '_' label '.png']));
            % Close both figures AND waitbars
            close(fig);
        end
        
        disp('Figures saved and closed, continuing with the script.');
        close all

        save(fullfile(outputDynDir, ['Results_State' num2str(k_i) '_' label '.mat']), 'nbs_outputs');      
        clear nbs_outputs nbs;
        
    end
end



%% Time averaged Functional Connectome (stationary FC estimated by the model)

outputStatDir = fullfile(outputDir, 'NBS', sizes{1}, 'staticFC');
mkdir(outputStatDir);

for sub_i = 1:N

    % Extract state correlation matrix
    for k_i = 1:(k+1)
        [~, FC(:,:,k_i)] = getFuncConn(HMMs_dualregr{k,sub_i},k_i); 
    end

    num_states = size(FC, 3);
    weighted_sum = zeros(size(FC,1));
    FO = FOdual{rep_i,k}(sub_i, :);

    for s = 1:num_states
        weighted_sum = weighted_sum + FO(s) * FC(:,:,s);
    end

    time_averaged_FC{sub_i} = weighted_sum / sum(FO);
    matrices_static(:,:,sub_i) = atanh(time_averaged_FC{sub_i}); % from r to z-fisher       

end
save(fullfile(outputStatDir, 'matrices_TimeAveragedFC.mat'), 'matrices_static');


for cont_i = 1:2
    
    clear global nbs
        
    if cont_i == 1
        contrastFile = contrastFile1;
        label = 'EBvsSC';
    else
        contrastFile = contrastFile2;
        label = 'SCvsEB';
    end

    % settings nbs
    UI = [];
    UI.test.ui = 't-test';
    UI.size.ui=sizes{1};
    UI.matrices.ui = fullfile(outputStatDir, 'matrices_TimeAveragedFC.mat');
    UI.design.ui = fullfile(outputDir, 'NBS/designMatrix.mat');
    UI.contrast.ui = contrastFile;
    UI.exchange.ui = '';
    UI.node_label.ui = fullfile(outputDir, 'NBS/nodeLabels.txt');
    UI.node_coor.ui = fullfile(outputDir, 'NBS/COG.txt');
    UI.perms.ui = num_perm;
    UI.method.ui = 'Run NBS';
    UI.thresh.ui = num2str(thresholds(1));
    UI.alpha.ui = '0.05';

    % run nbs
    NBSrun(UI,S);   
    global nbs
    nbs_outputs = nbs;    
    
    % Pause briefly to ensure figure is created
    pause(1);  % adjust if needed

    % Find all figure windows
    figHandles = findall(0,'Type','figure');

    for f = 1:length(figHandles)
        fig = figHandles(f);
        % Save only real NBS result figures
        saveas(fig, fullfile(outputStatDir, ...
            ['Results_TimeAveragedFC_' label '.png']));
        % Close both figures AND waitbars
        close(fig);
    end

    disp('Figures saved and closed, continuing with the script.');
    close all  

    save(fullfile(outputStatDir, ['Results_TimeAveragedFC_' label '.mat']), 'nbs_outputs');    
    clear nbs;

end






%% Static Functional connectome (traditional way of calculation)

derivativesDir = fullfile(projectDir, 'derivatives/A1');
DirData = fullfile(derivativesDir, 'HMM');

for sub_i=1:N
    my_string = ['sub-' subjects{sub_i}];        
    f = fullfile(DirData, ['sub-' subjects{sub_i}], ['sub-' subjects{sub_i} '_data_Atlas_reduced.txt']);           
    timecourses = load(f); 
    
    % Extract the stationary Functional Connectome
    static_FC = corr(timecourses);    % returns 100 x 100 parcel-by-parcel FC matrix
    static_FCz(:,:,sub_i) = atanh(static_FC);
    save(fullfile(outputStatDir, 'matrices_staticFC.mat'), 'static_FCz');

end 


% Run NBS for the stationary FC (estimated from the empirical data)
    
for cont_i = 1:2
    
    clear global nbs

    if cont_i == 1
        contrastFile = contrastFile1;
        label = 'EBvsSC';
    else
        contrastFile = contrastFile2;
        label = 'SCvsEB';
    end

    % settings nbs
    UI = [];
    UI.test.ui = 't-test';
    UI.size.ui=sizes{1};
    UI.matrices.ui = fullfile(outputStatDir, 'matrices_staticFC.mat');
    UI.design.ui = fullfile(outputDir, 'NBS/designMatrix.mat');
    UI.contrast.ui = contrastFile;
    UI.exchange.ui = '';
    UI.node_label.ui = fullfile(outputDir, 'NBS/nodeLabels.txt');
    UI.node_coor.ui = fullfile(outputDir, 'NBS/COG.txt');
    UI.perms.ui = num_perm;
    UI.method.ui = 'Run NBS';
    UI.thresh.ui = num2str(thresholds(1));
    UI.alpha.ui = '0.05';

    % run nbs
    NBSrun(UI,S);
    
    global nbs
    nbs_outputs = nbs;    
    
    % Pause briefly to ensure figure is created
    pause(1);  % adjust if needed

    % Find all figure windows
    figHandles = findall(0,'Type','figure');

    for f = 1:length(figHandles)
        fig = figHandles(f);
        % Save only real NBS result figures
        saveas(fig, fullfile(outputStatDir, ...
            ['Results_staticFC_' label '.png']));
        % Close both figures AND waitbars
        close(fig);
    end

    disp('Figures saved and closed, continuing with the script.');
    close all  
    
    save(fullfile(outputStatDir, ['Results_staticFC_' label '.mat']), 'nbs_outputs');
    clear nbs nbs_outputs;
end



