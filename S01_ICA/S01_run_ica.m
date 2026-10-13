% Caterina Pedersini, October 2024, cpedersini86@gmail.com

% Script to Perform single-subject ICA denoise

clearvars; 
close all;


% set paths
baseDir = uigetdir(pwd, 'Select the base directory');
projectDir=fullfile(baseDir, 'osf');

subjects = importdata(fullfile(projectDir, 'list_subjects_dataset2023_2017.txt'));
path_fmri = fullfile(baseDir, 'derivatives/fmriprep');
path_output = fullfile(baseDir, 'derivatives/ica');
task = 'rest'; % define the task name
bet_f = 0.4;
mkdir(path_output)



%% bet, registration and ICA

disp('registration and ICA');

for sub_i = 1:length(subjects)
        
    % set subject name
    subject = ['sub-' subjects{sub_i}];
    disp(subject);
    
    try 
    
        % Organize the data:
        mkdir(fullfile(path_output, subject));
        
        % copy and rename the functional and structural preprocessed image from the fmriprep output directory
        system(['cp ' fullfile(path_fmri, subject, 'func', [subject '_task-' task '_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz']) ' ' ...
            fullfile(path_output, subject, [subject '_func.nii.gz'])])
        system(['cp ' fullfile(path_fmri, subject, 'anat', [subject '_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz']) ' ' ...
            fullfile(path_output, subject, [subject '_T1w.nii.gz'])])                
        
        
        if isfolder(fullfile(path_output, subject))

            if ~isfile(fullfile(path_output, subject, [subject '_filtered_func.ica'], 'mean.nii.gz')) && ~isfile(fullfile(path_output, subject, [subject '_func_bet.ica'], 'mean.nii.gz'))
                            
                disp(['Prepare data of subject '  subject])
                
                % load data
                functional_raw = [path_output '/' subject];
                functional = [functional_raw '/' [subject '_func.nii.gz']]; % fMRI data in .nii format
                gunzip(functional, functional_raw)
                system(['rm ' functional])

                % bet on functional 4D
                functional_unzip = fullfile(path_output, subject, [subject '_func.nii']); % fMRI data in .nii format
                functional_bet = fullfile(path_output, subject, [subject '_func_bet.nii']);
               
                if ~isfile(functional_bet)
                    
                    system(['bet ' functional_unzip ' ' functional_bet ' -f ' num2str(bet_f) ' -R -S -F -v']);
                    % erode functional
                    system(['fslmaths ' functional_bet ' -ero ' functional_bet]);
                    
                end

                % bet on structural t1
                t1 = fullfile(path_output, subject, [subject '_T1w.nii.gz']);
                t1_bet3 = fullfile(path_output, subject, [subject '_T1w_bet.nii.gz']);

                system(['bet ' fullfile(path_output, subject, [subject '_T1w.nii']) ' ' t1_bet3 ' -f ' num2str(bet_f) ' -R -S -v']);

                
                %% ICA melodic

                disp(['Run Melodic ICA for subject '  subject])
                system(['melodic -i ' functional_bet ' --nobet --report --vn --bgthreshold=10 --Oall -v']);

            else
                
                disp(['ICA for subject '  subject ' already done'])
                
            end
            
            
            if ~isfile(fullfile(path_output, subject, [subject '_filtered_func.ica'], 'mean_func.nii.gz'))

                            
                % rename folder
                system(['mv ' fullfile(path_output, subject, [subject '_func_bet.ica/']) ' ' fullfile(path_output, subject, [subject '_filtered_func.ica/'])])
                system(['cp ' fullfile(path_output, subject, [subject '_func_bet.nii.gz']) ' ' fullfile(path_output, subject, [subject '_filtered_bet.nii.gz'])])
                system(['cp ' fullfile(path_output, subject, [subject '_filtered_bet.nii.gz']) ' ' fullfile(path_output, subject, [subject '_filtered_func.ica'], [subject '_filtered_func.nii.gz'])])
                
                % mask and mean files                 
                system(['cp ' fullfile(path_output, subject, [subject '_filtered_func.ica'], 'mean.nii.gz') ' ' fullfile(path_output, subject, [subject '_filtered_func.ica'], 'mean_func.nii.gz')]);
            
            end
        end
        
    catch
        
    end

end


%% Visualize melodic-ICA

%subj="sub-ctrl18"

%echo subject ${subj}

%fsleyes --scene melodic -ad ${baseDir}/derivatives/ica/${subj}/${subj}_filtered_func.ica/