
% Function to get the first valid repetition directory for a given state
function [firstRepetitionDir, net_mean] = getFirstRepetition(resultsDir, states)
    % Create the path to the NumStates_XX directory
    numStatesDir = fullfile(resultsDir, ['NumStates_' num2str(states)]);
    
    % List all directories in the NumStates_XX directory
    repetitionsDirs = dir(numStatesDir);
    repetitionsDirs = repetitionsDirs([repetitionsDirs.isdir]); % Filter out non-directory entries
    repetitionsDirs = repetitionsDirs(~ismember({repetitionsDirs.name}, {'.', '..'})); % Exclude '.' and '..'
    
    % Get the first valid repetition directory
    if ~isempty(repetitionsDirs)
        firstRepetitionDir = repetitionsDirs(1).name; % Get the name of the first repetition directory
        net_mean = fullfile(numStatesDir, firstRepetitionDir, 'MeanActivation.mat');
    else
        error(['No Repetitions_XX directory found for state ' num2str(states)]);
    end
end
