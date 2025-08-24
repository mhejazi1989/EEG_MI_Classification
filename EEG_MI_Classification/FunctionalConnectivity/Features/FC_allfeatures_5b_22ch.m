clear all;
close all;
clc;
warning off;

% Define the total number of subjects
totalSubjects = 29;

% List of subjects to exclude
excludeSubjects = [6, 8, 10, 13, 21, 23, 25, 27];

% Turn off the specific warning
warning('off', 'MATLAB:colon:operandsNotInteger');
fbands = 5;
Fs = 256;
for numSubjects = 1 :totalSubjects
    % Skip excluded subjects
    if any(numSubjects == excludeSubjects)
        continue;
    end

    % Start measuring time
    tic;
    subjectFolder = sprintf('Sub-%03d', numSubjects);
    mypath = fullfile('/lustre03/project/6067835/mhejazi/EEG-DATA', subjectFolder);

    % Load ALLEEG from the saved .mat file
    load(fullfile(mypath, sprintf('categorized_5b_values_sub%03d.mat', numSubjects)))
    % Assign your variables to fields in the struct
    conditions = {'REST', 'ME_DB', 'MI_DB', 'ME_FT', 'MI_FT'};

    selected_conditions = {categorized_band_values.BandValues_REST, ...
        categorized_band_values.BandValues_ME_DB, ...
        categorized_band_values.BandValues_MI_DB, ...
        categorized_band_values.BandValues_ME_FT, ...
        categorized_band_values.BandValues_MI_FT};

    % Initialize a cell array to hold the concatenated data for each frequency band
    numChannels2 = 22;
    numTrials = size(categorized_band_values.BandValues_ME_DB.PrefrontalRegion.Right{1,1},1);

    % Loop through each subject (assuming you have a numSubjects variable)
    for condIdx = 1:length(conditions)
        selected_condition = selected_conditions{condIdx}; % Current condition
        concatenated = cell(1, fbands); % Reset for each condition

        for bandIdx = 1:fbands
            % Extract data from each region and hemisphere for the current band
            prefrontalLeft = selected_condition.PrefrontalRegion.Left{bandIdx};    % Trials x 4 x 1025
            prefrontalRight = selected_condition.PrefrontalRegion.Right{bandIdx};  % Trials x 4 x 1025
            premotorLeft = selected_condition.PremotorCortex.Left{bandIdx};        % Trials x 3 x 1025
            premotorRight = selected_condition.PremotorCortex.Right{bandIdx};      % Trials x 3 x 1025
            motorLeft = selected_condition.MotorCortex.Left{bandIdx};              % Trials x 4 x 1025
            motorRight = selected_condition.MotorCortex.Right{bandIdx};            % Trials x 4 x 1025

            % Concatenate across the second dimension (channels) to form Trials x 14 x 1025
            concatenated{bandIdx} = cat(2, prefrontalLeft, prefrontalRight, ...
                premotorLeft, premotorRight, ...
                motorLeft, motorRight);
        end

        % Initialize a structure to hold the correlation ,coherence and phase-locking values for each band
        correlationMatrices_condition = cell(1, fbands); % Reset for each condition
        coherenceMatrices_condition = cell(1, fbands);
        phaseLockingValues_condition = cell(1, fbands);

        % Loop through each frequency band
        for bandIdx = 1:fbands
            % Get the data for the current band
            data = concatenated{bandIdx}; % Size: Trials x 22 x 1025

            % Initialize a 3D matrix to hold correlation,coherence and phase-locking values for each trial
            correlationMatrices_condition{bandIdx} = zeros(numChannels2, numChannels2, numTrials); % Size: 22 x 22 x Trials
            coherenceMatrices_condition{bandIdx} = zeros(numChannels2, numChannels2, numTrials); % Size: 22 x 22 x Trials
            phaseLockingValues_condition{bandIdx} = zeros(numChannels2, numChannels2, numTrials); % Size: 22 x 22 x Trials

            % Loop through each trial
            for trialIdx = 1:numTrials
                % Extract the data for the current trial
                trialData = squeeze(data(trialIdx, :, :)); % Size: 22 x 1025

                % Calculate the correlation matrix for the current trial
                correlationMatrices_condition{bandIdx}(:, :, trialIdx) = corr(trialData'); % Transpose to have channels as rows

                % Calculate coherence and phase-locking values
                for ch1 = 1:numChannels2
                    for ch2 = ch1+1:numChannels2 % Only compute for upper triangle to save computation

                        % Coherence calculation between channels ch1 and ch2
                        coherenceMatrices_condition{bandIdx}(ch1, ch2, trialIdx) = mean(mscohere(trialData(ch1, :), trialData(ch2, :), [], [], [], Fs));

                        % Phase-locking value (PLV) calculation
                        % Compute the instantaneous phase using the Hilbert transform
                         % Extract the instantaneous phases
                        phaseCh1 = angle(hilbert(trialData(ch1, :)));
                        phaseCh2 = angle(hilbert(trialData(ch2, :)));

                        % Calculate PLV between channels ch1 and ch2
                        phaseLockingValues_condition{bandIdx}(ch1, ch2, trialIdx) = abs(mean(exp(1i * (phaseCh1 - phaseCh2))));
                    end
                end
            end
        end

        % Initialize cell arrays to store upper triangular values for each band and each trial
        correlationUpper_condition = cell(1, fbands);
        coherenceUpper_condition = cell(1, fbands);
        phaseLockingUpper_condition = cell(1, fbands);

        for bandIdx = 1:fbands
            numTrials = size(correlationMatrices_condition{bandIdx}, 3);
            numChannels = size(correlationMatrices_condition{bandIdx}, 1);

            % Pre-allocate matrices to store upper triangular values (excluding diagonal)
            upperTriSize = (numChannels * (numChannels - 1)) / 2; % Number of upper triangle elements excluding diagonal
            correlationUpper_condition{bandIdx} = zeros(numTrials, upperTriSize);
            coherenceUpper_condition{bandIdx} = zeros(numTrials, upperTriSize);
            phaseLockingUpper_condition{bandIdx} = zeros(numTrials, upperTriSize);

            % Logical index for the upper triangle (ignoring diagonal)
            upperTriIdx = find(triu(ones(numChannels), 1));

            for trialIdx = 1:numTrials
                % Extract the upper triangle (excluding diagonal) for correlation
                corrMatrix = correlationMatrices_condition{bandIdx}(:, :, trialIdx);
                correlationUpper_condition{bandIdx}(trialIdx,:) = corrMatrix(upperTriIdx);

                % Extract the upper triangle (excluding diagonal) for coherence
                coherenceMatrix = coherenceMatrices_condition{bandIdx}(:, :, trialIdx);
                coherenceUpper_condition{bandIdx}(trialIdx,:) = coherenceMatrix(upperTriIdx);

                % Extract the upper triangle (excluding diagonal) for phase-locking values
                plvMatrix = phaseLockingValues_condition{bandIdx}(:, :, trialIdx);
                phaseLockingUpper_condition{bandIdx}(trialIdx,:) = plvMatrix(upperTriIdx)';
            end
        end

        % Save the upper triangular values for each band and each trial
        save(fullfile(mypath, sprintf('upper_triangle_connectivity_5b_22ch_sub%03d_%s.mat', numSubjects, conditions{condIdx})), ...
           'correlationUpper_condition', 'coherenceUpper_condition', 'phaseLockingUpper_condition');
    end

    toc;  % End measuring time for this subject
end
% Stop profiling
profile off;

% Restore the figure visibility setting
set(0, 'DefaultFigureVisible', 'on');

% Suppress warning
warning off;
