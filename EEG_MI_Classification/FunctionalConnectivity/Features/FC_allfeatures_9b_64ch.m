clear all;
close all;
clc;
warning off;

%% Define Parameters
% Total number of subjects and excluded subjects
totalSubjects = 29;
%excludeSubjects = [6, 8, 10, 13, 21, 23, 25, 27];
warning('off', 'MATLAB:colon:operandsNotInteger');

% Frequency bands definition
bands = struct('name', {'band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'band7', 'band8', 'band9'}, ...
    'min', {4, 8, 12, 16, 20, 24, 28, 32, 36}, ...
    'max', {8, 12, 16, 20, 24, 28, 32, 36, 40});
numChannels = 64;
Fs = 256;

%% Loop Through Subjects
%for numSubjects = 1:totalSubjects
excludeSubjects = [18,19,21,22,23];

for numSubjects = 24

    if any(numSubjects == excludeSubjects)
        continue;
    end

    tic; % Start measuring time
    subjectID = sprintf('%03d', numSubjects);
    subjectFolder = sprintf('Sub-%s', subjectID);
    mypath = fullfile('/lustre03/project/6067835/mhejazi/EEG-DATA', subjectFolder);

    % Load subject data
    load(fullfile(mypath, sprintf('sub%s_9BandValues.mat', subjectID)), ...
        'BandValues_ME_DB', 'BandValues_MI_DB', 'BandValues_ME_FT', 'BandValues_MI_FT', 'BandValues_REST');

    % Task and Band Matrices
    taskNames = {'ME_DB', 'MI_DB', 'ME_FT', 'MI_FT', 'REST'};
    BandMatrices = {BandValues_ME_DB, BandValues_MI_DB, BandValues_ME_FT, BandValues_MI_FT, BandValues_REST};
    minAvailableTrials = size(BandValues_ME_DB{1}, 1);

    % Initialize Connectivity and Feature Matrices
    ConnectivityPLV = cell(numel(taskNames), numel(bands));
    ConnectivityCOH = cell(numel(taskNames), numel(bands));
    ConnectivityCORR = cell(numel(taskNames), numel(bands));
    Feature_MatrixPLV = cell(numel(taskNames), numel(bands));
    Feature_MatrixCOH = cell(numel(taskNames), numel(bands));
    Feature_MatrixCORR = cell(numel(taskNames), numel(bands));

    %% Compute Connectivity for Each Task and Band
    for taskIdx = 1:numel(taskNames)
        for bandIdx = 1:numel(bands)
            currentBandData = BandMatrices{taskIdx}{bandIdx};

            % Initialize Connectivity Matrices
            plvMatrix = zeros(numChannels, numChannels, minAvailableTrials);
            cohMatrix = zeros(numChannels, numChannels, minAvailableTrials);
            corrMatrix = zeros(numChannels, numChannels, minAvailableTrials);

            % Compute connectivity measures for each trial
            for trial = 1:minAvailableTrials
                trialData = squeeze(currentBandData(trial, :, :));
                normalizedTrialData = (trialData - mean(trialData, 2)) ./ std(trialData, 0, 2);

                for ch1 = 1:numChannels
                    for ch2 = ch1+1:numChannels
                        phaseCh1 = angle(hilbert(normalizedTrialData(ch1, :)));
                        phaseCh2 = angle(hilbert(normalizedTrialData(ch2, :)));
                        plvMatrix(ch1, ch2, trial) = abs(mean(exp(1i * (phaseCh1 - phaseCh2))));

                        % Coherence Calculation
                        [Cxy, ~] = mscohere(normalizedTrialData(ch1, :), normalizedTrialData(ch2, :), [], [], [], Fs);
                        cohMatrix(ch1, ch2, trial) = mean(Cxy);

                        % Correlation Calculation
                        corrMatrix(ch1, ch2, trial) = corr(normalizedTrialData(ch1, :)', normalizedTrialData(ch2, :)');
                    end
                end
            end

            % Store Connectivity Data
            ConnectivityPLV{taskIdx, bandIdx} = plvMatrix;
            ConnectivityCOH{taskIdx, bandIdx} = cohMatrix;
            ConnectivityCORR{taskIdx, bandIdx} = corrMatrix;

            % Extract Features
            numConns = (numChannels * (numChannels - 1)) / 2;
            FeatureMatrixPLV = zeros(minAvailableTrials, numConns);
            FeatureMatrixCOH = zeros(minAvailableTrials, numConns);
            FeatureMatrixCORR = zeros(minAvailableTrials, numConns);

            for trial = 1:minAvailableTrials
                PLVMatrix = squeeze(plvMatrix(:, :, trial));
                COHMatrix = squeeze(cohMatrix(:, :, trial));
                CORRMatrix = squeeze(corrMatrix(:, :, trial));

                plvUpperTri = PLVMatrix(triu(true(numChannels), 1));
                cohUpperTri = COHMatrix(triu(true(numChannels), 1));
                corrUpperTri = CORRMatrix(triu(true(numChannels), 1));

                FeatureMatrixPLV(trial, :) = plvUpperTri;
                FeatureMatrixCOH(trial, :) = cohUpperTri;
                FeatureMatrixCORR(trial, :) = corrUpperTri;

            end

            Feature_MatrixPLV{taskIdx, bandIdx} = FeatureMatrixPLV;
            Feature_MatrixCOH{taskIdx, bandIdx} = FeatureMatrixCOH;
            Feature_MatrixCORR{taskIdx, bandIdx} = FeatureMatrixCORR;
        end
    end

% Initialize variables to store concatenated feature matrices for each condition
metrics = {'plv', 'coh', 'corr'};
AllFeatureMatrices = struct();

% Loop through each condition and each metric, then concatenate the feature matrices
for cIdx = 1:numel(taskNames)
    for mIdx = 1:numel(metrics)
        % Initialize the corresponding feature matrix as empty
        fieldName = sprintf('AllFeatureMatrix_%s_%s', taskNames{cIdx}, metrics{mIdx});
        AllFeatureMatrices.(fieldName) = [];

        % Concatenate for each frequency band
        for bandIdx = 1:numel(bands)
            matrixField = sprintf('Feature_Matrix%s{cIdx, bandIdx}', upper(metrics{mIdx}));
            AllFeatureMatrices.(fieldName) = [AllFeatureMatrices.(fieldName), eval(matrixField)];
        end
    end
end

% Save all feature matrices in one file for the current subject
savePath = fullfile('/lustre03/project/6067835/mhejazi/EEG-DATA', subjectFolder);
save(fullfile(savePath, sprintf('sub%03d_9b_AllFeatureMatrix.mat', numSubjects)), ...
    '-struct', 'AllFeatureMatrices');
   % Display completion message, elapsed time, and memory usage
    elapsedTime = toc; % Stop the timer
    fprintf('Subject %d completed in %.2f seconds.\n', numSubjects, elapsedTime);
end
