%% Initialization
clear all;
close all;
clc;
warning off;
% Connec for 22 channels
% Define the total number of subjects
totalSubjects = 29;

% List of subjects to exclude
excludeSubjects = [6, 8, 10, 13, 21,22, 23,24, 25,26, 27,28,29];

classes = 1:6; % Define your classes
numIteration = 15;
K = 5; % Number of folds for cross-validation

fBands = {'Band1', 'Band2', 'Band3', 'Band4', 'Band5', 'Band6', 'Band7', 'Band8', 'Band9'};
MaxFeatureNum = 10*length(fBands);
classProbs = {[1, 2], [1, 3], [1, 4], [1, 5], [2, 3], [4, 5]}; % Define class pairs

% Define the number of features to evaluate
numFeaturesList = [10, 20, 30, 40, MaxFeatureNum];

localAccuracies_lda_corr = zeros(length(classProbs),K, numIteration, length(numFeaturesList));
localAccuracies_svm_corr = zeros(length(classProbs),K, numIteration, length(numFeaturesList));
localAccuracies_dt_corr = zeros(length(classProbs),K, numIteration, length(numFeaturesList));

localAccuracies_lda_coh = zeros(length(classProbs),K, numIteration, length(numFeaturesList));
localAccuracies_svm_coh = zeros(length(classProbs),K, numIteration, length(numFeaturesList));
localAccuracies_dt_coh = zeros(length(classProbs),K, numIteration, length(numFeaturesList));

localAccuracies_lda_plv = zeros(length(classProbs),K, numIteration, length(numFeaturesList));
localAccuracies_svm_plv = zeros(length(classProbs),K, numIteration, length(numFeaturesList));
localAccuracies_dt_plv = zeros(length(classProbs),K, numIteration, length(numFeaturesList));

% Initialize a counter for valid subjects
validSubjectsCount = 0;

% Assign your variables to fields in the struct
conditions = {'REST', 'ME_DB', 'MI_DB', 'ME_FT', 'MI_FT'};

for numSubjects = 19
    % Skip excluded subjects
    if any(numSubjects == excludeSubjects)
        continue;
    end
    % Increment the valid subjects count
    validSubjectsCount = validSubjectsCount + 1;
    % Start measuring time
    tic;

    % Define the subject folder path
    subjectFolder = sprintf('Sub-%03d', numSubjects);
    myPath = fullfile('/project/def-b09sdp/mhejazi/EEG-DATA', subjectFolder);

    % Initialize struct to store the band values for this subject
    BandValues = struct();

    % Loop through each condition
    for condIdx = 1:length(conditions)
        % Define the path to the .mat file for each condition
        matFilePath = fullfile(myPath, sprintf('upper_triangle_connectivity_9b_22ch_sub%03d_%s.mat', numSubjects, conditions{condIdx}));

        % Check if the file exists
        if exist(matFilePath, 'file') ~= 2
            fprintf('Error loading file for Subject %d: Unable to find file or directory "%s".\n', numSubjects, matFilePath);
            continue;  % Skip this condition if the file is missing
        end

        % Load the .mat file containing the correlation, coherence, and plv locking matrices
        load(matFilePath, 'correlationUpper_condition', 'coherenceUpper_condition', 'phaseLockingUpper_condition');

        % Check if 'correlationUpper_condition', 'coherenceUpper_condition', and 'phaseLockingUpper_condition' have 5 cells for each frequency band
        if length(correlationUpper_condition) ~= length(fBands) || ...
                length(coherenceUpper_condition) ~= length(fBands) || ...
                length(phaseLockingUpper_condition) ~= length(fBands)
            fprintf('Warning: data for Subject %d, Condition %s does not match expected number of frequency bands.\n', numSubjects, conditions{condIdx});
            continue;
        end

        % Assign each frequency band to its respective field in the struct for correlation
        for bandIdx = 1:length(fBands)
            BandValues.(conditions{condIdx}).(fBands{bandIdx}).correlation = correlationUpper_condition{bandIdx};
            BandValues.(conditions{condIdx}).(fBands{bandIdx}).coherence = coherenceUpper_condition{bandIdx};
            BandValues.(conditions{condIdx}).(fBands{bandIdx}).phaseLocking = phaseLockingUpper_condition{bandIdx};
        end
    end

    % Prepare band values for classification based on correlation
    BandValues_REST_corr = [BandValues.REST.Band1.correlation, BandValues.REST.Band2.correlation, BandValues.REST.Band3.correlation, BandValues.REST.Band4.correlation, BandValues.REST.Band5.correlation, BandValues.REST.Band6.correlation, BandValues.REST.Band7.correlation, BandValues.REST.Band8.correlation, BandValues.REST.Band9.correlation];
    BandValues_MI_FT_corr = [BandValues.MI_FT.Band1.correlation, BandValues.MI_FT.Band2.correlation, BandValues.MI_FT.Band3.correlation, BandValues.MI_FT.Band4.correlation, BandValues.MI_FT.Band5.correlation, BandValues.MI_FT.Band6.correlation, BandValues.MI_FT.Band7.correlation, BandValues.MI_FT.Band8.correlation, BandValues.MI_FT.Band9.correlation];
    BandValues_MI_DB_corr = [BandValues.MI_DB.Band1.correlation, BandValues.MI_DB.Band2.correlation, BandValues.MI_DB.Band3.correlation, BandValues.MI_DB.Band4.correlation, BandValues.MI_DB.Band5.correlation, BandValues.MI_DB.Band6.correlation, BandValues.MI_DB.Band7.correlation, BandValues.MI_DB.Band8.correlation, BandValues.MI_DB.Band9.correlation];
    BandValues_ME_DB_corr = [BandValues.ME_DB.Band1.correlation, BandValues.ME_DB.Band2.correlation, BandValues.ME_DB.Band3.correlation, BandValues.ME_DB.Band4.correlation, BandValues.ME_DB.Band5.correlation, BandValues.ME_DB.Band6.correlation, BandValues.ME_DB.Band7.correlation, BandValues.ME_DB.Band8.correlation, BandValues.ME_DB.Band9.correlation];
    BandValues_ME_FT_corr = [BandValues.ME_FT.Band1.correlation, BandValues.ME_FT.Band2.correlation, BandValues.ME_FT.Band3.correlation, BandValues.ME_FT.Band4.correlation, BandValues.ME_FT.Band5.correlation, BandValues.ME_FT.Band6.correlation, BandValues.ME_FT.Band7.correlation, BandValues.ME_FT.Band8.correlation, BandValues.ME_FT.Band9.correlation];

    % Prepare band values for classification based on coherence
    BandValues_REST_coh = [BandValues.REST.Band1.coherence, BandValues.REST.Band2.coherence, BandValues.REST.Band3.coherence, BandValues.REST.Band4.coherence, BandValues.REST.Band5.coherence, BandValues.REST.Band6.coherence, BandValues.REST.Band7.coherence, BandValues.REST.Band8.coherence, BandValues.REST.Band9.coherence];
    BandValues_MI_FT_coh = [BandValues.MI_FT.Band1.coherence, BandValues.MI_FT.Band2.coherence, BandValues.MI_FT.Band3.coherence, BandValues.MI_FT.Band4.coherence, BandValues.MI_FT.Band5.coherence, BandValues.MI_FT.Band6.coherence, BandValues.MI_FT.Band7.coherence, BandValues.MI_FT.Band8.coherence, BandValues.MI_FT.Band9.coherence];
    BandValues_MI_DB_coh = [BandValues.MI_DB.Band1.coherence, BandValues.MI_DB.Band2.coherence, BandValues.MI_DB.Band3.coherence, BandValues.MI_DB.Band4.coherence, BandValues.MI_DB.Band5.coherence, BandValues.MI_DB.Band6.coherence, BandValues.MI_DB.Band7.coherence, BandValues.MI_DB.Band8.coherence, BandValues.MI_DB.Band9.coherence];
    BandValues_ME_DB_coh = [BandValues.ME_DB.Band1.coherence, BandValues.ME_DB.Band2.coherence, BandValues.ME_DB.Band3.coherence, BandValues.ME_DB.Band4.coherence, BandValues.ME_DB.Band5.coherence, BandValues.ME_DB.Band6.coherence, BandValues.ME_DB.Band7.coherence, BandValues.ME_DB.Band8.coherence, BandValues.ME_DB.Band9.coherence];
    BandValues_ME_FT_coh = [BandValues.ME_FT.Band1.coherence, BandValues.ME_FT.Band2.coherence, BandValues.ME_FT.Band3.coherence, BandValues.ME_FT.Band4.coherence, BandValues.ME_FT.Band5.coherence, BandValues.ME_FT.Band6.coherence, BandValues.ME_FT.Band7.coherence, BandValues.ME_FT.Band8.coherence, BandValues.ME_FT.Band9.coherence];

    % Prepare band values for classification based on PLV phase locking
    BandValues_REST_plv = [BandValues.REST.Band1.phaseLocking, BandValues.REST.Band2.phaseLocking, BandValues.REST.Band3.phaseLocking, BandValues.REST.Band4.phaseLocking, BandValues.REST.Band5.phaseLocking, BandValues.REST.Band6.phaseLocking, BandValues.REST.Band7.phaseLocking, BandValues.REST.Band8.phaseLocking, BandValues.REST.Band9.phaseLocking];
    BandValues_MI_FT_plv = [BandValues.MI_FT.Band1.phaseLocking, BandValues.MI_FT.Band2.phaseLocking, BandValues.MI_FT.Band3.phaseLocking, BandValues.MI_FT.Band4.phaseLocking, BandValues.MI_FT.Band5.phaseLocking, BandValues.MI_FT.Band6.phaseLocking, BandValues.MI_FT.Band7.phaseLocking, BandValues.MI_FT.Band8.phaseLocking, BandValues.MI_FT.Band9.phaseLocking];
    BandValues_MI_DB_plv = [BandValues.MI_DB.Band1.phaseLocking, BandValues.MI_DB.Band2.phaseLocking, BandValues.MI_DB.Band3.phaseLocking, BandValues.MI_DB.Band4.phaseLocking, BandValues.MI_DB.Band5.phaseLocking, BandValues.MI_DB.Band6.phaseLocking, BandValues.MI_DB.Band7.phaseLocking, BandValues.MI_DB.Band8.phaseLocking, BandValues.MI_DB.Band9.phaseLocking];
    BandValues_ME_DB_plv = [BandValues.ME_DB.Band1.phaseLocking, BandValues.ME_DB.Band2.phaseLocking, BandValues.ME_DB.Band3.phaseLocking, BandValues.ME_DB.Band4.phaseLocking, BandValues.ME_DB.Band5.phaseLocking, BandValues.ME_DB.Band6.phaseLocking, BandValues.ME_DB.Band7.phaseLocking, BandValues.ME_DB.Band8.phaseLocking, BandValues.ME_DB.Band9.phaseLocking];
    BandValues_ME_FT_plv = [BandValues.ME_FT.Band1.phaseLocking, BandValues.ME_FT.Band2.phaseLocking, BandValues.ME_FT.Band3.phaseLocking, BandValues.ME_FT.Band4.phaseLocking, BandValues.ME_FT.Band5.phaseLocking, BandValues.ME_FT.Band6.phaseLocking, BandValues.ME_FT.Band7.phaseLocking, BandValues.ME_FT.Band8.phaseLocking, BandValues.ME_FT.Band9.phaseLocking];

    % Now you can continue with your classification code
    % Create the Labels cell array with matrices for correlation
    classLabels_corr = {BandValues_REST_corr, BandValues_MI_FT_corr, BandValues_MI_DB_corr, ...
        BandValues_ME_DB_corr, BandValues_ME_FT_corr};

    % Create the Labels cell array with matrices for coherence
    classLabels_coh = {BandValues_REST_coh, BandValues_MI_FT_coh, BandValues_MI_DB_coh, ...
        BandValues_ME_DB_coh, BandValues_ME_FT_coh};

    % Create the Labels cell array with matrices for plv locking
    classLabels_plv = {BandValues_REST_plv, BandValues_MI_FT_plv, BandValues_MI_DB_plv, ...
        BandValues_ME_DB_plv, BandValues_ME_FT_plv};

    % Initialize storage for accuracies for each category
    meanAccuraciesLDA_corr = zeros(length(classProbs), length(numFeaturesList));
    meanAccuraciesSVM_corr = zeros(length(classProbs), length(numFeaturesList));
    meanAccuraciesDT_corr = zeros(length(classProbs), length(numFeaturesList));

    meanAccuraciesLDA_coh = zeros(length(classProbs), length(numFeaturesList));
    meanAccuraciesSVM_coh = zeros(length(classProbs), length(numFeaturesList));
    meanAccuraciesDT_coh = zeros(length(classProbs), length(numFeaturesList));

    meanAccuraciesLDA_plv = zeros(length(classProbs), length(numFeaturesList));
    meanAccuraciesSVM_plv = zeros(length(classProbs), length(numFeaturesList));
    meanAccuraciesDT_plv = zeros(length(classProbs), length(numFeaturesList));

    % Example: Use classLabels_corr for your classification
    for classPairIdx = 1:length(classProbs)

        % Define the data for each category
        categories = {'correlation', 'coherence', 'phaseLocking'};
        for categoryIdx = 1:length(categories)
            category = categories{categoryIdx};

            switch category
                case 'correlation'
                    classLabels = {BandValues_REST_corr, BandValues_MI_FT_corr, BandValues_MI_DB_corr, ...
                        BandValues_ME_DB_corr, BandValues_ME_FT_corr};
                case 'coherence'
                    classLabels = {BandValues_REST_coh, BandValues_MI_FT_coh, BandValues_MI_DB_coh, ...
                        BandValues_ME_DB_coh, BandValues_ME_FT_coh};
                case 'phaseLocking'
                    classLabels = {BandValues_REST_plv, BandValues_MI_FT_plv, BandValues_MI_DB_plv, ...
                        BandValues_ME_DB_plv, BandValues_ME_FT_plv};
            end

            % Main loop for cross-validation for the current category
            for iteration = 1:numIteration
                for c = 1:numel(classProbs)
                    tasks = classProbs{c};
                    taskA = classLabels{tasks(1)};
                    taskB = classLabels{tasks(2)};

                    % Labels: 1 for task A, -1 for task B
                    labels = [ones(size(taskA, 1), 1); -1 * ones(size(taskB, 1), 1)];

                    % Combine data from both tasks
                    combinedData = cat(1, taskA, taskB);

                    % Stratified cross-validation
                    cv = cvpartition(labels, 'KFold', K, 'Stratify', true);

                    for fold = 1:K
                        testIdx = cv.test(fold);
                        trainIdx = cv.training(fold);

                        % Train and test sets
                        trainData = combinedData(trainIdx, :);
                        testData = combinedData(testIdx, :);

                        trainLabels = labels(trainIdx);
                        testLabels = labels(testIdx);

                        % Select features using MRMR
                        selectedFeaturesIdxTrain = fscmrmr(trainData, trainLabels);

                        % Store selected features and perform classification
                        for fIdx = 1:length(numFeaturesList)
                            numFeatures = numFeaturesList(fIdx);
                            if size(trainData, 2) < numFeatures
                                continue;
                            end

                            selectedFeaturesTrain = trainData(:, selectedFeaturesIdxTrain(1:numFeatures));
                            selectedFeaturesTest = testData(:, selectedFeaturesIdxTrain(1:numFeatures));

                            % Perform LDA classification
                            ClassifierLDA = fitcdiscr(selectedFeaturesTrain, trainLabels);
                            accuracyLDA = mean(predict(ClassifierLDA, selectedFeaturesTest) == testLabels) * 100;

                            % Perform SVM classification
                            ClassifierSVM = fitcsvm(selectedFeaturesTrain, trainLabels);
                            accuracySVM = mean(predict(ClassifierSVM, selectedFeaturesTest) == testLabels) * 100;

                            % Perform Decision DT classification
                            ClassifierDT = fitctree(selectedFeaturesTrain, trainLabels);
                            accuracyDT = mean(predict(ClassifierDT, selectedFeaturesTest) == testLabels) * 100;

                            % Store accuracies for each category
                            switch category
                                case 'correlation'
                                    localAccuracies_lda_corr(c, fold, iteration, fIdx) = accuracyLDA;
                                    localAccuracies_svm_corr(c, fold, iteration, fIdx) = accuracySVM;
                                    localAccuracies_dt_corr(c, fold, iteration, fIdx) = accuracyDT;
                                case 'coherence'
                                    localAccuracies_lda_coh(c, fold, iteration, fIdx) = accuracyLDA;
                                    localAccuracies_svm_coh(c, fold, iteration, fIdx) = accuracySVM;
                                    localAccuracies_dt_coh(c, fold, iteration, fIdx) = accuracyDT;
                                case 'phaseLocking'
                                    localAccuracies_lda_plv(c, fold, iteration, fIdx) = accuracyLDA;
                                    localAccuracies_svm_plv(c, fold, iteration, fIdx) = accuracySVM;
                                    localAccuracies_dt_plv(c, fold, iteration, fIdx) = accuracyDT;
                            end
                        end
                    end
                end

                % Compute the mean accuracy for each condition
                for c = 1:numel(classProbs)
                    for fIdx = 1:length(numFeaturesList)
                        switch category
                            case 'correlation'
                                meanAccuraciesLDA_corr(c, fIdx) = mean(localAccuracies_lda_corr(c, :, :, fIdx), 'all');
                                meanAccuraciesSVM_corr(c, fIdx) = mean(localAccuracies_svm_corr(c, :, :, fIdx), 'all');
                                meanAccuraciesDT_corr(c, fIdx) = mean(localAccuracies_dt_corr(c, :, :, fIdx), 'all');
                            case 'coherence'
                                meanAccuraciesLDA_coh(c, fIdx) = mean(localAccuracies_lda_coh(c, :, :, fIdx), 'all');
                                meanAccuraciesSVM_coh(c, fIdx) = mean(localAccuracies_svm_coh(c, :, :, fIdx), 'all');
                                meanAccuraciesDT_coh(c, fIdx) = mean(localAccuracies_dt_coh(c, :, :, fIdx), 'all');
                            case 'phaseLocking'
                                meanAccuraciesLDA_plv(c, fIdx) = mean(localAccuracies_lda_plv(c, :, :, fIdx), 'all');
                                meanAccuraciesSVM_plv(c, fIdx) = mean(localAccuracies_svm_plv(c, :, :, fIdx), 'all');
                                meanAccuraciesDT_plv(c, fIdx) = mean(localAccuracies_dt_plv(c, :, :, fIdx), 'all');
                        end
                    end
                end
            end
        end
    end
    % Save mean accuracies for each category
    save(fullfile(myPath, sprintf('meanAccuracies_corr_9b_22ch_Sub-%03d.mat', numSubjects)), 'meanAccuraciesLDA_corr', 'meanAccuraciesSVM_corr', 'meanAccuraciesDT_corr');
    save(fullfile(myPath, sprintf('meanAccuracies_coh_9b_22ch_Sub-%03d.mat', numSubjects)), 'meanAccuraciesLDA_coh', 'meanAccuraciesSVM_coh', 'meanAccuraciesDT_coh');
    save(fullfile(myPath, sprintf('meanAccuracies_pLV_9b_22ch_Sub-%03d.mat', numSubjects)), 'meanAccuraciesLDA_plv', 'meanAccuraciesSVM_plv', 'meanAccuraciesDT_plv');

    %% Save Local Accuracies (LDA, SVM, DT)
    save(fullfile(myPath, sprintf('localAccuracies_lda_corr_9b_22ch_Sub-%03d.mat', numSubjects)), 'localAccuracies_lda_corr');
    save(fullfile(myPath, sprintf('localAccuracies_lda_coh_9b_22ch_Sub-%03d.mat', numSubjects)), 'localAccuracies_lda_coh');
    save(fullfile(myPath, sprintf('localAccuracies_lda_plv_9b_22ch_Sub-%03d.mat', numSubjects)), 'localAccuracies_lda_plv');

    save(fullfile(myPath, sprintf('localAccuracies_svm_corr_9b_22ch_Sub-%03d.mat', numSubjects)), 'localAccuracies_svm_corr');
    save(fullfile(myPath, sprintf('localAccuracies_svm_coh_9b_22ch_Sub-%03d.mat', numSubjects)), 'localAccuracies_svm_coh');
    save(fullfile(myPath, sprintf('localAccuracies_svm_plv_9b_22ch_Sub-%03d.mat', numSubjects)), 'localAccuracies_svm_plv');

    save(fullfile(myPath, sprintf('localAccuracies_dt_corr_9b_22ch_Sub-%03d.mat', numSubjects)), 'localAccuracies_dt_corr');
    save(fullfile(myPath, sprintf('localAccuracies_dt_coh_9b_22ch_Sub-%03d.mat', numSubjects)), 'localAccuracies_dt_coh');
    save(fullfile(myPath, sprintf('localAccuracies_dt_plv_9b_22ch_Sub-%03d.mat', numSubjects)), 'localAccuracies_dt_plv');

    % Display completion message, elapsed time, and memory usage
    elapsedTime = toc; % Stop the timer
    fprintf('Subject %d completed in %.2f seconds.\n', numSubjects, elapsedTime);
end

% Initialize accumulators for overall mean accuracies
overallMeanAccuraciesLDA_corr = zeros(length(classProbs), length(numFeaturesList));
overallMeanAccuraciesSVM_corr = zeros(length(classProbs), length(numFeaturesList));
overallMeanAccuraciesDT_corr = zeros(length(classProbs), length(numFeaturesList));

overallMeanAccuraciesLDA_coh = zeros(length(classProbs), length(numFeaturesList));
overallMeanAccuraciesSVM_coh = zeros(length(classProbs), length(numFeaturesList));
overallMeanAccuraciesDT_coh = zeros(length(classProbs), length(numFeaturesList));

overallMeanAccuraciesLDA_plv = zeros(length(classProbs), length(numFeaturesList));
overallMeanAccuraciesSVM_plv = zeros(length(classProbs), length(numFeaturesList));
overallMeanAccuraciesDT_plv = zeros(length(classProbs), length(numFeaturesList));

% After subject loop, calculate overall mean accuracies by averaging across valid subjects
for subjectIdx = 1:totalSubjects
    % Check if the subject was excluded
    if ismember(subjectIdx, excludeSubjects)
        continue;
    end

    % Accumulate mean accuracies for each subject
    overallMeanAccuraciesLDA_corr = overallMeanAccuraciesLDA_corr + meanAccuraciesLDA_corr;
    overallMeanAccuraciesSVM_corr = overallMeanAccuraciesSVM_corr + meanAccuraciesSVM_corr;
    overallMeanAccuraciesDT_corr = overallMeanAccuraciesDT_corr + meanAccuraciesDT_corr;

    overallMeanAccuraciesLDA_coh = overallMeanAccuraciesLDA_coh + meanAccuraciesLDA_coh;
    overallMeanAccuraciesSVM_coh = overallMeanAccuraciesSVM_coh + meanAccuraciesSVM_coh;
    overallMeanAccuraciesDT_coh = overallMeanAccuraciesDT_coh + meanAccuraciesDT_coh;

    overallMeanAccuraciesLDA_plv = overallMeanAccuraciesLDA_plv + meanAccuraciesLDA_plv;
    overallMeanAccuraciesSVM_plv = overallMeanAccuraciesSVM_plv + meanAccuraciesSVM_plv;
    overallMeanAccuraciesDT_plv = overallMeanAccuraciesDT_plv + meanAccuraciesDT_plv;
end

% Calculate the overall mean accuracies by dividing by the number of valid subjects
validSubjectCount = totalSubjects - length(excludeSubjects);

% Divide by the number of valid subjects to get the average
overallMeanAccuraciesLDA_corr = overallMeanAccuraciesLDA_corr / validSubjectsCount;
overallMeanAccuraciesSVM_corr = overallMeanAccuraciesSVM_corr / validSubjectsCount;
overallMeanAccuraciesDT_corr = overallMeanAccuraciesDT_corr / validSubjectsCount;

overallMeanAccuraciesLDA_coh = overallMeanAccuraciesLDA_coh / validSubjectsCount;
overallMeanAccuraciesSVM_coh = overallMeanAccuraciesSVM_coh / validSubjectsCount;
overallMeanAccuraciesDT_coh = overallMeanAccuraciesDT_coh / validSubjectsCount;

overallMeanAccuraciesLDA_plv = overallMeanAccuraciesLDA_plv / validSubjectsCount;
overallMeanAccuraciesSVM_plv = overallMeanAccuraciesSVM_plv / validSubjectsCount;
overallMeanAccuraciesDT_plv = overallMeanAccuraciesDT_plv / validSubjectsCount;


% Create a structure to store the overall mean accuracies for each measure and classifier type
overallMeanAccuracies.correlation.LDA = overallMeanAccuraciesLDA_corr;
overallMeanAccuracies.correlation.SVM = overallMeanAccuraciesSVM_corr;
overallMeanAccuracies.correlation.DT = overallMeanAccuraciesDT_corr;

overallMeanAccuracies.coherence.LDA = overallMeanAccuraciesLDA_coh;
overallMeanAccuracies.coherence.SVM = overallMeanAccuraciesSVM_coh;
overallMeanAccuracies.coherence.DT = overallMeanAccuraciesDT_coh;

overallMeanAccuracies.phaseLocking.LDA = overallMeanAccuraciesLDA_plv;
overallMeanAccuracies.phaseLocking.SVM = overallMeanAccuraciesSVM_plv;
overallMeanAccuracies.phaseLocking.DT = overallMeanAccuraciesDT_plv;

% Save the structure to a .mat file
save(fullfile(myPath, sprintf('overall_FC_9b_22ch.mat')), 'overallMeanAccuracies');

fprintf('Overall mean accuracies have been saved to overallMeanAccuracies.mat.\n');


% Display overall mean accuracies
fprintf('Overall Mean Accuracies (LDA - Correlation):\n');
disp(overallMeanAccuraciesLDA_corr);
fprintf('Overall Mean Accuracies (SVM - Correlation):\n');
disp(overallMeanAccuraciesSVM_corr);
fprintf('Overall Mean Accuracies (DT - Correlation):\n');
disp(overallMeanAccuraciesDT_corr);

fprintf('Overall Mean Accuracies (LDA - Coherence):\n');
disp(overallMeanAccuraciesLDA_coh);
fprintf('Overall Mean Accuracies (SVM - Coherence):\n');
disp(overallMeanAccuraciesSVM_coh);
fprintf('Overall Mean Accuracies (DT - Coherence):\n');
disp(overallMeanAccuraciesDT_coh);

fprintf('Overall Mean Accuracies (LDA - plv Locking):\n');
disp(overallMeanAccuraciesLDA_plv);
fprintf('Overall Mean Accuracies (SVM - plv Locking):\n');
disp(overallMeanAccuraciesSVM_plv);
fprintf('Overall Mean Accuracies (DT - plv Locking):\n');
disp(overallMeanAccuraciesDT_plv);
