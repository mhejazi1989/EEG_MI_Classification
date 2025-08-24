clear all;
close all;
clc;

% Define the total number of subjects
totalSubjects = 29;

% List of subjects to exclude
excludeSubjects = [6, 8, 10, 13, 21, 23, 25, 27];

% Turn off the specific warning
warning('off', 'MATLAB:colon:operandsNotInteger');
numChannels = 22;
fbands = 5;        % Adjust according to the number of frequency bands you're using

for numSubjects = 1:totalSubjects
    % Skip excluded subjects
    if any(numSubjects == excludeSubjects)
        continue;
    end

    % Start measuring time
    tic;

    % Define a cell array of participant numbers
    subjectIDs = cell(1, totalSubjects);
    subjectIDs{numSubjects} = sprintf('%03d', numSubjects);

    subjectFolder = sprintf('Sub-%03d', numSubjects);
    mypath = fullfile('/lustre03/project/6067835/mhejazi/EEG-DATA', subjectFolder);

    % Load ALLEEG from the saved .mat file
    load(fullfile(mypath, sprintf('concatenated_5b_values_sub%03d.mat', numSubjects)), ...
        'BandValues_REST_new', 'BandValues_MI_FT_new', 'BandValues_MI_DB_new', 'BandValues_ME_DB_new', 'BandValues_ME_FT_new');

    BandValues_REST = BandValues_REST_new;
    BandValues_MI_FT = BandValues_MI_FT_new;
    BandValues_MI_DB = BandValues_MI_DB_new;
    BandValues_ME_DB = BandValues_ME_DB_new;
    BandValues_ME_FT = BandValues_ME_FT_new;

    % Classification based on Band values
    classes = [1 2 3 4 5 6]; % Define your classes
    numIteration = 15;
    K = 5;
    numSamples = 1025;
    classProbs = {
        [1, 2], % Condition 1: Task 1 vs Task 2
        [1, 3], % Condition 2: Task 1 vs Task 3
        [1, 4], % Condition 3: Task 1 vs Task 4
        [1, 5], % Condition 4: Task 1 vs Task 5
        [2, 3], % Condition 5: Task 2 vs Task 3
        [4, 5], % Condition 5: Task 2 vs Task 3
        };
    % Create the Labels cell array with matrices
    classLabels = {
        BandValues_REST;
        BandValues_MI_FT;
        BandValues_MI_DB;
        BandValues_ME_DB;
        BandValues_ME_FT;
        };
    % Define condition names
    conditionNames = {
        'Condition 1: rest vs MI fingertips';
        'Condition 2: rest vs MI duckbill';
        'Condition 3: rest vs ME fingeratips';
        'Condition 4: rest vs ME duckbill';
        'Condition 5: MI fingertips vs MI duckbill';
        'Condition 6: ME fingertips vs ME duckbill';
        };

    % Define the number of features to evaluate
    numFeaturesList = [10, 20, 30, 40, 50];

    % Initialize cells for storing accuracies separately for LDA, SVM, and Decision Tree
    localAccuraciesLDA = cell(length(classProbs), K, numIteration, length(numFeaturesList));
    localAccuraciesSVM = cell(length(classProbs), K, numIteration, length(numFeaturesList));
    localAccuraciesDT = cell(length(classProbs), K, numIteration, length(numFeaturesList));

    subjectMeanAccuraciesSVM = cell(numSubjects, 1);
    subjectMeanAccuraciesDT = cell(numSubjects, 1);
    subjectMeanAccuraciesLDA = cell(numSubjects, 1);
    % CSP bfore the classification loop
    numCSPComponents = 5; % Choose the number of CSP components

    % Initialize cell arrays to store CSP filters for each band
    W_bands = cell(fbands, 1);

    % Initialize cell arrays to store transformed data for each frequency band
    transformedTrainData = cell(fbands, 1);
    transformedTestData = cell(fbands, 1);
    % Define the classifiers
    classifiers = {'LDA', 'SVM', 'DecisionTree'};

    for iteration = 1:numIteration
        for c = 1:numel(classProbs)
            tasks = classProbs{c};
            taskA = classLabels{tasks(1)};
            taskB = classLabels{tasks(2)};
            labels = [ones(size(taskA{1,1}, 1), 1); -1 * ones(size(taskB{1,1}, 1), 1)];

            % Stratified cross-validation
            cv = cvpartition(labels, 'KFold', K, 'Stratify', true);

            for fold = 1:K
                allTrainFeatures = [];
                allTestFeatures = [];

                for band = 1:fbands
                    taskA_bands = taskA{band};
                    taskB_bands = taskB{band};

                    % Prepare the data for k-fold cross-validation
                    combinedData = cat(1, taskA_bands, taskB_bands);

                    testIdx = cv.test(fold);
                    trainIdx = cv.training(fold);

                    trainData = combinedData(trainIdx, :, :);
                    testData = combinedData(testIdx, :, :);

                    trainLabels = labels(trainIdx);
                    testLabels = labels(testIdx);

                    if isempty(trainData) || isempty(trainLabels)
                        continue; % Skip empty training data
                    end

                    % Separate data for taskA and taskB
                    taskA_train = trainData(trainLabels == 1, :, :);
                    taskB_train = trainData(trainLabels == -1, :, :);

                    % Initialize the covariance matrices
                    S1 = zeros(numChannels, numChannels);
                    S2 = zeros(numChannels, numChannels);

                    % Accumulate covariance matrices for class A (taskA_train)
                    for i = 1:size(taskA_train, 1)
                        S1 = S1 + cov(squeeze(taskA_train(i, :, :))');
                    end
                    % Accumulate covariance matrices for class B (taskB_train)
                    for i = 1:size(taskB_train, 1)
                        S2 = S2 + cov(squeeze(taskB_train(i, :, :))');
                    end

                    % Normalize by the number of trials
                    S1 = S1 / size(taskA_train, 1);
                    S2 = S2 / size(taskB_train, 1);

                    % Compute CSP filters using the combined covariance matrices
                    % Solve the eigenvalue problem S1Â·W = lÂ·S2Â·W
                    [W, L] = eig(S1, S1 + S2);   % Mixing matrix W (spatial filters are columns)
                    lambda = diag(L);            % Eigenvalues
                    A = inv(W)';                 % Demixing matrix

                    % Select the first and last five components
                    W_bands{band} = W;

                    % Transform the training data for the current frequency band
                    transformedTrainData{band} = arrayfun(@(i) computeCSPFeatures(squeeze(trainData(i, :, :)), W_bands{band}, numCSPComponents), 1:size(trainData, 1), 'UniformOutput', false);
                    transformedTrainData{band} = cell2mat(transformedTrainData{band}');

                    % Transform the testing data for the current frequency band
                    transformedTestData{band} = arrayfun(@(i) computeCSPFeatures(squeeze(testData(i, :, :)), W_bands{band}, numCSPComponents), 1:size(testData, 1), 'UniformOutput', false);
                    transformedTestData{band} = cell2mat(transformedTestData{band}');

                    % Concatenate transformed train data from all frequency bands
                    allTrainFeatures = [allTrainFeatures, transformedTrainData{band}];

                    % Concatenate transformed test data from all frequency bands
                    allTestFeatures = [allTestFeatures, transformedTestData{band}];
                end

                selectedFeaturesIdxTrain = fscmrmr(allTrainFeatures, trainLabels);

                % Evaluate accuracy for each numFeatures in numFeaturesList
                for fIdx = 1:length(numFeaturesList)
                    numFeatures = numFeaturesList(fIdx);
                    if size(allTrainFeatures, 2) < numFeatures
                        continue; % Skip if there are fewer features than numFeatures
                    end

                    % Select top numFeatures from training data
                    selectedFeaturesIdxTrain1 = selectedFeaturesIdxTrain(1:numFeatures);

                    % Select top numFeatures from test data using indices from training selection
                    selectedFeaturesTrain = allTrainFeatures(:, selectedFeaturesIdxTrain1);
                    selectedFeaturesTest = allTestFeatures(:, selectedFeaturesIdxTrain1);

                    % Loop through each classifier
                    for clfIdx = 1:numel(classifiers)
                        switch classifiers{clfIdx}
                            case 'LDA'
                                % Train LDA classifier
                                LDAModel = fitcdiscr(selectedFeaturesTrain, trainLabels, 'DiscrimType', 'linear');
                                % Test LDA classifier
                                predictedLabels = predict(LDAModel, selectedFeaturesTest);
                                % Calculate accuracy
                                ldaAccuracy = mean(predictedLabels == testLabels);
                                localAccuraciesLDA{c, fold, iteration, fIdx} = ldaAccuracy;

                            case 'SVM'
                                % Train SVM classifier
                                SVMModel = fitcsvm(selectedFeaturesTrain, trainLabels, 'KernelFunction', 'linear', 'Standardize', true);
                                % Predict using SVM
                                predictedLabels = predict(SVMModel, selectedFeaturesTest);

                                % Calculate accuracy for SVM
                                accuracySVM = sum(predictedLabels == testLabels) / length(testLabels);
                                localAccuraciesSVM{c, fold, iteration, fIdx} = accuracySVM;

                            case 'DecisionTree'
                                % Train Decision Tree classifier
                                treeModel = fitctree(selectedFeaturesTrain, trainLabels);
                                % Predict using Decision Tree
                                predictedLabels = predict(treeModel, selectedFeaturesTest);

                                % Calculate accuracy for Decision Tree
                                accuracyDT = sum(predictedLabels == testLabels) / length(testLabels);
                                localAccuraciesDT{c, fold, iteration, fIdx} = accuracyDT;
                        end
                    end
                end
            end
        end
    end
    % Initialize arrays to store mean accuracies for SVM, Decision Tree, and LDA
    meanAccuracies_svm = cell(length(classProbs), length(numFeaturesList));
    meanAccuracies_dt = cell(length(classProbs), length(numFeaturesList));
    meanAccuracies_lda = cell(length(classProbs), length(numFeaturesList));

    % Compute the mean accuracy for each condition, feature index, and classifier
    for c = 1:numel(classProbs)
        for fIdx = 1:length(numFeaturesList)
            % Extract accuracies for the current condition and feature index
            accuracies_svm = localAccuraciesSVM(c, :, :, fIdx); % SVM accuracies
            accuracies_dt = localAccuraciesDT(c, :, :, fIdx);   % Decision Tree accuracies
            accuracies_lda = localAccuraciesLDA(c, :, :, fIdx); % LDA accuracies

            % Compute the mean accuracy over folds and iterations
            meanAccuracies_svm{c, fIdx} = mean(cell2mat(accuracies_svm), 'all') * 100;
            meanAccuracies_dt{c, fIdx} = mean(cell2mat(accuracies_dt), 'all') * 100;
            meanAccuracies_lda{c, fIdx} = mean(cell2mat(accuracies_lda), 'all') * 100;
        end
    end

    % Display mean accuracies for each classifier
    for c = 1:numel(classProbs)
        fprintf('Condition %d:\n', c);
        for fIdx = 1:length(numFeaturesList)
            fprintf('  Number of features: %d\n', numFeaturesList(fIdx));
            fprintf('    SVM Mean Accuracy: %.4f\n', meanAccuracies_svm{c, fIdx});
            fprintf('    Decision Tree Mean Accuracy: %.4f\n', meanAccuracies_dt{c, fIdx});
            fprintf('    LDA Mean Accuracy: %.4f\n', meanAccuracies_lda{c, fIdx});
        end
    end

    % Save meanAccuracies for this subject
    meanAccuraciesFileSVM = fullfile(mypath, sprintf('CSP_Accuracies_svm_5b_22ch_Sub-%03d.mat', numSubjects));
    meanAccuraciesFileDT = fullfile(mypath, sprintf('CSP_Accuracies_dt_5b_22ch_Sub-%03d.mat', numSubjects));
    meanAccuraciesFileLDA = fullfile(mypath, sprintf('CSP_Accuracies_lda_5b_22ch_Sub-%03d.mat', numSubjects));
     save(meanAccuraciesFileSVM, 'meanAccuracies_svm', 'localAccuraciesSVM');
    save(meanAccuraciesFileDT, 'meanAccuracies_dt', 'localAccuraciesDT');
    save(meanAccuraciesFileLDA, 'meanAccuracies_lda', 'localAccuraciesLDA');

    % Store the results for this subject
    subjectMeanAccuraciesSVM{numSubjects} = meanAccuracies_svm;
    subjectMeanAccuraciesDT{numSubjects} = meanAccuracies_dt;
    subjectMeanAccuraciesLDA{numSubjects} = meanAccuracies_lda;

    toc;
    % Display elapsed time for this subject
    fprintf('Subject %d completed in %.2f seconds.\n', numSubjects, toc);

end
% Define the final save path
finalSavePath = fullfile(mypath, 'AllSubjects_AccuracyResults_22ch_5b.mat');

% Save the results in a .mat file
save(finalSavePath, 'subjectMeanAccuraciesSVM', 'subjectMeanAccuraciesDT', 'subjectMeanAccuraciesLDA');

% Confirmation message
fprintf('All subjects\\'' accuracy results saved successfully in %s\n', finalSavePath);


% Restore the figure visibility setting
set(0, 'DefaultFigureVisible', 'on');

% Suppress warning
warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function transformedData = computeCSPFeatures(trialData, W, numCSPComponents)
% Initialize array to store CSP features for the trial data
numChannels = size(trialData, 1);
numSamples = size(trialData, 2);
transformedData = zeros(1, numCSPComponents);

% Iterate over CSP components
for com = 1:numCSPComponents
    % Extract the corresponding spatial filter from W
    Wb = [W(:, com), W(:, end - com + 1)];

    % Compute the covariance matrix Zb_i using the spatial filters
    Zb_i = Wb' * trialData * trialData' * Wb;

    % Compute the diagonal elements of Zb_i
    diag_Zb = diag(Zb_i);

    % Compute the trace of Zb_i
    trace_Zb = trace(Zb_i);

    % Compute the log-variance features (vb,i)
    vb_i = log(diag_Zb / trace_Zb);

    % Store vb_i in transformedData
    transformedData((com - 1) * 2 + 1:com * 2) = vb_i';
end
end




