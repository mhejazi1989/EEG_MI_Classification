clear all;
close all;
clc;
warning off;
% Conn for 64 channels

% Define the total number of subjects
totalSubjects = 29;

% List of subjects to exclude
excludeSubjects = [6, 8, 10, 13, 21, 23, 25, 27];

for numSubjects = 5:8
    % Skip excluded subjects
    if any(numSubjects == excludeSubjects)
        continue;
    end

    % Start measuring time
    tic;
    % Define a cell array of participant numbers
    subjectIDs = sprintf('%03d', numSubjects);
    subjectFolder = sprintf('Sub-%03d', numSubjects);
    myPath = fullfile('/lustre03/project/6067835/mhejazi/EEG-DATA', subjectFolder);

    % Check if the directory exists
    if ~exist(myPath, 'dir')
        fprintf('Error: The directory %s does not exist.\n', myPath);
        continue; % Skip this subject if the directory does not exist
    end

    % Define the path to the .mat file for 'coh' features
    matFilePath = fullfile(myPath, sprintf('sub%03d_5b_AllFeatureMatrix.mat', numSubjects));
    if exist(matFilePath, 'file') ~= 2
        fprintf('Error loading coh file for Subject %d: Unable to find file "%s".\n', numSubjects, matFilePath);
        continue;  % Skip this subject if the file does not exist
    end

    % Load the structure from the .mat file
    loadedData = load(matFilePath);  % This loads the entire structure

    % Store the 'corr' features with unique names
    corr_AllFeatureMatrix_ME_DB = loadedData.AllFeatureMatrix_ME_DB_corr;
    corr_AllFeatureMatrix_MI_DB = loadedData.AllFeatureMatrix_MI_DB_corr;
    corr_AllFeatureMatrix_ME_FT = loadedData.AllFeatureMatrix_ME_FT_corr;
    corr_AllFeatureMatrix_MI_FT = loadedData.AllFeatureMatrix_MI_FT_corr;
    corr_AllFeatureMatrix_REST = loadedData.AllFeatureMatrix_REST_corr;


    % Store the 'coh' features with unique names
    coh_AllFeatureMatrix_ME_DB = loadedData.AllFeatureMatrix_ME_DB_coh;
    coh_AllFeatureMatrix_MI_DB = loadedData.AllFeatureMatrix_ME_FT_coh;
    coh_AllFeatureMatrix_ME_FT = loadedData.AllFeatureMatrix_MI_DB_coh;
    coh_AllFeatureMatrix_MI_FT = loadedData.AllFeatureMatrix_MI_FT_coh;
    coh_AllFeatureMatrix_REST = loadedData.AllFeatureMatrix_REST_coh;

    % Store the 'plv' features with unique names
    plv_AllFeatureMatrix_ME_DB = loadedData.AllFeatureMatrix_ME_DB_plv;
    plv_AllFeatureMatrix_MI_DB = loadedData.AllFeatureMatrix_ME_FT_plv;
    plv_AllFeatureMatrix_ME_FT = loadedData.AllFeatureMatrix_MI_DB_plv;
    plv_AllFeatureMatrix_MI_FT = loadedData.AllFeatureMatrix_MI_FT_plv;
    plv_AllFeatureMatrix_REST = loadedData.AllFeatureMatrix_REST_plv;

    % Classification based on Band values
    classes = 1:6; % Define your classes
    numIteration = 15;
    K = 5;
    classProbs = {[1, 2], [1, 3], [1, 4], [1, 5], [2, 3], [4, 5]}; % Adjusted to remove duplicates
    fbands = 5;
    MaxFeatureNum = 10*fbands;

    % Define the classifiers
    classifiers = {'LDA', 'SVM', 'DecisionTree'};

    % Create the Labels cell array with matrices
    AllFeatureMatrix_corr = {corr_AllFeatureMatrix_REST, corr_AllFeatureMatrix_MI_FT, corr_AllFeatureMatrix_MI_DB, ...
        corr_AllFeatureMatrix_ME_DB, corr_AllFeatureMatrix_ME_FT};
    AllFeatureMatrix_coh = {coh_AllFeatureMatrix_REST, coh_AllFeatureMatrix_MI_FT, coh_AllFeatureMatrix_MI_DB, ...
        coh_AllFeatureMatrix_ME_DB, coh_AllFeatureMatrix_ME_FT};
    AllFeatureMatrix_plv = {plv_AllFeatureMatrix_REST, plv_AllFeatureMatrix_MI_FT, plv_AllFeatureMatrix_MI_DB, ...
        plv_AllFeatureMatrix_ME_DB, plv_AllFeatureMatrix_ME_FT};

    % Define condition names
    conditionNames = {'Condition 1: rest vs MI fingertips', 'Condition 2: rest vs MI duckbill', ...
        'Condition 3: rest vs ME fingertips', 'Condition 4: rest vs ME duckbill', ...
        'Condition 5: MI fingertips vs MI duckbill', 'Condition 6: ME fingertips vs ME duckbill'};

    % Define the number of features to evaluate
    numFeaturesList = [10, 20, 30, 40, MaxFeatureNum];

    % Initialize storage for mean accuracies for each classifier
    meanAccuracies_corr_SVM = zeros(length(classProbs), length(numFeaturesList));
    meanAccuracies_coh_SVM = zeros(length(classProbs), length(numFeaturesList));
    meanAccuracies_plv_SVM = zeros(length(classProbs), length(numFeaturesList));

    meanAccuracies_corr_DT = zeros(length(classProbs), length(numFeaturesList));
    meanAccuracies_coh_DT = zeros(length(classProbs), length(numFeaturesList));
    meanAccuracies_plv_DT = zeros(length(classProbs), length(numFeaturesList));

    meanAccuracies_corr_LDA = zeros(length(classProbs), length(numFeaturesList));
    meanAccuracies_coh_LDA = zeros(length(classProbs), length(numFeaturesList));
    meanAccuracies_plv_LDA = zeros(length(classProbs), length(numFeaturesList));

    % Initialize 4D arrays to hold local accuracies for each classifier
    localAccuracies_svm_corr = zeros(length(classProbs), K, numIteration, length(numFeaturesList));
    localAccuracies_svm_coh = zeros(length(classProbs), K, numIteration, length(numFeaturesList));
    localAccuracies_svm_plv = zeros(length(classProbs), K, numIteration, length(numFeaturesList));

    localAccuracies_dt_corr = zeros(length(classProbs), K, numIteration, length(numFeaturesList));
    localAccuracies_dt_coh = zeros(length(classProbs), K, numIteration, length(numFeaturesList));
    localAccuracies_dt_plv = zeros(length(classProbs), K, numIteration, length(numFeaturesList));

    localAccuracies_lda_corr = zeros(length(classProbs), K, numIteration, length(numFeaturesList));
    localAccuracies_lda_coh = zeros(length(classProbs), K, numIteration, length(numFeaturesList));
    localAccuracies_lda_plv = zeros(length(classProbs), K, numIteration, length(numFeaturesList));

    % Initialize a new cell array to hold the selected top 50 feature indices for each iteration and fold
    selectedTopFeatureIndices_corr = cell(numIteration, K);
    selectedTopFeatureIndices_coh = cell(numIteration, K);
    selectedTopFeatureIndices_plv = cell(numIteration, K);

    % Preallocate cell arrays for selected features for the current subject
    selectedFeaturesTrainSubject_corr = cell(numIteration, 1);
    selectedFeaturesTestSubject_corr = cell(numIteration, 1);
    selectedFeaturesTrainSubject_coh = cell(numIteration, 1);
    selectedFeaturesTestSubject_coh = cell(numIteration, 1);
    selectedFeaturesTrainSubject_plv = cell(numIteration, 1);
    selectedFeaturesTestSubject_plv = cell(numIteration, 1);

    selectedFeatureIndices_corr = cell(numIteration, K);
    selectedFeatureIndices_coh = cell(numIteration, K);
    selectedFeatureIndices_plv = cell(numIteration, K);

    for iteration = 1:numIteration
        for c = 1:numel(classProbs)
            tasks = classProbs{c};

            % Load the feature matrices for both coh
            taskA_corr = AllFeatureMatrix_corr{tasks(1)};
            taskB_corr = AllFeatureMatrix_corr{tasks(2)};
            labels_corr = [ones(size(taskA_corr, 1), 1); -1 * ones(size(taskB_corr, 1), 1)];

            taskA_coh = AllFeatureMatrix_coh{tasks(1)};
            taskB_coh = AllFeatureMatrix_coh{tasks(2)};
            labels_coh = [ones(size(taskA_coh, 1), 1); -1 * ones(size(taskB_coh, 1), 1)];

            taskA_plv = AllFeatureMatrix_plv{tasks(1)};
            taskB_plv = AllFeatureMatrix_plv{tasks(2)};
            labels_plv = [ones(size(taskA_plv, 1), 1); -1 * ones(size(taskB_plv, 1), 1)];


            % Stratified cross-validation
            cv_corr = cvpartition(labels_corr, 'KFold', K, 'Stratify', true);
            combinedData_corr = cat(1, taskA_corr, taskB_corr); % Prepare combined data for coherence

            cv_coh = cvpartition(labels_coh, 'KFold', K, 'Stratify', true);
            combinedData_coh = cat(1, taskA_coh, taskB_coh); % Prepare combined data for coherence

            cv_plv = cvpartition(labels_plv, 'KFold', K, 'Stratify', true);
            combinedData_plv = cat(1, taskA_plv, taskB_plv); % Prepare combined data for coherence


            for fold = 1:K
                testIdx_corr = cv_corr.test(fold);
                trainIdx_corr = cv_corr.training(fold);

                testIdx_coh = cv_coh.test(fold);
                trainIdx_coh = cv_coh.training(fold);

                testIdx_plv = cv_plv.test(fold);
                trainIdx_plv = cv_plv.training(fold);

                % Split the data into training and test sets for both coh
                trainData_corr = combinedData_corr(trainIdx_corr, :, :);
                testData_corr = combinedData_corr(testIdx_corr, :, :);

                trainData_coh = combinedData_coh(trainIdx_coh, :, :);
                testData_coh = combinedData_coh(testIdx_coh, :, :);

                trainData_plv = combinedData_plv(trainIdx_plv, :, :);
                testData_plv = combinedData_plv(testIdx_plv, :, :);

                trainLabels_corr = labels_corr(trainIdx_corr);
                testLabels_corr = labels_corr(testIdx_corr);

                trainLabels_coh = labels_coh(trainIdx_coh);
                testLabels_coh = labels_coh(testIdx_coh);

                trainLabels_plv = labels_plv(trainIdx_plv);
                testLabels_plv = labels_plv(testIdx_plv);

                % Select features using MRMR for both coh
                selectedFeaturesIdxTrain_corr = fscmrmr(trainData_corr, trainLabels_corr); % For coherence
                selectedFeaturesIdxTrain_coh = fscmrmr(trainData_coh, trainLabels_coh);
                selectedFeaturesIdxTrain_plv = fscmrmr(trainData_plv, trainLabels_plv);


                % Store the top 50 feature indices for both coh
                selectedTopFeatureIndices_corr{iteration, fold} = selectedFeaturesIdxTrain_corr(1:MaxFeatureNum);
                selectedTopFeatureIndices_coh{iteration, fold} = selectedFeaturesIdxTrain_coh(1:MaxFeatureNum);
                selectedTopFeatureIndices_plv{iteration, fold} = selectedFeaturesIdxTrain_plv(1:MaxFeatureNum);

                % Store the selected feature indices
                selectedFeatureIndices_corr{iteration, fold} = selectedFeaturesIdxTrain_corr;
                selectedFeatureIndices_coh{iteration, fold} = selectedFeaturesIdxTrain_coh;
                selectedFeatureIndices_plv{iteration, fold} = selectedFeaturesIdxTrain_plv;

                for fIdx = 1:length(numFeaturesList)
                    numFeatures = numFeaturesList(fIdx);
                    if size(trainData_corr, 2) < numFeatures
                        continue; % Skip if there are fewer features than numFeatures
                    end
                    if size(trainData_coh, 2) < numFeatures
                        continue; % Skip if there are fewer features than numFeatures
                    end
                    if size(trainData_plv, 2) < numFeatures
                        continue; % Skip if there are fewer features than numFeatures
                    end

                    % Select top numFeatures from training data for both coh
                    selectedFeaturesIdxTrain_corr1 = selectedFeaturesIdxTrain_corr(1:numFeatures);
                    selectedFeaturesIdxTrain_coh1 = selectedFeaturesIdxTrain_coh(1:numFeatures);
                    selectedFeaturesIdxTrain_plv1 = selectedFeaturesIdxTrain_plv(1:numFeatures);

                    % Select the features from training and testing data for both coh
                    selectedFeaturesTrain_corr = trainData_corr(:, selectedFeaturesIdxTrain_corr1);
                    selectedFeaturesTest_corr = testData_corr(:, selectedFeaturesIdxTrain_corr1);

                    selectedFeaturesTrain_coh = trainData_coh(:, selectedFeaturesIdxTrain_coh1);
                    selectedFeaturesTest_coh = testData_coh(:, selectedFeaturesIdxTrain_coh1);

                    selectedFeaturesTrain_plv = trainData_plv(:, selectedFeaturesIdxTrain_plv1);
                    selectedFeaturesTest_plv = testData_plv(:, selectedFeaturesIdxTrain_plv1);

                    % Save selected features for the current iteration
                    selectedFeaturesTrainSubject_corr{iteration} = selectedFeaturesTrain_corr;
                    selectedFeaturesTestSubject_corr{iteration} = selectedFeaturesTest_corr;

                    selectedFeaturesTrainSubject_coh{iteration} = selectedFeaturesTrain_coh;
                    selectedFeaturesTestSubject_coh{iteration} = selectedFeaturesTest_coh;

                    selectedFeaturesTrainSubject_plv{iteration} = selectedFeaturesTrain_plv;
                    selectedFeaturesTestSubject_plv{iteration} = selectedFeaturesTest_plv;

                    % Loop through each classifier
                    for clfIdx = 1:numel(classifiers)
                        switch classifiers{clfIdx}
                            case 'LDA'
                                % Perform LDA classification for coherence (coh)
                                LDA_Model_corr = fitcdiscr(selectedFeaturesTrain_corr, trainLabels_corr);
                                Predictions_corr = predict(LDA_Model_corr, selectedFeaturesTest_corr);
                                accuracy_corr = sum(Predictions_corr == testLabels_corr) / numel(testLabels_corr);
                                localAccuracies_lda_corr(c, fold, iteration, fIdx) = accuracy_corr; % Store accuracy for corr

                                LDA_Model_coh = fitcdiscr(selectedFeaturesTrain_coh, trainLabels_coh);
                                Predictions_coh = predict(LDA_Model_coh, selectedFeaturesTest_coh);
                                accuracy_coh = sum(Predictions_coh == testLabels_coh) / numel(testLabels_coh);
                                localAccuracies_lda_coh(c, fold, iteration, fIdx) = accuracy_coh; % Store accuracy for coh

                                LDA_Model_plv = fitcdiscr(selectedFeaturesTrain_plv, trainLabels_plv);
                                Predictions_plv = predict(LDA_Model_plv, selectedFeaturesTest_plv);
                                accuracy_plv = sum(Predictions_plv == testLabels_plv) / numel(testLabels_plv);
                                localAccuracies_lda_plv(c, fold, iteration, fIdx) = accuracy_plv; % Store accuracy for plv
                            case 'SVM'
                                % Perform SVM classification for coherence (coh)
                                SVM_Model_corr = fitcsvm(selectedFeaturesTrain_corr, trainLabels_corr, 'KernelFunction', 'linear');
                                Predictions_SVM_corr = predict(SVM_Model_corr, selectedFeaturesTest_corr);
                                accuracy_SVM_corr = sum(Predictions_SVM_corr == testLabels_corr) / numel(testLabels_corr);
                                localAccuracies_svm_corr(c, fold, iteration, fIdx) = accuracy_SVM_corr; % Store accuracy for corr

                                SVM_Model_coh = fitcsvm(selectedFeaturesTrain_coh, trainLabels_coh, 'KernelFunction', 'linear');
                                Predictions_SVM_coh = predict(SVM_Model_coh, selectedFeaturesTest_coh);
                                accuracy_SVM_coh = sum(Predictions_SVM_coh == testLabels_coh) / numel(testLabels_coh);
                                localAccuracies_svm_coh(c, fold, iteration, fIdx) = accuracy_SVM_coh; % Store accuracy for coh

                                SVM_Model_plv = fitcsvm(selectedFeaturesTrain_plv, trainLabels_plv, 'KernelFunction', 'linear');
                                Predictions_SVM_plv = predict(SVM_Model_plv, selectedFeaturesTest_plv);
                                accuracy_SVM_plv = sum(Predictions_SVM_plv == testLabels_plv) / numel(testLabels_plv);
                                localAccuracies_svm_plv(c, fold, iteration, fIdx) = accuracy_SVM_plv; % Store accuracy for plv
                            case 'DecisionTree'
                                % Perform Decision Tree classification for coherence (coh)
                                DT_Model_corr = fitctree(selectedFeaturesTrain_corr, trainLabels_corr);
                                Predictions_DT_corr = predict(DT_Model_corr, selectedFeaturesTest_corr);
                                accuracy_DT_corr = sum(Predictions_DT_corr == testLabels_corr) / numel(testLabels_corr);
                                localAccuracies_dt_corr(c, fold, iteration, fIdx) = accuracy_DT_corr; % Store accuracy for corr

                                DT_Model_coh = fitctree(selectedFeaturesTrain_coh, trainLabels_coh);
                                Predictions_DT_coh = predict(DT_Model_coh, selectedFeaturesTest_coh);
                                accuracy_DT_coh = sum(Predictions_DT_coh == testLabels_coh) / numel(testLabels_coh);
                                localAccuracies_dt_coh(c, fold, iteration, fIdx) = accuracy_DT_coh; % Store accuracy for coh

                                DT_Model_plv = fitctree(selectedFeaturesTrain_plv, trainLabels_plv);
                                Predictions_DT_plv = predict(DT_Model_plv, selectedFeaturesTest_plv);
                                accuracy_DT_plv = sum(Predictions_DT_plv == testLabels_plv) / numel(testLabels_plv);
                                localAccuracies_dt_plv(c, fold, iteration, fIdx) = accuracy_DT_plv; % Store accuracy for plv
                        end
                    end
                end
            end
        end
    end


    % Compute the mean accuracy for each condition and each fIdx
    for c = 1:numel(classProbs)
        for fIdx = 1:length(numFeaturesList)
            % For coherence (coh)
            accuracies_coh_SVM = localAccuracies_svm_coh(c , :, :, fIdx);
            accuracies_coh_DT = localAccuracies_dt_coh(c , :, :, fIdx);
            accuracies_coh_LDA = localAccuracies_lda_coh(c , :, :, fIdx);

            % Compute the mean accuracy over folds and iterations for SVM and Decision Tree
            meanAccuracies_coh_SVM(c, fIdx) = mean(accuracies_coh_SVM(:)) * 100;
            meanAccuracies_coh_DT(c, fIdx) = mean(accuracies_coh_DT(:)) * 100;
            meanAccuracies_coh_LDA(c, fIdx) = mean(accuracies_coh_LDA(:)) * 100;

            % For PLV (plv)
            accuracies_plv_SVM = localAccuracies_svm_plv(c , :, :, fIdx);
            accuracies_plv_DT = localAccuracies_dt_plv(c , :, :, fIdx);
            accuracies_plv_LDA = localAccuracies_lda_plv(c , :, :, fIdx);

            % Compute the mean accuracy over folds and iterations for SVM and Decision Tree
            meanAccuracies_plv_SVM(c, fIdx) = mean(accuracies_plv_SVM(:)) * 100;
            meanAccuracies_plv_DT(c, fIdx) = mean(accuracies_plv_DT(:)) * 100;
            meanAccuracies_plv_LDA(c, fIdx) = mean(accuracies_plv_LDA(:)) * 100;

            % For Correlation (corr)
            accuracies_corr_SVM = localAccuracies_svm_corr(c , :, :, fIdx);
            accuracies_corr_DT = localAccuracies_dt_corr(c , :, :, fIdx);
            accuracies_corr_LDA = localAccuracies_lda_corr(c , :, :, fIdx);

            % Compute the mean accuracy over folds and iterations for SVM and Decision Tree
            meanAccuracies_corr_SVM(c, fIdx) = mean(accuracies_corr_SVM(:)) * 100;
            meanAccuracies_corr_DT(c, fIdx) = mean(accuracies_corr_DT(:)) * 100;
            meanAccuracies_corr_LDA(c, fIdx) = mean(accuracies_corr_LDA(:)) * 100;
        end
    end


    % Define the save path
    savePath = fullfile(myPath, sprintf('Subject-%03d', numSubjects));
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end


    % Save the top 50 selected feature indices for each method
    save(fullfile(savePath, sprintf('topSelectedFeatureIndices_corr_5b_Sub-%03d.mat', numSubjects)), 'selectedTopFeatureIndices_corr');
    save(fullfile(savePath, sprintf('topSelectedFeatureIndices_coh_5b_Sub-%03d.mat', numSubjects)), 'selectedTopFeatureIndices_coh');
    save(fullfile(savePath, sprintf('topSelectedFeatureIndices_plv_5b_Sub-%03d.mat', numSubjects)), 'selectedTopFeatureIndices_plv');

    fprintf('Saved top 50 selected feature indices for Subject %d (corr, coh, plv).\n', numSubjects);

    % Save the full list of selected feature indices
    save(fullfile(savePath, sprintf('selectedFeatureIndices_corr_5b_Sub-%03d.mat', numSubjects)), 'selectedFeatureIndices_corr');
    save(fullfile(savePath, sprintf('selectedFeatureIndices_coh_5b_Sub-%03d.mat', numSubjects)), 'selectedFeatureIndices_coh');
    save(fullfile(savePath, sprintf('selectedFeatureIndices_plv_5b_Sub-%03d.mat', numSubjects)), 'selectedFeatureIndices_plv');

    fprintf('Saved all selected feature indices for Subject %d (corr, coh, plv).\n', numSubjects);

    %% Save Local Accuracies (LDA, SVM, DT)
    save(fullfile(savePath, sprintf('localAccuracies_lda_corr_5b_Sub-%03d.mat', numSubjects)), 'localAccuracies_lda_corr');
    save(fullfile(savePath, sprintf('localAccuracies_lda_coh_5b_Sub-%03d.mat', numSubjects)), 'localAccuracies_lda_coh');
    save(fullfile(savePath, sprintf('localAccuracies_lda_plv_5b_Sub-%03d.mat', numSubjects)), 'localAccuracies_lda_plv');

    save(fullfile(savePath, sprintf('localAccuracies_svm_corr_5b_Sub-%03d.mat', numSubjects)), 'localAccuracies_svm_corr');
    save(fullfile(savePath, sprintf('localAccuracies_svm_coh_5b_Sub-%03d.mat', numSubjects)), 'localAccuracies_svm_coh');
    save(fullfile(savePath, sprintf('localAccuracies_svm_plv_5b_Sub-%03d.mat', numSubjects)), 'localAccuracies_svm_plv');

    save(fullfile(savePath, sprintf('localAccuracies_dt_corr_5b_Sub-%03d.mat', numSubjects)), 'localAccuracies_dt_corr');
    save(fullfile(savePath, sprintf('localAccuracies_dt_coh_5b_Sub-%03d.mat', numSubjects)), 'localAccuracies_dt_coh');
    save(fullfile(savePath, sprintf('localAccuracies_dt_plv_5b_Sub-%03d.mat', numSubjects)), 'localAccuracies_dt_plv');

    fprintf('Saved local accuracies for Subject %d.\n', numSubjects);

    %% Save Mean Accuracies (SVM, DT, LDA) for coh, plv, corr
    if exist('meanAccuracies_coh_SVM', 'var') && exist('meanAccuracies_coh_DT', 'var') && exist('meanAccuracies_coh_LDA', 'var')
        save(fullfile(savePath, sprintf('FC_Accuracies_coh_5b_Sub-%03d.mat', numSubjects)), 'meanAccuracies_coh_SVM', 'meanAccuracies_coh_DT', 'meanAccuracies_coh_LDA');
        fprintf('Saved mean accuracies for Subject %d (coh).\n', numSubjects);
    else
        fprintf('Warning: Mean accuracies for coh are missing for Subject %d.\n', numSubjects);
    end

    if exist('meanAccuracies_plv_SVM', 'var') && exist('meanAccuracies_plv_DT', 'var') && exist('meanAccuracies_plv_LDA', 'var')
        save(fullfile(savePath, sprintf('FC_Accuracies_plv_5b_Sub-%03d.mat', numSubjects)), 'meanAccuracies_plv_SVM', 'meanAccuracies_plv_DT', 'meanAccuracies_plv_LDA');
        fprintf('Saved mean accuracies for Subject %d (plv).\n', numSubjects);
    else
        fprintf('Warning: Mean accuracies for plv are missing for Subject %d.\n', numSubjects);
    end

    if exist('meanAccuracies_corr_SVM', 'var') && exist('meanAccuracies_corr_DT', 'var') && exist('meanAccuracies_corr_LDA', 'var')
        save(fullfile(savePath, sprintf('FC_Accuracies_corr_5b_Sub-%03d.mat', numSubjects)), 'meanAccuracies_corr_SVM', 'meanAccuracies_corr_DT', 'meanAccuracies_corr_LDA');
        fprintf('Saved mean accuracies for Subject %d (corr).\n', numSubjects);
    else
        fprintf('Warning: Mean accuracies for corr are missing for Subject %d.\n', numSubjects);
    end

    % Display completion message, elapsed time, and memory usage
    elapsedTime = toc; % Stop the timer
    fprintf('Subject %d completed in %.2f seconds.\n', numSubjects, elapsedTime);

end
