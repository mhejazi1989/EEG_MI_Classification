
clear; close all; clc;

%% Define Parameters
totalSubjects       = 4;
excludeSubjects     = [6, 8, 10, 13, 21, 23, 25, 27];
nConditions = 6;
classifiers         = {'LDA', 'DT', 'SVM'};
connectivity_measures = {'coh', 'corr', 'plv'};
bands               = {'5b', '9b'};
electrodes          = {'64', '22'};
conditions          = 1:nConditions;
feature_number_5b   = {'10','20','30','40','50'};
feature_number_9b   = {'10','20','30','40','90'};

%% Dynamic variable creation function (using eval)
createVars = @(prefix, cond, clf, elec, band, feat) ...
    sprintf('%s_condition%d_%s_CSP_%s_%s_feat%s', prefix, cond, clf, elec, band, feat);


%% Create data containers for connectivity measures (with feature numbers)
for elec = electrodes
    for band = bands
        features = getFeatures(band{1});
        for conn = connectivity_measures
            for feat = features
                for cond = conditions
                    % Create container for data
                    eval(sprintf('%s_CSP_%s_%s_feat%s_condition%d = [];', 'dummy', elec{1}, band{1}, feat{1}, cond)); %#ok<NASGU>
                end
            end
        end
    end
end

%% Create mean accuracy containers for connectivity measures (with feature numbers)
for elec = electrodes
    for band = bands
        features = getFeatures(band{1});
        for conn = connectivity_measures
            for feat = features
                for cond = conditions
                    eval(sprintf('meanAccuracy_condition%d_%s_CSP_%s_%s_feat%s = [];', ...
                        cond, 'dummy', elec{1}, band{1}, feat{1})); %#ok<NASGU>
                end
            end
        end
    end
end

%% Initialize subject numbers container if needed
subjectNumbers = [];

%% Loop over subjects to load and extract SVM, DT, and LDA data
for iF = 1:length(feature_number_5b)
    for cond = 1:6
        eval(sprintf('SVM_CSP_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_5b{iF}, cond));
        eval(sprintf('DT_CSP_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_5b{iF}, cond));
        eval(sprintf('LDA_CSP_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_5b{iF}, cond));
    end
end
%% Loop over subjects to load and extract SVM, DT, and LDA data
for iF = 1:length(feature_number_5b)
    for cond = 1:6
        eval(sprintf('SVM_CSP_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_5b{iF}, cond));
        eval(sprintf('DT_CSP_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_5b{iF}, cond));
        eval(sprintf('LDA_CSP_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_5b{iF}, cond));
    end
end

%% Loop over subjects to load and extract SVM, DT, and LDA data
for iF = 1:length(feature_number_5b)
    for cond = 1:6
        eval(sprintf('SVM_CSP_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_9b{iF}, cond));
        eval(sprintf('DT_CSP_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_9b{iF}, cond));
        eval(sprintf('LDA_CSP_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_9b{iF}, cond));
    end
end
%% Loop over subjects to load and extract SVM, DT, and LDA data
for iF = 1:length(feature_number_9b)
    for cond = 1:6
        eval(sprintf('SVM_CSP_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_9b{iF}, cond));
        eval(sprintf('DT_CSP_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_9b{iF}, cond));
        eval(sprintf('LDA_CSP_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feature_number_9b{iF}, cond));
    end
end


for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile('/lustre03/project/6067835/mhejazi/EEG-DATA', subjectFolder);

    % Define the Compute Canada save path
    savePath = fullfile('/lustre03/project/6067835/mhejazi/EEG-DATA');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CSP 64 channnel 5 bands
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define file names (only SVM, DT, LDA files needed here)
    fileSVM = fullfile(myPath, sprintf('CSP_Accuracies_svm_5b_Sub-%03d.mat', numSub));
    fileDT  = fullfile(myPath, sprintf('CSP_Accuracies_dt_5b_Sub-%03d.mat', numSub));
    fileLDA = fullfile(myPath, sprintf('CSP_Accuracies_lda_5b_Sub-%03d.mat', numSub));

    % Load SVM data
    if exist(fileSVM, 'file') == 2
        load(fileSVM, 'meanAccuracies_svm');
        if exist('meanAccuracies_svm', 'var') && ~isempty(meanAccuracies_svm)
            tempSVM = cell2mat(meanAccuracies_svm); % Expected: [5 x 5]
        else
            fprintf('Subject %d file not found or is empty: %s\n', numSub, fileSVM);
            continue;
        end
    else
        fprintf('Subject %d file not found: %s\n', numSub, fileSVM);
        continue;
    end

    % Load DT data
    if exist(fileDT, 'file') == 2
        load(fileDT, 'meanAccuracies_dt');
        if exist('meanAccuracies_dt', 'var') && ~isempty(meanAccuracies_dt)
            tempDT = cell2mat(meanAccuracies_dt); % Expected: [6 x 5]
        else
            fprintf('Subject %d file not found or is empty: %s\n', numSub, fileDT);
            continue;
        end
    else
        fprintf('Subject %d file not found: %s\n', numSub, fileDT);
        continue;
    end

    % Load LDA data
    tempLDA = [];
    if exist(fileLDA, 'file') == 2
        load(fileLDA, 'meanAccuracies_lda');
        if exist('meanAccuracies_lda', 'var') && ~isempty(meanAccuracies_lda)
            tempLDA = cell2mat(meanAccuracies_lda); % Expected: [6 x 5]
        end
    else
        fprintf('Subject %d file not found: %s\n', numSub, fileLDA);
    end

    % Store SVM, DT, and LDA data in containers
    if ~isempty(tempSVM) && size(tempSVM,1) >= nConditions && size(tempSVM,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('SVM_CSP_64_5b_feat%s_condition%d{numSub} = tempSVM(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT) && size(tempDT,1) >= nConditions && size(tempDT,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('DT_CSP_64_5b_feat%s_condition%d{numSub} = tempDT(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempLDA) && size(tempLDA,1) >= nConditions && size(tempLDA,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('LDA_CSP_64_5b_feat%s_condition%d{numSub} = tempLDA(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_5b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_SVM_CSP_64_5b_feat%s_condition%d = mean(cell2mat(SVM_CSP_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_DT_CSP_64_5b_feat%s_condition%d = mean(cell2mat(DT_CSP_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_LDA_CSP_64_5b_feat%s_condition%d = mean(cell2mat(LDA_CSP_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
    end
end

%% Display the computed means for each feature
fprintf('Computed Mean Accuracies for SVM, DT, and LDA CSP 64 5b:\n');
for iF = 1:length(feature_number_5b)
    fprintf('\nFeature %s:\n', feature_number_5b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for SVM, DT, and LDA
        mean_SVM = eval(sprintf('mean_SVM_CSP_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_CSP_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_LDA = eval(sprintf('mean_LDA_CSP_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        fprintf('  Condition %d: SVM %.2f, DT %.2f, LDA %.2f\n', cond, mean_SVM, mean_DT, mean_LDA);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CSP 22 channel 5 bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Define file names (only SVM, DT, LDA files needed here)
    fileSVM = fullfile(myPath, sprintf('CSP_Accuracies_svm_5b_22ch_Sub-%03d.mat', numSub));
    fileDT  = fullfile(myPath, sprintf('CSP_Accuracies_dt_5b_22ch_Sub-%03d.mat', numSub));
    fileLDA = fullfile(myPath, sprintf('CSP_Accuracies_lda_5b_22ch_Sub-%03d.mat', numSub));

    % Load SVM data
    if exist(fileSVM, 'file') == 2
        load(fileSVM, 'meanAccuracies_svm');
        if exist('meanAccuracies_svm', 'var') && ~isempty(meanAccuracies_svm)
            tempSVM = cell2mat(meanAccuracies_svm); % Expected: [5 x 5]
        else
            fprintf('Subject %d file not found or is empty: %s\n', numSub, fileSVM);
            continue;
        end
    else
        fprintf('Subject %d file not found: %s\n', numSub, fileSVM);
        continue;
    end

    % Load DT data
    if exist(fileDT, 'file') == 2
        load(fileDT, 'meanAccuracies_dt');
        if exist('meanAccuracies_dt', 'var') && ~isempty(meanAccuracies_dt)
            tempDT = cell2mat(meanAccuracies_dt); % Expected: [6 x 5]
        else
            fprintf('Subject %d file not found or is empty: %s\n', numSub, fileDT);
            continue;
        end
    else
        fprintf('Subject %d file not found: %s\n', numSub, fileDT);
        continue;
    end

    % Load LDA data
    tempLDA = [];
    if exist(fileLDA, 'file') == 2
        load(fileLDA, 'meanAccuracies_lda');
        if exist('meanAccuracies_lda', 'var') && ~isempty(meanAccuracies_lda)
            tempLDA = cell2mat(meanAccuracies_lda); % Expected: [6 x 5]
        end
    else
        fprintf('Subject %d file not found: %s\n', numSub, fileLDA);
    end

    % Store SVM, DT, and LDA data in containers
    if ~isempty(tempSVM) && size(tempSVM,1) >= nConditions && size(tempSVM,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('SVM_CSP_22_5b_feat%s_condition%d{numSub} = tempSVM(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT) && size(tempDT,1) >= nConditions && size(tempDT,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('DT_CSP_22_5b_feat%s_condition%d{numSub} = tempDT(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempLDA) && size(tempLDA,1) >= nConditions && size(tempLDA,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('LDA_CSP_22_5b_feat%s_condition%d{numSub} = tempLDA(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_5b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_SVM_CSP_22_5b_feat%s_condition%d = mean(cell2mat(SVM_CSP_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_DT_CSP_22_5b_feat%s_condition%d = mean(cell2mat(DT_CSP_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_LDA_CSP_22_5b_feat%s_condition%d = mean(cell2mat(LDA_CSP_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
    end
end

%% Display the computed means for each feature
fprintf('Computed Mean Accuracies for SVM, DT, and LDA CSP 22 5b:\n');
for iF = 1:length(feature_number_5b)
    fprintf('\nFeature %s:\n', feature_number_5b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for SVM, DT, and LDA
        mean_SVM = eval(sprintf('mean_SVM_CSP_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_CSP_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_LDA = eval(sprintf('mean_LDA_CSP_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        fprintf('  Condition %d: SVM %.2f, DT %.2f, LDA %.2f\n', cond, mean_SVM, mean_DT, mean_LDA);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CSP 64 channnel 9 bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % CSP 64 channel 9 bands
    % Define file names (only SVM, DT, LDA files needed here)
    fileSVM_64 = fullfile(myPath, sprintf('CSP_Accuracies_svm_9b_Sub-%03d.mat', numSub));
    fileDT_64  = fullfile(myPath, sprintf('CSP_Accuracies_dt_9b_Sub-%03d.mat', numSub));
    fileLDA_64 = fullfile(myPath, sprintf('CSP_Accuracies_lda_9b_Sub-%03d.mat', numSub));

    % Initialize containers for CSP 64 channel data
    tempSVM_64 = []; tempDT_64 = []; tempLDA_64 = [];

    % Load CSP 64 channel SVM data
    if exist(fileSVM_64, 'file') == 2
        load(fileSVM_64, 'meanAccuracies_svm');
        if exist('meanAccuracies_svm', 'var') && ~isempty(meanAccuracies_svm)
            tempSVM_64 = cell2mat(meanAccuracies_svm); % Expected: [5 x 5]
        else
            fprintf('Subject %d file not found or is empty: %s\n', numSub, fileSVM_6nConditions);
        end
    else
        fprintf('Subject %d file not found: %s\n', numSub, fileSVM_64);
    end

    % Load CSP 64 channel DT data
    if exist(fileDT_64, 'file') == 2
        load(fileDT_64, 'meanAccuracies_dt');
        if exist('meanAccuracies_dt', 'var') && ~isempty(meanAccuracies_dt)
            tempDT_64 = cell2mat(meanAccuracies_dt); % Expected: [6 x 5]
        else
            fprintf('Subject %d file not found or is empty: %s\n', numSub, fileDT_64);
        end
    else
        fprintf('Subject %d file not found: %s\n', numSub, fileDT_64);
    end

    % Load CSP 64 channel LDA data
    if exist(fileLDA_64, 'file') == 2
        load(fileLDA_64, 'meanAccuracies_lda');
        if exist('meanAccuracies_lda', 'var') && ~isempty(meanAccuracies_lda)
            tempLDA_64 = cell2mat(meanAccuracies_lda); % Expected: [6 x 5]
        end
    else
        fprintf('Subject %d file not found: %s\n', numSub, fileLDA_64);
    end

    % Process and store CSP 64 channel data if valid
    if ~isempty(tempSVM_64) && size(tempSVM_64, 1) >= nConditions && size(tempSVM_64, 2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('SVM_CSP_64_9b_feat%s_condition%d{numSub} = tempSVM_64(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
    if ~isempty(tempDT_64) && size(tempDT_64, 1) >= nConditions && size(tempDT_64, 2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('DT_CSP_64_9b_feat%s_condition%d{numSub} = tempDT_64(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
    if ~isempty(tempLDA_64) && size(tempLDA_64, 1) >= nConditions && size(tempLDA_64, 2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('LDA_CSP_64_9b_feat%s_condition%d{numSub} = tempLDA_64(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    % CSP 22 channel 9 bands
    meanAccuraciesFileSVM_22_9b = fullfile(myPath, sprintf('CSP_Accuracies_svm_9b_22ch_Sub-%03d.mat', numSub));
    meanAccuraciesFileDT_22_9b = fullfile(myPath, sprintf('CSP_Accuracies_dt_9b_22ch_Sub-%03d.mat', numSub));
    meanAccuraciesFileLDA_22_9b = fullfile(myPath, sprintf('CSP_Accuracies_lda_9b_22ch_Sub-%03d.mat', numSub));

    % Initialize containers for CSP 22 channel data
    tempSVM_22 = []; tempDT_22 = []; tempLDA_22 = [];

    % Load CSP 22 channel SVM data
    if exist(meanAccuraciesFileSVM_22_9b, 'file') == 2
        load(meanAccuraciesFileSVM_22_9b, 'meanAccuracies_svm');
        if exist('meanAccuracies_svm', 'var') && ~isempty(meanAccuracies_svm)
            tempSVM_22 = cell2mat(meanAccuracies_svm); % Expected: [6 x 5]
        end
    end

    % Load CSP 22 channel DT data
    if exist(meanAccuraciesFileDT_22_9b, 'file') == 2
        load(meanAccuraciesFileDT_22_9b, 'meanAccuracies_dt');
        if exist('meanAccuracies_dt', 'var') && ~isempty(meanAccuracies_dt)
            tempDT_22 = cell2mat(meanAccuracies_dt); % Expected: [6 x 5]
        end
    end

    % Load CSP 22 channel LDA data
    if exist(meanAccuraciesFileLDA_22_9b, 'file') == 2
        load(meanAccuraciesFileLDA_22_9b, 'meanAccuracies_lda');
        if exist('meanAccuracies_lda', 'var') && ~isempty(meanAccuracies_lda)
            tempLDA_22 = cell2mat(meanAccuracies_lda); % Expected: [6 x 5]
        end
    end

    % Process and store CSP 22 channel data if valid
    if ~isempty(tempSVM_22) && size(tempSVM_22, 1) >= nConditions && size(tempSVM_22, 2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('SVM_CSP_22_9b_feat%s_condition%d{numSub} = tempSVM_22(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
    if ~isempty(tempDT_22) && size(tempDT_22, 1) >= nConditions && size(tempDT_22, 2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('DT_CSP_22_9b_feat%s_condition%d{numSub} = tempDT_22(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
    if ~isempty(tempLDA_22) && size(tempLDA_22, 1) >= nConditions && size(tempLDA_22, 2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('LDA_CSP_22_9b_feat%s_condition%d{numSub} = tempLDA_22(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
end

% Compute mean accuracy across subjects for 9b features & conditions (64 and 22 channels)
for iF = 1:length(feature_number_9b)
    for cond = 1:nConditions
        % CSP 64 channel mean computation
        eval(sprintf('mean_SVM_CSP_64_9b_feat%s_condition%d = mean(cell2mat(SVM_CSP_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_DT_CSP_64_9b_feat%s_condition%d = mean(cell2mat(DT_CSP_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_LDA_CSP_64_9b_feat%s_condition%d = mean(cell2mat(LDA_CSP_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));

        % CSP 22 channel mean computation
        eval(sprintf('mean_SVM_CSP_22_9b_feat%s_condition%d = mean(cell2mat(SVM_CSP_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_DT_CSP_22_9b_feat%s_condition%d = mean(cell2mat(DT_CSP_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_LDA_CSP_22_9b_feat%s_condition%d = mean(cell2mat(LDA_CSP_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
    end
end

fprintf('Computed Mean Accuracies for SVM, DT, and LDA CSP 64 9b:\n');
for iF = 1:length(feature_number_9b)
    fprintf('\nFeature %s:\n', feature_number_9b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for SVM, DT, and LDA
        mean_SVM = eval(sprintf('mean_SVM_CSP_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_CSP_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_LDA = eval(sprintf('mean_LDA_CSP_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        fprintf('  Condition %d: SVM %.2f, DT %.2f, LDA %.2f\n', cond, mean_SVM, mean_DT, mean_LDA);
    end
end

fprintf('Computed Mean Accuracies for SVM, DT, and LDA CSP 22 9b:\n');
for iF = 1:length(feature_number_9b)
    fprintf('\nFeature %s:\n', feature_number_9b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for SVM, DT, and LDA
        mean_SVM = eval(sprintf('mean_SVM_CSP_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_CSP_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_LDA = eval(sprintf('mean_LDA_CSP_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        fprintf('  Condition %d: SVM %.2f, DT %.2f, LDA %.2f\n', cond, mean_SVM, mean_DT, mean_LDA);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% OVERALL CSP BASED ON FEATURE NUMBER, CLASSIFIERS, FREQUENCY BANDS AND CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%
% Define feature numbers for 5b and 9b
feature_number_5b = {'10', '20', '30', '40', '50'};
feature_number_9b = {'10', '20', '30', '40', '90'};

% Define electrodes 22 and 64
electrodes = {'22', '64'};

% Initialize containers to store mean values
meanSVM = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA = zeros(totalSubjects, length(feature_number_5b), nConditions);

overallMean_5b = zeros(totalSubjects, length(feature_number_5b));
overallMean_9b = zeros(totalSubjects, length(feature_number_9b));

overallMean = zeros(totalSubjects, 1);  % Final overall mean per subject
overallMean_SVM = zeros(totalSubjects, 1);
overallMean_DT = zeros(totalSubjects, 1);
overallMean_LDA = zeros(totalSubjects, 1);

% Initialize containers for storing results for electrodes 22 and 64
meanSVM_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);

meanSVM_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);

% Initialize variables to store final overall mean across all classifiers for each electrode
overallMean_22 = zeros(totalSubjects, 1);
overallMean_64 = zeros(totalSubjects, 1);

% Loop through subjects
for numSub = 1:totalSubjects
    % Initialize accumulators for all classifiers
    svmAll = [];
    dtAll = [];
    ldaAll = [];
    svmAll_22 = [];
    dtAll_22 = [];
    ldaAll_22 = [];

    svmAll_64 = [];
    dtAll_64 = [];
    ldaAll_64 = [];
    % Loop through different configurations
    configIndex = 0;
    for el = electrodes
        for band = bands
            configIndex = configIndex + 1;
            feature_number = eval(['feature_number_', band{1}]);

            % Check if the 9b band is being used
            if strcmp(band{1}, '9b') && isempty(feature_number_9b)
                continue;  % Skip 9b calculations if there are no features for 9b
            end
            for iF = 1:length(feature_number)
                % Loop through conditions
                for cond = 1:nConditions
                    % Dynamically construct variable names
                    svmVar = sprintf('SVM_CSP_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);
                    dtVar = sprintf('DT_CSP_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);
                    ldaVar = sprintf('LDA_CSP_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);

                    % Extract values only if they are not empty
                    svmValue = eval(svmVar);
                    dtValue = eval(dtVar);
                    ldaValue = eval(ldaVar);

                    if ~isempty(svmValue{numSub})
                        svmAll = [svmAll, svmValue{numSub}];
                    end
                    if ~isempty(dtValue{numSub})
                        dtAll = [dtAll, dtValue{numSub}];
                    end
                    if ~isempty(ldaValue{numSub})
                        ldaAll = [ldaAll, ldaValue{numSub}];
                    end
                    if ~isempty(svmValue{numSub})
                        if strcmp(el{1}, '22')
                            svmAll_22 = [svmAll_22, svmValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            svmAll_64 = [svmAll_64, svmValue{numSub}];
                        end
                    end

                    if ~isempty(dtValue{numSub})
                        if strcmp(el{1}, '22')
                            dtAll_22 = [dtAll_22, dtValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            dtAll_64 = [dtAll_64, dtValue{numSub}];
                        end
                    end

                    if ~isempty(ldaValue{numSub})
                        if strcmp(el{1}, '22')
                            ldaAll_22 = [ldaAll_22, ldaValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            ldaAll_64 = [ldaAll_64, ldaValue{numSub}];
                        end
                    end
                    % Compute mean for each classifier per feature number, ignoring NaN
                    meanSVM(numSub, iF, cond) = mean(svmAll, 'omitnan');
                    meanDT(numSub, iF, cond) = mean(dtAll, 'omitnan');
                    meanLDA(numSub, iF, cond) = mean(ldaAll, 'omitnan');
                end
                % Compute overall mean per feature number for each subject
                allAccuracies = [svmAll, dtAll, ldaAll];
                if strcmp(band{1}, '5b')
                    overallMean_5b(numSub, iF) = mean(allAccuracies, 'omitnan');
                else
                    overallMean_9b(numSub, iF) = mean(allAccuracies, 'omitnan');
                end
            end

        end
    end

    % Compute overall mean across all configurations, classifiers, and feature indices
    overallMean(numSub) = mean([meanSVM(numSub, :, :), meanDT(numSub, :, :), meanLDA(numSub, :, :)], 'all', 'omitnan');
    % Compute the overall CSP mean per classifier
    overallMean_SVM(numSub) = mean(svmAll, 'omitnan');
    overallMean_DT(numSub) = mean(dtAll, 'omitnan');
    overallMean_LDA(numSub) = mean(ldaAll, 'omitnan');


    % Calculate the mean values for each classifier (SVM, DT, LDA) per feature number for electrode 22 and 64
    for iF = 1:length(feature_number_5b)
        for cond = 1:nConditions
            meanSVM_22(numSub, iF, cond) = mean(svmAll_22, 'omitnan');
            meanDT_22(numSub, iF, cond) = mean(dtAll_22, 'omitnan');
            meanLDA_22(numSub, iF, cond) = mean(ldaAll_22, 'omitnan');

            meanSVM_64(numSub, iF, cond) = mean(svmAll_64, 'omitnan');
            meanDT_64(numSub, iF, cond) = mean(dtAll_64, 'omitnan');
            meanLDA_64(numSub, iF, cond) = mean(ldaAll_64, 'omitnan');
        end
    end

    % Compute the overall mean across all classifiers for each electrode (22 and 64)
    overallMean_22(numSub) = mean([meanSVM_22(numSub, :, :), meanDT_22(numSub, :, :), meanLDA_22(numSub, :, :)], 'all', 'omitnan');
    overallMean_64(numSub) = mean([meanSVM_64(numSub, :, :), meanDT_64(numSub, :, :), meanLDA_64(numSub, :, :)], 'all', 'omitnan');
end


% Display results for all subjects
fprintf('\n=== Summary of Mean Values for All Subjects (CSP) ===\n');
for numSub = 1:totalSubjects
    fprintf('Subject %d:\n', numSub);
    fprintf('  Overall Mean - CSP: %.4f\n\n', overallMean(numSub));
    fprintf('  Overall Mean - SVM: %.4f\n', overallMean_SVM(numSub));
    fprintf('  Overall Mean - DT: %.4f\n', overallMean_DT(numSub));
    fprintf('  Overall Mean - LDA: %.4f\n', overallMean_LDA(numSub));
    % Print overall mean per feature for 5b
    for iF = 1:length(feature_number_5b)
        fprintf('  Overall Mean - 5b (Feature %s): %.4f\n', feature_number_5b{iF}, overallMean_5b(numSub, iF));
    end

    % Print overall mean per feature for 9b
    for iF = 1:length(feature_number_9b)
        fprintf('  Overall Mean - 9b (Feature %s): %.4f\n', feature_number_9b{iF}, overallMean_9b(numSub, iF));
    end
    % Print overall mean for electrode 22
    fprintf('  Overall Mean - Electrode 22: %.4f\n', overallMean_22(numSub));

    % Print overall mean for electrode 64
    fprintf('  Overall Mean - Electrode 64: %.4f\n', overallMean_64(numSub));

    fprintf('\n');  % Print a newline after each subject's results
end
% Saving to a specific path
save(fullfile(myPath, 'overal_CSP.mat'), 'overallMean', 'meanSVM', 'meanDT', 'meanLDA');
save(fullfile(myPath, 'overall_CSP_classifiers.mat'), 'overallMean_SVM', 'overallMean_DT', 'overallMean_LDA');
save(fullfile(myPath, 'overall_CSP_features.mat'), ...
    'overallMean_5b', 'overallMean_9b');
save(fullfile(myPath, 'overall_CSP_channels.mat'), 'overallMean_22', 'overallMean_64');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COHERENCE 64 channel 5 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for Coherence 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA coherence measures
            eval(sprintf('SVM_COH_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_COH_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_COH_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for Coherence 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_COH_64_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_COH_64_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_COH_64_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for Coherence
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load coherence data for LDA
    meanAccuraciesFileLDA_coh_64_5b = fullfile(myPath, sprintf('meanAccuracies_coh_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_coh_64_5b, 'file') == 2
        load(meanAccuraciesFileLDA_coh_64_5b, 'meanAccuracies_coh');
        if exist('meanAccuracies_coh', 'var') && ~isempty(meanAccuracies_coh)
            tempLDA_coh = meanAccuracies_coh; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_coh_64_5b);
            continue;
        end
    else
        fprintf('Subject %d LDA coherence file not found: %s\n', numSub, meanAccuraciesFileLDA_coh_64_5b);
        continue;
    end

    % Load coherence data for DT
    meanAccuraciesFileDT_coh_64_5b = fullfile(myPath, sprintf('meanAccuracies_coh_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_coh_64_5b, 'file') == 2
        load(meanAccuraciesFileDT_coh_64_5b, 'meanAccuracies_coh_DT');
        if exist('meanAccuracies_coh_DT', 'var') && ~isempty(meanAccuracies_coh_DT)
            tempDT_coh = meanAccuracies_coh_DT; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_coh_64_5b);
            continue;
        end
    else
        fprintf('Subject %d DT coherence file not found: %s\n', numSub, meanAccuraciesFileDT_coh_64_5b);
        continue;
    end

    % Load coherence data for SVM
    meanAccuraciesFileSVM_coh_64_5b = fullfile(myPath, sprintf('meanAccuracies_coh_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_coh_64_5b, 'file') == 2
        load(meanAccuraciesFileSVM_coh_64_5b, 'meanAccuracies_coh_SVM');
        if exist('meanAccuracies_coh_SVM', 'var') && ~isempty(meanAccuracies_coh_SVM)
            tempSVM_coh = meanAccuracies_coh_SVM; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_coh_64_5b);
            continue;
        end
    else
        fprintf('Subject %d SVM coherence file not found: %s\n', numSub, meanAccuraciesFileSVM_coh_64_5b);
        continue;
    end

    % Store coherence data for LDA, DT, and SVM
    if ~isempty(tempLDA_coh) && size(tempLDA_coh,1) >= nConditions && size(tempLDA_coh,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('LDA_COH_64_5b_feat%s_condition%d{numSub} = tempLDA_coh(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_coh) && size(tempDT_coh,1) >= nConditions && size(tempDT_coh,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('DT_COH_64_5b_feat%s_condition%d{numSub} = tempDT_coh(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_coh) && size(tempSVM_coh,1) >= nConditions && size(tempSVM_coh,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('SVM_COH_64_5b_feat%s_condition%d{numSub} = tempSVM_coh(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean coherence accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_5b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_COH_64_5b_feat%s_condition%d = mean(cell2mat(LDA_COH_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_DT_COH_64_5b_feat%s_condition%d = mean(cell2mat(DT_COH_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_SVM_COH_64_5b_feat%s_condition%d = mean(cell2mat(SVM_COH_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
    end
end

%% Display the computed means for coherence features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM COHERENCE 64 5b:\n');
for iF = 1:length(feature_number_5b)
    fprintf('\nFeature %s:\n', feature_number_5b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_COH_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_COH_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_COH_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COHERENCE 22 channel 5 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for Coherence 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA coherence measures
            eval(sprintf('SVM_COH_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_COH_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_COH_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for Coherence 22-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_COH_22_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_COH_22_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_COH_22_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for Coherence
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load coherence data for LDA
    meanAccuraciesFileLDA_coh_22_5b = fullfile(myPath, sprintf('meanAccuracies_coh_5b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_coh_22_5b, 'file') == 2
        load(meanAccuraciesFileLDA_coh_22_5b, 'meanAccuraciesLDA_coh');
        if exist('meanAccuraciesLDA_coh', 'var') && ~isempty(meanAccuraciesLDA_coh)
            tempLDA_coh = meanAccuraciesLDA_coh; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_coh_22_5b);
            continue;
        end
    else
        fprintf('Subject %d LDA coherence file not found: %s\n', numSub, meanAccuraciesFileLDA_coh_22_5b);
        continue;
    end

    % Load coherence data for DT
    meanAccuraciesFileDT_coh_22_5b = fullfile(myPath, sprintf('meanAccuracies_coh_5b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_coh_22_5b, 'file') == 2
        load(meanAccuraciesFileDT_coh_22_5b, 'meanAccuraciesDT_coh');
        if exist('meanAccuraciesDT_coh', 'var') && ~isempty(meanAccuraciesDT_coh)
            tempDT_coh = meanAccuraciesDT_coh; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_coh_22_5b);
            continue;
        end
    else
        fprintf('Subject %d DT coherence file not found: %s\n', numSub, meanAccuraciesFileDT_coh_22_5b);
        continue;
    end

    % Load coherence data for SVM
    meanAccuraciesFileSVM_coh_22_5b = fullfile(myPath, sprintf('meanAccuracies_coh_5b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_coh_22_5b, 'file') == 2
        load(meanAccuraciesFileSVM_coh_22_5b, 'meanAccuraciesSVM_coh');
        if exist('meanAccuraciesSVM_coh', 'var') && ~isempty(meanAccuraciesSVM_coh)
            tempSVM_coh = meanAccuraciesSVM_coh; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_coh_22_5b);
            continue;
        end
    else
        fprintf('Subject %d SVM coherence file not found: %s\n', numSub, meanAccuraciesFileSVM_coh_22_5b);
        continue;
    end

    % Store coherence data for LDA, DT, and SVM
    if ~isempty(tempLDA_coh) && size(tempLDA_coh,1) >= nConditions && size(tempLDA_coh,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('LDA_COH_22_5b_feat%s_condition%d{numSub} = tempLDA_coh(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_coh) && size(tempDT_coh,1) >= nConditions && size(tempDT_coh,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('DT_COH_22_5b_feat%s_condition%d{numSub} = tempDT_coh(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_coh) && size(tempSVM_coh,1) >= nConditions && size(tempSVM_coh,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('SVM_COH_22_5b_feat%s_condition%d{numSub} = tempSVM_coh(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean coherence accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_5b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_COH_22_5b_feat%s_condition%d = mean(cell2mat(LDA_COH_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_DT_COH_22_5b_feat%s_condition%d = mean(cell2mat(DT_COH_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_SVM_COH_22_5b_feat%s_condition%d = mean(cell2mat(SVM_COH_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
    end
end

%% Display the computed means for coherence features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM COHERENCE 22 5b:\n');
for iF = 1:length(feature_number_5b)
    fprintf('\nFeature %s:\n', feature_number_5b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_COH_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_COH_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_COH_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COHERENCE 64 channel 9 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for Coherence 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA coherence measures
            eval(sprintf('SVM_COH_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_COH_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_COH_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for Coherence 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_COH_64_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_COH_64_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_COH_64_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for Coherence
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load coherence data for LDA
    meanAccuraciesFileLDA_coh_64_9b = fullfile(myPath, sprintf('meanAccuracies_coh_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_coh_64_9b, 'file') == 2
        load(meanAccuraciesFileLDA_coh_64_9b, 'meanAccuracies_coh');
        if exist('meanAccuracies_coh', 'var') && ~isempty(meanAccuracies_coh)
            tempLDA_coh = meanAccuracies_coh; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_coh_64_9b);
            continue;
        end
    else
        fprintf('Subject %d LDA coherence file not found: %s\n', numSub, meanAccuraciesFileLDA_coh_64_9b);
        continue;
    end

    % Load coherence data for DT
    meanAccuraciesFileDT_coh_64_9b = fullfile(myPath, sprintf('meanAccuracies_coh_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_coh_64_9b, 'file') == 2
        load(meanAccuraciesFileDT_coh_64_9b, 'meanAccuracies_coh_DT');
        if exist('meanAccuracies_coh_DT', 'var') && ~isempty(meanAccuracies_coh_DT)
            tempDT_coh = meanAccuracies_coh_DT; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_coh_64_9b);
            continue;
        end
    else
        fprintf('Subject %d DT coherence file not found: %s\n', numSub, meanAccuraciesFileDT_coh_64_9b);
        continue;
    end

    % Load coherence data for SVM
    meanAccuraciesFileSVM_coh_64_9b = fullfile(myPath, sprintf('meanAccuracies_coh_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_coh_64_9b, 'file') == 2
        load(meanAccuraciesFileSVM_coh_64_9b, 'meanAccuracies_coh_SVM');
        if exist('meanAccuracies_coh_SVM', 'var') && ~isempty(meanAccuracies_coh_SVM)
            tempSVM_coh = meanAccuracies_coh_SVM; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_coh_64_9b);
            continue;
        end
    else
        fprintf('Subject %d SVM coherence file not found: %s\n', numSub, meanAccuraciesFileSVM_coh_64_9b);
        continue;
    end

    % Store coherence data for LDA, DT, and SVM
    if ~isempty(tempLDA_coh) && size(tempLDA_coh,1) >= nConditions && size(tempLDA_coh,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('LDA_COH_64_9b_feat%s_condition%d{numSub} = tempLDA_coh(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_coh) && size(tempDT_coh,1) >= nConditions && size(tempDT_coh,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('DT_COH_64_9b_feat%s_condition%d{numSub} = tempDT_coh(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_coh) && size(tempSVM_coh,1) >= nConditions && size(tempSVM_coh,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('SVM_COH_64_9b_feat%s_condition%d{numSub} = tempSVM_coh(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean coherence accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_9b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_COH_64_9b_feat%s_condition%d = mean(cell2mat(LDA_COH_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_DT_COH_64_9b_feat%s_condition%d = mean(cell2mat(DT_COH_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_SVM_COH_64_9b_feat%s_condition%d = mean(cell2mat(SVM_COH_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
    end
end

%% Display the computed means for coherence features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM COHERENCE 64 9b:\n');
for iF = 1:length(feature_number_9b)
    fprintf('\nFeature %s:\n', feature_number_9b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_COH_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_COH_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_COH_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COHERENCE 22 channel 9 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for Coherence 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA coherence measures
            eval(sprintf('SVM_COH_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_COH_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_COH_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for Coherence 22-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_COH_22_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_COH_22_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_COH_22_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for Coherence
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load coherence data for LDA
    meanAccuraciesFileLDA_coh_22_9b = fullfile(myPath, sprintf('meanAccuracies_coh_9b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_coh_22_9b, 'file') == 2
        load(meanAccuraciesFileLDA_coh_22_9b, 'meanAccuraciesLDA_coh');
        if exist('meanAccuraciesLDA_coh', 'var') && ~isempty(meanAccuraciesLDA_coh)
            tempLDA_coh = meanAccuraciesLDA_coh; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_coh_22_9b);
            continue;
        end
    else
        fprintf('Subject %d LDA coherence file not found: %s\n', numSub, meanAccuraciesFileLDA_coh_22_9b);
        continue;
    end

    % Load coherence data for DT
    meanAccuraciesFileDT_coh_22_9b = fullfile(myPath, sprintf('meanAccuracies_coh_9b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_coh_22_9b, 'file') == 2
        load(meanAccuraciesFileDT_coh_22_9b, 'meanAccuraciesDT_coh');
        if exist('meanAccuraciesDT_coh', 'var') && ~isempty(meanAccuraciesDT_coh)
            tempDT_coh = meanAccuraciesDT_coh; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_coh_22_9b);
            continue;
        end
    else
        fprintf('Subject %d DT coherence file not found: %s\n', numSub, meanAccuraciesFileDT_coh_22_9b);
        continue;
    end

    % Load coherence data for SVM
    meanAccuraciesFileSVM_coh_22_9b = fullfile(myPath, sprintf('meanAccuracies_coh_9b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_coh_22_9b, 'file') == 2
        load(meanAccuraciesFileSVM_coh_22_9b, 'meanAccuraciesSVM_coh');
        if exist('meanAccuraciesSVM_coh', 'var') && ~isempty(meanAccuraciesSVM_coh)
            tempSVM_coh = meanAccuraciesSVM_coh; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM coherence file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_coh_22_9b);
            continue;
        end
    else
        fprintf('Subject %d SVM coherence file not found: %s\n', numSub, meanAccuraciesFileSVM_coh_22_9b);
        continue;
    end

    % Store coherence data for LDA, DT, and SVM
    if ~isempty(tempLDA_coh) && size(tempLDA_coh,1) >= nConditions && size(tempLDA_coh,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('LDA_COH_22_9b_feat%s_condition%d{numSub} = tempLDA_coh(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_coh) && size(tempDT_coh,1) >= nConditions && size(tempDT_coh,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('DT_COH_22_9b_feat%s_condition%d{numSub} = tempDT_coh(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_coh) && size(tempSVM_coh,1) >= nConditions && size(tempSVM_coh,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('SVM_COH_22_9b_feat%s_condition%d{numSub} = tempSVM_coh(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean coherence accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_9b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_COH_22_9b_feat%s_condition%d = mean(cell2mat(LDA_COH_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_DT_COH_22_9b_feat%s_condition%d = mean(cell2mat(DT_COH_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_SVM_COH_22_9b_feat%s_condition%d = mean(cell2mat(SVM_COH_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
    end
end

%% Display the computed means for coherence features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM COHERENCE 22 9b:\n');
for iF = 1:length(feature_number_9b)
    fprintf('\nFeature %s:\n', feature_number_9b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_COH_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_COH_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_COH_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%
%% OVERALL COH BASED ON FEATURE NUMBER, CLASSIFIERS, FREQUENCY BANDS AND CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%
% Define feature numbers for 5b and 9b
feature_number_5b = {'10', '20', '30', '40', '50'};
feature_number_9b = {'10', '20', '30', '40', '90'};

% Define electrodes 22 and 64
electrodes = {'22', '64'};

% Initialize containers to store mean values
meanSVM = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA = zeros(totalSubjects, length(feature_number_5b), nConditions);

overallMean_5b = zeros(totalSubjects, length(feature_number_5b));
overallMean_9b = zeros(totalSubjects, length(feature_number_9b));

overallMean = zeros(totalSubjects, 1);  % Final overall mean per subject
overallMean_SVM = zeros(totalSubjects, 1);
overallMean_DT = zeros(totalSubjects, 1);
overallMean_LDA = zeros(totalSubjects, 1);

% Initialize containers for storing results for electrodes 22 and 64
meanSVM_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);

meanSVM_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);

% Initialize variables to store final overall mean across all classifiers for each electrode
overallMean_22 = zeros(totalSubjects, 1);
overallMean_64 = zeros(totalSubjects, 1);

% Loop through subjects
for numSub = 1:totalSubjects
    % Initialize accumulators for all classifiers
    svmAll = [];
    dtAll = [];
    ldaAll = [];
    svmAll_22 = [];
    dtAll_22 = [];
    ldaAll_22 = [];

    svmAll_64 = [];
    dtAll_64 = [];
    ldaAll_64 = [];
    % Loop through different configurations
    configIndex = 0;
    for el = electrodes
        for band = bands
            % Check if the 9b band is being used
            if strcmp(band{1}, '9b') && isempty(feature_number_9b)
                continue;  % Skip 9b calculations if there are no features for 9b
            end
            configIndex = configIndex + 1;
            feature_number = eval(['feature_number_', band{1}]);

            for iF = 1:length(feature_number)
                % Loop through conditions
                for cond = 1:nConditions
                    % Dynamically construct variable names
                    svmVar = sprintf('SVM_COH_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);
                    dtVar = sprintf('DT_COH_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);
                    ldaVar = sprintf('LDA_COH_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);

                    % Extract values only if they are not empty
                    svmValue = eval(svmVar);
                    dtValue = eval(dtVar);
                    ldaValue = eval(ldaVar);

                    if ~isempty(svmValue{numSub})
                        svmAll = [svmAll, svmValue{numSub}];
                    end
                    if ~isempty(dtValue{numSub})
                        dtAll = [dtAll, dtValue{numSub}];
                    end
                    if ~isempty(ldaValue{numSub})
                        ldaAll = [ldaAll, ldaValue{numSub}];
                    end
                    if ~isempty(svmValue{numSub})
                        if strcmp(el{1}, '22')
                            svmAll_22 = [svmAll_22, svmValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            svmAll_64 = [svmAll_64, svmValue{numSub}];
                        end
                    end

                    if ~isempty(dtValue{numSub})
                        if strcmp(el{1}, '22')
                            dtAll_22 = [dtAll_22, dtValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            dtAll_64 = [dtAll_64, dtValue{numSub}];
                        end
                    end

                    if ~isempty(ldaValue{numSub})
                        if strcmp(el{1}, '22')
                            ldaAll_22 = [ldaAll_22, ldaValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            ldaAll_64 = [ldaAll_64, ldaValue{numSub}];
                        end
                    end
                    % Compute mean for each classifier per feature number, ignoring NaN
                    meanSVM(numSub, iF, cond) = mean(svmAll, 'omitnan');
                    meanDT(numSub, iF, cond) = mean(dtAll, 'omitnan');
                    meanLDA(numSub, iF, cond) = mean(ldaAll, 'omitnan');
                end
                % Compute overall mean per feature number for each subject
                allAccuracies = [svmAll, dtAll, ldaAll];
                if strcmp(band{1}, '5b')
                    overallMean_5b(numSub, iF) = mean(allAccuracies, 'omitnan');
                else
                    overallMean_9b(numSub, iF) = mean(allAccuracies, 'omitnan');
                end
            end

        end
    end

    % Compute overall mean across all configurations, classifiers, and feature indices
    overallMean(numSub) = mean([meanSVM(numSub, :, :), meanDT(numSub, :, :), meanLDA(numSub, :, :)], 'all', 'omitnan');
    % Compute the overall COH mean per classifier
    overallMean_SVM(numSub) = mean(svmAll, 'omitnan');
    overallMean_DT(numSub) = mean(dtAll, 'omitnan');
    overallMean_LDA(numSub) = mean(ldaAll, 'omitnan');


    % Calculate the mean values for each classifier (SVM, DT, LDA) per feature number for electrode 22 and 64
    for iF = 1:length(feature_number_5b)
        for cond = 1:nConditions
            meanSVM_22(numSub, iF, cond) = mean(svmAll_22, 'omitnan');
            meanDT_22(numSub, iF, cond) = mean(dtAll_22, 'omitnan');
            meanLDA_22(numSub, iF, cond) = mean(ldaAll_22, 'omitnan');

            meanSVM_64(numSub, iF, cond) = mean(svmAll_64, 'omitnan');
            meanDT_64(numSub, iF, cond) = mean(dtAll_64, 'omitnan');
            meanLDA_64(numSub, iF, cond) = mean(ldaAll_64, 'omitnan');
        end
    end

    % Compute the overall mean across all classifiers for each electrode (22 and 64)
    overallMean_22(numSub) = mean([meanSVM_22(numSub, :, :), meanDT_22(numSub, :, :), meanLDA_22(numSub, :, :)], 'all', 'omitnan');
    overallMean_64(numSub) = mean([meanSVM_64(numSub, :, :), meanDT_64(numSub, :, :), meanLDA_64(numSub, :, :)], 'all', 'omitnan');
end


% Display results for all subjects
fprintf('\n=== Summary of Mean Values for All Subjects (COH) ===\n');
for numSub = 1:totalSubjects
    fprintf('Subject %d:\n', numSub);
    fprintf('  Overall Mean - COH: %.4f\n\n', overallMean(numSub));
    fprintf('  Overall Mean - SVM: %.4f\n', overallMean_SVM(numSub));
    fprintf('  Overall Mean - DT: %.4f\n', overallMean_DT(numSub));
    fprintf('  Overall Mean - LDA: %.4f\n', overallMean_LDA(numSub));
    % Print overall mean per feature for 5b
    for iF = 1:length(feature_number_5b)
        fprintf('  Overall Mean - 5b (Feature %s): %.4f\n', feature_number_5b{iF}, overallMean_5b(numSub, iF));
    end

    % Print overall mean per feature for 9b
    for iF = 1:length(feature_number_9b)
        fprintf('  Overall Mean - 9b (Feature %s): %.4f\n', feature_number_9b{iF}, overallMean_9b(numSub, iF));
    end
    % Print overall mean for electrode 22
    fprintf('  Overall Mean - Electrode 22: %.4f\n', overallMean_22(numSub));

    % Print overall mean for electrode 64
    fprintf('  Overall Mean - Electrode 64: %.4f\n', overallMean_64(numSub));

    fprintf('\n');  % Print a newline after each subject's results
end
% Saving to a specific path
save(fullfile(myPath, 'overal_COH.mat'), 'overallMean', 'meanSVM', 'meanDT', 'meanLDA');
save(fullfile(myPath, 'overall_COH_classifiers.mat'), 'overallMean_SVM', 'overallMean_DT', 'overallMean_LDA');
save(fullfile(myPath, 'overall_COH_features.mat'), ...
    'overallMean_5b', 'overallMean_9b');
save(fullfile(myPath, 'overall_COH_channels.mat'), 'overallMean_22', 'overallMean_64');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORRELATION 64 channel 5 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for Correlation 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA Correlation measures
            eval(sprintf('SVM_CORR_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_CORR_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_CORR_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for Correlation 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_CORR_64_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_CORR_64_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_CORR_64_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for Correlation
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load Correlation data for LDA
    meanAccuraciesFileLDA_corr_64_5b = fullfile(myPath, sprintf('meanAccuraciesFC_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_corr_64_5b, 'file') == 2
        load(meanAccuraciesFileLDA_corr_64_5b, 'meanAccuracies');
        if exist('meanAccuracies', 'var') && ~isempty(meanAccuracies)
            tempLDA_corr = meanAccuracies; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_corr_64_5b);
            continue;
        end
    else
        fprintf('Subject %d LDA Correlation file not found: %s\n', numSub, meanAccuraciesFileLDA_corr_64_5b);
        continue;
    end

    % Load Correlation data for DT
    meanAccuraciesFileDT_corr_64_5b = fullfile(myPath, sprintf('meanAccuracies_corr_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_corr_64_5b, 'file') == 2
        load(meanAccuraciesFileDT_corr_64_5b, 'meanAccuracies_DT');
        if exist('meanAccuracies_DT', 'var') && ~isempty(meanAccuracies_DT)
            tempDT_corr = meanAccuracies_DT; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_corr_64_5b);
            continue;
        end
    else
        fprintf('Subject %d DT Correlation file not found: %s\n', numSub, meanAccuraciesFileDT_corr_64_5b);
        continue;
    end

    % Load Correlation data for SVM
    meanAccuraciesFileSVM_corr_64_5b = fullfile(myPath, sprintf('meanAccuracies_corr_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_corr_64_5b, 'file') == 2
        load(meanAccuraciesFileSVM_corr_64_5b, 'meanAccuracies_SVM');
        if exist('meanAccuracies_SVM', 'var') && ~isempty(meanAccuracies_SVM)
            tempSVM_corr = meanAccuracies_SVM; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_corr_64_5b);
            continue;
        end
    else
        fprintf('Subject %d SVM Correlation file not found: %s\n', numSub, meanAccuraciesFileSVM_corr_64_5b);
        continue;
    end

    % Store Correlation data for LDA, DT, and SVM
    if ~isempty(tempLDA_corr) && size(tempLDA_corr,1) >= nConditions && size(tempLDA_corr,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('LDA_CORR_64_5b_feat%s_condition%d{numSub} = tempLDA_corr(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_corr) && size(tempDT_corr,1) >= nConditions && size(tempDT_corr,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('DT_CORR_64_5b_feat%s_condition%d{numSub} = tempDT_corr(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_corr) && size(tempSVM_corr,1) >= nConditions && size(tempSVM_corr,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('SVM_CORR_64_5b_feat%s_condition%d{numSub} = tempSVM_corr(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean Correlation accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_5b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_CORR_64_5b_feat%s_condition%d = mean(cell2mat(LDA_CORR_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_DT_CORR_64_5b_feat%s_condition%d = mean(cell2mat(DT_CORR_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_SVM_CORR_64_5b_feat%s_condition%d = mean(cell2mat(SVM_CORR_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
    end
end

%% Display the computed means for Correlation features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM Correlation 64 5b:\n');
for iF = 1:length(feature_number_5b)
    fprintf('\nFeature %s:\n', feature_number_5b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_CORR_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_CORR_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_CORR_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlation 22 channel 5 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for Correlation 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA Correlation measures
            eval(sprintf('SVM_CORR_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_CORR_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_CORR_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for Correlation 22-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_CORR_22_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_CORR_22_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_CORR_22_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for Correlation
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load Correlation data for LDA
    meanAccuraciesFileLDA_corr_22_5b = fullfile(myPath, sprintf('meanAccuracies_corr_5b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_corr_22_5b, 'file') == 2
        load(meanAccuraciesFileLDA_corr_22_5b, 'meanAccuraciesLDA_corr');
        if exist('meanAccuraciesLDA_corr', 'var') && ~isempty(meanAccuraciesLDA_corr)
            tempLDA_corr = meanAccuraciesLDA_corr; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_corr_22_5b);
            continue;
        end
    else
        fprintf('Subject %d LDA Correlation file not found: %s\n', numSub, meanAccuraciesFileLDA_corr_22_5b);
        continue;
    end

    % Load Correlation data for DT
    meanAccuraciesFileDT_corr_22_5b = fullfile(myPath, sprintf('meanAccuracies_corr_5b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_corr_22_5b, 'file') == 2
        load(meanAccuraciesFileDT_corr_22_5b, 'meanAccuraciesDT_corr');
        if exist('meanAccuraciesDT_corr', 'var') && ~isempty(meanAccuraciesDT_corr)
            tempDT_corr = meanAccuraciesDT_corr; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_corr_22_5b);
            continue;
        end
    else
        fprintf('Subject %d DT Correlation file not found: %s\n', numSub, meanAccuraciesFileDT_corr_22_5b);
        continue;
    end

    % Load Correlation data for SVM
    meanAccuraciesFileSVM_corr_22_5b = fullfile(myPath, sprintf('meanAccuracies_corr_5b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_corr_22_5b, 'file') == 2
        load(meanAccuraciesFileSVM_corr_22_5b, 'meanAccuraciesSVM_corr');
        if exist('meanAccuraciesSVM_corr', 'var') && ~isempty(meanAccuraciesSVM_corr)
            tempSVM_corr = meanAccuraciesSVM_corr; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_corr_22_5b);
            continue;
        end
    else
        fprintf('Subject %d SVM Correlation file not found: %s\n', numSub, meanAccuraciesFileSVM_corr_22_5b);
        continue;
    end

    % Store Correlation data for LDA, DT, and SVM
    if ~isempty(tempLDA_corr) && size(tempLDA_corr,1) >= nConditions && size(tempLDA_corr,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('LDA_CORR_22_5b_feat%s_condition%d{numSub} = tempLDA_corr(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_corr) && size(tempDT_corr,1) >= nConditions && size(tempDT_corr,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('DT_CORR_22_5b_feat%s_condition%d{numSub} = tempDT_corr(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_corr) && size(tempSVM_corr,1) >= nConditions && size(tempSVM_corr,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('SVM_CORR_22_5b_feat%s_condition%d{numSub} = tempSVM_corr(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean Correlation accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_5b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_CORR_22_5b_feat%s_condition%d = mean(cell2mat(LDA_CORR_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_DT_CORR_22_5b_feat%s_condition%d = mean(cell2mat(DT_CORR_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_SVM_CORR_22_5b_feat%s_condition%d = mean(cell2mat(SVM_CORR_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
    end
end

%% Display the computed means for Correlation features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM Correlation 22 5b:\n');
for iF = 1:length(feature_number_5b)
    fprintf('\nFeature %s:\n', feature_number_5b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_CORR_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_CORR_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_CORR_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORRELATION 64 channel 9 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for Correlation 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA Correlation measures
            eval(sprintf('SVM_CORR_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_CORR_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_CORR_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for Correlation 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_CORR_64_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_CORR_64_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_CORR_64_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for Correlation
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load Correlation data for LDA
    meanAccuraciesFileLDA_corr_64_9b = fullfile(myPath, sprintf('meanAccuraciesFC_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_corr_64_9b, 'file') == 2
        load(meanAccuraciesFileLDA_corr_64_9b, 'meanAccuracies');
        if exist('meanAccuracies', 'var') && ~isempty(meanAccuracies)
            tempLDA_corr = meanAccuracies; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_corr_64_9b);
            continue;
        end
    else
        fprintf('Subject %d LDA Correlation file not found: %s\n', numSub, meanAccuraciesFileLDA_corr_64_9b);
        continue;
    end

    % Load Correlation data for DT
    meanAccuraciesFileDT_corr_64_9b = fullfile(myPath, sprintf('meanAccuracies_corr_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_corr_64_9b, 'file') == 2
        load(meanAccuraciesFileDT_corr_64_9b, 'meanAccuracies_DT');
        if exist('meanAccuracies_DT', 'var') && ~isempty(meanAccuracies_DT)
            tempDT_corr = meanAccuracies_DT; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_corr_64_9b);
            continue;
        end
    else
        fprintf('Subject %d DT Correlation file not found: %s\n', numSub, meanAccuraciesFileDT_corr_64_9b);
        continue;
    end

    % Load Correlation data for SVM
    meanAccuraciesFileSVM_corr_64_9b = fullfile(myPath, sprintf('meanAccuracies_corr_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_corr_64_9b, 'file') == 2
        load(meanAccuraciesFileSVM_corr_64_9b, 'meanAccuracies_SVM');
        if exist('meanAccuracies_SVM', 'var') && ~isempty(meanAccuracies_SVM)
            tempSVM_corr = meanAccuracies_SVM; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_corr_64_9b);
            continue;
        end
    else
        fprintf('Subject %d SVM Correlation file not found: %s\n', numSub, meanAccuraciesFileSVM_corr_64_9b);
        continue;
    end

    % Store Correlation data for LDA, DT, and SVM
    if ~isempty(tempLDA_corr) && size(tempLDA_corr,1) >= nConditions && size(tempLDA_corr,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('LDA_CORR_64_9b_feat%s_condition%d{numSub} = tempLDA_corr(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_corr) && size(tempDT_corr,1) >= nConditions && size(tempDT_corr,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('DT_CORR_64_9b_feat%s_condition%d{numSub} = tempDT_corr(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_corr) && size(tempSVM_corr,1) >= nConditions && size(tempSVM_corr,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('SVM_CORR_64_9b_feat%s_condition%d{numSub} = tempSVM_corr(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean Correlation accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_9b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_CORR_64_9b_feat%s_condition%d = mean(cell2mat(LDA_CORR_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_DT_CORR_64_9b_feat%s_condition%d = mean(cell2mat(DT_CORR_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_SVM_CORR_64_9b_feat%s_condition%d = mean(cell2mat(SVM_CORR_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
    end
end

%% Display the computed means for Correlation features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM Correlation 64 9b:\n');
for iF = 1:length(feature_number_9b)
    fprintf('\nFeature %s:\n', feature_number_9b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_CORR_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_CORR_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_CORR_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlation 22 channel 9 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for Correlation 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA Correlation measures
            eval(sprintf('SVM_CORR_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_CORR_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_CORR_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for Correlation 22-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_CORR_22_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_CORR_22_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_CORR_22_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for Correlation
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load Correlation data for LDA
    meanAccuraciesFileLDA_corr_22_9b = fullfile(myPath, sprintf('meanAccuracies_corr_9b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_corr_22_9b, 'file') == 2
        load(meanAccuraciesFileLDA_corr_22_9b, 'meanAccuraciesLDA_corr');
        if exist('meanAccuraciesLDA_corr', 'var') && ~isempty(meanAccuraciesLDA_corr)
            tempLDA_corr = meanAccuraciesLDA_corr; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_corr_22_9b);
            continue;
        end
    else
        fprintf('Subject %d LDA Correlation file not found: %s\n', numSub, meanAccuraciesFileLDA_corr_22_9b);
        continue;
    end

    % Load Correlation data for DT
    meanAccuraciesFileDT_corr_22_9b = fullfile(myPath, sprintf('meanAccuracies_corr_9b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_corr_22_9b, 'file') == 2
        load(meanAccuraciesFileDT_corr_22_9b, 'meanAccuraciesDT_corr');
        if exist('meanAccuraciesDT_corr', 'var') && ~isempty(meanAccuraciesDT_corr)
            tempDT_corr = meanAccuraciesDT_corr; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_corr_22_9b);
            continue;
        end
    else
        fprintf('Subject %d DT Correlation file not found: %s\n', numSub, meanAccuraciesFileDT_corr_22_9b);
        continue;
    end

    % Load Correlation data for SVM
    meanAccuraciesFileSVM_corr_22_9b = fullfile(myPath, sprintf('meanAccuracies_corr_9b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_corr_22_9b, 'file') == 2
        load(meanAccuraciesFileSVM_corr_22_9b, 'meanAccuraciesSVM_corr');
        if exist('meanAccuraciesSVM_corr', 'var') && ~isempty(meanAccuraciesSVM_corr)
            tempSVM_corr = meanAccuraciesSVM_corr; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM Correlation file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_corr_22_9b);
            continue;
        end
    else
        fprintf('Subject %d SVM Correlation file not found: %s\n', numSub, meanAccuraciesFileSVM_corr_22_9b);
        continue;
    end

    % Store Correlation data for LDA, DT, and SVM
    if ~isempty(tempLDA_corr) && size(tempLDA_corr,1) >= nConditions && size(tempLDA_corr,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('LDA_CORR_22_9b_feat%s_condition%d{numSub} = tempLDA_corr(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_corr) && size(tempDT_corr,1) >= nConditions && size(tempDT_corr,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('DT_CORR_22_9b_feat%s_condition%d{numSub} = tempDT_corr(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_corr) && size(tempSVM_corr,1) >= nConditions && size(tempSVM_corr,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('SVM_CORR_22_9b_feat%s_condition%d{numSub} = tempSVM_corr(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean Correlation accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_9b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_CORR_22_9b_feat%s_condition%d = mean(cell2mat(LDA_CORR_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_DT_CORR_22_9b_feat%s_condition%d = mean(cell2mat(DT_CORR_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_SVM_CORR_22_9b_feat%s_condition%d = mean(cell2mat(SVM_CORR_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
    end
end

%% Display the computed means for Correlation features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM Correlation 22 9b:\n');
for iF = 1:length(feature_number_9b)
    fprintf('\nFeature %s:\n', feature_number_9b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_CORR_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_CORR_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_CORR_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end





%%%%%%%%%%%%%%%%%%%%%%%%%
%% OVERALL CORR BASED ON FEATURE NUMBER, CLASSIFIERS, FREQUENCY BANDS AND CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%
% Define feature numbers for 5b and 9b
feature_number_5b = {'10', '20', '30', '40', '50'};
feature_number_9b = {'10', '20', '30', '40', '90'};

% Define electrodes 22 and 64
electrodes = {'22', '64'};

% Initialize containers to store mean values
meanSVM = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA = zeros(totalSubjects, length(feature_number_5b), nConditions);

overallMean_5b = zeros(totalSubjects, length(feature_number_5b));
overallMean_9b = zeros(totalSubjects, length(feature_number_9b));

overallMean = zeros(totalSubjects, 1);  % Final overall mean per subject
overallMean_SVM = zeros(totalSubjects, 1);
overallMean_DT = zeros(totalSubjects, 1);
overallMean_LDA = zeros(totalSubjects, 1);

% Initialize containers for storing results for electrodes 22 and 64
meanSVM_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);

meanSVM_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);

% Initialize variables to store final overall mean across all classifiers for each electrode
overallMean_22 = zeros(totalSubjects, 1);
overallMean_64 = zeros(totalSubjects, 1);

% Loop through subjects
for numSub = 1:totalSubjects
    % Initialize accumulators for all classifiers
    svmAll = [];
    dtAll = [];
    ldaAll = [];
    svmAll_22 = [];
    dtAll_22 = [];
    ldaAll_22 = [];

    svmAll_64 = [];
    dtAll_64 = [];
    ldaAll_64 = [];
    % Loop through different configurations
    configIndex = 0;
    for el = electrodes
        for band = bands
            % Check if the 9b band is being used
            if strcmp(band{1}, '9b') && isempty(feature_number_9b)
                continue;  % Skip 9b calculations if there are no features for 9b
            end
            configIndex = configIndex + 1;
            feature_number = eval(['feature_number_', band{1}]);

            for iF = 1:length(feature_number)
                % Loop through conditions
                for cond = 1:nConditions
                    % Dynamically construct variable names
                    svmVar = sprintf('SVM_CORR_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);
                    dtVar = sprintf('DT_CORR_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);
                    ldaVar = sprintf('LDA_CORR_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);

                    % Extract values only if they are not empty
                    svmValue = eval(svmVar);
                    dtValue = eval(dtVar);
                    ldaValue = eval(ldaVar);

                    if ~isempty(svmValue{numSub})
                        svmAll = [svmAll, svmValue{numSub}];
                    end
                    if ~isempty(dtValue{numSub})
                        dtAll = [dtAll, dtValue{numSub}];
                    end
                    if ~isempty(ldaValue{numSub})
                        ldaAll = [ldaAll, ldaValue{numSub}];
                    end
                    if ~isempty(svmValue{numSub})
                        if strcmp(el{1}, '22')
                            svmAll_22 = [svmAll_22, svmValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            svmAll_64 = [svmAll_64, svmValue{numSub}];
                        end
                    end

                    if ~isempty(dtValue{numSub})
                        if strcmp(el{1}, '22')
                            dtAll_22 = [dtAll_22, dtValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            dtAll_64 = [dtAll_64, dtValue{numSub}];
                        end
                    end

                    if ~isempty(ldaValue{numSub})
                        if strcmp(el{1}, '22')
                            ldaAll_22 = [ldaAll_22, ldaValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            ldaAll_64 = [ldaAll_64, ldaValue{numSub}];
                        end
                    end
                    % Compute mean for each classifier per feature number, ignoring NaN
                    meanSVM(numSub, iF, cond) = mean(svmAll, 'omitnan');
                    meanDT(numSub, iF, cond) = mean(dtAll, 'omitnan');
                    meanLDA(numSub, iF, cond) = mean(ldaAll, 'omitnan');
                end
                % Compute overall mean per feature number for each subject
                allAccuracies = [svmAll, dtAll, ldaAll];
                if strcmp(band{1}, '5b')
                    overallMean_5b(numSub, iF) = mean(allAccuracies, 'omitnan');
                else
                    overallMean_9b(numSub, iF) = mean(allAccuracies, 'omitnan');
                end
            end

        end
    end

    % Compute overall mean across all configurations, classifiers, and feature indices
    overallMean(numSub) = mean([meanSVM(numSub, :, :), meanDT(numSub, :, :), meanLDA(numSub, :, :)], 'all', 'omitnan');
    % Compute the overall CORR mean per classifier
    overallMean_SVM(numSub) = mean(svmAll, 'omitnan');
    overallMean_DT(numSub) = mean(dtAll, 'omitnan');
    overallMean_LDA(numSub) = mean(ldaAll, 'omitnan');


    % Calculate the mean values for each classifier (SVM, DT, LDA) per feature number for electrode 22 and 64
    for iF = 1:length(feature_number_5b)
        for cond = 1:nConditions
            meanSVM_22(numSub, iF, cond) = mean(svmAll_22, 'omitnan');
            meanDT_22(numSub, iF, cond) = mean(dtAll_22, 'omitnan');
            meanLDA_22(numSub, iF, cond) = mean(ldaAll_22, 'omitnan');

            meanSVM_64(numSub, iF, cond) = mean(svmAll_64, 'omitnan');
            meanDT_64(numSub, iF, cond) = mean(dtAll_64, 'omitnan');
            meanLDA_64(numSub, iF, cond) = mean(ldaAll_64, 'omitnan');
        end
    end

    % Compute the overall mean across all classifiers for each electrode (22 and 64)
    overallMean_22(numSub) = mean([meanSVM_22(numSub, :, :), meanDT_22(numSub, :, :), meanLDA_22(numSub, :, :)], 'all', 'omitnan');
    overallMean_64(numSub) = mean([meanSVM_64(numSub, :, :), meanDT_64(numSub, :, :), meanLDA_64(numSub, :, :)], 'all', 'omitnan');
end


% Display results for all subjects
fprintf('\n=== Summary of Mean Values for All Subjects (CORR) ===\n');
for numSub = 1:totalSubjects
    fprintf('Subject %d:\n', numSub);
    fprintf('  Overall Mean - CORR: %.4f\n\n', overallMean(numSub));
    fprintf('  Overall Mean - SVM: %.4f\n', overallMean_SVM(numSub));
    fprintf('  Overall Mean - DT: %.4f\n', overallMean_DT(numSub));
    fprintf('  Overall Mean - LDA: %.4f\n', overallMean_LDA(numSub));
    % Print overall mean per feature for 5b
    for iF = 1:length(feature_number_5b)
        fprintf('  Overall Mean - 5b (Feature %s): %.4f\n', feature_number_5b{iF}, overallMean_5b(numSub, iF));
    end

    % Print overall mean per feature for 9b
    for iF = 1:length(feature_number_9b)
        fprintf('  Overall Mean - 9b (Feature %s): %.4f\n', feature_number_9b{iF}, overallMean_9b(numSub, iF));
    end
    % Print overall mean for electrode 22
    fprintf('  Overall Mean - Electrode 22: %.4f\n', overallMean_22(numSub));

    % Print overall mean for electrode 64
    fprintf('  Overall Mean - Electrode 64: %.4f\n', overallMean_64(numSub));

    fprintf('\n');  % Print a newline after each subject's results
end
% Saving to a specific path
save(fullfile(myPath, 'overal_CORR.mat'), 'overallMean', 'meanSVM', 'meanDT', 'meanLDA');
save(fullfile(myPath, 'overall_CORR_classifiers.mat'), 'overallMean_SVM', 'overallMean_DT', 'overallMean_LDA');
save(fullfile(myPath, 'overall_CORR_features.mat'), ...
    'overallMean_5b', 'overallMean_9b');
save(fullfile(myPath, 'overall_CORR_channels.mat'), 'overallMean_22', 'overallMean_64');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHASE LOCKING 64 channel 5 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for PLV 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA PLV measures
            eval(sprintf('SVM_PLV_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_PLV_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_PLV_64_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for PLV 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_PLV_64_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_PLV_64_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_PLV_64_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for PLV
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load PLV data for LDA
    meanAccuraciesFileLDA_plv_64_5b = fullfile(myPath, sprintf('meanAccuracies_plv_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_plv_64_5b, 'file') == 2
        load(meanAccuraciesFileLDA_plv_64_5b, 'meanAccuracies_plv');
        if exist('meanAccuracies_plv', 'var') && ~isempty(meanAccuracies_plv)
            tempLDA_plv = meanAccuracies_plv; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_plv_64_5b);
            continue;
        end
    else
        fprintf('Subject %d LDA PLV file not found: %s\n', numSub, meanAccuraciesFileLDA_plv_64_5b);
        continue;
    end

    % Load PLV data for DT
    meanAccuraciesFileDT_plv_64_5b = fullfile(myPath, sprintf('meanAccuracies_plv_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_plv_64_5b, 'file') == 2
        load(meanAccuraciesFileDT_plv_64_5b, 'meanAccuracies_plv_DT');
        if exist('meanAccuracies_plv_DT', 'var') && ~isempty(meanAccuracies_plv_DT)
            tempDT_plv = meanAccuracies_plv_DT; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_plv_64_5b);
            continue;
        end
    else
        fprintf('Subject %d DT PLV file not found: %s\n', numSub, meanAccuraciesFileDT_plv_64_5b);
        continue;
    end

    % Load PLV data for SVM
    meanAccuraciesFileSVM_plv_64_5b = fullfile(myPath, sprintf('meanAccuracies_plv_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_plv_64_5b, 'file') == 2
        load(meanAccuraciesFileSVM_plv_64_5b, 'meanAccuracies_plv_SVM');
        if exist('meanAccuracies_plv_SVM', 'var') && ~isempty(meanAccuracies_plv_SVM)
            tempSVM_plv = meanAccuracies_plv_SVM; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_plv_64_5b);
            continue;
        end
    else
        fprintf('Subject %d SVM PLV file not found: %s\n', numSub, meanAccuraciesFileSVM_plv_64_5b);
        continue;
    end

    % Store PLV data for LDA, DT, and SVM
    if ~isempty(tempLDA_plv) && size(tempLDA_plv,1) >= nConditions && size(tempLDA_plv,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('LDA_PLV_64_5b_feat%s_condition%d{numSub} = tempLDA_plv(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_plv) && size(tempDT_plv,1) >= nConditions && size(tempDT_plv,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('DT_PLV_64_5b_feat%s_condition%d{numSub} = tempDT_plv(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_plv) && size(tempSVM_plv,1) >= nConditions && size(tempSVM_plv,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('SVM_PLV_64_5b_feat%s_condition%d{numSub} = tempSVM_plv(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean PLV accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_5b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_PLV_64_5b_feat%s_condition%d = mean(cell2mat(LDA_PLV_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_DT_PLV_64_5b_feat%s_condition%d = mean(cell2mat(DT_PLV_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_SVM_PLV_64_5b_feat%s_condition%d = mean(cell2mat(SVM_PLV_64_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
    end
end

%% Display the computed means for PLV features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM PLV 64 5b:\n');
for iF = 1:length(feature_number_5b)
    fprintf('\nFeature %s:\n', feature_number_5b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_PLV_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_PLV_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_PLV_64_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHASE LOCKING 22 channel 5 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for PLV 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA PLV measures
            eval(sprintf('SVM_PLV_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_PLV_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_PLV_22_5b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for PLV 22-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_PLV_22_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_PLV_22_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_PLV_22_5b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for PLV
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load PLV data for LDA
    meanAccuraciesFileLDA_plv_22_5b = fullfile(myPath, sprintf('meanAccuracies_pLV_5b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_plv_22_5b, 'file') == 2
        load(meanAccuraciesFileLDA_plv_22_5b, 'meanAccuraciesLDA_plv');
        if exist('meanAccuraciesLDA_plv', 'var') && ~isempty(meanAccuraciesLDA_plv)
            tempLDA_plv = meanAccuraciesLDA_plv; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_plv_22_5b);
            continue;
        end
    else
        fprintf('Subject %d LDA PLV file not found: %s\n', numSub, meanAccuraciesFileLDA_plv_22_5b);
        continue;
    end

    % Load PLV data for DT
    meanAccuraciesFileDT_plv_22_5b = fullfile(myPath, sprintf('meanAccuracies_pLV_5b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_plv_22_5b, 'file') == 2
        load(meanAccuraciesFileDT_plv_22_5b, 'meanAccuraciesDT_plv');
        if exist('meanAccuraciesDT_plv', 'var') && ~isempty(meanAccuraciesDT_plv)
            tempDT_plv = meanAccuraciesDT_plv; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_plv_22_5b);
            continue;
        end
    else
        fprintf('Subject %d DT PLV file not found: %s\n', numSub, meanAccuraciesFileDT_plv_22_5b);
        continue;
    end

    % Load PLV data for SVM
    meanAccuraciesFileSVM_plv_22_5b = fullfile(myPath, sprintf('meanAccuracies_pLV_5b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_plv_22_5b, 'file') == 2
        load(meanAccuraciesFileSVM_plv_22_5b, 'meanAccuraciesSVM_plv');
        if exist('meanAccuraciesSVM_plv', 'var') && ~isempty(meanAccuraciesSVM_plv)
            tempSVM_plv = meanAccuraciesSVM_plv; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_plv_22_5b);
            continue;
        end
    else
        fprintf('Subject %d SVM PLV file not found: %s\n', numSub, meanAccuraciesFileSVM_plv_22_5b);
        continue;
    end

    % Store PLV data for LDA, DT, and SVM
    if ~isempty(tempLDA_plv) && size(tempLDA_plv,1) >= nConditions && size(tempLDA_plv,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('LDA_PLV_22_5b_feat%s_condition%d{numSub} = tempLDA_plv(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_plv) && size(tempDT_plv,1) >= nConditions && size(tempDT_plv,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('DT_PLV_22_5b_feat%s_condition%d{numSub} = tempDT_plv(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_plv) && size(tempSVM_plv,1) >= nConditions && size(tempSVM_plv,2) >= length(feature_number_5b)
        for iF = 1:length(feature_number_5b)
            for cond = 1:nConditions
                eval(sprintf('SVM_PLV_22_5b_feat%s_condition%d{numSub} = tempSVM_plv(%d,%d);', ...
                    feature_number_5b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean PLV accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_5b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_PLV_22_5b_feat%s_condition%d = mean(cell2mat(LDA_PLV_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_DT_PLV_22_5b_feat%s_condition%d = mean(cell2mat(DT_PLV_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
        eval(sprintf('mean_SVM_PLV_22_5b_feat%s_condition%d = mean(cell2mat(SVM_PLV_22_5b_feat%s_condition%d));', ...
            feature_number_5b{iF}, cond, feature_number_5b{iF}, cond));
    end
end

%% Display the computed means for PLV features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM PLV 22 5b:\n');
for iF = 1:length(feature_number_5b)
    fprintf('\nFeature %s:\n', feature_number_5b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_PLV_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_PLV_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_PLV_22_5b_feat%s_condition%d', feature_number_5b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHASE LOCKING 64 channel 9 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for PLV 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA PLV measures
            eval(sprintf('SVM_PLV_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_PLV_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_PLV_64_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for PLV 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_PLV_64_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_PLV_64_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_PLV_64_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for PLV
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load PLV data for LDA
    meanAccuraciesFileLDA_plv_64_9b = fullfile(myPath, sprintf('meanAccuracies_plv_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_plv_64_9b, 'file') == 2
        load(meanAccuraciesFileLDA_plv_64_9b, 'meanAccuracies_plv');
        if exist('meanAccuracies_plv', 'var') && ~isempty(meanAccuracies_plv)
            tempLDA_plv = meanAccuracies_plv; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_plv_64_9b);
            continue;
        end
    else
        fprintf('Subject %d LDA PLV file not found: %s\n', numSub, meanAccuraciesFileLDA_plv_64_9b);
        continue;
    end

    % Load PLV data for DT
    meanAccuraciesFileDT_plv_64_9b = fullfile(myPath, sprintf('meanAccuracies_plv_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_plv_64_9b, 'file') == 2
        load(meanAccuraciesFileDT_plv_64_9b, 'meanAccuracies_plv_DT');
        if exist('meanAccuracies_plv_DT', 'var') && ~isempty(meanAccuracies_plv_DT)
            tempDT_plv = meanAccuracies_plv_DT; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_plv_64_9b);
            continue;
        end
    else
        fprintf('Subject %d DT PLV file not found: %s\n', numSub, meanAccuraciesFileDT_plv_64_9b);
        continue;
    end

    % Load PLV data for SVM
    meanAccuraciesFileSVM_plv_64_9b = fullfile(myPath, sprintf('meanAccuracies_plv_svmdt_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_plv_64_9b, 'file') == 2
        load(meanAccuraciesFileSVM_plv_64_9b, 'meanAccuracies_plv_SVM');
        if exist('meanAccuracies_plv_SVM', 'var') && ~isempty(meanAccuracies_plv_SVM)
            tempSVM_plv = meanAccuracies_plv_SVM; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_plv_64_9b);
            continue;
        end
    else
        fprintf('Subject %d SVM PLV file not found: %s\n', numSub, meanAccuraciesFileSVM_plv_64_9b);
        continue;
    end

    % Store PLV data for LDA, DT, and SVM
    if ~isempty(tempLDA_plv) && size(tempLDA_plv,1) >= nConditions && size(tempLDA_plv,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('LDA_PLV_64_9b_feat%s_condition%d{numSub} = tempLDA_plv(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_plv) && size(tempDT_plv,1) >= nConditions && size(tempDT_plv,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('DT_PLV_64_9b_feat%s_condition%d{numSub} = tempDT_plv(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_plv) && size(tempSVM_plv,1) >= nConditions && size(tempSVM_plv,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('SVM_PLV_64_9b_feat%s_condition%d{numSub} = tempSVM_plv(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean PLV accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_9b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_PLV_64_9b_feat%s_condition%d = mean(cell2mat(LDA_PLV_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_DT_PLV_64_9b_feat%s_condition%d = mean(cell2mat(DT_PLV_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_SVM_PLV_64_9b_feat%s_condition%d = mean(cell2mat(SVM_PLV_64_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
    end
end

%% Display the computed means for PLV features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM PLV 64 9b:\n');
for iF = 1:length(feature_number_9b)
    fprintf('\nFeature %s:\n', feature_number_9b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_PLV_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_PLV_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_PLV_64_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PHASE LOCKING 22 channel 9 band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create empty containers for PLV 64-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            % Create containers for SVM, DT, and LDA PLV measures
            eval(sprintf('SVM_PLV_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('DT_PLV_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
            eval(sprintf('LDA_PLV_22_9b_feat%s_condition%d = cell(totalSubjects,1);', feat{1}, cond));
        end
    end
end
%% Create mean accuracy containers for PLV 22-channel 5-band
for band = bands
    features = getFeatures(band{1});
    for feat = features
        for cond = 1:nConditions
            eval(sprintf('mean_SVM_PLV_22_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_DT_PLV_22_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
            eval(sprintf('mean_LDA_PLV_22_9b_feat%s_condition%d = [];', feat{1}, cond)); %#ok<NASGU>
        end
    end
end
%% Initialize subjectNumbers container for PLV
subjectNumbers = [];

% Loop through all subjects
for numSub = 1:totalSubjects
    if any(numSub == excludeSubjects), continue; end
    subjectFolder = sprintf('Sub-%03d', numSub);
    myPath = fullfile(myPath, subjectFolder);

    % Load PLV data for LDA
    meanAccuraciesFileLDA_plv_22_9b = fullfile(myPath, sprintf('meanAccuracies_pLV_9b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileLDA_plv_22_9b, 'file') == 2
        load(meanAccuraciesFileLDA_plv_22_9b, 'meanAccuraciesLDA_plv');
        if exist('meanAccuraciesLDA_plv', 'var') && ~isempty(meanAccuraciesLDA_plv)
            tempLDA_plv = meanAccuraciesLDA_plv; % Expected: [6 x 5]
        else
            fprintf('Subject %d LDA PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileLDA_plv_22_9b);
            continue;
        end
    else
        fprintf('Subject %d LDA PLV file not found: %s\n', numSub, meanAccuraciesFileLDA_plv_22_9b);
        continue;
    end

    % Load PLV data for DT
    meanAccuraciesFileDT_plv_22_9b = fullfile(myPath, sprintf('meanAccuracies_pLV_9b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileDT_plv_22_9b, 'file') == 2
        load(meanAccuraciesFileDT_plv_22_9b, 'meanAccuraciesDT_plv');
        if exist('meanAccuraciesDT_plv', 'var') && ~isempty(meanAccuraciesDT_plv)
            tempDT_plv = meanAccuraciesDT_plv; % Expected: [6 x 5]
        else
            fprintf('Subject %d DT PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileDT_plv_22_9b);
            continue;
        end
    else
        fprintf('Subject %d DT PLV file not found: %s\n', numSub, meanAccuraciesFileDT_plv_22_9b);
        continue;
    end

    % Load PLV data for SVM
    meanAccuraciesFileSVM_plv_22_9b = fullfile(myPath, sprintf('meanAccuracies_pLV_9b_22ch_Sub-%03d.mat', numSub));
    if exist(meanAccuraciesFileSVM_plv_22_9b, 'file') == 2
        load(meanAccuraciesFileSVM_plv_22_9b, 'meanAccuraciesSVM_plv');
        if exist('meanAccuraciesSVM_plv', 'var') && ~isempty(meanAccuraciesSVM_plv)
            tempSVM_plv = meanAccuraciesSVM_plv; % Expected: [6 x 5]
        else
            fprintf('Subject %d SVM PLV file is empty or missing: %s\n', numSub, meanAccuraciesFileSVM_plv_22_9b);
            continue;
        end
    else
        fprintf('Subject %d SVM PLV file not found: %s\n', numSub, meanAccuraciesFileSVM_plv_22_9b);
        continue;
    end

    % Store PLV data for LDA, DT, and SVM
    if ~isempty(tempLDA_plv) && size(tempLDA_plv,1) >= nConditions && size(tempLDA_plv,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('LDA_PLV_22_9b_feat%s_condition%d{numSub} = tempLDA_plv(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempDT_plv) && size(tempDT_plv,1) >= nConditions && size(tempDT_plv,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('DT_PLV_22_9b_feat%s_condition%d{numSub} = tempDT_plv(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end

    if ~isempty(tempSVM_plv) && size(tempSVM_plv,1) >= nConditions && size(tempSVM_plv,2) >= length(feature_number_9b)
        for iF = 1:length(feature_number_9b)
            for cond = 1:nConditions
                eval(sprintf('SVM_PLV_22_9b_feat%s_condition%d{numSub} = tempSVM_plv(%d,%d);', ...
                    feature_number_9b{iF}, cond, cond, iF));
            end
        end
    end
end

%% Compute mean PLV accuracy across subjects for each feature & condition
for iF = 1:length(feature_number_9b)
    for cond = 1:nConditions
        % Compute the mean using cell2mat conversion (ignore empty cells)
        eval(sprintf('mean_LDA_PLV_22_9b_feat%s_condition%d = mean(cell2mat(LDA_PLV_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_DT_PLV_22_9b_feat%s_condition%d = mean(cell2mat(DT_PLV_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
        eval(sprintf('mean_SVM_PLV_22_9b_feat%s_condition%d = mean(cell2mat(SVM_PLV_22_9b_feat%s_condition%d));', ...
            feature_number_9b{iF}, cond, feature_number_9b{iF}, cond));
    end
end

%% Display the computed means for PLV features
fprintf('Computed Mean Accuracies for LDA, DT, and SVM PLV 22 9b:\n');
for iF = 1:length(feature_number_9b)
    fprintf('\nFeature %s:\n', feature_number_9b{iF});
    for cond = 1:nConditions
        % Retrieve the computed mean values for LDA, DT, and SVM
        mean_LDA = eval(sprintf('mean_LDA_PLV_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_DT  = eval(sprintf('mean_DT_PLV_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        mean_SVM = eval(sprintf('mean_SVM_PLV_22_9b_feat%s_condition%d', feature_number_9b{iF}, cond));
        fprintf('  Condition %d: LDA %.2f, DT %.2f, SVM %.2f\n', cond, mean_LDA, mean_DT, mean_SVM);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% OVERALL PLV BASED ON FEATURE NUMBER, CLASSIFIERS, FREQUENCY BANDS AND CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%
% Define feature numbers for 5b and 9b
feature_number_5b = {'10', '20', '30', '40', '50'};
feature_number_9b = {'10', '20', '30', '40', '90'};
% Define electrodes 22 and 64
electrodes = {'22', '64'};

% Initialize containers to store mean values
meanSVM = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA = zeros(totalSubjects, length(feature_number_5b), nConditions);

overallMean_5b = zeros(totalSubjects, length(feature_number_5b));
overallMean_9b = zeros(totalSubjects, length(feature_number_9b));

overallMean = zeros(totalSubjects, 1);  % Final overall mean per subject
overallMean_SVM = zeros(totalSubjects, 1);
overallMean_DT = zeros(totalSubjects, 1);
overallMean_LDA = zeros(totalSubjects, 1);

% Initialize containers for storing results for electrodes 22 and 64
meanSVM_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA_22 = zeros(totalSubjects, length(feature_number_5b), nConditions);

meanSVM_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanDT_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);
meanLDA_64 = zeros(totalSubjects, length(feature_number_5b), nConditions);

% Initialize variables to store final overall mean across all classifiers for each electrode
overallMean_22 = zeros(totalSubjects, 1);
overallMean_64 = zeros(totalSubjects, 1);

% Loop through subjects
for numSub = 1:totalSubjects
    % Initialize accumulators for all classifiers
    svmAll = [];
    dtAll = [];
    ldaAll = [];
    svmAll_22 = [];
    dtAll_22 = [];
    ldaAll_22 = [];

    svmAll_64 = [];
    dtAll_64 = [];
    ldaAll_64 = [];
    % Loop through different configurations
    configIndex = 0;
    for el = electrodes
        for band = bands
            % Check if the 9b band is being used
            if strcmp(band{1}, '9b') && isempty(feature_number_9b)
                continue;  % Skip 9b calculations if there are no features for 9b
            end
            configIndex = configIndex + 1;
            feature_number = eval(['feature_number_', band{1}]);

            for iF = 1:length(feature_number)
                % Loop through conditions
                for cond = 1:nConditions
                    % Dynamically construct variable names
                    svmVar = sprintf('SVM_PLV_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);
                    dtVar = sprintf('DT_PLV_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);
                    ldaVar = sprintf('LDA_PLV_%s_%s_feat%s_condition%d', el{1}, band{1}, feature_number{iF}, cond);

                    % Extract values only if they are not empty
                    svmValue = eval(svmVar);
                    dtValue = eval(dtVar);
                    ldaValue = eval(ldaVar);

                    if ~isempty(svmValue{numSub})
                        svmAll = [svmAll, svmValue{numSub}];
                    end
                    if ~isempty(dtValue{numSub})
                        dtAll = [dtAll, dtValue{numSub}];
                    end
                    if ~isempty(ldaValue{numSub})
                        ldaAll = [ldaAll, ldaValue{numSub}];
                    end
                    if ~isempty(svmValue{numSub})
                        if strcmp(el{1}, '22')
                            svmAll_22 = [svmAll_22, svmValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            svmAll_64 = [svmAll_64, svmValue{numSub}];
                        end
                    end

                    if ~isempty(dtValue{numSub})
                        if strcmp(el{1}, '22')
                            dtAll_22 = [dtAll_22, dtValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            dtAll_64 = [dtAll_64, dtValue{numSub}];
                        end
                    end

                    if ~isempty(ldaValue{numSub})
                        if strcmp(el{1}, '22')
                            ldaAll_22 = [ldaAll_22, ldaValue{numSub}];
                        elseif strcmp(el{1}, '64')
                            ldaAll_64 = [ldaAll_64, ldaValue{numSub}];
                        end
                    end
                    % Compute mean for each classifier per feature number, ignoring NaN
                    meanSVM(numSub, iF, cond) = mean(svmAll, 'omitnan');
                    meanDT(numSub, iF, cond) = mean(dtAll, 'omitnan');
                    meanLDA(numSub, iF, cond) = mean(ldaAll, 'omitnan');
                end
                % Compute overall mean per feature number for each subject
                allAccuracies = [svmAll, dtAll, ldaAll];
                if strcmp(band{1}, '5b')
                    overallMean_5b(numSub, iF) = mean(allAccuracies, 'omitnan');
                else
                    overallMean_9b(numSub, iF) = mean(allAccuracies, 'omitnan');
                end
            end

        end
    end

    % Compute overall mean across all configurations, classifiers, and feature indices
    overallMean(numSub) = mean([meanSVM(numSub, :, :), meanDT(numSub, :, :), meanLDA(numSub, :, :)], 'all', 'omitnan');
    % Compute the overall PLV mean per classifier
    overallMean_SVM(numSub) = mean(svmAll, 'omitnan');
    overallMean_DT(numSub) = mean(dtAll, 'omitnan');
    overallMean_LDA(numSub) = mean(ldaAll, 'omitnan');


    % Calculate the mean values for each classifier (SVM, DT, LDA) per feature number for electrode 22 and 64
    for iF = 1:length(feature_number_5b)
        for cond = 1:nConditions
            meanSVM_22(numSub, iF, cond) = mean(svmAll_22, 'omitnan');
            meanDT_22(numSub, iF, cond) = mean(dtAll_22, 'omitnan');
            meanLDA_22(numSub, iF, cond) = mean(ldaAll_22, 'omitnan');

            meanSVM_64(numSub, iF, cond) = mean(svmAll_64, 'omitnan');
            meanDT_64(numSub, iF, cond) = mean(dtAll_64, 'omitnan');
            meanLDA_64(numSub, iF, cond) = mean(ldaAll_64, 'omitnan');
        end
    end

    % Compute the overall mean across all classifiers for each electrode (22 and 64)
    overallMean_22(numSub) = mean([meanSVM_22(numSub, :, :), meanDT_22(numSub, :, :), meanLDA_22(numSub, :, :)], 'all', 'omitnan');
    overallMean_64(numSub) = mean([meanSVM_64(numSub, :, :), meanDT_64(numSub, :, :), meanLDA_64(numSub, :, :)], 'all', 'omitnan');
end


% Display results for all subjects
fprintf('\n=== Summary of Mean Values for All Subjects (PLV) ===\n');
for numSub = 1:totalSubjects
    fprintf('Subject %d:\n', numSub);
    fprintf('  Overall Mean - PLV: %.4f\n\n', overallMean(numSub));
    fprintf('  Overall Mean - SVM: %.4f\n', overallMean_SVM(numSub));
    fprintf('  Overall Mean - DT: %.4f\n', overallMean_DT(numSub));
    fprintf('  Overall Mean - LDA: %.4f\n', overallMean_LDA(numSub));
    % Print overall mean per feature for 5b
    for iF = 1:length(feature_number_5b)
        fprintf('  Overall Mean - 5b (Feature %s): %.4f\n', feature_number_5b{iF}, overallMean_5b(numSub, iF));
    end

    % Print overall mean per feature for 9b
    for iF = 1:length(feature_number_9b)
        fprintf('  Overall Mean - 9b (Feature %s): %.4f\n', feature_number_9b{iF}, overallMean_9b(numSub, iF));
    end
    % Print overall mean for electrode 22
    fprintf('  Overall Mean - Electrode 22: %.4f\n', overallMean_22(numSub));

    % Print overall mean for electrode 64
    fprintf('  Overall Mean - Electrode 64: %.4f\n', overallMean_64(numSub));

    fprintf('\n');  % Print a newline after each subject's results
end
% Saving to a specific path
save(fullfile(myPath, 'overal_PLV.mat'), 'overallMean', 'meanSVM', 'meanDT', 'meanLDA');
save(fullfile(myPath, 'overall_PLV_classifiers.mat'), 'overallMean_SVM', 'overallMean_DT', 'overallMean_LDA');
save(fullfile(myPath, 'overall_PLV_features.mat'), ...
    'overallMean_5b', 'overallMean_9b');
save(fullfile(myPath, 'overall_PLV_channels.mat'), 'overallMean_22', 'overallMean_64');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPH based on FEATURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data from the saved .mat files, ensuring the variable names are different
data_CSP = load(fullfile(myPath, 'overal_CSP.mat'));
data_COH = load(fullfile(myPath, 'overal_COH.mat'));
data_CORR = load(fullfile(myPath, 'overal_CORR.mat'));
data_PLV = load(fullfile(myPath, 'overal_PLV.mat'));

% Extracting the 'overallMean' from each loaded structure
CSP_data = data_CSP.overallMean;  % Example for CSP
COH_data = data_COH.overallMean;  % Example for COH
CORR_data = data_CORR.overallMean; % Example for CORR
PLV_data = data_PLV.overallMean;  % Example for PLV

% Combine the data into a single matrix (each column will correspond to a category)
data = [CSP_data(:), COH_data(:), CORR_data(:), PLV_data(:)];

% Create the box plot
figure;
boxplot(data, 'Labels', {'CSP', 'COH', 'CORR', 'PLV'});
hold on;  % Keep the box plot active for adding scatter

% Overlay scatter plot with individual data points
scatter(ones(size(CSP_data)), CSP_data, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
scatter(2*ones(size(COH_data)), COH_data, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
scatter(3*ones(size(CORR_data)), CORR_data, 'filled', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
scatter(4*ones(size(PLV_data)), PLV_data, 'filled', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm');

% Calculate and display mean and std for each category
categories = {'CSP', 'COH', 'CORR', 'PLV'};
for i = 1:4
    % Calculate mean and std
    category_data = data(:, i);
    category_mean = mean(category_data);
    category_std = std(category_data);

    % Display the mean and std beside each box
    text(i, category_mean, sprintf('Mean: %.2f\nSTD: %.2f', category_mean, category_std), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'black');
end

% Hold off to stop adding to the current plot
hold off;

% Title and labels
title('Comparison Between Features', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Features', 'FontSize', 12);
ylabel('Mean Accuracy (%)', 'FontSize', 12);
saveas(gcf, fullfile(savePath, 'features_comparison.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistical Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all data into one column and create group labels
all_data = [CSP_data(:); COH_data(:); CORR_data(:); PLV_data(:)];
group_labels = [ones(size(CSP_data)); 2*ones(size(COH_data)); 3*ones(size(CORR_data)); 4*ones(size(PLV_data))];

% Perform Kruskal-Wallis test (non-parametric ANOVA alternative)
[p_kw, ~, stats] = kruskalwallis(all_data, group_labels, 'off');
fprintf('Kruskal-Wallis test p-value: %.4f\n', p_kw);

% Perform post-hoc multiple comparisons using Dunn’s test (always performed)
comparisons = multcompare(stats, 'CType', 'dunn-sidak', 'Display', 'off');

% Store p-values for each pairwise comparison
p_values = zeros(4,4);
for i = 1:size(comparisons, 1)
    idx1 = comparisons(i, 1);
    idx2 = comparisons(i, 2);
    p_values(idx1, idx2) = comparisons(i, 6);
    p_values(idx2, idx1) = comparisons(i, 6); % Ensure symmetry
end

% Print formatted p-values for each pair
fprintf('\nPairwise comparison p-values for features:\n');
fprintf('P value between CSP vs COH: %.4f\n', p_values(1,2));
fprintf('P value between CSP vs CORR: %.4f\n', p_values(1,3));
fprintf('P value between CSP vs PLV: %.4f\n', p_values(1,4));
fprintf('P value between COH vs CORR: %.4f\n', p_values(2,3));
fprintf('P value between COH vs PLV: %.4f\n', p_values(2,4));
fprintf('P value between CORR vs PLV: %.4f\n', p_values(3,4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPH based on CLASSIFIERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-allocate arrays to store means for each classifier across subjects
mean_LDA = zeros(totalSubjects, 4);  % 4 feature categories (CSP, COH, CORR, PLV)
mean_SVM = zeros(totalSubjects, 4);
mean_DT = zeros(totalSubjects, 4);
% Initialize arrays to hold the LDA values for each subject
LDA_values = zeros(totalSubjects, 4);  % 4 feature categories (CSP, COH, CORR, PLV)

% Loop through each subject and calculate the means
for subj = 1:totalSubjects
    % Load the data for the current subject
    data_CSP_classifiers = load(fullfile(myPath, 'overall_CSP_classifiers.mat'));
    data_COH_classifiers = load(fullfile(myPath, 'overall_COH_classifiers.mat'));
    data_CORR_classifiers = load(fullfile(myPath, 'overall_CORR_classifiers.mat'));
    data_PLV_classifiers = load(fullfile(myPath, 'overall_PLV_classifiers.mat'));

    % Extract the means for each classifier (LDA, SVM, DT)
    mean_LDA(subj, 1) = data_CSP_classifiers.overallMean_LDA(subj);
    mean_LDA(subj, 2) = data_COH_classifiers.overallMean_LDA(subj);
    mean_LDA(subj, 3) = data_CORR_classifiers.overallMean_LDA(subj);
    mean_LDA(subj, 4) = data_PLV_classifiers.overallMean_LDA(subj);

    mean_SVM(subj, 1) = data_CSP_classifiers.overallMean_SVM(subj);
    mean_SVM(subj, 2) = data_COH_classifiers.overallMean_SVM(subj);
    mean_SVM(subj, 3) = data_CORR_classifiers.overallMean_SVM(subj);
    mean_SVM(subj, 4) = data_PLV_classifiers.overallMean_SVM(subj);

    mean_DT(subj, 1) = data_CSP_classifiers.overallMean_DT(subj);
    mean_DT(subj, 2) = data_COH_classifiers.overallMean_DT(subj);
    mean_DT(subj, 3) = data_CORR_classifiers.overallMean_DT(subj);
    mean_DT(subj, 4) = data_PLV_classifiers.overallMean_DT(subj);
end
% Calculate overall mean for each classifier
all_LDA = mean(mean_LDA, 2);  % Mean across subjects
all_SVM = mean(mean_SVM, 2);
all_DT = mean(mean_DT, 2);

% Create boxplot for each classifier
figure;
data = [all_LDA, all_SVM, all_DT];  % Combine classifier data
boxplot(data, 'Labels', {'LDA', 'SVM', 'DT'}, 'Whisker', Inf);

% Hold on to add scatter plots for individual points
hold on;

% Overlay scatter plot with individual data points for LDA, SVM, and DT
scatter(ones(size(all_LDA)), all_LDA, 100, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
scatter(2*ones(size(all_SVM)), all_SVM, 100, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
scatter(3*ones(size(all_DT)), all_DT, 100, 'filled', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');

% Calculate and display mean and std for each classifier
classifiers = {'LDA', 'SVM', 'DT'};
for i = 1:3
    category_data = data(:, i);  % Data for each classifier
    category_mean = mean(category_data);
    category_std = std(category_data);

    % Display the mean and std beside each box
    text(i, category_mean + category_std + 0.02, ...
        sprintf('Mean: %.2f\nSTD: %.2f', category_mean, category_std), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'black');
end

% Hold off to stop adding to the current plot
hold off;

% Title and labels
title('Comparison Between Classifiers', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Classifiers', 'FontSize', 12);
ylabel('Mean Accuracy (%)', 'FontSize', 12);
grid on;

% Save the figure as an image file
saveas(gcf, fullfile(savePath, 'classifier_comparison.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistical Analysis for Classifiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all data into one column and create group labels
all_LDA_values = all_LDA(:);
all_SVM_values = all_SVM(:);
all_DT_values = all_DT(:);

all_data_classifiers = [all_LDA_values; all_SVM_values; all_DT_values];
group_labels_classifiers = [ones(size(all_LDA_values)); 2*ones(size(all_SVM_values)); 3*ones(size(all_DT_values))];

% Perform Kruskal-Wallis test (non-parametric ANOVA alternative)
[p_kw_classifiers, ~, stats_classifiers] = kruskalwallis(all_data_classifiers, group_labels_classifiers, 'off');
fprintf('Kruskal-Wallis test p-value: %.4f\n', p_kw_classifiers);

% Perform post-hoc multiple comparisons using Dunn’s test (always performed)
comparisons_classifiers = multcompare(stats_classifiers, 'CType', 'dunn-sidak', 'Display', 'off');

% Store p-values for each pairwise comparison
p_values_classifiers = zeros(3,3);
for i = 1:size(comparisons_classifiers, 1)
    idx1 = comparisons_classifiers(i, 1);
    idx2 = comparisons_classifiers(i, 2);
    p_values_classifiers(idx1, idx2) = comparisons_classifiers(i, 6);
    p_values_classifiers(idx2, idx1) = comparisons_classifiers(i, 6); % Ensure symmetry
end

% Print formatted p-values for each pair
fprintf('\nPairwise comparison p-values for Classifiers:\n');
fprintf('P value between LDA vs SVM: %.4f\n', p_values_classifiers(1,2));
fprintf('P value between LDA vs DT: %.4f\n', p_values_classifiers(1,3));
fprintf('P value between SVM vs DT: %.4f\n', p_values_classifiers(2,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPH based on CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data for different feature extraction methods (CSP, COH, CORR, PLV)
data_CSP = load(fullfile(myPath, 'overall_CSP_channels.mat'));
data_COH = load(fullfile(myPath, 'overall_COH_channels.mat'));
data_CORR = load(fullfile(myPath, 'overall_CORR_channels.mat'));
data_PLV = load(fullfile(myPath, 'overall_PLV_channels.mat'));

% Combine data for 22 and 64 channels across all feature methods
overall_22 = [data_CSP.overallMean_22, data_COH.overallMean_22, data_CORR.overallMean_22, data_PLV.overallMean_22];
overall_64 = [data_CSP.overallMean_64, data_COH.overallMean_64, data_CORR.overallMean_64, data_PLV.overallMean_64];

% Compute mean across subjects
mean_22 = mean(overall_22, 2);
mean_64 = mean(overall_64, 2);

% Combine data for boxplot
data = [mean_22, mean_64];

% Create boxplot
figure;
boxplot(data, 'Labels', {'22 Channels', '64 Channels'}, 'Whisker', Inf);

% Hold on to add scatter plots for individual points
hold on;

% Overlay scatter plot with individual data points
scatter(ones(size(mean_22)), mean_22, 100, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
scatter(2*ones(size(mean_64)), mean_64, 100, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

% Display mean and standard deviation
categories = {'22 Channels', '64 Channels'};
for i = 1:2
    category_data = data(:, i);
    category_mean = mean(category_data);
    category_std = std(category_data);

    text(i, category_mean + category_std + 0.02, ...
        sprintf('Mean: %.2f\nSTD: %.2f', category_mean, category_std), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'black');
end

% Hold off to stop adding to the current plot
hold off;

% Title and labels
title('Comparison Between Channels', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Number of Channels', 'FontSize', 12);
ylabel('Mean Accuracy (%)', 'FontSize', 12);
grid on;

% Save the figure
saveas(gcf, fullfile(savePath, 'channel_comparison.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistical Analysis for Channels (22 vs 64 Channels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all data for 22 and 64 channels
all_22_values = mean_22(:);  % Mean values for 22 Channels
all_64_values = mean_64(:);  % Mean values for 64 Channels

% Combine data for channels into one column and create group labels
all_data_channels = [all_22_values; all_64_values];
group_labels_channels = [ones(size(all_22_values)); 2*ones(size(all_64_values))];

% Perform Kruskal-Wallis test (non-parametric ANOVA alternative)
[p_kw_channels, ~, stats_channels] = kruskalwallis(all_data_channels, group_labels_channels, 'off');
fprintf('Kruskal-Wallis test p-value for channels: %.4f\n', p_kw_channels);

% Perform post-hoc multiple comparisons using Dunn’s test (always performed)
comparisons_channels = multcompare(stats_channels, 'CType', 'dunn-sidak', 'Display', 'off');

% Store p-values for pairwise comparison (though there is only one comparison here)
p_values_channels = zeros(2,2);
for i = 1:size(comparisons_channels, 1)
    idx1 = comparisons_channels(i, 1);
    idx2 = comparisons_channels(i, 2);
    p_values_channels(idx1, idx2) = comparisons_channels(i, 6);
    p_values_channels(idx2, idx1) = comparisons_channels(i, 6); % Ensure symmetry
end

% Print formatted p-values for pairwise comparisons
fprintf('\nPairwise comparison p-values for Channels (22 vs 64):\n');
fprintf('P value between 22 Channels vs 64 Channels: %.4f\n', p_values_channels(1,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPH based on Frequency Bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the feature data
data_CSP = load(fullfile(myPath, 'overall_CSP_features.mat'));
data_COH = load(fullfile(myPath, 'overall_COH_features.mat'));
data_CORR = load(fullfile(myPath, 'overall_CORR_features.mat'));
data_PLV = load(fullfile(myPath, 'overall_PLV_features.mat'));

% Concatenate all feature categories
all_features_5b = [data_CSP.overallMean_5b, data_COH.overallMean_5b, data_CORR.overallMean_5b, data_PLV.overallMean_5b];
all_features_9b = [data_CSP.overallMean_9b, data_COH.overallMean_9b, data_CORR.overallMean_9b, data_PLV.overallMean_9b];

% Compute the mean across all five feature values for each subject
subject_mean_5b = mean(all_features_5b, 2);
subject_mean_9b = mean(all_features_9b, 2);

% Create boxplot
figure('Position', [100, 100, 600, 400]); % Adjust figure size
data = [subject_mean_5b, subject_mean_9b];  % Combine data
boxplot(data, 'Labels', {'5 Bands', '9 Bands'}, 'Whisker', Inf, 'Widths', 0.5);

% Hold on to add scatter plots for individual points
hold on;
scatter(ones(size(subject_mean_5b)), subject_mean_5b, 100, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
scatter(2 * ones(size(subject_mean_9b)), subject_mean_9b, 100, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

% Calculate and display mean and std for each category
categories = {'5 Bands', '9 Bands'};
for i = 1:2
    category_data = data(:, i);  % Data for each band set
    category_mean = mean(category_data);
    category_std = std(category_data);

    % Display the mean and std beside each box
    text(i, category_mean + category_std + 0.02, ...
        sprintf('Mean: %.2f\nSTD: %.2f', category_mean, category_std), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'black');
end

% Hold off to stop adding to the current plot
hold off;

% Title and labels with consistent font
title('Comparison Between Frequency Bands', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequency Bands', 'FontSize', 12);
ylabel('Mean Accuracy (%)', 'FontSize', 12);
grid on;

% Save the figure
saveas(gcf, fullfile(savePath, 'bands_comparison.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistical Analysis for Frequency Bands (5 vs 9 Bands)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all data for 5 bands and 9 bands
all_5b_values = subject_mean_5b(:);  % Mean values for 5 Bands
all_9b_values = subject_mean_9b(:);  % Mean values for 9 Bands

% Combine data for frequency bands into one column and create group labels
all_data_bands = [all_5b_values; all_9b_values];
group_labels_bands = [ones(size(all_5b_values)); 2*ones(size(all_9b_values))];

% Perform Kruskal-Wallis test (non-parametric ANOVA alternative)
[p_kw_bands, ~, stats_bands] = kruskalwallis(all_data_bands, group_labels_bands, 'off');
fprintf('Kruskal-Wallis test p-value for frequency bands: %.4f\n', p_kw_bands);

% Perform post-hoc multiple comparisons using Dunn’s test (always performed)
comparisons_bands = multcompare(stats_bands, 'CType', 'dunn-sidak', 'Display', 'off');

% Store p-values for pairwise comparison (though there is only one comparison here)
p_values_bands = zeros(2,2);
for i = 1:size(comparisons_bands, 1)
    idx1 = comparisons_bands(i, 1);
    idx2 = comparisons_bands(i, 2);
    p_values_bands(idx1, idx2) = comparisons_bands(i, 6);
    p_values_bands(idx2, idx1) = comparisons_bands(i, 6); % Ensure symmetry
end

% Print formatted p-values for pairwise comparisons
fprintf('\nPairwise comparison p-values for Frequency Bands (5 vs 9):\n');
fprintf('P value between 5 Bands vs 9 Bands: %.4f\n', p_values_bands(1,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPH based on Feature Number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the feature data
data_CSP = load(fullfile(myPath, 'overall_CSP_features.mat'));
data_COH = load(fullfile(myPath, 'overall_COH_features.mat'));
data_CORR = load(fullfile(myPath, 'overall_CORR_features.mat'));
data_PLV = load(fullfile(myPath, 'overall_PLV_features.mat'));

% Concatenate all feature categories
all_features_5b = cat(3, data_CSP.overallMean_5b, data_COH.overallMean_5b, data_CORR.overallMean_5b, data_PLV.overallMean_5b);
all_features_9b = cat(3, data_CSP.overallMean_9b, data_COH.overallMean_9b, data_CORR.overallMean_9b, data_PLV.overallMean_9b);

% Compute the mean per feature number for each subject
mean_features_5b = mean(all_features_5b, 3); % Average across the 4 feature types
mean_features_9b = mean(all_features_9b, 3);
% Compute the overall feature mean for columns 1 to 4
overall_feature = [(mean_features_5b(:,1:4) + mean_features_9b(:,1:4)) / 2, mean_features_5b(:,5), mean_features_9b(:,5)];

% Define the corresponding feature numbers
feature_numbers = [10, 20, 30, 40, 50, 90];

% Create the boxplot
figure('Position', [100, 100, 600, 400]); % Adjust figure size
boxplot(overall_feature, 'Labels', arrayfun(@num2str, feature_numbers, 'UniformOutput', false), 'Widths', 0.5);

% Hold the plot to overlay scatter points
hold on;

% Add scatter plot: scatter each feature's data points for each subject
[~, num_features] = size(overall_feature); % Get the number of features
for i = 1:num_features
    scatter(repmat(i, size(overall_feature, 1), 1), overall_feature(:,i), 100, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
end

% Compute and display the mean and std for each feature number
for i = 1:num_features
    feature_data = overall_feature(:, i);
    mean_val = mean(feature_data);
    std_val = std(feature_data);

    % Display the mean and std beside each box
    text(i, mean_val + std_val + 0.02, ...
        sprintf('Mean: %.2f\nSTD: %.2f', mean_val, std_val), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'black');
end

% Customize the plot
title('Comparison BetweenFeature Numbers', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Feature Number', 'FontSize', 12);
ylabel('Mean Accuracy (%)', 'FontSize', 12);
grid on;

% Release hold to finish the plot
hold off;
% Save the figure
saveas(gcf, fullfile(savePath, 'feature_number_comparison.png'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistical Analysis for Feature Numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine all data for feature numbers
all_feature_values = reshape(overall_feature, [], 1);  % Flatten the matrix for analysis
feature_labels = repmat(feature_numbers, size(overall_feature, 1), 1);  % Corresponding feature number labels
feature_labels = feature_labels(:);  % Reshape for analysis

% Perform Kruskal-Wallis test (non-parametric ANOVA alternative)
[p_kw_features, ~, stats_features] = kruskalwallis(all_feature_values, feature_labels, 'off');
fprintf('Kruskal-Wallis test p-value for feature numbers: %.4f\n', p_kw_features);

% Perform post-hoc multiple comparisons using Dunn’s test (always performed)
comparisons_features = multcompare(stats_features, 'CType', 'dunn-sidak', 'Display', 'off');

% Store p-values for pairwise comparison
p_values_features = zeros(length(feature_numbers), length(feature_numbers));
for i = 1:size(comparisons_features, 1)
    idx1 = comparisons_features(i, 1);
    idx2 = comparisons_features(i, 2);
    p_values_features(idx1, idx2) = comparisons_features(i, 6);
    p_values_features(idx2, idx1) = comparisons_features(i, 6);  % Ensure symmetry
end

% Print formatted p-values for pairwise comparisons
fprintf('\nPairwise comparison p-values for Feature Numbers:\n');
for i = 1:length(feature_numbers)
    for j = i+1:length(feature_numbers)
        fprintf('P value between Feature %d vs Feature %d: %.4f\n', feature_numbers(i), feature_numbers(j), p_values_features(i,j));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function features = getFeatures(band)
feature_number_5b   = {'10','20','30','40','50'};
feature_number_9b   = {'10','20','30','40','90'};
if strcmp(band, '5b')
    features = feature_number_5b;
else
    features = feature_number_9b;
end
end
