%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing before extracting events
% Mona Hejazi
% mhejazi@mun.ca+
% MUN - Memorial University of Newfoundland
% August 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Step 3: Pre-processing of EEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
warning off;


% Add the EEGLAB path
addpath('E:\1. PhD project\2022 Files\1. ProjectMaterials\EEG\EEG-DATA\example\eeglab2022.1');

% Load EEGLAB
eeglab;
% Define the total number of subjects
totalSubjects = 29;

% List of subjects to exclude
excludeSubjects = [23,25];
warning('off', 'MATLAB:colon:operandsNotInteger');

for numSubjects = 1
%     :totalSubjects
    % Skip excluded subjects
    if any(numSubjects == excludeSubjects)
        continue;
    end

    % start EEGLAB
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

    subjectFolder = sprintf('Sub-%03d', numSubjects);
    mypath = fullfile('E:\1. PhD project\2022 Files\1. ProjectMaterials\EEG\EEG-DATA', subjectFolder);
    EEGpath = 'E:\1. PhD project\2022 Files\1. ProjectMaterials\EEG\EEG-DATA\example\eeglab2023.1\eeglab.m';
    addpath(EEGpath)
    STUDY = [];
    ALLEEG = [];
    EEG = [];
    CURRENTSET = [];
    allfiles = dir(fullfile(mypath, '*.vhdr'));
    for n = 1:length(allfiles)
        loadName = allfiles(n).name;
        dataName = loadName(1:end-5);
        % Step 2: Import data.
        EEG = pop_loadbv(mypath, loadName);
        EEG.setname = dataName;
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', dataName);
        % Downsampling
        EEG = pop_resample(EEG, 256);
        % Filtering
        EEG = pop_eegfiltnew(EEG, 'locutoff', 59, 'hicutoff', 61, 'revfilt', 1, 'plotfreqz', 0);

        % Add reference electrode (FCz)
        % Add the reference channel (FCz) using pop_chanedit()
        locpath = 'E:\1. PhD project\2022 Files\1. ProjectMaterials\EEG\EEG-DATA\example\eeglab2023.0_old\plugins\dipfit\standard_BEM\elec\standard_1005.elc';
        EEG = pop_chanedit(EEG, 'append', 63, 'changefield', {64, 'labels', 'FCz'}, 'lookup', locpath, 'setref', {'1:63', 'FCz'});
        % Re-read the electrode locations from the elc file
        EEG.chanlocs = readlocs(locpath);
        % Specify the reference electrode label
        refLabel = 'FCz';
        % Find the index of FCz channel
        refIndex = find(strcmp({EEG.chanlocs.labels}, refLabel));
        % Re-reference to the average reference
        EEG = pop_reref(EEG, [], 'refloc', EEG.chanlocs(refIndex), 'keepref', 'on');
        % Store the data
        [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
    end
    eeglab redraw
    % Save the latest version of ALLEEG in mypath
    save(fullfile(mypath, 'ALLEEG_latest.mat'), 'ALLEEG');

end
