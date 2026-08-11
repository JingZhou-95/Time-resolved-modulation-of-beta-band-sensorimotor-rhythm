%% EEG Preprocessing Pipeline
% Author: Jing
%
% Description:
% This script preprocesses EEG data recorded during a repeated motor task.
%
% Main processing steps:
%   1. Extract and filter EEG signals
%   2. Extract torque and EMG signals
%   3. Create task-onset event markers
%   4. Retain the analysis periods from each trial
%   5. Import the data into EEGLAB
%   6. Manually reject large artifacts
%   7. Perform Independent Component Analysis (ICA)
%   8. Manually remove artifact-related ICA components
%   9. Extract task-locked epochs
%  10. Save the preprocessed data
%
% Required variables:
%   data : Raw data matrix [samples x channels]
%   fs   : Sampling frequency [Hz]
%
% Expected data structure:
%   Column 1      : Torque / MVC signal
%   Columns 2-30  : EEG channels (29 channels)
%   Column 31     : Tibialis anterior (TA) EMG
%
% Required software:
%   EEGLAB
%
% Required custom functions:
%   bp.m : Band-pass filter
%   bs.m : Band-stop filter


%% Configuration

% Set the EEGLAB directory before running this script.
eeglabPath = '/path/to/eeglab';
addpath(genpath(eeglabPath));

% Channel indices in the original data matrix.
TORQUE_CH = 1;
EEG_CH    = 2:30;
TA_CH     = 31;

% Experimental timing.
trialDuration = 25;      % Trial duration [s]
taskOnsetTime = 18;      % Task onset from the beginning of each trial [s]

% Segment timing used in the original analysis.
segmentStartTime = 9;    % Segment start from trial onset [s]
segmentEndOffset = 4;    % Offset from the beginning of the next trial [s]

% Epoch window relative to task onset.
epochWindow = [-8 12];   % [s]


%% Check input data

assert(exist('data', 'var') == 1, ...
    'Variable "data" was not found in the workspace.');

assert(exist('fs', 'var') == 1, ...
    'Variable "fs" was not found in the workspace.');

assert(size(data, 2) >= 32, ...
    'The input data must contain at least 32 channels.');

numSamples = size(data, 1);

% Number of complete 25-s trials.
trialnum = floor(numSamples / (trialDuration * fs));

fprintf('Number of trials detected: %d\n', trialnum);


%% EEG preprocessing

% Extract EEG channels.
EEGdata = data(:, EEG_CH);

% Convert EEG signal amplitude.
EEGdata = EEGdata * 100;

% Remove linear trends.
EEGdata = detrend(EEGdata);

% Band-pass filter: 1-100 Hz.
EEGdata = bp(EEGdata, fs, 3, 1, 100);

% Remove 50-Hz power-line noise.
EEGdata = bs(EEGdata, fs, 3, 49, 51);

% Remove the 100-Hz harmonic.
EEGdata = bs(EEGdata, fs, 3, 99, 101);

% EEGLAB data format: channels x samples.
EEGdata = EEGdata';


%% Extract torque and EMG signals

torque = data(:, TORQUE_CH)';
TA     = data(:, TA_CH)';

%% Create task-onset event markers

% A single-sample pulse is placed at task onset in each trial.

eventData = zeros(1, numSamples);

for trial = 1:trialnum

    eventSample = round( ...
        taskOnsetTime * fs + ...
        (trial - 1) * trialDuration * fs + 1);

    eventData(eventSample) = 1;

end


%% Combine signals

% Channel structure used for EEGLAB:
%
%   Channel 1      : Torque
%   Channels 2-30  : EEG
%   Channel 31     : TA EMG
%   Channel 33     : Task-onset event marker

data4eeglabFull = [
    torque;
    EEGdata;
    TA;
    eventData
];


%% Retain analysis periods

% Extract the same trial segments used in the original analysis.
%
% For each trial:
%   start = 9 s after the beginning of the current trial
%
% For all trials except the final trial:
%   end = 4 s after the beginning of the next trial
%
% For the final trial:
%   all remaining samples are retained.

segments = cell(1, trialnum);

for trial = 1:trialnum

    startSample = round( ...
        segmentStartTime * fs + 1 + ...
        (trial - 1) * trialDuration * fs);

    if trial == trialnum

        endSample = size(data4eeglabFull, 2);

    else

        endSample = round( ...
            segmentEndOffset * fs + ...
            trial * trialDuration * fs);

    end

    segments{trial} = ...
        data4eeglabFull(:, startSample:endSample);

end

data4eeglab = [segments{:}];


%% Import data into EEGLAB

EEG = pop_importdata( ...
    'setname',    'EEG_preprocessing', ...
    'data',       'data4eeglab', ...
    'dataformat', 'float32', ...
    'srate',      fs);

EEG = eeg_checkset(EEG);


%% Import task-onset events

% Channel 33 contains single-sample task-onset pulses.
%
% Only the leading edge is extracted so that each pulse produces
% one TaskOnset event.
%
% The event-marker channel is removed after event extraction.

EEG = pop_chanevent( ...
    EEG, ...
    33, ...
    'edge',     'leading', ...
    'edgelen',  1, ...
    'delchan',  'on', ...
    'delevent', 'on', ...
    'nbtype',   1, ...
    'typename', 'TaskOnset');

EEG = eeg_checkset(EEG, 'eventconsistency');


%% Prepare EEGLAB dataset

% These variables are required for interactive EEGLAB operations.

ALLEEG = EEG;
CURRENTSET = 1;


%% Manual artifact rejection

% Display the continuous data.
%
% Manually mark periods containing large movement artifacts,
% electrode artifacts, or other obvious non-physiological signals.
%
% Close the EEGLAB plotting window after completing rejection.

pop_eegplot(EEG, 1, 0, 1);

fig = gcf;
uiwait(fig);

EEG = eeg_checkset(EEG);


%% Independent Component Analysis (ICA)

% Run extended Infomax ICA using EEG channels only.
%
% Channel 1:
%   Torque
%
% Channels 2-30:
%   EEG
%
% Channels 31-32:
%   EMG
%
% Therefore, only channels 2-30 are included in ICA.

EEG = pop_runica( ...
    EEG, ...
    'icatype',   'runica', ...
    'chanind',   2:30, ...
    'extended',  1, ...
    'interrupt', 'on');

EEG = eeg_checkset(EEG);


%% Inspect ICA components

% Inspect ICA component activations before selecting components
% associated with eye movements, muscle activity, or other artifacts.

pop_eegplot(EEG, 0, 1, 1);

fig = gcf;
uiwait(fig);


%% Remove artifact-related ICA components

componentInput = inputdlg( ...
    'Enter ICA components to reject (e.g., [1 3 5]):', ...
    'ICA Component Rejection');

if ~isempty(componentInput) && ~isempty(componentInput{1})

    rejectedComponents = str2num(componentInput{1}); %#ok<ST2NM>

    EEG = pop_subcomp( ...
        EEG, ...
        rejectedComponents, ...
        1);

end

EEG = eeg_checkset(EEG);


%% Epoch extraction

% Extract all epochs from -8 s to +12 s relative to task onset.
%
% No odd/even epoch selection is performed.

[EEG_epoch, epochidx] = pop_epoch( ...
    EEG, ...
    {'TaskOnset'}, ...
    epochWindow, ...
    'newname',   'Task_Epochs', ...
    'epochinfo', 'yes');

EEG_epoch = eeg_checkset(EEG_epoch);


%% Store preprocessed data

% Continuous data after manual artifact rejection and ICA.
EEGdata = EEG.data;

% All task-locked epochs after ICA.
% Dimensions:
%   channels x time points x epochs
EEGEpoch = EEG_epoch.data;


%% Concatenate epochs

% Concatenate all epochs along the time dimension.
%
% Output dimensions:
%   channels x (time points * number of epochs)

dataarej = reshape( ...
    EEGEpoch, ...
    size(EEGEpoch, 1), ...
    []);


%% Data description

ExplanationText = [ ...
    "EEGdata contains continuous data after manual artifact rejection " ...
    "and ICA. EEGEpoch contains all task-locked epochs extracted from " ...
    "-8 s to +12 s relative to task onset. No odd/even epoch selection " ...
    "is applied."
];


%% Save preprocessed data

% Specify the output directory and filename as needed.

outputDir = 'CleanData';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

outputFile = fullfile( ...
    outputDir, ...
    'EEG_preprocessed.mat');

save( ...
    outputFile, ...
    'EEGdata', ...
    'EEGEpoch', ...
    'dataarej', ...
    'fs', ...
    'epochidx', ...
    'ExplanationText', ...
    '-v7.3');

fprintf('Preprocessed data saved to:\n%s\n', outputFile);
