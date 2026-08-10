%% SMR/ CMC / Force Analysis
%
% Main analyses:
%   1. Cz surface-Laplacian EEG
%   2. EEG time-frequency power
%   3. Beta-band power change: MRBD, SRBS, and Delta SRBS-MRBD
%   4. Corticomuscular coherence (CMC) using TA EMG
%   5. CMC-selected 3-Hz EEG power change
%   6. Torque performance using RMSE
%
% Supported task types:
%   Experiment 1:
%       '10pct'  : constant target at 10 %MVC
%       '15pct'  : constant target at 15 %MVC
%
%   Experiment 2:
%       'St'     : constant target at 10 %MVC
%       '1sin'   : one sinusoidal target
%       '3sin'   : three sinusoidal segments
%       '4sin'   : four sinusoidal segments
%
% Required variables in the workspace:
%   data        : preprocessed data matrix
%   IEMG_ratio  : torque-to-%MVC conversion factor
%
% Expected data columns used here:
%   Column 1      : Torque
%   Column 2      : Cz
%   Column 3      : FCz
%   Column 4      : CPz
%   Column 5      : C1
%   Column 6      : C2
%   Column 31     : Tibialis anterior (TA) EMG
%
% Required custom function:
%   bs.m         : band-stop filter

%% Task condition
% Change this value for the dataset currently being analyzed:
% '10pct', '15pct', 'St', '1sin', '3sin', or '4sin'

taskType = '10pct';

validTaskTypes = {'10pct', '15pct', 'St', '1sin', '3sin', '4sin'};
assert(any(strcmp(taskType, validTaskTypes)), ...
    'Invalid taskType. Use: 10pct, 15pct, St, 1sin, 3sin, or 4sin.');

%% Parameters
fs = 1000;
nfft = 1000;
overlap = 0;

trialDuration = 20;      % Duration of one trial [s]
analysisDuration = 18;   % Analyzed duration within each trial [s]
maxTrials = 24;

% Middle 5-s interval used for torque-performance analysis.
% With the current time axis, this corresponds to +1 to +6 s
% relative to task onset.
analysisStartSec = 9;
analysisEndSec = 14;
analysisIdx = analysisStartSec * fs + 1 : analysisEndSec * fs;

% Reproducible random trial selection when more than 24 trials exist.
rng(42);

%% Keep at most 24 trials
numSamples = size(data, 1);
trialnumInitial = floor(numSamples / (trialDuration * fs));

if trialnumInitial > maxTrials
    trialsToRemove = sort( ...
        randperm(trialnumInitial, trialnumInitial - maxTrials), ...
        'descend');

    for i = trialsToRemove
        data(trialDuration * fs * (i-1) + 1 : ...
             trialDuration * fs * i, :) = [];
    end
end

%% Extract EEG channels and calculate Cz surface Laplacian
EEG_Cz  = data(:, 2);
EEG_FCz = data(:, 3);
EEG_CPz = data(:, 4);
EEG_C1  = data(:, 5);
EEG_C2  = data(:, 6);

lapEEG_Cz = EEG_Cz - ...
    (EEG_FCz + EEG_CPz + EEG_C1 + EEG_C2) / 4;

trialnum = floor(length(lapEEG_Cz) / (trialDuration * fs));

%% EEG time-frequency analysis
% A 1-s window is shifted in 50-ms steps.
% The 18-s analysis period produces 360 time windows.

numTimeWindows = 360;
windowStep = round(0.05 * fs);
windowLength = fs;

ConCz = zeros(trialnum * fs, numTimeWindows);

for trial = 1:trialnum
    for i = 1:numTimeWindows
        startSample = 1 + (i-1) * windowStep + ...
            (trial-1) * trialDuration * fs;
        endSample = startSample + windowLength - 1;

        ConCz(1 + (trial-1)*fs : trial*fs, i) = ...
            lapEEG_Cz(startSample:endSample);
    end
end

pEEGConCz = zeros(nfft/2 + 1, numTimeWindows);

for i = 1:numTimeWindows
    pEEGConCz(:, i) = pwelch( ...
        ConCz(:, i), hanning(nfft), overlap, nfft, fs);
end

%% EEG time-frequency map
figure(1)
imagesc(pEEGConCz(1:51, :), [0 0.1]);
colormap(parula);
colorbar;

set(gca, ...
    'YTick', [1 10 20 30 40 50], ...
    'YTickLabel', [0 10 20 30 40 50], ...
    'XTick', [1 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [1 50], ...
    'XLim', [1 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('Frequency (Hz)', 'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', 'FontName', 'Arial', 'FontSize', 20);

%% Torque preprocessing
Force = data(:, 1);

Force = bs(Force, fs, 3, 49, 51);
Force = bs(Force, fs, 3, 99, 101);
Force = bs(Force, fs, 3, 149, 151);
Force = bs(Force, fs, 3, 199, 201);
Force = bs(Force, fs, 3, 249, 251);
Force = bs(Force, fs, 3, 299, 301);
Force = bs(Force, fs, 3, 349, 351);
Force = bs(Force, fs, 3, 399, 401);
Force = bs(Force, fs, 3, 449, 451);

tormvc = abs(Force) * abs(IEMG_ratio(1,1));
tormvc = lowpass(tormvc, 30, fs);
tormvc = movmean(tormvc, 200);

%% TA EMG preprocessing
TA = (data(:, 31) - mean(data(:, 31))) * 1000;

[b50,  a50]  = butter(3, [49  51] / (fs/2), 'stop');
[b100, a100] = butter(3, [99  101] / (fs/2), 'stop');
[b150, a150] = butter(3, [149 151] / (fs/2), 'stop');
[b200, a200] = butter(3, [199 201] / (fs/2), 'stop');
[b250, a250] = butter(3, [249 251] / (fs/2), 'stop');
[b300, a300] = butter(3, [299 301] / (fs/2), 'stop');
[b350, a350] = butter(3, [349 351] / (fs/2), 'stop');
[b400, a400] = butter(3, [399 401] / (fs/2), 'stop');
[b450, a450] = butter(3, [449 451] / (fs/2), 'stop');

TA = filtfilt(b50,  a50,  TA);
TA = filtfilt(b100, a100, TA);
TA = filtfilt(b150, a150, TA);
TA = filtfilt(b200, a200, TA);
TA = filtfilt(b250, a250, TA);
TA = filtfilt(b300, a300, TA);
TA = filtfilt(b350, a350, TA);
TA = filtfilt(b400, a400, TA);
TA = filtfilt(b450, a450, TA);

TA = highpass(TA, 2, fs);
% Rectified TA EMG.
rTA = abs(TA);

%% Divide continuous signals into 18-s trials
dtor = zeros(analysisDuration * fs, trialnum);
drTA = zeros(analysisDuration * fs, trialnum);
dTA  = zeros(analysisDuration * fs, trialnum);
dCz  = zeros(analysisDuration * fs, trialnum);

for trial = 1:trialnum
    startSample = 1 + (trial-1) * trialDuration * fs;
    endSample = analysisDuration * fs + ...
        (trial-1) * trialDuration * fs;

    dtor(:, trial) = tormvc(startSample:endSample);
    drTA(:, trial) = rTA(startSample:endSample);
    dTA(:, trial)  = TA(startSample:endSample);
    dCz(:, trial)  = lapEEG_Cz(startSample:endSample);
end

%% Extract the fixed middle 5-s interval
% The previous "most stable 3-s" search is not used.
% The same fixed 5-s interval is extracted from every trial.

exdtor = dtor(analysisIdx, :);
exrTA  = drTA(analysisIdx, :);
exTA   = dTA(analysisIdx, :);
exCz   = dCz(analysisIdx, :);

% Concatenated 5-s data across trials.
coexdtor = exdtor(:);
coexrTA  = exrTA(:);
coexTA   = exTA(:);
coexCz   = exCz(:);

%% Beta-band power change: MRBD and SRBS
% Baseline power is calculated from time windows 81:101.
% Beta-band power is integrated over frequency bins 16:36.
%
% MRBD:
%   Minimum beta-band power change during the predefined task interval.
%
% SRBS:
%   Maximum beta-band power change during the predefined post-MRBD interval.
%
% Delta SRBS-MRBD:
%   SRBS - MRBD.

betaBaselineValues = zeros(21, 1);

for n = 81:101
    betaBaselineValues(n-80) = ...
        trapz(pEEGConCz(16:36, n));
end

betaBaseline = mean(betaBaselineValues);

betaPower = zeros(numTimeWindows, 1);

for n = 1:numTimeWindows
    betaPower(n) = trapz(pEEGConCz(16:36, n));
end

betaPowerChange = betaPower / betaBaseline * 100 - 100;

[MRBD, mrbdLocalIdx] = min(betaPowerChange(101:181));
MRBD_idx = mrbdLocalIdx + 100;

[SRBS, srbsLocalIdx] = max(betaPowerChange(150:275));
SRBS_idx = srbsLocalIdx + 149;

Delta_SRBS_MRBD = SRBS - MRBD;

figure(2)
plot(betaPowerChange, 'LineWidth', 3);

set(gca, ...
    'YTick', [-100 -50 0 50 100 150], ...
    'XTick', [0 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [-100 150], ...
    'XLim', [0 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('Beta power change (%)', ...
    'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 20);

%% Time-frequency CMC
% Rectified and non-rectified TA EMG are analyzed separately.

ConTA  = zeros(trialnum * fs, numTimeWindows);
ConrTA = zeros(trialnum * fs, numTimeWindows);

for trial = 1:trialnum
    for i = 1:numTimeWindows
        startSample = 1 + (i-1) * windowStep + ...
            (trial-1) * trialDuration * fs;
        endSample = startSample + windowLength - 1;

        ConrTA(1 + (trial-1)*fs : trial*fs, i) = ...
            rTA(startSample:endSample);

        ConTA(1 + (trial-1)*fs : trial*fs, i) = ...
            TA(startSample:endSample);
    end
end

rtCoh = zeros(nfft/2 + 1, numTimeWindows);
tCoh  = zeros(nfft/2 + 1, numTimeWindows);

for i = 1:numTimeWindows
    [rtCoh(:, i), rF] = mscohere( ...
        ConCz(:, i), ConrTA(:, i), ...
        hanning(nfft), overlap, nfft, fs);

    [tCoh(:, i), F] = mscohere( ...
        ConCz(:, i), ConTA(:, i), ...
        hanning(nfft), overlap, nfft, fs);
end

%% Rectified-EMG CMC time-frequency map
figure(3)
imagesc(rtCoh(1:51, :), [0 0.5]);
colormap(parula);
colorbar;

set(gca, ...
    'YTick', [1 10 20 30 40 50], ...
    'YTickLabel', [0 10 20 30 40 50], ...
    'XTick', [1 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [1 50], ...
    'XLim', [1 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('Frequency (Hz)', ...
    'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 20);

%% Non-rectified-EMG CMC time-frequency map
figure(4)
imagesc(tCoh(1:51, :), [0 0.5]);
colormap(parula);
colorbar;

set(gca, ...
    'YTick', [1 10 20 30 40 50], ...
    'YTickLabel', [0 10 20 30 40 50], ...
    'XTick', [1 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [1 50], ...
    'XLim', [1 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('Frequency (Hz)', ...
    'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 20);

%% Beta-band CMC
rbeCMC = zeros(numTimeWindows, 1);
beCMC  = zeros(numTimeWindows, 1);

for n = 1:numTimeWindows
    rbeCMC(n) = trapz(rtCoh(16:36, n)) / 21;
    beCMC(n)  = trapz(tCoh(16:36, n)) / 21;
end

rbetaCMCinc = max(rbeCMC(161:300));
betaCMCinc  = max(beCMC(161:300));

figure(5)
plot(rbeCMC, 'LineWidth', 3);

set(gca, ...
    'YTick', [0 0.05 0.1 0.2 0.3], ...
    'XTick', [0 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [0 0.3], ...
    'XLim', [0 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('Beta CMC', ...
    'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 20);

figure(6)
plot(beCMC, 'LineWidth', 3);

set(gca, ...
    'YTick', [0 0.05 0.1 0.2 0.3], ...
    'XTick', [0 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [0 0.3], ...
    'XLim', [0 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('Beta CMC', ...
    'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 20);

%% Maximum 3-Hz CMC band: rectified EMG
meanCMC = zeros(21, numTimeWindows);

for i = 16:36
    meanCMC(i-15, :) = mean(rtCoh(i-1:i+1, :), 1);
end

[inc, FmaxRow] = ...
    max(max(meanCMC(:, 161:300), [], 2));

sanHzCMC = meanCMC(FmaxRow, :)';

figure(7)
plot(meanCMC(FmaxRow, :), 'LineWidth', 3);

set(gca, ...
    'YTick', [0 0.1 0.2 0.3 0.4 0.5], ...
    'XTick', [0 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [0 0.5], ...
    'XLim', [0 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('CMC', ...
    'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 20);

% Convert the selected row to the first frequency-bin index
% of the corresponding 3-Hz band.
Fmax = FmaxRow + 14;

%% Maximum 3-Hz CMC band: non-rectified EMG
meannCMC = zeros(21, numTimeWindows);

for i = 16:36
    meannCMC(i-15, :) = mean(tCoh(i-1:i+1, :), 1);
end

[ninc, nFmaxRow] = ...
    max(max(meannCMC(:, 161:300), [], 2));

nsanHznCMC = meannCMC(nFmaxRow, :)';

figure(8)
plot(meannCMC(nFmaxRow, :), 'LineWidth', 3);

set(gca, ...
    'YTick', [0 0.1 0.2 0.3 0.4 0.5], ...
    'XTick', [0 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [0 0.5], ...
    'XLim', [0 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('CMC', ...
    'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 20);

nFmax = nFmaxRow + 14;

%% Power change in the rectified-EMG CMC-selected 3-Hz band
powerBaseline = mean(pEEGConCz(1:51, 81:101), 2);

powerChangeMap = ...
    pEEGConCz(1:51, :) ./ powerBaseline * 100 - 100;

cmcPowerChange = ...
    mean(powerChangeMap(Fmax:Fmax+2, :), 1);

cmcMRBD = min(cmcPowerChange(101:201));
cmcSRBS = max(cmcPowerChange(161:280));
cmcDelta_SRBS_MRBD = cmcSRBS - cmcMRBD;

figure(9)
plot(cmcPowerChange, 'LineWidth', 3);

set(gca, ...
    'YTick', [-150 -100 -50 0 50 100], ...
    'XTick', [0 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [-150 100], ...
    'XLim', [0 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('Power change (%)', ...
    'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 20);

cmcPowerChange = cmcPowerChange';

%% Power change in the non-rectified-EMG CMC-selected 3-Hz band
ncmcPowerChange = ...
    mean(powerChangeMap(nFmax:nFmax+2, :), 1);

ncmcMRBD = min(ncmcPowerChange(101:201));
ncmcSRBS = max(ncmcPowerChange(161:280));
ncmcDelta_SRBS_MRBD = ncmcSRBS - ncmcMRBD;

figure(10)
plot(ncmcPowerChange, 'LineWidth', 3);

set(gca, ...
    'YTick', [-150 -100 -50 0 50 100], ...
    'XTick', [0 100 160 300 360], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [-150 100], ...
    'XLim', [0 360], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('Power change (%)', ...
    'FontName', 'Arial', 'FontSize', 20);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 20);

ncmcPowerChange = ncmcPowerChange';

%% Torque plot
meanTorque = mean(dtor, 2);
trialTime = 1:analysisDuration*fs;

figure(11)
plot(trialTime, dtor, 'LineWidth', 0.5);
hold on;
plot(trialTime, meanTorque, 'LineWidth', 3);

set(gca, ...
    'XTick', [0 5000 8000 15000 18000], ...
    'XTickLabel', [-8 -3 0 7 10], ...
    'YLim', [0 30], ...
    'XLim', [0 18000], ...
    'YDir', 'normal', ...
    'FontSize', 20);

ylabel('Torque (%MVC)', ...
    'FontName', 'Arial', 'FontSize', 24);
xlabel('Time (s)', ...
    'FontName', 'Arial', 'FontSize', 24);

%% Torque performance: RMSE
% A less-smoothed torque signal is used for performance evaluation,
% consistent with the original performance section.

ptormvc = abs(Force) * abs(IEMG_ratio(1,1));
ptormvc = lowpass(ptormvc, 50, fs);
ptormvc = movmean(ptormvc, 5);

pdtor = zeros(analysisDuration * fs, trialnum);

for trial = 1:trialnum
    startSample = 1 + (trial-1) * trialDuration * fs;
    endSample = analysisDuration * fs + ...
        (trial-1) * trialDuration * fs;

    pdtor(:, trial) = ptormvc(startSample:endSample);
end

% Use the same fixed middle 5-s interval for all tasks.
performanceTorque = pdtor(analysisIdx, :);

% Generate the task-specific 7-s target trajectory.
fullTargetTorque = generateTargetTorque(taskType, fs);

% The analyzed torque interval corresponds to 1-6 s
% of the 7-s target trajectory.
target5s = fullTargetTorque(fs+1 : 6*fs);

% RMSE is calculated separately for each trial.
targetMatrix = repmat(target5s, 1, trialnum);
rmse_each = sqrt(mean((performanceTorque - targetMatrix).^2, 1));

% Mean RMSE across trials for this participant/condition.
rmse_mean = mean(rmse_each);

% Structure compatible with the downstream RMSE analysis workflow.
RMSE_raw = struct;
RMSE_raw.taskType = taskType;
RMSE_raw.rmse_each = rmse_each(:);
RMSE_raw.rmse_mean = rmse_mean;
RMSE_raw.target5s = target5s;
RMSE_raw.performanceTorque = performanceTorque;

%% Collect main outcome variables
results = struct;

results.taskType = taskType;

results.rbetaCMCinc = rbetaCMCinc;
results.betaCMCinc = betaCMCinc;

results.MRBD = MRBD;
results.SRBS = SRBS;
results.Delta_SRBS_MRBD = Delta_SRBS_MRBD;
results.MRBD_idx = MRBD_idx;
results.SRBS_idx = SRBS_idx;

results.Fmax = Fmax;
results.maxRectifiedCMC = inc;
results.cmcMRBD = cmcMRBD;
results.cmcSRBS = cmcSRBS;
results.cmcDelta_SRBS_MRBD = cmcDelta_SRBS_MRBD;

results.nFmax = nFmax;
results.maxNonRectifiedCMC = ninc;
results.ncmcMRBD = ncmcMRBD;
results.ncmcSRBS = ncmcSRBS;
results.ncmcDelta_SRBS_MRBD = ncmcDelta_SRBS_MRBD;

results.rmse_each = rmse_each;
results.rmse_mean = rmse_mean;

% Compact numeric vector for legacy downstream code.
c = [ ...
    rbetaCMCinc, ...
    betaCMCinc, ...
    MRBD, ...
    SRBS, ...
    Delta_SRBS_MRBD, ...
    Fmax, ...
    inc, ...
    cmcMRBD, ...
    cmcSRBS, ...
    cmcDelta_SRBS_MRBD, ...
    nFmax, ...
    ninc, ...
    ncmcMRBD, ...
    ncmcSRBS, ...
    ncmcDelta_SRBS_MRBD, ...
    rmse_mean];

%% Save results and figures
% Results are saved relative to this script.
% Each task has its own output folder to avoid overwriting other tasks.

scriptDir = fileparts(mfilename('fullpath'));

if isempty(scriptDir)
    scriptDir = pwd;
end

outputRoot = fullfile(scriptDir, 'analysis_output');
outputDir = fullfile(outputRoot, taskType);
figureDir = fullfile(outputDir, 'figures');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

% Save the complete analysis output.
save(fullfile(outputDir, 'analysis_results.mat'), ...
    'results', ...
    'c', ...
    'RMSE_raw', ...
    'betaPowerChange', ...
    'cmcPowerChange', ...
    'ncmcPowerChange', ...
    'rbeCMC', ...
    'beCMC', ...
    'rtCoh', ...
    'tCoh', ...
    'pEEGConCz', ...
    'target5s');

% Save RMSE separately for direct use in the downstream RMSE script.
save(fullfile(outputDir, ['RMSE_' taskType '.mat']), ...
    'RMSE_raw');

% Figure files have no internal plot titles.
figureNames = { ...
    'EEG_TimeFrequency', ...
    'MRBD_SRBS_BetaPower', ...
    'Time_CMC_RectifiedEMG', ...
    'Time_CMC_NonRectifiedEMG', ...
    'Beta_CMC_RectifiedEMG', ...
    'Beta_CMC_NonRectifiedEMG', ...
    'Max_CMC_3Hz_RectifiedEMG', ...
    'Max_CMC_3Hz_NonRectifiedEMG', ...
    'CMCBand_MRBD_SRBS_RectifiedEMG', ...
    'CMCBand_MRBD_SRBS_NonRectifiedEMG', ...
    'Torque'};

for figNum = 1:numel(figureNames)
    saveas(figure(figNum), ...
        fullfile(figureDir, [figureNames{figNum} '.jpg']));

    saveas(figure(figNum), ...
        fullfile(figureDir, [figureNames{figNum} '.fig']));
end

%% Local function: task-specific target trajectory
function targetTorque = generateTargetTorque(taskType, fs)
%GENERATETARGETTORQUE Generate the 7-s torque target for each task.
%
% All outputs are column vectors with a total duration of 7 s.
%
% 10pct:
%   Constant 10 %MVC.
%
% 15pct:
%   Constant 15 %MVC.
%
% St:
%   Constant 10 %MVC.
%
% 1sin:
%   7-s sinusoid with amplitude 2 %MVC, offset 10 %MVC,
%   and period 7 s.
%
% 3sin:
%   A 2.3-s sinusoidal segment repeated three times,
%   followed by the first 0.1 s of the same waveform.
%
% 4sin:
%   Four consecutive sinusoidal segments with periods
%   2.3 s, 1.2 s, 2.3 s, and 1.2 s.
%
% The sine-wave amplitude is 2 %MVC and the offset is 10 %MVC.

    switch taskType

        case '10pct'
            targetTorque = 10 * ones(7*fs, 1);

        case '15pct'
            targetTorque = 15 * ones(7*fs, 1);

        case 'St'
            targetTorque = 10 * ones(7*fs, 1);

        case '1sin'
            t = (0 : 1/fs : 7 - 1/fs)';
            targetTorque = ...
                2 * sin(2*pi*(1/7) * t) + 10;

        case '3sin'
            period = 2.3;

            tSegment = (0 : 1/fs : period - 1/fs)';
            oneSegment = ...
                2 * sin(2*pi/period * tSegment) + 10;

            tTail = (0 : 1/fs : 0.1 - 1/fs)';
            tail = ...
                2 * sin(2*pi/period * tTail) + 10;

            targetTorque = [ ...
                oneSegment;
                oneSegment;
                oneSegment;
                tail];

        case '4sin'
            period1 = 2.3;
            period2 = 1.2;

            t1 = (0 : 1/fs : period1 - 1/fs)';
            t2 = (0 : 1/fs : period2 - 1/fs)';

            segment1 = ...
                2 * sin(2*pi/period1 * t1) + 10;

            segment2 = ...
                2 * sin(2*pi/period2 * t2) + 10;

            targetTorque = [ ...
                segment1;
                segment2;
                segment1;
                segment2];

    end

    assert(length(targetTorque) == 7*fs, ...
        'Generated target trajectory must be exactly 7 seconds long.');
end
