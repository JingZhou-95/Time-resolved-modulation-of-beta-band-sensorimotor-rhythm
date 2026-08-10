# Inter- and intra-individual differences in β-band sensorimotor rhythm modulation and its relation to corticomuscular coherence in force-tracking tasks
 
# EEG–EMG Analysis Pipeline

This repository contains MATLAB code used for EEG, EMG, force, corticomuscular coherence (CMC), and SMR–CMC cross-correlation analyses.

The repository provides the complete analysis workflow used in the study. Raw experimental data are not included.

## Overview

The analysis pipeline includes:

* EEG preprocessing using EEGLAB
* EEG time-frequency analysis
* Motor-related beta desynchronization (MRBD)
* Post-movement beta synchronization (SRBS)
* ΔSRBS–MRBD calculation
* EEG–EMG corticomuscular coherence (CMC)
* Identification of the CMC-related 3-Hz frequency band
* Force performance analysis using RMSE
* SMR–CMC cross-correlation analysis
* Cross-task comparison of SMR–CMC temporal relationships

## Software Requirements

The code was developed in MATLAB R2023a.

Required software and toolboxes include:

* MATLAB
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox
* EEGLAB

Some preprocessing scripts also use custom filtering functions:

```text
bp.m    Band-pass filter
bs.m    Band-stop filter
```

EEGLAB must be installed separately and added to the MATLAB path before running the preprocessing scripts.

## Repository Structure

### `EEG_preprocessing.m`

Preprocesses the continuous EEG recordings.

Main steps:

1. Extract EEG, torque, and EMG signals.
2. Remove linear trends from EEG signals.
3. Apply band-pass filtering.
4. Remove 50-Hz power-line noise and harmonics.
5. Generate task-onset event markers.
6. Import the data into EEGLAB.
7. Perform manual artifact rejection.
8. Perform Independent Component Analysis (ICA).
9. Remove artifact-related ICA components.
10. Extract task-locked EEG epochs.

The resulting continuous and epoched EEG data are saved for subsequent analyses.

---

### `SMR_CMC_Force_analysis.m`

Performs the main EEG, EMG, CMC, and torque-performance analyses.

Supported task conditions are:

```text
10pct
15pct
St
1sin
3sin
4sin
```

The task condition is specified using:

```matlab
taskType = '10pct';
```

#### Experiment 1

* `10pct`: constant target at 10 %MVC
* `15pct`: constant target at 15 %MVC

#### Experiment 2

* `St`: constant target at 10 %MVC
* `1sin`: sinusoidal target with a 7-s period
* `3sin`: repeated sinusoidal target with a 2.3-s period
* `4sin`: consecutive sinusoidal segments with periods of 2.3, 1.2, 2.3, and 1.2 s

Main analyses include:

1. Cz surface-Laplacian EEG calculation.
2. EEG time-frequency power analysis using 1-s windows shifted every 50 ms.
3. Beta-band power analysis.
4. Calculation of MRBD.
5. Calculation of SRBS.
6. Calculation of ΔSRBS–MRBD:

```text
ΔSRBS–MRBD = SRBS − MRBD
```

7. EEG–TA corticomuscular coherence using both rectified and non-rectified EMG.
8. Identification of the 3-Hz frequency band showing maximum CMC.
9. EEG power analysis within the CMC-selected frequency band.
10. Force-performance analysis.

Force performance is evaluated using the middle 5 s of each trial.

RMSE is calculated as:

```text
RMSE = sqrt(mean((observed force − target force)^2))
```

RMSE is first calculated separately for each trial and then averaged across trials for each participant and task condition.

---

### `Exp1_SMR_CMC_crosscorrelation.m`

Performs pooled SMR–CMC cross-correlation analysis for the 10% and 15% tasks.

Input matrices contain SMR and CMC time courses for each task observation.

The analysis is restricted to the 0–7 s contraction period.

Cross-correlation is calculated as:

```matlab
xcorr(CMC, SMR)
```

The lag convention is:

```text
Positive lag: SMR leads CMC
Negative lag: CMC leads SMR
```

Main steps include:

1. Calculate an SMR–CMC cross-correlation function for each observation.
2. Perform observation-level permutation tests.
3. Identify observations with significant cross-correlations.
4. Calculate the pooled group-average cross-correlation function.
5. Perform cluster-based permutation testing across lag values.
6. Determine the group-level peak lag.
7. Compare CMC-positive and CMC-negative observations at the peak lag using a linear mixed-effects model.

---

### `Exp2_SMR_CMC_crosscorrelation.m`

Compares SMR–CMC temporal relationships across the four Experiment 2 tasks:

```text
St
1sin
3sin
4sin
```

The same participants are analyzed across all four conditions.

Main steps include:

1. Calculate the SMR–CMC cross-correlation function for each participant and task.
2. Calculate the group-average cross-correlation function for each task.
3. Perform one-tailed cluster-based permutation tests for positive cross-correlations.
4. Determine the task-specific peak lag.
5. Extract each participant's cross-correlation coefficient at the corresponding task-specific peak lag.
6. Compare tasks using the Friedman test.
7. Perform Dunn's post-hoc multiple comparisons with correction for multiple comparisons.
8. Examine individual peak-lag distributions across tasks.

---


## Data Format

Raw data are not distributed with this repository.

The original recordings contain synchronized EEG, EMG, and torque signals. The main analysis scripts assume MATLAB matrices with channels organized according to the corresponding script documentation.

For example, the main EEG/CMC analysis uses:

```text
Column 1       Force
Column 2       Cz
Column 3       FCz
Column 4       CPz
Column 5       C1
Column 6       C2
Column 31      Tibialis anterior EMG
```

Users applying the code to their own datasets should adapt the channel indices and experimental parameters as necessary.

## Data Availability

The raw data used in this study are not publicly available due to ethical restrictions.

The analysis code is provided to document the complete signal-processing and statistical-analysis workflow.

Researchers interested in reproducing the analyses may apply the code to appropriately formatted EEG, EMG, and force datasets.
