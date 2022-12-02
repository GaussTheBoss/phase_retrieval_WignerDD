% This script generates all figures from Section 5 of the manuscript
%
% Note: Some of these scripts (especially those comparing the performance 
% of the proposed algorithms with other methods) may take a long time to 
% complete executing

clear; close all; clc

%% Section 5.1 - Algorithm 1 - Bandlimited Masks
cd bandlimited_masks
HIOER_ElbowPlot_Fig1
PostProcessing_Fig2a
StudyShifts_Fig2b
NoiseRobustness_Fig3a
ExecutionTime_Fig3b


%% Section 5.2 - Lemma 11 - Compactly Supported Masks
cd ../compactly_supported_masks
Comp_PrevPaper_Fig4a
StudyModes_Fig4b
NoiseRobustness_Fig5a
ExecutionTime_Fig5b


%% Section 5.3 - Algorithm 2 - Bandlimited Signals
cd ../bandlimited_signals
NoiseRobustness_Fig6a
ExecutionTime_Fig6b

