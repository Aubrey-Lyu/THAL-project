clear
home_dir = '/data/dian';%getuserdir; %'C:\Users\lvdia';

addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/Cohort_Organization'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/ThalamocoricalLoop-project/CCEP'))
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'))

dir_base = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL');

result_folder = fullfile(dir_base, 'CCEP' , 'results', 'explore2_redoPreProc_dataSort_Decomposition');
cd(result_folder)
plot_folder = fullfile(dir_base,'CCEP', 'Plots', 'explore2');
if ~exist(plot_folder, 'dir'); mkdir(plot_folder); end

% load data
metaT = readtable(fullfile(dir_base, 'COHORT', 'table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'));
load CCEP_all_flat_meanTr_cleaned
vars = metaT.Properties.VariableNames;
ttpIdx = find(contains(vars, 'pks_time_'));

% clean up data
metaT(badChan,:) = [];
zccep_clean(badChan,:) = [];
metaT.JP_label_out1 = ListSortAnatLabel_THAL(metaT.JP_label_out1, 1);
metaT.JP_label_in1 = ListSortAnatLabel_THAL(metaT.JP_label_in1, 1);
metaT.JP_label_out2 = ListSortAnatLabel_THAL(metaT.JP_label_out2, 1);
metaT.JP_label_in2 = ListSortAnatLabel_THAL(metaT.JP_label_in2, 1);
metaT.JP_label_out = ListSortAnatLabel_THAL(metaT.JP_label_out, 1);
metaT.JP_label_in = ListSortAnatLabel_THAL(metaT.JP_label_in, 1);

%% sort stim pairs
ROIs = {'antTH', 'midTH', 'pstTH'};
fromTHAL = ismember(metaT.JP_label_out, ROIs)  | metaT.sCrossBorder == 0;
toTHAL = ismember(metaT.JP_label_in, ROIs)  | metaT.rCrossBorder == 0;
% from THAL, exclude self-connection (to THAL)
%-----------------------------------------------
filterIdx = fromTHAL & ~toTHAL;

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_fromTHAL, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)
figure; 
n = 2000; %random select 2000 signals to check
zccep_ = zccep_fromTHAL;
randIdx = randsample(size(zccep_,1), n);
plot(zccep_(randIdx,:)')
%-----------------------------------------------
% to THAL, exclude self-connection (from THAL)
%-----------------------------------------------
filterIdx = ~fromTHAL & toTHAL;

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_toTHAL, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)
figure; 
n = 2000; %random select 2000 signals to check
zccep_ = zccep_toTHAL;
randIdx = randsample(size(zccep_,1), n);
plot(zccep_(randIdx,:)')
%-----------------------------------------------
% no THAL, mainly among cortical areas
%-----------------------------------------------
filterIdx = (~fromTHAL & ~toTHAL) & ...
    (~strcmp(metaT.JP_label_in1, 'EXCLUDE') | ~strcmp(metaT.JP_label_out1, 'EXCLUDE') | ...
    ~strcmp(metaT.JP_label_in2, 'EXCLUDE') | ~strcmp(metaT.JP_label_out2, 'EXCLUDE')); % stronger condition on cortex

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_noTHAL, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)
figure; 
n = 2000; %random select 2000 signals to check
zccep_ = zccep_noTHAL;
randIdx = randsample(size(zccep_,1), n);
plot(zccep_(randIdx,:)')
%-----------------------------------------------
% THAL self
%-----------------------------------------------
filterIdx = fromTHAL & toTHAL;

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_onlyTHAL, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)
figure; 
n = 2000; %random select 2000 signals to check
zccep_ = zccep_onlyTHAL;
randIdx = randsample(size(zccep_,1), n);
plot(zccep_(randIdx,:)')
%-----------------------------------------------

%% plot for sanity check
figure; 
n = 2000; %random select 2000 signals to check
zccep_ = zccep_toTHAL;
randIdx = randsample(size(zccep_,1), n);
plot(zccep_(randIdx,:)')


