clear
home_dir = '/data/dian';%getuserdir; %'C:\Users\lvdia';

addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/Cohort_Organization'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/ThalamocoricalLoop-project/CCEP'))
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/tools'));

dir_base = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL');

result_folder = fullfile(dir_base, 'CCEP' , 'results', 'explore5_locked');
filpsign_folder = fullfile(result_folder , 'spectCCEP');
cd(result_folder)
plot_folder = fullfile(dir_base,'CCEP', 'Plots', 'explore5');
if ~exist(plot_folder, 'dir'); mkdir(plot_folder); end

%parpool(4)
% load data
metaT = readtable(fullfile(result_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone2.csv'));
load(fullfile(result_folder, 'CCEP_all_flat_meanTr_cleaned.mat'));
vars = metaT.Properties.VariableNames;
ttpIdx = find(contains(vars, 'pks_time_'));
% goodChan = ones(size(metaT,1),1);
% goodChan(badChan) = 0;
% prefilterIdx =  find(goodChan & ...
%         ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'NA'}) & ...
%         ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'NA'}))  ;
%
% % flip CCEP sign
% sublist = unique(metaT.subject);
% for ss = 1:length(sublist)
%     sb = sublist{ss};
%     sbID = find(strcmp(metaT.subject, sb));
%     [filterIdx, ~, ~] = intersect(sbID, prefilterIdx);
%
%     l = load(fullfile(filpsign_folder, ['flipCCEP_' sb, '.mat']));
%     zccep_ = l.toflip .* zccep_clean(filterIdx,:);
%
% end

% clean up data
metaT(badChan,:) = [];
zccep_clean(badChan,:) = [];
% metaT.JP_label_out1 = ListSortAnatLabel_THAL(metaT.JP_label_out1, 1);
% metaT.JP_label_in1 = ListSortAnatLabel_THAL(metaT.JP_label_in1, 1);
% metaT.JP_label_out2 = ListSortAnatLabel_THAL(metaT.JP_label_out2, 1);
% metaT.JP_label_in2 = ListSortAnatLabel_THAL(metaT.JP_label_in2, 1);
% metaT.JP_label_out = ListSortAnatLabel_THAL(metaT.JP_label_out, 1);
% metaT.JP_label_in = ListSortAnatLabel_THAL(metaT.JP_label_in, 1);

%% sort stim pairs
% prefilter
prefilterIdx =  ...
    ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'NA'}) & ...
    ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'NA'}) & ...
    metaT.sCrossBorder == 0 & metaT.rCrossBorder == 0;

ROIs = {'antTH', 'midTH', 'pstTH'};
fromTHAL = prefilterIdx & (ismember(metaT.JP_label_out, ROIs) );
toTHAL = prefilterIdx & (ismember(metaT.JP_label_in, ROIs) );
ipsi   = (metaT.MNIout_coord_1 .* metaT.MNIin_coord_1) >=0;
% from THAL, exclude self-connection (to THAL)
%-----------------------------------------------
% - ipsilateal
filterIdx = fromTHAL & ~toTHAL & ipsi;

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_fromTHAL_ipsi,~, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)

% - contralateral
filterIdx = fromTHAL & ~toTHAL & ~ipsi;

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_fromTHAL_contr,~, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)
% figure;
% n = 2000; %random select 2000 signals to check
% zccep_ = zccep_fromTHAL;
% randIdx = randsample(size(zccep_,1), n);
% plot(zccep_(randIdx,:)')
%-----------------------------------------------
% to THAL, exclude self-connection (from THAL)
%-----------------------------------------------
filterIdx = ~fromTHAL & toTHAL;

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_toTHAL,~, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)
% figure;
% n = 2000; %random select 2000 signals to check
% zccep_ = zccep_toTHAL;
% randIdx = randsample(size(zccep_,1), n);
% plot(zccep_(randIdx,:)')
%-----------------------------------------------
% no THAL, mainly among cortical areas
%-----------------------------------------------
ROIs =  {'CLT','AMY','BG'};

% - ipsilateal
filterIdx = (~fromTHAL & ~toTHAL) & ...
    ~ismember(metaT.JP_label_in, ROIs) & ipsi; % stronger condition on cortex

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_noTHAL_ipsi,~, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)

% - contralateral
filterIdx = (~fromTHAL & ~toTHAL) & ...
    ~ismember(metaT.JP_label_in, ROIs) & ~ipsi; % stronger condition on cortex

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_noTHAL_contr,~, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)
% figure;
% n = 2000; %random select 2000 signals to check
% zccep_ = zccep_noTHAL;
% randIdx = randsample(size(zccep_,1), n);
% plot(zccep_(randIdx,:)')
%-----------------------------------------------
% THAL self
%-----------------------------------------------
filterIdx = fromTHAL & toTHAL;

peakTimeMat = table2array(metaT(filterIdx ,ttpIdx));
[zccep_onlyTHAL,~, C] = flipSignROI(zccep_clean(filterIdx, :), peakTimeMat);
close;figure;plot(C)
% figure;
% n = 2000; %random select 2000 signals to check
% zccep_ = zccep_onlyTHAL;
% randIdx = randsample(size(zccep_,1), n);
% plot(zccep_(randIdx,:)')
%-----------------------------------------------

%% plot for sanity check
% figure;
% n = 2000; %random select 2000 signals to check
% zccep_ = zccep_toTHAL;
% randIdx = randsample(size(zccep_,1), n);
% plot(zccep_(randIdx,:)')
%
%% plot

clrs = [254 186 7; 92 179 204; 220 145 35; 18 110 130]/256; %hupohuang, biqing, fengfanhuang, yuqinlan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close
figure
%figure('Position', [3934         449        1695         667]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1)

dat = zccep_fromTHAL_ipsi(:, 200+[-20:800]);
time = [1:size(dat,2)] -20;
mn = mean(dat, 1, 'omitnan');
se = std(dat, 1,  'omitnan')./sqrt(sum(~isnan(dat))); % standard error;

shadedErrorBar(time, mn, se,...
    'lineProps',{'Color',clrs(1,:),'LineWidth', 1.35});
% xline(0); xlim([-20 800]);
% xlabel('time (ms)'); ylabel('z-score')
% title('thalamus-stim evoked potential in ipsilateral cortex')
%--------------------------------------
hold on

dat = zccep_noTHAL_ipsi(:, 200+[-20:800]);
time = [1:size(dat,2)] -20;
mn = mean(dat, 1, 'omitnan');
se = std(dat, 1,  'omitnan')./sqrt(sum(~isnan(dat))); % standard error;

shadedErrorBar(time, mn, se,...
    'lineProps',{'Color',clrs(2,:),'LineWidth', 1.35});
xline(0); xlim([-20 800]);
xlabel('time (ms)'); ylabel('z-score')
% title('cortex-stim evoked potential in ipsilateral cortex')
title('SPES-evoked potential in ipsilateral cortex')
legend({'from THAL (n=16056)', 'from COR (n=125126)'})
%--------------------------------------
subplot(2,1,2)

dat = zccep_fromTHAL_contr(:, 200+[-20:800]);
time = [1:size(dat,2)] -20;
mn = mean(dat, 1, 'omitnan');
se = std(dat, 1,  'omitnan')./sqrt(sum(~isnan(dat))); % standard error;

shadedErrorBar(time, mn, se,...
    'lineProps',{'Color',clrs(3,:),'LineWidth', 1.35});
% xline(0); xlim([-20 800]);
% xlabel('time (ms)'); ylabel('z-score')
% title('thalamus-stim evoked potential in contralateral cortex')
%--------------------------------------
% subplot(2,2,4)
hold on
dat = zccep_noTHAL_contr(:, 200+[-20:800]);
time = [1:size(dat,2)] -20;
mn = mean(dat, 1, 'omitnan');
se = std(dat, 1,  'omitnan')./sqrt(sum(~isnan(dat))); % standard error;

shadedErrorBar(time, mn, se,...
    'lineProps',{'Color',clrs(4,:),'LineWidth', 1.35});
xline(0); xlim([-20 800]);
xlabel('time (ms)'); ylabel('z-score')
% title('cortex-stim evoked potential in contralateral cortex')
title('SPES-evoked potential in contralateral cortex')
legend({'from THAL (n=13800)', 'from COR (n=78926)'})